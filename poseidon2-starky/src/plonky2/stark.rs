use crate::plonky2::columns::{
    COL_1ST_FULLROUND_STATE_START, COL_2ND_FULLROUND_STATE_START,
    COL_PARTIAL_ROUND_END_STATE_START, COL_PARTIAL_ROUND_STATE_START, NUM_COLS, ROUNDS_F, ROUNDS_P,
    SBOX_DEGREE, STATE_SIZE,
};
use plonky2::field::extension::{Extendable, FieldExtension};
use plonky2::field::packed::PackedField;
use plonky2::field::polynomial::PolynomialValues;
use plonky2::field::types::Field;
use plonky2::hash::hash_types::RichField;
use plonky2::hash::poseidon2::Poseidon2;
use plonky2::plonk::circuit_builder::CircuitBuilder;
use starky::constraint_consumer::{ConstraintConsumer, RecursiveConstraintConsumer};
use starky::stark::Stark;
use starky::vars::{StarkEvaluationTargets, StarkEvaluationVars};
use std::fmt::Display;
use std::marker::PhantomData;

// degree: 1
fn add_rc_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
    const STATE_SIZE: usize,
>(
    state: &[P; STATE_SIZE],
    r: usize,
) -> [P; STATE_SIZE]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 12);
    let mut out = [P::ZEROS; STATE_SIZE];

    for i in 0..STATE_SIZE {
        out[i] =
            state[i] + FE::from_basefield(F::from_canonical_u64(<F as Poseidon2>::RC12[r + i]));
    }

    out
}

// degree: SBOX_DEGREE (7)
fn sbox_p_constraints<F: RichField + Extendable<D>, const D: usize, FE, P, const D2: usize>(
    state: &P,
) -> P
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    let mut out = P::ONES;

    for _ in 0..SBOX_DEGREE {
        out = out.mul(*state);
    }

    out
}

fn matmul_m4_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
    const STATE_SIZE: usize,
>(
    state: &[P; STATE_SIZE],
) -> [P; STATE_SIZE]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    // input x = (x0, x1, x2, x3)
    assert_eq!(STATE_SIZE, 12);
    let mut out = [P::ZEROS; STATE_SIZE];
    let t4 = STATE_SIZE / 4;
    for i in 0..t4 {
        let start_index = i * 4;
        // t0 = x0 + x1
        let t_0 = state[start_index] + state[start_index + 1];

        // t1 = x2 + x3
        let t_1 = state[start_index + 2] + state[start_index + 3];

        // 2x1
        let x1_2 = state[start_index + 1].mul(FE::TWO);
        // 2x3
        let x3_2 = state[start_index + 3].mul(FE::TWO);
        let four = FE::TWO + FE::TWO;

        // t2 = 2x1 + t1
        let t_2 = x1_2 + t_1;

        // t3 = 2x3 + t0
        let t_3 = x3_2 + t_0;

        // t4 = 4t1 + t3
        let t_4 = t_3 + t_1.mul(four);

        // t5 = 4t0 + t2
        let t_5 = t_2 + t_0.mul(four);

        // t6 = t3 + t5
        let t_6 = t_3 + t_5;

        // t7 = t2 + t4
        let t_7 = t_2 + t_4;

        out[start_index] = t_6;
        out[start_index + 1] = t_5;
        out[start_index + 2] = t_7;
        out[start_index + 3] = t_4;
    }
    out
}

fn matmul_external12_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
    const STATE_SIZE: usize,
>(
    state: &[P; STATE_SIZE],
) -> [P; STATE_SIZE]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 12);
    let mut out = [P::ZEROS; STATE_SIZE];
    let updated_state = matmul_m4_constraints(state);

    let t4 = STATE_SIZE / 4;
    let mut stored = [P::ZEROS; 4];

    for l in 0..4 {
        stored[l] = updated_state[l];
        for j in 1..t4 {
            stored[l] += updated_state[4 * j + l];
        }
    }
    for i in 0..STATE_SIZE {
        out[i] = updated_state[i].add(stored[i % 4]);
    }
    out
}

// degree: 1
fn matmul_internal12_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
    const STATE_SIZE: usize,
>(
    state: &[P; STATE_SIZE],
) -> [P; STATE_SIZE]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 12);
    let mut out = [P::ZEROS; STATE_SIZE];
    let mut sum = P::ZEROS;

    for item in state {
        sum += *item;
    }

    for i in 0..STATE_SIZE {
        out[i] = state[i]
            * FE::from_basefield(F::from_canonical_u64(
                <F as Poseidon2>::MAT_DIAG12_M_1[i] - 1,
            ));
        out[i] += sum;
    }

    out
}

#[derive(Copy, Clone, Default)]
#[allow(clippy::module_name_repetitions)]
pub struct Poseidon2_12Stark<F, const D: usize> {
    pub _f: PhantomData<F>,
}

impl<F, const D: usize> Display for Poseidon2_12Stark<F, D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Poseidon2_12Stark")
    }
}

impl<F: RichField + Extendable<D>, const D: usize> Stark<F, D> for Poseidon2_12Stark<F, D> {
    const COLUMNS: usize = NUM_COLS;
    const PUBLIC_INPUTS: usize = 0;
    fn eval_packed_generic<FE, P, const D2: usize>(
        &self,
        vars: StarkEvaluationVars<FE, P, { Self::COLUMNS }, { Self::PUBLIC_INPUTS }>,
        yield_constr: &mut ConstraintConsumer<P>,
    ) where
        FE: FieldExtension<D2, BaseField = F>,
        P: PackedField<Scalar = FE>,
    {
        let lv = vars.local_values;
        let mut state: [P; STATE_SIZE] =
            matmul_external12_constraints(lv[0..STATE_SIZE].try_into().unwrap());
        // first full rounds
        for r in 0..(ROUNDS_F / 2) {
            state = add_rc_constraints(&state, r);
            #[allow(clippy::needless_range_loop)]
            for i in 0..STATE_SIZE {
                state[i] = sbox_p_constraints(&state[i]);
            }
            state = matmul_external12_constraints(&state);
            for i in 0..STATE_SIZE {
                yield_constr
                    .constraint(state[i] - lv[COL_1ST_FULLROUND_STATE_START + r * STATE_SIZE + i]);
                state[i] = lv[COL_1ST_FULLROUND_STATE_START + r * STATE_SIZE + i];
            }
        }

        // partial rounds
        for i in 0..ROUNDS_P {
            state[0] += FE::from_basefield(F::from_canonical_u64(<F as Poseidon2>::RC12_MID[i]));
            state[0] = sbox_p_constraints(&state[0]);
            state = matmul_internal12_constraints(&state);
            yield_constr.constraint(state[0] - lv[COL_PARTIAL_ROUND_STATE_START + i]);
            state[0] = lv[COL_PARTIAL_ROUND_STATE_START + i];
        }

        // // the state before last full rounds
        for i in 0..STATE_SIZE {
            yield_constr.constraint(state[i] - lv[COL_PARTIAL_ROUND_END_STATE_START + i]);
            state[i] = lv[COL_PARTIAL_ROUND_END_STATE_START + i];
        }

        // // last full rounds
        for i in 0..(ROUNDS_F / 2) {
            let r = (ROUNDS_F / 2) + i;
            state = add_rc_constraints(&state, r);
            #[allow(clippy::needless_range_loop)]
            for j in 0..STATE_SIZE {
                state[j] = sbox_p_constraints(&state[j]);
            }
            state = matmul_external12_constraints(&state);
            for j in 0..STATE_SIZE {
                yield_constr
                    .constraint(state[j] - lv[COL_2ND_FULLROUND_STATE_START + i * STATE_SIZE + j]);
                state[j] = lv[COL_2ND_FULLROUND_STATE_START + i * STATE_SIZE + j];
            }
        }
    }
    fn constraint_degree(&self) -> usize {
        7
    }

    fn eval_ext_circuit(
        &self,
        _builder: &mut CircuitBuilder<F, D>,
        _vars: StarkEvaluationTargets<D, { Self::COLUMNS }, { Self::PUBLIC_INPUTS }>,
        _yield_constr: &mut RecursiveConstraintConsumer<F, D>,
    ) {
        unimplemented!()
    }
}

pub fn trace_to_poly_values<F: Field, const COLUMNS: usize>(
    trace: [Vec<F>; COLUMNS],
) -> Vec<PolynomialValues<F>> {
    trace.into_iter().map(PolynomialValues::new).collect()
}

#[cfg(test)]
mod tests {
    use crate::plonky2::columns::STATE_SIZE;
    use crate::plonky2::generation::{generate_poseidon2_trace, Row};
    use crate::plonky2::stark::{trace_to_poly_values, Poseidon2_12Stark};
    use anyhow::Result;
    use plonky2::field::types::Sample;
    use plonky2::plonk::config::{GenericConfig, Poseidon2GoldilocksConfig};
    use plonky2::util::timing::TimingTree;
    use starky::config::StarkConfig;
    use starky::prover::prove;
    use starky::stark_testing::test_stark_low_degree;
    use starky::verifier::verify_stark_proof;

    const D: usize = 2;
    type C = Poseidon2GoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;
    type S = Poseidon2_12Stark<F, D>;

    #[test]
    fn poseidon2_constraints() -> Result<()> {
        let mut config = StarkConfig::standard_fast_config();
        config.fri_config.cap_height = 0;
        config.fri_config.rate_bits = 3; // to meet the constraint degree bound

        let num_rows = 12;
        let mut step_rows = Vec::with_capacity(num_rows);
        for _ in 0..num_rows {
            let preimage = (0..STATE_SIZE).map(|_| F::rand()).collect::<Vec<_>>();
            step_rows.push(Row {
                preimage: preimage.try_into().unwrap(),
            });
        }

        let stark = S::default();
        let trace = generate_poseidon2_trace(&step_rows);
        let trace_poly_values = trace_to_poly_values(trace);

        let proof = prove::<F, C, S, D>(
            stark,
            &config,
            trace_poly_values,
            [],
            &mut TimingTree::default(),
        )?;
        verify_stark_proof(stark, proof, &config)
    }

    #[test]
    fn poseidon2_stark_degree() -> Result<()> {
        let stark = S::default();
        test_stark_low_degree(stark)
    }
}
