use crate::columns::{NUM_COLS, SBOX_DEGREE, STATE_SIZE};
use crate::poseidon2::{MAT_DIAG8_M_1, RC8};
use ark_ff::{BigInteger, PrimeField};
use num::BigUint;
use plonky2::field::extension::{Extendable, FieldExtension};
use plonky2::field::packed::PackedField;
use plonky2::field::polynomial::PolynomialValues;
use plonky2::field::types::Field;
use plonky2::hash::hash_types::RichField;
use plonky2::plonk::circuit_builder::CircuitBuilder;
use starky::constraint_consumer::{ConstraintConsumer, RecursiveConstraintConsumer};
use starky::stark::Stark;
use starky::vars::{StarkEvaluationTargets, StarkEvaluationVars};
use std::marker::PhantomData;
use zkhash::fields::goldilocks::FpGoldiLocks;

const M4: [[usize; 4]; 4] = [
    [5, 7, 1, 3], //
    [4, 6, 1, 1], //
    [1, 3, 5, 7], //
    [1, 1, 4, 6],
];

fn scalar_to_fe<F: RichField + Extendable<D>, const D: usize, FE, const D2: usize, PF: PrimeField>(
    scalar: PF,
) -> FE
where
    FE: FieldExtension<D2, BaseField = F>,
{
    FE::from_canonical_u64(
        F::from_noncanonical_biguint(BigUint::from_bytes_le(&scalar.into_bigint().to_bytes_le()))
            .to_canonical_u64(),
    )
}

// linear layer (degree = 1)
fn matmul_external8_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
>(
    state: &[P; 8],
) -> [P; 8]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 8);
    let mut out = [P::ZEROS; 8];

    for i in 0..8 {
        for j in 0..4 {
            out[i] = if i < 4 {
                out[i] + state[j] * FE::from_canonical_usize(M4[i][j])
            } else {
                out[i] + state[j + 4] * FE::from_canonical_usize(M4[i - 4][j])
            }
        }
    }

    let mut stored = [P::ZEROS; 4];
    for i in 0..4 {
        stored[i] = out[i] + out[4 + i];
    }
    for i in 0..8 {
        out[i] = out[i] + stored[i % 4];
    }

    out
}

fn matmul_internal8_constraints<
    F: RichField + Extendable<D>,
    const D: usize,
    FE,
    P,
    const D2: usize,
>(
    state: &[P; 8],
) -> [P; 8]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 8);
    let mut out = [P::ZEROS; 8];
    let mut sum = P::ZEROS;

    for i in 0..STATE_SIZE {
        sum = sum + state[i];
    }

    for i in 0..STATE_SIZE {
        out[i] = state[i] * scalar_to_fe::<F, D, FE, D2, FpGoldiLocks>(MAT_DIAG8_M_1[i]);
        out[i] = out[i] + sum;
    }

    out
}

// degree: 1
fn add_rc_constraints<F: RichField + Extendable<D>, const D: usize, FE, P, const D2: usize>(
    state: &[P; 8],
    r: usize,
) -> [P; 8]
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 8);
    let mut out = [P::ZEROS; 8];

    for i in 0..8 {
        out[i] = state[i] + scalar_to_fe::<F, D, FE, D2, FpGoldiLocks>(RC8[r][i]);
    }

    out
}

// degree: 7
fn sbox_p_constraints<F: RichField + Extendable<D>, const D: usize, FE, P, const D2: usize>(
    state: &P,
) -> P
where
    FE: FieldExtension<D2, BaseField = F>,
    P: PackedField<Scalar = FE>,
{
    assert_eq!(STATE_SIZE, 8);
    let mut out = P::ONES;

    for _ in 0..SBOX_DEGREE {
        out = out * *state;
    }

    out
}

#[derive(Copy, Clone, Default)]
pub struct Poseidon2Stark<F, const D: usize> {
    pub _f: PhantomData<F>,
}

impl<F: RichField + Extendable<D>, const D: usize> Stark<F, D> for Poseidon2Stark<F, D> {
    const COLUMNS: usize = NUM_COLS;
    const PUBLIC_INPUTS: usize = 0;

    fn eval_packed_generic<FE, P, const D2: usize>(
        &self,
        vars: StarkEvaluationVars<FE, P, { Self::COLUMNS }, { Self::PUBLIC_INPUTS }>,
        _yield_constr: &mut ConstraintConsumer<P>,
    ) where
        FE: FieldExtension<D2, BaseField = F>,
        P: PackedField<Scalar = FE>,
    {
        let lv = vars.local_values;
        let _ = matmul_external8_constraints(lv[0..8].try_into().unwrap());
    }

    fn constraint_degree(&self) -> usize {
        unimplemented!()
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
    use crate::columns::{
        COL_1ST_FULLROUND_STATE_START, COL_PARTIAL_ROUND_STATE_START, NUM_COLS, ROUNDS_F, ROUNDS_P,
        STATE_SIZE,
    };
    use crate::generation::{generate_poseidon2_trace, Row};
    use crate::poseidon2::RC8;
    use crate::stark::FpGoldiLocks;
    use crate::stark::{scalar_to_fe, trace_to_poly_values, Poseidon2Stark};
    use anyhow::Result;
    use plonky2::field::extension::{Extendable, FieldExtension};
    use plonky2::field::packed::PackedField;
    use plonky2::field::types::Field;
    use plonky2::hash::hash_types::RichField;
    use plonky2::plonk::circuit_builder::CircuitBuilder;
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    use plonky2::util::timing::TimingTree;
    use starky::config::StarkConfig;
    use starky::constraint_consumer::{ConstraintConsumer, RecursiveConstraintConsumer};
    use starky::prover::prove;
    use starky::stark::Stark;
    use starky::stark_testing::test_stark_low_degree;
    use starky::vars::{StarkEvaluationTargets, StarkEvaluationVars};
    use starky::verifier::verify_stark_proof;
    use std::marker::PhantomData;

    #[derive(Copy, Clone, Default)]
    pub struct PoseidonTestStark<F, const D: usize> {
        pub _f: PhantomData<F>,
    }
    impl<F: RichField + Extendable<D>, const D: usize> Stark<F, D> for PoseidonTestStark<F, D> {
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
            let mut state = super::matmul_external8_constraints(lv[0..8].try_into().unwrap());

            // first full rounds
            for r in 0..ROUNDS_F {
                state = super::add_rc_constraints(&state, r);
                for i in 0..STATE_SIZE {
                    state[i] = super::sbox_p_constraints(&state[i]);
                }
                state = super::matmul_external8_constraints(&state);
                for i in 0..STATE_SIZE {
                    yield_constr.constraint(
                        state[i] - lv[COL_1ST_FULLROUND_STATE_START + r * STATE_SIZE + i],
                    );
                    state[i] = lv[COL_1ST_FULLROUND_STATE_START + r * STATE_SIZE + i];
                }
            }

            // partial rounds
            for i in 0..ROUNDS_P {
                let r = ROUNDS_F + i;
                state[0] = state[0] + scalar_to_fe::<F, D, FE, D2, FpGoldiLocks>(RC8[r][0]);
                state[0] = super::sbox_p_constraints(&state[0]);
                state = super::matmul_internal8_constraints(&state);
                yield_constr.constraint(state[0] - lv[COL_PARTIAL_ROUND_STATE_START + i]);
                state[0] = lv[COL_PARTIAL_ROUND_STATE_START + i];
            }

            for i in 0..STATE_SIZE {
                yield_constr.constraint(state[i] - lv[8 + i]);
            }
        }

        fn eval_ext_circuit(
            &self,
            _builder: &mut CircuitBuilder<F, D>,
            _vars: StarkEvaluationTargets<D, { Self::COLUMNS }, { Self::PUBLIC_INPUTS }>,
            _yield_constr: &mut RecursiveConstraintConsumer<F, D>,
        ) {
            unimplemented!()
        }

        fn constraint_degree(&self) -> usize {
            7
        }
    }

    #[test]
    fn poseidon2_constraints() -> Result<()> {
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;
        type S = PoseidonTestStark<F, D>;
        let mut config = StarkConfig::standard_fast_config();
        config.fri_config.cap_height = 0;
        config.fri_config.rate_bits = 3;

        let num_rows = 1;
        let mut step_rows = Vec::with_capacity(num_rows);
        for _ in 0..num_rows {
            let preimage = (0..STATE_SIZE)
                .map(|i| F::from_canonical_usize(i + 1))
                .collect::<Vec<_>>();
            step_rows.push(Row {
                preimage: preimage.try_into().unwrap(),
            });
        }

        let stark = S::default();
        let mut trace = generate_poseidon2_trace(step_rows);
        for row in 0..num_rows {
            trace[8][row] = F::from_canonical_u64(8819548283525844653);
            trace[9][row] = F::from_canonical_u64(12228992105858652070);
            trace[10][row] = F::from_canonical_u64(4861132991556502187);
            trace[11][row] = F::from_canonical_u64(5587681080442376794);
            trace[12][row] = F::from_canonical_u64(1339891673634421267);
            trace[13][row] = F::from_canonical_u64(13568472719586988937);
            trace[14][row] = F::from_canonical_u64(16342764561928125140);
            trace[15][row] = F::from_canonical_u64(13652518956506335588);
        }
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
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;
        type S = PoseidonTestStark<F, D>;

        let num_rows = 1 << 5;
        let stark = S::default();
        test_stark_low_degree(stark)
    }
}
