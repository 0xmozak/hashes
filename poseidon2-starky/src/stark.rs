use crate::columns::{NUM_COLS, SBOX_DEGREE, STATE_SIZE};
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
use zkhash::ark_ff::{BigInteger, PrimeField};
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::poseidon2::poseidon2_instance_goldilocks::RC8;

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
    use crate::columns::{NUM_COLS, STATE_SIZE};
    use crate::generation::{generate_poseidon2_trace, Row};
    use crate::stark::trace_to_poly_values;
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
            let mut out = super::matmul_external8_constraints(lv[0..8].try_into().unwrap());
            out = super::add_rc_constraints(&out, 0);
            for i in 0..8 {
                out[i] = super::sbox_p_constraints(&out[i]);
            }
            for i in 0..8 {
                yield_constr.constraint(out[i] - lv[8 + i]);
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
            trace[8][row] = F::from_canonical_u64(4924480480635527917);
            trace[9][row] = F::from_canonical_u64(7328260706673212971);
            trace[10][row] = F::from_canonical_u64(12769208979081050669);
            trace[11][row] = F::from_canonical_u64(226507612748720379);
            trace[12][row] = F::from_canonical_u64(1452606391015137905);
            trace[13][row] = F::from_canonical_u64(13954454011491297451);
            trace[14][row] = F::from_canonical_u64(12425154049474751493);
            trace[15][row] = F::from_canonical_u64(14960278433617001970);
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
}
