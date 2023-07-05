use crate::columns::{
    COL_1ST_FULLROUND_STATE_START, COL_2ND_FULLROUND_STATE_START, COL_INPUT_START,
    COL_OUTPUT_START, COL_PARTIAL_ROUND_END_STATE_START, COL_PARTIAL_ROUND_STATE_START, NUM_COLS,
    ROUNDS_F, ROUNDS_P, STATE_SIZE,
};
use crate::poseidon2::Poseidon2;
use crate::poseidon2::POSEIDON2_GOLDILOCKS_8_PARAMS;
use ark_ff::{BigInteger, PrimeField};
use num::bigint::BigUint;
use plonky2::hash::hash_types::RichField;
use std::convert::TryInto;
use std::ops::AddAssign;
use zkhash::fields::goldilocks::FpGoldiLocks;

pub(crate) fn scalar_to_field<RF: RichField, PF: PrimeField>(scalar: &PF) -> RF {
    RF::from_noncanonical_biguint(BigUint::from_bytes_le(&scalar.into_bigint().to_bytes_le()))
}

pub(crate) fn scalar_to_field_vec<RF: RichField, PF: PrimeField>(scalar: &[PF]) -> Vec<RF> {
    scalar.iter().map(|s| scalar_to_field(s)).collect()
}

pub(crate) fn field_to_scalar<RF: RichField, PF: PrimeField>(field: &RF) -> PF {
    PF::from_le_bytes_mod_order(&field.to_canonical_biguint().to_bytes_le())
}

pub(crate) fn field_to_scalar_vec<RF: RichField, PF: PrimeField>(field: &[RF]) -> Vec<PF> {
    field.iter().map(|f| field_to_scalar(f)).collect()
}

// Represent a row of the preimage
#[derive(Debug, Clone, Default)]
pub struct Row<Field: RichField> {
    pub(crate) preimage: [Field; STATE_SIZE],
}

/// Pad the trace to a power of 2.
#[must_use]
fn pad_trace<F: RichField>(mut trace: Vec<Vec<F>>) -> Vec<Vec<F>> {
    let ext_trace_len = trace[0].len().next_power_of_two();

    // All columns have their last value duplicated.
    for row in &mut trace {
        row.resize(ext_trace_len, *row.last().unwrap());
    }

    trace
}

fn generate_1st_full_round_state<Field: RichField>(
    preimage: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);

    let mut current_state = field_to_scalar_vec(preimage);

    // Linear layer at beginning
    instance.matmul_external(&mut current_state);

    for r in 0..instance.params.rounds_f_beginning {
        current_state = instance.add_rc(&current_state, &instance.params.round_constants[r]);
        current_state = instance.sbox(&current_state);
        instance.matmul_external(&mut current_state);
        outputs.push(scalar_to_field_vec(&current_state).try_into().unwrap());
    }

    outputs
}

fn generate_partial_round_state<Field: RichField>(
    last_rount_output: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);

    let mut current_state: Vec<FpGoldiLocks> = field_to_scalar_vec(last_rount_output);

    let p_end = instance.params.rounds_f_beginning + instance.params.rounds_p;
    for r in instance.params.rounds_f_beginning..p_end {
        current_state[0].add_assign(&instance.params.round_constants[r][0]);
        current_state[0] = instance.sbox_p(&current_state[0]);
        instance.matmul_internal(&mut current_state, &instance.params.mat_internal_diag_m_1);
        outputs.push(scalar_to_field_vec(&current_state).try_into().unwrap());
    }

    outputs
}

fn generate_2st_full_round_state<Field: RichField>(
    last_rount_output: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);

    let mut current_state = field_to_scalar_vec(last_rount_output);

    let p_end = instance.params.rounds_f_beginning + instance.params.rounds_p;
    for r in p_end..instance.params.rounds {
        current_state = instance.add_rc(&current_state, &instance.params.round_constants[r]);
        current_state = instance.sbox(&current_state);
        instance.matmul_external(&mut current_state);
        outputs.push(scalar_to_field_vec(&current_state).try_into().unwrap());
    }

    outputs
}

// Generate the outputs for a given preimage
fn generate_outputs<Field: RichField>(preimage: &[Field; STATE_SIZE]) -> [Field; STATE_SIZE] {
    let mut outputs = [Field::ZERO; STATE_SIZE];
    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);

    let input = field_to_scalar_vec(preimage);
    let perm = instance.permutation(&input);

    for i in 0..STATE_SIZE {
        outputs[i] = scalar_to_field(&perm[i]);
    }
    outputs
}

// Function to generate the Poseidon2 trace
pub fn generate_poseidon2_trace<F: RichField>(step_rows: Vec<Row<F>>) -> [Vec<F>; NUM_COLS] {
    let trace_len = step_rows.len();
    let mut trace: Vec<Vec<F>> = vec![vec![F::ZERO; trace_len]; NUM_COLS];

    for (i, row) in step_rows.iter().enumerate() {
        for j in 0..STATE_SIZE {
            trace[COL_INPUT_START + j][i] = row.preimage[j];
        }
        let outputs = generate_outputs(&row.preimage);
        for j in 0..STATE_SIZE {
            trace[COL_OUTPUT_START + j][i] = outputs[j];
        }

        // Generate the full round states
        let first_full_round_state = generate_1st_full_round_state(&row.preimage);
        let partial_round_state = generate_partial_round_state(
            first_full_round_state.last().unwrap().try_into().unwrap(),
        );
        let second_full_round_state =
            generate_2st_full_round_state(partial_round_state.last().unwrap().try_into().unwrap());
        for j in 0..ROUNDS_F {
            for k in 0..STATE_SIZE {
                trace[COL_1ST_FULLROUND_STATE_START + j * STATE_SIZE + k][i] =
                    first_full_round_state[j][k];
                trace[COL_2ND_FULLROUND_STATE_START + j * STATE_SIZE + k][i] =
                    second_full_round_state[j][k];
            }
        }
        for j in 0..ROUNDS_P {
            trace[COL_PARTIAL_ROUND_STATE_START + j][i] = partial_round_state[j][0];
        }
        for j in 0..STATE_SIZE {
            trace[COL_PARTIAL_ROUND_END_STATE_START + j][i] = partial_round_state[ROUNDS_P - 1][j];
        }
    }

    trace = pad_trace(trace);
    trace.try_into().unwrap_or_else(|v: Vec<Vec<F>>| {
        panic!(
            "Expected a Vec of length {} but it was {}",
            NUM_COLS,
            v.len()
        )
    })
}

#[cfg(test)]
mod test {
    use crate::columns::{COL_OUTPUT_START, STATE_SIZE};
    use crate::generation::{
        field_to_scalar_vec, generate_1st_full_round_state, generate_2st_full_round_state,
        generate_outputs, generate_partial_round_state, scalar_to_field_vec, Row,
    };
    use plonky2::field::goldilocks_field::GoldilocksField;
    use plonky2::field::types::{PrimeField64, Sample};
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    use zkhash::fields::goldilocks::FpGoldiLocks;
    use zkhash::poseidon2::poseidon2::Poseidon2;
    use zkhash::poseidon2::poseidon2_instance_goldilocks::POSEIDON2_GOLDILOCKS_8_PARAMS;

    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    #[test]
    fn field_vec_conversion() {
        let fields: Vec<GoldilocksField> = (0..16).map(|_| GoldilocksField::rand()).collect();
        let scalars: Vec<FpGoldiLocks> = field_to_scalar_vec(&fields);
        let fields2: Vec<GoldilocksField> = scalar_to_field_vec(&scalars);
        assert_eq!(fields, fields2);
    }

    #[test]
    fn rounds_generation() {
        let preimage = (0..STATE_SIZE).map(|_| F::rand()).collect::<Vec<_>>();
        let output0: Vec<[F; STATE_SIZE]> =
            generate_1st_full_round_state(&preimage.clone().try_into().unwrap());
        let output1: Vec<[F; STATE_SIZE]> =
            generate_partial_round_state(output0.last().unwrap().try_into().unwrap());
        let output2: Vec<[F; STATE_SIZE]> =
            generate_2st_full_round_state(output1.last().unwrap().try_into().unwrap());
        let expected_output = generate_outputs(&preimage.try_into().unwrap());
        assert_eq!(expected_output, *output2.last().unwrap());
    }

    #[test]
    fn generate_poseidon2_trace() {
        let num_rows = 12;
        let mut step_rows = Vec::with_capacity(num_rows);

        for _ in 0..num_rows {
            let preimage = (0..STATE_SIZE).map(|_| F::rand()).collect::<Vec<_>>();
            step_rows.push(Row {
                preimage: preimage.try_into().unwrap(),
            });
        }

        let trace = super::generate_poseidon2_trace(step_rows.clone());
        for i in 0..trace.len() {
            assert_eq!(trace[i].len(), 16);
        }

        let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
        for i in 0..num_rows {
            let input = step_rows[i]
                .preimage
                .iter()
                .map(|x| FpGoldiLocks::from(x.to_canonical_u64()))
                .collect::<Vec<_>>();
            let perm = instance.permutation(&input);

            for j in 0..STATE_SIZE {
                let expected_val =
                    FpGoldiLocks::from(trace[COL_OUTPUT_START + j][i].to_canonical_u64());
                assert_eq!(
                    perm[j], expected_val,
                    "Mismatch at row {}, position {}",
                    i, j
                );
            }
        }
    }
}
