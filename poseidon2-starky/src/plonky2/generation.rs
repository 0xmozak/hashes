use crate::plonky2::columns::{
    COL_1ST_FULLROUND_STATE_START, COL_2ND_FULLROUND_STATE_START, COL_INPUT_START,
    COL_OUTPUT_START, COL_PARTIAL_ROUND_END_STATE_START, COL_PARTIAL_ROUND_STATE_START, NUM_COLS,
    ROUNDS_F, ROUNDS_P, STATE_SIZE,
};
use plonky2::hash::{
    hash_types::RichField,
    poseidon2::{Poseidon2, WIDTH},
};

use super::columns::COL_IS_EXE;

// Represent a row of the preimage
#[derive(Debug, Clone, Default)]
pub struct Row<Field: RichField> {
    pub preimage: [Field; STATE_SIZE],
}

/// Pad the trace to a power of 2.
#[must_use]
fn pad_trace<F: RichField>(mut trace: Vec<Vec<F>>) -> Vec<Vec<F>> {
    let original_len = trace[0].len();
    let ext_trace_len = original_len.next_power_of_two();

    // All columns have their last value duplicated.
    for row in &mut trace {
        row.resize(ext_trace_len, *row.last().unwrap());
    }
    // Set COL_IS_EXE to ZERO
    for i in original_len..ext_trace_len {
        trace[COL_IS_EXE][i] = F::ZERO;
    }

    trace
}

fn generate_1st_full_round_state<Field: RichField>(
    preimage: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    assert_eq!(STATE_SIZE, WIDTH);
    let mut current_state = *preimage;

    // Linear layer at start
    Field::matmul_external(&mut current_state);

    for r in 0..(ROUNDS_F / 2) {
        <Field as Poseidon2>::constant_layer(&mut current_state, r);
        <Field as Poseidon2>::sbox_layer(&mut current_state);
        Field::matmul_external(&mut current_state);
        outputs.push(current_state);
    }

    outputs
}

fn generate_partial_round_state<Field: RichField>(
    last_rount_output: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    assert_eq!(STATE_SIZE, WIDTH);
    let mut current_state = *last_rount_output;

    for r in 0..ROUNDS_P {
        current_state[0] += Field::from_canonical_u64(<Field as Poseidon2>::RC12_MID[r]);
        current_state[0] = <Field as Poseidon2>::sbox_monomial(current_state[0]);
        Field::matmul_internal(&mut current_state, &<Field as Poseidon2>::MAT_DIAG12_M_1);
        outputs.push(current_state);
    }

    outputs
}

fn generate_2st_full_round_state<Field: RichField>(
    last_rount_output: &[Field; STATE_SIZE],
) -> Vec<[Field; STATE_SIZE]> {
    let mut outputs = Vec::new();
    assert_eq!(STATE_SIZE, WIDTH);
    let mut current_state = *last_rount_output;

    for r in (ROUNDS_F / 2)..ROUNDS_F {
        <Field as Poseidon2>::constant_layer(&mut current_state, r);
        <Field as Poseidon2>::sbox_layer(&mut current_state);
        Field::matmul_external(&mut current_state);
        outputs.push(current_state);
    }

    outputs
}

/// Generate the outputs for a given preimage
fn generate_outputs<Field: RichField>(preimage: &[Field; STATE_SIZE]) -> [Field; STATE_SIZE] {
    assert_eq!(STATE_SIZE, WIDTH);
    <Field as Poseidon2>::poseidon2(*preimage)
}

/// Function to generate the Poseidon2 trace
pub fn generate_poseidon2_trace<F: RichField>(step_rows: &Vec<Row<F>>) -> [Vec<F>; NUM_COLS] {
    let mut trace_len = step_rows.len();
    if trace_len == 0 {
        trace_len = 1;
    }
    let mut trace: Vec<Vec<F>> = vec![vec![F::ZERO; trace_len]; NUM_COLS];

    let mut add_rows = |step_rows: &Vec<Row<F>>, is_exe: bool| {
        for (i, row) in step_rows.iter().enumerate() {
            if is_exe {
                trace[COL_IS_EXE][i] = F::ONE;
            } else {
                trace[COL_IS_EXE][i] = F::ZERO;
            }

            for j in 0..STATE_SIZE {
                trace[COL_INPUT_START + j][i] = row.preimage[j];
            }
            let outputs = generate_outputs(&row.preimage);
            for j in 0..STATE_SIZE {
                trace[COL_OUTPUT_START + j][i] = outputs[j];
            }

            let first_full_round_state = generate_1st_full_round_state(&row.preimage);
            let partial_round_state = generate_partial_round_state(
                first_full_round_state.last().unwrap().try_into().unwrap(),
            );
            let second_full_round_state = generate_2st_full_round_state(
                partial_round_state.last().unwrap().try_into().unwrap(),
            );
            for j in 0..(ROUNDS_F / 2) {
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
                trace[COL_PARTIAL_ROUND_END_STATE_START + j][i] =
                    partial_round_state[ROUNDS_P - 1][j];
            }
        }
    };

    add_rows(step_rows, true);

    if step_rows.len() == 0 {
        let preimage = (0..STATE_SIZE).map(|_| F::rand()).collect::<Vec<_>>();
        let dummy_rows = vec![Row {
            preimage: preimage.try_into().expect("can't fail"),
        }];
        add_rows(&dummy_rows, false);
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
    use crate::plonky2::columns::{COL_OUTPUT_START, NUM_COLS, STATE_SIZE};
    use crate::plonky2::generation::{
        generate_1st_full_round_state, generate_2st_full_round_state, generate_outputs,
        generate_partial_round_state, Row,
    };
    use plonky2::field::types::Sample;
    use plonky2::hash::poseidon2::Poseidon2;
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

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

        let trace = super::generate_poseidon2_trace(&step_rows);
        for (i, step_row) in step_rows.iter().enumerate().take(num_rows) {
            let expected_hash = <F as Poseidon2>::poseidon2(step_row.preimage);
            for j in 0..STATE_SIZE {
                assert_eq!(
                    expected_hash[j],
                    trace[COL_OUTPUT_START + j][i],
                    "Mismatch at row {i}, position {j}"
                );
            }
        }
    }
    #[test]
    fn generate_poseidon2_trace_with_dummy() {
        let step_rows = vec![];
        let trace: [Vec<F>; NUM_COLS] = super::generate_poseidon2_trace(&step_rows);
        assert_eq!(trace[0].len(), 1);
    }
}
