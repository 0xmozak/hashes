use crate::columns::{COL_INPUT_START, COL_OUTPUT_START, NUM_COLS, STATE_SIZE};
use num::bigint::BigUint;
use plonky2::hash::hash_types::RichField;
use std::convert::TryInto;
use zkhash::ark_ff::{BigInteger, PrimeField};
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::poseidon2::poseidon2::Poseidon2;
use zkhash::poseidon2::poseidon2_instance_goldilocks::POSEIDON2_GOLDILOCKS_8_PARAMS;

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

// Generate the outputs for a given preimage
fn generate_outputs<Field: RichField>(preimage: &[Field; STATE_SIZE]) -> [Field; STATE_SIZE] {
    let mut outputs = [Field::ZERO; STATE_SIZE];
    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);

    let input = preimage
        .iter()
        .map(|i| FpGoldiLocks::from(i.to_canonical_u64()))
        .collect::<Vec<_>>();
    let perm = instance.permutation(&input);

    for i in 0..STATE_SIZE {
        outputs[i] = Field::from_noncanonical_biguint(BigUint::from_bytes_le(
            &perm[i].into_bigint().to_bytes_le(),
        ));
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
    use crate::generation::Row;
    use plonky2::field::types::{PrimeField64, Sample};
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    use zkhash::fields::goldilocks::FpGoldiLocks;
    use zkhash::poseidon2::poseidon2::Poseidon2;

    #[test]
    fn generate_poseidon2_trace() {
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;

        let num_rows = 12;
        let mut step_rows = Vec::with_capacity(num_rows);

        for _ in 0..num_rows {
            let preimage = (0..STATE_SIZE).map(|_| F::rand()).collect::<Vec<_>>();
            step_rows.push(Row {
                preimage: preimage.try_into().unwrap(),
            });
        }

        let trace = super::generate_poseidon2_trace(step_rows.clone());
        assert_eq!(trace.len(), 16);

        let instance = Poseidon2::new(&super::POSEIDON2_GOLDILOCKS_8_PARAMS);
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
