use crate::columns::{COL_INPUT_START, NUM_COLS, STATE_SIZE};
use num::bigint::BigUint;
use plonky2::hash::hash_types::RichField;
use std::convert::TryInto;
use zkhash::ark_ff::{BigInteger, PrimeField};
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::poseidon2::poseidon2::Poseidon2;
use zkhash::poseidon2::poseidon2_instance_goldilocks::POSEIDON2_GOLDILOCKS_8_PARAMS;

#[derive(Debug, Clone, Default)]
// A row of the preimage
pub struct Row<Field: RichField> {
    preimage: [Field; STATE_SIZE],
}

fn generate_outputs<Field: RichField>(preimage: &[Field; STATE_SIZE]) -> [Field; STATE_SIZE] {
    let mut outputs = [Field::ZERO; STATE_SIZE];

    let instance = Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS);
    assert_eq!(instance.get_t(), STATE_SIZE);
    let mut input = Vec::with_capacity(STATE_SIZE);
    for i in 0..STATE_SIZE {
        input.push(FpGoldiLocks::from(preimage[i].to_canonical_u64()));
    }
    let perm = instance.permutation(&input);
    for i in 0..STATE_SIZE {
        outputs[i] = Field::from_noncanonical_biguint(BigUint::from_bytes_le(
            &perm[i].into_bigint().to_bytes_le(),
        ));
    }

    outputs
}

// Function to generate Poseidon2 trace
pub fn generate_poseidon2_trace<Field: RichField>(
    step_rows: Vec<Row<Field>>,
) -> Result<[Vec<Field>; NUM_COLS], String> {
    let trace_len = step_rows.len();
    let mut trace: Vec<Vec<Field>> = vec![vec![Field::ZERO; trace_len]; NUM_COLS];

    for (i, row) in step_rows.iter().enumerate() {
        for j in 0..STATE_SIZE {
            trace[i][COL_INPUT_START + j] = row.preimage[j];
        }
        let outputs = generate_outputs(&row.preimage);
        for j in 0..STATE_SIZE {
            trace[i][COL_INPUT_START + STATE_SIZE + j] = outputs[j];
        }
    }

    trace.try_into().map_err(|trace_vector: Vec<Vec<Field>>| {
        format!(
            "Expected a Vec of length {} but it was {}",
            NUM_COLS,
            trace_vector.len()
        )
    })
}

#[cfg(test)]
mod test {
    use crate::columns::{COL_OUTPUT_START, STATE_SIZE};
    use crate::generation::Row;
    use plonky2::field::types::{Field, PrimeField64, Sample};
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    use zkhash::fields::goldilocks::FpGoldiLocks;
    use zkhash::poseidon2::poseidon2::Poseidon2;

    #[test]
    fn generate_poseidon2_trace() {
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;

        let num_rows = 16;
        let mut step_rows = Vec::with_capacity(num_rows);

        for _ in 0..num_rows {
            let mut preimage = [F::ZERO; STATE_SIZE];
            for i in 0..STATE_SIZE {
                preimage[i] = F::rand();
            }

            step_rows.push(Row { preimage });
        }

        let result = super::generate_poseidon2_trace(step_rows.clone()).unwrap();

        let instance = Poseidon2::new(&super::POSEIDON2_GOLDILOCKS_8_PARAMS);
        for i in 0..num_rows {
            let mut input = Vec::with_capacity(STATE_SIZE);
            for j in 0..STATE_SIZE {
                input.push(FpGoldiLocks::from(
                    step_rows[i].preimage[j].to_canonical_u64(),
                ));
            }
            let perm = instance.permutation(&input);
            for j in 0..STATE_SIZE {
                assert_eq!(
                    perm[j],
                    FpGoldiLocks::from(result[i][COL_OUTPUT_START + j].to_canonical_u64())
                );
            }
        }
    }
}
