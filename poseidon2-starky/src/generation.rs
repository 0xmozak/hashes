use crate::columns::{NUM_COLS, STATE_SIZE};
use plonky2::hash::hash_types::RichField;
use std::convert::TryInto;

#[derive(Debug, Clone, Default)]
// A row of the preimage
pub struct Row<Field: RichField> {
    preimage: [Field; STATE_SIZE],
}

// Function to generate Poseidon2 trace
pub fn generate_poseidon2_trace<Field: RichField>(
    step_rows: Vec<Row<Field>>,
) -> Result<[Vec<Field>; NUM_COLS], String> {
    let trace_len = step_rows.len();
    let mut trace: Vec<Vec<Field>> = vec![vec![Field::ZERO; trace_len]; NUM_COLS];

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
    use crate::columns::STATE_SIZE;
    use crate::generation::Row;
    use plonky2::field::types::{Field, Sample};
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};

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

        let result = super::generate_poseidon2_trace(step_rows);
        assert!(result.is_ok());
    }
}
