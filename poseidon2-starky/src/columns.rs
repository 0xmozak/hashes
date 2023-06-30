// The size of the state
pub(crate) const STATE_SIZE: usize = 8;

// Poseidon constants
pub(crate) const ROUNDS_F: usize = 4;
pub(crate) const ROUNDS_P: usize = 22;
pub(crate) const ROUNDS: usize = ROUNDS_F * 2 + ROUNDS_P; // 30

// The starting point of the column input
pub(crate) const COL_INPUT_START: usize = 0;

// The starting point of the column output
pub(crate) const COL_OUTPUT_START: usize = COL_INPUT_START + STATE_SIZE;

// The total number of columns
pub(crate) const NUM_COLS: usize = COL_OUTPUT_START + STATE_SIZE;
