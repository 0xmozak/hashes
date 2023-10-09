/// The size of the state
pub const STATE_SIZE: usize = 12;
pub(crate) const SBOX_DEGREE: usize = 7;

/// Poseidon2 constants
pub(crate) const ROUNDS_F: usize = 8;
pub(crate) const ROUNDS_P: usize = 22;

/// If the row is execution or padding.
pub(crate) const COL_IS_EXE: usize = 0;

/// The starting point of the column input
pub(crate) const COL_INPUT_START: usize = 1;

/// The starting point of the state after each 1st full round
pub(crate) const COL_1ST_FULLROUND_STATE_START: usize = COL_INPUT_START + STATE_SIZE; // 1 + 12

/// The value of state[0] after each partial round
pub(crate) const COL_PARTIAL_ROUND_STATE_START: usize =
    COL_1ST_FULLROUND_STATE_START + STATE_SIZE * (ROUNDS_F / 2); // 1 + 12 + 48

/// The starting point of the state after the partial round
pub(crate) const COL_PARTIAL_ROUND_END_STATE_START: usize =
    COL_PARTIAL_ROUND_STATE_START + ROUNDS_P - 1; // 1 + 12 + 48 + 22

/// The starting point of the state after each 2nd full round
pub(crate) const COL_2ND_FULLROUND_STATE_START: usize =
    COL_PARTIAL_ROUND_END_STATE_START + STATE_SIZE; // 1 + 12 + 48 + 22 + 12

/// The starting point of the column output
/// This is the same as the last state after the 2nd full round
pub(crate) const COL_OUTPUT_START: usize =
    COL_2ND_FULLROUND_STATE_START + STATE_SIZE * ((ROUNDS_F / 2) - 1); // 1 + 12 + 48 + 22 + 12 + 36

/// The total number of columns
pub(crate) const NUM_COLS: usize = COL_2ND_FULLROUND_STATE_START + STATE_SIZE * (ROUNDS_F / 2); // 1 + 12 + 96 + 22 + 12 + 36 + 12
