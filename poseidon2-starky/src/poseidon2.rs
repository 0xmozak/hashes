use ark_ff::PrimeField;
use lazy_static::lazy_static;
use std::sync::Arc;
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::fields::utils::from_hex;

// This file is lifted from the following library:
// https://github.com/HorizenLabs/poseidon2
//
// Local modifications:
// Some originally private functions have been made public in order to generate STARK traces.

#[derive(Clone, Debug)]
pub struct Poseidon2Params<F: PrimeField> {
    pub(crate) t: usize, // statesize
    pub(crate) d: usize, // sbox degree
    pub(crate) rounds_f_beginning: usize,
    pub(crate) rounds_p: usize,
    #[allow(dead_code)]
    pub(crate) rounds_f_end: usize,
    pub(crate) rounds: usize,
    pub(crate) mat_internal_diag_m_1: Vec<F>,
    pub(crate) _mat_internal: Vec<Vec<F>>,
    pub(crate) round_constants: Vec<Vec<F>>,
}

impl<F: PrimeField> Poseidon2Params<F> {
    #[allow(clippy::too_many_arguments)]

    pub const INIT_SHAKE: &'static str = "Poseidon2";

    pub fn new(
        t: usize,
        d: usize,
        rounds_f: usize,
        rounds_p: usize,
        mat_internal_diag_m_1: &[F],
        mat_internal: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert!(d == 3 || d == 5 || d == 7 || d == 11);
        assert_eq!(rounds_f % 2, 0);
        let r = rounds_f / 2;
        let rounds = rounds_f + rounds_p;

        Poseidon2Params {
            t,
            d,
            rounds_f_beginning: r,
            rounds_p,
            rounds_f_end: r,
            rounds,
            mat_internal_diag_m_1: mat_internal_diag_m_1.to_owned(),
            _mat_internal: mat_internal.to_owned(),
            round_constants: round_constants.to_owned(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Poseidon2<F: PrimeField> {
    pub(crate) params: Arc<Poseidon2Params<F>>,
}

impl<F: PrimeField> Poseidon2<F> {
    pub fn new(params: &Arc<Poseidon2Params<F>>) -> Self {
        Poseidon2 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();

        // Linear layer at beginning
        self.matmul_external(&mut current_state);

        for r in 0..self.params.rounds_f_beginning {
            let r = 0;
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            self.matmul_external(&mut current_state);
        }

        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        for r in self.params.rounds_f_beginning..p_end {
            current_state[0].add_assign(&self.params.round_constants[r][0]);
            current_state[0] = self.sbox_p(&current_state[0]);
            self.matmul_internal(&mut current_state, &self.params.mat_internal_diag_m_1);
        }

        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            self.matmul_external(&mut current_state);
        }
        current_state
    }

    fn sbox(&self, input: &[F]) -> Vec<F> {
        input.iter().map(|el| self.sbox_p(el)).collect()
    }

    fn sbox_p(&self, input: &F) -> F {
        let mut input2 = *input;
        input2.square_in_place();

        match self.params.d {
            3 => {
                let mut out = input2;
                out.mul_assign(input);
                out
            }
            5 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(input);
                out
            }
            7 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(&input2);
                out.mul_assign(input);
                out
            }
            _ => {
                panic!()
            }
        }
    }

    fn matmul_m4(&self, input: &mut [F]) {
        let t = self.params.t;
        let t4 = t / 4;
        for i in 0..t4 {
            let start_index = i * 4;
            let mut t_0 = input[start_index];
            t_0.add_assign(&input[start_index + 1]);
            let mut t_1 = input[start_index + 2];
            t_1.add_assign(&input[start_index + 3]);
            let mut t_2 = input[start_index + 1];
            t_2.double_in_place();
            t_2.add_assign(&t_1);
            let mut t_3 = input[start_index + 3];
            t_3.double_in_place();
            t_3.add_assign(&t_0);
            let mut t_4 = t_1;
            t_4.double_in_place();
            t_4.double_in_place();
            t_4.add_assign(&t_3);
            let mut t_5 = t_0;
            t_5.double_in_place();
            t_5.double_in_place();
            t_5.add_assign(&t_2);
            let mut t_6 = t_3;
            t_6.add_assign(&t_5);
            let mut t_7 = t_2;
            t_7.add_assign(&t_4);
            input[start_index] = t_6;
            input[start_index + 1] = t_5;
            input[start_index + 2] = t_7;
            input[start_index + 3] = t_4;
        }
    }

    fn matmul_external(&self, input: &mut [F]) {
        let t = self.params.t;
        match t {
            2 => {
                // Matrix circ(2, 1)
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
            }
            3 => {
                // Matrix circ(2, 1, 1)
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].add_assign(&sum);
            }
            4 => {
                // Applying cheap 4x4 MDS matrix to each 4-element part of the state
                self.matmul_m4(input);
            }
            8 | 12 | 16 | 20 | 24 => {
                // Applying cheap 4x4 MDS matrix to each 4-element part of the state
                self.matmul_m4(input);

                // Applying second cheap matrix for t > 4
                let t4 = t / 4;
                let mut stored = [F::zero(); 4];
                for l in 0..4 {
                    stored[l] = input[l];
                    for j in 1..t4 {
                        stored[l].add_assign(&input[4 * j + l]);
                    }
                }
                for i in 0..input.len() {
                    input[i].add_assign(&stored[i % 4]);
                }
            }
            _ => {
                panic!()
            }
        }
    }

    fn matmul_internal(&self, input: &mut [F], mat_internal_diag_m_1: &[F]) {
        let t = self.params.t;

        match t {
            2 => {
                // [2, 1]
                // [1, 3]
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                input[0].add_assign(&sum);
                input[1].double_in_place();
                input[1].add_assign(&sum);
            }
            3 => {
                // [2, 1, 1]
                // [1, 2, 1]
                // [1, 1, 3]
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].double_in_place();
                input[2].add_assign(&sum);
            }
            4 | 8 | 12 | 16 | 20 | 24 => {
                // Compute input sum
                let mut sum = input[0];
                input
                    .iter()
                    .skip(1)
                    .take(t - 1)
                    .for_each(|el| sum.add_assign(el));
                // Add sum + diag entry * element to each element
                for i in 0..input.len() {
                    input[i].mul_assign(&mat_internal_diag_m_1[i]);
                    input[i].add_assign(&sum);
                }
            }
            _ => {
                panic!()
            }
        }
    }

    fn add_rc(&self, input: &[F], rc: &[F]) -> Vec<F> {
        input
            .iter()
            .zip(rc.iter())
            .map(|(a, b)| {
                let mut r = *a;
                r.add_assign(b);
                r
            })
            .collect()
    }
}

type Scalar = FpGoldiLocks;

lazy_static! {
    pub static ref MAT_DIAG8_M_1: Vec<Scalar> = vec![
        from_hex("0xa98811a1fed4e3a5"),
        from_hex("0x1cc48b54f377e2a0"),
        from_hex("0xe40cd4f6c5609a26"),
        from_hex("0x11de79ebca97a4a3"),
        from_hex("0x9177c73d8b7e929c"),
        from_hex("0x2a6fe8085797e791"),
        from_hex("0x3de6e93329f8d5ad"),
        from_hex("0x3f7af9125da962fe"),
    ];
    pub static ref MAT_INTERNAL8: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0xa98811a1fed4e3a6"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x1cc48b54f377e2a1"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0xe40cd4f6c5609a27"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x11de79ebca97a4a4"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x9177c73d8b7e929d"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x2a6fe8085797e792"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x3de6e93329f8d5ae"),
            from_hex("0x0000000000000001"),
        ],
        vec![
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x0000000000000001"),
            from_hex("0x3f7af9125da962ff"),
        ],
    ];
    pub static ref RC8: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0xdd5743e7f2a5a5d9"),
            from_hex("0xcb3a864e58ada44b"),
            from_hex("0xffa2449ed32f8cdc"),
            from_hex("0x42025f65d6bd13ee"),
            from_hex("0x7889175e25506323"),
            from_hex("0x34b98bb03d24b737"),
            from_hex("0xbdcc535ecc4faa2a"),
            from_hex("0x5b20ad869fc0d033"),
        ],
        vec![
            from_hex("0xf1dda5b9259dfcb4"),
            from_hex("0x27515210be112d59"),
            from_hex("0x4227d1718c766c3f"),
            from_hex("0x26d333161a5bd794"),
            from_hex("0x49b938957bf4b026"),
            from_hex("0x4a56b5938b213669"),
            from_hex("0x1120426b48c8353d"),
            from_hex("0x6b323c3f10a56cad"),
        ],
        vec![
            from_hex("0xce57d6245ddca6b2"),
            from_hex("0xb1fc8d402bba1eb1"),
            from_hex("0xb5c5096ca959bd04"),
            from_hex("0x6db55cd306d31f7f"),
            from_hex("0xc49d293a81cb9641"),
            from_hex("0x1ce55a4fe979719f"),
            from_hex("0xa92e60a9d178a4d1"),
            from_hex("0x002cc64973bcfd8c"),
        ],
        vec![
            from_hex("0xcea721cce82fb11b"),
            from_hex("0xe5b55eb8098ece81"),
            from_hex("0x4e30525c6f1ddd66"),
            from_hex("0x43c6702827070987"),
            from_hex("0xaca68430a7b5762a"),
            from_hex("0x3674238634df9c93"),
            from_hex("0x88cee1c825e33433"),
            from_hex("0xde99ae8d74b57176"),
        ],
        vec![
            from_hex("0x488897d85ff51f56"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x1140737ccb162218"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xa7eeb9215866ed35"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x9bd2976fee49fcc9"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xc0c8f0de580a3fcc"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x4fb2dae6ee8fc793"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x343a89f35f37395b"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x223b525a77ca72c8"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x56ccb62574aaa918"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xc4d507d8027af9ed"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xa080673cf0b7e95c"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xf0184884eb70dcf8"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x044f10b0cb3d5c69"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xe9e3f7993938f186"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x1b761c80e772f459"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x606cec607a1b5fac"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x14a0c2e1d45f03cd"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x4eace8855398574f"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xf905ca7103eff3e6"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xf8c8f8d20862c059"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xb524fe8bdd678e5a"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0xfbb7865901a1ec41"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
            from_hex("0x0000000000000000"),
        ],
        vec![
            from_hex("0x014ef1197d341346"),
            from_hex("0x9725e20825d07394"),
            from_hex("0xfdb25aef2c5bae3b"),
            from_hex("0xbe5402dc598c971e"),
            from_hex("0x93a5711f04cdca3d"),
            from_hex("0xc45a9a5b2f8fb97b"),
            from_hex("0xfe8946a924933545"),
            from_hex("0x2af997a27369091c"),
        ],
        vec![
            from_hex("0xaa62c88e0b294011"),
            from_hex("0x058eb9d810ce9f74"),
            from_hex("0xb3cb23eced349ae4"),
            from_hex("0xa3648177a77b4a84"),
            from_hex("0x43153d905992d95d"),
            from_hex("0xf4e2a97cda44aa4b"),
            from_hex("0x5baa2702b908682f"),
            from_hex("0x082923bdf4f750d1"),
        ],
        vec![
            from_hex("0x98ae09a325893803"),
            from_hex("0xf8a6475077968838"),
            from_hex("0xceb0735bf00b2c5f"),
            from_hex("0x0a1a5d953888e072"),
            from_hex("0x2fcb190489f94475"),
            from_hex("0xb5be06270dec69fc"),
            from_hex("0x739cb934b09acf8b"),
            from_hex("0x537750b75ec7f25b"),
        ],
        vec![
            from_hex("0xe9dd318bae1f3961"),
            from_hex("0xf7462137299efe1a"),
            from_hex("0xb1f6b8eee9adb940"),
            from_hex("0xbdebcc8a809dfe6b"),
            from_hex("0x40fc1f791b178113"),
            from_hex("0x3ac1c3362d014864"),
            from_hex("0x9a016184bdb8aeba"),
            from_hex("0x95f2394459fbc25e"),
        ],
    ];
    pub static ref POSEIDON2_GOLDILOCKS_8_PARAMS: Arc<Poseidon2Params<Scalar>> = Arc::new(
        Poseidon2Params::new(8, 7, 8, 22, &MAT_DIAG8_M_1, &MAT_INTERNAL8, &RC8)
    );
}
