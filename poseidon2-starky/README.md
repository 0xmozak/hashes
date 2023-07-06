# Poseidon2-Starky
This repository contains STARK circuits of a hash function called Poseidon2. The implementation is based on Starky, a powerful STARK library. You can find the Poseidon2 hash function repository [here](https://github.com/HorizenLabs/poseidon2) and the Starky library [here](https://github.com/mir-protocol/plonky2/tree/main/starky).

## Benchmark Results
Below are some benchmarking results for Poseidon2-Starky, including table generation and proving. These tests were conducted on a Macbook M1 Pro, with a STARK table of 102 columns and constraints at a degree of 7.

| Row    | Time |
| ------ | ---- |
| 2^16   | 7s   |
| 2^17   | 14s  |
| 2^18   | 28s  |
| 2^19   | 58s  |

## Future Improvements
1. Speed Enhancements: Future versions could potentially speed up the process. For example, we can reduce the constraints degree with more STARK table columns.

2. Support for more STATE_SIZE: The current version supports a fixed STATE_SIZE of 8. It should be relatively straightforward to adjust the code if your use case requires a larger STATE_SIZE.

## Contributing
Feel free to fork this repository and submit your pull requests for review. Any contribution that helps improve the performance of Poseidon2-Starky is welcome.