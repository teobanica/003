# On Fourier and Hadamard matrices of prime order

Companion code for the paper by Teo B and Chat GPT. Section numbers below refer to the manuscript. The original order-11 inputs, programs and certificates are kept together in order11/.

## Running the checks

From the extracted code003 directory, with Python 3.10 or later and a C++17 compiler:

    python3 -m venv .venv
    . .venv/bin/activate
    python -m pip install -r requirements.txt
    python run.py

To run one part, use e.g. `python run.py order11` or `python run.py corners`. The names accepted are `w19`, `order11`, `minors`, `orbits` and `corners`. Each group prints PASS after all its checks have passed. The orbits group includes an exhaustive comparison of two generators at order 31 and takes the longest. The groups run serially. On a Mac the C++ compiler supplied by Command Line Tools is sufficient. The default compiler is `c++`; to select another, set `CXX`. Python assertions must be enabled.

The runner makes a fresh temporary directory for compiled programs and regenerated order-11 results. The archived data stay in place. To keep a transcript, redirect standard output to a file. Individual historical order-11 programs write outputs next to themselves when invoked directly.

| Folder | Section | What to check |
| --- | --- | --- |
| w19 | 2 | exponents.txt gives the displayed 19-by-19 matrix. completions.py enumerates the sign choices by exact polynomial pseudo-division over the Eisenstein integers, checks the four surviving representatives, and computes the minors 19, 98, 46 and 58 modulo 101. |
| order11/a_branch | 4 | The three-row incidence graphs, cycle colors, normalized labels and the Laurent coefficients of the affine example. brute_labels.cpp independently counts the labels. The six counts are 40, 368, 480, 3744, 2592 and 2880. |
| order11/a_square | 4 | Regenerate all 3441 dephased rows and the 2016 compatibility edges. The C++ clique search and the Python neighborhood check recover the unique Fourier clique. exclusion_certificate.json gives the missing neighbor edges. |
| minors | 6 | Integer counts of minor magnitudes, including 144 and 225 in the two blocks at order 19. A separate floating-point illustration compares polar factors and varies the free phase. |
| orbits | 7 | Two C++ programs enumerate and compare the complete word lists through 31. small.py checks complementary and mixed differences. odd.py verifies the character-polynomial calculations, even.py and parity.py verify the printed half-word and rank tables, and prime109.py checks the lift and the two squared-gap matrices at 109. |
| corners | 8–9 | deletions.cpp enumerates the normalized columns through 103. anchors.py reproduces all twelve selected rows at 503. strips.py reconstructs four test examples over exact cyclotomic fields at 5, 7 and 11. inverse.py illustrates the one-deletion argument at small primes. |

The symmetry range through 109 and the doubly deleted corner range through 103 are different tests. The rank and recursion figures are checked against the manuscript. The general theorems, including the border theorem of Section 10, have proofs in the paper. The strip examples use Fourier-equivalent inputs with permutations and phases applied outside the prescribed rectangle. They do not search for counterexamples of exponent rank five.
