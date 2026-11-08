# On Fourier and Hadamard matrices of prime order

**Teo B and Chat GPT**, November 2026

How much of the Fourier matrix $F_p$ must one see to recognize it inside a complex Hadamard matrix of prime order $p$? For regular and Butson matrices, a single block of $F_p$ suffices at every odd prime. Throughout, matrices are considered equivalent up to permutations of rows and columns and multiplication of rows and columns by scalars of modulus one.

## Main results

- **Entrywise powers.** Raising each entry of a regular or Butson matrix of order $p$ to the $p$th power gives either a Hadamard matrix or a matrix of rank one. Rank one occurs exactly when the original matrix is equivalent to a Butson matrix of level $p$.
- **Fourier rows.** For every odd prime $p$, a regular or Butson matrix of order $p$ is equivalent to $F_p$ if it contains $(p-1)/2$ distinct Fourier rows in a common column ordering.
- **Two-orbit symmetry.** For every odd prime $p\leq 109$, a Butson matrix of order and level $p$ is equivalent to $F_p$ if it admits a simultaneous row and column symmetry with one fixed point and two cycles of length $(p-1)/2$ on each side.
- **The missing border.** Take the rows and columns of $F_p$ indexed by the nonzero squares modulo $p$. For every odd prime $p$, any regular or Butson matrix of order $p$ containing a submatrix equivalent to this block is itself equivalent to $F_p$.

## Contents

- [003.pdf](003.pdf) — the paper (42 pages).
- [003.tex](003.tex) — LaTeX source.
- [code003/](code003/) — Python and C++ programs, inputs, and certificates for the finite calculations.
- [code003/README.md](code003/README.md) — instructions for reproducing the checks, organized by section of the paper.

## Running the code

The checks use Python 3.10 or later, NumPy, SymPy, and a C++17 compiler. From the repository root, run

```sh
cd code003
python3 -m venv ../code003_env
. ../code003_env/bin/activate
python -m pip install -r requirements.txt
python run.py
```

To run just the order-11 enumeration and verify its certificates, use `python run.py order11`. The other groups are `w19`, `minors`, `orbits`, and `corners`. Each group prints `PASS` when its checks finish successfully.

