# Quantum Factoring

Implementations and simulations of quantum factoring algorithms in **Qiskit** — Shor's algorithm, the Quantum Fourier Transform, and Regev's 2023 multidimensional algorithm — alongside a full mathematical writeup.

Factoring large integers is believed to be hard for classical computers, and that hardness is what secures RSA. Quantum computers change the picture: this repository works through *how*, building each algorithm from elementary gates and simulating it end to end.

---

## What's inside

| Project | What it does |
|---|---|---|
| [`qft/`](qft/) | The Quantum Fourier Transform built from scratch (Hadamards, controlled-phase rotations, swaps), for 2–5 qubits. |
| [`shor_n21/`](shor_n21/) | Shor's algorithm factoring **N = 21** via quantum order-finding, with a comparison of two bases. |
| [`regev/`](regev/) | Regev's 2023 algorithm — a multidimensional generalization of Shor's that uses fewer gates. Theory written up; implementation in progress. |
| [`quantum_utils/`](quantum_utils/) | Shared code (the QFT / inverse-QFT routines) imported by the projects above. |

A companion LaTeX writeup develops the mathematics behind each implementation. *(Link it here when public.)*

---

## Highlights

### Shor's algorithm factoring N = 21

The order-finding circuit encodes the period of `a^x mod 21` into the phase of a counting register, then reads it out with an inverse QFT:

![Order-finding circuit](shor_n21/figures/circuit_diagram.png)

Both `a = 8` and `a = 2` correctly recover `21 = 3 × 7`, but with very different efficiency — the number of measurement peaks equals the order `r`:

| Metric | a = 8 | a = 2 |
|---|---|---|
| Order *r* | 2 | 6 |
| Counting qubits | 3 | 5 |
| Transpiled depth | 1240 | 6199 |
| Distinct outcomes (peaks) | 2 | 32 |

![Measured phase distributions](shor_n21/figures/compare_bases.png)

### The Quantum Fourier Transform, from scratch

Built directly from Hadamard and controlled-phase gates, with the bit-reversal swap layer — shown here scaling from 2 to 5 qubits. The growing "triangle" of phase gates is the `O(n²)` gate count made visible:

![QFT circuits 2–5 qubits](qft/figures/qft_all.png)

---

## Quick start

```bash
git clone https://github.com/<your-username>/quantum-factoring.git
cd quantum-factoring
pip install -r requirements.txt
```

Then open any project's notebook, e.g.:

```bash
jupyter notebook shor_n21/shor_factoring.ipynb
```

Each subfolder has its own README with details on running that piece.

---

## Background

Shor's algorithm (1994) factors an integer `N` by reducing the problem to **order-finding**: for a base `a`, find the smallest `r` with `a^r ≡ 1 (mod N)`. A quantum circuit finds `r` using phase estimation and the inverse QFT, after which the factors follow classically from `gcd(a^(r/2) ± 1, N)`. The QFT is the engine that detects the hidden period. Regev's algorithm (2023) generalizes this to several bases at once, turning the single period into a **multidimensional lattice** and cutting the gate count from roughly `n²` to `n^1.5` — the first asymptotic improvement to Shor's circuit in nearly three decades.

---

## Built with

[Qiskit](https://www.ibm.com/quantum/qiskit) · [Qiskit Aer](https://github.com/Qiskit/qiskit-aer) (simulator) · NumPy · Matplotlib · Python 3

---

## Author

Carolina Amorim — built while studying quantum computing during a summer research program at UIUC (2026).

*A study and reproduction project; not original research. Algorithms are due to Peter Shor and Oded Regev as cited in the writeup.*