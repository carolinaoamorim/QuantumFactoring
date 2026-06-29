# Base comparison

Order `r` of each base `a`, recovered from the counting register on
`AerSimulator`. For an order-`r` base the phase peaks land at multiples of
`2^n_count / r`.

## N = 15 (hand-built swap gates, `n_count = 8`)

| base `a` | true order `r` | recovered `r` | clean phase peaks (measured / 256) |
|----------|----------------|---------------|------------------------------------|
| 2  | 4 | 4 | 0, 64, 128, 192 |
| 4  | 2 | 2 | 0, 128 |
| 7  | 4 | 4 | 0, 64, 128, 192 |
| 8  | 4 | 4 | 0, 64, 128, 192 |
| 11 | 2 | 2 | 0, 128 |
| 13 | 4 | 4 | 0, 64, 128, 192 |

All recovered orders match the true multiplicative orders. Across
`random.seed(0..4)`, `shor(AerSimulator(), N=15)` returns `(3, 5)` every run.

## N = 21 (general permutation-unitary gates)

Two coprime bases that both factor `21 = 3 × 7`, showing why `a = 8` is the more
efficient choice (`analyze_base(a, backend, N=21, ...)`):

| metric            | a = 8 (order 2) | a = 2 (order 6) |
|-------------------|-----------------|-----------------|
| order `r`         | 2               | 6               |
| factors found     | (7, 3)          | (7, 3)          |
| counting qubits   | 3               | 5               |
| total qubits      | 8               | 10              |
| circuit depth     | 1240            | 6205            |
| distinct outcomes | 2               | 32              |

- **a = 8 (order 2):** minimal — fewer qubits, shallower circuit, two clean peaks.
- **a = 2 (order 6):** illustrative — the full six-peak period spectrum.

Both recover nontrivial factors of 21; the small-order base reaches the answer
with a fraction of the circuit depth.
