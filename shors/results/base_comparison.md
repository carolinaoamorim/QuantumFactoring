# Base comparison — order-finding on AerSimulator

Order `r` of each base `a` mod 15, recovered from the 8-qubit counting register
(`build_order_finding_circuit(a, n_count=8)`, 1024 shots). For an order-`r` base
the phase peaks land at multiples of `256 / r`.

| base `a` | true order `r` | recovered `r` | clean phase peaks (measured / 256) |
|----------|----------------|---------------|------------------------------------|
| 2  | 4 | 4 | 0, 64, 128, 192 |
| 4  | 2 | 2 | 0, 128 |
| 7  | 4 | 4 | 0, 64, 128, 192 |
| 8  | 4 | 4 | 0, 64, 128, 192 |
| 11 | 2 | 2 | 0, 128 |
| 13 | 4 | 4 | 0, 64, 128, 192 |

All recovered orders match the true multiplicative orders. Every base with even
order and `a^(r/2) ≢ −1 mod 15` yields the factors `(3, 5)`; `shor()` retries the
others automatically.

Across `random.seed(0..4)`, `shor(AerSimulator())` returns `(3, 5)` every run.
