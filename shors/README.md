# Shor's Algorithm — N = 15 and N = 21

An implementation of Shor's factoring algorithm for N = 15 and N = 21, built on a
Qiskit Runtime harness that runs **both** on the local `AerSimulator` and on
**real IBM Quantum hardware** through the same code path.

The inverse QFT, the controlled modular-multiplication gates, the order-finding
circuit, and the classical wrapper are all built from scratch.

For **N = 15** the modular multiplication collapses to a tidy swap network
(`c_amod15`). For **N = 21** (and any other N) the residues no longer permute so
neatly, so the controlled map `|y> -> |a^(2^k)·y mod N>` is built directly as a
permutation unitary (`c_amod`). Everything downstream — the order-finding circuit,
inverse QFT, and classical wrapper — is shared.

---

## File organization

| File | Purpose |
|------|---------|
| `shor_factoring.py` | Main program: runner harness + order-finding circuit + classical wrapper. |
| `shor_factoring.ipynb` | Runnable walkthrough: builds the circuit, draws it, plots the phase peaks, factors 15. |
| `results/base_comparison.md` | Per-base order-finding results (orders and phase peaks). |

---

## How Shor's algorithm is structured here

Shor's = **one quantum subroutine (order-finding) wrapped in classical code.**
Five of the six steps are classical; only order-finding is quantum.

**Classical wrapper** (`shor` / `find_order` in `shor_factoring.py`):
1. If N is even → return 2. (Handled up front.)
2. Pick a random base `a` coprime to N; if `gcd(a, N) ≠ 1` you've already found a factor.
3. **Quantum:** find the order `r` of `a` mod N (smallest `r` with `a^r ≡ 1 mod N`).
4. Recover `r` from the measured phase via continued fractions
   (`Fraction(...).limit_denominator(N)`).
5. If `r` is even and `a^(r/2) ≢ −1 mod N`, the factors are `gcd(a^(r/2) ± 1, N)`.
   Otherwise retry with a new `a`.

**Quantum order-finding circuit** (`build_order_finding_circuit`) — three blocks,
`n_count` counting qubits + `ceil(log2 N)` work qubits (12 qubits for N=15 with 8
counting):
- **Counting register** put into uniform superposition with Hadamards.
- **Controlled-Uₐ blocks**: `count[k]` controls multiplication by `a^(2^k) mod N`
  on the work register (`c_amod15` for N=15, `c_amod` otherwise).
- **Inverse QFT** (`inverse_qft`, hand-built) on the counting register, then measure.

This is standard quantum phase estimation: the counting register reads out a
phase ≈ s/r, from which the classical step extracts `r`. For a correct order-4
base (e.g. `a=2`, N=15), the phase peaks land cleanly at multiples of 256/4 = 64
(i.e. 0, 64, 128, 192).

### The N=15 modular-multiplication trick

For N = 15 the residues fit in 4 qubits and `a^4 ≡ 1`, so the controlled-Uₐ
collapses to a small swap network. Pinning the convention **work qubit 0 = LSB**
of the residue:

- The swap chain `swap(2,3); swap(1,2); swap(0,1)` cyclically shifts the bits and
  realizes `|y> → |2y mod 15>`; reversing the chain gives the `×8` inverse.
- Bases `7, 11, 13` need an extra `y → 15 − y` step, supplied by an X gate on all
  four work qubits.

(The work register starts at `|1>` and only ever visits residues in the
multiplicative orbit, so the unused `|0>↔|15>` swap that the X-all flip introduces
is harmless.)

### The general case (N=21 and beyond)

When the swap trick no longer applies, `c_amod(a, power, N)` builds the controlled
map straight from its permutation matrix: `|y> → |a^power · y mod N>` is a
permutation of the residues, so we reduce `a^power mod N` classically and emit a
single permutation unitary on `ceil(log2 N)` work qubits. Because the exponent is
reduced first, one gate covers even large `2^k` powers. `analyze_base` and the
`results/base_comparison.md` table use this path to compare bases for N=21
(e.g. `a=8` order 2 vs `a=2` order 6), both factoring `21 = 3 × 7`.

---

## Setup

```bash
pip install -r ../requirements.txt
```

### Run on the simulator

```bash
python3 shor_factoring.py
```

The `__main__` block factors **15** (hand-built gates) and **21** (general gates)
on `AerSimulator` and prints the recovered factors and peaks. Expected output:

```text
factors: (3, 5)
order r for a=7 mod 15: 4
top peaks: [('01000000', ...), ('00000000', ...), ('10000000', ...), ('11000000', ...)]

factors: (3, 7)
a=8 mod 21 -> order 2 factors (7, 3)
```

i.e. for N=15 clean spikes at measured values 0, 64, 128, 192, and N=21 factored
via the order-2 base `a=8`.

### Run on the IBM hardware

1. Make a free account at the IBM Quantum Platform and copy your API token.
2. Save credentials once (note: the old `ibm_quantum` channel was removed — use
   `ibm_quantum_platform`):
   ```python
   from qiskit_ibm_runtime import QiskitRuntimeService
   QiskitRuntimeService.save_account(
       channel="ibm_quantum_platform",
       token="YOUR_API_TOKEN",
       overwrite=True,
   )
   ```
3. Pick a backend and pass it instead of `AerSimulator()`:
   ```python
   from qiskit_ibm_runtime import QiskitRuntimeService
   service = QiskitRuntimeService()
   backend = service.least_busy(operational=True, simulator=False)
   print(shor(backend))
   ```

The harness (`run_circuit_and_get_counts`) is backend-agnostic, so nothing else
changes between simulator and hardware.

> **Expect noise on real hardware.** Shor's circuit for N=15 is deep; on current
> (NISQ) devices the phase peaks smear and factoring is unreliable. This is
> physics, not a code bug. The intended workflow is: **validate correctness on
> the simulator, then run on hardware and report the degradation** as a result —
> the sim-vs-hardware histogram comparison is the strongest part of a writeup.
> Also: real jobs sit in a queue (minutes to hours), so set up your account and
> test it early, not on the deadline.

---

## Notes

- Built and tested with Qiskit 2.x and `qiskit-ibm-runtime` 0.47.
- The execution harness is adapted from a shared base; the order-finding circuit,
  inverse QFT, and the controlled modular-multiplication gates are the
  from-scratch work.
