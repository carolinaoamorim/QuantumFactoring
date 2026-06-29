# Manual Quantum Fourier Transform and Noise Analysis

A from-scratch QFT / inverse-QFT built out of elementary gates, plus the tooling
to study how noise in the inverse QFT degrades period recovery in Shor's
algorithm.

**Central question:** how does noise in a manually implemented QFT affect the
Fourier-peak structure, period recovery, multiplicative-order recovery, and the
overall factoring success of Shor's algorithm?

---

## File organization

| File | Purpose |
|------|---------|
| `qft.py` | Library: manual `qft` / `inverse_qft`, approximate variants, validation against Qiskit, noise models, distribution metrics, period reconstruction, and the manual inverse QFT inside an N=15 Shor circuit. |
| `qft_constructions.ipynb` | Runnable walkthrough: builds and draws the circuits, validates them, runs the noise sweep, and plots the results. |

---

## How it is built

**Manual QFT** (`qft(n)`) — for each qubit from the top down: a Hadamard, then a
ladder of controlled-phase rotations `cp(pi / 2^(j-k))` from the lower qubits,
finished with a bit-reversal swap network.

**Manual inverse QFT** (`inverse_qft(n)`) — the exact reverse: bit-reversal swaps
first, then for each qubit the controlled-phase angles are *negated* and the
Hadamard is applied after the rotations. This is the piece that replaces the
built-in inverse QFT inside Shor's order-finding circuit.

**Approximate QFT** (`approximate_qft`, `approximate_inverse_qft`) — the same
construction with controlled rotations below an angle threshold dropped, trading
precision for a shallower circuit.

---

## What the analysis measures

- **Validation** (`validate_qft_against_builtin`, `round_trip_fidelity`): the
  manual QFT matches Qiskit's built-in QFT up to global phase, and QFT followed
  by inverse QFT returns the original state with fidelity ~ 1.
- **Known-period test** (`qft_distribution_for_periodic_state`): feeding a state
  with known period `r` into the inverse QFT produces peaks near multiples of
  `2^n / r`, isolating QFT behavior from modular-exponentiation errors.
- **Noise models** (`get_noise_model`): depolarizing (general gate error), phase
  damping (loss of phase coherence), and readout error.
- **Distribution metrics**: total-variation distance and Hellinger fidelity
  measure global histogram change; `peak_probability` and
  `estimated_period_success_from_distribution` measure whether the output is
  still *useful* for order recovery.
- **Period reconstruction** (`candidate_period_from_measurement`): continued
  fractions turn a measured `y / 2^n` into a candidate order.

A recurring finding: period recovery often survives well past the point where the
full distribution has visibly degraded.

---

## Setup

```bash
pip install -r ../requirements.txt
```

```bash
jupyter notebook qft_constructions.ipynb
```

The notebook imports the library from `qft.py`, so run it from inside this
directory.

---

## Notes

- Built and tested with Qiskit 2.x and `qiskit-aer`.
- The noise sweep organizes its results with Polars.
- The manual inverse QFT here is the same construction used (in-place) by the
  [Shor implementation](../shors/shor_factoring.py).
