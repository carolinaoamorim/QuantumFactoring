"""
Manual Quantum Fourier Transform and Noise Analysis.

A from-scratch QFT / inverse-QFT built out of elementary gates (Hadamards,
controlled-phase rotations, and bit-reversal swaps), plus the tooling used to
study how noise in the inverse QFT degrades period recovery in Shor's algorithm:

  * manual `qft` / `inverse_qft` (and approximate variants),
  * validation against Qiskit's built-in QFT (unitary equality up to global phase),
  * noise models (depolarizing, phase damping, readout, and a combined model
    stacking all three), with a decompose-first runner for the combined case,
  * distribution metrics (total-variation distance, Hellinger fidelity, peak mass),
  * period reconstruction via continued fractions,
  * the manual inverse QFT dropped into an N=15 Shor order-finding circuit.

This module holds the library functions; `qft_constructions.ipynb` imports them
and runs the experiments and plots.
"""

from fractions import Fraction
from math import pi

import numpy as np

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate
from qiskit.quantum_info import Operator, Statevector, state_fidelity
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel,
    ReadoutError,
    depolarizing_error,
    phase_damping_error,
)

SHOTS = 4096
RANDOM_SEED = 42

# STAGE 1 & 2: MANUAL QFT AND INVERSE QFT
def qft(n):
    """Manual QFT on n qubits: Hadamards + controlled-phase ladder + swaps."""
    circuit = QuantumCircuit(n)
    for j in reversed(range(n)):
        circuit.h(j)                               # Step 1: Hadamard on qubit j
        for k in reversed(range(j)):
            angle = pi / (2 ** (j - k))            # Step 2: controlled rotations
            circuit.cp(angle, k, j)
        circuit.barrier()                          # barriers are for visual only
    for i in range(n // 2):                        # Step 3: bit-reversal swaps
        circuit.swap(i, n - 1 - i)
    return circuit


def inverse_qft(n):
    """Manual inverse QFT on n qubits: the reverse of `qft`, with negated angles."""
    circuit = QuantumCircuit(n)
    for i in range(n // 2):                        # Step 1: bit-reversal swaps first
        circuit.swap(i, n - 1 - i)
        circuit.barrier()
    for j in range(n):                             # Step 2: undo phase kicks, then H
        for m in range(j):
            circuit.cp(-pi / (2 ** (j - m)), m, j)
        circuit.h(j)
        circuit.barrier()
    return circuit


# STAGE 9: APPROXIMATE QFT (drop controlled rotations below a threshold)
def approximate_qft(n, min_angle=pi / 16):
    """QFT with controlled rotations smaller than `min_angle` removed."""
    circuit = QuantumCircuit(n)
    for j in reversed(range(n)):
        circuit.h(j)
        for k in reversed(range(j)):
            angle = pi / (2 ** (j - k))
            if abs(angle) >= min_angle:
                circuit.cp(angle, k, j)
        circuit.barrier()
    for i in range(n // 2):
        circuit.swap(i, n - 1 - i)
    return circuit


def approximate_inverse_qft(n, min_angle=pi / 16):
    """Inverse QFT with controlled rotations smaller than `min_angle` removed."""
    circuit = QuantumCircuit(n)
    for i in range(n // 2):
        circuit.swap(i, n - 1 - i)
    circuit.barrier()
    for j in range(n):
        for m in range(j):
            angle = -pi / (2 ** (j - m))
            if abs(angle) >= min_angle:
                circuit.cp(angle, m, j)
        circuit.h(j)
        circuit.barrier()
    return circuit


# STAGE 3: IDEAL VALIDATION AGAINST QISKIT'S BUILT-IN QFT
def strip_barriers(circuit):
    """Return a copy of `circuit` without barriers (so Operator/Statevector
    comparisons work cleanly). Barriers are visual only."""
    clean = QuantumCircuit(circuit.num_qubits, circuit.num_clbits)
    for instruction, qargs, cargs in circuit.data:
        if instruction.name == "barrier":
            continue
        q_indices = [circuit.find_bit(q).index for q in qargs]
        c_indices = [circuit.find_bit(c).index for c in cargs]
        clean.append(instruction.copy(), q_indices, c_indices)
    return clean


def equal_up_to_global_phase(matrix_a, matrix_b, tolerance=1e-8):
    """Two unitaries that differ only by a global phase are physically equal."""
    if matrix_a.shape != matrix_b.shape:
        return False
    nonzero = np.argwhere(np.abs(matrix_b) > tolerance)
    if len(nonzero) == 0:
        return np.allclose(matrix_a, matrix_b, atol=tolerance)
    row, col = nonzero[0]
    phase = matrix_a[row, col] / matrix_b[row, col]
    if abs(phase) < tolerance:
        return False
    phase = phase / abs(phase)
    return np.allclose(matrix_a, phase * matrix_b, atol=tolerance)


def builtin_qft(n):
    """Qiskit's built-in QFT circuit, for comparison."""
    circuit = QuantumCircuit(n)
    circuit.append(QFTGate(n), range(n))
    return circuit


def validate_qft_against_builtin(max_qubits=5):
    """Print whether the manual QFT matches Qiskit's built-in QFT for 2..max_qubits."""
    print("Manual QFT versus Qiskit built-in QFT")
    print("-" * 55)
    for n in range(2, max_qubits + 1):
        manual = Operator(strip_barriers(qft(n))).data
        built = Operator(builtin_qft(n)).data
        equivalent = equal_up_to_global_phase(manual, built)
        print(f"{n} qubits:", "equivalent up to global phase =", equivalent)


def round_trip_fidelity(n, input_integer):
    """Apply QFT then inverse QFT to |input_integer>; fidelity should be ~1."""
    dimension = 2 ** n
    initial = Statevector.from_int(input_integer, dimension)
    final = initial.evolve(strip_barriers(qft(n)))
    final = final.evolve(strip_barriers(inverse_qft(n)))
    return state_fidelity(initial, final)


# STAGE 5: KNOWN-PERIOD STATE (isolate QFT behavior from mod-exp errors)
def periodic_statevector(n, period, offset=0):
    """Uniform superposition over |offset>, |offset+period>, ... — known period."""
    dimension = 2 ** n
    if period <= 0:
        raise ValueError("period must be positive.")
    if offset < 0 or offset >= dimension:
        raise ValueError("offset must be inside the register.")
    selected = list(range(offset, dimension, period))
    amplitudes = np.zeros(dimension, dtype=complex)
    for state in selected:
        amplitudes[state] = 1 / np.sqrt(len(selected))
    return Statevector(amplitudes)


def expected_peak_locations(n, period):
    """Ideal inverse-QFT peaks fall near multiples of 2^n / period."""
    dimension = 2 ** n
    spacing = dimension / period
    peaks = [round(k * spacing) % dimension for k in range(period)]
    return sorted(set(peaks))


def qft_distribution_for_periodic_state(n, period, offset=0):
    """Probabilities after applying the manual inverse QFT to a known-period state."""
    state = periodic_statevector(n, period, offset)
    output = state.evolve(strip_barriers(inverse_qft(n)))
    return np.abs(output.data) ** 2


# STAGE 6: NOISE MODELS
def create_depolarizing_noise_model(one_qubit_error, two_qubit_error):
    """Depolarizing noise for general gate imperfections."""
    noise_model = NoiseModel()
    one = depolarizing_error(one_qubit_error, 1)
    two = depolarizing_error(two_qubit_error, 2)
    noise_model.add_all_qubit_quantum_error(one, ["h"])
    noise_model.add_all_qubit_quantum_error(two, ["cp", "swap", "cx"])
    return noise_model


def create_phase_damping_noise_model(damping_probability):
    """Phase damping reduces coherence -- especially relevant to the phase-based QFT."""
    noise_model = NoiseModel()
    one = phase_damping_error(damping_probability)
    two = one.tensor(one)
    noise_model.add_all_qubit_quantum_error(one, ["h"])
    noise_model.add_all_qubit_quantum_error(two, ["cp", "swap", "cx"])
    return noise_model


def create_readout_noise_model(p_0_to_1, p_1_to_0):
    """Readout error: the measured classical bit flips with the given probabilities."""
    noise_model = NoiseModel()
    matrix = [[1 - p_0_to_1, p_0_to_1], [p_1_to_0, 1 - p_1_to_0]]
    noise_model.add_all_qubit_readout_error(ReadoutError(matrix))
    return noise_model


def create_combined_noise_model(
    one_qubit_error, two_qubit_error, damping_probability, p_0_to_1, p_1_to_0
):
    """Stack depolarizing + phase damping + readout into a single NoiseModel.

    The isolated builders above each model one channel; a real device suffers all
    three at once. On every gate we *compose* the depolarizing and phase-damping
    channels (depolarizing first, then phase damping) so both act on the same gate,
    and we add readout error on top for the measurement. This is genuinely new
    behavior -- you cannot get it by passing different arguments to
    `create_depolarizing_noise_model`, which only knows how to build depolarizing
    noise.
    """
    noise_model = NoiseModel()

    depolarizing_1 = depolarizing_error(one_qubit_error, 1)
    phase_1 = phase_damping_error(damping_probability)
    combined_1 = depolarizing_1.compose(phase_1)

    depolarizing_2 = depolarizing_error(two_qubit_error, 2)
    phase_2 = phase_1.tensor(phase_1)
    combined_2 = depolarizing_2.compose(phase_2)

    noise_model.add_all_qubit_quantum_error(combined_1, ["h"])
    noise_model.add_all_qubit_quantum_error(combined_2, ["cp", "swap", "cx"])

    matrix = [[1 - p_0_to_1, p_0_to_1], [p_1_to_0, 1 - p_1_to_0]]
    noise_model.add_all_qubit_readout_error(ReadoutError(matrix))
    return noise_model


def get_noise_model(noise_type, noise_level):
    """Dispatch wrapper so experiment sweeps stay clean. 'ideal' -> None."""
    if noise_type == "ideal":
        return None
    if noise_type == "depolarizing":
        return create_depolarizing_noise_model(noise_level, 2 * noise_level)
    if noise_type == "phase_damping":
        return create_phase_damping_noise_model(noise_level)
    if noise_type == "readout":
        return create_readout_noise_model(noise_level, noise_level)
    if noise_type == "combined":
        # A single sweep level drives every channel, so combined and isolated
        # runs sit on the same x-axis for a like-for-like comparison.
        return create_combined_noise_model(
            one_qubit_error=noise_level,
            two_qubit_error=2 * noise_level,
            damping_probability=noise_level,
            p_0_to_1=noise_level,
            p_1_to_0=noise_level,
        )
    raise ValueError(f"Unknown noise type: {noise_type}")


# DISTRIBUTION METRICS
def counts_to_probability_vector(counts, n):
    """Convert Qiskit counts into a length-2^n probability vector."""
    dimension = 2 ** n
    probabilities = np.zeros(dimension, dtype=float)
    total = sum(counts.values())
    for bitstring, count in counts.items():
        integer = int(bitstring.replace(" ", ""), 2)
        probabilities[integer] = count / total
    return probabilities


def total_variation_distance(p, q):
    """D_TV(P, Q) = 1/2 sum_y |P(y) - Q(y)|. Zero means identical."""
    return 0.5 * np.sum(np.abs(p - q))


def hellinger_fidelity(p, q):
    """Hellinger fidelity; near 1 indicates strong similarity."""
    return (np.sum(np.sqrt(p * q))) ** 2


def peak_probability(probability_vector, peaks, window=0):
    """Probability mass within `window` of the expected peaks.

    Uses a set of indices so overlapping windows are not double-counted.
    """
    dimension = len(probability_vector)
    included = set()
    for peak in peaks:
        for shift in range(-window, window + 1):
            included.add((peak + shift) % dimension)
    return float(sum(probability_vector[i] for i in included))


def top_measurement(probability_vector):
    """Index of the most likely measurement outcome."""
    return int(np.argmax(probability_vector))


def closest_expected_peak_distance(measured_peak, expected_peaks, dimension):
    """Circular distance from a measured peak to the nearest expected peak."""
    distances = []
    for peak in expected_peaks:
        direct = abs(measured_peak - peak)
        distances.append(min(direct, dimension - direct))
    return min(distances)


def summarize_peak_behavior(probability_vector, expected_peaks):
    """Simple Fourier-peak diagnostics for a measured distribution."""
    dimension = len(probability_vector)
    measured_peak = top_measurement(probability_vector)
    return {
        "measured_peak": measured_peak,
        "peak_distance": closest_expected_peak_distance(
            measured_peak, expected_peaks, dimension
        ),
        "expected_peak_mass_window_1": peak_probability(
            probability_vector, expected_peaks, window=1
        ),
        "max_probability": float(probability_vector[measured_peak]),
    }


# STAGE 7: PERIOD RECONSTRUCTION (continued fractions)
def candidate_period_from_measurement(measured_integer, n, max_denominator):
    """Approximate y / 2^n by a rational and return its denominator (candidate r)."""
    if measured_integer == 0:
        return None
    phase = Fraction(measured_integer, 2 ** n)
    return phase.limit_denominator(max_denominator).denominator


def measurement_recovers_period(measured_integer, n, true_period):
    """Whether a measurement yields the true period or a useful divisor of it."""
    candidate = candidate_period_from_measurement(measured_integer, n, 2 ** n)
    if candidate is None:
        return False
    if candidate == true_period:
        return True
    return candidate != 0 and true_period % candidate == 0


def estimated_period_success_from_distribution(probability_vector, n, true_period):
    """Total probability of measurements that would recover the known period."""
    success = 0.0
    for measured_integer, probability in enumerate(probability_vector):
        if measurement_recovers_period(measured_integer, n, true_period):
            success += probability
    return float(success)

# KNOWN-PERIOD INVERSE-QFT EXPERIMENT (ideal vs noisy)
def build_periodic_iqft_circuit(n, period, offset=0):
    """Prepare a known-period state, apply the manual inverse QFT, and measure."""
    return build_periodic_iqft_circuit_with_custom_iqft(n, period, inverse_qft, offset)


def build_periodic_iqft_circuit_with_custom_iqft(n, period, iqft_builder, offset=0):
    """Known-period circuit using whichever inverse-QFT builder is passed in."""
    state = periodic_statevector(n, period, offset)
    circuit = QuantumCircuit(n, n)
    circuit.initialize(state.data, range(n))
    circuit.barrier()
    circuit.compose(strip_barriers(iqft_builder(n)), qubits=range(n), inplace=True)
    circuit.barrier()
    circuit.measure(range(n), range(n))
    return circuit


# Basis gates for the "safe" runner. The four noise-targeted names (h, cp, swap,
# cx) are kept so the noise model still attaches; the rest let the opaque
# `initialize` instruction unroll into reset + rotations instead of staying a
# single black-box gate.
SAFE_BASIS_GATES = ["h", "cp", "swap", "cx", "rz", "ry", "rx", "x", "id", "reset"]


def run_measured_circuit(circuit, noise_model=None, shots=SHOTS, seed=RANDOM_SEED):
    """Transpile and run a measured circuit on Aer; return (counts, compiled).

    `seed` drives both transpilation and simulation; vary it across repetitions to
    get independent noise realizations (and therefore error bars).
    """
    backend = AerSimulator(noise_model=noise_model)
    compiled = transpile(
        circuit, backend, optimization_level=1, seed_transpiler=seed
    )
    result = backend.run(compiled, shots=shots, seed_simulator=seed).result()
    return result.get_counts(), compiled


def run_measured_circuit_safe(circuit, noise_model=None, shots=SHOTS, seed=RANDOM_SEED):
    """Like `run_measured_circuit`, but decomposes to basis gates first.

    Necessary, not cosmetic: the combined noise model crashes on the opaque
    `initialize` instruction (an "empty Kraus" / non-hermitian eigensystem error)
    because Aer tries to attach noise to a black-box state-prep gate. Transpiling
    to an explicit basis first (see `SAFE_BASIS_GATES`) unrolls `initialize` into
    reset + rotations while preserving the noise-targeted gate names, so the model
    both runs and still applies its errors.
    """
    backend = AerSimulator(noise_model=noise_model)
    compiled = transpile(
        circuit,
        backend,
        basis_gates=SAFE_BASIS_GATES,
        optimization_level=1,
        seed_transpiler=seed,
    )
    result = backend.run(compiled, shots=shots, seed_simulator=seed).result()
    return result.get_counts(), compiled


def run_periodic_qft_experiment(n, period, noise_type, noise_level, seed=RANDOM_SEED):
    """Run one known-period inverse-QFT experiment and return a metrics row.

    The combined noise model is routed through `run_measured_circuit_safe`, which
    it needs to run at all. `seed` lets a caller repeat the same configuration with
    independent noise draws for error bars.
    """
    runner = run_measured_circuit_safe if noise_type == "combined" else run_measured_circuit
    circuit = build_periodic_iqft_circuit(n, period)
    ideal_counts, _ = runner(circuit, noise_model=None, seed=seed)
    noise_model = get_noise_model(noise_type, noise_level)
    noisy_counts, noisy_compiled = runner(circuit, noise_model=noise_model, seed=seed)

    ideal = counts_to_probability_vector(ideal_counts, n)
    noisy = counts_to_probability_vector(noisy_counts, n)
    peaks = expected_peak_locations(n, period)

    return {
        "num_qubits": n,
        "period": period,
        "noise_type": noise_type,
        "noise_level": noise_level,
        "seed": seed,
        "shots": SHOTS,
        "logical_depth": circuit.depth(),
        "transpiled_depth": noisy_compiled.depth(),
        "total_gates": sum(noisy_compiled.count_ops().values()),
        "tv_distance": float(total_variation_distance(ideal, noisy)),
        "hellinger_fidelity": float(hellinger_fidelity(ideal, noisy)),
        "peak_probability_window_0": peak_probability(noisy, peaks, window=0),
        "peak_probability_window_1": peak_probability(noisy, peaks, window=1),
        "period_success_probability": estimated_period_success_from_distribution(
            noisy, n, period
        ),
        **summarize_peak_behavior(noisy, peaks),
    }


def run_exact_vs_approx_experiment(n, period, qft_label, iqft_builder, noise_type, noise_level):
    """Compare an exact vs approximate inverse QFT under the same noise."""
    circuit = build_periodic_iqft_circuit_with_custom_iqft(n, period, iqft_builder)
    ideal_circuit = build_periodic_iqft_circuit_with_custom_iqft(n, period, inverse_qft)

    ideal_counts, _ = run_measured_circuit(ideal_circuit, noise_model=None)
    noise_model = get_noise_model(noise_type, noise_level)
    noisy_counts, noisy_compiled = run_measured_circuit(circuit, noise_model=noise_model)

    ideal = counts_to_probability_vector(ideal_counts, n)
    noisy = counts_to_probability_vector(noisy_counts, n)
    peaks = expected_peak_locations(n, period)

    return {
        "qft_version": qft_label,
        "num_qubits": n,
        "period": period,
        "noise_type": noise_type,
        "noise_level": noise_level,
        "transpiled_depth": noisy_compiled.depth(),
        "total_gates": sum(noisy_compiled.count_ops().values()),
        "tv_distance": float(total_variation_distance(ideal, noisy)),
        "hellinger_fidelity": float(hellinger_fidelity(ideal, noisy)),
        "peak_probability_window_1": peak_probability(noisy, peaks, window=1),
        "period_success_probability": estimated_period_success_from_distribution(
            noisy, n, period
        ),
    }

# STAGE 8: MANUAL INVERSE QFT INSIDE SHOR (N = 15, a = 2, r = 4)
# Expected 8-qubit control peaks: 0, 64, 128, 192.
def multiply_by_2_mod_15():
    """|x> -> |2x mod 15> on a 4-qubit register (swap network)."""
    circuit = QuantumCircuit(4, name="M_2 mod 15")
    circuit.swap(2, 3)
    circuit.swap(1, 2)
    circuit.swap(0, 1)
    return circuit.to_gate()


def multiply_by_4_mod_15():
    """|x> -> |4x mod 15> on a 4-qubit register (swap network)."""
    circuit = QuantumCircuit(4, name="M_4 mod 15")
    circuit.swap(1, 3)
    circuit.swap(0, 2)
    return circuit.to_gate()


def build_shor_n15_manual_iqft():
    """N=15, a=2 order-finding circuit using the manual inverse QFT."""
    from qiskit import ClassicalRegister, QuantumRegister

    num_control, num_target = 8, 4
    control = QuantumRegister(num_control, name="control")
    target = QuantumRegister(num_target, name="target")
    classical = ClassicalRegister(num_control, name="measurement")
    circuit = QuantumCircuit(control, target, classical)

    circuit.x(target[0])                           # target := |1>
    circuit.barrier()
    circuit.h(control)                             # control into superposition
    circuit.barrier()

    circuit.append(multiply_by_2_mod_15().control(1), [control[0], *target])
    circuit.append(multiply_by_4_mod_15().control(1), [control[1], *target])
    circuit.barrier()

    circuit.compose(                               # manual inverse QFT on control
        strip_barriers(inverse_qft(num_control)), qubits=control, inplace=True
    )
    circuit.barrier()
    circuit.measure(control, classical)
    return circuit


def run_shor_manual_iqft(noise_type="ideal", noise_level=0.0):
    """Run the N=15 Shor circuit with the manual inverse QFT; return (counts, compiled)."""
    circuit = build_shor_n15_manual_iqft()
    noise_model = get_noise_model(noise_type, noise_level)
    backend = AerSimulator(noise_model=noise_model)
    compiled = transpile(
        circuit, backend, optimization_level=1, seed_transpiler=RANDOM_SEED
    )
    result = backend.run(compiled, shots=SHOTS, seed_simulator=RANDOM_SEED).result()
    return result.get_counts(), compiled
