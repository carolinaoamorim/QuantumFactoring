from fractions import Fraction
from math import gcd, pi
import random

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2 as Sampler

"""
Shor's algorithm for N = 15.

One quantum subroutine (order-finding via quantum phase estimation) wrapped in a
classical control loop. Built on a backend-agnostic Qiskit Runtime harness so the
exact same code path runs on the local AerSimulator and on real IBM hardware.

  * develop + validate on AerSimulator (noise-free -> clean factoring)
  * final run on real IBM hardware via the SAME code path (expect noisy peaks)

The inverse QFT, the controlled modular-multiplication gates, the order-finding
circuit, and the classical wrapper are all built from scratch.

----------------------------------------------------------------------
SETUP for real hardware (do this early -- the queue can be slow):

    from qiskit_ibm_runtime import QiskitRuntimeService
    QiskitRuntimeService.save_account(
        channel="ibm_quantum_platform",   # NOTE: the old "ibm_quantum" channel was removed
        token="YOUR_API_TOKEN",
        overwrite=True,
    )
Then get a backend with:
    service = QiskitRuntimeService()
    backend = service.least_busy(operational=True, simulator=False)
----------------------------------------------------------------------
"""

# ======================================================================
# HARNESS  --  backend-agnostic. Pass AerSimulator() for dev, or an IBM
# backend for the real run; nothing else changes.
# ======================================================================
def run_circuit_and_get_counts(circuit, backend, shots=1024):
    """Run a circuit on any backend (sim or real HW) and return counts."""
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)
    sampler = Sampler(mode=backend)
    job = sampler.run([isa_circuit], shots=shots)
    result = job.result()
    return result[0].data.c.get_counts()   # 'c' = our ClassicalRegister name


# ======================================================================
# 1. INVERSE QFT  --  hand-built. The core "from scratch" piece.
# ======================================================================
def inverse_qft(circuit, qubits):
    """In-place inverse QFT on the list `qubits`."""
    n = len(qubits)
    for i in range(n // 2):                       # bit-reversal swap network
        circuit.swap(qubits[i], qubits[n - 1 - i])
    for j in range(n):                            # reversed rotation ladder
        for m in range(j):
            circuit.cp(-pi / (2 ** (j - m)), qubits[m], qubits[j])
        circuit.h(qubits[j])
    return circuit


# ======================================================================
# 2. CONTROLLED MODULAR MULTIPLICATION  a^power mod 15
#    For N=15 the controlled-U reduces to a swap network plus X gates.
#    Work-qubit convention: qubit 0 is the LSB of the 4-bit residue.
#    The swaps alone permute |y> -> |2y mod 15> (or the |8y> inverse);
#    the X-all flip on {7,11,13} supplies the extra y -> 15-y those need.
# ======================================================================
def c_amod15(a, power):
    if a not in (2, 4, 7, 8, 11, 13):
        raise ValueError("a must be one of 2,4,7,8,11,13 for N=15")
    U = QuantumCircuit(4)
    for _ in range(power):                         # one pass = one multiply by a
        if a in (2, 13):
            U.swap(2, 3); U.swap(1, 2); U.swap(0, 1)
        if a in (7, 8):
            U.swap(0, 1); U.swap(1, 2); U.swap(2, 3)
        if a in (4, 11):
            U.swap(1, 3); U.swap(0, 2)
        if a in (7, 11, 13):                        # the bit-flips the swaps alone miss
            for q in range(4):
                U.x(q)
    gate = U.to_gate()
    gate.name = f"{a}^{power} mod 15"
    return gate.control()


# ======================================================================
# 3. ORDER-FINDING CIRCUIT (phase estimation).
#    8 counting + 4 work = 12 qubits.
# ======================================================================
def build_order_finding_circuit(a, n_count=8):
    count = QuantumRegister(n_count, "count")
    work = QuantumRegister(4, "work")
    creg = ClassicalRegister(n_count, "c")
    qc = QuantumCircuit(count, work, creg)

    qc.h(count)
    qc.x(work[0])                                  # work register := |1>
    for k in range(n_count):                       # count[k] controls U^(2^k)
        qc.append(c_amod15(a, 2 ** k), [count[k]] + list(work))
    inverse_qft(qc, list(count))
    qc.measure(count, creg)
    return qc


# ======================================================================
# 4. CLASSICAL WRAPPER
# ======================================================================
def find_order(a, backend, N=15, n_count=8, shots=1024):
    """Find the multiplicative order r of a mod N from the measured phase."""
    qc = build_order_finding_circuit(a, n_count)
    counts = run_circuit_and_get_counts(qc, backend, shots)
    for bitstring, _ in sorted(counts.items(), key=lambda kv: -kv[1]):
        measured = int(bitstring, 2)
        phase = measured / (2 ** n_count)
        if phase == 0:
            continue
        r = Fraction(phase).limit_denominator(N).denominator
        if r > 0 and pow(a, r, N) == 1:
            return r, counts
    return None, counts


def shor(backend, N=15, max_tries=12):
    """Factor N. Returns a (p, q) tuple, or None if every attempt failed."""
    if N % 2 == 0:
        return (2, N // 2)
    for _ in range(max_tries):
        a = random.choice([2, 7, 8, 11, 13])
        g = gcd(a, N)
        if g != 1:                                 # lucky guess: a shares a factor
            return (g, N // g)
        r, _ = find_order(a, backend, N)
        if r is None or r % 2 != 0:
            continue
        x = pow(a, r // 2, N)
        if x == N - 1:
            continue
        for cand in (gcd(x - 1, N), gcd(x + 1, N)):
            if cand not in (1, N):
                return (cand, N // cand)
    return None


if __name__ == "__main__":
    random.seed(2)
    backend = AerSimulator()                       # swap for an IBM backend later
    print("Factoring 15 on", backend, "...")
    print("factors:", shor(backend))

    r, counts = find_order(7, backend)
    print("order r for a=7:", r)
    print("top peaks:", sorted(counts.items(), key=lambda kv: -kv[1])[:6])
