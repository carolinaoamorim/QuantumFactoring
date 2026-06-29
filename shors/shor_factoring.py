from fractions import Fraction
from math import ceil, gcd, log2, pi
import random

import numpy as np

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import Operator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2 as Sampler

"""
Shor's algorithm for N = 15 and N = 21.

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

# HARNESS  --  backend-agnostic. Pass AerSimulator() for dev, or an IBM
def run_circuit_and_get_counts(circuit, backend, shots=1024):
    """Run a circuit on any backend (sim or real HW) and return counts."""
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)
    sampler = Sampler(mode=backend)
    job = sampler.run([isa_circuit], shots=shots)
    result = job.result()
    return result[0].data.c.get_counts()   # 'c' = our ClassicalRegister name


# 1. INVERSE QFT
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


# 2. CONTROLLED MODULAR MULTIPLICATION  a^power mod 15
#    For N=15 the controlled-U reduces to a swap network plus X gates.
#    Work-qubit convention: qubit 0 is the LSB of the 4-bit residue.
#    The swaps alone permute |y> -> |2y mod 15> (or the |8y> inverse);
#    the X-all flip on {7,11,13} supplies the extra y -> 15-y those need.
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


def n_work_qubits(N):
    """Work-register width: enough qubits to hold residues 0..N-1."""
    return max(1, ceil(log2(N)))

# 2b. GENERAL CONTROLLED MODULAR MULTIPLICATION  a^power mod N
#    Works for any N (e.g. N=21). The map |y> -> |a^power * y mod N> is a
#    permutation of the residues, so we build it directly as a permutation
#    unitary. a^power is reduced classically first, so a single gate covers
#    even large powers -- no need to repeat the multiply 2^k times.

def c_amod(a, power, N, n_work=None):
    if n_work is None:
        n_work = n_work_qubits(N)
    mult = pow(a, power, N)
    dim = 2 ** n_work
    M = np.zeros((dim, dim))
    for y in range(dim):
        fy = (mult * y) % N if y < N else y        # leave unused basis states fixed
        M[fy, y] = 1.0
    U = QuantumCircuit(n_work)
    U.unitary(Operator(M), range(n_work), label=f"{a}^{power} mod {N}")
    return U.to_gate().control()


# 3. ORDER-FINDING CIRCUIT (phase estimation).
#    n_count counting qubits + ceil(log2 N) work qubits.
#    N=15 uses the hand-built swap network; other N uses the general unitary.
def build_order_finding_circuit(a, N=15, n_count=8):
    n_work = n_work_qubits(N)
    count = QuantumRegister(n_count, "count")
    work = QuantumRegister(n_work, "work")
    creg = ClassicalRegister(n_count, "c")
    qc = QuantumCircuit(count, work, creg)

    qc.h(count)
    qc.x(work[0])                                  # work register := |1>
    for k in range(n_count):                       # count[k] controls U^(2^k)
        if N == 15:
            gate = c_amod15(a, 2 ** k)
        else:
            gate = c_amod(a, 2 ** k, N, n_work)
        qc.append(gate, [count[k]] + list(work))
    inverse_qft(qc, list(count))
    qc.measure(count, creg)
    return qc


# 4. CLASSICAL WRAPPER
def order_from_counts(counts, a, N, n_count):
    """Recover the order r of a mod N from measured phases via continued fractions.
    Tries successively larger denominator caps so divisors of r are caught too."""
    for bits, _ in sorted(counts.items(), key=lambda kv: -kv[1]):
        phase = int(bits, 2) / (2 ** n_count)
        if phase == 0:
            continue
        for cap in range(2, N + 1):
            r = Fraction(phase).limit_denominator(cap).denominator
            if r > 0 and pow(a, r, N) == 1:
                return r
    return None


def find_order(a, backend, N=15, n_count=8, shots=1024):
    """Find the multiplicative order r of a mod N from the measured phase."""
    qc = build_order_finding_circuit(a, N=N, n_count=n_count)
    counts = run_circuit_and_get_counts(qc, backend, shots)
    return order_from_counts(counts, a, N, n_count), counts


def shor(backend, N=15, n_count=8, max_tries=12):
    """Factor N. Returns a (p, q) tuple, or None if every attempt failed."""
    if N % 2 == 0:
        return (2, N // 2)
    for _ in range(max_tries):
        # N=15 only has hand-built gates for these coprime bases
        a = random.choice([2, 7, 8, 11, 13]) if N == 15 else random.randrange(2, N)
        g = gcd(a, N)
        if g != 1:
            return (g, N // g)
        r, _ = find_order(a, backend, N=N, n_count=n_count)
        if r is None or r % 2 != 0:
            continue
        x = pow(a, r // 2, N)
        if x == N - 1:
            continue
        for cand in (gcd(x - 1, N), gcd(x + 1, N)):
            if cand not in (1, N):
                return (cand, N // cand)
    return None


def analyze_base(a, backend, N=21, n_count=3, shots=2000):
    """Order-find and factor with a single base. Returns a summary dict including
    the recovered order, the factors, the circuit depth, and the raw counts."""
    from qiskit import transpile

    qc = build_order_finding_circuit(a, N=N, n_count=n_count)
    depth = transpile(qc, backend).depth()
    counts = run_circuit_and_get_counts(qc, backend, shots)
    r = order_from_counts(counts, a, N, n_count)
    factors = None
    if r and r % 2 == 0:
        x = pow(a, r // 2, N)
        for f in (gcd(x - 1, N), gcd(x + 1, N)):
            if f not in (1, N):
                factors = (f, N // f)
                break
    return {
        "a": a,
        "order": r,
        "factors": factors,
        "n_count": n_count,
        "num_qubits": qc.num_qubits,
        "depth": depth,
        "distinct_outcomes": len(counts),
        "counts": counts,
    }


if __name__ == "__main__":
    random.seed(2)
    backend = AerSimulator()                       # TODO swap for an IBM backend later

    print("Factoring 15 on", backend, "...")
    print("factors:", shor(backend, N=15))
    r, counts = find_order(7, backend, N=15)
    print("order r for a=7 mod 15:", r)
    print("top peaks:", sorted(counts.items(), key=lambda kv: -kv[1])[:6])

    print("\nFactoring 21 on", backend, "...")
    print("factors:", shor(backend, N=21, n_count=8))
    res = analyze_base(8, backend, N=21, n_count=3)
    print("a=8 mod 21 -> order", res["order"], "factors", res["factors"])
