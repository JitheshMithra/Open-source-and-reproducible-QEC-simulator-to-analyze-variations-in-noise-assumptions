from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error
from QECops.analytical import analytical_logical_error
import math

def build_repetition_circuit(n, logicalbit=0):
    qc = QuantumCircuit(n, n)
    if logicalbit == 1:
        qc.x(range(n))
    qc.measure(range(n), range(n))
    return qc

def majority_vote(bitstring, n):
    # qiskit returns bitstring reversed
    bits = [int(b) for b in bitstring[::-1]]
    return 1 if sum(bits) > n // 2 else 0

def qiskit_ler(n, p, shots=100000):
    noise_model = NoiseModel()
    error = pauli_error([('X', p), ('I', 1 - p)])
    noise_model.add_all_qubit_quantum_error(error, ['measure'])

    qc = build_repetition_circuit(n, logicalbit=0)
    sim = AerSimulator(noise_model=noise_model)
    result = sim.run(qc, shots=shots).result()
    counts = result.get_counts()

    failures = sum(
        count for bitstring, count in counts.items()
        if majority_vote(bitstring, n) != 0
    )
    ler = failures / shots
    stderr = math.sqrt(ler * (1 - ler) / shots)
    return ler, stderr

def run_validation():
    distances = [3, 5, 7]
    pvalues = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    shots = 100000

    print("Qiskit Aer Validation vs Analytical vs QECops")
    print("=" * 70)

    for n in distances:
        print(f"\nCode distance d={n}")
        print(f"{'p':>6}  {'Qiskit LER':>12}  {'Analytical':>12}  {'Abs Error':>12}")
        print("-" * 50)

        for p in pvalues:
            qiskit_l, qiskit_se = qiskit_ler(n, p, shots)
            analytical_l = analytical_logical_error(n, p)
            abs_error = abs(qiskit_l - analytical_l)
            print(f"{p:>6.2f}  {qiskit_l:>12.6f}  {analytical_l:>12.6f}  {abs_error:>12.6f}")

def run_qecops_vs_qiskit():
    from QECops.simulation import runtrials
    distances = [3, 5, 7]
    pvalues = [0.10, 0.20, 0.30, 0.40, 0.50]
    shots = 100000

    print("\nQECops vs Qiskit Direct Comparison")
    print("=" * 70)

    for n in distances:
        print(f"\nCode distance d={n}")
        print(f"{'p':>6}  {'QECops LER':>12}  {'Qiskit LER':>12}  {'Abs Diff':>12}")
        print("-" * 50)

        for idx, p in enumerate(pvalues):
            pt_seed = 42 + idx * 99991 + n * 7
            qecops = runtrials(n=n, trials=shots, seed=pt_seed,
                             noisetype="bitflip", p=p)
            qiskit_l, _ = qiskit_ler(n, p, shots)
            diff = abs(qecops["LER"] - qiskit_l)
            print(f"{p:>6.2f}  {qecops['LER']:>12.6f}  {qiskit_l:>12.6f}  {diff:>12.6f}")
def run_validation():
    from QECops.simulation import runtrials
    distances = [3, 5, 7]
    pvalues = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    shots = 100000

    print("Qiskit Aer Validation vs Analytical vs QECops")
    print("=" * 70)

    for n in distances:
        print(f"\nCode distance d={n}")
        print(f"{'p':>6}  {'QECops LER':>12}  {'Qiskit LER':>12}  {'Analytical':>12}  {'Max Abs Err':>12}")
        print("-" * 65)

        for idx, p in enumerate(pvalues):
            pt_seed = 42 + idx * 99991 + n * 7
            qecops = runtrials(n=n, trials=shots, seed=pt_seed,
                              noisetype="bitflip", p=p)
            qiskit_l, _ = qiskit_ler(n, p, shots)
            analytical_l = analytical_logical_error(n, p)
            max_err = max(
                abs(qecops["LER"] - analytical_l),
                abs(qiskit_l - analytical_l),
                abs(qecops["LER"] - qiskit_l)
            )
            print(f"{p:>6.2f}  {qecops['LER']:>12.6f}  {qiskit_l:>12.6f}  "
                  f"{analytical_l:>12.6f}  {max_err:>12.6f}")

if __name__ == "__main__":
    run_validation()
    run_qecops_vs_qiskit()