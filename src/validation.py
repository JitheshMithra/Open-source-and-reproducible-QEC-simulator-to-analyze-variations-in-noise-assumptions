from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error
from QECops.analytical import analyticalLogicalError
import math

def buildRepetitionCircuit(n, logicalBit=0):
    qc = QuantumCircuit(n, n)
    if logicalBit == 1:
        qc.x(range(n))
    qc.measure(range(n), range(n))
    return qc

def majorityVote(bitstring, n):
    bits = [int(b) for b in bitstring[::-1]]
    return 1 if sum(bits) > n // 2 else 0

def qiskitLer(n, p, shots=100000):
    noiseModel = NoiseModel()
    error = pauli_error([('X', p), ('I', 1 - p)])
    noiseModel.add_all_qubit_quantum_error(error, ['measure'])
    qc = buildRepetitionCircuit(n, logicalBit=0)
    sim = AerSimulator(noise_model=noiseModel)
    result = sim.run(qc, shots=shots, seed_simulator=42).result()
    counts = result.get_counts()
    failures = sum(
        count for bitstring, count in counts.items()
        if majorityVote(bitstring, n) != 0
    )
    ler = failures / shots
    stderr = math.sqrt(ler * (1 - ler) / shots)
    return ler, stderr

def runValidation():
    from QECops.simulation import runTrials
    distances = [3, 5, 7]
    pValues = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    shots = 100000

    print("Qiskit Aer Validation vs Analytical vs QECops")
    print("=" * 70)

    for n in distances:
        print(f"\nCode distance d={n}")
        print(f"{'p':>6}  {'QECops LER':>12}  {'Qiskit LER':>12}  {'Analytical':>12}  {'Max Abs Err':>12}")
        print("-" * 65)

        for idx, p in enumerate(pValues):
            ptSeed = 42 + idx * 99991 + n * 7
            qecops = runTrials(n=n, trials=shots, seed=ptSeed,
                              noiseType="bitflip", p=p)
            qiskitL, _ = qiskitLer(n, p, shots)
            analyticalL = analyticalLogicalError(n, p)
            maxErr = max(
                abs(qecops["LER"] - analyticalL),
                abs(qiskitL - analyticalL),
                abs(qecops["LER"] - qiskitL)
            )
            print(f"{p:>6.2f}  {qecops['LER']:>12.6f}  {qiskitL:>12.6f}  "
                  f"{analyticalL:>12.6f}  {maxErr:>12.6f}")

def runQecopsVsQiskit():
    from QECops.simulation import runTrials
    distances = [3, 5, 7]
    pValues = [0.10, 0.20, 0.30, 0.40, 0.50]
    shots = 100000

    print("\nQECops vs Qiskit Direct Comparison")
    print("=" * 70)

    for n in distances:
        print(f"\nCode distance d={n}")
        print(f"{'p':>6}  {'QECops LER':>12}  {'Qiskit LER':>12}  {'Abs Diff':>12}")
        print("-" * 50)

        for idx, p in enumerate(pValues):
            ptSeed = 42 + idx * 99991 + n * 7
            qecops = runTrials(n=n, trials=shots, seed=ptSeed,
                             noiseType="bitflip", p=p)
            qiskitL, _ = qiskitLer(n, p, shots)
            diff = abs(qecops["LER"] - qiskitL)
            print(f"{p:>6.2f}  {qecops['LER']:>12.6f}  {qiskitL:>12.6f}  {diff:>12.6f}")

if __name__ == "__main__":
    runValidation()
    runQecopsVsQiskit()
