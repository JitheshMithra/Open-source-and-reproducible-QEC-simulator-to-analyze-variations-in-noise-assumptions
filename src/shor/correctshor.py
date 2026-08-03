#Encode -> inject error -> get syndrome -> correct -> check we recovered the logical state. Uses the analytic syndrome extraction from verifyshor.py (statevector expectation values) instead of the physical ancilla circuit in shorcode.py; that one's for later when we actually want to simulate noisy/sampled measurement.

# run from src/shor/
# python correctshor.py

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
import numpy as np

from shorcode import buildEncodingCircuit, nQubits
from verifyshor import knownLogicalStates, getSyndromeAnalytic


#syndrome -> correction lookup. same tables as shorcode.py's decode functions but keeping them local here since decodeAndCorrect needs both together
bitflipSyndromeTable = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}
phaseSyndromeTable = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}

blockQubits = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]


def decodeAndCorrect(syndrome):
    #syndrome = [z01,z12, z34,z45, z67,z78, xGrp1, xGrp2] -> a circuit of X/Z gates on the 9 data qubits that undoes whatever error caused it. Worth noting the phase correction is degenerate: any Z error within a block gives the same syndrome as any other Z error in that block, since the outer code only tells us which block, not which of the 3 qubits. Doesn't matter - fixing the block's first qubit corrects the logical state either way, since in-block Z errors are all equivalent up to a stabilizer.
    correction = QuantumCircuit(nQubits, name="correction")

    z01, z12, z34, z45, z67, z78, xGrp1, xGrp2 = syndrome

    blockSyndromes = [(z01, z12), (z34, z45), (z67, z78)]
    for blockIdx, (s0, s1) in enumerate(blockSyndromes):
        localQ = bitflipSyndromeTable[(s0, s1)]
        if localQ is not None:
            physicalQ = blockQubits[blockIdx][localQ]
            correction.x(physicalQ)

    blockWithPhaseError = phaseSyndromeTable[(xGrp1, xGrp2)]
    if blockWithPhaseError is not None:
        canonicalQ = blockQubits[blockWithPhaseError][0]
        correction.z(canonicalQ)

    return correction


def runFullPipeline(qubit, errorType, verbose=True):
    #Encode |0_L>, hit it with a single-qubit error, correct, return fidelity vs the original.
    qc = QuantumCircuit(nQubits)
    qc.initialize([1, 0], 0)
    qc.compose(buildEncodingCircuit(), inplace=True)

    if errorType == "X":
        qc.x(qubit)
    elif errorType == "Z":
        qc.z(qubit)
    elif errorType == "Y":
        qc.y(qubit)
    else:
        raise ValueError(errorType)

    svWithError = Statevector.from_instruction(qc)
    syndrome = getSyndromeAnalytic(svWithError)

    correction = decodeAndCorrect(syndrome)
    qc.compose(correction, inplace=True)

    svCorrected = Statevector.from_instruction(qc)
    zeroL, _ = knownLogicalStates()
    fidelity = abs(np.vdot(zeroL, svCorrected.data)) ** 2

    if verbose:
        print(f"qubit={qubit} error={errorType:1} syndrome={syndrome} "
              f"-> fidelity after correction = {fidelity:.6f}")

    return fidelity


def testFullCorrection():
    print("=== Full encode -> error -> syndrome -> correct pipeline ===\n")
    fidelities = []
    for q in range(nQubits):
        for err in ["X", "Y", "Z"]:
            fid = runFullPipeline(q, err)
            fidelities.append(fid)

    minFid = min(fidelities)
    print(f"\nMinimum fidelity across all 27 single-qubit errors: {minFid:.6f}")

    if minFid > 1 - 1e-8:
        print("PASS: all 27 single-qubit errors (X, Y, Z on all 9 qubits) "
              "are corrected with fidelity 1.0.")
    else:
        print("FAIL: at least one error was not fully corrected. "
              "Check the case(s) with fidelity < 1 above.")


if __name__ == "__main__":
    testFullCorrection()
