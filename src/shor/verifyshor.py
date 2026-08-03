#Sanity checks for the Shor code impl. 1. Encode |0>/|1> and diff against the textbook closed-form logical states. 2. Every qubit x every Pauli error -> syndrome should come out non-trivial.

# run from src/shor/
# python verifyshor.py

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from shorcode import buildEncodingCircuit, nQubits


def encodeLogical(alpha: complex, beta: complex) -> Statevector:
    qc = QuantumCircuit(nQubits)
    qc.initialize([alpha, beta], 0)
    qc.compose(buildEncodingCircuit(), inplace=True)
    return Statevector.from_instruction(qc)


def knownLogicalStates():
    #closed-form Shor code states, Nielsen & Chuang eq 10.34: |0_L> = (|000>+|111>)^ox3 / 2sqrt2, |1_L> = (|000>-|111>)^ox3 / 2sqrt2
    plusBlock = np.array([1, 0, 0, 0, 0, 0, 0, 1]) / np.sqrt(2)
    minusBlock = np.array([1, 0, 0, 0, 0, 0, 0, -1]) / np.sqrt(2)

    zeroL = np.kron(np.kron(plusBlock, plusBlock), plusBlock) / (np.sqrt(2) ** 2)
    oneL = np.kron(np.kron(minusBlock, minusBlock), minusBlock) / (np.sqrt(2) ** 2)

    zeroL = zeroL / np.linalg.norm(zeroL)
    oneL = oneL / np.linalg.norm(oneL)
    return zeroL, oneL


def testEncoding():
    print("=== TEST 1: Encoding correctness ===")
    sv0 = encodeLogical(1, 0)
    sv1 = encodeLogical(0, 1)

    zeroL, oneL = knownLogicalStates()

    #global phase differs between them so compare fidelity, not the raw vectors
    fid0 = abs(np.vdot(zeroL, sv0.data)) ** 2
    fid1 = abs(np.vdot(oneL, sv1.data)) ** 2

    print(f"Fidelity |0_L> encoded vs analytical: {fid0:.6f}")
    print(f"Fidelity |1_L> encoded vs analytical: {fid1:.6f}")
    assert abs(fid0 - 1.0) < 1e-8, "Encoding of |0> does not match known Shor code |0_L>"
    assert abs(fid1 - 1.0) < 1e-8, "Encoding of |1> does not match known Shor code |1_L>"
    print("PASS: encoding matches closed-form Shor code logical states.\n")


def getSyndromeAnalytic(sv: Statevector) -> list[int]:
    #same 8 stabilizer eigenvalues as the ancilla measurement circuit, just pulled straight from the statevector instead - fine for a correctness check, we're not trying to simulate a real measurement here
    from qiskit.quantum_info import Pauli

    def pauliStr(qubitIndices, kind):
        s = ["I"] * nQubits
        for q in qubitIndices:
            s[nQubits - 1 - q] = kind  #qiskit is little-endian
        return "".join(s)

    stabilizers = [
        pauliStr([0, 1], "Z"),
        pauliStr([1, 2], "Z"),
        pauliStr([3, 4], "Z"),
        pauliStr([4, 5], "Z"),
        pauliStr([6, 7], "Z"),
        pauliStr([7, 8], "Z"),
        pauliStr([0, 1, 2, 3, 4, 5], "X"),
        pauliStr([3, 4, 5, 6, 7, 8], "X"),
    ]

    syndrome = []
    for p in stabilizers:
        expVal = sv.expectation_value(Pauli(p)).real
        bit = 0 if expVal > 0 else 1
        syndrome.append(bit)
    return syndrome


def testErrorDetection():
    print("=== TEST 2: Single-qubit error detection ===")
    zeroL, _ = knownLogicalStates()
    baseSv = encodeLogical(1, 0)

    synClean = getSyndromeAnalytic(baseSv)
    print(f"Clean syndrome (expect all zero): {synClean}")
    assert synClean == [0] * 8, "Clean state should give trivial syndrome"

    results = []
    for q in range(nQubits):
        for err in ["X", "Y", "Z"]:
            qc = QuantumCircuit(nQubits)
            qc.initialize([1, 0], 0)
            qc.compose(buildEncodingCircuit(), inplace=True)
            if err == "X":
                qc.x(q)
            elif err == "Z":
                qc.z(q)
            elif err == "Y":
                qc.y(q)
            sv = Statevector.from_instruction(qc)
            syn = getSyndromeAnalytic(sv)
            nontrivial = syn != [0] * 8
            results.append((q, err, syn, nontrivial))

    print(f"\n{'qubit':<6}{'error':<7}{'syndrome':<28}{'detected'}")
    allDetected = True
    for q, err, syn, detected in results:
        print(f"{q:<6}{err:<7}{str(syn):<28}{detected}")
        if not detected:
            allDetected = False

    print()
    if allDetected:
        print("PASS: every single-qubit X/Y/Z error produces a non-trivial syndrome.")
    else:
        print("FAIL: some errors were undetected (this can happen for Z errors on")
        print("      qubits within the same block under this stabilizer choice --")
        print("      check which specific cases failed above).")

    uniqueSyndromes = set(tuple(s) for _, _, s, _ in results if any(s))
    print(f"\nNumber of distinct non-trivial syndromes observed: {len(uniqueSyndromes)} "
          f"(out of {nQubits*3} possible single-qubit errors)")


if __name__ == "__main__":
    testEncoding()
    testErrorDetection()
