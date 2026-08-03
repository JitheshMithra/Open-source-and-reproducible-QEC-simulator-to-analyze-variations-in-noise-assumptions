#Shor's [[9,1,3]] code. 1 logical qubit -> 9 physical qubits, corrects any single-qubit Pauli error. Two 3-qubit repetition codes stacked: outer one catches phase flips (qubit 0 repeated into 0,3,6, then rotated to the +/- basis), inner one catches bit flips (each of those 3 gets its own block of 3). blocks: (0,1,2), (3,4,5), (6,7,8) stabilizers: Z0Z1, Z1Z2, Z3Z4, Z4Z5, Z6Z7, Z7Z8 (bit-flip, per block), X0X1X2X3X4X5, X3X4X5X6X7X8 (phase-flip, between blocks)

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import Statevector, Pauli
import numpy as np


nQubits = 9


def buildEncodingCircuit() -> QuantumCircuit:
    #Encodes whatever's on qubit 0 into the full 9-qubit state. Qubits 1-8 must start at |0>.
    qc = QuantumCircuit(nQubits, name="shor_encode")

    #outer repetition: 0 -> 0,3,6
    qc.cx(0, 3)
    qc.cx(0, 6)

    #rotate the three block leaders before splitting them out - this is what gives us phase protection
    qc.h(0)
    qc.h(3)
    qc.h(6)

    #inner repetition, each leader -> its own block of 3
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.cx(3, 4)
    qc.cx(3, 5)
    qc.cx(6, 7)
    qc.cx(6, 8)

    return qc


def buildSyndromeCircuit() -> QuantumCircuit:
    #9 data qubits + 8 ancillas, measures all 8 stabilizers into an 8-bit classical register. anc 0/1 -> block 0 bit-flip syndrome (Z0Z1, Z1Z2) anc 2/3 -> block 1 anc 4/5 -> block 2 anc 6 -> phase syndrome, blocks 0 vs 1 (X0..X5) anc 7 -> phase syndrome, blocks 1 vs 2 (X3..X8)
    data = QuantumRegister(nQubits, "d")
    anc = QuantumRegister(8, "a")
    creg = ClassicalRegister(8, "syn")
    qc = QuantumCircuit(data, anc, creg, name="shor_syndrome")

    zPairs = [(0, 1, 0), (1, 2, 1), (3, 4, 2), (4, 5, 3), (6, 7, 4), (7, 8, 5)]
    for qi, qj, a in zPairs:
        qc.cx(data[qi], anc[a])
        qc.cx(data[qj], anc[a])

    #X-type stabilizers need the ancilla flipped to the X basis first
    xGroup1 = [0, 1, 2, 3, 4, 5]
    xGroup2 = [3, 4, 5, 6, 7, 8]

    qc.h(anc[6])
    for q in xGroup1:
        qc.cx(anc[6], data[q])
    qc.h(anc[6])

    qc.h(anc[7])
    for q in xGroup2:
        qc.cx(anc[7], data[q])
    qc.h(anc[7])

    qc.measure(anc, creg)

    return qc


def applyError(qc: QuantumCircuit, qubit: int, errorType: str) -> None:
    if errorType == "X":
        qc.x(qubit)
    elif errorType == "Z":
        qc.z(qubit)
    elif errorType == "Y":
        qc.y(qubit)
    else:
        raise ValueError(f"Unknown error type {errorType}")


def decodeBitflipSyndrome(s0: int, s1: int) -> int | None:
    #(s0,s1) from a block's Z_iZ_{i+1} pairs -> which local qubit (0,1,2) flipped
    table = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}
    return table[(s0, s1)]


def decodePhaseSyndrome(s6: int, s7: int) -> int | None:
    #same idea but for the two 6-qubit X stabilizers -> which block had the phase flip
    table = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}
    return table[(s6, s7)]


if __name__ == "__main__":
    enc = buildEncodingCircuit()
    syn = buildSyndromeCircuit()
    print("Encoding circuit depth:", enc.depth(), "gates:", len(enc.data))
    print("Syndrome circuit depth:", syn.depth(), "gates:", len(syn.data))
