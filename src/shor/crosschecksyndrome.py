#Cross-check: does the PHYSICAL ancilla-measurement circuit give the same syndrome as the ANALYTIC (statevector expectation value) method, for the same injected errors? This validates that montecarloShor.py's use of the physical circuit is actually correct before we trust it for the threshold sweep.

# run from src/shor/
# python crosschecksyndrome.py

import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

from shorcode import buildEncodingCircuit, nQubits
from verifyshor import getSyndromeAnalytic
from qiskit.quantum_info import Statevector


simulator = AerSimulator(method="statevector")


def getSyndromePhysical(injectedErrors):
    #Build the circuit with SPECIFIC injected errors and measure syndrome via the physical ancilla circuit (deterministic errors, not random).
    qc = QuantumCircuit(nQubits + 8, 8)
    data = list(range(nQubits))
    anc = list(range(nQubits, nQubits + 8))

    qc.initialize([1, 0], 0)
    qc.compose(buildEncodingCircuit(), qubits=data, inplace=True)

    for q, err in injectedErrors:
        if err == "X":
            qc.x(q)
        elif err == "Y":
            qc.y(q)
        elif err == "Z":
            qc.z(q)

    zPairs = [(0, 1, 0), (1, 2, 1), (3, 4, 2), (4, 5, 3), (6, 7, 4), (7, 8, 5)]
    for qi, qj, a in zPairs:
        qc.cx(data[qi], anc[a])
        qc.cx(data[qj], anc[a])

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

    qc.measure(anc, list(range(8)))

    result = simulator.run(qc, shots=1, memory=True).result()
    bitstring = result.get_memory()[0]
    return [int(b) for b in reversed(bitstring)]


def getSyndromeAnalyticForErrors(injectedErrors):
    #Same errors, but via the analytic statevector expectation method.
    qc = QuantumCircuit(nQubits)
    qc.initialize([1, 0], 0)
    qc.compose(buildEncodingCircuit(), inplace=True)
    for q, err in injectedErrors:
        if err == "X":
            qc.x(q)
        elif err == "Y":
            qc.y(q)
        elif err == "Z":
            qc.z(q)
    sv = Statevector.from_instruction(qc)
    return getSyndromeAnalytic(sv)


def crossCheck():
    testCases = [
        [],
        [(0, "X")],
        [(3, "Z")],
        [(5, "Y")],
        [(0, "X"), (4, "Z")],
        [(1, "X"), (2, "X")],  #two X errors in same block -- degenerate case
        [(0, "X"), (3, "X"), (6, "X")],  #one X error per block
        [(2, "Z"), (5, "Z"), (8, "Z")],
    ]

    print(f"{'injected errors':<35}{'physical':<28}{'analytic':<28}{'match'}")
    allMatch = True
    for errs in testCases:
        phys = getSyndromePhysical(errs)
        anal = getSyndromeAnalyticForErrors(errs)
        match = phys == anal
        allMatch = allMatch and match
        print(f"{str(errs):<35}{str(phys):<28}{str(anal):<28}{match}")

    print()
    if allMatch:
        print("PASS: physical (sampled ancilla circuit) syndrome matches "
              "analytic syndrome for all test cases. Monte Carlo pipeline "
              "is using a validated syndrome extraction method.")
    else:
        print("FAIL: mismatch found. Do not trust montecarloShor.py "
              "results until this is fixed.")


if __name__ == "__main__":
    crossCheck()
