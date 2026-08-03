#Monte Carlo simulation of the Shor code under uniform stochastic Pauli noise. Unlike verifyshor.py (which used exact statevector expectation values to verify correctness deterministically), this module simulates the FULL physical pipeline stochastically: 1. Encode logical |0> 2. Independently apply a random Pauli error to EACH of the 9 qubits with probability p (the "physical error rate"), using a depolarizing-style model: with probability p, one of {X, Y, Z} is applied uniformly at random; with probability 1-p, no error. 3. Measure the 8-qubit syndrome using the PHYSICAL ancilla circuit (qiskit_aer sampling), exactly as a real device would. 4. Decode the syndrome -> apply correction. 5. Check whether the corrected state matches the original logical state (a "success") or not (a "logical failure"). Repeat N times per physical error rate p, estimate the logical error rate as failures/N, following the same methodology as QECops (binomial estimator, bootstrap confidence intervals to follow in the next step).

# run from src/shor/ -- quick smoke test only, thresholdsweep.py is the real entry point
# python montecarloShor.py

import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator

from shorcode import buildEncodingCircuit, nQubits
from correctshor import decodeAndCorrect, blockQubits
from verifyshor import knownLogicalStates


simulator = AerSimulator(method="statevector")

#Precomputed once -- avoid rebuilding this every trial (was ~18ms of pure waste)
_zeroL, _ = knownLogicalStates()


def buildFullCircuit(rng: np.random.Generator, p):
    #Build one full trial circuit: encode -> stochastic noise -> physical syndrome measurement -> (classical decode happens after, outside the circuit, since Qiskit doesn't support mid-circuit classical branching easily in this simple setup) -> return circuit + which errors were actually injected (for optional debugging/validation). p may be a scalar (uniform noise) or a length-9 array (per-qubit error rate, e.g. for spatial inhomogeneity studies).
    qc = QuantumCircuit(nQubits + 8, 8)  #9 data + 8 ancilla, 8 classical bits
    data = list(range(nQubits))
    anc = list(range(nQubits, nQubits + 8))

    #Prepare logical |0>
    qc.initialize([1, 0], 0)
    qc.compose(buildEncodingCircuit(), qubits=data, inplace=True)

    #--- Apply independent stochastic noise to each of the 9 data qubits ---
    pVec = np.broadcast_to(np.asarray(p, dtype=float), (nQubits,))
    injected = []
    for i, q in enumerate(data):
        if rng.random() < pVec[i]:
            err = rng.choice(["X", "Y", "Z"])
            if err == "X":
                qc.x(q)
            elif err == "Y":
                qc.y(q)
            elif err == "Z":
                qc.z(q)
            injected.append((q, err))

    #--- Physical syndrome measurement circuit (same as shorcode.py) ---
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

    return qc, injected


def runSingleTrial(rng: np.random.Generator, p: float) -> bool:
    #Run one Monte Carlo trial. Returns True if the logical state was correctly recovered after correction, False if there was a logical failure.
    qc, injected = buildFullCircuit(rng, p)

    #Run the circuit with measurement to get one shot of the syndrome
    result = simulator.run(qc, shots=1, memory=True).result()
    bitstring = result.get_memory()[0]  #e.g. '10110010', classical bit order

    #Qiskit classical bit order: rightmost char = classical bit 0 We measured anc[0..7] -> creg[0..7], so reverse to get [c0, c1, ..., c7]
    syndrome = [int(b) for b in reversed(bitstring)]

    #Rebuild the pre-measurement state to apply correction and check fidelity. (Re-simulate deterministically with the SAME injected errors, then apply the correction derived from the SAMPLED syndrome -- this checks whether the sampled syndrome, if wrong due to measurement noise model differences, would lead to a correct or incorrect recovery. Here the ancilla measurement is exact/noiseless in the circuit itself, so the sampled syndrome should match the deterministic one.)
    checkQc = QuantumCircuit(nQubits)
    checkQc.initialize([1, 0], 0)
    checkQc.compose(buildEncodingCircuit(), inplace=True)
    for q, err in injected:
        if err == "X":
            checkQc.x(q)
        elif err == "Y":
            checkQc.y(q)
        elif err == "Z":
            checkQc.z(q)

    correction = decodeAndCorrect(syndrome)
    checkQc.compose(correction, inplace=True)

    sv = Statevector.from_instruction(checkQc)
    fidelity = abs(np.vdot(_zeroL, sv.data)) ** 2

    return fidelity > 0.999  #success if recovered (allowing float tolerance)


def estimateLogicalErrorRate(p: float, nTrials: int, seed: int) -> float:
    #Run nTrials Monte Carlo trials at physical error rate p, return LER.
    rng = np.random.default_rng(seed)
    failures = 0
    for _ in range(nTrials):
        success = runSingleTrial(rng, p)
        if not success:
            failures += 1
    return failures / nTrials


if __name__ == "__main__":
    #Quick smoke test at a few error rates with a small trial count first, to confirm the pipeline runs end-to-end before scaling up.
    print("Smoke test: Shor code logical error rate at small trial counts\n")
    testPs = [0.01, 0.05, 0.10, 0.20]
    for p in testPs:
        ler = estimateLogicalErrorRate(p, nTrials=200, seed=42)
        print(f"p = {p:.2f}  ->  logical error rate estimate = {ler:.4f}  (n=200)")
