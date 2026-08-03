#Quenched spatial disorder variant of the Shor code Monte Carlo pipeline. Instead of one physical error rate p applied uniformly to all 9 data qubits (montecarloShor.py), each qubit gets its own fixed rate p_i drawn once per "disorder realization" from Uniform(p*(1-w), p*(1+w)) (clipped to [0,1]), mean(p_i) = p by construction. That realization's p_i vector is then held fixed while many stochastic noise trials are run against it -- this is the standard quenched-disorder protocol: average over noise (fast/annealed) within a realization, then average over disorder realizations (slow/quenched) on top. Reuses buildFullCircuit/runSingleTrial from montecarloShor.py unchanged (they already accept a length-9 p array), so the circuit, decoder, and success criterion are byte-for-byte identical to the uniform-noise pipeline -- only the noise-injection rates differ. That's what makes the comparison to the uniform threshold meaningful.

# run from src/shor/ -- quick smoke test only, disordersweep.py is the real entry point
# python montecarloShorinhomogeneous.py

import numpy as np

from montecarloShor import runSingleTrial
from shorcode import nQubits


def drawDisorderedRates(rng: np.random.Generator, pMean: float, disorderStrength: float) -> np.ndarray:
    #p_i ~ Uniform(pMean*(1-w), pMean*(1+w)), clipped to [0,1]. w=0 recovers the uniform-noise case exactly.
    lo = pMean * (1 - disorderStrength)
    hi = pMean * (1 + disorderStrength)
    pVec = rng.uniform(lo, hi, size=nQubits)
    return np.clip(pVec, 0.0, 1.0)


def estimateLogicalErrorRateDisordered(
    pMean: float,
    disorderStrength: float,
    nTrials: int,
    nRealizations: int,
    seed: int,
):
    #Split nTrials evenly across nRealizations disorder draws. Returns (overallLer, perRealizationLers, perRealizationPVecs) so callers can separate disorder-to-disorder variance from binomial trial noise.
    if nTrials % nRealizations != 0:
        raise ValueError(f"nTrials ({nTrials}) must be divisible by nRealizations ({nRealizations})")
    trialsPerRealization = nTrials // nRealizations

    disorderRng = np.random.default_rng(seed)
    totalFailures = 0
    perRealizationLers = []
    perRealizationPVecs = []

    for r in range(nRealizations):
        pVec = drawDisorderedRates(disorderRng, pMean, disorderStrength)
        #independent trial-noise stream per realization, derived from the master seed so the whole sweep stays reproducible
        trialRng = np.random.default_rng(seed * 1_000_003 + r)

        failures = 0
        for _ in range(trialsPerRealization):
            if not runSingleTrial(trialRng, pVec):
                failures += 1

        totalFailures += failures
        perRealizationLers.append(failures / trialsPerRealization)
        perRealizationPVecs.append(pVec.tolist())

    overallLer = totalFailures / nTrials
    return overallLer, perRealizationLers, perRealizationPVecs


if __name__ == "__main__":
    print("Smoke test: disordered-noise Shor code LER at small trial counts\n")
    for p in [0.05, 0.10, 0.20]:
        ler, perReal, _ = estimateLogicalErrorRateDisordered(
            p, disorderStrength=0.5, nTrials=200, nRealizations=10, seed=42
        )
        spread = f"[{min(perReal):.3f}, {max(perReal):.3f}]"
        print(f"p_mean={p:.2f}  overall LER={ler:.4f}  per-realization range={spread}")
