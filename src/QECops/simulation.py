# estimate logical error rate vs p
import random
import math
from .noise import applyNoise
from .decode import majorityDecoder
import numpy as np

def bootstrapThreshold(resultsByDistance, nBootstrap=1000, confidence=0.95):
    # bootstrap confidence intervals on pseudo-threshold estimates, resamples failure counts from Binomial(N, LER) each iteration, finds the crossing each time, and returns a CI across crossings
    rng = np.random.default_rng(0)
    distances = sorted(resultsByDistance.keys())
    out = {}

    for i in range(len(distances) - 1):
        d1, d2 = distances[i], distances[i + 1]
        r1 = resultsByDistance[d1]
        r2 = resultsByDistance[d2]

        pVals = [r["physical_error_rate"] for r in r1]
        l1 = np.array([r["LER"] for r in r1])
        l2 = np.array([r["LER"] for r in r2])
        n1 = np.array([r["trials"] for r in r1])
        n2 = np.array([r["trials"] for r in r2])

        hits = []
        for _ in range(nBootstrap):
            s1 = rng.binomial(n1, l1) / n1
            s2 = rng.binomial(n2, l2) / n2
            diff = s2 - s1
            for j in range(len(diff) - 1):
                if diff[j] * diff[j + 1] < 0:
                    p1, p2 = pVals[j], pVals[j + 1]
                    hits.append(p1 - diff[j] * (p2 - p1) / (diff[j + 1] - diff[j]))
                    break

        if not hits:
            out[(d1, d2)] = None
        else:
            h = np.array(hits)
            a = (1 - confidence) / 2
            out[(d1, d2)] = {
                "mean": float(h.mean()),
                "lower": float(np.quantile(h, a)),
                "upper": float(np.quantile(h, 1 - a)),
                "std": float(h.std())
            }

    return out

# turns one bit into a repetition code
def encodeRepetition(logicalBit, n):
    if logicalBit not in (0, 1):
        raise ValueError("logicalbit must be 0 or 1")
    if n <= 0:
        raise ValueError("n must be positive")
    if n % 2 == 0:
        raise ValueError("n must be odd for majority decoding")

    return [logicalBit] * n

# runs pSweep for multiple code distances
def distanceSweep(distances, pValues, trials, seed, logicalBit=0, noiseType="depolarizing", sweepParam="p", **extraNoiseParams):
    # returns a dict keyed by code distance, values are the pSweep result lists for that distance
    results = {}
    for d in distances:
        results[d] = pSweep(n=d, pValues=pValues, trials=trials, seed=seed, logicalBit=logicalBit, noiseType=noiseType, sweepParam=sweepParam, **extraNoiseParams)
    return results

# monte carlo
def runTrials(n, trials, seed, logicalBit=0, noiseType="depolarizing", **noiseParams):
    # encodes a logical bit into a length-n repetition code, applies noise, decodes with majority vote, and counts logical failures across all trials. LER is failures/trials, stderr is the binomial standard error sqrt(LER*(1-LER)/trials)
    if trials <= 0:
        raise ValueError("trials must be positive")
    rng = random.Random(seed)
    failures = 0

    for _ in range(trials):
        encodedBit = encodeRepetition(logicalBit, n)
        noisedBit = applyNoise(encodedBit, rng, noiseType=noiseType, **noiseParams)
        decodedBit = majorityDecoder(noisedBit)

        if decodedBit != logicalBit:
            failures += 1

    ler = failures / trials
    stderr = math.sqrt(ler * (1 - ler) / trials)

    return {"distance": n, "trials": trials, "failures": failures, "LER": ler, "stderr": stderr, "logicalbit": logicalBit, "noise_type": noiseType, "noise_params": noiseParams}

# special sweep for correlated noise
def correlationSweep(n, p, correlations, trials, seed, logicalBit=0):
    results = []

    for c in correlations:
        result = runTrials(n=n, trials=trials, seed=seed, logicalBit=logicalBit, noiseType="correlated", p=p, correlation=c)
        result["physical_error_rate"] = p
        result["correlation"] = c
        results.append(result)

    return results

# finds where 2 distance curves cross
def estimatePseudoThreshold(resultsByDistance):
    # the pseudo-threshold is the physical error rate below which increasing code distance reduces the logical error rate, estimated by linear interpolation between the two points straddling the crossing of adjacent distance curves. keys of the returned dict are (d_low, d_high) tuples, values are the crossing p (or None if no crossing was found). only meaningful when sweeping p
    thresholds = {}
    distances = sorted(resultsByDistance.keys())

    for i in range(len(distances) - 1):
        dLow = distances[i]
        dHigh = distances[i + 1]
        low = resultsByDistance[dLow]
        high = resultsByDistance[dHigh]
        crossing = None
        for j in range(len(low) - 1):
            p1 = low[j]["physical_error_rate"]
            p2 = low[j + 1]["physical_error_rate"]
            diff1 = high[j]["LER"] - low[j]["LER"]
            diff2 = high[j + 1]["LER"] - low[j + 1]["LER"]
            if diff1 * diff2 < 0:
                crossing = p1 - diff1 * (p2 - p1) / (diff2 - diff1)
                break
        thresholds[(dLow, dHigh)] = crossing
    return thresholds

# sweeps one noise parameter across many values
def pSweep(n, pValues, trials, seed, logicalBit=0, noiseType="depolarizing", sweepParam="p", **extraNoiseParams):
    # each result dict contains LER, stderr, failures, trials, distance, noise_type, noise_params, sweep_param, sweep_value, physical_error_rate
    results = []
    for idx, value in enumerate(pValues):
        noiseParams = extraNoiseParams.copy()
        noiseParams[sweepParam] = value
        ptSeed = seed + idx * 99991 + n * 7
        result = runTrials(n=n, trials=trials, seed=ptSeed, logicalBit=logicalBit, noiseType=noiseType, **noiseParams)
        result["sweep_param"] = sweepParam
        result["sweep_value"] = value

        if "p" in noiseParams:
            result["physical_error_rate"] = noiseParams["p"]
        elif sweepParam in ["px", "pz"]:
            result["physical_error_rate"] = value
        else:
            result["physical_error_rate"] = None
        results.append(result)
    return results

# checks whether larger distance reduces LER
def thresholdScalingSummary(resultsByDistance):
    # for each swept p value, reports the LER at every distance and whether error suppression held across the full distance sweep
    summary = []
    distances = sorted(resultsByDistance.keys())
    pValues = [
        r["physical_error_rate"]
        for r in resultsByDistance[distances[0]]
    ]

    for idx, p in enumerate(pValues):
        lers = {
            d: resultsByDistance[d][idx]["LER"]
            for d in distances
        }

        suppression = all(
            lers[distances[i + 1]] <= lers[distances[i]]
            for i in range(len(distances) - 1)
        )

        summary.append({"physical_error_rate": p, "LERs_by_distance": lers, "error_suppression_with_distance": suppression})

    return summary

def robustnessMetric(thresholdByC, ciByC=None):
    # S = d(threshold)/d(correlation_strength) across swept c values, numerical derivative with optional propagated CI uncertainty
    cVals = sorted(k for k, v in thresholdByC.items() if v is not None)
    out = []

    for i in range(len(cVals) - 1):
        c1, c2 = cVals[i], cVals[i + 1]
        t1, t2 = thresholdByC[c1], thresholdByC[c2]
        S = (t2 - t1) / (c2 - c1)
        entry = {"c_mid": (c1 + c2) / 2, "c1": c1, "c2": c2, "S": S}

        if ciByC:
            s1 = ciByC.get(c1, {}).get("std", 0) if ciByC.get(c1) else 0
            s2 = ciByC.get(c2, {}).get("std", 0) if ciByC.get(c2) else 0
            entry["S_uncertainty"] = float(np.sqrt(s1**2 + s2**2) / abs(c2 - c1))

        out.append(entry)
    return out


def failureBoundary(summary):
    # first p where distance stops helping
    for row in summary:
        if not row["error_suppression_with_distance"]:
            return row["physical_error_rate"]
    return None


def crossingConsistency(thresholds, ci=None):
    pairs = sorted(thresholds.keys())
    out = {}

    for pair in pairs:
        t = thresholds[pair]
        if t is None:
            out[pair] = {"consistent": False, "reason": "no crossing found"}
            continue
        if ci and ci.get(pair) is None:
            out[pair] = {"consistent": False, "reason": "no bootstrap crossing"}
            continue
        out[pair] = {"consistent": True, "reason": "crossing found with CI" if ci else "crossing found"}

    vals = [(p, thresholds[p]) for p in pairs if thresholds[p] is not None]
    if len(vals) >= 2:
        for idx in range(len(vals) - 1):
            p1, t1 = vals[idx]
            p2, t2 = vals[idx + 1]
            if ci and ci.get(p1) and ci.get(p2):
                combinedUncertainty = 1.96 * np.sqrt(ci[p1]["std"]**2 + ci[p2]["std"]**2)
            else:
                combinedUncertainty = 0.05
            if abs(t1 - t2) > combinedUncertainty:
                out[p1]["consistent"] = False
                out[p1]["reason"] = "estimates diverge beyond combined bootstrap uncertainty"
                out[p2]["consistent"] = False
                out[p2]["reason"] = "estimates diverge beyond combined bootstrap uncertainty"
    return out
