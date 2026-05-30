#estimate logical error rate vs p 
import random
import math
from .noise import applynoise
from .decode import majoritydecoder
import numpy as np

def bootstrapthreshold(results_by_distance, nbootstrap=1000, confidence=0.95):
    """
    Bootstrap confidence intervals on pseudo-threshold estimates.
    Resamples LER values within stderr bounds nbootstrap times,
    finds crossing each time, returns CI across crossings.
    """
    rng = np.random.default_rng(0)
    distances = sorted(results_by_distance.keys())
    out = {}

    for i in range(len(distances) - 1):
        d1, d2 = distances[i], distances[i+1]
        r1 = results_by_distance[d1]
        r2 = results_by_distance[d2]

        pvals = [r["physical_error_rate"] for r in r1]
        l1 = np.array([r["LER"] for r in r1])
        l2 = np.array([r["LER"] for r in r2])
        e1 = np.array([r["stderr"] for r in r1])
        e2 = np.array([r["stderr"] for r in r2])

        hits = []
        for _ in range(nbootstrap):
            s1 = np.clip(rng.normal(l1, e1), 0, 1)
            s2 = np.clip(rng.normal(l2, e2), 0, 1)
            diff = s2 - s1
            for j in range(len(diff) - 1):
                if diff[j] * diff[j+1] < 0:
                    p1, p2 = pvals[j], pvals[j+1]
                    hits.append(p1 - diff[j] * (p2-p1) / (diff[j+1]-diff[j]))
                    break

        if not hits:
            out[(d1, d2)] = None
        else:
            h = np.array(hits)
            a = (1 - confidence) / 2
            out[(d1, d2)] = {
                "mean": float(h.mean()),
                "lower": float(np.quantile(h, a)),
                "upper": float(np.quantile(h, 1-a)),
                "std": float(h.std())
            }

    return out
#turns one bit into a repetition code
def encoderepetition(logicalbit, n):
    if logicalbit not in (0, 1):
        raise ValueError("logicalbit must be 0 or 1")
    if n<= 0:
        raise ValueError("n must be positive")
    if n% 2 == 0:
        raise ValueError("n must be odd for majority decoding")

    return [logicalbit] * n
#runs psweep for multiple code distances
def distancesweep(distances,pvalues,trials,seed,logicalbit=0,noisetype="depolarizing",sweepparam="p",**extranoiseparams):
    """
    Run psweep for multiple code distances.
    
    Returns
    -------
    dict
        Keys are code distances, values are psweep result lists.
    """    
    results = {}
    for d in distances:
        results[d] = psweep(n=d,pvalues=pvalues,trials=trials,seed=seed,logicalbit=logicalbit,noisetype=noisetype,sweepparam=sweepparam,**extranoiseparams)
    return results
#monte carlo
def runtrials(n,trials,seed,logicalbit=0,noisetype="depolarizing",**noiseparams):
    """
    Run Monte Carlo simulation of a repetition code under a specified noise model.

    Encodes a logical bit into a length-n repetition code, applies noise,
    decodes with majority vote, and counts logical failures across all trials.

    Parameters
    ----------
    n : int
        Code distance (must be odd positive integer)
    trials : int
        Number of Monte Carlo trials to run
    seed : int
        Random seed for reproducibility
    logicalbit : int, optional
        Logical bit to encode, 0 or 1 (default 0)
    noisetype : str, optional
        Noise model: 'bitflip', 'depolarizing', 'biased', 'correlated' (default 'depolarizing')
    **noiseparams
        Parameters for the noise model, e.g. p=0.1 for bitflip

    Returns
    -------
    dict
        Keys: distance, trials, failures, LER, stderr, logicalbit, noise_type, noise_params
        LER is logical error rate (failures/trials)
        stderr is binomial standard error sqrt(LER*(1-LER)/trials)
    """
    if trials <= 0:
        raise ValueError("trials must be positive")
    rng = random.Random(seed)
    failures = 0

    for _ in range(trials):
        encodedbit = encoderepetition(logicalbit, n)
        noisedbit = applynoise(encodedbit,rng,noisetype=noisetype,**noiseparams)
        decodedbit = majoritydecoder(noisedbit)

        if decodedbit != logicalbit:
            failures += 1

    ler =failures/trials
    stderr= math.sqrt(ler*(1 - ler) /trials)

    return {"distance": n,"trials": trials,"failures": failures,"LER": ler,"stderr": stderr,"logicalbit": logicalbit,"noise_type": noisetype,"noise_params": noiseparams}
#special sweep for correlated noise
def correlationsweep(n,p,correlations,trials,seed,logicalbit=0):
    results = []

    for c in correlations:
        result = runtrials(n=n,trials=trials,seed=seed,logicalbit=logicalbit,noisetype="correlated",p=p,correlation=c)
        result["physical_error_rate"] = p
        result["correlation"] = c
        results.append(result)

    return results

#finds where 2 distance curves cross
def estimatepseudothreshold(results_by_distance):
    """
    Estimate pseudo-threshold by finding where LER curves for adjacent distances cross.

    The pseudo-threshold is the physical error rate below which increasing code
    distance reduces the logical error rate. It is estimated by linear interpolation
    between the two data points straddling the crossing of adjacent distance curves.

    Parameters
    ----------
    results_by_distance : dict
        Keys are code distances (int), values are lists of result dicts
        as returned by psweep. Each result dict must contain
        'physical_error_rate' and 'LER'.

    Returns
    -------
    dict
        Keys are tuples (d_low, d_high) of adjacent distances.
        Values are estimated crossing p values (float), or None if no
        crossing was found in the swept range.

    Notes
    -----
    Only meaningful when sweeping physical error rate p.
    If sweepparam is not 'p', physical_error_rate will be None
    and this function should not be called.
    """
    thresholds = {}
    distances = sorted(results_by_distance.keys())

    for i in range(len(distances) - 1):
        dlow = distances[i]
        dhigh = distances[i + 1]
        low = results_by_distance[dlow]
        high = results_by_distance[dhigh]
        crossing = None
        for j in range(len(low) - 1):
            p1 =low[j]["physical_error_rate"]
            p2= low[j + 1]["physical_error_rate"]
            diff1= high[j]["LER"] - low[j]["LER"]
            diff2 =high[j + 1]["LER"] - low[j + 1]["LER"]
            if diff1*diff2 < 0:
                crossing=p1 -diff1* (p2 -p1)/(diff2- diff1)
                break
        thresholds[(dlow, dhigh)] = crossing
    return thresholds
#Sweeps one noise parameter across many values
def psweep(n,pvalues,trials,seed,logicalbit=0,noisetype="depolarizing",sweepparam="p",**extranoiseparams):
    """
    Sweep one noise parameter across a range of values for a single code distance.

    Parameters
    ----------
    n : int
        Code distance (must be odd positive integer)
    pvalues : list of float
        Values of the sweep parameter to evaluate
    trials : int
        Number of Monte Carlo trials per data point
    seed : int
        Random seed for reproducibility
    logicalbit : int, optional
        Logical bit to encode, 0 or 1 (default 0)
    noisetype : str, optional
        Noise model to apply: 'bitflip', 'depolarizing', 'biased', 'correlated' (default 'depolarizing')
    sweepparam : str, optional
        Name of the parameter being swept, e.g. 'p', 'px', 'correlation' (default 'p')
    **extranoiseparams
        Additional noise parameters passed to the noise model

    Returns
    -------
    list of dict
        Each dict contains: LER, stderr, failures, trials, distance, 
        noise_type, noise_params, sweep_param, sweep_value, physical_error_rate
    """
    results = []
    for value in pvalues:
        noiseparams = extranoiseparams.copy()
        noiseparams[sweepparam] = value
        result = runtrials(n=n,trials=trials,seed=seed,logicalbit=logicalbit,noisetype=noisetype,**noiseparams)
        result["sweep_param"] = sweepparam
        result["sweep_value"] = value

        if "p" in noiseparams:
            result["physical_error_rate"] = noiseparams["p"]
        elif sweepparam in ["px", "pz"]:
            result["physical_error_rate"] = value
        else:
            result["physical_error_rate"] = None
        results.append(result)
    return results
#checks whether larger distance reduces LER
def thresholdscalingsummary(results_by_distance):
    """
    Check whether increasing code distance suppresses logical errors at each p value.
    
    Returns
    -------
    list of dict
        Each dict contains physical_error_rate, LERs_by_distance, 
        and error_suppression_with_distance (bool).
    """
    summary =[]
    distances= sorted(results_by_distance.keys())
    pvalues = [
        r["physical_error_rate"]
        for r in results_by_distance[distances[0]]
    ]

    for idx, p in enumerate(pvalues):
        lers = {
            d: results_by_distance[d][idx]["LER"]
            for d in distances
        }

        suppression = all(
            lers[distances[i + 1]] <= lers[distances[i]]
            for i in range(len(distances) - 1)
        )

        summary.append({"physical_error_rate": p,"LERs_by_distance": lers,"error_suppression_with_distance": suppression})

    return summary

def robustnessmetric(threshold_by_c, ci_by_c=None):
    """
    S = d(threshold)/d(correlation_strength) across swept c values.
    Numerical derivative with optional propagated CI uncertainty.
    """
    cvals = sorted(k for k, v in threshold_by_c.items() if v is not None)
    out = []

    for i in range(len(cvals) - 1):
        c1, c2 = cvals[i], cvals[i+1]
        t1, t2 = threshold_by_c[c1], threshold_by_c[c2]
        S = (t2 - t1) / (c2 - c1)
        entry = {"c_mid": (c1+c2)/2, "c1": c1, "c2": c2, "S": S}

        if ci_by_c:
            s1 = ci_by_c.get(c1, {}).get("std", 0) if ci_by_c.get(c1) else 0
            s2 = ci_by_c.get(c2, {}).get("std", 0) if ci_by_c.get(c2) else 0
            entry["S_uncertainty"] = float(np.sqrt(s1**2 + s2**2) / abs(c2-c1))

        out.append(entry)
    return out


def failureboundary(summary):
    #first p where distance stops helping
    for row in summary:
        if not row["error_suppression_with_distance"]:
            return row["physical_error_rate"]
    return None


def crossingconsistency(thresholds, ci=None):
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

    #flag if pair estimates diverge too much
    vals = [thresholds[p] for p in pairs if thresholds[p] is not None]
    if len(vals) >= 2 and max(vals) - min(vals) > 0.05:
        for pair in pairs:
            if thresholds[pair] is not None:
                out[pair]["consistent"] = False
                out[pair]["reason"] = "estimates diverge across distance pairs"

    return out