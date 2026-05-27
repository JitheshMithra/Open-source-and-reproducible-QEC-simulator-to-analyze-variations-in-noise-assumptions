#estimate logical error rate vs p 
import random
import math
from .noise import applynoise
from .decode import majoritydecoder

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