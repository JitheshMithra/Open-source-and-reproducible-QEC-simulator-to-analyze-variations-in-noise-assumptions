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
    results = {}
    for d in distances:
        results[d] = psweep(n=d,pvalues=pvalues,trials=trials,seed=seed,logicalbit=logicalbit,noisetype=noisetype,sweepparam=sweepparam,**extranoiseparams)
    return results
#monte carlo
def runtrials(n,trials,seed,logicalbit=0,noisetype="depolarizing",**noiseparams):
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
            if diff1== 0:
                crossing =p1
                break
            if diff1*diff2 < 0:
                crossing=p1 -diff1* (p2 -p1)/(diff2- diff1)
                break
        thresholds[(dlow, dhigh)] = crossing
    return thresholds
#Sweeps one noise parameter across many values
def psweep(n,pvalues,trials,seed,logicalbit=0,noisetype="depolarizing",sweepparam="p",**extranoiseparams):
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