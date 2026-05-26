#estimate logical error rate vs p 
import random
import math
from .noise import applynoise
from .decode import majoritydecoder

def encoderepetition(logicalbit, n):
    if logicalbit not in (0, 1):
        raise ValueError("logicalbit must be 0 or 1")
    if n<= 0:
        raise ValueError("n must be positive")
    if n% 2 == 0:
        raise ValueError("n must be odd for majority decoding")

    return [logicalbit] * n


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

def psweep(n,pvalues,trials,seed,logicalbit=0,noisetype="depolarizing",**extranoiseparams):
    results = []
    for p in pvalues:
        noiseparams = extranoiseparams.copy()
        noiseparams["p"] = p
        result = runtrials(n=n,trials=trials,seed=seed,logicalbit=logicalbit,noisetype=noisetype,**noiseparams)
        result["physical_error_rate"]=p
        results.append(result)
    return results


def distancesweep(distances,pvalues,trials,seed,logicalbit=0,noisetype="depolarizing",**extranoiseparams):
    results = {}
    for d in distances:
        results[d] = psweep(n=d,pvalues=pvalues,trials=trials,seed=seed,logicalbit=logicalbit,noisetype=noisetype,**extranoiseparams)
    return results

def correlationsweep(n,p,correlations,trials,seed,logicalbit=0):
    results = []

    for c in correlations:
        result = runtrials(n=n,trials=trials,seed=seed,logicalbit=logicalbit,noisetype="correlated",p=p,correlation=c)
        result["physical_error_rate"] = p
        result["correlation"] = c
        results.append(result)

    return results


def estimatepseudothreshold(results_by_distance):
    thresholds = {}
    distances = sorted(results_by_distance.keys())
    for i in range(len(distances) - 1):
        dlow = distances[i]
        dhigh = distances[i + 1]
        curvelow = results_by_distance[dlow]
        curvehigh = results_by_distance[dhigh]
        crossing = None
        for rlow, rhigh in zip(curvelow, curvehigh):
            p = rlow["physical_error_rate"]

            if rhigh["LER"] <= rlow["LER"]:
                crossing = p
                break

        thresholds[(dlow, dhigh)] = crossing
    return thresholds


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