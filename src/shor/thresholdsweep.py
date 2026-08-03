#Threshold characterization sweep for the Shor [[9,1,3]] code under uniform independent Pauli noise, with bootstrap confidence intervals. Methodology matches QECops (repetition code): sweep physical error rate p over a range, estimate logical error rate (LER) via Monte Carlo at each p, bootstrap resample the binomial failure counts to get 95% CIs on LER. (Threshold crossing requires a second code distance for comparison -- the Shor code doesn't have a natural "d=9 vs d=?" comparison the way repetition codes do, since [[9,1,3]] is a fixed single code. For this first sweep we characterize the LER curve for the Shor code alone. Comparing against distance scaling will require either concatenating the Shor code with itself, or comparing against the d=3 sub-codes, which is a decision to flag with Pan.) Given the ~36ms/trial cost of 17-qubit circuit simulation, N=3000 trials per point is the default (vs. QECops' N=100,000 for the much cheaper classical repetition code simulation). This gives wider bootstrap CIs, which is an honest and expected tradeoff given the higher per-trial compute cost of genuine quantum circuit simulation.

# run from src/shor/
# python thresholdsweep.py
# python thresholdsweep.py --trials 5000 --pmin 0.05 --pmax 0.40 --pstep 0.05 --seed 7

import argparse
import json
import time

import numpy as np

from montecarloShor import estimateLogicalErrorRate


def argParser():
    parser = argparse.ArgumentParser(description="Shor code threshold sweep")

    parser.add_argument("--trials", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pmin", type=float, default=0.05)
    parser.add_argument("--pmax", type=float, default=0.35)
    parser.add_argument("--pstep", type=float, default=0.05)
    parser.add_argument("--nbootstrap", type=int, default=1000)
    parser.add_argument("--bootstrapseed", type=int, default=0)
    parser.add_argument("--out", default="shor_threshold_sweep_results.json")

    return parser.parse_args()


def pValuesGen(pmin, pmax, pstep):
    if pstep <= 0:
        raise ValueError("pstep must be positive")
    if pmin < 0 or pmax > 1 or pmin > pmax:
        raise ValueError("pmin and pmax must be between 0 and 1, with pmin <= pmax")

    pValues = []
    p = pmin
    while p <= pmax + 1e-12:
        pValues.append(round(p, 12))
        p += pstep
    return pValues


def runSweep(pValues, nTrials, baseSeed, outPath):
    #Run the Monte Carlo sweep, save incrementally after each point so progress is never lost regardless of execution time limits.
    results = []

    #resume from existing partial results if present
    try:
        with open(outPath) as f:
            results = json.load(f)
        donePs = {r["p"] for r in results}
        print(f"Resuming: {len(results)} points already complete: {sorted(donePs)}")
    except FileNotFoundError:
        donePs = set()

    for idx, p in enumerate(pValues):
        if p in donePs:
            continue
        seed = baseSeed + idx * 99991  #same seeding philosophy as QECops
        t0 = time.time()
        ler = estimateLogicalErrorRate(p, nTrials=nTrials, seed=seed)
        t1 = time.time()
        failures = round(ler * nTrials)  #recover raw failure count
        results.append({
            "p": p,
            "n_trials": nTrials,
            "failures": failures,
            "ler": ler,
            "seed": seed,
            "wall_time_s": round(t1 - t0, 1),
        })
        print(f"p={p:.2f}  LER={ler:.4f}  ({failures}/{nTrials} failures)  "
              f"[{t1-t0:.1f}s]")

        #save after EVERY point, not just at the end
        with open(outPath, "w") as f:
            json.dump(results, f, indent=2)

    return results


def bootstrapCi(failures, trials, bootN, seed):
    #Parametric bootstrap CI on the logical error rate, following the same method as QECops: resample failure counts from Binomial(trials, p_hat), compute empirical 2.5th/97.5th percentiles.
    pHat = failures / trials
    rng = np.random.default_rng(seed)
    bootFailures = rng.binomial(trials, pHat, size=bootN)
    bootLer = bootFailures / trials

    ciLow = np.percentile(bootLer, 2.5)
    ciHigh = np.percentile(bootLer, 97.5)
    bootStd = np.std(bootLer)

    return ciLow, ciHigh, bootStd


def addBootstrapCis(results, nBootstrap, bootstrapSeed):
    for r in results:
        ciLow, ciHigh, bootStd = bootstrapCi(
            r["failures"], r["n_trials"], nBootstrap, bootstrapSeed
        )
        r["ci_low"] = round(ciLow, 4)
        r["ci_high"] = round(ciHigh, 4)
        r["boot_std"] = round(bootStd, 4)
    return results


def main():
    args = argParser()
    pValues = pValuesGen(args.pmin, args.pmax, args.pstep)

    print(f"Running Shor code threshold sweep: N={args.trials} trials/point, "
          f"{len(pValues)} points\n")

    tStart = time.time()
    results = runSweep(pValues, args.trials, args.seed, args.out)
    tEnd = time.time()

    print(f"\nTotal sweep time: {tEnd - tStart:.1f}s "
          f"({(tEnd - tStart)/60:.1f} min)")

    results = addBootstrapCis(results, args.nbootstrap, args.bootstrapseed)

    print("\n=== Results with 95% bootstrap CIs ===")
    print(f"{'p':<8}{'LER':<10}{'95% CI':<22}{'boot_std'}")
    for r in results:
        ciStr = f"[{r['ci_low']:.4f}, {r['ci_high']:.4f}]"
        print(f"{r['p']:<8}{r['ler']:<10.4f}{ciStr:<22}{r['boot_std']:.4f}")

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved raw results to {args.out}")


if __name__ == "__main__":
    main()
