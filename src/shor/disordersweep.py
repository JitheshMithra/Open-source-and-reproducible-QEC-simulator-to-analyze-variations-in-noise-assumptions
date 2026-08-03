#Disorder comparison sweep: does spatially inhomogeneous noise change the Shor code's logical error rate relative to uniform noise at the same mean p? Same bootstrap CI methodology as thresholdsweep.py, but here the x-axis is disorderStrength at fixed p, run at a few representative p values by default (0.10, 0.15, 0.20, 0.25). N=2000/point since we're now running 3x as many points as thresholdsweep.py for a similar wall-clock budget.

# run from src/shor/
# python disordersweep.py
# python disordersweep.py --trials 4000 --realizations 40 --p 0.10 0.20 0.30 --delta 0.0 0.10 0.20

import argparse
import json
import time

import numpy as np

from montecarloShorinhomogeneous import estimateLogicalErrorRateDisordered


def argParser():
    parser = argparse.ArgumentParser(description="Shor code disorder sweep")

    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--realizations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--p", nargs="+", type=float, default=[0.10, 0.15, 0.20, 0.25])
    parser.add_argument("--delta", nargs="+", type=float, default=[0.0, 0.05, 0.10])
    parser.add_argument("--nbootstrap", type=int, default=1000)
    parser.add_argument("--bootstrapseed", type=int, default=0)
    parser.add_argument("--out", default="disorder_sweep_results.json")

    return parser.parse_args()


def runSweep(pValues, disorderStrengths, nTrials, nRealizations, baseSeed, outPath):
    #same incremental-save/resume pattern as thresholdsweep.py
    results = []

    try:
        with open(outPath) as f:
            results = json.load(f)
        done = {(r["p"], r["disorder_strength"]) for r in results}
        print(f"Resuming: {len(results)} points already complete")
    except FileNotFoundError:
        done = set()

    idx = 0
    for p in pValues:
        for delta in disorderStrengths:
            idx += 1
            if (p, delta) in done:
                continue
            seed = baseSeed + idx * 99991
            t0 = time.time()
            ler, perReal, _ = estimateLogicalErrorRateDisordered(
                p, delta, nTrials=nTrials, nRealizations=nRealizations, seed=seed
            )
            t1 = time.time()
            failures = round(ler * nTrials)
            results.append({
                "p": p,
                "disorder_strength": delta,
                "n_trials": nTrials,
                "n_realizations": nRealizations,
                "failures": failures,
                "ler": ler,
                "seed": seed,
                "wall_time_s": round(t1 - t0, 1),
            })
            print(f"p={p:.2f} delta={delta:.2f}  LER={ler:.4f}  ({failures}/{nTrials})  [{t1-t0:.1f}s]")

            with open(outPath, "w") as f:
                json.dump(results, f, indent=2)

    return results


def bootstrapCi(failures, trials, bootN, seed):
    pHat = failures / trials
    rng = np.random.default_rng(seed)
    bootFailures = rng.binomial(trials, pHat, size=bootN)
    bootLer = bootFailures / trials
    return np.percentile(bootLer, 2.5), np.percentile(bootLer, 97.5)


def addBootstrapCis(results, nBootstrap, bootstrapSeed):
    for r in results:
        lo, hi = bootstrapCi(r["failures"], r["n_trials"], nBootstrap, bootstrapSeed)
        r["ci_low"] = round(lo, 4)
        r["ci_high"] = round(hi, 4)
    return results


def overlaps(a, b):
    return a["ci_low"] <= b["ci_high"] and b["ci_low"] <= a["ci_high"]


def summarize(results):
    print("\n=== Disorder comparison (LER vs disorder strength, at fixed mean p) ===")
    byP = {}
    for r in results:
        byP.setdefault(r["p"], []).append(r)

    for p, rows in sorted(byP.items()):
        rows.sort(key=lambda r: r["disorder_strength"])
        line = " | ".join(
            f"delta={r['disorder_strength']:.2f} LER={r['ler']:.4f} [{r['ci_low']:.4f},{r['ci_high']:.4f}]"
            for r in rows
        )
        print(f"p={p:.2f}: {line}")

        baseline = rows[0]
        for r in rows[1:]:
            tag = "overlap" if overlaps(baseline, r) else "NO OVERLAP"
            print(f"    delta={r['disorder_strength']:.2f} vs baseline: {tag}")


def main():
    args = argParser()

    print(f"Running disorder sweep: {len(args.p)} p-values x {len(args.delta)} disorder strengths, "
          f"N={args.trials} trials/point across {args.realizations} disorder realizations\n")

    tStart = time.time()
    results = runSweep(args.p, args.delta, args.trials, args.realizations, args.seed, args.out)
    tEnd = time.time()
    print(f"\nTotal sweep time: {(tEnd - tStart)/60:.1f} min")

    results = addBootstrapCis(results, args.nbootstrap, args.bootstrapseed)
    summarize(results)

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved raw results to {args.out}")


if __name__ == "__main__":
    main()
