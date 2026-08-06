#Superseded, see concatenation_recursion_v2.py: brute force (rigorous_model_check.py) shows F'(p*)=0, not the ~2 this toy model gives, because the Z-channel is really a parity/XOR check not a majority vote, so kept only as a labeled negative-control comparison.
#Classical (no Qiskit, no circuit noise) recursion of hard-decision logical error rate through concatenation levels of a Shor-type 9-qubit block, modeled as two levels of 3-way majority vote, the standard toy recursion for threshold estimates, deliberately ignoring the code's real Z-error degeneracy since folding that in is what drives F'(p*) to 0.

#run from src/shor/
#python concatenation_recursion.py
#python concatenation_recursion.py --deltas 0.05 0.10 0.20 --n-trees 3000

import argparse
import json
import os
import time

#OpenBLAS's threaded allocator has thrown spurious ArrayMemoryError on modest-sized arrays
#(tens of MB) under this workload. Nothing here is BLAS-bound (pure numpy ufuncs), so forcing
#single-threaded BLAS costs nothing and avoids the crash. Must be set before numpy is imported.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.stats import binom


def majorityFail3(a, b, c):
    #P(at least 2 of 3 independent Bernoulli events), elementwise, works on scalars or arrays
    return a * b * (1 - c) + a * (1 - b) * c + (1 - a) * b * c + a * b * c


def f(p):
    #one 3-way majority-vote level, iid inputs
    return majorityFail3(p, p, p)


def F(p):
    #Shor-type 9-qubit block: majority-of-majority, two concatenated 3-way levels
    return f(f(p))


def combineBlockOf9(children):
    #children: array (...,9) grouped as 3 sub-blocks of 3, matching the Shor code's physical
    #layout 0-2,3-5,6-8. Same majority-of-majority rule as F, but allows the 9 inputs to
    #differ, needed once disorder makes siblings non-identical.
    qA = majorityFail3(children[..., 0], children[..., 1], children[..., 2])
    qB = majorityFail3(children[..., 3], children[..., 4], children[..., 5])
    qC = majorityFail3(children[..., 6], children[..., 7], children[..., 8])
    return majorityFail3(qA, qB, qC)


def findFixedPoint(lo=0.3, hi=0.7, tol=1e-13):
    #bisection on F(p)-p. 0, 1, and 0.5 are all fixed points; bracket around 0.5 to land on
    #the nontrivial threshold instead of the trivial ones at 0/1.
    def g(p):
        return F(p) - p

    gLo, gHi = g(lo), g(hi)
    if gLo * gHi > 0:
        raise ValueError(f"no sign change in [{lo},{hi}]: g(lo)={gLo}, g(hi)={gHi}")
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def numericalDerivative(func, x, h=1e-6):
    return (func(x + h) - func(x - h)) / (2 * h)


def simulateLevelStats(pStar, delta, levels, nTrees, batchTrees, seed):
    #Vectorized bottom-up evaluation: draw all 9**levels leaves for a batch of trees at once,
    #then repeatedly fold groups of 9 into their parent via combineBlockOf9, accumulating a
    #streaming sum/sumsq/count per level. Storing every node at low levels (up to
    #nTrees*9**5 values at level 1) would blow up memory, so only the running moments are kept.
    rng = np.random.default_rng(seed)
    nLeaves = 9 ** levels
    lo = np.clip(pStar * (1 - delta), 0.0, 1.0)
    hi = np.clip(pStar * (1 + delta), 0.0, 1.0)

    sums = {k: 0.0 for k in range(1, levels + 1)}
    sumsqs = {k: 0.0 for k in range(1, levels + 1)}
    counts = {k: 0 for k in range(1, levels + 1)}

    nDone = 0
    while nDone < nTrees:
        b = min(batchTrees, nTrees - nDone)
        current = rng.uniform(lo, hi, size=(b, nLeaves))
        for lvl in range(1, levels + 1):
            nNodes = current.shape[1] // 9
            children = current.reshape(b, nNodes, 9)
            current = combineBlockOf9(children)
            flat = current.ravel()
            sums[lvl] += flat.sum()
            sumsqs[lvl] += (flat ** 2).sum()
            counts[lvl] += flat.size
        nDone += b

    stats = {}
    for lvl in range(1, levels + 1):
        mean = sums[lvl] / counts[lvl]
        var = max(sumsqs[lvl] / counts[lvl] - mean ** 2, 0.0)
        std = var ** 0.5
        stats[lvl] = {
            "mean": mean,
            "std": std,
            "relWidth": std / mean if mean > 0 else float("nan"),
            "nSamples": counts[lvl],
        }
    return stats


def addRatios(stats, levels):
    #ratio of relative width at level k over level k-1, to compare against Pan's guessed 2/3
    #self-averaging factor
    ratios = {}
    for lvl in range(2, levels + 1):
        prev, cur = stats[lvl - 1]["relWidth"], stats[lvl]["relWidth"]
        ratios[lvl] = cur / prev if prev > 0 else float("nan")
    return ratios


def majorityFailN(n, p):
    #P(more than n/2 of n iid Bernoulli(p) trials), n odd, exact via the binomial survival
    #function (safe for large n where direct summation would overflow/underflow)
    thresh = n // 2
    return binom.sf(thresh, n, p)


def squareFamilyCheck(pValues, mValues):
    #Reproduces Pan's formula for the "square" (2n+1)^2-qubit family under a single
    #non-concatenated parity-style decoder: q_m = [1-(1-2p)^m]/2. q_m -> 1/2 as m grows for
    #any p in (0,1), so no p* separates "growing helps" from "growing hurts": this family
    #has no genuine threshold.
    rows = []
    for p in pValues:
        qm = {m: (1 - (1 - 2 * p) ** m) / 2 for m in mValues}
        rows.append({"p": p, "qm": qm})
    return rows


def rectangularFamilyCheck(pValues, nInner, nOuterValues):
    #Generalized (possibly asymmetric) Shor-type block [[nInner*nOuter,1,min(nInner,nOuter)]]
    #with genuine majority-vote decoding at both the inner and outer level (unlike the square
    #family's parity-style q_m). nInner is held fixed while nOuter grows unevenly (9,81,729
    #total qubits via nOuter=3,27,243), the "rectangular" case since nInner != nOuter once
    #nOuter > 3. Confirms a real threshold at p*=0.5 survives as the block grows.
    rows = []
    for p in pValues:
        failProb = {}
        innerFail = majorityFailN(nInner, p)
        for nOuter in nOuterValues:
            outerFail = majorityFailN(nOuter, innerFail)
            failProb[nInner * nOuter] = outerFail
        rows.append({"p": p, "failProb": failProb})
    return rows


def argParser():
    parser = argparse.ArgumentParser(description="Classical concatenation recursion for Shor-type 9-qubit blocks")
    parser.add_argument("--deltas", nargs="+", type=float, default=[0.05, 0.10, 0.20])
    parser.add_argument("--p-near", type=float, default=0.45, help="physical p near threshold p*")
    parser.add_argument("--p-deep", type=float, default=0.10, help="physical p deep in the coding phase")
    parser.add_argument("--levels", type=int, default=6)
    parser.add_argument("--n-trees", type=int, default=2000)
    parser.add_argument("--batch-trees", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="concatenation_recursion_results.json")
    return parser.parse_args()


def printLevelTable(regimeName, pStar, delta, stats, ratios, levels):
    print(f"\n  {regimeName}: p={pStar:.3f}  delta={delta:.2f}")
    print(f"    {'level':>5} {'mean(p_k)':>12} {'std(p_k)':>12} {'rel.width':>10} {'ratio vs k-1':>13}")
    for lvl in range(1, levels + 1):
        s = stats[lvl]
        ratioStr = f"{ratios[lvl]:.4f}" if lvl in ratios else "--"
        print(f"    {lvl:>5} {s['mean']:>12.6f} {s['std']:>12.6f} {s['relWidth']:>10.4f} {ratioStr:>13}")


def main():
    args = argParser()
    tStart = time.time()

    print("Task 1: level-by-level hard-decision concatenation recursion\n")

    pStar = findFixedPoint()
    deriv = numericalDerivative(F, pStar)
    print(f"Nontrivial fixed point p* = {pStar:.10f}")
    print(f"F'(p*) (central finite difference, h=1e-6) = {deriv:.6f}  (Pan recalled ~2)")

    regimes = {"near_threshold": args.p_near, "deep_coding_phase": args.p_deep}

    deltaResults = []
    for regimeName, pRegime in regimes.items():
        print(f"\nRegime: {regimeName} (p={pRegime})")
        for delta in args.deltas:
            seed = args.seed + hash((regimeName, delta)) % 10_000
            stats = simulateLevelStats(pRegime, delta, args.levels, args.n_trees, args.batch_trees, seed)
            ratios = addRatios(stats, args.levels)
            printLevelTable(regimeName, pRegime, delta, stats, ratios, args.levels)

            meanRatioLate = float(np.mean([ratios[l] for l in range(4, args.levels + 1) if l in ratios]))
            deltaResults.append({
                "regime": regimeName,
                "p": pRegime,
                "delta": delta,
                "levels": {str(k): v for k, v in stats.items()},
                "ratios": {str(k): v for k, v in ratios.items()},
                "meanRatioLevels4to6": meanRatioLate,
                "seed": seed,
            })

    print("\n\nTask 1.6: square vs rectangular family threshold comparison")

    squarePs = [0.1, 0.3, 0.5, 0.7, 0.9]
    squareMs = [9, 25, 49, 81, 169, 361, 729]
    squareRows = squareFamilyCheck(squarePs, squareMs)

    print("\n  Square family q_m = [1-(1-2p)^m]/2  (parity-style decoder, single non-concatenated block)")
    header = "    p    | " + " | ".join(f"m={m:>4}" for m in squareMs)
    print(header)
    for row in squareRows:
        line = f"    {row['p']:.2f} | " + " | ".join(f"{row['qm'][m]:.4f}" for m in squareMs)
        print(line)
    print("  all rows drift to 0.5 as m grows, for every p tested: no p* separates 'grows helps' from 'grows hurts'. No genuine threshold.")

    rectPs = [0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70]
    nInner = 3
    nOuterValues = [3, 27, 243]
    rectRows = rectangularFamilyCheck(rectPs, nInner, nOuterValues)

    print(f"\n  Rectangular family (nInner={nInner} fixed, nOuter grows 3->27->243; genuine majority-vote decoding both levels)")
    totalQubits = [nInner * n for n in nOuterValues]
    header = "    p    | " + " | ".join(f"n={n:>4}" for n in totalQubits)
    print(header)
    for row in rectRows:
        line = f"    {row['p']:.2f} | " + " | ".join(f"{row['failProb'][n]:.4e}" for n in totalQubits)
        print(line)
    print("  below p=0.5 the failure probability collapses toward 0 as the block grows; above p=0.5 it grows toward 1. Genuine threshold at p*~0.5 survives block growth.")

    results = {
        "fixedPoint": pStar,
        "derivativeAtFixedPoint": deriv,
        "deltaSweep": deltaResults,
        "squareFamily": {"pValues": squarePs, "mValues": squareMs, "rows": squareRows},
        "rectangularFamily": {
            "pValues": rectPs, "nInner": nInner, "nOuterValues": nOuterValues, "rows": rectRows,
        },
        "config": vars(args),
    }

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)

    tEnd = time.time()
    print(f"\nSaved results to {args.out}")
    print(f"Total run time: {tEnd - tStart:.1f}s")


if __name__ == "__main__":
    main()
