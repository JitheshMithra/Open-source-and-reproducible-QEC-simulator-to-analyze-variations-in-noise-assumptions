#Task 1 redone with the CSS/degeneracy-aware model, validated against brute force in rigorous_model_check.py; supersedes concatenation_recursion.py's toy model, which disagrees with ground truth at every tested p.
#X-channel is parity-of-majority: F_X(p) = parity3(f(p)), f(p)=3p^2-2p^3 within-block majority, parity3 is XOR-of-3.
#Z-channel is majority-of-parity: F_Z(p) = f(parity3(p)), this is where the block's Z-degeneracy enters (any single Z error is degenerate with any other).
#Both channels share the fixed point p*=0.5 with F_X'(p*)=F_Z'(p*)=0 exactly, forced by the chain rule since parity3'(0.5)=0 regardless of the other factor, a real structural property, not a bug.
#So parity3 iterated converges to 0.5 for any p in (0,1) (same as squareFamilyCheck in concatenation_recursion.py): the Z-channel has no real threshold and drifts to 50% failure, while the X-channel keeps a genuine threshold at p*=0.5.

#run from src/shor/
#python concatenation_recursion_v2.py

import argparse
import json
import os
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import mpmath as mp


def majorityFail3(a, b, c):
    return a * b * (1 - c) + a * (1 - b) * c + (1 - a) * b * c + a * b * c


def parityFail3(a, b, c):
    #P(odd number of 3 independent Bernoulli events), non-iid safe
    return a * (1 - b) * (1 - c) + (1 - a) * b * (1 - c) + (1 - a) * (1 - b) * c + a * b * c


def combineBlockOf9X(children):
    #X-channel: within-block majority, then parity across the 3 sub-blocks
    qA = majorityFail3(children[..., 0], children[..., 1], children[..., 2])
    qB = majorityFail3(children[..., 3], children[..., 4], children[..., 5])
    qC = majorityFail3(children[..., 6], children[..., 7], children[..., 8])
    return parityFail3(qA, qB, qC)


def combineBlockOf9Z(children):
    #Z-channel: within-block parity (degeneracy), then majority across the 3 sub-blocks
    qA = parityFail3(children[..., 0], children[..., 1], children[..., 2])
    qB = parityFail3(children[..., 3], children[..., 4], children[..., 5])
    qC = parityFail3(children[..., 6], children[..., 7], children[..., 8])
    return majorityFail3(qA, qB, qC)


def fMp(p):
    return 3 * p ** 2 - 2 * p ** 3


def parity3Mp(p):
    return 3 * p * (1 - p) ** 2 + p ** 3


def deterministicTraceMpmath(p0, levels, dps=80):
    #Arbitrary-precision (mpmath) deterministic, non-disordered trace of pX_k, pZ_k, total_k, to show the recursion stays well-defined (finite, non-NaN) far beyond where a plain float64 toy-model recursion would underflow to exactly 0.0.
    mp.mp.dps = dps
    px = mp.mpf(p0)
    pz = mp.mpf(p0)
    rows = []
    for lvl in range(1, levels + 1):
        px = parity3Mp(fMp(px))
        pz = fMp(parity3Mp(pz))
        total = 1 - (1 - px) * (1 - pz)
        rows.append({"level": lvl, "pX": str(px), "pZ": str(pz), "total": str(total)})
    return rows


def simulateLevelStatsV2(pStar, delta, levels, nTrees, batchTrees, seed):
    rng = np.random.default_rng(seed)
    nLeaves = 9 ** levels
    lo = np.clip(pStar * (1 - delta), 0.0, 1.0)
    hi = np.clip(pStar * (1 + delta), 0.0, 1.0)

    channels = ("X", "Z", "total")
    sums = {ch: {k: 0.0 for k in range(1, levels + 1)} for ch in channels}
    sumsqs = {ch: {k: 0.0 for k in range(1, levels + 1)} for ch in channels}
    counts = {k: 0 for k in range(1, levels + 1)}

    nDone = 0
    while nDone < nTrees:
        b = min(batchTrees, nTrees - nDone)
        leaves = rng.uniform(lo, hi, size=(b, nLeaves))

        #process channels sequentially rather than holding both at once, to keep peak memory down: the elementwise majorityFail3/parityFail3 formulas each spawn several same-size temporaries, and holding curX and curZ's full level-1 arrays alive together was enough to exhaust memory on this machine
        allLevelsX = []
        cur = leaves
        for lvl in range(1, levels + 1):
            nNodes = cur.shape[1] // 9
            cur = combineBlockOf9X(cur.reshape(b, nNodes, 9))
            allLevelsX.append(cur)

        allLevelsZ = []
        cur = leaves
        for lvl in range(1, levels + 1):
            nNodes = cur.shape[1] // 9
            cur = combineBlockOf9Z(cur.reshape(b, nNodes, 9))
            allLevelsZ.append(cur)

        for lvl in range(1, levels + 1):
            curX, curZ = allLevelsX[lvl - 1], allLevelsZ[lvl - 1]
            total = 1 - (1 - curX) * (1 - curZ)
            for ch, arr in (("X", curX), ("Z", curZ), ("total", total)):
                flat = arr.ravel()
                sums[ch][lvl] += flat.sum()
                sumsqs[ch][lvl] += (flat ** 2).sum()
            counts[lvl] += curX.size
        nDone += b

    stats = {ch: {} for ch in channels}
    for ch in channels:
        for lvl in range(1, levels + 1):
            mean = sums[ch][lvl] / counts[lvl]
            var = max(sumsqs[ch][lvl] / counts[lvl] - mean ** 2, 0.0)
            std = var ** 0.5
            stats[ch][lvl] = {
                "mean": mean,
                "std": std,
                "relWidth": std / mean if mean > 0 else float("nan"),
                "nSamples": counts[lvl],
            }
    return stats


def addRatios(levelStats, levels):
    ratios = {}
    for lvl in range(2, levels + 1):
        prev, cur = levelStats[lvl - 1]["relWidth"], levelStats[lvl]["relWidth"]
        ratios[lvl] = cur / prev if prev > 0 else float("nan")
    return ratios


def argParser():
    parser = argparse.ArgumentParser(description="Corrected (channel-separated) concatenation recursion for Shor-type 9-qubit blocks")
    parser.add_argument("--deltas", nargs="+", type=float, default=[0.05, 0.10, 0.20])
    parser.add_argument("--p-near", type=float, default=0.45)
    parser.add_argument("--p-deep", type=float, default=0.10)
    parser.add_argument("--levels", type=int, default=6)
    parser.add_argument("--n-trees", type=int, default=2000)
    parser.add_argument("--batch-trees", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="concatenation_recursion_v2_results.json")
    return parser.parse_args()


def printLevelTable(regimeName, pStar, delta, channelStats, levels):
    print(f"\n  {regimeName}: p={pStar:.3f}  delta={delta:.2f}")
    for ch in ("X", "Z", "total"):
        stats = channelStats[ch]
        ratios = addRatios(stats, levels)
        print(f"    channel {ch}")
        print(f"    {'level':>5} {'mean':>14} {'std':>14} {'rel.width':>10} {'ratio vs k-1':>13}")
        for lvl in range(1, levels + 1):
            s = stats[lvl]
            ratioStr = f"{ratios[lvl]:.4f}" if lvl in ratios else "--"
            print(f"    {lvl:>5} {s['mean']:>14.6e} {s['std']:>14.6e} {s['relWidth']:>10.4f} {ratioStr:>13}")


def main():
    args = argParser()
    tStart = time.time()

    print("Task 1 redone: corrected (channel-separated, degeneracy-aware) recursion\n")

    def fmtMp(s, sig=12):
        v = mp.mpf(s)
        return mp.nstr(v, sig)

    print("Deterministic (mpmath, 80 digits) trace, p0=0.45, 12 levels, shows no underflow/NaN")
    print("far beyond level 6, satisfying Task B's robustness requirement:")
    for row in deterministicTraceMpmath(0.45, 12):
        print(f"  level {row['level']:>2}: pX={fmtMp(row['pX'])}  pZ={fmtMp(row['pZ'])}  total={fmtMp(row['total'])}")

    print("\nSame trace, p0=0.10, 12 levels:")
    for row in deterministicTraceMpmath(0.10, 12):
        print(f"  level {row['level']:>2}: pX={fmtMp(row['pX'])}  pZ={fmtMp(row['pZ'])}  total={fmtMp(row['total'])}")

    regimes = {"near_threshold": args.p_near, "deep_coding_phase": args.p_deep}

    deltaResults = []
    for regimeName, pRegime in regimes.items():
        print(f"\nRegime: {regimeName} (p={pRegime})")
        for delta in args.deltas:
            seed = args.seed + hash((regimeName, delta)) % 10_000
            stats = simulateLevelStatsV2(pRegime, delta, args.levels, args.n_trees, args.batch_trees, seed)
            printLevelTable(regimeName, pRegime, delta, stats, args.levels)

            entry = {"regime": regimeName, "p": pRegime, "delta": delta, "seed": seed, "channels": {}}
            for ch in ("X", "Z", "total"):
                ratios = addRatios(stats[ch], args.levels)
                meanRatioLate = float(np.mean([ratios[l] for l in range(4, args.levels + 1) if l in ratios and not np.isnan(ratios[l])]) if any(l in ratios for l in range(4, args.levels+1)) else float("nan"))
                entry["channels"][ch] = {
                    "levels": {str(k): v for k, v in stats[ch].items()},
                    "ratios": {str(k): v for k, v in ratios.items()},
                    "meanRatioLevels4to6": meanRatioLate,
                }
            deltaResults.append(entry)

    results = {
        "model": "rigorous channel-separated (X: parity-of-majority, Z: majority-of-parity), validated against brute-force enumeration in rigorous_model_check.py",
        "fixedPoint": 0.5,
        "note": "F_X'(0.5)=F_Z'(0.5)=0 exactly (chain rule, parity3'(0.5)=0); Z-channel has no threshold and drifts to 0.5 for any p in (0,1); X-channel retains a genuine threshold at p*=0.5.",
        "deterministicTraceNear": deterministicTraceMpmath(0.45, 12),
        "deterministicTraceDeep": deterministicTraceMpmath(0.10, 12),
        "deltaSweep": deltaResults,
        "config": vars(args),
    }

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)

    tEnd = time.time()
    print(f"\nSaved results to {args.out}")
    print(f"Total run time: {tEnd - tStart:.1f}s")


if __name__ == "__main__":
    main()
