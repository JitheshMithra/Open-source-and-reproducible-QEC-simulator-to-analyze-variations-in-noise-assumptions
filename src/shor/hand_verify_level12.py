#Task C: independent hand-verification of concatenation_recursion_v2.py at levels 1-2, p=0.45, delta=0.10.
#First attempt (single-shot Bernoulli sampling) was abandoned: the quenched-disorder std the main script measures (~1e-4 at level 1) is way smaller than a single win/lose trial's sampling noise (std~0.5), so a raw MC would need an infeasible number of samples.
#Instead computes the exact conditional failure probability per disorder realization via brute-force enumeration over all 2^9 error-pattern subsets, generalizing rigorous_model_check.py's bruteForceX/Z to non-iid per-qubit probabilities, exactly where an indexing bug in the main script's combine functions could hide.
#Fully separate code path from the main script, doesn't import majorityFail3/parityFail3/combineBlockOf9X/combineBlockOf9Z, just checks the same physics independently.

#run from src/shor/
#python hand_verify_level12.py

import argparse
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import itertools
import numpy as np

_PATTERNS = np.array(list(itertools.product([0, 1], repeat=9)), dtype=np.int64)  #(512,9)


def _failMask(channel):
    #channel 'X': majority within each 3-group, parity across the 3 groups
    #channel 'Z': parity within each 3-group, majority across the 3 groups
    mask = np.zeros(512, dtype=bool)
    for idx, pat in enumerate(_PATTERNS):
        groupVals = []
        for g in range(3):
            cnt = pat[3 * g] + pat[3 * g + 1] + pat[3 * g + 2]
            if channel == "X":
                groupVals.append(1 if cnt >= 2 else 0)
            else:
                groupVals.append(1 if cnt % 2 == 1 else 0)
        total = sum(groupVals)
        if channel == "X":
            mask[idx] = (total % 2 == 1)
        else:
            mask[idx] = (total >= 2)
    return mask


_MASK_X = _failMask("X")
_MASK_Z = _failMask("Z")


def exactBlockFailProb(p9, mask):
    #p9: (..., 9) array of per-qubit probabilities. Returns (...,) exact failure probability via brute-force weighted sum over all 512 error-pattern subsets (vectorized).
    p9 = np.asarray(p9)
    shape = p9.shape[:-1]
    #subsetProb: (..., 512), built by broadcasting each qubit's p or (1-p) per subset
    subsetProb = np.ones(shape + (512,), dtype=np.float64)
    for i in range(9):
        bit = _PATTERNS[:, i]  #(512,)
        pi = p9[..., i:i+1]    #(...,1)
        factor = np.where(bit == 1, pi, 1 - pi)  #(...,512)
        subsetProb *= factor
    return subsetProb[..., mask].sum(axis=-1)


_CHUNK = 2000  #keep single allocations small, this environment OOMs on modest arrays


def simulateLevel1(pStar, delta, nTrees, seed):
    rng = np.random.default_rng(seed)
    lo, hi = pStar * (1 - delta), pStar * (1 + delta)
    xChunks, zChunks = [], []
    done = 0
    while done < nTrees:
        b = min(_CHUNK, nTrees - done)
        p9 = rng.uniform(lo, hi, size=(b, 9))
        xChunks.append(exactBlockFailProb(p9, _MASK_X))
        zChunks.append(exactBlockFailProb(p9, _MASK_Z))
        done += b
    xProb, zProb = np.concatenate(xChunks), np.concatenate(zChunks)
    total = 1 - (1 - xProb) * (1 - zProb)
    return xProb, zProb, total


def simulateLevel2(pStar, delta, nTrees, seed):
    rng = np.random.default_rng(seed)
    lo, hi = pStar * (1 - delta), pStar * (1 + delta)
    xChunks, zChunks = [], []
    done = 0
    while done < nTrees:
        b = min(_CHUNK, nTrees - done)
        p81 = rng.uniform(lo, hi, size=(b, 9, 9))  #9 level-1 blocks of 9 leaves each
        childXProb = exactBlockFailProb(p81, _MASK_X)  #(b, 9)
        childZProb = exactBlockFailProb(p81, _MASK_Z)  #(b, 9)
        xChunks.append(exactBlockFailProb(childXProb, _MASK_X))  #(b,)
        zChunks.append(exactBlockFailProb(childZProb, _MASK_Z))  #(b,)
        done += b
    xProb, zProb = np.concatenate(xChunks), np.concatenate(zChunks)
    total = 1 - (1 - xProb) * (1 - zProb)
    return xProb, zProb, total


def meanStd(arr):
    return arr.mean(), arr.std()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=0.45)
    parser.add_argument("--delta", type=float, default=0.10)
    parser.add_argument("--n-trees", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    print(f"Task C: independent hand-verification (exact subset enumeration), "
          f"p={args.p}, delta={args.delta}, nTrees={args.n_trees}\n")

    for lvl, simFn in ((1, simulateLevel1), (2, simulateLevel2)):
        xProb, zProb, total = simFn(args.p, args.delta, args.n_trees, args.seed + lvl)
        n = args.n_trees
        print(f"Level {lvl} (n={n} disorder realizations, exact per-realization probability):")
        for name, arr in (("X", xProb), ("Z", zProb), ("total", total)):
            mean, std = meanStd(arr)
            se = std / np.sqrt(n)  #Monte Carlo SE of the mean, across disorder draws
            print(f"  {name:>5}: mean={mean:.8f}  std={std:.8f}  relWidth={std/mean:.4f}  (SE of mean ~= {se:.2e})")
        print()


if __name__ == "__main__":
    main()
