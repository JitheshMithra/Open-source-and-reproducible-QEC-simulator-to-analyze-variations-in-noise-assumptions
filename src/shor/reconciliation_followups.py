#Follow-ups #2 and #3 on the rectangular-family reconciliation, see rectangular_family_check_v2.py. Follow-up #1 (literature check) is a writing task, reported directly to Jithesh, not here.
#Follow-up #2: does growing nInner and nOuter together relabel squareFamilyCheck, or is it genuinely different? squareFamilyCheck's q_m(p)=parityN(m,p) is single-level with no inner majority stage, while fxSq/fzSq is a real two-level composition using n at both levels.
#They're algebraically distinct and diverge asymptotically: q_m(p)->0.5 for every p, while fxSq(p,n)->0 as n->inf for p<0.5, so X keeps (even sharpens) its threshold growing both dimensions together, unlike growing nOuter alone or the pure square family.
#fzSq behaves oppositely: its inner stage drifts to 0.5 regardless of p, so the outer majority vote sees a coin flip and the Z-channel still has no threshold, same failure as before.
#Combined total still saturates near 0.5 as n->inf for p<0.5, the same endpoint as the other families, just reached via a different per-channel split (X protected now, Z still degenerate, instead of the reverse under nOuter-only growth).
#Follow-up #3: dense failure-probability table at p=0.10 for nOuter in {9,27,81,243,729} (nInner=3 fixed), X and Z reported separately alongside the combined total, plus a plot (figures/rectangular_family_p010_trend.png).

#run from src/shor/
#python reconciliation_followups.py

from pathlib import Path

import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

from rectangular_family_check_v2 import fxRect, fzRect, parityN, majorityFailN, f, parity3


#Follow-up #2: "grow both together" (n x n) family

def fxSq(p, n):
    return parityN(n, majorityFailN(n, p))


def fzSq(p, n):
    return majorityFailN(n, parityN(n, p))


def squareFamilyOld(p, m):
    #the already-tested squareFamilyCheck formula (concatenation_recursion.py), reproduced here
    #for direct side-by-side comparison. Single-level parity, no inner majority stage.
    return parityN(m, p)


def consistencyCheckAtN3():
    #At n=3, "grow both together" must reduce exactly to the already brute-force-validated base
    #Shor-block formulas (rigorous_model_check.py / rectangular_family_check_v2.py's nOuter=3
    #case), since nInner=nOuter=3 is the base 9-qubit block. This is a free, exact analytic check.
    print("Consistency check: fxSq/fzSq at n=3 must equal the base Shor-block formulas")
    for p in (0.1, 0.3, 0.45, 0.5):
        fxSqVal, fzSqVal = fxSq(p, 3), fzSq(p, 3)
        fxBase, fzBase = fxRect(p, 3), fzRect(p, 3)  #fxRect/fzRect at nOuter=3 is the base block
        matchX = abs(fxSqVal - fxBase) < 1e-12
        matchZ = abs(fzSqVal - fzBase) < 1e-12
        print(f"  p={p:.2f}  fxSq={fxSqVal:.10f} fxBase={fxBase:.10f} match={matchX}   "
              f"fzSq={fzSqVal:.10f} fzBase={fzBase:.10f} match={matchZ}")


def bruteForceSquareFamily(n, pValues, batchSize=1 << 20):
    #Full physical-qubit enumeration for the n x n family, n=5 (25 qubits, 2^25=33554432 patterns), re-deriving both the inner and outer decode from raw error patterns, batched to keep memory bounded. Same decode logic as fullPhysicalBruteForce in rectangular_family_check_v2.py, just with block size n instead of fixed at 3.
    nQubits = n * n
    total = 1 << nQubits
    assert nQubits <= 26, "full enumeration only tractable up to ~26 qubits"
    bitPositions = np.arange(nQubits, dtype=np.int64)

    xProbs = {p: 0.0 for p in pValues}
    zProbs = {p: 0.0 for p in pValues}

    idxStart = 0
    while idxStart < total:
        idxEnd = min(idxStart + batchSize, total)
        idx = np.arange(idxStart, idxEnd, dtype=np.int64)
        bits = ((idx[:, None] >> bitPositions) & 1).astype(np.int8)          #(batch, nQubits)
        bits3d = bits.reshape(-1, n, n)                                      #(batch, nBlocks=n, blockSize=n)
        blockWeight = bits3d.sum(axis=2)                                     #(batch, nBlocks)

        blockXFail = (blockWeight > n // 2)                                  #per-block majority-vote failure
        xNumFail = blockXFail.sum(axis=1)
        xIsFail = (xNumFail % 2 == 1)                                        #outer: parity across blocks

        blockZOdd = (blockWeight % 2 == 1)                                   #per-block parity (degenerate)
        zWeight = blockZOdd.sum(axis=1)
        zIsFail = zWeight > n // 2                                           #outer: majority-vote decode

        w = bits.sum(axis=1)
        for p in pValues:
            probs = p ** w * (1 - p) ** (nQubits - w)
            xProbs[p] += probs[xIsFail].sum()
            zProbs[p] += probs[zIsFail].sum()
        idxStart = idxEnd

    return xProbs, zProbs


def validateSquareFamilyBruteForce():
    print("\nBrute-force validation: n=5 'grow both together' family (25 qubits, 2^25=33,554,432 patterns)")
    pTest = [0.1, 0.3, 0.45]
    bfX, bfZ = bruteForceSquareFamily(5, pTest)
    print(f"\n{'p':>6} {'bruteForce F_X':>15} {'closedForm F_X':>15} {'match':>6}   "
          f"{'bruteForce F_Z':>15} {'closedForm F_Z':>15} {'match':>6}")
    for p in pTest:
        cfX, cfZ = fxSq(p, 5), fzSq(p, 5)
        mX = abs(bfX[p] - cfX) < 1e-9
        mZ = abs(bfZ[p] - cfZ) < 1e-9
        print(f"{p:>6.3f} {bfX[p]:>15.8f} {cfX:>15.8f} {str(mX):>6}   {bfZ[p]:>15.8f} {cfZ:>15.8f} {str(mZ):>6}")


def compareFamilies():
    print("\n'Grow both together' (n x n) vs already-tested pure square family (q_m)")
    print("(pure square family q_m(p)=parityN(m,p) has no inner majority stage at all, shown at m=n^2,")
    print(" i.e. matched to the same total qubit count as the n x n family, for a fair side-by-side)\n")
    pValues = [0.10, 0.30, 0.45, 0.50, 0.55, 0.70]
    nValues = [3, 5, 9, 15, 27, 51, 99]

    header = f"{'p':>6}  " + "  ".join(f"{'n='+str(n)+' (N='+str(n*n)+')':>22}" for n in nValues)
    print(header)
    print("       " + "  ".join(f"{'fxSq/fzSq/qm':>22}" for _ in nValues))
    for p in pValues:
        row = f"{p:>6.3f}  "
        for n in nValues:
            fx, fz = fxSq(p, n), fzSq(p, n)
            qm = squareFamilyOld(p, n * n)
            row += f"{fx:>6.3f}/{fz:>6.3f}/{qm:>5.3f}  "
        print(row)

    print("\nfxSq drops toward 0 for p<0.5 as n grows (X-channel keeps a threshold when both")
    print("dimensions grow together), while qm (pure parity family) drifts to 0.5 regardless of p.")
    print("These are not the same curve, 'grow both together' is a genuinely different, not yet")
    print("previously tested, scaling regime, not a relabeling of squareFamilyCheck. fzSq still")
    print("drifts to 0.5 like qm does, the Z-channel mechanism is the same degenerate-parity drift")
    print("as the pure square family, just relocated to one channel instead of the whole lumped")
    print("code. Combined total F_total_sq still saturates near 0.5 as n->inf for p<0.5")
    print("(1-(1-0)(1-0.5)=0.5), the same endpoint as every other family tested so far, reached")
    print("via a different per-channel split.")


#Follow-up #3: dense p=0.10 table across nOuter in {9,27,81,243,729}, nInner=3 fixed

def denseTrendTable():
    print("\nFollow-up #3: rectangular family (nInner=3 fixed) at p=0.10, dense nOuter sweep")
    p = 0.10
    nOuterValues = [9, 27, 81, 243, 729]
    rows = []
    print(f"\n{'nOuter':>8} {'totalQubits':>12} {'F_X':>14} {'F_Z':>14} {'F_total':>14}")
    for nOuter in nOuterValues:
        fx, fz = fxRect(p, nOuter), fzRect(p, nOuter)
        total = 1 - (1 - fx) * (1 - fz)
        rows.append((nOuter, 3 * nOuter, fx, fz, total))
        print(f"{nOuter:>8} {3*nOuter:>12} {fx:>14.6e} {fz:>14.6e} {total:>14.6f}")

    print("\nF_X (a genuine parity-across-more-blocks check) climbs toward 0.5 as nOuter grows,")
    print("the X-channel degrades monotonically. F_Z (a genuine majority vote across more blocks)")
    print("shrinks toward 0, the Z-channel improves monotonically. F_total is dominated by F_X")
    print("once F_X becomes non-negligible, so total failure gets worse then flattens near 0.5,")
    print("not monotonically better, confirming the reversal is real and visible across all 5 points.")
    return rows


def plotDenseTrend(rows):
    outDir = Path("figures")
    outDir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.labelsize": 12,
                          "legend.fontsize": 10, "figure.dpi": 300})

    nOuterValues = [r[0] for r in rows]
    fxValues = [r[2] for r in rows]
    fzValues = [r[3] for r in rows]
    totalValues = [r[4] for r in rows]

    colors = {"X": "#1f77b4", "Z": "#ff7f0e", "total": "#2ca02c"}

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(nOuterValues, fxValues, marker="o", color=colors["X"], linewidth=1.5, markersize=6, label=r"$F_X$ (degrades)")
    ax.plot(nOuterValues, fzValues, marker="s", color=colors["Z"], linewidth=1.5, markersize=6, label=r"$F_Z$ (improves)")
    ax.plot(nOuterValues, totalValues, marker="^", color=colors["total"], linewidth=1.5, markersize=6, label=r"$F_{total}$ (combined)")
    ax.set_xscale("log")
    ax.set_xticks(nOuterValues)
    ax.set_xticklabels([str(n) for n in nOuterValues])
    ax.set_xlabel(r"$n_{outer}$ (blocks; $n_{inner}=3$ fixed)")
    ax.set_ylabel("Failure probability")
    ax.set_ylim(-0.02, 0.55)
    ax.set_title("Rectangular family, $p=0.10$:\nX degrades, Z improves, total worsens then flattens")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.grid(alpha=0.3, linewidth=0.5)
    fig.tight_layout(rect=(0, 0, 1, 1))
    outPath = outDir / "rectangular_family_p010_trend.png"
    fig.savefig(outPath)
    plt.close(fig)
    print(f"\nSaved plot to {outPath.resolve()}")


def main():
    consistencyCheckAtN3()
    validateSquareFamilyBruteForce()
    compareFamilies()
    rows = denseTrendTable()
    plotDenseTrend(rows)


if __name__ == "__main__":
    main()
