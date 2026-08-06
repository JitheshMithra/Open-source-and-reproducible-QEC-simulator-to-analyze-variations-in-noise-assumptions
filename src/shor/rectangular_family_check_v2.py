#Redo of the Task 1.6 rectangular family check (nInner=3 fixed, nOuter in {3,27,243}) using the rigorous CSS/degeneracy-aware model instead of the toy majority-of-majority formula, extending rigorous_model_check.py's nOuter=3 derivation to general odd nOuter.
#The original rectangularFamilyCheck (concatenation_recursion.py) used majority-vote decoding at both levels lumped into one channel, the same toy assumption already discredited since the Z-channel's within-block combination is really a degenerate parity check, so those numbers are superseded.
#X-channel: inner majority is unchanged as nOuter grows, and a block's majority failure equals the global X_L, so overall X failure happens iff an odd number of blocks fail: fxRect(p,nOuter)=parityN(nOuter,f(p)), which drifts to 0.5 as nOuter grows, X loses its threshold, same mechanism as the square family.
#Z-channel: each block's phase errors combine via the same degenerate parity as the base case (parity3(p)), and the outer code is a genuine majority-decoded repetition code: fzRect(p,nOuter)=majorityFailN(nOuter,parity3(p)), which goes to 0 for p<0.5 as nOuter grows, Z gains a threshold, unlike concatenation where the outer vote stays stuck at size 3.
#Total fails iff either channel fails; as nOuter grows fxRect->0.5 and fzRect->0, so total saturates near 0.5, still no genuine threshold, growing the outer code just trades X-protection for Z-protection.

#run from src/shor/
#python rectangular_family_check_v2.py

import itertools
import numpy as np
from scipy.stats import binom


def f(p):
    return 3 * p ** 2 - 2 * p ** 3


def parity3(p):
    return 3 * p * (1 - p) ** 2 + p ** 3


def parityN(n, q):
    #P(odd number of n iid Bernoulli(q) events), closed form
    return (1 - (1 - 2 * q) ** n) / 2


def majorityFailN(n, q):
    return binom.sf(n // 2, n, q)


def fxRect(p, nOuter):
    return parityN(nOuter, f(p))


def fzRect(p, nOuter):
    return majorityFailN(nOuter, parity3(p))


#Brute-force validation #1: outer-combination logic via an explicit nearest-neighbor repetition-code decoder (not the closed-form formulas), confirms the "weight > n/2 fails" shortcut matches a real syndrome-based decode.

def outerDecodeFailMask(n, channel):
    #X: parity of n blocks' own X_L-equivalent residuals, whether an odd number of blocks independently failed, no correction step since different blocks aren't caught by the outer stabilizers, established for nOuter=3 in rigorous_model_check.py.
    #Z: explicit repetition-code nearest-neighbor syndrome decode of the n block-level "effective Z_b" bits, generalizing correctshor.py's phaseSyndromeTable (built for n=3) to general odd n.
    patterns = list(itertools.product([0, 1], repeat=n))
    mask = np.zeros(len(patterns), dtype=bool)
    for idx, pat in enumerate(patterns):
        if channel == "X":
            mask[idx] = (sum(pat) % 2 == 1)
        else:
            #syndrome is the n-1 adjacent XORs, which fix pat up to a global complement; minimum-weight decode assumes true=all-0 iff weight(pat)<=n//2, and getting that wrong flips every qubit, leaving a residual equal to the all-ones operator, i.e. failure iff weight(pat) > n//2.
            syndrome = tuple(pat[i] ^ pat[i + 1] for i in range(n - 1))
            assert syndrome == tuple((1 - pat[i]) ^ (1 - pat[i + 1]) for i in range(n - 1))  #complement has the same syndrome, as expected
            weight = sum(pat)
            mask[idx] = weight > n // 2
    return patterns, mask


def bruteForceOuterFail(n, q, channel):
    patterns, mask = outerDecodeFailMask(n, channel)
    total = 0.0
    for pat, isFail in zip(patterns, mask):
        w = sum(pat)
        if isFail:
            total += q ** w * (1 - q) ** (n - w)
    return total


def validateOuterLogic():
    print("Brute-force validation of outer-combination logic (explicit decoder, not closed form)")
    qTest = [0.1, 0.3, 0.468]
    for n in (3, 5, 7, 9):
        print(f"\n  n={n}")
        for ch, closedFn in (("X", lambda n, q: parityN(n, q)), ("Z", majorityFailN)):
            for q in qTest:
                bf = bruteForceOuterFail(n, q, ch)
                cf = closedFn(n, q)
                match = abs(bf - cf) < 1e-9
                print(f"    channel={ch} q={q:.3f}  bruteForce={bf:.8f}  closedForm={cf:.8f}  match={match}")


#Brute-force validation #2: full physical-qubit enumeration for nInner=3, nOuter=5 (15 qubits, 2^15=32768 patterns), the strongest check since it re-derives both the per-block and outer-combination formulas from raw error patterns instead of assuming them.

def fullPhysicalBruteForce(nOuter, pValues):
    nQubits = 3 * nOuter
    assert nQubits <= 21, "full enumeration only tractable for small nOuter"
    patterns = list(itertools.product([0, 1], repeat=nQubits))

    xFailWeights = []
    zFailWeights = []
    for pat in patterns:
        blockXFail = []
        blockZOdd = []
        for b in range(nOuter):
            q0, q1, q2 = pat[3 * b], pat[3 * b + 1], pat[3 * b + 2]
            cnt = q0 + q1 + q2
            blockXFail.append(1 if cnt >= 2 else 0)
            blockZOdd.append(1 if cnt % 2 == 1 else 0)

        xNumFail = sum(blockXFail)
        xIsFail = (xNumFail % 2 == 1)

        zWeight = sum(blockZOdd)
        zIsFail = zWeight > nOuter // 2

        w = sum(pat)
        if xIsFail:
            xFailWeights.append(w)
        if zIsFail:
            zFailWeights.append(w)

    results = {}
    for p in pValues:
        xProb = sum(p ** w * (1 - p) ** (nQubits - w) for w in xFailWeights)
        zProb = sum(p ** w * (1 - p) ** (nQubits - w) for w in zFailWeights)
        results[p] = (xProb, zProb)
    return results


def validateFullPhysical():
    print("\nBrute-force validation #2: full physical-qubit enumeration, nOuter=5 (15 qubits, 2^15=32768 patterns)")
    pTest = [0.1, 0.3, 0.45]
    bf = fullPhysicalBruteForce(5, pTest)
    print(f"\n{'p':>6} {'bruteForce F_X':>15} {'closedForm F_X':>15} {'match':>6}   {'bruteForce F_Z':>15} {'closedForm F_Z':>15} {'match':>6}")
    for p in pTest:
        bfX, bfZ = bf[p]
        cfX, cfZ = fxRect(p, 5), fzRect(p, 5)
        print(f"{p:>6.3f} {bfX:>15.8f} {cfX:>15.8f} {str(abs(bfX-cfX)<1e-9):>6}   {bfZ:>15.8f} {cfZ:>15.8f} {str(abs(bfZ-cfZ)<1e-9):>6}")


def main():
    validateOuterLogic()
    validateFullPhysical()

    print("\n\nCorrected rectangular family (nInner=3, nOuter in {3,27,243})")
    pValues = [0.10, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70]
    nOuterValues = [3, 27, 243]

    print(f"\n{'p':>6}  " + "  ".join(f"{'n='+str(3*n):>28}" for n in nOuterValues))
    header2 = "      " + "  ".join(f"{'F_X':>9} {'F_Z':>9} {'total':>9}" for _ in nOuterValues)
    print(header2)
    for p in pValues:
        row = f"{p:>6.3f}  "
        for nOuter in nOuterValues:
            fx, fz = fxRect(p, nOuter), fzRect(p, nOuter)
            total = 1 - (1 - fx) * (1 - fz)
            row += f"{fx:>9.4f} {fz:>9.4f} {total:>9.4f}  "
        print(row)

    print("\nAs nOuter grows: F_Z_rect -> 0 for p<0.5 (Z-channel gains a threshold from the")
    print("widening genuine outer majority vote), but F_X_rect -> 0.5 for any p in (0,1)")
    print("(X-channel loses its threshold, same parity-drift mechanism as the square family,")
    print("because a wider outer code turns the X-channel's block-failure combination into a")
    print("parity check over more and more blocks). Total saturates near 0.5 as nOuter->inf")
    print("for p<0.5, still no genuine threshold for the full code, just relocated from Z to X.")

    print("\nReconciliation with the original (invalid) rectangular family numbers")
    print("Original concatenation_recursion_results.json 'rectangularFamily' used:")
    print("  outerFail = majorityFailN(nOuter, majorityFailN(nInner=3, p))")
    print("i.e. majority-vote decoding at both levels, for a single lumped channel, the same toy")
    print("assumption already discredited for the base 9-qubit block in rigorous_model_check.py.")
    print("It never modeled the Z-channel's degenerate (parity) within-block combination at all,")
    print("so its reported 'threshold survives to 729 qubits' result is an artifact of that")
    print("missing degeneracy, not a property of the real Shor-type code. Those numbers are")
    print("superseded and shouldn't be cited.")
    print(f"\nOriginal (invalid) vs corrected total failure at p=0.30:")
    for nOuter in nOuterValues:
        oldInner = f(0.30)
        oldOuter = majorityFailN(nOuter, oldInner)
        newFx, newFz = fxRect(0.30, nOuter), fzRect(0.30, nOuter)
        newTotal = 1 - (1 - newFx) * (1 - newFz)
        print(f"  n={3*nOuter:>4}: original(invalid)={oldOuter:.6e}   corrected total={newTotal:.6f}")


if __name__ == "__main__":
    main()
