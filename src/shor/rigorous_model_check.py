#Task A: brute-force ground truth for a single Shor 9-qubit block, cross-checked against the toy majority-of-majority model f(f(p)) and the CSS/degeneracy-aware closed form derived here, using the actual decoder tables imported from correctshor.py.
#run from src/shor/
#python rigorous_model_check.py

import itertools
from correctshor import bitflipSyndromeTable, phaseSyndromeTable


def f(p):
    #3-way majority-vote failure prob (at least 2 of 3), iid, the toy model's single level
    return 3 * p ** 2 - 2 * p ** 3


def parity3(p):
    #P(odd number of 3 iid Bernoulli(p) events), the XOR/degenerate combination
    return 3 * p * (1 - p) ** 2 + p ** 3


def toyF(p):
    return f(f(p))


def rigorousFX(p):
    #X-channel: if at least 2 real X errors land in a block, the corrected residual is a 3-qubit "all flipped" pattern equivalent to the global logical X_L, and two block failures cancel, so overall X failure happens iff an odd number of the 3 blocks individually fail, parity of majority.
    return parity3(f(p))


def rigorousFZ(p):
    #Z-channel: no protection within a block (any single Z equals block-level Z_b, pairs cancel), so the per-block "effective Z_b" probability is parity3(p), and the outer code majority-decodes across the 3 blocks, majority of parity.
    return f(parity3(p))


def rigorousFtotal(p):
    #Independent X/Z channels, each qubit has per-channel error prob p; the block fails overall if either channel's logical operator survives.
    fx, fz = rigorousFX(p), rigorousFZ(p)
    return 1 - (1 - fx) * (1 - fz)


def bruteForceX(pValues):
    #Enumerates all 2^9 bit-flip-only patterns, computes each block's syndrome the same way correctshor.py's decoder does, applies the correction, and checks whether the residual is the "all three flipped" pattern; overall X failure happens iff an odd number of the 3 blocks show that residual.
    patterns = list(itertools.product([0, 1], repeat=9))
    #precompute, per pattern, (weight, isFailure)
    weightAndFail = []
    for pat in patterns:
        blockFails = 0
        for b in range(3):
            q0, q1, q2 = pat[3 * b], pat[3 * b + 1], pat[3 * b + 2]
            z01, z12 = q0 ^ q1, q1 ^ q2
            outlier = bitflipSyndromeTable[(z01, z12)]
            corr = [0, 0, 0]
            if outlier is not None:
                corr[outlier] = 1
            residual = (q0 ^ corr[0], q1 ^ corr[1], q2 ^ corr[2])
            if residual == (1, 1, 1):
                blockFails += 1
            elif residual != (0, 0, 0):
                raise AssertionError(f"unexpected residual {residual} for block {pat[3*b:3*b+3]}")
        isFail = (blockFails % 2) == 1
        weightAndFail.append((sum(pat), isFail))

    results = {}
    for p in pValues:
        total = 0.0
        for w, isFail in weightAndFail:
            if isFail:
                total += p ** w * (1 - p) ** (9 - w)
        results[p] = total
    return results


def bruteForceZ(pValues):
    #Enumerates all 2^9 phase-only patterns; per block the "effective Z_b" bit is the parity of its 3 physical Z-bits (what the outer X-type stabilizers actually measure), then decodes via phaseSyndromeTable the same way decodeAndCorrect does and checks whether at least 2 of the 3 blocks end up with residual parity 1.
    patterns = list(itertools.product([0, 1], repeat=9))
    weightAndFail = []
    for pat in patterns:
        b = [pat[0] ^ pat[1] ^ pat[2], pat[3] ^ pat[4] ^ pat[5], pat[6] ^ pat[7] ^ pat[8]]
        xGrp1, xGrp2 = b[0] ^ b[1], b[1] ^ b[2]
        blockWithPhaseError = phaseSyndromeTable[(xGrp1, xGrp2)]
        residual = list(b)
        if blockWithPhaseError is not None:
            residual[blockWithPhaseError] ^= 1
        numOdd = sum(residual)
        isFail = numOdd >= 2
        weightAndFail.append((sum(pat), isFail))

    results = {}
    for p in pValues:
        total = 0.0
        for w, isFail in weightAndFail:
            if isFail:
                total += p ** w * (1 - p) ** (9 - w)
        results[p] = total
    return results


def main():
    pTest = [0.1, 0.3, 0.45, 0.5]

    print("Task A.1/A.2: closed-form comparison")
    print(f"{'p':>6} {'toy F(p)=f(f(p))':>18} {'rigorous F_X(p)':>16} {'rigorous F_Z(p)':>16} {'rigorous F_total(p)':>20}")
    for p in pTest:
        print(f"{p:>6.3f} {toyF(p):>18.6f} {rigorousFX(p):>16.6f} {rigorousFZ(p):>16.6f} {rigorousFtotal(p):>20.6f}")

    print("\nDerivative check at p*=0.5 (analytic, exact via chain rule):")
    print("  f'(p)      = 6p-6p^2            -> f'(0.5)      =", 6*0.5-6*0.25)
    print("  parity3'(p)= 3-12p+12p^2        -> parity3'(0.5)=", 3-12*0.5+12*0.25)
    print("  F_X'(0.5) = f'(parity3(0.5))*parity3'(0.5) ... wait order is parity3(f(p)):")
    print("  F_X'(0.5) = parity3'(f(0.5)) * f'(0.5) =", (3-12*0.5+12*0.25) * (6*0.5-6*0.25))
    print("  F_Z'(0.5) = f'(parity3(0.5)) * parity3'(0.5) =", (6*0.5-6*0.25) * (3-12*0.5+12*0.25))
    print("  Both are exactly 0 because parity3'(0.5)=0 is forced by the chain rule,")
    print("  regardless of the other factor. This is a structural consequence of the")
    print("  degenerate (XOR/parity) combination appearing anywhere in the chain, not a bug.")

    print("\nTask A.3: brute-force ground truth (exact enumeration of 2^9=512 patterns)")
    bfX = bruteForceX(pTest)
    bfZ = bruteForceZ(pTest)
    print(f"\n{'p':>6} {'bruteForce F_X':>15} {'closedForm F_X':>15} {'match?':>7}")
    for p in pTest:
        match = abs(bfX[p] - rigorousFX(p)) < 1e-9
        print(f"{p:>6.3f} {bfX[p]:>15.8f} {rigorousFX(p):>15.8f} {str(match):>7}")

    print(f"\n{'p':>6} {'bruteForce F_Z':>15} {'closedForm F_Z':>15} {'match?':>7}")
    for p in pTest:
        match = abs(bfZ[p] - rigorousFZ(p)) < 1e-9
        print(f"{p:>6.3f} {bfZ[p]:>15.8f} {rigorousFZ(p):>15.8f} {str(match):>7}")

    print(f"\n{'p':>6} {'bruteForce F_total':>18} {'closedForm F_total':>18} {'toy F(p)':>10} {'which matches?':>15}")
    for p in pTest:
        bfTotal = 1 - (1 - bfX[p]) * (1 - bfZ[p])
        closed = rigorousFtotal(p)
        toy = toyF(p)
        matchesRigorous = abs(bfTotal - closed) < 1e-9
        matchesToy = abs(bfTotal - toy) < 1e-6
        verdict = "RIGOROUS" if matchesRigorous and not matchesToy else ("TOY" if matchesToy else "NEITHER")
        print(f"{p:>6.3f} {bfTotal:>18.8f} {closed:>18.8f} {toy:>10.6f} {verdict:>15}")


if __name__ == "__main__":
    main()
