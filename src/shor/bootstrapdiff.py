#Bootstrap-of-difference test for the disorder sweep: is LER(delta) - LER(0) at fixed p distinguishable from zero once we resample at the realization level instead of pretending all N trials are iid? Each (p, delta) point's 20 realizations come from an independently seeded RNG stream (disordersweep.py: seed = baseSeed + idx*99991), so there's no shared realization structure to pair against the delta=0 baseline, making this an unpaired two-sample bootstrap over each group's 20 per-realization LERs.

import argparse
import json

import numpy as np

ALPHA = 0.05
N_COMPARISONS = 8
N_RESAMPLES = 10000
SEED = 0


def argParser():
    parser = argparse.ArgumentParser(description="Bootstrap-of-difference test for the disorder sweep")
    parser.add_argument("--in", dest="inPath", default="disorder_sweep_results.json")
    parser.add_argument("--out", default="disorder_bootstrap_diff_results.json")
    return parser.parse_args()


def bootstrapDiffCi(baselineLers, deltaLers, nResamples, seed, alpha):
    rng = np.random.default_rng(seed)
    baselineLers = np.asarray(baselineLers)
    deltaLers = np.asarray(deltaLers)
    nBase, nDelta = len(baselineLers), len(deltaLers)

    baseSamples = rng.choice(baselineLers, size=(nResamples, nBase), replace=True).mean(axis=1)
    deltaSamples = rng.choice(deltaLers, size=(nResamples, nDelta), replace=True).mean(axis=1)
    diffSamples = deltaSamples - baseSamples

    ciLow, ciHigh = np.percentile(diffSamples, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return diffSamples, ciLow, ciHigh


def main():
    args = argParser()

    with open(args.inPath) as f:
        results = json.load(f)

    if not all("per_realization_ler" in r for r in results):
        raise SystemExit("per_realization_ler missing from results -- rerun disordersweep.py first")

    byP = {}
    for r in results:
        byP.setdefault(r["p"], {})[r["disorder_strength"]] = r

    bonferroniAlpha = ALPHA / N_COMPARISONS

    print(f"Unpaired cluster bootstrap over realizations, N={N_RESAMPLES} resamples/comparison")
    print(f"alpha={ALPHA} (uncorrected) vs Bonferroni alpha={bonferroniAlpha:.5f} across {N_COMPARISONS} comparisons\n")
    header = f"{'p':<6}{'delta':<8}{'obs diff':<12}{'95% CI':<24}{'sig (raw)':<12}{f'{100*(1-bonferroniAlpha):.3f}% CI':<26}{'sig (Bonf.)'}"
    print(header)

    nComparisons = 0
    comparisons = []
    for p, byDelta in sorted(byP.items()):
        baseline = byDelta[0.0]
        baselineLers = baseline["per_realization_ler"]
        for delta, r in sorted(byDelta.items()):
            if delta == 0.0:
                continue
            nComparisons += 1
            deltaLers = r["per_realization_ler"]
            observedDiff = np.mean(deltaLers) - np.mean(baselineLers)

            _, ciLow, ciHigh = bootstrapDiffCi(baselineLers, deltaLers, N_RESAMPLES, SEED, ALPHA)
            sigRaw = not (ciLow <= 0 <= ciHigh)

            _, bfLow, bfHigh = bootstrapDiffCi(baselineLers, deltaLers, N_RESAMPLES, SEED, bonferroniAlpha)
            sigBonf = not (bfLow <= 0 <= bfHigh)

            ciStr = f"[{ciLow:+.4f},{ciHigh:+.4f}]"
            bfStr = f"[{bfLow:+.4f},{bfHigh:+.4f}]"
            diffStr = f"{observedDiff:+.4f}"
            print(f"{p:<6.2f}{delta:<8.2f}{diffStr:<12}{ciStr:<24}{str(sigRaw):<12}{bfStr:<26}{sigBonf}")

            comparisons.append({
                "p": p,
                "disorder_strength": delta,
                "observed_diff": round(float(observedDiff), 6),
                "ci_low": round(float(ciLow), 6),
                "ci_high": round(float(ciHigh), 6),
                "significant_raw_alpha0.05": sigRaw,
                "bonferroni_ci_low": round(float(bfLow), 6),
                "bonferroni_ci_high": round(float(bfHigh), 6),
                "significant_bonferroni": sigBonf,
            })

    if nComparisons != N_COMPARISONS:
        print(f"\nWARNING: expected {N_COMPARISONS} comparisons, found {nComparisons} "
              f"-- Bonferroni denominator in this script no longer matches the data")

    output = {
        "method": "unpaired cluster bootstrap over disorder realizations",
        "n_resamples": N_RESAMPLES,
        "bootstrap_seed": SEED,
        "alpha_raw": ALPHA,
        "n_comparisons": nComparisons,
        "alpha_bonferroni": bonferroniAlpha,
        "comparisons": comparisons,
    }
    with open(args.out, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved bootstrap-of-difference results to {args.out}")


if __name__ == "__main__":
    main()
