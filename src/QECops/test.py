from .simulation import (
    distancesweep,
    estimatepseudothreshold,
    thresholdscalingsummary,
)

pvalues = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
distances = [3, 5, 7]

results = distancesweep(
    distances=distances,
    pvalues=pvalues,
    trials=10000,
    seed=42,
    logicalbit=0,
    noisetype="depolarizing",
)

print("\n=== RESULTS ===")
for d, curve in results.items():
    print(f"\nDistance d={d}")
    for r in curve:
        print(
            f"p={r['physical_error_rate']:.3f} "
            f"LER={r['LER']:.6f} "
            f"stderr={r['stderr']:.6f} "
            f"failures={r['failures']}/{r['trials']}"
        )

print("\n=== PSEUDO-THRESHOLDS ===")
thresholds = estimatepseudothreshold(results)
for pair, threshold in thresholds.items():
    print(f"{pair}: {threshold}")

print("\n=== SCALING SUMMARY ===")
summary = thresholdscalingsummary(results)
for row in summary:
    print(row)