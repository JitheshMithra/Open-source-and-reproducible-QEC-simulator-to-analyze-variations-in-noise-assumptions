from QECops.simulation import robustnessMetric

# (d=3,5) pair
thresholdByC35 = {
    0.1: 0.4609,
    0.2: 0.4220,
    0.3: 0.3663,
    0.4: 0.3172,
    0.5: 0.2372
}

ciByC35 = {
    0.1: {"std": 0.0046},
    0.2: {"std": 0.0054},
    0.3: {"std": 0.0056},
    0.4: {"std": 0.0055},
    0.5: {"std": 0.0088}
}

# (d=5,7) pair
thresholdByC57 = {
    0.1: 0.4752,
    0.2: 0.4346,
    0.3: 0.4026,
    0.4: 0.3351,
    0.5: 0.2844
}

ciByC57 = {
    0.1: {"std": 0.0058},
    0.2: {"std": 0.0066},
    0.3: {"std": 0.0071},
    0.4: {"std": 0.0063},
    0.5: {"std": 0.0076}
}

s35 = robustnessMetric(thresholdByC35, ciByC35)
s57 = robustnessMetric(thresholdByC57, ciByC57)

print("Robustness metric S = d(threshold)/d(correlation_strength)")
print("------------------------------------------------------------")
print("(d=3,5) pair:")
for entry in s35:
    print(f"c = {entry['c1']} to {entry['c2']}: S = {entry['S']:.4f} ± {entry.get('S_uncertainty', 'N/A'):.4f}")

print()
print("(d=5,7) pair:")
for entry in s57:
    print(f"c = {entry['c1']} to {entry['c2']}: S = {entry['S']:.4f} ± {entry.get('S_uncertainty', 'N/A'):.4f}")
