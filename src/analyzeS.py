from QECops.simulation import robustnessmetric

#using mean threshold from d=3v5 pair across c values
threshold_by_c = {
    0.1: 0.4609,
    0.2: 0.4220,
    0.3: 0.3663,
    0.4: 0.3172,
    0.5: 0.2372
}

ci_by_c = {
    0.1: {"std": 0.0046},
    0.2: {"std": 0.0054},
    0.3: {"std": 0.0056},
    0.4: {"std": 0.0055},
    0.5: {"std": 0.0088}
}

S = robustnessmetric(threshold_by_c, ci_by_c)

print("Robustness metric S = d(threshold)/d(correlation_strength)")
print("------------------------------------------------------------")
for entry in S:
    print(f"c = {entry['c1']} to {entry['c2']}: S = {entry['S']:.4f} ± {entry.get('S_uncertainty', 'N/A'):.4f}")