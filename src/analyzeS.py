from QECops.simulation import robustnessmetric

#using mean threshold from d=3v5 pair across c values
threshold_by_c = {
    0.1: 0.4635,
    0.2: 0.4274,
    0.3: 0.3523,
    0.4: 0.2976,
    0.5: 0.2381
}

ci_by_c = {
    0.1: {"std": 0.0155},
    0.2: {"std": 0.0210},
    0.3: {"std": 0.0202},
    0.4: {"std": 0.0161},
    0.5: {"std": 0.0301}
}

S = robustnessmetric(threshold_by_c, ci_by_c)

print("Robustness metric S = d(threshold)/d(correlation_strength)")
print("------------------------------------------------------------")
for entry in S:
    print(f"c = {entry['c1']} to {entry['c2']}: S = {entry['S']:.4f} ± {entry.get('S_uncertainty', 'N/A'):.4f}")