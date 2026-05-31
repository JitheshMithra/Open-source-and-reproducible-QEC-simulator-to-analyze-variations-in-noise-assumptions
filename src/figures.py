#figures.py
#run from src/
#python figures.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from QECops.simulation import distancesweep, bootstrapthreshold

plt.rcParams.update({"font.family": "serif","font.size": 11,"axes.labelsize": 12,"legend.fontsize": 10,"figure.dpi": 300,})

outdir =Path("figures")
outdir.mkdir(exist_ok=True)

distances= [3, 5, 7]
trials =100000
seed =42
colors ={3: "#1f77b4", 5: "#ff7f0e", 7: "#2ca02c"}
markers= {3: "o", 5: "s", 7: "^"}
pvalues= [round(0.05 + i*0.05, 2) for i in range(11)]


#fig 1: LER curves at c=0.1, 0.3, 0.5
print("Running fig 1 simulations...")
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, c in zip(axes, [0.1, 0.3, 0.5]):
    results = distancesweep(
        distances=distances, pvalues=pvalues, trials=trials,
        seed=seed, noisetype="correlated", sweepparam="p", correlation=c)
    ci = bootstrapthreshold(results, nbootstrap=1000, confidence=0.95)
    for d in distances:
        x = [r["physical_error_rate"] for r in results[d]]
        y = [r["LER"] for r in results[d]]
        err = [r["stderr"] for r in results[d]]
        ax.errorbar(x, y, yerr=err, marker=markers[d], color=colors[d],
                   capsize=3, label=f"$d={d}$", linewidth=1.5)
    for pair, val in ci.items():
        if val is not None and pair == (3, 5):
            ax.axvspan(val["lower"], val["upper"], alpha=0.12, color="gray")
    ax.set_title(f"$c = {c}$")
    ax.set_xlabel("Physical error rate $p$")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper left")
axes[0].set_ylabel("Logical error rate $P_L$")
fig.suptitle("Logical error rate curves under correlated noise", fontsize=13)
plt.tight_layout()
plt.savefig(outdir / "fig1_correlated_ler_curves.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 1 done")

#fig 2: sensitivity metric S
c_mids = [0.15, 0.25, 0.35, 0.45]
S_vals = [0.389, 0.557, 0.491, 0.800]  # magnitudes
S_errs = [0.071, 0.078, 0.079, 0.104]
labels = ["0.1→0.2", "0.2→0.3", "0.3→0.4", "0.4→0.5"]
fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(c_mids, S_vals, yerr=S_errs, fmt="o", color="#d62728",
           capsize=5, capthick=1.5, linewidth=1.5, markersize=7)
#shade endpoints to show non-overlap
ax.axhspan(S_vals[0]-S_errs[0], S_vals[0]+S_errs[0], alpha=0.1,
          color="#1f77b4", label="$c=0.1\\to0.2$ range")
ax.axhspan(S_vals[3]-S_errs[3], S_vals[3]+S_errs[3], alpha=0.1,
          color="#ff7f0e", label="$c=0.4\\to0.5$ range")
ax.set_xticks(c_mids)
ax.set_xticklabels(labels, fontsize=10)
ax.set_xlabel("Correlation strength interval")
ax.set_ylabel(r"$|S|$ (sensitivity metric magnitude)")
ax.set_title(r"Sensitivity metric $|S|$ across correlation intervals")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "fig2_sensitivity_S.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 2 done")
# fig 3: inversion at p=0.400
bf = {3: (0.3515, 0.001510), 5: (0.3172, 0.001472), 7: (0.2885, 0.001433)}
co = {3: (0.4388, 0.001569), 5: (0.4496, 0.001573), 7: (0.4481, 0.001573)}
d_vals = [3, 5, 7]

fig, ax = plt.subplots(figsize=(5, 4))
ax.errorbar(d_vals, [bf[d][0] for d in d_vals], yerr=[bf[d][1] for d in d_vals],
           marker="o", color="#1f77b4", capsize=4, linewidth=1.5, label="Bitflip")
ax.errorbar(d_vals, [co[d][0] for d in d_vals], yerr=[co[d][1] for d in d_vals],
           marker="s", color="#d62728", capsize=4, linewidth=1.5,
           label=r"Correlated ($c=0.3$)")
ax.set_xticks(d_vals)
ax.set_xticklabels([f"$d={d}$" for d in d_vals])
ax.set_xlabel("Code distance $d$")
ax.set_ylabel("Logical error rate $P_L$")
ax.set_title(r"Distance scaling at $p=0.400$")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend()
plt.tight_layout()
plt.savefig(outdir / "fig3_inversion_p400.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 3 done")

#fig 5:pseudo-threshold vs correlation strength with 95% CI
#uses confirmed 100k trial threshold estimates

c_vals = [0.1, 0.2, 0.3, 0.4, 0.5]
#d=3v5 thresholds and CIs
t35 = [0.4609, 0.4220, 0.3663, 0.3172, 0.2372]
lo35 = [0.4504, 0.4108, 0.3543, 0.3053, 0.2215]
hi35 = [0.4690, 0.4325, 0.3764, 0.3273, 0.2543]

#d=5v7 thresholds and CIs
t57 = [0.4752, 0.4346, 0.4026, 0.3351, 0.2844]
lo57 = [0.4637, 0.4236, 0.3853, 0.3240, 0.2710]
hi57 = [0.4864, 0.4495, 0.4133, 0.3488, 0.3013]

# bitflip baseline for reference
bf35 = 0.4961
bf57 = 0.5030
fig, ax = plt.subplots(figsize=(6, 4))
#d=3v5
err35_lo = [t - l for t, l in zip(t35, lo35)]
err35_hi = [h - t for t, h in zip(t35, hi35)]
ax.errorbar(c_vals, t35, yerr=[err35_lo, err35_hi],
           fmt="o-", color="#1f77b4", capsize=4, linewidth=1.5,
           label=r"$d=3$ vs $d=5$")
#d=5v7
err57_lo = [t - l for t, l in zip(t57, lo57)]
err57_hi = [h - t for t, h in zip(t57, hi57)]
ax.errorbar(c_vals, t57, yerr=[err57_lo, err57_hi],
           fmt="s--", color="#ff7f0e", capsize=4, linewidth=1.5,
           label=r"$d=5$ vs $d=7$")
#bitflip baseline
ax.axhline(bf35, color="#1f77b4", linestyle=":", linewidth=1,
          alpha=0.6, label="Bitflip baseline (d=3v5)")
ax.axhline(bf57, color="#ff7f0e", linestyle=":", linewidth=1,
          alpha=0.6, label="Bitflip baseline (d=5v7)")
ax.set_xlabel("Correlation strength $c$")
ax.set_ylabel(r"Pseudo-threshold $p^*$")
ax.set_title("Pseudo-threshold vs correlation strength")
ax.set_xticks(c_vals)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "fig5_threshold_vs_c.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 5 done")

print("\nall figures saved to figures/")