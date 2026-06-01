#figures.py
#run from src/
#python figures.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from QECops.simulation import distancesweep, bootstrapthreshold
from QECops.analytical import analytical_logical_error

plt.rcParams.update({"font.family": "serif","font.size": 11,"axes.labelsize": 12,"legend.fontsize": 10,"figure.dpi": 300,})
outdir= Path("figures")
outdir.mkdir(exist_ok=True)
distances = [3, 5, 7]
trials =100000
seed =42
colors ={3: "#1f77b4", 5: "#ff7f0e", 7: "#2ca02c"}
markers ={3: "o", 5: "s", 7: "^"}
pvalues =[round(0.05 + i*0.05, 2) for i in range(11)]

#fig 1: bitflip validation, single panel, no residuals
print("Running fig 1...")
results_bf = distancesweep(distances=distances, pvalues=pvalues, trials=trials,seed=seed, noisetype="bitflip", sweepparam="p")
ci_bf= bootstrapthreshold(results_bf, nbootstrap=1000, confidence=0.95)
fig, ax= plt.subplots(figsize=(6, 4))
for d in distances:
    x = [r["physical_error_rate"] for r in results_bf[d]]
    y = [r["LER"] for r in results_bf[d]]
    err = [r["stderr"] for r in results_bf[d]]
    analytical = [analytical_logical_error(d, p) for p in x]
    ax.errorbar(x, y, yerr=err, marker=markers[d], color=colors[d],
               capsize=3, label=f"$d={d}$ MC", linewidth=1.5, markersize=6)
    ax.plot(x, analytical, linestyle="--", color=colors[d], alpha=0.6,
        label=f"$d={d}$ analytical")
for pair, val in ci_bf.items():
    if val is not None and pair == (3, 5):
        ax.axvspan(val["lower"], val["upper"], alpha=0.1, color="gray",
                  label=r"95% CI $(d=3,d=5)$")
ax.set_xlabel("Physical error rate $p$")
ax.set_ylabel("Logical error rate $P_L$")
ax.set_title("Bitflip noise validation")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "fig1_bitflip_validation.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 1 done")

#fig 2: correlated LER curves at c=0.1, 0.3, 0.5
print("Running fig 2 simulations...")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
c_list = [0.1, 0.3, 0.5]
# confirmed crossing locations from 100k runs
crossings = {0.1: 0.4609, 0.3: 0.3663, 0.5: 0.2372}

for ax, c in zip(axes, c_list):
    results = distancesweep(distances=distances, pvalues=pvalues, trials=trials,seed=seed, noisetype="correlated", sweepparam="p", correlation=c)
    ci = bootstrapthreshold(results, nbootstrap=1000, confidence=0.95)

    for d in distances:
        x = [r["physical_error_rate"] for r in results[d]]
        y = [r["LER"] for r in results[d]]
        err = [r["stderr"] for r in results[d]]
        ax.errorbar(x, y, yerr=err, marker=markers[d], color=colors[d],
                   capsize=3, label=f"$d={d}$", linewidth=2.0, markersize=7)
    for pair, val in ci.items():
        if val is not None and pair == (3, 5):
            ax.axvspan(val["lower"], val["upper"], alpha=0.12, color="gray")
#annotate crossing in each correlated-noise panel
    ann_y = {0.1: 0.43, 0.3: 0.38, 0.5: 0.35}
    text_y = {0.1: 0.32, 0.3: 0.30, 0.5: 0.28}

    ax.axvline(crossings[c], color="black", linestyle=":", linewidth=1.2, alpha=0.7)
    ax.annotate(
        f"$p^* = {crossings[c]}$",
        xy=(crossings[c], ann_y[c]),
        xytext=(crossings[c] + 0.04, text_y[c]),
        fontsize=9,
        color="black",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    ax.set_title(f"$c = {c}$")
    ax.set_xlabel("Physical error rate $p$")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize=9)
axes[0].set_ylabel("Logical error rate $P_L$")
fig.suptitle("Logical error rate curves under correlated noise", fontsize=13)
plt.tight_layout()
plt.savefig(outdir / "fig2_correlated_ler_curves.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 2 done")
#fig 3: pseudo-threshold vs correlation strength (main result figure)
c_vals = [0.1, 0.2, 0.3, 0.4, 0.5]
t35 =  [0.4609, 0.4220, 0.3663, 0.3172, 0.2372]
lo35 = [0.4504, 0.4108, 0.3543, 0.3053, 0.2215]
hi35 = [0.4690, 0.4325, 0.3764, 0.3273, 0.2543]
t57 =  [0.4752, 0.4346, 0.4026, 0.3351, 0.2844]
lo57 = [0.4637, 0.4236, 0.3853, 0.3240, 0.2710]
hi57 = [0.4864, 0.4495, 0.4133, 0.3488, 0.3013]
bf35, bf57 = 0.4961, 0.5030
fig, ax = plt.subplots(figsize=(6, 4))
err35 = [[t-l for t,l in zip(t35,lo35)], [h-t for t,h in zip(t35,hi35)]]
err57 = [[t-l for t,l in zip(t57,lo57)], [h-t for t,h in zip(t57,hi57)]]
ax.errorbar(c_vals, t35, yerr=err35, fmt="o-", color="#1f77b4",
           capsize=4, linewidth=2.0, markersize=7, label=r"$d=3$ vs $d=5$")
ax.errorbar(c_vals, t57, yerr=err57, fmt="s--", color="#ff7f0e",
           capsize=4, linewidth=2.0, markersize=7, label=r"$d=5$ vs $d=7$")
ax.axhline(bf35, color="#1f77b4", linestyle=":", linewidth=1.2,
          alpha=0.5, label=r"Bitflip baseline $(d=3,d=5)$")
ax.axhline(bf57, color="#ff7f0e", linestyle=":", linewidth=1.2,
          alpha=0.5, label=r"Bitflip baseline $(d=5,d=7)$")
ax.set_xlabel("Correlation strength $c$")
ax.set_ylabel(r"Pseudo-threshold $p^*$")
ax.set_title("Pseudo-threshold vs correlation strength")
ax.set_xticks(c_vals)
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "fig3_threshold_vs_c.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 3 done")

#fig 4: sensitivity metric S, faded middle intervals
c_mids =  [0.15, 0.25, 0.35, 0.45]
S_vals =  [0.389, 0.557, 0.491, 0.800]
S_errs =  [0.071, 0.078, 0.079, 0.104]
labels =  ["0.1→0.2", "0.2→0.3", "0.3→0.4", "0.4→0.5"]
#endpoints solid, middle faded
alphas = [1.0, 0.35, 0.35, 1.0]
point_colors = ["#d62728", "#aaaaaa", "#aaaaaa", "#d62728"]
fig, ax = plt.subplots(figsize=(6, 4))
for i in range(4):
    ax.errorbar(c_mids[i], S_vals[i], yerr=S_errs[i],
               fmt="o", color=point_colors[i], capsize=5,
               capthick=1.5, linewidth=1.5, markersize=8,
               alpha=alphas[i])
#shade endpoint ranges
ax.axhspan(S_vals[0]-S_errs[0], S_vals[0]+S_errs[0],
          alpha=0.08, color="#1f77b4", label="$c=0.1\\to0.2$ range")
ax.axhspan(S_vals[3]-S_errs[3], S_vals[3]+S_errs[3],
          alpha=0.08, color="#ff7f0e", label="$c=0.4\\to0.5$ range")
ax.set_xticks(c_mids)
ax.set_xticklabels(labels, fontsize=10)
ax.set_xlabel("Correlation strength interval")
ax.set_ylabel(r"$|S|$ (sensitivity metric magnitude)")
ax.set_title(r"Sensitivity metric $|S|$ across correlation intervals")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "fig4_sensitivity_S.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 4 done")

#fig 5: inversion at p=0.400
bf =  {3: (0.3515, 0.001510), 5: (0.3172, 0.001472), 7: (0.2885, 0.001433)}
co =  {3: (0.4388, 0.001569), 5: (0.4496, 0.001573), 7: (0.4481, 0.001573)}
d_vals = [3, 5, 7]

fig, ax = plt.subplots(figsize=(5, 4))
ax.errorbar(d_vals, [bf[d][0] for d in d_vals],
           yerr=[bf[d][1] for d in d_vals],
           marker="o", color="#1f77b4", capsize=4,
           linewidth=2.0, markersize=7, label="Bitflip")
ax.errorbar(d_vals, [co[d][0] for d in d_vals],
           yerr=[co[d][1] for d in d_vals],
           marker="s", color="#d62728", capsize=4,
           linewidth=2.0, markersize=7,
           label=r"Correlated ($c=0.3$)")
ax.set_xticks(d_vals)
ax.set_xticklabels([f"$d={d}$" for d in d_vals])
ax.set_xlabel("Code distance $d$")
ax.set_ylabel("Logical error rate $P_L$")
ax.set_title(r"Distance scaling at $p=0.400$")
ax.grid(True, linestyle="--", alpha=0.4)
ax.legend()
plt.tight_layout()
plt.savefig(outdir / "fig5_inversion_p400.png", dpi=300, bbox_inches="tight")
plt.close()
print("fig 5 done")

print("\nall figures saved to figures/")
