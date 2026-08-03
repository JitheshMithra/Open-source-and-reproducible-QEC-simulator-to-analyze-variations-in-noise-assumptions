# figuresdisorder.py
# run from src/shor/
# python figuresdisorder.py

import json
from pathlib import Path

import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.labelsize": 12, "legend.fontsize": 10, "figure.dpi": 300})
outDir = Path("figures")
outDir.mkdir(exist_ok=True)

with open("disorder_sweep_results.json") as f:
    results = json.load(f)

byP = {}
for r in results:
    byP.setdefault(r["p"], []).append(r)

colors = {0.10: "#1f77b4", 0.15: "#ff7f0e", 0.20: "#2ca02c", 0.25: "#d62728"}
markers = {0.10: "o", 0.15: "s", 0.20: "^", 0.25: "D"}

# fig 1: LER vs disorder strength, one line per physical error rate
fig, ax = plt.subplots(figsize=(6, 4.5))
for p, rows in sorted(byP.items()):
    rows.sort(key=lambda r: r["disorder_strength"])
    x = [r["disorder_strength"] for r in rows]
    y = [r["ler"] for r in rows]
    yerrLow = [r["ler"] - r["ci_low"] for r in rows]
    yerrHigh = [r["ci_high"] - r["ler"] for r in rows]
    ax.errorbar(x, y, yerr=[yerrLow, yerrHigh], marker=markers[p], color=colors[p],
                capsize=3, linewidth=1.5, markersize=6, label=f"$p={p}$")
ax.set_xlabel(r"Disorder strength $\delta$")
ax.set_ylabel("Logical error rate")
ax.set_title("Shor code LER vs. spatial disorder strength")
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()
fig.savefig(outDir / "disorder_ler_vs_delta.png")
plt.close(fig)

# fig 2: bar chart, grouped by p, one bar per delta, with CI whiskers -- easier to eyeball CI overlap
deltas = sorted({r["disorder_strength"] for r in results})
ps = sorted(byP.keys())
width = 0.25
fig, ax = plt.subplots(figsize=(7, 4.5))
for i, delta in enumerate(deltas):
    y = []
    yerrLow = []
    yerrHigh = []
    for p in ps:
        row = next(r for r in byP[p] if r["disorder_strength"] == delta)
        y.append(row["ler"])
        yerrLow.append(row["ler"] - row["ci_low"])
        yerrHigh.append(row["ci_high"] - row["ler"])
    xPos = [j + (i - 1) * width for j in range(len(ps))]
    ax.bar(xPos, y, width=width, yerr=[yerrLow, yerrHigh], capsize=3, label=rf"$\delta={delta}$")
ax.set_xticks(range(len(ps)))
ax.set_xticklabels([f"$p={p}$" for p in ps])
ax.set_ylabel("Logical error rate")
ax.set_title("Shor code LER by disorder strength (error bars = 95% bootstrap CI)")
ax.legend()
fig.tight_layout()
fig.savefig(outDir / "disorder_ler_bars.png")
plt.close(fig)

print(f"Saved figures to {outDir.resolve()}")
