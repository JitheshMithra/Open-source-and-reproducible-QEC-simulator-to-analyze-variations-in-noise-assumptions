import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import plotly.graph_objects as go

from .analytical import analyticalLogicalError
from .simulation import distanceSweep, estimatePseudoThreshold, bootstrapThreshold, thresholdScalingSummary, robustnessMetric, failureBoundary, crossingConsistency


def argParser():
    parser = argparse.ArgumentParser(description="QECops experiment runner")

    parser.add_argument("--n", nargs="+", type=int, required=True)
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--logicalbit", type=int, choices=[0, 1], default=0)

    parser.add_argument("--pmin", type=float, default=0.05)
    parser.add_argument("--pmax", type=float, default=0.4)
    parser.add_argument("--pstep", type=float, default=0.05)

    parser.add_argument("--noise", choices=["bitflip", "depolarizing", "biased", "correlated"], default="bitflip")
    parser.add_argument("--sweepparam", default="p")

    parser.add_argument("--px", type=float, default=0.05)
    parser.add_argument("--pz", type=float, default=0.05)
    parser.add_argument("--correlation", type=float, default=0.3)
    parser.add_argument("--fixedp", type=float, default=0.1)

    parser.add_argument("--plotmode", choices=["validation", "threshold"], default="threshold")
    parser.add_argument("--logscale", action="store_true")
    parser.add_argument("--showthresholds", action="store_true")
    parser.add_argument("--export", choices=["txt", "csv", "json", "all"], default="all")
    parser.add_argument("--bootstrap", action="store_true")
    parser.add_argument("--nbootstrap", type=int, default=1000)
    parser.add_argument("--confidence", type=float, default=0.95)

    parser.add_argument("--mode", choices=["phenomenological", "circuit"], default="phenomenological")
    parser.add_argument("--cyclic", action="store_true")
    parser.add_argument("--rounds", type=int, default=1)

    return parser.parse_args()


def pValuesGen(pmin, pmax, pstep):
    if pstep <= 0:
        raise ValueError("pstep must be positive")
    if pmin < 0 or pmax > 1 or pmin > pmax:
        raise ValueError("pmin and pmax must be between 0 and 1, with pmin <= pmax")

    pValues = []
    p = pmin
    while p <= pmax + 1e-12:
        pValues.append(round(p, 12))
        p += pstep
    return pValues


def validateDistances(nValues):
    for n in nValues:
        if n <= 0:
            raise ValueError("n values must be positive")
        if n % 2 == 0:
            raise ValueError("n values must be odd")


def validateArgs(args):
    if args.mode == "circuit":
        raise NotImplementedError("Circuit-level noise is not implemented yet.")

    if args.cyclic:
        raise NotImplementedError("Cyclic QEC is not implemented yet.")

    if args.plotmode == "validation" and args.noise not in ["bitflip", "depolarizing"]:
        raise ValueError("Validation mode only supports bitflip or depolarizing noise.")

    if args.plotmode == "validation" and args.sweepparam != "p":
        raise ValueError("Validation mode requires --sweepparam p.")

    if args.noise == "biased" and args.sweepparam == "p":
        raise ValueError("For biased noise, use --sweepparam px or --sweepparam pz.")

    if args.noise == "correlated" and args.sweepparam not in ["p", "correlation"]:
        raise ValueError("For correlated noise, use --sweepparam p or --sweepparam correlation.")

    validateDistances(args.n)


def buildNoiseParams(args):
    if args.noise == "biased":
        return {"px": args.px, "pz": args.pz}

    if args.noise == "correlated":
        if args.sweepparam == "correlation":
            return {"p": args.fixedp, "correlation": args.correlation}
        return {"correlation": args.correlation}

    return {}


def printReadableSummary(results, thresholds, ci, summary, args):
    print("\nQECops Simulation Summary")
    print("--------------------------")

    for d, curve in sorted(results.items()):
        print(f"\nCode distance d={d}")
        print(f"{args.sweepparam}-value    Logical Error Rate    Uncertainty      Failures")

        for r in curve:
            print(
                f"{r['sweep_value']:.3f}      "
                f"{r['LER']:.6f}              "
                f"±{r['stderr']:.6f}      "
                f"{r['failures']}/{r['trials']}"
            )

    print("\nPseudo-threshold estimates")
    print("--------------------------")

    if not thresholds:
        print("Not computed because sweep parameter is not p.")
    else:
        for pair, value in thresholds.items():
            if value is None:
                print(f"d={pair[0]} vs d={pair[1]}: no crossing found")
            else:
                print(f"d={pair[0]} vs d={pair[1]}: p ~= {value:.4f}")

    if ci:
        print(f"\nPseudo-threshold confidence intervals ({int(args.confidence*100)}%)")
        print("--------------------------------------------")
        for pair, val in ci.items():
            if val is None:
                print(f"d={pair[0]} vs d={pair[1]}: no crossing found")
            else:
                print(f"d={pair[0]} vs d={pair[1]}: {val['mean']:.4f} [{val['lower']:.4f}, {val['upper']:.4f}] ±{val['std']:.4f}")

    print("\nScaling behavior")
    print("----------------")
    for row in summary:
        status = "improves with distance" if row["error_suppression_with_distance"] else "does not improve cleanly"
        print(f"p={row['physical_error_rate']:.3f}: {status}")


def exportResults(resultsDir, results, thresholds, ci, summary, exportType, args):
    rows = []

    for d, curve in sorted(results.items()):
        for r in curve:
            rows.append({
                "distance": d,
                "physical_error_rate": r["physical_error_rate"],
                "sweep_param": r["sweep_param"],
                "sweep_value": r["sweep_value"],
                "LER": r["LER"],
                "stderr": r["stderr"],
                "failures": r["failures"],
                "trials": r["trials"],
                "noise_type": r["noise_type"],
            })

    if exportType in ["txt", "all"]:
        with open(resultsDir / "summary.txt", "w") as f:
            f.write("QECops Simulation Summary\n")
            f.write("--------------------------\n\n")

            for d, curve in sorted(results.items()):
                f.write(f"Code distance d={d}\n")
                f.write(f"{args.sweepparam}-value    Logical Error Rate    Uncertainty      Failures\n")

                for r in curve:
                    f.write(
                        f"{r['sweep_value']:.3f}      "
                        f"{r['LER']:.6f}              "
                        f"±{r['stderr']:.6f}      "
                        f"{r['failures']}/{r['trials']}\n"
                    )
                f.write("\n")

            f.write("Pseudo-threshold estimates\n")
            f.write("--------------------------\n")

            if not thresholds:
                f.write("Not computed because sweep parameter is not p.\n")
            else:
                for pair, value in thresholds.items():
                    if value is None:
                        f.write(f"d={pair[0]} vs d={pair[1]}: no crossing found\n")
                    else:
                        f.write(f"d={pair[0]} vs d={pair[1]}: p ~= {value:.4f}\n")

            if ci:
                f.write(f"\nConfidence intervals ({int(args.confidence*100)}%)\n")
                f.write("--------------------------\n")
                for pair, val in ci.items():
                    if val is None:
                        f.write(f"d={pair[0]} vs d={pair[1]}: no crossing found\n")
                    else:
                        f.write(f"d={pair[0]} vs d={pair[1]}: {val['mean']:.4f} [{val['lower']:.4f}, {val['upper']:.4f}]\n")

    if exportType in ["csv", "all"]:
        with open(resultsDir / "raw_results.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    if exportType in ["json", "all"]:
        ciSerializable = {}
        if ci:
            for k, v in ci.items():
                ciSerializable[str(k)] = v

        with open(resultsDir / "raw_results.json", "w") as f:
            json.dump(
                {
                    "results": rows,
                    "thresholds": {str(k): v for k, v in thresholds.items()},
                    "confidence_intervals": ciSerializable,
                    "scaling_summary": summary,
                },
                f,
                indent=4,
            )


def makeThresholdPlot(resultsDir, results, thresholds, ci, args):
    fig, ax = plt.subplots(figsize=(8, 6))

    for d, curve in sorted(results.items()):
        x = [r["sweep_value"] for r in curve]
        y = [max(r["LER"], 1 / r["trials"]) for r in curve]
        err = [max(r["stderr"], 1 / r["trials"]) for r in curve]
        ax.errorbar(x, y, yerr=err, marker="o", capsize=3, label=f"d={d}")

    if args.showthresholds and thresholds:
        for pair, threshold in thresholds.items():
            if threshold is not None:
                ax.axvline(threshold, linestyle=":", label=f"threshold {pair} ~= {threshold:.3f}")

    if ci:
        for pair, val in ci.items():
            if val is not None:
                ax.axvspan(val["lower"], val["upper"], alpha=0.1,
                          label=f"{int(args.confidence*100)}% CI {pair}")

    ax.set_xlabel(args.sweepparam)
    ax.set_ylabel("Logical error rate")
    ax.set_title(f"Threshold plot: LER vs {args.sweepparam}, noise={args.noise}")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()

    if args.logscale:
        ax.set_yscale("log")

    plt.tight_layout()
    plt.savefig(resultsDir / "threshold_plot.png", dpi=300, bbox_inches="tight")
    plt.close()


def makeValidationPlot(resultsDir, results, args):
    if args.noise not in ["bitflip", "depolarizing"]:
        raise ValueError("Validation plot only supports bitflip/depolarizing.")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)

    for d, curve in sorted(results.items()):
        x = [r["sweep_value"] for r in curve]
        ysimRaw = [r["LER"] for r in curve]
        ysimPlot = [max(r["LER"], 1 / r["trials"]) for r in curve]
        err = [max(r["stderr"], 1 / r["trials"]) for r in curve]

        yanalytic = [analyticalLogicalError(d, p) for p in x]
        yabsError = [
            max(abs(s - a), 1 / curve[0]["trials"])
            for s, a in zip(ysimRaw, yanalytic)
        ]

        ax1.errorbar(x, ysimPlot, yerr=err, marker="o", capsize=3, label=f"d={d} Monte Carlo")
        ax1.plot(x, yanalytic, linestyle="--", label=f"d={d} Analytical")
        ax2.plot(x, yabsError, marker="x", label=f"d={d}")

    ax1.set_ylabel("Logical error rate")
    ax1.set_title(f"Validation plot, noise={args.noise}")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend()

    ax2.set_xlabel(args.sweepparam)
    ax2.set_ylabel("|Simulation - Analytical|")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", linestyle="--", alpha=0.5)
    ax2.legend(title="Distance")

    if args.logscale:
        ax1.set_yscale("log")

    plt.tight_layout()
    plt.savefig(resultsDir / "validation_plot.png", dpi=300, bbox_inches="tight")
    plt.close()


def makeInteractivePlot(resultsDir, results, args):
    fig = go.Figure()

    for d, curve in sorted(results.items()):
        x = [r["sweep_value"] for r in curve]
        y = [max(r["LER"], 1 / r["trials"]) for r in curve]

        fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=f"d={d} Monte Carlo"))

    fig.update_layout(
        title=f"QECops LER vs {args.sweepparam}, noise={args.noise}",
        xaxis_title=args.sweepparam,
        yaxis_title="Logical error rate",
        hovermode="x unified",
    )

    if args.logscale:
        fig.update_yaxes(type="log")

    fig.write_html(resultsDir / "interactive_plot.html")


def plotRun(args):
    validateArgs(args)

    pValues = pValuesGen(args.pmin, args.pmax, args.pstep)

    runId = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    resultsDir = Path.cwd() / "results" / runId
    resultsDir.mkdir(parents=True, exist_ok=True)

    noiseParams = buildNoiseParams(args)

    results = distanceSweep(
        distances=args.n,
        pValues=pValues,
        trials=args.trials,
        seed=args.seed,
        logicalBit=args.logicalbit,
        noiseType=args.noise,
        sweepParam=args.sweepparam,
        **noiseParams,
    )

    thresholds = estimatePseudoThreshold(results)

    if args.sweepparam != "p":
        print("Warning: pseudo-thresholds are only standard when sweeping p.")
        thresholds = {}

    ci = None
    if args.bootstrap and thresholds:
        ci = bootstrapThreshold(results, nBootstrap=args.nbootstrap, confidence=args.confidence)

    summary = thresholdScalingSummary(results)

    printReadableSummary(results, thresholds, ci, summary, args)
    exportResults(resultsDir, results, thresholds, ci, summary, args.export, args)

    if args.plotmode == "validation":
        makeValidationPlot(resultsDir, results, args)
    else:
        makeThresholdPlot(resultsDir, results, thresholds, ci, args)

    makeInteractivePlot(resultsDir, results, args)
    boundary = failureBoundary(summary)
    consistency = crossingConsistency(thresholds, ci)
    if boundary is not None:
        print(f"\nFailure boundary detected at p = {boundary:.3f}")
    else:
        print("\nNo failure boundary detected in swept range")

    print("\nCrossing consistency")
    print("--------------------")
    for pair, result in consistency.items():
        status = "consistent" if result["consistent"] else "inconsistent"
        print(f"d={pair[0]} vs d={pair[1]}: {status} ({result['reason']})")
    print("\nSaved files:")
    print(resultsDir)

def main():
    args = argParser()
    plotRun(args)


if __name__ == "__main__":
    main()
