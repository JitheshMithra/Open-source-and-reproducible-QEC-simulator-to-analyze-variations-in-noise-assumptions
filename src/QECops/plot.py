import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import plotly.graph_objects as go

from .analytical import analytical_logical_error
from .simulation import distancesweep, estimatepseudothreshold, bootstrapthreshold, thresholdscalingsummary


def argparser():
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


def pvaluesgen(pmin, pmax, pstep):
    if pstep <= 0:
        raise ValueError("pstep must be positive")
    if pmin < 0 or pmax > 1 or pmin > pmax:
        raise ValueError("pmin and pmax must be between 0 and 1, with pmin <= pmax")

    pvalues = []
    p = pmin
    while p <= pmax + 1e-12:
        pvalues.append(round(p, 12))
        p += pstep
    return pvalues


def validate_distances(nvalues):
    for n in nvalues:
        if n <= 0:
            raise ValueError("n values must be positive")
        if n % 2 == 0:
            raise ValueError("n values must be odd")


def validate_args(args):
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

    validate_distances(args.n)


def build_noise_params(args):
    if args.noise == "biased":
        return {"px": args.px, "pz": args.pz}

    if args.noise == "correlated":
        if args.sweepparam == "correlation":
            return {"p": args.fixedp, "correlation": args.correlation}
        return {"correlation": args.correlation}

    return {}


def print_readable_summary(results, thresholds, ci, summary, args):
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


def export_results(resultsdir, results, thresholds, ci, summary, export_type, args):
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

    if export_type in ["txt", "all"]:
        with open(resultsdir / "summary.txt", "w") as f:
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

    if export_type in ["csv", "all"]:
        with open(resultsdir / "raw_results.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    if export_type in ["json", "all"]:
        ci_serializable = {}
        if ci:
            for k, v in ci.items():
                ci_serializable[str(k)] = v

        with open(resultsdir / "raw_results.json", "w") as f:
            json.dump(
                {
                    "results": rows,
                    "thresholds": {str(k): v for k, v in thresholds.items()},
                    "confidence_intervals": ci_serializable,
                    "scaling_summary": summary,
                },
                f,
                indent=4,
            )


def make_threshold_plot(resultsdir, results, thresholds, ci, args):
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
    plt.savefig(resultsdir / "threshold_plot.png", dpi=300, bbox_inches="tight")
    plt.close()


def make_validation_plot(resultsdir, results, args):
    if args.noise not in ["bitflip", "depolarizing"]:
        raise ValueError("Validation plot only supports bitflip/depolarizing.")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)

    for d, curve in sorted(results.items()):
        x = [r["sweep_value"] for r in curve]
        ysim_raw = [r["LER"] for r in curve]
        ysim_plot = [max(r["LER"], 1 / r["trials"]) for r in curve]
        err = [max(r["stderr"], 1 / r["trials"]) for r in curve]

        yanalytic = [analytical_logical_error(d, p) for p in x]
        yabs_error = [
            max(abs(s - a), 1 / curve[0]["trials"])
            for s, a in zip(ysim_raw, yanalytic)
        ]

        ax1.errorbar(x, ysim_plot, yerr=err, marker="o", capsize=3, label=f"d={d} Monte Carlo")
        ax1.plot(x, yanalytic, linestyle="--", label=f"d={d} Analytical")
        ax2.plot(x, yabs_error, marker="x", label=f"d={d}")

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
    plt.savefig(resultsdir / "validation_plot.png", dpi=300, bbox_inches="tight")
    plt.close()


def make_interactive_plot(resultsdir, results, args):
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

    fig.write_html(resultsdir / "interactive_plot.html")


def plotrun(args):
    validate_args(args)

    pvalues = pvaluesgen(args.pmin, args.pmax, args.pstep)

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    resultsdir = Path.cwd() / "results" / run_id
    resultsdir.mkdir(parents=True, exist_ok=True)

    noise_params = build_noise_params(args)

    results = distancesweep(
        distances=args.n,
        pvalues=pvalues,
        trials=args.trials,
        seed=args.seed,
        logicalbit=args.logicalbit,
        noisetype=args.noise,
        sweepparam=args.sweepparam,
        **noise_params,
    )

    thresholds = estimatepseudothreshold(results)

    if args.sweepparam != "p":
        print("Warning: pseudo-thresholds are only standard when sweeping p.")
        thresholds = {}

    ci = None
    if args.bootstrap and thresholds:
        ci = bootstrapthreshold(results, nbootstrap=args.nbootstrap, confidence=args.confidence)

    summary = thresholdscalingsummary(results)

    print_readable_summary(results, thresholds, ci, summary, args)
    export_results(resultsdir, results, thresholds, ci, summary, args.export, args)

    if args.plotmode == "validation":
        make_validation_plot(resultsdir, results, args)
    else:
        make_threshold_plot(resultsdir, results, thresholds, ci, args)

    make_interactive_plot(resultsdir, results, args)

    print("\nSaved files:")
    print(resultsdir)


def main():
    args = argparser()
    plotrun(args)


if __name__ == "__main__":
    main()