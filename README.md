**Current Version:** v2.5

**Status:** Comparative noise sensitivity analysis framework with bootstrap threshold estimation, extended to real quantum codes (ACTIVE)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19410365.svg)](https://doi.org/10.5281/zenodo.19410365)
[![License](https://img.shields.io/badge/License-MIT-green)](https://github.com/JitheshMithra/QECops/blob/main/LICENSE)
[![Qiskit](https://img.shields.io/badge/Qiskit-Aer-6929C4)](https://qiskit.github.io/qiskit-aer/)
![Field](https://img.shields.io/badge/Field-quant--ph-purple) 

[![Featured in The Quantum Insider](https://img.shields.io/badge/Featured-The%20Quantum%20Insider-00CED1)](https://thequantuminsider.com/2026/07/22/guest-post-what-a-high-school-student-found-when-he-stress-tested-a-quantum-benchmark/)

<p align="center">
   <img width="1024" height="291" alt="logos (3)" src="https://github.com/user-attachments/assets/f0e8e979-5e60-4ae4-9f42-d738679a8899" />
</p>



QECops is a lightweight, open-source Monte Carlo simulation framework for studying how noise assumptions influence logical error behavior in quantum error correction (QEC); Python. The question it seeks to answer is: How sensitive are QEC performance conclusions to the choice of noise model assumptions?

This tool simulates how physical noise models translate into logical error rates in repetition codes, compares behavior across four noise models with consistent methodology, and estimates pseudo-thresholds with bootstrap confidence intervals. It has since been extended to the Shor [[9,1,3]] code, simulating logical error rates under both uniform and spatially disordered noise with the same bootstrap methodology. Everything runs locally from the command line with no institutional access required.

### What it does

- Simulates repetition codes of distance d=3, 5, 7 (and beyond) under four noise models
- Simulates the Shor [[9,1,3]] code under uniform and spatially disordered noise
- Decodes using majority-vote (repetition code) or full syndrome-based correction (Shor code)
- Sweeps physical error rate p (or other noise parameters) across a configurable range
- Estimates logical error rate (LER) with Monte Carlo simulation
- Computes exact analytical baseline using the binomial distribution for bitflip noise
- Estimates pseudo-thresholds by finding where LER curves for adjacent distances cross
- Computes bootstrap confidence intervals on threshold and error rate estimates
- Exports results as PNG plots, interactive HTML, CSV, and JSON
- Validates simulation against analytical predictions with absolute error subplots

Most QEC tools need institutional access or require complex environments and setup. This tool runs locally from the command line and has no other dependencies, being fully reproducible with fixed seed.

### Noise Models:
**Bitflip**: independent per-qubit X errors at rate p. Analytical solution exists via binomial distribution. Used as baseline.

**Depolarizing**: symmetric X/Y/Z errors. Effective flip probability is 2p/3, accounting for the three error channels. Threshold is suppressed relative to bitflip.

**Biased**: asymmetric X and Z error rates via separate px and pz parameters. Z errors are invisible to the classical repetition code. Useful for studying hardware with asymmetric noise.

**Correlated**: spatially propagating errors. A flip at qubit i propagates to qubit i+1 with probability correlation. Models crosstalk and physically realistic error spreading.

**Disordered** (Shor code): each of the 9 physical qubits independently draws its own error rate from Uniform[p−δ, p+δ], where δ controls disorder strength while keeping the mean error rate equal to p. Models device-level gate fidelity variation across physical qubits.

**Current Project structure:**
```
src/
  QECops/
    __init__.py
    noise.py
    decode.py
    simulation.py
    plot.py
    analytical.py
  shor/
    shorcode.py
    verifyshor.py
    correctshor.py
    crosschecksyndrome.py
    montecarloShor.py
    montecarloShorinhomogeneous.py
    thresholdsweep.py
    disordersweep.py
    bootstrapdiff.py
    figuresdisorder.py
    figures/
    disorder_sweep_results.json
    disorder_bootstrap_diff_results.json
  analyzeS.py
  validation.py
  figures.py
  figures/
  results/
REPORTS/
  EXPERIMENTS.md
  reportlinks.md
requirements.txt
README.md
LICENSE
reproducibility.md
.gitignore
```
## Technical Reports:
A detailed explanation of the v1/v2 methodology, experiments, and results is available in the REPORTS folder.

**DOI**: [https://doi.org/10.5281/zenodo.19410365](https://doi.org/10.5281/zenodo.19410365) 

## Getting Started:
### Installation:
Make sure to download/update latest versions of pip, python, git, and related packages prior to running this simulation for best/optimal results. Once complete, proceed with installation instructions.

Clone the repository:
```bash
git clone https://github.com/JitheshMithra/QECops.git
cd QECops
```
Install required dependencies:
```bash
pip install -r requirements.txt
```
### Running the simulation:
All simulations are executed from the src directory. Paper-quality results use 100k trials.
```bash
cd src
```
**Basic bitflip threshold sweep:**
```bash
python -m QECops.plot --n 3 5 7 --trials 10000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise bitflip --showthresholds
```

**With bootstrap confidence intervals:**
```bash
python -m QECops.plot --n 3 5 7 --trials 10000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise bitflip --bootstrap --nbootstrap 1000 --showthresholds
```

**Validation mode (Monte Carlo vs analytical):**
```bash
python -m QECops.plot --n 3 5 7 --trials 20000 --seed 42 --pmin 0.01 --pmax 0.4 --pstep 0.02 --noise bitflip --plotmode validation
```

**Correlated noise sweep:**
```bash
python -m QECops.plot --n 3 5 7 --trials 10000 --seed 42 --pmin 0.05 --pmax 0.45 --pstep 0.05 --noise correlated --correlation 0.3 --bootstrap --nbootstrap 1000
```

**Depolarizing noise:**
```bash
python -m QECops.plot --n 3 5 7 --trials 10000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise depolarizing --bootstrap --nbootstrap 1000 --showthresholds
```
**Biased noise (sweeping X error rate, fixed Z):**
```bash
python -m QECops.plot --n 3 5 7 --trials 10000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise biased --sweepparam px --pz 0.05 --bootstrap --nbootstrap 1000 --showthresholds
```

### Running the Shor code study
Requires `qiskit` and `qiskit-aer` (installed via requirements.txt above). All commands run from `src/shor`.
```bash
cd shor
```
**Verify the implementation (encoding, syndrome detection, correction):**
```bash
python3 verifyshor.py
python3 correctshor.py
python3 crosschecksyndrome.py
```
**Uniform-noise threshold sweep with bootstrap CIs:**
```bash
python3 thresholdsweep.py
```
**Disorder comparison sweep (uniform vs. spatially inhomogeneous noise):**
```bash
python3 disordersweep.py
```
**Significance test on disorder vs. uniform (bootstrap-of-difference):**
```bash
python3 bootstrapdiff.py
```

### Command line arguments (repetition code: `QECops.plot`)
| Argument | Description | Default |
|---|---|---|
| `--n` | One or more odd code distances | required |
| `--trials` | Monte Carlo trials per data point | 10000 |
| `--seed` | Random seed for reproducibility | 0 |
| `--logicalbit` | Logical bit to encode, 0 or 1 | 0 |
| `--pmin` | Minimum physical error rate | 0.05 |
| `--pmax` | Maximum physical error rate | 0.4 |
| `--pstep` | Step size for error rate sweep | 0.05 |
| `--noise` | Noise model: bitflip, depolarizing, biased, correlated | bitflip |
| `--sweepparam` | Parameter to sweep: p, px, pz, correlation | p |
| `--px` | X error rate for biased noise | 0.05 |
| `--pz` | Z error rate for biased noise | 0.05 |
| `--correlation` | Correlation strength for correlated noise | 0.3 |
| `--fixedp` | Fixed p when sweeping correlation | 0.1 |
| `--plotmode` | validation or threshold | threshold |
| `--logscale` | Log scale y-axis | False |
| `--showthresholds` | Show threshold markers on plot | False |
| `--bootstrap` | Run bootstrap CI on threshold estimates | False |
| `--nbootstrap` | Number of bootstrap samples | 1000 |
| `--confidence` | Confidence level for CI | 0.95 |
| `--export` | Export format: txt, csv, json, all | all |

### Command line arguments (Shor code: `src/shor`)
| Script | Arguments | Notes |
|---|---|---|
| `thresholdsweep.py` | `--trials --seed --pmin --pmax --pstep --nbootstrap --bootstrapseed --out` | defaults: trials=3000, seed=42, pmin=0.05, pmax=0.35, pstep=0.05 |
| `disordersweep.py` | `--trials --realizations --seed --p (nargs) --delta (nargs) --nbootstrap --bootstrapseed --out` | defaults: p=[0.10, 0.15, 0.20, 0.25], delta=[0.0, 0.05, 0.10] |
| `bootstrapdiff.py` | `--in --out` | reads disordersweep.py output |

Trial counts are much lower by default than the repetition code tool (~3000 vs. 10000+) due to the per-trial cost of simulating actual Qiskit circuits; override with `--trials` for larger runs.

Example:
```bash
python thresholdsweep.py --trials 5000 --pmin 0.05 --pmax 0.40 --pstep 0.05 --seed 7
python disordersweep.py --trials 4000 --realizations 40 --p 0.10 0.20 0.30 --delta 0.0 0.10 0.20
```

`verifyshor.py`, `correctshor.py`, `crosschecksyndrome.py`, `montecarloShor.py`, and `montecarloShorinhomogeneous.py` are unparameterized correctness/smoke-test scripts, run directly with no flags.

## Limitations
- Phenomenological noise only: no circuit-level gate or measurement noise
- Repetition code: majority-vote decoding only, no surface codes or other stabilizer codes, classical simulation with no quantum state representation
- Shor code: single fixed code distance (no natural second distance for a pseudo-threshold crossing); disorder model is a simple i.i.d. per-qubit uniform draw, not spatially correlated disorder; lower Monte Carlo trial counts (~3000-4000 vs. 10000+) due to per-trial circuit simulation cost
- No hardware integration

## Future Work
- Circuit-level noise modeling
- Cyclic QEC with multiple syndrome rounds
- Surface code support
- Extended disorder strength range and code concatenation for a genuine second Shor code distance
- [Stim](https://github.com/quantumlib/stim) comparison layer for cross-validation
- pip installable package
- Relative Error Subplots

### Acknowledgements:
- Special thanks to _Dr. Haining Pan_, University of Florida, for discussion about possible directions and feedback on assumptions. 
- Special thanks to _Daniel Strano_, developer of [qrack](https://github.com/unitaryfoundation/qrack), from the Unitary Foundation for external review and consistent feedback and mentoring on my methodology
- Special thanks to _Dr. Zebo Yang_, Florida Atlantic University, for external review and feedback

### License:
If used or mentioned in published works please cite in the recommended format and reference this repository.

Copyright (c) [2025] [Jithesh Mithra].
It is licensed under the MIT License, available at [https://github.com/JitheshMithra/QECops].

**Contact**: 
- _Email_: jitheshmithra412 [at] gmail [dot] com
- _Linkedin_: https://www.linkedin.com/in/jitheshmithra/ 
