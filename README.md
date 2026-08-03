**Current Version:** v2.5

**Status:** Shor [[9,1,3]] code robustness to spatially inhomogeneous noise disorder & repetition-code noise sensitivity analysis with bootstrap threshold estimation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19410365.svg)](https://doi.org/10.5281/zenodo.19410365) ![License](https://img.shields.io/badge/License-MIT-green) 
![Python](https://img.shields.io/badge/Python-3.10-blue) 
![Field](https://img.shields.io/badge/Field-quant--ph-purple) 
![Qiskit](https://img.shields.io/badge/Qiskit-Aer-6929C4)

<p align="center">
   <img width="734" height="236" alt="image" src="https://github.com/user-attachments/assets/50658cb1-b35d-450b-af73-4a9f9f6be833" />
</p>



QECops is a lightweight, open-source Monte Carlo simulation framework for studying how noise assumptions influence logical error behavior in quantum error correction (QEC); Python. The question it seeks to answer is: _How sensitive are QEC performance conclusions to the choice of noise model assumptions?_

This tool simulates how physical noise models translate into logical error rates in repetition codes, compares behavior across four noise models with consistent methodology, and estimates pseudo-thresholds with bootstrap confidence intervals. Everything runs locally from the command line with no institutional access required.

### What it does
   - Simulates repetition codes of distance d=3, 5, 7 (and beyond) under four noise models
   - Decodes using majority-vote
   - Sweeps physical error rate p (or other noise parameters) across a configurable range
   - Estimates logical error rate (LER) with Monte Carlo simulation
   - Computes exact analytical baseline using the binomial distribution for bitflip noise
   - Estimates pseudo-thresholds by finding where LER curves for adjacent distances cross
   - Computes bootstrap confidence intervals on threshold estimates
   - Exports results as PNG plots, interactive HTML, CSV, and JSON
   - Validates simulation against analytical predictions with absolute error subplots

Most QEC tools need instituitional access or require complex environments and setup. This tool runs locally from the command line and has no other dependencies, being fully reproducible with fixed seed.

### Noise Models:
**Bitflip**: independent per-qubit X errors at rate p. Analytical solution exists via binomial distribution. Used as baseline.

**Depolarizing**: symmetric X/Y/Z errors. Effective flip probability is 2p/3, accounting for the three error channels. Threshold is suppressed relative to bitflip.

**Biased**: asymmetric X and Z error rates via separate px and pz parameters. Z errors are invisible to the classical repetition code. Useful for studying hardware with asymmetric noise.

**Correlated**: spatially propagating errors. A flip at qubit i propagates to qubit i+1 with probability correlation. Models crosstalk and physically realistic error spreading.

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
A detailed explanation of the v1 methodology, experiments, and results is available in the REPORTS folder.

Simulation results show strong agreement with analytical predictions across code distances n=3 to n=11.

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
### Command line arguments

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

### Output

Each run generates a timestamped results directory containing:

   - `threshold_plot.png`: LER vs sweep parameter for all distances with error bars
   - `validation_plot.png`: Monte Carlo vs analytical with absolute error subplot (validation mode)
   - `interactive_plot.html`: interactive Plotly visualization
   - `summary.txt`: readable numerical summary with threshold estimates and CI
   - `raw_results.csv`: raw results table
   - `raw_results.json`: full results including threshold estimates and scaling summary

## Limitations

- Phenomenological noise only: no circuit-level gate or measurement noise
- Repetition code only: no surface codes or other stabilizer codes
- Majority vote decoding only: no minimum weight perfect matching
- Classical simulation: no quantum state representation
- No hardware integration

## Future Work

- Circuit-level noise modeling
- Cyclic QEC with multiple syndrome rounds
- Surface code support
- Density Matrices/Quantum State simulation/Superposition states
- [Qiskit](https://github.com/Qiskit/qiskit)/[Stim](https://github.com/quantumlib/stim) comparison layer for cross-validation
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
