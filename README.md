**Current Version:** v1.0

**Status:** Baseline code/implementation
# Open source and reproducible QEC simulator to analyze variations in noise assumptions
## Overview:
This project is meant to be an open‑source, reproducible, research-grade, and reusable simulation tool made to analyze and simulate how noise assumptions affect the performance of quantum error correction.
Quantum error correction depends strongly on assumed noise models, yet many existing simulations are held to specific experiments or papers, and are not open to the public and far-reaching educated use. This tool gives a simple and clear framework for measuring and simulating how error correction performance changes as noise changes and varies.
Simply stated: This tool aims to analyze error correction performance under explicit noise assumptions than measuring noise directly from hardware coded with the python language.

**Why does it matter:**
Error correction suppresses noise only under certain conditions, so small changes in noise assumptions can lead to very large differences in outputs/conclusions.
This tool makes those assumptions inspectable, transparent and reproducible, thus allowing users to:
  - test the strength of results
  - compare code sizes
  - understand failure behaviors
  - reproduce and resimulate baseline behaviors easily

### This tool aims to:
  - Simulate bit-flip noise
  - Will apply a repetition code
  - Will decode using majority-vote decoding
  - Measures logical error as a function of physical error rate

Primary output is essentially a standard diagnostic plot used in many quantum error correction research - being logical error rate vs physical error rate, but for different code sizes.

**Noise Model:**
This project models independent bit-flip noise, where each bit is flipped with a probability known as _p_. This represents assumed hardware error rate, which is supplied as an input to the simulation.

**Error correction:**
Logical information is protected with a reptition code, which encodes one logical bit into multiple physical bits. Decoding is done with majority-vote decoding, which selects the most likely logical value after the noise is applied.

**Output:**
For a given physical error rate _p_ and the size of the code _n_, the tool will estimate the logical error rate, which is the probability that decoding will fail, and it uses the Monte Carlo simulation by repeatedly sampling random noise.
So by producing _p_ across a range of values and comparing the different code sizes, the tool will make a standard diagnostic plot that shows logical error rate vs physical error rate. The plot quantifies the data and shows how effectively the redundancy suppresses the errors and identifies where error correction will succeed or fail.

### Limitations:
  - The tool does not measure hardware noise, it analyzes the error correction performance of a given assumed noise parameters
  - Noise is assumed to be independent and uncorrelated
  - The current version mode focuses on repetition codes and majority-vote decoding

The goal of this project is to give a simple, transparent, and reproducible simulation framework open to the public for studying how noise assumptions can change and impact error correction and its results.

**Assumptions:**
Currently:
  - Independent bit flip noise
  - no correlated errors
  - classical repetition code
  - majority vote decoding

This project does not attempt to create or show new error-correction codes or to prove theoretical thresholds and boundaries. It is made and intended to be used as a measurement and analysis tool, **NOT** a theory checker or unfound breakthrough, just an easy tool for all to use.

**Project status:**
v1 has been successfully pushed! v2 planning/coding has begun.
Future extensions/ideas are seen down below under v2 plans.

## Getting Started:
### Installation:
Clone the repository:
```bash
git clone https://github.com/JitheshMithra/Open-source-reproducible-simulator-to-analyze-noise-assumptions-in-quantum-error-correction.git
cd Open-source-reproducible-simulator-to-analyze-noise-assumptions-in-quantum-error-correction
```
Install required dependencies:
```bash
pip install -r requirements.txt
```
Intended Project structure (v1):
```
src/
  QECsim/
    __init__.py
    noise.py
    decode.py
    simulation.py
    plot.py
requirements.txt
README.md
LICENSE
.gitignore
```
### Running the simulation:
All simulations are executed from the src directory
```bash
#example simulation:
cd src
python -m QECsim.plot --n 3 5 7 --trials 5000 --seed 0 --logicalbit 0 --pmin 0.0 --pmax 0.2 --pstep 0.02
```
Command Line arguments:
  - --n: one or more ODD repetition-code lengths
  - --trials: # of Monte Carlo trials per data point
  - --logicalbit: Logical bit to encode (0 or 1)
  - --pmin, --pmax, --pstep: Physical error rate sweep parameters
  - --seed: Random seed for reproducibility
### Output:
Each run will generate timestamped results directory which contains:
  - plot.png: Static plot image of logical error rate vs physical error rate
  - plot.html: interactive Plotly html graph, useful for visualization
  - QECsim_noise_results.txt: Table numerical values
These values/outputs will show direct comparison of logical error reduction across different code sizes under specified noise assumptions

### Reproducibility:
  - All simulations are fully reproducible using the fixed random seeds specified in the CLI
  - Noise models and decoding assumptions are expected to be explicitly and specifically defined
  - Results can be regenerated or extended by changing simulation parameters

Future additions (v2 plans):
  - Extended noise models (biased, correlated)
  - Compare analytic vs simulated results (overlaying analytic curves and explicit convergence)
  - Log scale plots (users can plot LER on a log scale, clear error suppression, effective distance scaling)
  - Seed control per-p (use rng across entire sweep, improves statistical independence and interpretation)
  - Parameter sweeps BEYOND p 
  - Optional saving results as JSON or CSV

### License:
If used or mentioned in published works please cite in the recommended format.

Copyright (c) [2025] [Jithesh Mithra].
It is licensed under the MIT License, available at [https://github.com/JitheshMithra/Open-source-reproducible-simulator-to-analyze-noise-assumptions-in-quantum-error-correction/].

Also check License page in files.
