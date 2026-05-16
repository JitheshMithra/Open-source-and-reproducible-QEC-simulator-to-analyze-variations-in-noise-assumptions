**Current Version:** v1.2

**Status:** Reproducible baseline framework with analytical validation (ACTIVE)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19410366.svg)](https://doi.org/10.5281/zenodo.19410366) 
![License](https://img.shields.io/badge/License-MIT-green) 
![Python](https://img.shields.io/badge/Python-3.10-blue) 
![Field](https://img.shields.io/badge/Field-quant--ph-purple) 

<p align="center">
  <img width="702" height="197" alt="image" src="https://github.com/user-attachments/assets/940d5ead-d342-42a9-8c86-4aefb890cbd0" />
</p>



QECops is a lightweight and open-source simulation framework for studying how noise assumptions influence logical error behavior in quantum error correction (QEC). The question that it seeks to answer is: _How sensitive are QEC results to specific noise models?_

This tool simulates how physical noise models translate into logical error rates in quantum error correction codes. It combines numerical simulation with analytical validation to study error behavior in repetition codes, and sticks to its core value of reproducibility among trials.

### What it does
  - This tool applies independent bit flip noise at rate _p_ to a repetition code of length _n_
  - It decodes using majority-vote
  - It sweeps _p_ across a range which you can configure and estimates logical error rate with Monte Carlo estimation model
  - This tool computes the exact analytical baseline using the binomial distribution model
  - It will output a plot with both curves and an absolute error subplot (log-scaled) to validate the agreement

Most QEC tools are held to large institutions or require significant setup. This tool runs from the command line and has no other dependencies, being fully reproducible with fixed seed.

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
REPORTS/
  EXPERIMENTS.md
  reportlinks.md
requirements.txt
README.md
LICENSE
ExampleResult.png (view an example result here)
.gitignore
```
## Technical Reports:
A detailed explanation of the methodology, experiments, and results is available in the REPORTS folder.

**DOI**: [https://doi.org/10.5281/zenodo.19410366](https://doi.org/10.5281/zenodo.19410366)

**See /REPORTS for excecuted experiments and more on the stored technical reports related to this project and its methodology**

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
All simulations are executed from the src directory
```bash
#example simulation:
cd src
python -m QECops.plot --n 3 5 7 --trials 20000 --seed 0 --logicalbit 0 --pmin 0.0 --pmax 0.2 --pstep 0.02
```
Command Line arguments:
  - --n: one or more ODD repetition-code lengths
  - --trials: # of Monte Carlo trials per data point
  - --logicalbit: Logical bit to encode (0 or 1)
  - --pmin, --pmax, --pstep: Physical error rate sweep parameters
  - --seed: Random seed for reproducibility
### Output:
Each run will generate timestamped results directory which contains:
  - plot.png: Static plot containing LER vs_p_, Monte Carlo + Analytical results, absolute error subplot
  - plot.html: interactive Plotly html graph, useful for visualization and analytical curve
  - QECops_noise_results.txt: Table numerical values

These outputs provide a direct comparison between simulation and analytical predictions, including quantitative validation through the absolute error subplot.

The framework reproduces expected scaling behavior and highlights how different noise assumptions impact logical error suppression.

### Noise Model and Assumptions:
- Independent, uncorrelated bit flip noise (for now)
- Classical repetition code with majority vote decoding
- Analytical baseline from the binomial distribution

This tool takes a noise assumption as input and shows you what error correction looks like under that assumption, it does not implement new codes or claim theoretical thresholds. It does not model hardware directly.

### Limitations:
- Noise is assumed to be independent, no correlated or biased errors YET
- Only repetition codes and majority-vote decoding are currently supported
- Does not connect with real hardware
- Assumes product-state initialization (|0...0⟩)
- Uses a simplified bit-flip noise model
- Does not yet model full quantum states
- Not a circuit-level noise model

**Future additions (v2 plans):**
  - Extended noise models (biased, correlated)
  - Threshold analysis for SURE
  - Implement thresholds and acknowledge them in the system
  - Log scale plots for Monte Carlo (users can plot LER on a log scale, clear error suppression, effective distance scaling)
  - Seed control per-p (use rng across entire sweep, improves statistical independence and interpretation)
  - Parameter sweeps BEYOND p 
  - Optional saving results as JSON or CSV
  - Switch from absolute error analysis to relative error analysis for better comparing in low probability regions
  - Superposition states (e.g., Bell, GHZ)
  - Density matrixes
  - Circuit-level noise modeling
  - Exploring more advanced QEC codes using tools like [Stim](https://github.com/quantumlib/Stim)
  - Front-end UI, like a website
  - Python library? make it pip installable/Installable Package
  - analysis layer built on top of Qiskit / existing tools

### Acknowledgements:
- Special thanks to _Daniel Strano_, developer of [Qrack](https://github.com/unitaryfoundation/qrack), from the Unitary Foundation for external review and consistent feedback and mentoring on my methodology
- Special thanks to Dr. Haining Pan, University of Florida, for discussion about possible directions and feedback on assumptions
- Special thanks to Dr. Zebo Yang, Florida Atlantic University, for external review and feedback

### License:
If used or mentioned in published works please cite in the recommended format and reference this repository.

Copyright (c) [2025] [Jithesh Mithra].
It is licensed under the MIT License, available at [https://github.com/JitheshMithra/QECops].

**Contact**: 
- _Email_: jitheshmithra412 [at] gmail [dot] com
- _Linkedin_: https://www.linkedin.com/in/jitheshmithra/ 
