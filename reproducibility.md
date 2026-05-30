# QECops v2 — Reproducibility Guide

All results in the paper are reproducible using QECops v2.0.
Repository: https://github.com/JitheshMithra/QECops

## Setup
```bash
cd QECops
pip install -r requirements.txt
```

## Primary Results
All commands are run from `QECops/src/`.

### Bitflip baseline

```bash
python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise bitflip --bootstrap --nbootstrap 1000 --showthresholds
```

### Depolarizing

```bash
python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.60 --pmax 0.85 --pstep 0.05 --noise depolarizing --bootstrap --nbootstrap 1000 --showthresholds
```

Note: depolarizing threshold falls near nominal p=0.75, corresponding to effective bitflip probability p_eff = 2p/3 = 0.50. The sweep range is extended to p=0.85 to capture the crossing cleanly.

### Biased noise (px sweep, pz=0.05 fixed)

```bash
python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise biased --sweepparam px --pz 0.05 --bootstrap --nbootstrap 1000 --showthresholds
```

Note: biased noise with fixed pz and swept px is equivalent to independent bitflip noise in px. This run serves as an internal consistency check.

### Correlated noise — correlation strength sweep

```bash
python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise correlated --correlation 0.1 --bootstrap --nbootstrap 1000 --showthresholds

python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise correlated --correlation 0.2 --bootstrap --nbootstrap 1000 --showthresholds

python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise correlated --correlation 0.3 --bootstrap --nbootstrap 1000 --showthresholds

python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise correlated --correlation 0.4 --bootstrap --nbootstrap 1000 --showthresholds

python -m QECops.plot --n 3 5 7 --trials 100000 --seed 42 --pmin 0.05 --pmax 0.55 --pstep 0.05 --noise correlated --correlation 0.5 --bootstrap --nbootstrap 1000 --showthresholds
```

## Sensitivity Metric S

After running the correlated noise sweep above, update threshold values in `src/analyzeS.py` and run:

```bash
cd src
python analyzeS.py
```

## Qiskit Validation

```bash
cd src
python validation.py
```

Validates QECops Monte Carlo against Qiskit Aer simulation and analytical bitflip prediction across d=3,5,7 and p=0.05 to 0.50. Maximum absolute error between QECops and Qiskit is 0.003, within statistical uncertainty for N=100,000 trials.

## Parameter Summary

| Parameter | Value |
|---|---|
| Code distances | d = 3, 5, 7 |
| Trials per data point | N = 100,000 |
| Base random seed | 42 |
| Seed scheme | Independent per (p, d, noise model) via pt_seed = seed + idx * 99991 + n * 7 |
| Physical error rate sweep | p = 0.05 to 0.55, step 0.05 (0.60 to 0.85 for depolarizing) |
| Bootstrap iterations | B = 1,000 |
| Bootstrap confidence level | 95% |
| Correlated noise strengths | c = 0.1, 0.2, 0.3, 0.4, 0.5 |
| Biased noise fixed pz | 0.05 |