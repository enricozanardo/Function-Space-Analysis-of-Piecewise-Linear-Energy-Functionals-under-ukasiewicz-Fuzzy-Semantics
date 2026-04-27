# Numerical validation — supplementary code

Companion code for Section 8 of *Function Space Analysis of Piecewise-Linear
Energy Functionals under Łukasiewicz Fuzzy Semantics*.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy
python numerical_validation.py --output-dir results
```

Wall-clock time: ≈ 60 seconds on a single CPU core.

## What is reproduced

| Manuscript section | Quantity | Reference value (n = 6, m = 10) |
|---|---|---|
| §8.2 Partition function | `log_Z` empirical (importance sampling) | `3.521 ± 0.002` |
| §8.2 | Bound `−W ≤ log Z ≤ W` (W ≈ 4.72) | satisfied |
| §8.3 Mixing time | Empirical (reflected Langevin, 120 chains) | `44 steps × h = 0.22 t.u.` |
| §8.3 | Theoretical upper bound (Cor. 7.5) | `≈ 3.85 × 10³` |
| §8.4 Laplace asymptotic | `log_Z(β·ŵ)` at β ∈ {1, 4, 16} | `0.73, 2.97, 12.51` |
| §8.4 | Leading order `β·E_max` | `1, 4, 16` |
| §8.5 Scaling | `mean W` for n ∈ {2, 4, 8, 16, 32, 64} | `1.0, 2.9, 8.6, 22.7, 55.4, 134.3` |

All numbers above are seeded for bit-exact reproducibility under
`numpy >= 1.24`, `scipy >= 1.10`, with the seed `20260426` declared at the top
of the script.

## File layout

```
experiments/
├── README.md                    (this file)
├── numerical_validation.py      (driver script)
└── results/                     (auto-created on first run)
    └── results.json             (machine-readable outputs)
```

## Citation

If you reuse this code, please cite the parent manuscript and the LIMEN-AI
companion papers listed in `paper6/references.bib`.
