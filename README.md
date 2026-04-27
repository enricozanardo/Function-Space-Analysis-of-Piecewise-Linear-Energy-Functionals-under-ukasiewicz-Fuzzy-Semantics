# Numerical validation — supplementary code

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/enricozanardo/Function-Space-Analysis-of-Piecewise-Linear-Energy-Functionals-under-ukasiewicz-Fuzzy-Semantics)](../../releases)

Companion code for Section 8 of the manuscript

> Zanardo, E. & Ragusa, M. A. (2026). *Function Space Analysis of
> Piecewise-Linear Energy Functionals under Łukasiewicz Fuzzy Semantics*.

## Quick start

```bash
git clone https://github.com/enricozanardo/Function-Space-Analysis-of-Piecewise-Linear-Energy-Functionals-under-ukasiewicz-Fuzzy-Semantics.git
cd Function-Space-Analysis-of-Piecewise-Linear-Energy-Functionals-under-ukasiewicz-Fuzzy-Semantics
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
of the script. The canonical reference output is committed to the repository
as `results/reference_results.json`; comparing the live `results/results.json`
against it (e.g. with `diff` or `jq`) verifies bit-for-bit reproducibility.

## File layout

```
experiments/
├── CITATION.cff                 (citation metadata, CFF v1.2.0)
├── LICENSE                      (MIT)
├── README.md                    (this file)
├── numerical_validation.py      (driver script)
└── results/
    ├── reference_results.json   (canonical output, committed)
    └── results.json             (live output, gitignored)
```

## Releases & DOI

This repository follows [semantic versioning](https://semver.org).

- **v1.0.2** (2026-04-27) — first Zenodo-archived release.  Metadata-only
  bump from v1.0.1 to fix `CITATION.cff` parsing (the previous file
  contained a placeholder DOI in its `preferred-citation` block, which
  caused Zenodo to reject the deposit with *Citation metadata load failed*).
  Code is byte-identical to v1.0.0/v1.0.1.
- **v1.0.1** (2026-04-27) — metadata bump after enabling the GitHub–Zenodo
  integration; not archived because of the malformed `CITATION.cff` noted above.
- **v1.0.0** (2026-04-27) — first public release accompanying the manuscript;
  not Zenodo-archived because the integration was enabled afterwards.

Tagged releases are archived on Zenodo through the
[GitHub–Zenodo integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content).
Once the deposit has been minted, the DOI badge below resolves to the
archived snapshot:

```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PLACEHOLDER.svg)](https://doi.org/10.5281/zenodo.PLACEHOLDER)
```

## Citation

If you use this software, please cite both the software and the manuscript:

```bibtex
@software{zanardo2026numericalvalidation,
  author    = {Zanardo, Enrico and Ragusa, Maria Alessandra},
  title     = {Numerical validation for ``Function Space Analysis of
               Piecewise-Linear Energy Functionals under {\L}ukasiewicz
               Fuzzy Semantics''},
  year      = {2026},
  version   = {v1.0.2},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.PLACEHOLDER},
  url       = {https://github.com/enricozanardo/Function-Space-Analysis-of-Piecewise-Linear-Energy-Functionals-under-ukasiewicz-Fuzzy-Semantics}
}
```

The accompanying manuscript should be cited from the same `references.bib`
distributed with the LaTeX source of the paper.

## License

MIT — see [LICENSE](LICENSE).
