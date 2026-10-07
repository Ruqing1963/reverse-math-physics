# How Much Infinity Does Physics Need?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23212848.svg)](https://doi.org/10.5281/zenodo.23212848)

Code, data, figures and manuscript for

> **R. Chen**, *How Much Infinity Does Physics Need? Reverse Mathematics of Physical Limits and the Empirical
> Boundary* (2026). DOI: [10.5281/zenodo.23212848](https://doi.org/10.5281/zenodo.23212848)

## Summary

- **Four equivalences over RCA₀.** Each is equivalent to arithmetical comprehension ACA₀:
  - Fekete's lemma for subadditive sequences (existence of thermodynamic limits);
  - the existence of the infrared limit of every bounded monotone c-function;
  - the existence of the scaling exponent of every Cantor construction with nondecreasing dimension parameters;
  - the existence of the set of indices whose computable fractal limit sets are connected.
- **Individual quantities.** A limit with a computable modulus that RCA₀ verifies is provably existent in RCA₀; a
  limit that is not a computable real is not, because the ω-model REC contains no real equal to it.
- **Moduli in the series.** The Ising free energy on an L×L torus has error 0.639913/L² at K_c, matching the Ising
  CFT torus amplitude ln[(θ₂+θ₃+θ₄)/2η](i) = 0.639912; the diamond-lattice critical-cluster mass converges geometrically
  with ratio Λ₂/Λ₁ = 0.3553; the embezzlement asymptotics has error O(1/(n ln n)). The Ω approximations have no
  computable modulus.
- **The empirical boundary.** The limits that entered explicit certification budgets are provable in RCA₀; those shown
  to be empirically inaccessible (Ω-dimensions, limit connectivity, finite-dimensional game values) require ACA₀. Since
  ACA₀ is conservative over Peano arithmetic, the infinity physical limits require is real but logically mild.

## Repository structure

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and compiled PDF of the manuscript |
| `code/reverse_math_moduli.py` | Script producing every number and the figure in the paper |
| `figures/` | Figure as vector PDF (used by the paper) and PNG |
| `data/` | Numerical data as CSV (metadata in `#` header lines) |
| `results/` | Console output of the script (the numbers quoted in the paper) |

## Reproducing the results

Requirements: Python ≥ 3.10 and the packages in `requirements.txt`
(tested with Python 3.12.4, NumPy 1.26.4, SciPy 1.13.1, Matplotlib 3.8.4). Runtime is about 15 seconds.

```bash
pip install -r requirements.txt
python code/reverse_math_moduli.py      # add --show to display the figure
```

Run from the repository root; the figure goes to `figures/` and data to `data/` (override with `RM_FIG_DIR`,
`RM_DATA_DIR`).

| Data file | Content | Paper |
|---|---|---|
| `ising_free_energy_error.csv` | free-energy error on the L×L torus at K_c and 0.8 K_c | §4, Fig. 1 |
| `dhl_mass_convergence.csv` | convergence of the normalised critical-cluster mass | §4, Fig. 1 |
| `embezzlement_error.csv` | error of the embezzlement asymptotics | §4, Fig. 1 |
| `omega_certified_interval.csv` | certified lower and upper bounds for Ω of binary lambda calculus | §4, Fig. 1 |

To rebuild the paper (pdfLaTeX, two passes):

```bash
cd paper
pdflatex Chen_2026_Reverse_Math_Physics.tex
pdflatex Chen_2026_Reverse_Math_Physics.tex
```

## Citation

```bibtex
@misc{Chen2026ReverseMathPhysics,
  author = {Chen, Ruqing},
  title  = {How Much Infinity Does Physics Need? Reverse Mathematics of Physical Limits and the Empirical Boundary},
  year   = {2026},
  doi    = {10.5281/zenodo.23212848},
  url    = {https://doi.org/10.5281/zenodo.23212848}
}
```

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Manuscript, figures, data and results** (`paper/`, `figures/`, `data/`, `results/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## Contact

Ruqing Chen — GUT Geoservice Inc., Montreal — ruqing@hotmail.com
