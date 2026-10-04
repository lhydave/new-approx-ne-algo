# Bound discovery

[`find_bound.py`](find_bound.py) reproduces the computation used to find the bound and its active comparisons in Section 5.3 and Appendix D of the paper. It directly evaluates the necessary scalar system `Sigma_r` in the coordinates `(r, rho, eta, xi, alpha, gamma)`, including the six payoff bounds, four weighted comparisons, two cross-call comparisons, and four segment conditions. Segment minima are evaluated at the endpoints and any intersection of the two affine regret bounds.

From the repository root, in the Python environment used for verification:

```sh
python -m pip install -r requirements-discovery.txt
python tools/find_bound.py --starts 32 --seed 20261004 --out results/find-bound.json
```

The pinned numerical environment is NumPy 2.0.2 and SciPy 1.13.1, tested with Python 3.9.6. Defaults match the paper: SLSQP minimizes `-r` with 32 starts, initial `r=0.305`, independent initial weights sampled uniformly from `[0.32, 0.68]`, seed `20261004`, at most 2,000 iterations per start, and `ftol=1e-12`. A returned point is counted as numerically feasible when every constraint residual is at least `-1e-8`. The report records all returned points, solver status, residuals, derived values, and the best numerically feasible value. Use `--starts`, `--seed`, or `--out` to change these options.

The report also includes a `paper_candidate` obtained by exchanging the players when needed to match the paper's orientation `S > r > T`; the original optimizer output is preserved in `best` and `runs`.

The 32-start calculation returns 22 numerically feasible runs and a best threshold of approximately **0.3095399654277576** in the pinned environment. In the paper's player orientation, the approximate returned coordinates are:

| Coordinate | Numerical candidate |
| --- | ---: |
| `r` | 0.3095399654 |
| `rho` | 0.4799488827 |
| `eta` | 0.5461144556 |
| `xi` | 0.4328331390 |
| `alpha` | 0.4572048381 |
| `gamma` | 0.5427951619 |
The `paper_candidate` provides the active cap, `H = K = r`, an interior segment comparison, and `E_1 = 0` used to derive the algebraic equations in Appendix D. The optimization imposes these as inequalities and discovers their near-equalities through the returned residuals. The [proof bundle](../proof/README.md) then isolates the exact root and establishes the global bound.

The local search uses floating-point arithmetic; feasibility is assessed with the stated residual tolerance. The complete generated report is written under `results/`, which is ignored by Git. The directory is created automatically. Fixed exact proof inputs remain versioned.
