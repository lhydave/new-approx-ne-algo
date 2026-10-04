# A Simple Algorithm Breaking the 1/3 Barrier for Approximate Nash Equilibria in Bimatrix Games

**Hanyu Li and Dongchen Li**

This repository contains the [paper](main.pdf), its [LaTeX source](main.tex), and the exact arithmetic checks supporting its proof. The algorithm achieves an approximation guarantee of approximately **0.30954 + δ**, breaking the 1/3 barrier, in polynomial time. It is deterministic and symmetric, with no tunable hyperparameters.

## Algorithm and proof architecture

The algorithm enriches the Tsaknakis–Spirakis framework with strategy-mixing building blocks. It obtains a stationary profile and dual strategies, solves four fixed-side SimpleMixing linear programs in two rounds, runs OptimalMixing on the original stationary/dual rectangle, and selects the best of five resulting profiles.

The analysis first instantiates the building-block properties at strategies from the same game. It then forgets payoff information in a sound, necessary direction to obtain a finite scalar failure system. Global reductions and exact arithmetic certificates exclude every failure threshold above the isolated root

```text
r* ≈ 0.30953996543144345
```

The final step restores the stationary-point and OptimalMixing error budgets. The repository provides proof verification software; the algorithm itself is specified in the paper.

The research extends the LegoNE building-block methodology. OpenAI Codex played a central role in discovering the algorithm and developing its proof, including the specialized exact verification tools; the paper discusses the broader role of LLM agents in this process.

| Path | Contents |
| --- | --- |
| `main.tex`, `main.pdf` | Manuscript and exported paper; the bibliography is embedded in the source |
| [`proof/README.md`](proof/README.md) | Exact verification setup, proof architecture, and scope |
| [`proof/bound/bundle/proof-manifest.json`](proof/bound/bundle/proof-manifest.json) | Authoritative unit registry, dependencies, scopes, and input hashes |
| [`tools/find_bound.py`](tools/find_bound.py) | Floating-point candidate discovery described in Section 5.3 and Appendix D |
| `requirements.txt`, `requirements-discovery.txt` | Pinned exact-verification and optional numerical dependencies |

## Verify the proof

From the repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python proof/bound/bundle/verify.py --check-metadata
python proof/bound/bundle/verify.py --replay
```

The first check validates all 63 registered input hashes and any locally available evidence reports. Missing generated reports are listed and can be created by `--replay`. The second freshly runs all 29 arithmetic verification units, including complete subdivision coverage. It uses only files in this repository and does not require the original Desktop workspace or a numerical search. Full replay can take several minutes. Use ordinary Python execution without `-O` or `PYTHONOPTIMIZE`, because the checkers use assertions.

The exact checks use Python's rational arithmetic and SymPy. They have been exercised with Python 3.9.6 and SymPy 1.14.0. See the [proof README](proof/README.md) for individual-unit commands and the distinction between the mathematical proof and its machine-checked arithmetic.

The [proof README](proof/README.md#registered-checkers) contains the correspondence for all 29 registered units.

## Reproduce the bound discovery

The numerical computation finds the candidate bound near 0.30954 and the active comparisons that lead to the algebraic equations in Section 5.3 and Appendix D. Exact root isolation and global exclusion then certify the bound.

```sh
python -m pip install -r requirements-discovery.txt
python tools/find_bound.py --starts 32 --seed 20261004 --out results/find-bound.json
```

The defaults reproduce the paper's 32 SLSQP starts, seed, iteration limit, and tolerances. The report saves every initial and returned point, solver status, and constraint residual, together with the candidate in the player orientation used by the paper. It is a generated experiment output and is ignored by Git. See [`tools/README.md`](tools/README.md).

## Build the paper

With a TeX distribution containing the packages used by `main.tex`, run from the repository root:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The source is standalone, includes its bibliography and TikZ drawings, and needs no external figure files. Experimental outputs, regenerated verification reports, and LaTeX intermediates are ignored; fixed proof inputs and PDF exports remain available for version control.

## License

The repository is released under the [MIT License](LICENSE).
