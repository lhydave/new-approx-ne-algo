# Exact arithmetic verification

This directory contains the arithmetic checkers and fixed certificates accompanying *A Simple Algorithm Breaking the 1/3 Barrier for Approximate Nash Equilibria in Bimatrix Games*. The paper supplies the strategy semantics and mathematical reductions; the registered checks verify the arithmetic obligations at the places identified in Appendices A–F.

The table below maps every unit to its paper location and script. Dependencies and precise scopes are recorded in the authoritative registry, [`bound/bundle/proof-manifest.json`](bound/bundle/proof-manifest.json).

## Setup and commands

From the repository root, create a Python environment and install the pinned dependency:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Then run:

```sh
# Inspect the 29 registered units.
python proof/bound/bundle/verify.py --list

# Check all 63 input hashes, sizes, and existing evidence predicates.
python proof/bound/bundle/verify.py --check-metadata

# Freshly execute all registered arithmetic checks.
python proof/bound/bundle/verify.py --replay
```

Metadata validation checks fixed input hashes and any available generated evidence. It lists absent reports; a fresh clone creates these with `--replay`. Metadata validation does not establish a fresh arithmetic replay. Full `--replay` reconstructs algebraic identities, polynomial signs, root enclosures, and subdivision coverage. It uses the final fixed certificates, including the canonical tree's `--deep` mode, and does not generate new trees or reuse historical prefix acceptance reports.

The runner resolves all file paths relative to its own directory, so it also works from another current directory. All required scripts, helper modules, polynomial data, and fixed certificate inputs are in this repository. Verification reports are regenerated locally. Python 3.9.6 with SymPy 1.14.0 is the tested exact-verification environment. Do not use optimized Python (`-O` or `PYTHONOPTIMIZE`): assertions are part of verification.

## Paper correspondence

| Paper location | Arithmetic obligation | Registered units |
| --- | --- | ---: |
| Appendix A | Typed, same-game instantiation and ideal payoff provenance | 2 |
| Appendix B | Signed Jensen-gap bridge and finite OptimalMixing comparisons | 2 |
| Appendix C | Necessary projection from the stationary regret to the comparison threshold | 3 |
| Appendix D | Exact local root isolation and scalar equality fixture | 2 |
| Appendix E.1–E.2 | Polynomial identities, source bounds, orientation, and shared-mass reductions | 7 |
| Appendix E.3 | Preclip continuation, stopping events, and postclip directions | 11 |
| Appendix E.4 | Full-range switch-face resultant and global exclusion | 1 |
| Appendix F | Transfer of finite stationary-point and mixing errors | 1 |
| **Total** | | **29** |

The three subdivision certificates independently verify **5,098**, **19,771**, and **4,945** terminal leaves, respectively, with full partition coverage and no open leaves. The checker scopes distinguish canonical and large stationary-weight regions, preclip and postclip branches, and the constraints preserved along each reduction.

The paper spells Greek symbols out in script identifiers (`rho`, `eta`, `xi`, `alpha`, `gamma`), with subscripts represented by underscores (`h_2`, `N_R`, `B_mix`, `N_1`).

Symbols that the paper writes with accents or Greek capitals use the following identifiers:

| Paper symbol | Script identifier | Meaning |
| --- | --- | --- |
| $\alpha_0$ | `alpha_0` | balancing value of $\alpha$, $(r-T)/(S-T)$ |
| $\hat g$ | `g_hat` | conservative stationary value $g-\delta_f$ |
| $\Theta$ | `Theta` | the bound $\alpha_0\le1-r$, multiplied by $m+n$ |
| $\Phi_1,\Phi_2$ | `Phi_1`, `Phi_2` | additive constants of the stationarity LP |
| $\Phi(\cdot)$ | `Phi(...)` | stationarity payoff expression, instantiated at the listed strategies |

In scripts, a trailing `_bar` denotes a complement, e.g. `alpha_bar` is $1-\alpha$ and `xi_bar` is $1-\xi$. The Stage-I contract text in `bound/quantified/stage-i-contracts.json` is preserved verbatim and keeps its original notation (for example `A_s`, `B_s` for $\Phi_1,\Phi_2$).

Some inherited directory and helper filenames describe an earlier construction. The registered algorithm is the final **four-LP + OptimalMixing + five-profile selection** algorithm; Appendices A–C derive the retained arithmetic from the first four LP calls and the OptimalMixing comparisons.

## Registered checkers

| Unit | Appendix | Script relative to `proof/` |
| --- | --- | --- |
| `quantified-interface-identities` | A | [replay_quantified_identities.py](bound/quantified/replay_quantified_identities.py) |
| `ideal-cap-provenance-ledger` | A | [replay_ideal_cap_provenance_ledger.py](bound/global/replay_ideal_cap_provenance_ledger.py) |
| `full-BR-rectangle-bridge` | B | [replay_fullbr_rectangle_bridge.py](bound/quantified/replay_fullbr_rectangle_bridge.py) |
| `parametric-OM-linearizer` | B | [replay_parametric_optimalmix.py](bound/quantified/replay_parametric_optimalmix.py) |
| `near-BR-delta-ledger` | F | [replay_delta_f_lift_ledger.py](bound/critical/replay_delta_f_lift_ledger.py) |
| `g-identities` | C | [cross_projection_symbolic_checks.py](analytic_scalar_bound/cross_projection_symbolic_checks.py) |
| `g-derivative-small` | C | [replay_cross_derivative_certificate.py](analytic_scalar_bound/replay_cross_derivative_certificate.py) |
| `g-derivative-large` | C | [replay_cross_derivative_large_mass_certificate.py](analytic_scalar_bound/replay_cross_derivative_large_mass_certificate.py) |
| `compact-identities` | E.1–E.2 | [verify_compact_identities.py](bound/simplify/verify_compact_identities.py) |
| `mass-identities` | E.1–E.2 | [check_reduction_identities.py](bound/global/check_reduction_identities.py) |
| `mass-M-small` | E.1–E.2 | [replay_m_chart_certificate.py](bound/global/replay_m_chart_certificate.py) |
| `mass-M` | E.1–E.2 | [replay_m_chart_full_mass_certificate.py](bound/global/replay_m_chart_full_mass_certificate.py) |
| `mass-h_2-xi` | E.1–E.2 | [replay-h_2-xi-certificates.py](bound/simplify/replay-h_2-xi-certificates.py) |
| `canonical-u-cap` | E.1–E.2 | [replay-relaxed-source-knot-u-cap.py](bound/simplify/replay-relaxed-source-knot-u-cap.py) |
| `full-mass-both-high` | E.1–E.2 | [replay-fullmass-both-high.py](bound/simplify/replay-fullmass-both-high.py) |
| `fixed-v-and-u1` | E.3 | [replay_preclip_left_projection_signs.py](bound/global/replay_preclip_left_projection_signs.py) |
| `canonical-alignment` | E.3 | [preclip_canonical_alignment.py](bound/preclip_canonical_alignment.py) |
| `canonical-alpha_0-stop` | E.3 | [replay-alpha_0-mass-stop-exclusion.py](bound/simplify/replay-alpha_0-mass-stop-exclusion.py) |
| `canonical-E1` | E.3 | [replay-arc-e1-source-budget.py](bound/simplify/replay-arc-e1-source-budget.py) |
| `canonical-boundary-stops` | E.3 | [replay-preclip-stop-events.py](bound/simplify/replay-preclip-stop-events.py) |
| `canonical-upper-arc` | E.3 | [replay_preclip_bernstein_interval.py](bound/kink/replay_preclip_bernstein_interval.py) |
| `high-rho-left` | E.3 | [replay_preclip_left_interval.py](bound/global/replay_preclip_left_interval.py) |
| `high-rho-below` | E.3 | [replay_preclip_below_interval.py](bound/global/replay_preclip_below_interval.py) |
| `postclip-envelope` | E.3 | [replay_postclip_crossing_certificate.py](bound/global/replay_postclip_crossing_certificate.py) |
| `postclip-orientation` | E.3 | [replay_postclip_orientation_certificate.py](bound/global/replay_postclip_orientation_certificate.py) |
| `postclip-knots` | E.3 | [replay_postclip_knot_direction.py](bound/global/replay_postclip_knot_direction.py) |
| `natural-root` | D | [local_root_certificate.py](bound/critical/local_root_certificate.py) |
| `full-rho-kink-face` | E.4 | [replay_full_rho_kink_face.py](bound/critical/replay_full_rho_kink_face.py) |
| `scalar-root-fixture` | D | [local_signs.py](bound/critical/local_signs.py) |

The exact Stage-I contract text read by the interface checker is preserved in [`bound/quantified/stage-i-contracts.json`](bound/quantified/stage-i-contracts.json). The companion registry specifies the typed building-block clauses and composition. Historical derivation notes are omitted; the paper contains the mathematical argument.

## Check the approximately 0.30954 root

To run the exact local root check:

```sh
python proof/bound/bundle/verify.py --replay --unit natural-root
```

It checks the paper's system `N_1 = N_2 = H_root = 0` in its rational box of radius `10^-30`. A rational Krawczyk enclosure lies strictly inside the box, and a rational interval contraction bound proves local uniqueness. Numerical values used to construct a preconditioner are rationalized before the certificate inequalities are tested.

To check the switch-face arithmetic together with its root dependency:

```sh
python proof/bound/bundle/verify.py --replay --unit full-rho-kink-face
```

This isolated command verifies that face's stated scope. The full global conclusion also needs the preceding reductions and all their certificates; use the full replay for those. The bound-discovery computation is [`../tools/find_bound.py`](../tools/find_bound.py). It supplies the numerical candidate and active comparisons used to derive Appendix D's algebraic equations; the exact root checks and global exclusion then certify the bound. Parameters, output coordinates, and reproduction instructions are in the [discovery README](../tools/README.md).

## Inputs, reports, and scope

The manifest registers SHA-256 hashes and sizes for 63 proof inputs. Its `requires` fields specify arithmetic dependencies; the complete mathematical composition is described by its `logical_chain`, the paper, and this README. `--unit ID` includes the selected unit's transitive registered dependencies; it does not constitute a full proof replay.

Runtime summaries and regenerated verification reports are ignored by Git. Fixed certificate data and all 63 registered inputs are retained, including `critical/local-root-certificate.json`, which is also read as an input by later checks. A clone without generated reports can validate input metadata and run the complete replay; metadata output explicitly lists reports that must first be generated. Do not update a hash merely to silence a mismatch: review the changed mathematics and certificate first.

The checkers use exact rational identities, interval bounds, Bernstein coefficients, and complete fixed subdivision certificates. Standard LP duality, building-block semantics, and sound instantiation/forgetting are justified by the paper's mathematical proof. This is not a proof-assistant kernel formalization. The scalar equality fixture does not establish realization by an actual game or worst-case tightness of the algorithm.
