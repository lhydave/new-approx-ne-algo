"""Discover the candidate bound and active comparisons in the manuscript's Sigma_r.

The coordinates are (r, rho, eta, xi, alpha, gamma). All constraints
are the relaxed scalar bounds and four segment predicates in Appendix C.
No observed equality or algebraic root is imposed on the optimizer.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize


def segment_minimum(a, b):
    """Minimum of the maximum of the two affine endpoint interpolants."""
    R_0, C_0 = a
    R_1, C_1 = b
    values = [max(a), max(b)]
    D = (R_1 - R_0) - (C_1 - C_0)
    if D != 0:
        t = (C_0 - R_0) / D
        if 0 <= t <= 1:
            values.append((1 - t) * R_0 + t * R_1)
    return min(values)


def evaluate(point):
    r, rho, eta, xi, alpha, gamma = point
    sigma = 1 - rho
    u, v = (1 - eta) / eta, (1 - xi) / xi
    s_1 = 1 - r - r / rho + u * (1 - r - r / sigma)
    s_2 = 1 - r - r / sigma + v * (1 - r - r / rho)
    affine_h_1 = 1 - r * (1 + sigma) / rho + sigma * (1 - r) * v / rho
    affine_h_2 = 1 - r * (1 + rho) / sigma + rho * (1 - r) * u / sigma
    h_1, h_2 = min(1, affine_h_1), min(1, affine_h_2)
    S, T = h_1 - r + u * (s_2 - r), h_2 - r + v * (s_1 - r)
    U, V = alpha * s_1 + (1 - alpha) * h_2, (1 - gamma) * h_1 + gamma * s_2
    H, K = alpha * S + (1 - alpha) * T, (1 - gamma) * S + gamma * T
    E_1 = (1 - gamma) * (H - r) + gamma * v * (U - r) + alpha * gamma * (1 - r - r * v) - alpha * r
    E_2 = (1 - alpha) * (K - r) + alpha * u * (V - r) + alpha * gamma * (1 - r - r * u) - gamma * r
    D_c, D_r = rho * (1 - eta) + sigma * eta, sigma * (1 - xi) + rho * xi
    p_dagger, q_dagger = rho * (1 - eta) / D_c, sigma * (1 - xi) / D_r
    ell_c, ell_r = r - eta * (1 - r / rho), r - xi * (1 - r / sigma)

    def regret_pair(p, q):
        j_c = ell_c * min((1 - p) * rho / (eta * sigma), p / (1 - eta))
        j_r = ell_r * min((1 - q) * sigma / (xi * rho), q / (1 - xi))
        R_r = (1 - p) * (1 - q) * r + q * (1 - p * r / rho) - j_r
        C_r = (1 - p) * (1 - q) * r + p * (1 - q * r / sigma) - j_c
        return R_r, C_r

    center = regret_pair(p_dagger, q_dagger)
    endpoints = [(p_dagger, 0), (p_dagger, 1), (0, q_dagger), (1, q_dagger)]
    residuals = {}
    for name, weight in zip(("rho", "eta", "xi", "alpha", "gamma"), point[1:]):
        residuals[name + "_lower"] = weight - r
        residuals[name + "_upper"] = 1 - r - weight
    for name, value in zip(("s_1", "s_2", "h_1", "h_2", "S", "T"), (s_1, s_2, h_1, h_2, S, T)):
        residuals[name + "_nonnegative"] = value
    for name, value in zip(("U", "V", "H", "K"), (U, V, H, K)):
        residuals[name + "_minus_r"] = value - r
    residuals.update(E_1=E_1, E_2=E_2)
    for j, endpoint in enumerate(endpoints, 1):
        residuals["segment_" + str(j)] = segment_minimum(center, regret_pair(*endpoint)) - r
    values = dict(s_1=s_1, s_2=s_2, h_1=h_1, h_2=h_2, S=S, T=T,
                  affine_h_1=affine_h_1, affine_h_2=affine_h_2,
                  p_dagger=p_dagger, q_dagger=q_dagger)
    return residuals, values


def run(starts=32, seed=20261004):
    rng = np.random.default_rng(seed)
    records = []
    for index in range(starts):
        initial = np.r_[0.305, rng.uniform(0.32, 0.68, 5)]
        result = minimize(lambda point: -point[0], initial, method="SLSQP",
                          bounds=[(0.3, 1 / 3)] + [(0.3, 0.7)] * 5,
                          constraints=[{"type": "ineq", "fun": lambda point: np.array(list(evaluate(point)[0].values()))}],
                          options={"ftol": 1e-12, "maxiter": 2000})
        residuals, values = evaluate(result.x)
        records.append({"start": index, "initial": initial.tolist(),
                        "point": result.x.tolist(), "solver_success": bool(result.success),
                        "message": str(result.message), "iterations": int(result.nit),
                        "minimum_residual": float(min(residuals.values())),
                        "residuals": {name: float(value) for name, value in residuals.items()},
                        "values": {name: float(value) for name, value in values.items()}})
    feasible = [item for item in records if item["minimum_residual"] >= -1e-8]
    if not feasible:
        raise RuntimeError("No numerically feasible point returned.")
    best = max(feasible, key=lambda item: item["point"][0])
    # Match the paper's S > r > T orientation while retaining every raw run.
    point = best["point"]
    exchanged = best["values"]["S"] < best["values"]["T"]
    if exchanged:
        r, rho, eta, xi, alpha, gamma = point
        point = [r, 1 - rho, xi, eta, gamma, alpha]
    residuals, values = evaluate(point)
    paper_candidate = {"players_exchanged": exchanged, "point": point,
                       "minimum_residual": float(min(residuals.values())),
                       "residuals": {name: float(value) for name, value in residuals.items()},
                       "values": {name: float(value) for name, value in values.items()}}
    return {"scope": "Floating-point local exploration of Sigma_r; no exact or global optimality claim.",
            "optimizer": "scipy.optimize.minimize(method='SLSQP')",
            "scipy_version": scipy.__version__, "seed": seed, "starts": starts,
            "coordinates": ["r", "rho", "eta", "xi", "alpha", "gamma"],
            "ftol": 1e-12, "maxiter": 2000, "feasibility_tolerance": 1e-8,
            "feasible_runs": len(feasible), "best": best,
            "paper_candidate": paper_candidate, "runs": records}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = run(args.starts, args.seed)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"starts": report["starts"], "feasible_runs": report["feasible_runs"],
                      "best_r": report["best"]["point"][0], "out": str(args.out)}))
