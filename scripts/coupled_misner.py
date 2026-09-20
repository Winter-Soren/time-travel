"""Coupled, homogeneous Einstein--massless-scalar benchmark (8*pi*G=c=1).

Run: .venv/Scripts/python.exe coupled_misner.py
This tests compact Bianchi I near its flat Misner limit, not Ori's curved core.
"""
from project_paths import project_path

import json
from pathlib import Path
import sys

import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import solve_ivp


FIELDS = ("A", "B", "H_A", "H_B", "Pi", "phi")


def initial_state(epsilon):
    """Use a rationalized root to preserve the constraint for tiny epsilon."""
    hb = -(epsilon * epsilon / 2) / (1 + np.sqrt(1 + epsilon * epsilon / 2))
    return np.array([1.0, 1.0, -1.0, hb, epsilon, 0.0])


def proper_time_rhs(_, y):
    a, b, ha, hb, pi, _ = y
    theta = ha + 2 * hb
    return np.array([ha * a, hb * b, -theta * ha, -theta * hb,
                     -theta * pi, pi])


def exact_state(q, y0):
    theta0 = y0[2] + 2 * y0[3]
    return np.array([q ** (y0[2] / theta0), q ** (y0[3] / theta0),
                     y0[2] / q, y0[3] / q, y0[4] / q,
                     (y0[4] / theta0) * np.log(q)])


def run(epsilon, rtol=1e-12, method="DOP853", q_end=1e-6):
    if not np.isfinite(epsilon) or epsilon < 0 or not 0 < q_end < 1:
        raise ValueError("Require finite epsilon >= 0 and 0 < q_end < 1")
    y0 = initial_state(epsilon)
    s = -(y0[2] + 2 * y0[3])
    lambdas = np.linspace(0, -np.log(q_end), 301)

    def log_time_rhs(lam, y):
        # q only prescribes the time coordinate; all six states are evolved.
        # d(tau)/d(lambda)=exp(-lambda)/s; no exact state is inserted in RHS.
        return np.exp(-lam) / s * proper_time_rhs(None, y)

    atol = rtol * 1e-3 * np.array([1, 1, 1, max(abs(y0[3]), 1e-16),
                                 max(epsilon, 1e-16), max(epsilon / s, 1e-16)])
    sol = solve_ivp(log_time_rhs, (0, lambdas[-1]), y0, method=method,
                    t_eval=lambdas, rtol=rtol, atol=atol)
    if not sol.success or sol.y.shape != (6, len(lambdas)):
        raise RuntimeError(f"{method}, epsilon={epsilon}: {sol.message}")
    q = np.exp(-lambdas)
    reference = exact_state(q, y0)
    errors = {}
    for index, field in enumerate(FIELDS):
        # Relative errors for nonzero fields; phi uses a small initial scale
        # only to avoid division by zero at its exactly zero initial value.
        floor = (epsilon / s) * 1e-6 if field == "phi" else 0.0
        scale = np.maximum(np.abs(reference[index]), floor)
        scale = np.where(scale == 0, 1, scale)
        errors[field] = float(np.max(np.abs(sol.y[index] - reference[index]) / scale))
    _, _, ha, hb, pi, _ = sol.y
    terms = np.array([2 * ha * hb, hb * hb, -pi * pi / 2])
    constraint = terms.sum(axis=0)
    denominator = np.abs(terms).sum(axis=0)
    normalized_constraint = np.divide(np.abs(constraint), denominator,
                                      out=np.zeros_like(constraint), where=denominator > 0)
    endpoint = sol.y[:, -1]
    return {
        "epsilon": epsilon, "method": method, "rtol": rtol,
        "q_end": q_end, "sample_count": len(lambdas), "nfev": sol.nfev,
        "solver_success": bool(sol.success),
        "tau_end": float((1 - q_end) / s), "q_zero_tau_exact": float(1 / s),
        "endpoint": dict(zip(FIELDS, map(float, endpoint))),
        "exact_endpoint": dict(zip(FIELDS, map(float, reference[:, -1]))),
        "max_component_scaled_error": errors,
        "max_normalized_constraint_residual": float(np.max(normalized_constraint)),
        "max_absolute_constraint_residual": float(np.max(np.abs(constraint))),
        "ricci_scalar_endpoint": float(-endpoint[4] ** 2),
        "ricci_tensor_squared_endpoint": float(endpoint[4] ** 4),
        "exact_ricci_scalar_endpoint": float(-(epsilon / q_end) ** 2),
    }


def main():
    amplitudes = [0, 1e-6, 1e-4, 0.01, 0.1, 0.5]
    runs = [run(epsilon) for epsilon in amplitudes]
    refinement = [run(0.1, rtol=tolerance) for tolerance in (1e-8, 1e-10, 1e-12)]
    independent = run(0.1, rtol=1e-11, method="Radau")
    refinement_errors = [max(item["max_component_scaled_error"].values())
                         for item in refinement]
    checks = {
        "all_finite": all(np.isfinite(list(item["endpoint"].values())).all() for item in runs),
        "reference_accuracy": all(max(item["max_component_scaled_error"].values()) < 1e-6
                                  for item in runs),
        "constraint_accuracy": all(item["max_normalized_constraint_residual"] < 1e-10
                                   for item in runs + [independent]),
        "vacuum_zero_ricci": runs[0]["ricci_scalar_endpoint"] == 0,
        "refinement_improves": refinement_errors[2] < refinement_errors[1] < refinement_errors[0],
        "independent_solver_accuracy": max(independent["max_component_scaled_error"].values()) < 1e-6,
    }
    result = {
        "model": "Compact, locally rotationally symmetric Bianchi I with canonical massless scalar",
        "software_versions": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy_version},
        "units": "8*pi*G=c=1; initial longitudinal contraction rate magnitude=1",
        "metric": "ds^2=-d(tau)^2+A(tau)^2*d(theta)^2+B(tau)^2*(dx^2+dy^2); theta periodic",
        "constraint": "2*H_A*H_B+H_B^2-Pi^2/2=0",
        "scope": "Fully coupled homogeneous sector. Not a stability analysis of the curved Ori core.",
        "time_coordinate": "lambda=-log(q), q=1+Theta_0*tau; analytic q sets lapse only",
        "error_normalization": "Relative to each exact component; phi denominator is floored at (epsilon/abs(Theta_0))*1e-6. Exactly zero reference components use absolute error.",
        "conclusion": "Finite-endpoint numerics reproduce exact solutions. For every nonzero tested epsilon, the exact solution has R=-epsilon^2/q^2 and a curvature singularity as q approaches zero. This asymptotic statement is analytic, not a numerical evolution through q=0.",
        "vacuum_limit": "epsilon=0 is flat compact Milne with a possible Misner extension; vanishing Ricci alone would not establish flatness for arbitrary metrics",
        "runs": runs, "refinement": refinement, "independent_solver": independent,
        "self_checks": checks,
    }
    output = project_path("coupled_misner_results.json")
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "self_checks": checks,
                      "refinement_errors": refinement_errors,
                      "run_errors": [max(item["max_component_scaled_error"].values()) for item in runs]}, indent=2))
    assert all(checks.values()), "Numerical check failed; inspect recorded errors before drawing conclusions"


if __name__ == "__main__":
    main()
