"""Israel layer required by a sharp synthetic AdS2 x S2 -> flat join.

Conventions: signature -+++, c=1, G left symbolic; K_ab=e_a^mu e_b^nu
nabla_mu n_nu; one normal goes from throat (-) to exterior (+) on BOTH
sides of each interface.  [K]=K_plus-K_minus.  The primary convention is
Poisson & Visser, gr-qc/9506083, equations (3), (6)-(8):
https://arxiv.org/pdf/gr-qc/9506083 .  The interior product metric is MMP
1807.04726, equation (5.18): https://arxiv.org/pdf/1807.04726 .

This evaluates the singular geometric stress of a prescribed sharp join.
It is not a smooth MMP collar, a matter realization, or a solution of the
bulk Einstein-Maxwell equations with the project's test magnetic field.
Run with --write to save this audit's JSON; importing or calling checks
does not change any file.  No SI length or physical value of G is assigned.
"""

from __future__ import annotations
from project_paths import project_path

import argparse
import json
import math
from pathlib import Path

import sympy as sp
from scipy.integrate import quad


def _connection_cut(metric, coordinates, normal):
    """Pull back nabla n to a static radial cut, then raise with its metric."""
    inverse = metric.inv()
    covector = metric * normal
    tangent = (0, 2, 3)
    connection = {}
    for a in tangent:
        for b in tangent:
            for c in range(4):
                connection[c, a, b] = sp.simplify(sum(
                    inverse[c, d] * (
                        sp.diff(metric[d, b], coordinates[a])
                        + sp.diff(metric[d, a], coordinates[b])
                        - sp.diff(metric[a, b], coordinates[d])
                    ) / 2 for d in range(4)
                ))
    curvature = sp.Matrix(3, 3, lambda i, j: sp.simplify(
        sp.diff(covector[tangent[j]], coordinates[tangent[i]])
        - sum(connection[c, tangent[i], tangent[j]] * covector[c]
              for c in range(4))
    ))
    induced = metric.extract(tangent, tangent)
    assert sp.simplify((normal.T * metric * normal)[0] - 1) == 0
    assert all(sp.simplify(covector[a]) == 0 for a in tangent)
    return induced, (induced.inv() * curvature).applyfunc(sp.simplify)


def _smooth_collar_check():
    """Independent 4D curvature and finite-collar affine-null integral check."""
    t, l, theta, phi = sp.symbols("t l theta phi", real=True)
    coordinates = (t, l, theta, phi)
    radius, lapse = sp.Function("r")(l), sp.Function("N")(l)
    metric = sp.diag(-lapse**2, 1, radius**2, radius**2 * sp.sin(theta)**2)
    inverse = metric.inv()
    gamma = {(a, b, c): sp.simplify(sum(inverse[a, d] * (
        sp.diff(metric[d, c], coordinates[b])
        + sp.diff(metric[d, b], coordinates[c])
        - sp.diff(metric[b, c], coordinates[d])) / 2 for d in range(4)))
        for a in range(4) for b in range(4) for c in range(4)}
    ricci = [sp.simplify(sum(
        sp.diff(gamma[k, i, i], coordinates[k])
        - sp.diff(gamma[k, i, k], coordinates[i])
        + sum(gamma[k, k, j] * gamma[j, i, i]
              - gamma[k, i, j] * gamma[j, i, k] for j in range(4))
        for k in range(4))) for i in range(4)]
    scalar = sp.simplify(sum(inverse[i, i] * ricci[i] for i in range(4)))
    energy = sp.simplify(ricci[0] / lapse**2 + scalar / 2)
    radial_pressure = sp.simplify(ricci[1] - scalar / 2)
    rp, rpp, np_ = sp.diff(radius, l), sp.diff(radius, l, 2), sp.diff(lapse, l)
    assert sp.simplify(energy - (1 - rp**2 - 2 * radius * rpp) / radius**2) == 0
    assert sp.simplify(radial_pressure - 2 * np_ * rp / (lapse * radius)
                       - (rp**2 - 1) / radius**2) == 0
    assert sp.simplify(energy + radial_pressure
                       + 2 * lapse / radius * sp.diff(rp / lapse, l)) == 0
    # Numerical control is a prescribed smooth geometry, not a matter model.
    b, width = 0.2, 0.1
    def profile(distance):
        s = distance / width
        return b + width * (s**3 - s**4 / 2), 3 * s**2 - 2 * s**3, (6 * s - 6 * s**2) / width
    def direct(distance):
        r, dr, ddr = profile(distance)
        g_energy = (1 - dr**2 - 2 * r * ddr) / r**2
        g_radial_pressure = (dr**2 - 1) / r**2
        return (g_energy + g_radial_pressure) / 2
    lhs, lhs_error = quad(direct, 0, width, epsabs=1e-13, epsrel=1e-13)
    squared, squared_error = quad(lambda distance: (profile(distance)[1] / profile(distance)[0])**2,
                                  0, width, epsabs=1e-13, epsrel=1e-13)
    r_out = profile(width)[0]
    rhs = -1 / r_out - squared
    assert profile(0) == (b, 0.0, 0.0)
    assert profile(width)[1:] == (1.0, 0.0)
    assert abs(lhs - rhs) < 1e-11 and lhs < 0 and squared >= 0
    return {
        "independent_symbolic_curvature": "Full 4D Christoffel/Ricci contraction reproduces both orthonormal Einstein components exactly.",
        "energy": "8*pi*G*rho=(1-r_prime^2-2*r*r_second)/r^2",
        "radial_pressure": "8*pi*G*p_l=2*N_prime*r_prime/(N*r)+(r_prime^2-1)/r^2",
        "radial_null_identity": "rho+p_l=-N/(4*pi*G*r)*(r_prime/N)_prime",
        "affine_null_vector": "k^t=E/N^2, k^l=E/N, E>0 constant",
        "finite_collar_integral": "I=int T_ab*k^a*k^b dLambda=-E/(4*pi*G)*([r_prime/(N*r)]_inner^outer+int r_prime^2/(N*r^2) dl)",
        "strict_sign_conditions": "Static spherical smooth collar, N>0, r>0, inner r_prime=0, outer r_prime=N=1: I<0.",
        "scope": "Smoothing with these endpoint conditions still requires a negative radial-null stress average inside the collar. This is a finite-segment identity, not a complete-geodesic ANEC statement, instability result or no-go for quantum-supported smooth collars.",
        "quadrature_control": {
            "profile": "N=1; r=b+w*(s^3-s^4/2), s=l/w in [0,1]",
            "scope": "Checks the integral identity and radial endpoint data, not the nonzero inner lapse derivative of the AdS throat.",
            "b": b, "width": width, "r_out": r_out,
            "four_pi_G_over_E_integral_direct": lhs,
            "four_pi_G_over_E_integral_boundary_formula": rhs,
            "four_pi_G_b_over_E_integral": b * lhs,
            "positive_squared_derivative_integral": squared,
            "absolute_difference": abs(lhs - rhs),
            "reported_quad_errors": [lhs_error, squared_error],
        },
    }


def throat_junction_checks():
    """Recompute exact geometry, Israel tensor and three dimensionless cases."""
    t, rho, theta, phi, r, proper = sp.symbols(
        "t rho theta phi r proper", real=True
    )
    b, ell, cut, G = sp.symbols("b ell rho_c G", positive=True)
    throat = sp.diag(-b**2 * (1 + rho**2) / ell**2,
                     b**2 / (1 + rho**2), b**2, b**2 * sp.sin(theta)**2)
    flat = sp.diag(-1, 1, r**2, r**2 * sp.sin(theta)**2)
    h_out, k_out = _connection_cut(
        flat, (t, r, theta, phi), sp.Matrix([0, 1, 0, 0])
    )
    h_out, k_out = h_out.subs(r, b), k_out.subs(r, b)
    calibration = {ell: b * sp.sqrt(1 + cut**2)}
    expected_inner = sp.diag(cut / (b * sp.sqrt(1 + cut**2)), 0, 0)
    expected_outer = sp.diag(0, 1 / b, 1 / b)
    assert (k_out - expected_outer).applyfunc(sp.simplify) == sp.zeros(3)
    for side in (-1, 1):
        h_in, k_in = _connection_cut(
            throat, (t, rho, theta, phi),
            sp.Matrix([0, side * sp.sqrt(1 + rho**2) / b, 0, 0])
        )
        h_in = h_in.subs(rho, side * cut).subs(calibration)
        k_in = k_in.subs(rho, side * cut).subs(calibration)
        assert (h_in - h_out).applyfunc(sp.simplify) == sp.zeros(3)
        assert (k_in - expected_inner).applyfunc(sp.simplify) == sp.zeros(3)
        # Independent proper-distance chart: rho=sinh(proper/b), g_ll=1.
        lapse = b * sp.cosh(proper / b) / ell
        gaussian_metric = sp.diag(-lapse**2, b**2, b**2 * sp.sin(theta)**2)
        gaussian_k = gaussian_metric.inv() * side * gaussian_metric.diff(proper) / 2
        gaussian_k = gaussian_k.subs(proper, side * b * sp.asinh(cut))
        assert (gaussian_k - k_in).applyfunc(sp.simplify) == sp.zeros(3)
    jump = k_out - expected_inner
    surface = -(jump - sp.eye(3) * sp.trace(jump)) / (8 * sp.pi * G)
    sigma, pressure = -surface[0, 0], surface[1, 1]
    sigma_hat = sp.simplify(4 * sp.pi * G * b * sigma)
    pressure_hat = sp.simplify(8 * sp.pi * G * b * pressure)
    nec_hat = sp.simplify(8 * sp.pi * G * b * (sigma + pressure))
    assert sigma_hat == -1
    assert sp.simplify(pressure_hat - 1 + cut / sp.sqrt(1 + cut**2)) == 0
    assert sp.simplify(nec_hat + 1 + cut / sp.sqrt(1 + cut**2)) == 0
    # Intrinsic null tangent verifies the NEC uses sigma+p, not sigma-p.
    null = sp.Matrix([1, 1 / b, 0])
    assert sp.simplify((null.T * h_out * null)[0]) == 0
    assert sp.simplify((null.T * h_out * surface * null)[0] - sigma - pressure) == 0
    # A smooth flat-to-flat identification with the same normal has no layer.
    flat_gaussian = sp.diag(-1, (b + proper)**2,
                            (b + proper)**2 * sp.sin(theta)**2)
    flat_inner = (flat_gaussian.inv() * flat_gaussian.diff(proper) / 2).subs(proper, 0)
    flat_jump = (k_out - flat_inner).applyfunc(sp.simplify)
    assert flat_jump - sp.eye(3) * sp.trace(flat_jump) == sp.zeros(3)
    examples = []
    for value in (2, 5, 10):
        radius = 0.2
        ell_value = radius * math.hypot(1, value)
        ratio = value / math.hypot(1, value)
        phat = 1 / (math.hypot(1, value) * (math.hypot(1, value) + value))
        assert abs(phat - (1 - ratio)) < 3e-16
        assert phat > 0 and -1 - ratio < 0
        examples.append({
            "rho_c": value, "b": radius, "ell": ell_value,
            "N_cut": radius / ell_value * math.hypot(1, value),
            "K_inner_mixed_times_b": [ratio, 0.0, 0.0],
            "K_outer_mixed_times_b": [0.0, 1.0, 1.0],
            "four_pi_G_b_sigma": -1.0,
            "eight_pi_G_b_pressure": phat,
            "eight_pi_G_b_sigma_plus_pressure": -1 - ratio,
        })
    return {
        "status": "PASS: geometric sharp-junction audit",
        "conventions": {
            "metric": "signature -+++, c=1, G symbolic",
            "normal": "throat toward exterior at each end; same oriented normal in both limiting geometries",
            "extrinsic": "K_ab=e_a^mu e_b^nu nabla_mu n_nu",
            "jump": "K_exterior-K_throat",
            "surface_tensor": "S^a_b=-(jumpK^a_b-delta^a_b trace(jumpK))/(8*pi*G)=diag(-sigma,p,p)",
            "units": "b=0.2 is in the existing chosen coordinate length units; no SI scale assigned",
        },
        "exact": {
            "calibration": "ell=b*sqrt(1+rho_c^2), N_cut=1, r_cut=b",
            "inner_K_mixed": "diag(rho_c/ell,0,0), at BOTH ends with the stated normals",
            "outer_K_mixed": "diag(0,1/b,1/b)",
            "sigma": "-1/(4*pi*G*b)",
            "pressure": "(1/b-rho_c/ell)/(8*pi*G)",
            "sigma_plus_pressure": "-(1/b+rho_c/ell)/(8*pi*G)",
            "dimensionless_sigma": "4*pi*G*b*sigma=-1",
            "dimensionless_pressure": "8*pi*G*b*p=1-rho_c/sqrt(1+rho_c^2)",
            "dimensionless_nec": "8*pi*G*b*(sigma+p)=-1-rho_c/sqrt(1+rho_c^2)<0",
        },
        "verification": {
            "connection_and_unit_normals": "exact symbolic; both ends",
            "first_fundamental_form": "exactly continuous after clock calibration",
            "independent_gaussian_normal_derivative": "exact symbolic; both ends",
            "israel_and_intrinsic_null_contraction": "exact symbolic",
            "smooth_flat_identification_control": "zero jump and zero layer",
            "numerical_examples": "same symbolic expressions; cancellation-safe positive pressure",
        },
        "examples": examples,
        "smooth_collar_gate": _smooth_collar_check(),
        "interpretation": {
            "required_layer": "Negative surface energy and tangential null-energy violation are required by this sharp prescribed metric join in Einstein gravity.",
            "not_established": [
                "A material or quantum state that provides this singular surface stress.",
                "The total stress of the nonbackreacting magnetic scattering exterior: flat space with a retained test magnetic field is not its Einstein-Maxwell solution.",
                "Bulk Einstein-Maxwell consistency, surface electromagnetic currents, or stability.",
                "A no-go theorem for a smooth MMP collar; smoothing changes the geometry and stress problem.",
            ],
        },
        "primary_sources": [
            {"title": "Poisson and Visser, Thin-shell wormholes: Linearization stability (1995)",
             "url": "https://arxiv.org/pdf/gr-qc/9506083",
             "equations": "(3), (6)-(8)",
             "use": "Extrinsic-curvature convention and Israel surface-stress relation; actual full text read."},
            {"title": "Maldacena, Milekhin and Popov, Traversable wormholes in four dimensions",
             "url": "https://arxiv.org/pdf/1807.04726",
             "equations": "(5.18)-(5.20)",
             "use": "Leading product throat geometry only; this sharp flat matching is our synthetic calculation."},
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="save the focused audit JSON")
    args = parser.parse_args()
    result = throat_junction_checks()
    if args.write:
        target = project_path("throat_junction_audit_results.json")
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, allow_nan=False))
