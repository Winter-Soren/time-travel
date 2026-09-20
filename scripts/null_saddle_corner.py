"""All-distance local corner estimate for the bilinear null-saddle Dirac model.

This proves a quadratic-form estimate on compactly supported corner fields.
It does not complete the four overlapping scattering arms or supply their
unknown common field, interface traces, or residual. Run only this module for
small symbolic and compact-spinor checks; no scattering solve is performed.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from null_saddle_leads import AU, AV, SY
from project_paths import project_path


def corner_gap_squared(radius, eta, product_sign=None):
    """Analytic comparison minimum, not the exact operator spectrum.

    For ``product_sign`` +/-1 return the individual corner value. None gives
    the uniform bound over all four corners. A positive result is required
    for the coercivity consequences. This helper uses ordinary floats; the
    proof and separate interval evaluations establish the example constants.
    """
    if (not math.isfinite(radius) or not math.isfinite(eta)
            or radius < 1 or abs(eta) > 1
            or product_sign not in (None, -1, 1)):
        raise ValueError('Require R>=1, |eta|<=1 and corner product sign +/-1 or None.')
    e = abs(eta) if product_sign is None else product_sign*eta
    return (2*radius*radius-e)**2-2*math.sqrt(2)*radius


def _symbolic_check():
    """Differentiate the first-order operator before comparing its square."""
    u, v, eta = sp.symbols('u v eta', real=True)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    au, av = (sx+sz)/sp.sqrt(2), (sx-sz)/sp.sqrt(2)
    mass = eta-2*u*v

    def operator(field):
        return -sp.I*au*field.diff(u)-sp.I*av*field.diff(v)+mass*sy*field

    field = sp.Matrix([(1+sp.I*u)*v**2+u**3, u*v+sp.I*(v**3+u)])
    gradient = 2*sp.I*(v*au+u*av)*sy
    rhs = -field.diff(u, 2)-field.diff(v, 2)+mass**2*field+gradient*field
    residual = (operator(operator(field))-rhs).applyfunc(sp.expand)
    gradient_square = (gradient**2-4*(u*u+v*v)*sp.eye(2)).applyfunc(sp.expand)
    matrix_error = max(np.max(abs(np.array(a.evalf(), complex)-b))
                       for a, b in ((au, AU), (av, AV), (sy, SY)))
    assert residual == sp.zeros(2, 1) and gradient_square == sp.zeros(2)
    assert gradient.conjugate().T == gradient and matrix_error < 3e-16
    # The derivative lower bound is exact once x,y>=R and 2R^2-e>0.
    x, y, e = sp.symbols('x y e', real=True, positive=True)
    f = (2*x*y-e)**2-2*sp.sqrt(x*x+y*y)
    dx = 4*y*(2*x*y-e)-2*x/sp.sqrt(x*x+y*y)
    assert sp.simplify(sp.diff(f, x)-dx) == 0
    # At R=6/5, eta=1 the uniform squared gap is strictly positive.
    threshold_gap = sp.Rational(2209, 625)-sp.Rational(12, 5)*sp.sqrt(2)
    assert sp.Rational(2209, 625)**2-2*sp.Rational(12, 5)**2 > 0
    return dict(first_order_square_residual=['0', '0'],
                gradient_square_residual='zero 2x2 matrix',
                actual_matrix_maximum_error=float(matrix_error),
                exact_uniform_gap_at_R_6_over_5=str(threshold_gap))


def _interval_constants():
    """Outward interval arithmetic for constants only; minimization is analytic."""
    old = mp.iv.dps
    mp.iv.dps = 70
    rows = []
    try:
        for radius in ('1.2', '2.5', '3'):
            r = mp.iv.mpf(radius)
            for eta in ('0', '1'):
                e = mp.iv.mpf(eta)
                gap2 = (2*r*r-e)**2-2*mp.iv.sqrt(2)*r
                assert gap2.a > 0
                rows.append(dict(R=radius, absolute_eta=eta,
                                 uniform_gap_squared_interval=str(gap2),
                                 uniform_gap_interval=str(mp.iv.sqrt(gap2))))
    finally:
        mp.iv.dps = old
    return rows


def _compact_spinor_check():
    """Direct H norm versus integrated H^2, using no spectral discretization.

    A C2 compact polynomial bump lies strictly inside each corner. The constant
    spinor minimizes the Pauli-gradient matrix at the bump center. Gauss
    quadrature is exact for these polynomial integrands in exact arithmetic.
    """
    nodes, weights = np.polynomial.legendre.leggauss(14)
    p, q = np.meshgrid(nodes, nodes, indexing='ij')
    quadrature = .4**2*np.outer(weights, weights)
    bump = (1-p*p)**3*(1-q*q)**3
    bp = -6*p*(1-p*p)**2*(1-q*q)**3/.4
    bq = -6*q*(1-q*q)**2*(1-p*p)**3/.4
    radius, eta, center = 2.5, 1., 3.1
    records = []
    for su, sv in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        u, v = su*(center+.4*p), sv*(center+.4*q)
        gcenter = 2j*(sv*center*AU+su*center*AV)@SY
        spin = np.linalg.eigh(gcenter)[1][:, 0]
        du, dv = su*bp, sv*bq
        mass = eta-2*u*v
        direct = (-1j*du[..., None]*(AU@spin)-1j*dv[..., None]*(AV@spin)
                  +(mass*bump)[..., None]*(SY@spin))
        norm2 = float(np.sum(quadrature*bump*bump))
        hnorm2 = float(np.sum(quadrature*np.sum(abs(direct)**2, axis=-1)))
        grad_expect = (2j*(v*(spin.conj()@AU@SY@spin)
                          +u*(spin.conj()@AV@SY@spin))).real
        form = float(np.sum(quadrature*(du*du+dv*dv
                             +(mass*mass+grad_expect)*bump*bump)))
        gap2 = corner_gap_squared(radius, eta, su*sv)
        error = abs(hnorm2-form)/max(hnorm2, form)
        assert error < 2e-13 and hnorm2 > gap2*norm2
        records.append(dict(corner_signs=[su, sv], R=radius, eta=eta,
                            support_absolute_coordinate_interval=[2.7, 3.5],
                            psi_norm_squared=norm2, direct_H_norm_squared=hnorm2,
                            integrated_square_form=form, relative_identity_error=error,
                            rayleigh_quotient=hnorm2/norm2,
                            corner_comparison_gap_squared=gap2))
    return records


def run_checks():
    symbolic = _symbolic_check()
    constants = _interval_constants()
    compact = _compact_spinor_check()
    tests = dict(exact_first_order_square=True, exact_gradient_eigenvalue_norm=True,
                 exact_unbounded_minimum_monotonicity=True,
                 positive_interval_constants=True,
                 all_four_compact_spinor_identities=True)
    return dict(passed=all(tests.values()), tests=tests,
        model=dict(operator='H=-i AU d_u-i AV d_v+SY*(eta-2uv)',
                   matrices='AU=(sigma_x+sigma_z)/sqrt(2), AV=(sigma_x-sigma_z)/sqrt(2)',
                   units='Dimensionless local saddle coordinates, mass and energy.',
                   domain='C_c^infinity in each open corner |u|>R, |v|>R; closure in the compact-support graph norm also allowed.'),
        proof=dict(
            square='H^2=-Delta+m^2+G, G=2i(v AU+u AV)SY, G^2=4(u^2+v^2)I.',
            sign='x=|u|, y=|v|, e=eta*sign(uv); m^2=(2xy-e)^2.',
            scalar_comparison='F(x,y)=(2xy-e)^2-2sqrt(x^2+y^2).',
            unbounded_minimum='For R>=1, |eta|<=1, each partial derivative of F is >=4R(2R^2-e)-2>=2. Thus inf_{x,y>R}F=F(R,R). No outer box or sampled minimization is used.',
            squared_gap='g_corner^2=(2R^2-e)^2-2sqrt(2)R; uniformly g^2=(2R^2-|eta|)^2-2sqrt(2)R.',
            positivity='R>=6/5 guarantees g^2>0 for |eta|<=1. R>=1 suffices for the minimum formula, but g^2 must separately be positive.',
            quadratic_form='||H phi||^2>=||grad phi||^2+g^2||phi||^2 for compact corner fields.',
            not_exact_spectrum='g is a coercivity lower bound, not the bottom eigenvalue of a chosen self-adjoint corner boundary problem.'),
        cutoff_consequences=dict(
            commutator='[H,chi]=-i(AU chi_u+AV chi_v), with pointwise matrix norm |grad chi|.',
            zero_energy='g ||chi psi|| <= ||chi Hpsi||+|| |grad chi| psi||.',
            shifted_energy='For E=xi+i*kappa, |xi|<g, set gamma_E=sqrt((g-|xi|)^2+kappa^2). Then gamma_E ||chi psi||<=||chi(H-E)psi||+|| |grad chi| psi||.',
            shifted_proof='For compact phi, symmetry gives ||(H-xi-i*kappa)phi||^2=||(H-xi)phi||^2+kappa^2||phi||^2; use ||(H-xi)phi||>=(g-|xi|)||phi||.',
            nested_corner='For a product linear ramp of width ell>0, and a field/residual admitting the infinity cutoff limit: ||psi||_{C_(R+ell)} <= [||r||_{C_R}+sqrt(2)/ell*||psi||_{C_R minus C_(R+ell)}]/gamma_E, r=(H-E)psi.',
            ramp_conditions='The ramp has zero trace at the two inner edges, equals one in the deeper corner, and |grad chi|<=sqrt(2)/ell. Its Lipschitz/graph approximation is permitted; finite corner L2 norm and residual justify removing an outer cutoff.',
            missing_data='The transition-strip norm or interface trace and a globally conforming common field with controlled residual are still needed. Independent whole-line arms do not supply them; a trace jump contributes an unaccounted distributional residual.'),
        interval_digits=70, constants=constants, symbolic=symbolic,
        compact_spinor_control=compact,
        scope=[
            'This is an all-distance estimate inside the four open corners of the stated local bilinear model.',
            'No boundary condition has been imposed on the physical scattering field at a corner edge. Compact cutoffs introduce the explicit gradient/transition term.',
            'Large corner mass alone gives no unconditional bound on corner amplitude, arm-overlap error or a scattering matrix. No exponential suppression or completed full-plane exterior is asserted.',
            'The local normal form, not the full curved dihole geometry, is being checked. No supporting quantum stress, stability, formation or time-travel conclusion follows.'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                      for name in ('null_saddle_leads.py', 'project_paths.py')})


if __name__ == '__main__':
    result = run_checks()
    project_path('null_saddle_corner_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'],
                         uniform_gap_squared_R2_5_eta1=corner_gap_squared(2.5, 1.),
                         largest_compact_identity_error=max(
                             x['relative_identity_error'] for x in result['compact_spinor_control'])), indent=2))
