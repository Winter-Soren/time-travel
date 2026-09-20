"""Analytic quadrature remainders for the exact central A4/N18 sine model.

This helper performs no physical integration and writes no result file.  It
bounds the operator-norm error of EXACT panelwise tensor Gauss--Legendre
quadrature for M=<H Phi,H Phi> and the full Gram of
Rtilde=A Phi-Phi Mhat, with Mhat the frozen123 central matrix.  Node, weight,
coefficient, contraction and summation errors must be enclosed separately.

Proof of the constants used below:

* Reuse analytic_geometry_bounds on its documented complex physical box.
  Positive real parts of the root radicands select analytic principal roots;
  rho and Sigma stay nonzero.  The combined mass numerator has no division
  by the lapse numerator F, so a zero of F creates no artificial pole.
* A width-h Bernstein ellipse with parameter2 has imaginary extent3h/8
  and extends h/8 beyond a real panel.  A u/v polydisc of radius1/512 around
  either one-coordinate ellipse fits in the old physical box.  Cauchy then
  gives |w_i|<=W/delta, |Delta w|<=4W/delta^2 and |m_i|<=U/delta.
* s_n(x)=sin(n*pi*(x+4)/8)/2 has unit L2 norm on[-4,4].  With only one
  coordinate complex, |sin|,|cos|<=exp(kmax*|Im x|).  The two-spin basis
  row Phi has Frobenius norm <=sqrt(2)*18/4*exp(kmax*3h/8).
* H=iK with real K on the real box.  The analytic matrix integrands are
  (K Phi)^T(K Phi) and Rtilde^T Rtilde, NOT conjugation at complex arguments.
  Their operator norms are at most the squared Frobenius bounds below.
* Analytic Chebyshev coefficients satisfy ||a_j||<=2M*2^-j.  Truncation at
  degree2n-1 and positivity/exactness of GLn give one-panel error <=h*delta_n*M,
  delta_n=4*2^(1-2n).  Tensorization and summing panel areas give
  ||I-Q||<=2*64*delta_n*M.  The same contour/triangle proof works for matrix
  operator norms, without another matrix-dimension factor.

Primary quadrature reference: Trefethen, SIAM Review50(2008), Theorems4.2--4.3,
https://people.maths.ox.ac.uk/~trefethen/publication/PDF/2008_127.pdf .
The geometry and differential expression are the exact-balanced Emparan
target already identified in dihole_probe_certificate.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

import dihole_interval_arithmetic as ia
from dihole_continuum_heat_tail import target
from dihole_quadrature_remainder import analytic_geometry_bounds, dyadic, outward
from project_paths import project_path


FROZEN_STATE_SHA256 = '218545a22a2f1d9d757f7d7535c5534c82a2233e1cdac0b938dd885ceed01bdd'


def analytic_control(iv):
    """GL2 tensor test with a known logarithmic integral, no fitted error."""
    p = iv.mpf
    node = 1/iv.sqrt(3)
    q1 = 1/(3-node)+1/(3+node)
    actual_error = iv.log(2)**2-q1**2
    # On ellipse2 x realpart<=5/4, hence |3-x|>=7/4; the other
    # coordinate is real and |3-y|>=2.  Each strip modulus is <=2/7.
    strip_modulus = p(2)/7
    rule_delta = 4*p(2)**(1-2*2)
    upper = 2*4*rule_delta*strip_modulus
    assert abs(q1-p(9)/13) < p('1e-45')
    assert actual_error > 0 and actual_error < upper
    return dict(passed=True, exact_integral='log(2)^2',
                exact_GL2_sum='(9/13)^2',
                error=outward(actual_error), upper_bound=outward(upper))


def remainder_bounds(order=64, panel_width=2, dps=75):
    """JSON-safe analytic remainder bounds for this fixed central sector.

    Returned lower/upper pairs enclose positive majorant constants; they are
    not a signed enclosure of the actual quadrature error.  The upper values
    may be added to an independently enclosed exact quadrature matrix sum.
    """
    if isinstance(order, bool) or not isinstance(order, int) or order < 1:
        raise ValueError('order must be a positive integer.')
    if (isinstance(panel_width, bool) or not isinstance(panel_width, int)
            or panel_width not in (1, 2)):
        raise ValueError('This bounded helper supports panel widths 1 or 2 only.')
    if isinstance(dps, bool) or not isinstance(dps, int) or dps < 50:
        raise ValueError('Use at least 50 interval decimal digits.')
    iv = MPIntervalContext(); iv.dps = dps
    p = iv.mpf
    count, radius, dimension = 18, 4, 648
    h, delta = p(panel_width), p(1)/512
    constants = target(iv)
    L, rs, vs, E = (constants[k] for k in ('L', 'rho_star', 'vstar', 'E'))

    # These are the exact balanced radicals, not independently exact stored
    # floats.  They verify applicability of the older conservative box.
    k, Q = p(1)/2, p(1)/5
    M = (iv.sqrt(k*k+4*Q*k)-k)/2
    ag = iv.sqrt(k*k-M*M)
    field = M/((k+ag)*(M+k))
    assert p('.153') < M < p('.154') and 0 < ag < p('.476')
    assert 0 < field < p('.241') and L < p('.011762')
    assert E > p('43.3') and vs > p('.509')
    mu, q = p(797711)/2, 512000
    assert mu < 400000 and q == 512000

    # Include BOTH complex Cauchy increments.  This is stronger than the
    # separate one-variable disks sufficient for pure second derivatives.
    imaginary = 3*h/8
    physical_real_extent = L*(2*radius+h/8+2*delta)/iv.sqrt(2)
    physical_imaginary_extent = L*(imaginary+2*delta)/iv.sqrt(2)
    assert rs-physical_real_extent > p('.515')
    assert rs+physical_real_extent < p('.820')
    assert physical_real_extent < p('.152')
    assert physical_imaginary_extent < p('.0063')
    geometry = analytic_geometry_bounds(iv)
    W, U = (geometry[k] for k in ('normalized_speed_upper', 'normalized_mass_upper'))
    assert W < 8 and U < 500000
    first_speed = W/delta
    lap_speed = 4*W/(delta*delta)
    first_mass = U/delta

    frozen = project_path('dihole_heat_angular_state.npz')
    state_hash = hashlib.sha256(frozen.read_bytes()).hexdigest()
    if state_hash != FROZEN_STATE_SHA256:
        raise ValueError('The frozen123 matrix state has changed.')
    with np.load(frozen) as saved:
        matrix = saved['M0'].copy()
    if (matrix.shape != (dimension, dimension) or matrix.dtype != np.float64
            or not np.isfinite(matrix).all() or not np.array_equal(matrix, matrix.T)):
        raise ValueError('Expected the finite real symmetric 648x648 binary64 M0.')
    # Directed pairwise summation, not an opaque floating sum or eigensolve.
    norm2 = ia.sum_interval(ia.square(ia.point(matrix)))
    norm2_iv = iv.mpf([dyadic(iv, float(norm2[0])).a,
                       dyadic(iv, float(norm2[1])).b])
    matrix_frobenius = iv.sqrt(norm2_iv)

    frequency = count*iv.pi/(2*radius)
    phi = iv.sqrt(2)*count/radius*iv.exp(frequency*imaginary)
    # For each coordinate ||partial_i Phi||F<=frequency*phi and
    # ||Delta Phi||F<=2*frequency^2*phi.
    K_factor = 2*W*frequency+first_speed+U
    # Exact strong A=-K^2: -w^2 Delta-2w gradw.grad plus
    # [m^2-w Delta(w)/2-|gradw|^2/4], the spin commutator, and
    # w(m_v AU-m_u AV).  Absolute values yield these six positive terms.
    A_factor = (2*W**2*frequency**2+6*W*first_speed*frequency+U**2
                +W*lap_speed/2+first_speed**2/2+2*W*first_mass)
    K_phi = K_factor*phi
    residual = (A_factor+matrix_frobenius)*phi
    M_modulus, G_modulus = K_phi**2, residual**2
    rule_delta = 4*p(2)**(1-2*order)
    area = p(4*radius*radius)
    M_error, G_error = 2*area*rule_delta*M_modulus, 2*area*rule_delta*G_modulus
    control = analytic_control(iv)
    files = ['dihole_heat_angular_state.npz', 'dihole_probe_certificate_results.json',
             'dihole_quadrature_remainder.py', 'dihole_continuum_heat_tail.py',
             'dihole_interval_arithmetic.py', 'project_paths.py']
    return dict(
        passed=True,
        target=dict(radius=radius, sine_count=count, dimension=dimension,
                    q=q, mu_center=398855.5, frozen_matrix='M0',
                    matrix_sha256=hashlib.sha256(matrix.tobytes(order='C')).hexdigest(),
                    coordinates={k:outward(constants[k]) for k in ('L', 'rho_star', 'vstar', 'E')},
                    exact_energy='E=vstar/L, with L,rho_star,vstar exact saved dyadics.',
                    exact_gravity='D=1,Q=1/5,k=1/2; M=(sqrt(k^2+4Qk)-k)/2; ag=sqrt(k^2-M^2); B=M/((k+ag)(M+k)).'),
        quadrature=dict(order_per_panel=order, panel_width=panel_width,
                        panel_edges=list(range(-radius, radius+1, panel_width)),
                        total_nodes_per_axis=(2*radius//panel_width)*order,
                        ellipse_parameter=2, cauchy_polydisc_radius='1/512',
                        analytic_rule_delta=outward(rule_delta)),
        geometry_box=dict(rho_real=[.515,.820], z_real=[-.152,.152],
                          imaginary_extent=.0063,
                          actual_physical_real_extent=outward(physical_real_extent),
                          actual_physical_imaginary_extent=outward(physical_imaginary_extent),
                          radical_parameters={k:outward(v) for k,v in
                                              (('M',M),('ag',ag),('B',field))},
                          branch_and_denominator_constants={k:outward(geometry[k]) for k in
                              ('radicand_real_lower','root_real_lower','sigma_modulus_lower','h_modulus_lower')}),
        analytic_bounds={k:outward(v) for k,v in
            (('speed',W),('mass',U),('first_speed',first_speed),('laplacian_speed',lap_speed),
             ('first_mass',first_mass),('basis_frobenius',phi),('matrix_frobenius',matrix_frobenius),
             ('K_basis_frobenius',K_phi),('strong_A_factor',A_factor),
             ('residual_frobenius',residual),('M_integrand_operator_modulus',M_modulus),
             ('G_integrand_operator_modulus',G_modulus))},
        operator_remainders=dict(M=outward(M_error), residual_Gram=outward(G_error)),
        requested_budgets_met=dict(M_below_1e_6=bool(M_error<p('1e-6')),
                                  residual_Gram_below_1e_4=bool(G_error<p('1e-4'))),
        controls=dict(exact_parameter_boxes=True, complex_polydiscs_inside_old_box=True,
                      frozen_matrix_provenance=True, directed_frobenius_norm=True,
                      independent_analytic_GL2_control=control),
        residual_interpretation=('The integrated Gram is of A Phi-Phi Mhat. It equals '
            'R*R+(Mexact-Mhat)^*(Mexact-Mhat), so it is a valid upper surrogate for '
            'the true complementary coupling Gram. This identity uses the exact orthonormal sine basis.'),
        scope=('Only analytic truncation of exact tensor Gauss--Legendre quadrature is bounded. '
               'No quadrature sum or arithmetic/integrand evaluation error is certified here. '
               'This is one central angular sector of the Dirichlet box, not pointwise heat or the natural exterior.'),
        interval_digits=dps,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in files})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, default=64)
    parser.add_argument('--panel-width', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    print(json.dumps(remainder_bounds(args.order, args.panel_width), indent=2))
