"""Local, constrained linear Einstein--Polyakov response on a null diamond.

This is a conditional radius-independent Polyakov model, not mouth matching.
Run with .venv/Scripts/python.exe -B polyakov_response.py.  No parent artifacts
are changed; the sole output is polyakov_response_results.json.
"""
from __future__ import annotations
from project_paths import project_path

import hashlib
import inspect
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from mmp_pulse import _polyakov_static_solution


def _bump_derivatives(z):
    z = np.asarray(z, dtype=float)
    bp = np.zeros_like(z)
    bpp = np.zeros_like(z)
    bppp = np.zeros_like(z)
    inside = (z > 0) & (z < 1)
    w = z[inside] * (1 - z[inside])
    b = np.exp(4 - 1 / w)
    bp[inside] = b * (1 - 2 * z[inside]) / w**2
    bpp[inside] = b * (1 - 6 * w + 6 * w**2) / w**4
    bppp[inside] = (b * (1 - 2 * z[inside])
                   * (1 - 10 * w + 24 * w**2 - 12 * w**3) / w**6)
    return bp, bpp, bppp


def _state_direction(u, omega2, start, width):
    """d t_u / d epsilon for V(u)=u+epsilon*width*B((u-start)/width)."""
    bp, _, bppp = _bump_derivatives((np.asarray(u) - start) / width)
    return -(bppp / width**2 + omega2 * bp) / 2


class StaticBackground:
    def __init__(self, alpha, omega2, endpoint):
        self.alpha, self.k = alpha, alpha / 3
        self.solution, self.rhs = _polyakov_static_solution(
            alpha, omega2, endpoint, refinement=2)

    def fields(self, x):
        x = np.asarray(x, dtype=float)
        raw = self.solution.sol(np.abs(x).ravel()).reshape((4,) + x.shape)
        y, yp, z, zp = raw
        sec2 = 1 / np.cos(x)**2
        sphere = 1 + self.alpha * y
        radial = np.sqrt(sphere)
        exponential = sec2 * np.exp(2 * self.alpha * z)
        # The normalized static equations avoid subtracting O(1) terms to
        # recover the physical O(alpha) radius and conformal corrections.
        a = (exponential * (y / sphere - 1 / (3 * sphere**2))
             - self.alpha * yp**2 / (4 * sphere)) / (sphere - self.k)
        ypp = 2 * sphere * a + self.alpha * yp**2 / (2 * sphere)
        zpp = (sec2 * np.expm1(2 * self.alpha * z
                              - 2 * np.log1p(self.alpha * y)) / self.alpha - a)
        radial_x = self.alpha * yp * np.sign(x) / (2 * radial)
        radial_xx = (self.alpha * ypp / (2 * radial)
                     - self.alpha**2 * yp**2 / (4 * radial**3))
        sigma_x = np.tan(x) + self.alpha * zp * np.sign(x)
        sigma_xx = sec2 + self.alpha * zpp
        return dict(R=radial, Rx=radial_x, Rxx=radial_xx,
                    sx=sigma_x, sxx=sigma_xx, e2s=exponential, Y=y,
                    sigma=-np.log(np.cos(x))+self.alpha*z)

    def linear_operator(self, x):
        """q_uv=M0*q+Mu*q_u+Mv*q_v for q=(delta R,delta sigma)."""
        b = self.fields(x)
        r, rx, e2, k = b['R'], b['Rx'], b['e2s'], self.k
        ru, rv = rx / 2, -rx / 2
        h = self.alpha * (b['Y'] / r**2 - 1 / (3 * r**4))
        factor = r / (r**2 - k)
        bracket = ru * rv + e2 * h / 4
        f = -factor * bracket
        ar = ((r**2 + k) / (r**2 - k)**2 * bracket
              - factor * e2 / 4 * (2 / r**3 + 4 * k / r**5))
        az = -factor * e2 * h / 2
        au, av = -factor * rv, -factor * ru
        shape = np.shape(r) + (2, 2)
        m0, mu, mv = (np.zeros(shape) for _ in range(3))
        m0[..., 0, 0], m0[..., 0, 1] = ar, az
        m0[..., 1, 0] = -ar / r + f / r**2 + e2 / r**5
        m0[..., 1, 1] = -az / r - e2 / (2 * r**4)
        mu[..., 0, 0], mu[..., 1, 0] = au, -au / r
        mv[..., 0, 0], mv[..., 1, 0] = av, -av / r
        return m0, mu, mv


def _independent_operator_checks(background, omega2, extent):
    """Check the Jacobian against exact complex jets and the static limit by
    differentiating a polynomial fit, independently of the static ODE RHS.
    """
    from numpy.polynomial import chebyshev as cheb

    def exact_bulk(jet):
        excess, sigma, ru, rv = jet
        r = 1+excess
        # Unlike 1-1/r**2, this preserves the very small radius excess.
        potential = excess*(2+excess)/r**2-background.k/r**4
        f = -r/(r*r-background.k)*(ru*rv+np.exp(2*sigma)*potential/4)
        return np.array([f, -f/r-np.exp(2*sigma)/(4*r**4)])

    errors = []
    for x in np.linspace(-.45*extent, .45*extent, 7):
        b = background.fields(x)
        excess = background.alpha*b['Y']/(b['R']+1)
        jet = np.array([excess, b['sigma'], b['Rx']/2, -b['Rx']/2], dtype=complex)
        reference = np.zeros((2, 4))
        for column in range(4):
            perturbed = jet.copy()
            perturbed[column] += 1j*1e-25
            reference[:, column] = exact_bulk(perturbed).imag/1e-25
        m0, mu, mv = background.linear_operator(x)
        jacobian = np.column_stack((m0, mu[:, 0], mv[:, 0]))
        errors.append(float(np.max(abs(reference-jacobian)/np.maximum(abs(jacobian), 1e-18))))
    assert max(errors) < 2e-11

    degree = 40
    nodes = np.cos(np.arange(degree+1)*math.pi/degree)
    x = nodes*extent/2
    raw = background.solution.sol(np.abs(x))
    points = np.linspace(-extent/2+.04, extent/2-.04, 81)
    values = []
    for row in (raw[0], raw[2]):
        coefficient = cheb.chebfit(nodes, row, degree)
        values.append([cheb.chebval(points*2/extent, cheb.chebder(coefficient, m))*(2/extent)**m
                       for m in (0, 1, 2)])
    y, yp, ypp = values[0]
    z, zp, zpp = values[1]
    alpha, k = background.alpha, background.k
    sphere = 1+alpha*y
    r = np.sqrt(sphere)
    rx_over_alpha = yp/(2*r)
    rxx_over_alpha = ypp/(2*r)-alpha*yp*yp/(4*r**3)
    sec2 = 1/np.cos(points)**2
    e2 = sec2*np.exp(2*alpha*z)
    sx, sxx = np.tan(points)+alpha*zp, sec2+alpha*zpp
    radial = rxx_over_alpha-r/(sphere-k)*(e2*(y/sphere-1/(3*sphere**2))-alpha*rx_over_alpha**2)
    angular = zpp-sec2*np.expm1(2*alpha*z-2*np.log1p(alpha*y))/alpha+rxx_over_alpha/r
    null = r*(rxx_over_alpha-2*sx*rx_over_alpha)+(sxx-sx*sx-omega2)/3
    static_errors = dict(radius_equation_over_alpha=float(np.max(abs(radial))),
                         angular_equation_over_alpha=float(np.max(abs(angular))),
                         four_times_null_constraint_over_alpha=float(np.max(abs(null))))
    assert max(static_errors.values()) < 2e-7
    return dict(complex_step_maximum_component_relative_error=max(errors),
                static_independent_chebyshev_degree=degree,
                static_independent_residuals=static_errors)


def _initial_null_edge(background, omega2, start, width, extent, strength=1.):
    """Fix delta sigma=0 as null-coordinate gauge; solve the uu constraint."""
    def rhs(u, state):
        b = background.fields(u / 2)
        r, rx, rxx, sx = (float(b[k]) for k in ('R', 'Rx', 'Rxx', 'sx'))
        source = strength * float(_state_direction(u, omega2, start, width))
        value, derivative = state
        second = (sx * derivative - (rxx / 4 - sx * rx / 2) / r * value
                  - source / (3 * r))
        return derivative, second
    sol = solve_ivp(rhs, (0, extent), (0., 0.), method='DOP853',
                    rtol=3e-12, atol=3e-13, max_step=width / 100,
                    dense_output=True)
    assert sol.success
    return sol


def _evolve(background, edge, extent, intervals):
    """Second-order implicit characteristic boxes, solved by antidiagonals."""
    n, h = intervals, extent / intervals
    coordinate = np.linspace(0, extent, n + 1)
    state = np.zeros((n + 1, n + 1, 2))
    state[:, 0, 0] = edge.sol(coordinate)[0]
    # u=0 has the unperturbed state and geometry.  On v=0 the uu
    # constraint supplies delta R; delta sigma=0 fixes conformal gauge.
    m0, mu, mv = background.linear_operator(np.arange(-n, n + 1) * h / 2)
    lhs = np.eye(2) - h*h*m0/4 - h*(mu+mv)/2
    smallest_determinant = math.inf
    for diagonal in range(2, 2*n + 1):
        i = np.arange(max(1, diagonal-n), min(n, diagonal-1) + 1)
        j = diagonal-i
        ix = i-j+n
        sw, se, nw = state[i-1, j-1], state[i, j-1], state[i-1, j]
        matvec = lambda a, v: np.einsum('nij,nj->ni', a, v)
        rhs = (se+nw-sw + h*h/4*matvec(m0[ix], sw+se+nw)
               + h/2*matvec(mu[ix], se-sw-nw)
               + h/2*matvec(mv[ix], nw-sw-se))
        matrix = lhs[ix]
        det = matrix[:, 0, 0]*matrix[:, 1, 1]-matrix[:, 0, 1]*matrix[:, 1, 0]
        smallest_determinant = min(smallest_determinant, float(np.min(abs(det))))
        state[i, j, 0] = (matrix[:, 1, 1]*rhs[:, 0]-matrix[:, 0, 1]*rhs[:, 1])/det
        state[i, j, 1] = (-matrix[:, 1, 0]*rhs[:, 0]+matrix[:, 0, 0]*rhs[:, 1])/det
    assert smallest_determinant > .5 and np.all(np.isfinite(state))
    return coordinate, state, smallest_determinant


def _constraint_check(background, state, coordinate, omega2, start, width):
    """Independent fourth-order Cartesian differences, fixed physical inset."""
    h = coordinate[1]-coordinate[0]
    f = state
    value = f[2:-2, 2:-2]
    du = (f[:-4, 2:-2]-8*f[1:-3, 2:-2]+8*f[3:-1, 2:-2]-f[4:, 2:-2])/(12*h)
    dv = (f[2:-2, :-4]-8*f[2:-2, 1:-3]+8*f[2:-2, 3:-1]-f[2:-2, 4:])/(12*h)
    duu = (-f[:-4, 2:-2]+16*f[1:-3, 2:-2]-30*value+16*f[3:-1, 2:-2]-f[4:, 2:-2])/(12*h*h)
    dvv = (-f[2:-2, :-4]+16*f[2:-2, 1:-3]-30*value+16*f[2:-2, 3:-1]-f[2:-2, 4:])/(12*h*h)
    u, v = np.meshgrid(coordinate[2:-2], coordinate[2:-2], indexing='ij')
    b = background.fields((u-v)/2)
    r, rx, rxx, sx = (b[k] for k in ('R', 'Rx', 'Rxx', 'sx'))
    k = background.k
    source = _state_direction(u, omega2, start, width)/3
    common = (rxx/4-sx*rx/2)*value[..., 0]
    cu = (r*duu[..., 0]-r*sx*du[..., 0]+common
          -r*rx*du[..., 1]+k*(duu[..., 1]-sx*du[..., 1])+source)
    cv = (r*dvv[..., 0]+r*sx*dv[..., 0]+common
          +r*rx*dv[..., 1]+k*(dvv[..., 1]+sx*dv[..., 1]))
    mask = (u >= .06) & (v >= .06) & (u <= coordinate[-1]-.06) & (v <= coordinate[-1]-.06)
    scale = 1+float(np.max(abs(source)))
    return dict(maximum_uu_absolute=float(np.max(abs(cu[mask]))),
                maximum_vv_absolute=float(np.max(abs(cv[mask]))),
                source_normalization=scale,
                maximum_source_scaled=float(max(np.max(abs(cu[mask])), np.max(abs(cv[mask])))/scale),
                physical_inset=.06,
                convention='Residual of linearized null constraint divided by alpha*epsilon; fourth-order independent differences')


def polyakov_response_checks(alpha=1e-9, exterior_delay=1.2, epsilon=1e-3,
                              extent=1.5, pulse_start=.2, pulse_width=.6,
                              intervals=(64, 128, 256, 512)):
    assert alpha > 0 and 0 < exterior_delay < math.pi
    assert 0 < pulse_start < pulse_start+pulse_width < extent < math.pi
    w0 = (3-math.sqrt(3))/6
    max_bp = math.exp(4-1/w0)*math.sqrt(1-4*w0)/w0**2
    assert 1-abs(epsilon)*max_bp > 0, 'Illustrative finite diffeomorphism must preserve orientation'
    omega2 = (2*math.pi/(math.pi+exterior_delay))**2
    background = StaticBackground(alpha, omega2, extent/2)
    operator_checks = _independent_operator_checks(background, omega2, extent)
    edge = _initial_null_edge(background, omega2, pulse_start, pulse_width, extent)
    zero_edge = _initial_null_edge(background, omega2, pulse_start, pulse_width, extent, 0.)
    _, zero, _ = _evolve(background, zero_edge, extent, 32)
    assert np.max(abs(zero)) == 0
    records, states = [], []
    for n in intervals:
        coordinates, state, determinant = _evolve(background, edge, extent, n)
        checks = _constraint_check(background, state, coordinates, omega2, pulse_start, pulse_width)
        records.append(dict(intervals=n, maximum_normalized_radius_response=float(np.max(abs(state[..., 0]))),
                            maximum_normalized_conformal_response=float(np.max(abs(state[..., 1]))),
                            final_corner=state[-1, -1].tolist(),
                            smallest_cell_determinant=determinant, constraints=checks))
        states.append(state)
    differences = [float(np.max(abs(coarse-fine[::2, ::2]))) for coarse, fine in zip(states[:-1], states[1:])]
    ratios = [differences[i]/differences[i+1] for i in range(len(differences)-1)]
    assert len(intervals) >= 3 and all(b == 2*a for a, b in zip(intervals[:-1], intervals[1:]))
    assert all(3.5 < ratio < 4.5 for ratio in ratios)
    constraints = [row['constraints']['maximum_source_scaled'] for row in records]
    assert all(b < a for a, b in zip(constraints[:-1], constraints[1:]))
    assert constraints[-1] < 1e-4
    assert differences[-1] < 1e-4
    final = states[-1]
    stride = intervals[-1]//32
    sample = final[::stride, ::stride]
    sample_coord = coordinates[::stride]
    r = background.fields((coordinates[:, None]-coordinates[None, :])/2)['R']
    area_fraction = 2*alpha*epsilon*final[..., 0]/r
    diag = final[np.arange(len(coordinates)), np.arange(len(coordinates))]
    result = dict(
        model='Linear response of spherical Einstein-Maxwell plus radius-independent Polyakov matter on a finite characteristic diamond; no physical mouth matching',
        parameters=dict(alpha=alpha, k=alpha/3, exterior_delay_used_only_to_set_static_state=exterior_delay,
                        omega_squared=omega2, infinitesimal_amplitude_epsilon=epsilon,
                        null_coordinate_interval=[0, extent], x_interval=[-extent/2, extent/2],
                        pulse_u_interval=[pulse_start, pulse_start+pulse_width],
                        minimum_finite_seed_Vprime=1-abs(epsilon)*max_bp),
        conventions=dict(metric='ds4^2=-exp(2sigma)du dv+R^2 dOmega^2, re=1, u=t+x, v=t-x',
            quantum_source='t2_uu=q/(12pi)*(sigma_uu-sigma_u^2+t_u(u)); t2_vv analogously; t2_uv=-q*sigma_uv/(12pi); p_perp=0; k=Gq/(12pi)=alpha/3',
            static_state='t_u=t_v=-Omega^2/4; calibrated static R,sigma from _polyakov_static_solution',
            bulk='R_uv=-R/(R^2-k)*(R_u*R_v+exp(2sigma)*(1-R^-2-kR^-4)/4); sigma_uv=-R_uv/R-exp(2sigma)/(4R^4)',
            constraints='R*(R_uu-2sigma_u*R_u)+k*(sigma_uu-sigma_u^2+t_u)=0, and vv counterpart',
            normalization='q_response=(delta R,delta sigma)/(alpha*epsilon); all evolution is linear in the amplitude direction'),
        initial_data=dict(state='V(u)=u+epsilon*b(u), b=width*exp(4-1/[z(1-z)]), z=(u-start)/width inside(0,1), zero elsewhere in the first conformal period. Extend by degree-one periodicity with period2pi/Omega outside this diamond; delta t_u/epsilon=-(b_third+Omega^2*b_prime)/2, delta t_v=0',
            geometry='delta sigma=0 on both initial null rays is a conformal coordinate choice. The uu constraint determines delta R along v=0 with zero corner value/slope. The u=0 edge is unperturbed; its vv constraint holds. No independent overprescription of delta R and its second null derivative',
            energy='This is an infinitesimal CFT diffeomorphism-state direction. Its integrated first-order excess is zero; finite positive packet energy starts at epsilon^2 and is not evolved here'),
        method='Second-order implicit characteristic boxes; independent fourth-order finite differences test null-constraint propagation on a fixed physical inset',
        refinement=records, common_grid_state_differences=differences, observed_refinement_ratios=ratios,
        independent_operator_and_static_checks=operator_checks,
        acceptance=dict(zero_source_exact=True, state_second_order_ratios_between=[3.5, 4.5],
                        final_common_grid_difference_below=1e-4,
                        source_scaled_null_constraints_monotonically_decrease=True,
                        final_source_scaled_constraint_below=1e-4,
                        complex_step_Jacobian_relative_tolerance=2e-11,
                        independently_differentiated_static_tolerance=2e-7,
                        all_pass=True),
        zero_source_maximum_response=float(np.max(abs(zero))),
        response=dict(maximum_fractional_sphere_area_change=float(np.max(abs(area_fraction))),
                      maximum_conformal_log_change=float(abs(alpha*epsilon)*np.max(abs(final[..., 1]))),
                      final_corner_radius_change=float(alpha*epsilon*final[-1,-1,0]),
                      final_corner_conformal_change=float(alpha*epsilon*final[-1,-1,1]),
                      interpretation='Linear predictions at the stated infinitesimal amplitude; maxima are sampled-grid maxima, not continuum bounds or a nonlinear finite-amplitude validation'),
        plot=dict(u=sample_coord.tolist(), v=sample_coord.tolist(), radius_response_over_alpha_epsilon=sample[...,0].tolist(),
                  conformal_response_over_alpha_epsilon=sample[...,1].tolist(),
                  center_time=coordinates[::stride].tolist(), center_radius_response=diag[::stride,0].tolist(),
                  center_conformal_response=diag[::stride,1].tolist()),
        provenance=dict(static_helper_sha256=hashlib.sha256(inspect.getsource(_polyakov_static_solution).encode()).hexdigest()),
        scope='A local linear-response benchmark, not the eta4 history, a finite-amplitude nonlinear evolution, a source-preparation recipe, exterior matching, a boundary-clock correction, stability proof, or chronology change. The null-edge metric choice fixes coordinates; no clock offset is imposed. O(epsilon^2) response and the full four-dimensional renormalized matter action are excluded',
        sources=['https://arxiv.org/pdf/hep-th/0006035 (Polyakov action and state-dependent stress)',
                 'https://arxiv.org/pdf/hep-th/0111091 (spherical Einstein-Maxwell reduction)',
                 'Existing polyakov_static_checks and angular_metric_checks; the characteristic benchmark and linearization are our application'])
    return result


# BEGIN_FINITE_SOURCE_RESPONSE
def _finite_background(background, x):
    b = background.fields(x)
    alpha = background.alpha
    b['S'] = 1+alpha*b['Y']
    b['Yu'] = b['R']*b['Rx']/alpha
    b['Yv'] = -b['Yu']
    b['Yuv'] = -(b['Rx']**2+b['R']*b['Rxx'])/(2*alpha)
    return b


def _finite_bulk(background, b, q, qu, qv):
    """Exact normalized nonlinear equations minus the stored static solution.

    q=(Y-Ybar,Z), R^2=1+alpha*Y, sigma=sigmabar+alpha*Z.
    Every difference vanishes exactly at q=qu=qv=0; no O(1) subtraction
    is subsequently divided by the physical alpha=1e-9.
    """
    alpha, k = background.alpha, background.k
    w, z = q[..., 0], q[..., 1]
    wu, wv = qu[..., 0], qv[..., 0]
    sb, yb = b['S'], b['Y']
    s = sb+alpha*w
    assert np.min(s-k) > 0
    ab = (yb-1/(3*sb))/(sb-k)
    da = w*(1+alpha/(3*s*sb)-alpha*ab)/(s-k)
    product = b['Yu']*b['Yv']
    dproduct = b['Yu']*wv+b['Yv']*wu+wu*wv
    db = (dproduct/(s*(s-k))
          - product*alpha*w*(s+sb-k)/(s*(s-k)*sb*(sb-k)))
    fw = (-b['e2s']/2*(np.expm1(2*alpha*z)*ab+np.exp(2*alpha*z)*da)
          - alpha*alpha*db/6)
    derivative_difference = fw/s-alpha*w*b['Yuv']/(s*sb)
    product_difference = (dproduct/s**2
                          - product*alpha*w*(s+sb)/(s**2*sb**2))
    exponent = 2*alpha*z-2*np.log1p(alpha*w/sb)
    fz = (-derivative_difference/2+alpha*product_difference/4
          -b['e2s']/(4*sb**2)*np.expm1(exponent)/alpha)
    return np.stack((fw, fz), axis=-1)


def _finite_bulk_algebra_check(background):
    """Compare normalized subtractions to the original equations at 70 digits."""
    import mpmath as mp
    errors = []
    with mp.workdps(70):
        alpha = mp.mpf(float(background.alpha))
        for x in (-.6, -.2, .25, .65):
            b = _finite_background(background, x)
            q, qu, qv = (np.array(pair) for pair in ((-.7, .2), (.31, -.13), (-.27, .18)))
            y0, yu0, yv0, e0 = (mp.mpf(float(b[key])) for key in ('Y', 'Yu', 'Yv', 'e2s'))

            def original(y, yu, yv, exponential):
                s = 1+alpha*y
                r = mp.sqrt(s)
                ru, rv = alpha*yu/(2*r), alpha*yv/(2*r)
                ruv = -r/(s-alpha/3)*(ru*rv+exponential*(1-1/s-alpha/(3*s*s))/4)
                return (2*(r*ruv+ru*rv)/alpha, -ruv/r-exponential/(4*s*s))

            old = original(y0, yu0, yv0, e0)
            new = original(y0+mp.mpf(float(q[0])), yu0+mp.mpf(float(qu[0])),
                           yv0+mp.mpf(float(qv[0])), e0*mp.exp(2*alpha*mp.mpf(float(q[1]))))
            exact = np.array([float(new[0]-old[0]), float((new[1]-old[1])/alpha)])
            calculated = _finite_bulk(background, b, q, qu, qv)
            errors.append(float(np.max(abs(exact-calculated)/np.maximum(abs(exact), 1e-12))))
    assert max(errors) < 2e-12
    return dict(original_equation_precision_digits=70, sampled_nonzero_jets=4,
                maximum_component_relative_error=max(errors))


def _finite_edges(background, source, extent, width):
    def nonlinear_rhs(u, state):
        b = _finite_background(background, u/2)
        w, derivative = state
        sb, yu, sx = (float(b[k]) for k in ('S', 'Yu', 'sx'))
        s = sb+background.alpha*w
        assert s-background.k > 0
        square_difference = ((2*yu*derivative+derivative**2)/s
                             -background.alpha*w*yu**2/(s*sb))
        return derivative, (background.alpha*square_difference/2
                            +sx*derivative-2*float(source(u))/3)

    def linear_rhs(u, state):
        b = background.fields(u/2)
        r, rx, rxx, sx = (float(b[k]) for k in ('R', 'Rx', 'Rxx', 'sx'))
        value, derivative = state
        return derivative, (sx*derivative-(rxx/4-sx*rx/2)/r*value-float(source(u))/(3*r))

    solutions = []
    for rhs in (nonlinear_rhs, linear_rhs):
        sol = solve_ivp(rhs, (0, extent), (0., 0.), method='DOP853',
                        rtol=2e-12, atol=2e-13, max_step=width/120,
                        dense_output=True)
        assert sol.success
        solutions.append(sol)
    return solutions


def _finite_evolve(background, edge, extent, intervals):
    n, h = intervals, extent/intervals
    coordinate = np.linspace(0, extent, n+1)
    state = np.zeros((n+1, n+1, 2))
    state[:, 0, 0] = edge.sol(coordinate)[0]
    cached = _finite_background(background, np.arange(-n, n+1)*h/2)
    maximum_iterations, maximum_defect = 0, 0.
    for diagonal in range(2, 2*n+1):
        i = np.arange(max(1, diagonal-n), min(n, diagonal-1)+1)
        j = diagonal-i
        b = {key:value[i-j+n] for key,value in cached.items()}
        sw, se, nw = state[i-1, j-1], state[i, j-1], state[i-1, j]
        free = se+nw-sw
        ne = free.copy()
        for iteration in range(1, 13):
            middle = (sw+se+nw+ne)/4
            du, dv = (se+ne-sw-nw)/(2*h), (nw+ne-sw-se)/(2*h)
            updated = free+h*h*_finite_bulk(background, b, middle, du, dv)
            change = float(np.max(abs(updated-ne)/(1+abs(updated))))
            ne = updated
            if change < 3e-13:
                break
        else:
            raise AssertionError(('Nonlinear characteristic box did not converge', intervals, diagonal, change))
        middle = (sw+se+nw+ne)/4
        du, dv = (se+ne-sw-nw)/(2*h), (nw+ne-sw-se)/(2*h)
        defect = float(np.max(abs(ne-free-h*h*_finite_bulk(background, b, middle, du, dv))/(1+abs(ne))))
        maximum_defect = max(maximum_defect, defect)
        maximum_iterations = max(maximum_iterations, iteration)
        state[i, j] = ne
    assert np.all(np.isfinite(state)) and maximum_defect < 2e-12
    return coordinate, state, dict(maximum_iterations=maximum_iterations,
                                   maximum_normalized_box_defect=maximum_defect)


def _finite_constraints(background, state, coordinate, source):
    h = coordinate[1]-coordinate[0]
    f, value = state, state[2:-2, 2:-2]
    du = (f[:-4, 2:-2]-8*f[1:-3, 2:-2]+8*f[3:-1, 2:-2]-f[4:, 2:-2])/(12*h)
    dv = (f[2:-2, :-4]-8*f[2:-2, 1:-3]+8*f[2:-2, 3:-1]-f[2:-2, 4:])/(12*h)
    duu = (-f[:-4, 2:-2]+16*f[1:-3, 2:-2]-30*value+16*f[3:-1, 2:-2]-f[4:, 2:-2])/(12*h*h)
    dvv = (-f[2:-2, :-4]+16*f[2:-2, 1:-3]-30*value+16*f[2:-2, 3:-1]-f[2:-2, 4:])/(12*h*h)
    u, v = np.meshgrid(coordinate[2:-2], coordinate[2:-2], indexing='ij')
    b = _finite_background(background, (u-v)/2)
    alpha, sb = background.alpha, b['S']
    s = sb+alpha*value[..., 0]
    statesource = np.asarray(source(u))
    residuals = []
    for derivative, second, ybase, sbase, pulse in (
            (du, duu, b['Yu'], b['sx']/2, statesource),
            (dv, dvv, b['Yv'], -b['sx']/2, 0.)):
        w1, z1 = derivative[...,0], derivative[...,1]
        square_difference = ((2*ybase*w1+w1*w1)/s-alpha*value[...,0]*ybase*ybase/(s*sb))
        c = (second[...,0]/2-alpha*square_difference/4-sbase*w1
             -alpha*z1*(ybase+w1)
             +alpha*(second[...,1]-2*sbase*z1-alpha*z1*z1)/3+pulse/3)
        residuals.append(c)
    mask = (u>=.06)&(v>=.06)&(u<=coordinate[-1]-.06)&(v<=coordinate[-1]-.06)
    scale = 1+float(np.max(abs(statesource)))/3
    norms = [float(np.max(abs(c[mask]))) for c in residuals]
    return dict(maximum_uu_over_alpha=norms[0], maximum_vv_over_alpha=norms[1],
                source_normalization=scale, maximum_source_scaled=max(norms)/scale,
                fixed_physical_inset=.06)


def _local_return_clock(background, coordinate, conformal, quadrature_order=5):
    """Proper elapsed time on the SAME prescribed x=+.2 timelike curve.

    A null ray launched t=.4 reflects ideally at x=-.2,t=.8 and returns
    t=1.2.  These reference trajectories are not physical wormhole mouths.
    """
    from scipy.interpolate import RectBivariateSpline
    x, initial, final = .2, .4, 1.2
    assert initial-x >= coordinate[0] and final+x <= coordinate[-1]
    interpolate = RectBivariateSpline(coordinate, coordinate, conformal, kx=3, ky=3)
    cuts = np.unique(np.r_[initial, final, coordinate-x, coordinate+x])
    cuts = cuts[(cuts>=initial)&(cuts<=final)]
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    half = np.diff(cuts)/2
    middle = (cuts[:-1]+cuts[1:])/2
    times = middle[:,None]+half[:,None]*nodes
    z = interpolate.ev((times+x).ravel(), (times-x).ravel()).reshape(times.shape)
    lapse = math.sqrt(float(background.fields(x)['e2s']))
    delta = lapse*float(np.sum(half[:,None]*weights*np.expm1(background.alpha*z)/background.alpha))
    return dict(static_proper_duration=(final-initial)*lapse,
                proper_duration_change_over_alpha=delta,
                proper_duration_change=background.alpha*delta)


def _finite_case(alpha, omega2, source, extent, width, intervals):
    background = StaticBackground(alpha, omega2, extent/2)
    static_checks = _independent_operator_checks(background, omega2, extent)
    algebra_checks = _finite_bulk_algebra_check(background)
    edge, linear_edge = _finite_edges(background, source, extent, width)
    records, nonlinear_states, linear_states = [], [], []
    for n in intervals:
        coord, state, iteration = _finite_evolve(background, edge, extent, n)
        _, linear, _ = _evolve(background, linear_edge, extent, n)
        rbar = background.fields((coord[:,None]-coord[None,:])/2)['R']
        full_r = np.sqrt(rbar*rbar+alpha*state[...,0])
        normalized = np.stack((state[...,0]/(full_r+rbar), state[...,1]),axis=-1)
        remainder = float(np.max(abs(normalized-linear)))
        constraints = _finite_constraints(background, state, coord, source)
        clock = _local_return_clock(background, coord, state[...,1])
        linear_clock = _local_return_clock(background, coord, linear[...,1])
        records.append(dict(intervals=n, iteration=iteration, constraints=constraints,
            minimum_R_squared=float(np.min(full_r*full_r)), minimum_R_squared_minus_k=float(np.min(full_r*full_r-background.k)),
            maximum_normalized_radius_response=float(np.max(abs(normalized[...,0]))),
            maximum_normalized_conformal_response=float(np.max(abs(normalized[...,1]))),
            maximum_fractional_sphere_area_change=float(np.max(abs(alpha*state[...,0]/(rbar*rbar)))),
            normalized_geometry_nonlinear_remainder=remainder,
            local_clock=clock, linear_geometry_same_source_local_clock=linear_clock,
            final_corner_radius_change=float(alpha*normalized[-1,-1,0]),
            final_corner_conformal_change=float(alpha*normalized[-1,-1,1])))
        nonlinear_states.append(normalized)
        linear_states.append(linear)
    differences = [float(np.max(abs(a-b[::2,::2]))) for a,b in zip(nonlinear_states[:-1],nonlinear_states[1:])]
    linear_differences = [float(np.max(abs(a-b[::2,::2]))) for a,b in zip(linear_states[:-1],linear_states[1:])]
    remainder_differences = [float(np.max(abs((a-la)-(b-lb)[::2,::2])))
                            for a,b,la,lb in zip(nonlinear_states[:-1],nonlinear_states[1:],linear_states[:-1],linear_states[1:])]
    ratios = [a/b for a,b in zip(differences[:-1],differences[1:])]
    clocks = [row['local_clock']['proper_duration_change_over_alpha'] for row in records]
    clock_differences = [abs(a-b) for a,b in zip(clocks[:-1],clocks[1:])]
    fine_clock_check = _local_return_clock(background, coord, nonlinear_states[-1][...,1],7)
    quadrature_error = abs(fine_clock_check['proper_duration_change_over_alpha']-clocks[-1])
    constraints = [row['constraints']['maximum_source_scaled'] for row in records]
    assert all(3.5 < ratio < 4.5 for ratio in ratios)
    assert all(b < a for a,b in zip(constraints[:-1], constraints[1:])) and constraints[-1] < 1e-4
    assert differences[-1] < 1e-4 and linear_differences[-1] < 1e-4
    clock_ratios = [a/b for a,b in zip(clock_differences[:-1],clock_differences[1:])]
    assert all(3.5 < ratio < 4.5 for ratio in clock_ratios)
    assert clock_differences[-1] < 1e-5 and quadrature_error < 1e-10
    stride = intervals[-1]//32
    sample = nonlinear_states[-1][::stride,::stride]
    return dict(alpha=alpha, refinement=records, common_grid_state_differences=differences,
                linear_geometry_grid_differences=linear_differences,
                nonlinear_remainder_grid_differences=remainder_differences,
                observed_refinement_ratios=ratios, normalized_local_clock_grid_differences=clock_differences,
                local_clock_refinement_ratios=clock_ratios,
                local_clock_quadrature_order_5_vs_7_error=quadrature_error,
                independent_static_and_linear_operator_checks=static_checks,
                independent_nonlinear_bulk_algebra_checks=algebra_checks,
                plot=dict(u=coord[::stride].tolist(),v=coord[::stride].tolist(),
                    radius_response_over_alpha=sample[...,0].tolist(),
                    conformal_response_over_alpha=sample[...,1].tolist(),
                    sphere_area_fraction_over_alpha=(state[::stride,::stride,0]/rbar[::stride,::stride]**2).tolist(),
                    center_time=coord[::stride].tolist(),
                    center_radius_response=nonlinear_states[-1][::stride,::stride,0].diagonal().tolist(),
                    center_conformal_response=nonlinear_states[-1][::stride,::stride,1].diagonal().tolist(),
                    linear_geometry_center_radius_response=linear_states[-1][::stride,::stride,0].diagonal().tolist(),
                    linear_geometry_center_conformal_response=linear_states[-1][::stride,::stride,1].diagonal().tolist(),
                    center_area_fraction_over_alpha=(state[::stride,::stride,0]/rbar[::stride,::stride]**2).diagonal().tolist(),
                    linearized_center_area_fraction_over_alpha=(2*linear_states[-1][::stride,::stride,0]/rbar[::stride,::stride]).diagonal().tolist()))


def polyakov_finite_checks(intervals=(128,256,512), control_alphas=(.01,.005,.0025)):
    """Focused finite-state nonlinear benchmark; old linear results are untouched."""
    from scipy.integrate import quad
    from mmp_pulse import _exterior_quantum_seed
    alpha, delay, amplitude, start, extent = 1e-9, 1.2, .15752457643370363, .15, 1.5
    assert len(intervals) >= 3 and all(b == 2*a for a,b in zip(intervals[:-1],intervals[1:]))
    assert intervals[-1] % 32 == 0
    assert len(control_alphas) >= 3 and all(b == a/2 for a,b in zip(control_alphas[:-1],control_alphas[1:]))
    omega2 = (2*math.pi/(math.pi+delay))**2
    b0 = math.pi*(omega2-1)/6
    w0 = (3-math.sqrt(3))/6
    max_bp = math.exp(4-1/w0)*math.sqrt(1-4*w0)/w0**2
    minimum_orientation = 1-amplitude*max_bp
    assert minimum_orientation > 0

    def source(u):
        a = np.asarray(u)
        values = np.fromiter((1.5*_exterior_quantum_seed(float(x)-start,amplitude,delay)
                              for x in a.ravel()), dtype=float, count=a.size)
        return values.reshape(a.shape)

    initial_energy = quad(lambda u:_exterior_quantum_seed(u,amplitude,delay),0,delay,
                          epsabs=2e-12,epsrel=2e-12)[0]
    assert initial_energy > 0 and abs(initial_energy/b0-4)<2e-12

    def integrated_state_energy(u):
        bp,bpp,_ = (float(x) for x in _bump_derivatives(np.asarray(u/delay)))
        return (-omega2*amplitude**2*bp**2+(amplitude*bpp/delay/(1+amplitude*bp))**2)/4

    independent_energy = quad(integrated_state_energy,0,delay,epsabs=2e-12,epsrel=2e-12)[0]
    energy_error = abs(independent_energy-1.5*initial_energy)
    assert independent_energy > 0 and energy_error < 2e-11
    background = StaticBackground(alpha,omega2,extent/2)
    zero_edge, zero_linear = _finite_edges(background,lambda u:0.,extent,delay)
    _,zero,_ = _finite_evolve(background,zero_edge,extent,32)
    _,zero_l,_ = _evolve(background,zero_linear,extent,32)
    assert np.max(abs(zero)) == np.max(abs(zero_l)) == 0
    physical = _finite_case(alpha,omega2,source,extent,delay,intervals)
    controls = [_finite_case(a,omega2,source,extent,delay,intervals) for a in control_alphas]
    nonlinear_remainders = [case['refinement'][-1]['normalized_geometry_nonlinear_remainder'] for case in controls]
    halving_ratios = [a/b for a,b in zip(nonlinear_remainders[:-1],nonlinear_remainders[1:])]
    error_margins = [remainder/(case['common_grid_state_differences'][-1]+case['linear_geometry_grid_differences'][-1])
                     for remainder,case in zip(nonlinear_remainders,controls)]
    assert all(1.8 < ratio < 2.2 for ratio in halving_ratios) and min(error_margins) > 10
    for case in controls:
        errors = case['nonlinear_remainder_grid_differences']
        assert all(3.5 < a/b < 4.5 for a,b in zip(errors[:-1],errors[1:]))
    source_grid = np.linspace(0,extent,401)
    return dict(model='Finite positive-energy CFT state with nonlinear local Einstein-Maxwell-Polyakov geometry on a characteristic diamond',
        parameters=dict(alpha=alpha,exterior_delay_used_only_to_set_state=delay,amplitude=amplitude,
            pulse_u_interval=[start,start+delay],null_coordinate_interval=[0,extent],
            omega_squared=omega2,b0=b0,initial_energy_helper_normalization=initial_energy,
            initial_energy_over_b0=initial_energy/b0,minimum_finite_map_derivative=minimum_orientation,
            integrated_delta_t=1.5*initial_energy,independent_integrated_delta_t=independent_energy,
            integrated_state_energy_identity_absolute_error=energy_error),
        conventions=dict(metric='ds4^2=-exp(2sigma)du dv+R^2 dOmega^2; re=1; u=t+x,v=t-x; durations in re/c',
            quantum_stress='t2_uu=q/(12pi)*(sigma_uu-sigma_u^2+t_u); t2_vv analogously; t2_uv=-q*sigma_uv/(12pi); p_perp=0; k=Gq/(12pi)=alpha/3',
            finite_state='V(u)=u+amplitude*delay*B((u-start)/delay), B=exp(4-1/[z*(1-z)]) inside(0,1), zero otherwise in the first conformal period; degree-one periodic extension of period2pi/Omega outside this local diamond',
            energy_identity='Integral(delta_t)=[Integral((Vsecond/Vprime)^2)-Omega^2*Integral((Vprime-1)^2)]/4; full finite-state energy retained'),
        source='Exact t_u=-Omega^2/4+(3/2)*_exterior_quantum_seed(u-start,amplitude,delay); t_v=-Omega^2/4. Finite initial CFT-state energy retained; state amplitude is not linearized',
        equations=dict(variables='S=R^2=1+alpha*Y; sigma=sigmabar+alpha*Z; evolved variables (Y-Ybar,Z)',
            radial='Y_uv=-exp(2sigma)*(Y-1/(3S))/(2*(S-alpha/3))-alpha^2*Y_u*Y_v/(6*S*(S-alpha/3))',
            angular='sigma_uv=-alpha*Y_uv/(2S)+alpha^2*Y_u*Y_v/(4*S^2)-exp(2sigma)/(4*S^2)',
            edge='Z=0 on both initial null rays is coordinate gauge; on v=0 solve Y_uu=alpha*Y_u^2/(2S)+2*sigmabar_u*Y_u-(2/3)*(sigmabar_uu-sigmabar_u^2+t_u), with static corner value/slope; u=0 remains static',
            stability='Subtract static equations algebraically using log1p/expm1; no subtraction of nearly equal unscaled lapse/radius fields before division by alpha',
            comparison='Geometry-linearized response uses the SAME exact finite t_u perturbation and the same static alpha background; this is not the earlier infinitesimal-state benchmark'),
        zero_source_maximum_response=float(np.max(abs(zero))), physical_alpha=physical,
        larger_alpha_numerical_controls=controls,
        nonlinear_scaling=dict(control_alpha_halving_ratios=halving_ratios,
            control_remainder_over_sum_of_full_and_linear_mesh_differences=error_margins,
            interpretation='Larger-alpha nonlinear geometry remainders exceed both solution discretization differences and halve approximately with alpha. The tiny physical-alpha remainder is below the total solution mesh difference; its paired-grid difference is recorded but is not a certified continuum nonlinear correction'),
        acceptance=dict(all_pass=True,zero_source_exact=True,state_and_clock_second_order_ratios_between=[3.5,4.5],
            final_full_and_linear_common_grid_difference_below=1e-4,
            source_scaled_constraints_decrease=True,final_source_scaled_null_constraint_below=1e-4,
            final_normalized_clock_grid_difference_below=1e-5,clock_quadrature_difference_below=1e-10,
            control_nonlinear_remainder_halving_ratios_between=[1.8,2.2],
            control_remainder_exceeds_sum_of_mesh_differences_by=10,
            exact_bulk_high_precision_relative_tolerance=2e-12,
            independently_differentiated_static_tolerance=2e-7),
        plot=dict(source_u=source_grid.tolist(),delta_t_u=source(source_grid).tolist()),
        local_clock_scope='Prescribed accelerated conformal-coordinate curves x=+.2 and x=-.2; null launch t=.4, ideal reflection t=.8, return t=1.2. Proper elapsed time is integrated on the SAME x=+.2 emitter/receiver curve, in re/c units. Causal by construction; no material mouths, exterior clock identification or time-travel test',
        local_clock_linear_comparison='The exact proper-time functional is evaluated on the geometry-linearized conformal field; this is not truncation of the clock functional itself to first order',
        scope='The initial data ASSUME an existing calibrated throat; this does not model throat formation from ordinary space. Local prescribed finite CFT state and conditional radius-independent Polyakov action. Not the eta4 circulating history: coordinate gauge, incoming-data prescription and absent exterior feedback differ. Larger-alpha cases are numerical controls, not validated designs. No physical state preparation, full4D renormalized matter completion, global matching, stability or chronology claim',
        provenance=dict(static_helper_sha256=hashlib.sha256(inspect.getsource(_polyakov_static_solution).encode()).hexdigest(),
                        exact_seed_helper_sha256=hashlib.sha256(inspect.getsource(_exterior_quantum_seed).encode()).hexdigest()),
        sources=['https://arxiv.org/pdf/math-ph/0412028 (Fewster-Hollands conformal stress transformation, Eqs2.9 and3.41)',
                 'https://arxiv.org/pdf/hep-th/0006035 (Polyakov action and state-dependent stress)',
                 'https://arxiv.org/pdf/hep-th/0111091 (spherical Einstein-Maxwell reduction)',
                 'The normalized nonlinear diamond and same-source comparison are our application of this conditional model'])


def save_finite_benchmark():
    result = polyakov_finite_checks()
    project_path('polyakov_finite_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_FINITE_SOURCE_RESPONSE


# BEGIN_RADIUS_TUBE_DIAGNOSTIC

def _radius_tube_field(background, coordinate, state):
    from scipy.interpolate import RectBivariateSpline
    y, z = (RectBivariateSpline(coordinate, coordinate, state[..., i]) for i in (0, 1))
    # FITPACK otherwise rebuilds whole-grid derivative coefficients at every scalar query.
    derivatives = {f: {(du,dv): f.partial_derivative(du,dv) if du+dv else f
        for du,dv in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))} for f in (y,z)}
    for f, cached in derivatives.items():
        for (du,dv), derivative in cached.items():
            assert np.isclose(float(derivative(.35,.55,grid=False)),
                float(f.ev(.35,.55,dx=du,dy=dv)),rtol=1e-12,atol=1e-12)

    def field(t, x, full=False):
        u, v = t+x, t-x
        assert 0 <= u <= coordinate[-1] and 0 <= v <= coordinate[-1], 'Outside solved diamond'
        b = _finite_background(background, x)
        a = background.alpha
        ev = lambda f, du=0, dv=0: float(derivatives[f][du,dv](u, v, grid=False))
        result = dict(Y=float(b['Y'])+ev(y), Yu=float(b['Yu'])+ev(y, 1),
            Yv=float(b['Yv'])+ev(y, 0, 1), sigma=float(b['sigma'])+a*ev(z))
        if full:
            result.update(Yuu=-float(b['Yuv'])+ev(y, 2),
            Yuv=float(b['Yuv'])+ev(y, 1, 1), Yvv=-float(b['Yuv'])+ev(y, 0, 2),
            su=float(b['sx'])/2+a*ev(z, 1),
            sv=-float(b['sx'])/2+a*ev(z, 0, 1), suu=float(b['sxx'])/4+a*ev(z, 2),
            suv=-float(b['sxx'])/4+a*ev(z, 1, 1), svv=float(b['sxx'])/4+a*ev(z, 0, 2))
        return result
    return field


def _trace_radius_tube(field, side, cutoff=.05, start=.075, finish=1.4, refinement=1):
    def rhs(t, x):
        f = field(t, float(x[0]))
        return [-(f['Yu']+f['Yv'])/(f['Yu']-f['Yv'])]

    def null(t, x):
        f = field(t, float(x[0]))
        return -4*f['Yu']*f['Yv']

    def chart(t, x):
        f = field(t, float(x[0]))
        return side*(f['Yu']-f['Yv'])

    for event in (null, chart):
        event.terminal, event.direction = True, -1
    sol = solve_ivp(rhs, (start, finish), [side*cutoff], method='DOP853',
        rtol=3e-10/refinement**2, atol=3e-12/refinement**2,
        max_step=.002/refinement, dense_output=True, events=(null, chart))
    assert sol.success
    stop = 'null_level_set' if len(sol.t_events[0]) else ('radius_chart_degenerate' if len(sol.t_events[1]) else 'test_interval_end')
    return sol, stop


def _tube_junction_values(field, source, alpha, omega2, side, t, x):
    f = field(t, x, full=True)
    yt, yx = f['Yu']+f['Yv'], f['Yu']-f['Yv']
    velocity = -yt/yx
    gap = 1-velocity**2
    assert gap > 0, 'Junction frame requires a timelike tube'
    gamma = 1/math.sqrt(gap)
    sigma = f['sigma']; s = 1+alpha*f['Y']; r = math.sqrt(s)
    ytt, ytx, yxx = f['Yuu']+2*f['Yuv']+f['Yvv'], f['Yuu']-f['Yvv'], f['Yuu']-2*f['Yuv']+f['Yvv']
    acceleration = -(ytt+2*velocity*ytx+velocity**2*yxx)/yx
    ktau = side*math.exp(-sigma)*gamma*(gamma**2*acceleration+f['su']-f['sv']+velocity*(f['su']+f['sv']))
    cu = f['suu']-f['su']**2-omega2/4+float(source(t+x))
    cv = f['svv']-f['sv']**2-omega2/4
    cuv = -f['suv']
    up, down = (1+velocity)/(1-velocity), (1-velocity)/(1+velocity)
    stress_scale = 2*math.exp(-2*sigma)/(3*s)
    return dict(time=t, x=x, velocity=velocity, proper_lapse=math.exp(sigma)/gamma,
        Ktheta_over_alpha=side*math.exp(-sigma)*math.sqrt(gap)*yx/(2*s),
        Ktau=ktau,
        mass_excess_over_alpha_squared=(4*f['Y']**2/(r+1)**2-math.exp(-2*sigma)*yx*yx*gap)/(8*r),
        quantum_normal_pressure_over_alpha=stress_scale*(cu*up+cv*down-2*cuv),
        quantum_TnU_over_alpha=side*stress_scale*(cu*up-cv*down))


def _radius_tube_case(background, source, omega2, extent, intervals, refinement=1):
    from scipy.integrate import quad
    from scipy.optimize import brentq
    edge, _ = _finite_edges(background, source, extent, 1.2)
    coordinate, state, _ = _finite_evolve(background, edge, extent, intervals)
    field = _radius_tube_field(background, coordinate, state)
    target = float(background.fields(.05)['Y'])
    solutions, records = {}, {}
    common_times = np.linspace(.075, .22, 41)
    for side in (-1, 1):
        sol, stop = _trace_radius_tube(field, side, refinement=refinement)
        solutions[side] = sol
        final, xfinal = float(sol.t[-1]), float(sol.y[0,-1])
        times = np.linspace(.075, final, 201)
        positions = sol.sol(times)[0]
        jets = [field(float(t), float(x)) for t,x in zip(times, positions)]
        speeds = [-(f['Yu']+f['Yv'])/(f['Yu']-f['Yv']) for f in jets]
        radius_error = max(abs(f['Y']-target) for f in jets)
        assert radius_error < 2e-8 and max(abs(np.asarray(speeds[:-1]))) < 1
        assert final > common_times[-1]
        surface = [_tube_junction_values(field, source, background.alpha, omega2, side,
                    float(t), float(sol.sol(t)[0])) for t in common_times]
        end = jets[-1]
        # These are derivatives along increasing t at fixed x, not beyond-null tube continuation.
        end_slice = [field(final+dt, xfinal) for dt in (-.0001, .0001)]
        records[str(side)] = dict(stop=stop, final_time=final, final_x=xfinal,
            final_Yu=end['Yu'], final_Yv=end['Yv'], final_velocity=speeds[-1],
            normalized_radius_constraint_maximum=radius_error,
            nearby_fixed_x_expansion_signs=[[int(np.sign(f['Yu'])),int(np.sign(f['Yv']))] for f in end_slice],
            boundary_data=surface,
            plot=dict(time=times.tolist(),x=positions.tolist(),velocity=speeds,
                proper_lapse=[math.exp(f['sigma'])*math.sqrt(max(0.,1-v*v)) for f,v in zip(jets,speeds)]))

    right, left = solutions[1], solutions[-1]
    launch = .08
    ulaunch = launch+float(right.sol(launch)[0])
    reflection = brentq(lambda t:t+float(left.sol(t)[0])-ulaunch, launch, left.t[-1],xtol=2e-13)
    vreturn = reflection-float(left.sol(reflection)[0])
    if right.t[-1]-float(right.sol(right.t[-1])[0]) <= vreturn:
        probe = dict(launch_time=launch, reflection_time=reflection, return_time=None,
            status='No return intersection before selected right tube ceases to be timelike; censored')
    else:
        arrival = brentq(lambda t:t-float(right.sol(t)[0])-vreturn, reflection, right.t[-1],xtol=2e-13)
        def lapse(t):
            x = float(right.sol(t)[0]); f=field(t,x)
            speed=-(f['Yu']+f['Yv'])/(f['Yu']-f['Yv'])
            assert abs(speed) < 1
            return math.exp(f['sigma'])*math.sqrt(1-speed*speed)
        duration, quad_error = quad(lapse, launch, arrival, epsabs=2e-10/refinement**2,
            epsrel=2e-10/refinement**2, points=right.t[(right.t>launch)&(right.t<arrival)],limit=1000)
        static_duration=.2*math.exp(float(background.fields(.05)['sigma']))
        probe=dict(launch_time=launch,launch_x=float(right.sol(launch)[0]),
            reflection_time=reflection,reflection_x=float(left.sol(reflection)[0]),
            return_time=arrival,return_x=float(right.sol(arrival)[0]),
            proper_round_trip_duration=duration, quadrature_error_estimate=quad_error,
            static_same_cutoff_round_trip_duration=static_duration,
            proper_duration_difference=duration-static_duration,
            return_before_right_null_gate=right.t[-1]-arrival,
            status='Completed on surviving timelike segments; ideal test-ray reflection only')
        assert duration > 0 and arrival > reflection > launch
    return dict(intervals=intervals, tube_integration_refinement=refinement,
        target_Y=target, sides=records, probe=probe)


def _boundary_verification(cases, zero, tight, background):
    keys=('x','velocity','proper_lapse','Ktheta_over_alpha','Ktau',
          'mass_excess_over_alpha_squared','quantum_normal_pressure_over_alpha','quantum_TnU_over_alpha')
    def differences(a,b):
        surface={key:max(abs(x[key]-y[key])/(1+abs(y[key]))
            for side in ('-1','1') for x,y in zip(a['sides'][side]['boundary_data'],b['sides'][side]['boundary_data'])) for key in keys}
        return dict(gate_time=max(abs(a['sides'][side]['final_time']-b['sides'][side]['final_time']) for side in ('-1','1')),
            return_time=abs(a['probe']['return_time']-b['probe']['return_time']),
            proper_duration=abs(a['probe']['proper_round_trip_duration']-b['probe']['proper_round_trip_duration']),
            boundary_component_scaled=surface)
    refinement=[differences(a,b) for a,b in zip(cases[:-1],cases[1:])]
    for field in ('gate_time','return_time','proper_duration'):
        assert all(b[field]<a[field] for a,b in zip(refinement[:-1],refinement[1:]))
    assert refinement[-1]['gate_time']<2e-6 and refinement[-1]['proper_duration']<2e-7
    assert max(refinement[-1]['boundary_component_scaled'].values())<.01
    integration=differences(cases[-1],tight)
    assert integration['gate_time']<2e-8 and integration['proper_duration']<2e-9
    assert max(integration['boundary_component_scaled'].values())<3e-6
    for case in cases+[tight]:
        for row in case['sides'].values():
            assert row['stop']=='null_level_set'
            assert min(abs(row['final_Yu']),abs(row['final_Yv']))<1e-10
            assert max(row['final_Yu'],row['final_Yv'])>.001
            assert abs(abs(row['final_velocity'])-1)<1e-8
    assert abs(zero['probe']['proper_duration_difference'])<1e-11
    assert abs(zero['probe']['return_time']-.28)<1e-11
    static_errors=[]; mass_errors=[]
    for side in ('-1','1'):
        row=zero['sides'][side]
        assert row['stop']=='test_interval_end'
        assert max(abs(np.asarray(row['plot']['x'])-int(side)*.05))<1e-11
        b=background.fields(int(side)*.05); s=float(b['R'])**2
        pressure=-(float(b['sx'])**2+(2*math.pi/(math.pi+1.2))**2)/(3*s*float(b['e2s']))
        ktau=int(side)*float(b['sx'])/math.sqrt(float(b['e2s']))
        for value in row['boundary_data']:
            static_errors.extend((abs(value['quantum_normal_pressure_over_alpha']-pressure),
                abs(value['quantum_TnU_over_alpha']),abs(value['Ktau']-ktau)))
    assert max(static_errors)<1e-11
    for case in cases+[zero,tight]:
        y=case['target_Y']; r=math.sqrt(1+background.alpha*y)
        for side in case['sides'].values():
            for value in side['boundary_data']:
                independent=y*y/(2*r*(r+1)**2)-r**3*value['Ktheta_over_alpha']**2/2
                mass_errors.append(abs(value['mass_excess_over_alpha_squared']-independent))
    assert max(mass_errors)<3e-9
    return dict(all_declared_checks_pass=True,mesh_differences=refinement,
        tighter_tube_integration_difference=integration,
        static_projection_and_acceleration_maximum_error=max(static_errors),
        charged_mass_extrinsic_curvature_identity_maximum_error=max(mass_errors),
        thresholds=dict(final_gate_time_difference=2e-6,final_proper_duration_difference=2e-7,
            final_boundary_component_difference_scaled_by_one_plus_abs=.01,
            tighter_integration_gate_time_difference=2e-8,tighter_integration_duration_difference=2e-9,
            tighter_integration_boundary_scaled_difference=3e-6),
        scope='Refinement comparisons, not certified continuum bounds. Mixed temporal curvature is the least accurate exported component. Boundary data are checked only on the fixed .075<=t<=.22 interval, away from the null endpoints.')


def polyakov_boundary_checks(intervals=(128,256,512,1024)):
    from mmp_pulse import _exterior_quantum_seed
    alpha, extent = 1e-9,1.5
    omega2=(2*math.pi/(math.pi+1.2))**2
    background=StaticBackground(alpha,omega2,extent/2)
    source=lambda u:1.5*_exterior_quantum_seed(float(u)-.15,.15752457643370363,1.2)
    cases=[_radius_tube_case(background,source,omega2,extent,n) for n in intervals]
    control=_radius_tube_case(background,lambda u:0.,omega2,extent,64)
    tight=_radius_tube_case(background,source,omega2,extent,intervals[-1],refinement=2)
    result=dict(model='Constant-areal-radius diagnostic tubes in the finite-source Einstein-Polyakov diamond; not material mouths or exterior matching',
        parameters=dict(alpha=alpha,extent=extent,cutoff_x=.05,initial_time=.075,
            target_Y=float(background.fields(.05)['Y']),launch_time=.08,
            source_amplitude=.15752457643370363,pulse_u_interval=[.15,1.35]),
        refinement=cases,zero_source_control=control,tight_tube_integration=tight,
        conventions=dict(velocity='dx/dt=-(Yu+Yv)/(Yu-Yv) on constant Y',
            timelike='Yu*Yv<0; null stop is a marginal sphere when exactly one derivative vanishes',
            proper_time='dproper/dt=exp(sigma)*sqrt(1-velocity^2); physical unit re/c',
            boundary_orientation='Normal orientation is initially outward: +x right, -x left; retained along each continuous branch',
            curvatures='Stored mixed Ktheta divided by alpha and mixed Ktau, both in inverse re units',
            mass='Charged mass aspect m: 1-2m/R+1/R²=(grad R)²; saved (m-1)/alpha² avoids catastrophic cancellation',
            stress='Quantum 8piG re² T_nn/alpha and T_nU/alpha; outward energy flux is minus T_nU; restore c^-4 for SI stress'),
        scope='Initial geometry assumes an existing throat. These are tiny local radius level sets initialized before the pulse, not entrance devices. The ideal reflected null probe has negligible stress. Finite-diamond marginal spheres do not establish an event horizon, complete trapping horizon, physical mirror, global stability, creation or backward travel. Nonzero boundary quantum stress cannot be silently removed to paste on RN vacuum.',
        sources=['https://arxiv.org/pdf/gr-qc/9303006 (Hayward SectionIII: marginal spheres and distinction from event horizons)',
                 'https://arxiv.org/pdf/gr-qc/0310107 (Visser-Wiltshire Eqs7-8 and22-24: junction and extrinsic-curvature conventions)',
                 'The moving level sets, local probe and normalized matching diagnostics are our conditional calculation'])
    result['verification']=_boundary_verification(cases,control,tight,background)
    return result


def save_boundary_benchmark():
    result=polyakov_boundary_checks()
    project_path('polyakov_boundary_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_RADIUS_TUBE_DIAGNOSTIC


# BEGIN_SOURCE_RETAINING_JUNCTION
def _single_null_candidates(cu, cv, cm, scale, side):
    """Necessary normal-stress matching for signed single-stream Vaidya."""
    candidates=[]
    for branch,ratio in ((1,cm/cu if cu else -1),(2,cv/cm if cm else -1)):
        if ratio>0:
            velocity=(ratio-1)/(ratio+1)
            p=scale*(ratio*cu+cv/ratio-2*cm); j=side*scale*(ratio*cu-cv/ratio)
            candidates.append(dict(branch=branch,velocity=velocity,pressure=p,nU=j,
                positive_energy=p>=0,residual=abs(p+side*j) if branch==1 else abs(p-side*j),
                invariant_pressure_error=abs(p-scale*(cu*cv-cm*cm)/cm)))
    return candidates


def _source_retaining_jet(field, source, alpha, omega2, side, t, x):
    """Boundary first derivatives in an areal/null chart, not a global exterior.

    Q=alpha*q, z=(R-1)/alpha keep the near-extremal chart well conditioned.
    The retarded q is normalized to psi=0 along the fixed-x boundary.
    """
    f=field(t,x,full=True)
    y,yu,yv=(f[k] for k in ('Y','Yu','Yv'))
    suu,suv,svv=(f[k] for k in ('suu','suv','svv'))
    su,sv=f['su'],f['sv']; sigma=f['sigma']
    yuu,yuv,yvv=(f[k] for k in ('Yuu','Yuv','Yvv'))
    s=1+alpha*y; r=math.sqrt(s); em=math.exp(-sigma)
    yt,yx=yu+yv,yu-yv; st,sx=su+sv,su-sv
    ytt,ytx=yuu+2*yuv+yvv,yuu-yvv
    d=em*yt/(2*r); b=side*em*yx/(2*r)
    assert b>0 and b+d>0, 'Chosen outward areal/null chart is not admissible'
    ddot=em**2*((ytt-st*yt)/(2*r)-alpha*yt*yt/(4*r**3))
    bdot=side*em**2*((ytx-st*yx)/(2*r)-alpha*yt*yx/(4*r**3))
    qdot=1/(b+d); qddot=-(bdot+ddot)/(b+d)**2
    h=y*y/(2*r*(r+1)**2)+em**2*yu*yv/(2*r)
    hu=y*yu/(4*r**3)+em**2*((yuu*yv+yu*yuv-2*su*yu*yv)/(2*r)-alpha*yu*yu*yv/(4*r**3))
    hv=y*yv/(4*r**3)+em**2*((yvv*yu+yv*yuv-2*sv*yu*yv)/(2*r)-alpha*yv*yv*yu/(4*r**3))
    yd,yo,ydd,sd,hd,ho=(yu,yv,yuu,su,hu,hv) if side==1 else (yv,yu,yvv,sv,hv,hu)
    mr=2*r*hd/yd                 # M_R/alpha
    mq=em**2*(yd*ho-yo*hd)/r     # M_q/alpha^3
    psir=2*r*(2*sd/yd-ydd/yd**2+alpha/(2*s))  # alpha*psi_R
    psiq=-psir*d/qdot            # psi_q/alpha, from psi|boundary=0
    fh=-em**2*yu*yv/s            # F/alpha^2
    fr=2*y/(r**3*(r+1))+2*alpha*h/s-2*mr/r
    fq=-2*mq/r
    metric=np.array([[-fh,-1.],[-1.,0.]])
    inverse=np.linalg.inv(metric)
    # Derivatives are ordered [coordinate, metric row, metric column].
    dg=np.array([[[-fq-2*psiq*fh,-psiq],[-psiq,0.]],
                 [[-fr-2*psir*fh,-psir],[-psir,0.]]])
    connection=np.zeros((2,2,2))
    for a in range(2):
        for i in range(2):
            for j in range(2):
                connection[a,i,j]=sum(inverse[a,k]*(dg[i,k,j]+dg[j,k,i]-dg[k,i,j])/2 for k in range(2))
    tangent=np.array([qdot,d]); normal=np.array([-qdot,b])
    acceleration=np.array([qddot,ddot])+np.einsum('aij,i,j->a',connection,tangent,tangent)
    ktau=float(normal@metric@acceleration)
    # Spherical warped-product Einstein tensor, after subtracting the same
    # monopole Maxwell stress. Algebraically subtract the O(1) background.
    trace=float(np.sum(inverse*connection[1]))
    einstein=2*connection[1]/r+metric*(-2*trace/r+alpha*fh/s-y/s**2)
    projected=dict(energy=float(tangent@einstein@tangent),
        pressure=float(normal@einstein@normal),nU=float(normal@einstein@tangent))
    cu=suu-su*su-omega2/4+float(source(t+x)); cv=svv-sv*sv-omega2/4; cm=-suv
    scale=2*em**2/(3*s)
    expected=dict(energy=scale*(cu+cv+2*cm),pressure=scale*(cu+cv-2*cm),nU=side*scale*(cu-cv))
    delta=expected['energy']+expected['pressure']+2*expected['nU']
    source_psir=r*delta/(2*(b+d)**2)
    source_mr=s/2*(-expected['pressure']-expected['nU']+b*delta/(b+d))
    source_mq=s/2*((b+d)**2*expected['nU']-b*d*delta)
    return dict(time=t,x=x,R_minus_one_over_alpha=y/(r+1),
        radius_proper_velocity_over_alpha=d,outward_radius_gradient_over_alpha=b,
        scaled_retarded_clock_rate=math.exp(sigma)*qdot,
        mass_excess_over_alpha_squared=h,mass_radial_derivative_over_alpha=mr,
        mass_null_derivative_over_alpha_cubed=mq,
        lapse_radial_derivative_times_alpha=psir,lapse_null_derivative_over_alpha=psiq,
        lapse_radial_derivative_from_quantum_source=source_psir,
        mass_radial_derivative_from_quantum_source=source_mr,
        mass_null_derivative_from_quantum_source=source_mq,
        quantum_components=dict(Cuu=cu,Cvv=cv,Cuv=cm),
        quantum_stress=expected,exterior_radial_einstein_stress=projected,
        missing_null_stress_if_radial_lapse_is_removed=expected['energy']+expected['pressure']+2*expected['nU'],
        single_null_candidates=_single_null_candidates(cu,cv,cm,scale,side),
        checks=dict(induced_metric=abs(float(tangent@metric@tangent)+1),
            normal_unit=abs(float(normal@metric@normal)-1),orthogonality=abs(float(normal@metric@tangent)),
            Ktau_difference=abs(ktau-side*em*sx),
            Ktheta_scaled_difference=abs((fh*qdot+d)/r-b/r),
            charged_mass_identity=abs(fh-(b*b-d*d)),
            mass_chain_rule=abs(hu+hv-(mq*math.exp(sigma)*qdot+mr*yt/(2*r))),
            quantum_einstein_scaled=max(abs(projected[k]-expected[k])/(1+abs(expected[k])) for k in expected),
            lapse_source_scaled=abs(psir-source_psir)/(1+abs(source_psir)),
            mass_radial_source_scaled=abs(mr-source_mr)/(1+abs(source_mr)),
            mass_null_source_scaled=abs(mq-source_mq)/(1+abs(source_mq))))


def polyakov_junction_checks(intervals=(128,256,512,1024)):
    from mmp_pulse import _exterior_quantum_seed
    alpha,extent=1e-9,1.5; omega2=(2*math.pi/(math.pi+1.2))**2
    background=StaticBackground(alpha,omega2,extent/2)
    source=lambda u:1.5*_exterior_quantum_seed(float(u)-.15,.15752457643370363,1.2)
    times=np.linspace(.075,.22,41)
    def run(n,forcing):
        edge,_=_finite_edges(background,forcing,extent,1.2)
        coordinate,state,_=_finite_evolve(background,edge,extent,n)
        field=_radius_tube_field(background,coordinate,state)
        rows={str(side):[_source_retaining_jet(field,forcing,alpha,omega2,side,float(t),side*.05)
                        for t in times] for side in (-1,1)}
        maxima={key:max(row['checks'][key] for values in rows.values() for row in values)
                for key in rows['1'][0]['checks']}
        scan=[]
        for side in (-1,1):
            for t in np.linspace(.075,1.4,161):
                x=side*.05; f=field(float(t),x,full=True)
                cu=f['suu']-f['su']**2-omega2/4+float(forcing(t+x))
                cv=f['svv']-f['sv']**2-omega2/4; cm=-f['suv']
                scale=2*math.exp(-2*f['sigma'])/(3*(1+alpha*f['Y']))
                scan.append(dict(time=float(t),side=side,Cuu=cu,Cvv=cv,Cuv=cm,
                    candidates=_single_null_candidates(cu,cv,cm,scale,side)))
        return dict(intervals=n,sides=rows,maximum_checks=maxima,single_stream_scan=scan)
    cases=[run(n,source) for n in intervals]; zero=run(64,lambda u:0.)
    differences=[]
    keys=('mass_excess_over_alpha_squared','mass_radial_derivative_over_alpha',
          'mass_null_derivative_over_alpha_cubed','lapse_radial_derivative_times_alpha',
          'lapse_null_derivative_over_alpha','scaled_retarded_clock_rate')
    for coarse,fine in zip(cases[:-1],cases[1:]):
        differences.append({key:max(abs(a[key]-b[key])/(1+abs(b[key])) for side in ('-1','1')
            for a,b in zip(coarse['sides'][side],fine['sides'][side])) for key in keys})
    result=dict(model='Source-retaining first derivatives at fixed-coordinate timelike boundary; not a global exterior solution',
        parameters=dict(alpha=alpha,omega_squared=omega2,extent=extent,x_cutoff=.05,time_interval=[.075,.22]),
        refinement=cases,zero_source_control=zero,mesh_scaled_differences=differences,
        conventions=dict(exterior_metric='-exp(2psi)F dq²-2exp(psi)dq dR+R²dOmega²; F=1-2M/R+1/R²',
            chart='Retarded q reparametrizes v on right or u on left. psi=0 along each boundary. Q=alpha*q,z=(R-1)/alpha are numerically scaled coordinates.',
            stress='Quantum stress normalized as 8piG re² T/alpha (c=1), common same-charge monopole Maxwell source subtracted. Outward energy flux=-T_nU.',
            derivatives='First partial derivatives in areal/null chart, all evaluated at the boundary. Retarded clock rate is d(alpha*q)/dt.'),
        scope='Existing throat and prescribed finite CFT state assumed. No actual material boundary, angular Einstein equation for an extended atmosphere, external CFT state, asymptotic solution, formation or stability is established. Matching first derivatives is boundary data, not a constructed medium or chronology change.',
        sources=['https://arxiv.org/pdf/gr-qc/0206076 (charged Vaidya and normal junction conditions)',
                 'https://arxiv.org/pdf/1301.1417 (generalized Vaidya mass derivatives and nonzero exterior stress)',
                 'Our coordinate transformation, first-derivative extraction and numerical Einstein checks'])
    geometric=('induced_metric','normal_unit','orthogonality','Ktau_difference',
        'Ktheta_scaled_difference','charged_mass_identity','mass_chain_rule')
    sourced=('quantum_einstein_scaled','lapse_source_scaled','mass_radial_source_scaled','mass_null_source_scaled')
    for case in cases+[zero]:
        assert max(case['maximum_checks'][k] for k in geometric)<1e-10
        scan=case['single_stream_scan']
        assert max(row['Cvv'] for row in scan)<-.27 and min(row['Cuv'] for row in scan)>.25
        for row in scan:
            assert all(not candidate['positive_energy'] for candidate in row['candidates'])
            assert all(max(candidate['residual'],candidate['invariant_pressure_error'])<1e-10 for candidate in row['candidates'])
    assert max(zero['maximum_checks'][k] for k in sourced)<1e-11
    for key in sourced:
        assert all(b['maximum_checks'][key]<a['maximum_checks'][key] for a,b in zip(cases[:-1],cases[1:]))
        assert cases[-1]['maximum_checks'][key]<.001
    assert max(differences[-1].values())<.005
    result['verification']=dict(all_declared_checks_pass=True,
        geometric_maximum=max(case['maximum_checks'][k] for case in cases+[zero] for k in geometric),
        source_maximum=cases[-1]['maximum_checks'],
        maximum_last_mesh_scaled_difference=max(differences[-1].values()),
        thresholds=dict(geometric=1e-10,zero_source=1e-11,final_source_scaled=.001,last_mesh_scaled=.005),
        single_stream_scope='322 sampled events per grid on fixed x=+/-.05,.075<=t<=1.4. Algebra excludes positive-energy single-null matching for every finite radial boost wherever Cvv<0<Cuv. Sampling is not a continuum sign bound or full junction sufficiency test.',
        scope='Geometric identities check the coordinate transformation; decreasing independent quantum-Einstein residuals test source consistency. Neither is a global exterior or continuum error bound.')
    return result


def save_junction_benchmark():
    result=polyakov_junction_checks()
    project_path('polyakov_junction_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_SOURCE_RETAINING_JUNCTION


# BEGIN_STATIC_QUANTUM_COLLAR
def _static_collar_fields(state, alpha, omega2):
    """Cancellation-safe static constraint and angular equation, re=1."""
    y,h,n,x=state
    r=1+alpha*y
    beta2=y*y/(r*r)-2*h/r
    assert beta2>0 and r*r>alpha/3
    beta=math.sqrt(beta2)
    bterm=y/r**2+alpha*h/r-omega2*math.exp(-2*n)/6
    discriminant=r*r*beta2+2*bterm/3
    assert discriminant>0
    acceleration=2*bterm/(r*beta+math.sqrt(discriminant))
    acceleration_prime=(1/r**4-alpha/(r*r-alpha/3)*
        (2*y/r**2+2*alpha*h/r-1/(3*r**4))-acceleration**2)
    return r,beta,acceleration,acceleration_prime


def _static_collar_solution(alpha, omega2, refinement=1, outer_radius=100.):
    background=StaticBackground(alpha,omega2,1.3)
    start_x=.5; f=background.fields(start_x); r=float(f['R'])
    y=float(f['Y'])/(r+1); beta=float(f['Rx'])/(alpha*math.sqrt(float(f['e2s'])))
    h=y*y/(2*r)-r*beta*beta/2
    initial=np.array([y,h,float(f['sigma']),start_x])
    def rhs(distance,state):
        r,beta,a,ap=_static_collar_fields(state,alpha,omega2)
        return beta,beta/6*(2*ap+a*a-omega2*math.exp(-2*state[2])),a,math.exp(-state[2])
    def finish(distance,state): return 1+alpha*state[0]-outer_radius
    finish.terminal,finish.direction=True,1
    sol=solve_ivp(rhs,(0,outer_radius+80),initial,method='DOP853',
        rtol=1e-11/refinement**2,atol=1e-12/refinement**2,
        max_step=.1/refinement,dense_output=True,events=finish)
    assert sol.success and len(sol.t_events[0])==1, 'Collar did not reach declared radius'
    return sol,background,initial


def _static_collar_record(alpha, omega2, refinement=1):
    from scipy.optimize import brentq
    sol,background,initial=_static_collar_solution(alpha,omega2,refinement)
    def row(radius):
        endpoint_radius=1+alpha*sol.y[0,-1]
        if abs(radius-endpoint_radius)<2e-12:
            distance=float(sol.t[-1])
        else:
            distance=brentq(lambda ell:1+alpha*sol.sol(ell)[0]-radius,0,sol.t[-1],xtol=2e-13)
        state=sol.sol(distance);y,h,n,x=state
        r,beta,a,ap=_static_collar_fields(state,alpha,omega2)
        b=alpha*beta; lapse=math.exp(n)
        energy_over_alpha=(2*ap+a*a-omega2/lapse**2)/(3*r*r)
        pressure_over_alpha=-(a*a+omega2/lapse**2)/(3*r*r)
        # The angular equation was used in ap; this is a constraint check,
        # not an independently differentiated angular residual.
        constraint=beta*(alpha*beta+2*r*a)+a*a/3+omega2/(3*lapse**2)-y*(2+alpha*y)/(r*r)
        return dict(radius=r,proper_distance_from_start=distance,optical_x=x,
            lapse_relative_to_center=lapse,normal_radius_derivative=b,
            proper_log_lapse_derivative=a,mass_excess=alpha*alpha*h,
            mixed_angular_curvature=b/r,mixed_temporal_curvature=a,
            mass_excess_over_alpha=alpha*h,quantum_energy_over_alpha=energy_over_alpha,
            quantum_radial_pressure_over_alpha=pressure_over_alpha,
            scaled_lapse_ratio=alpha*lapse/(1-1/r),
            remaining_optical_path=2*math.pi/math.sqrt(omega2)-2*x,
            localization_coefficient_over_sqrt_flux=2*math.sqrt(2)*abs(b),
            constraint_scaled=abs(constraint)/(1+abs(y*(2+alpha*y)/(r*r))))
    checkpoints=[row(r) for r in (1.01,1.1,1.25,2.,10.,30.,100.)]
    overlap=[]
    for x in (.6,.9,1.2):
        ell=brentq(lambda ell:sol.sol(ell)[3]-x,0,sol.t[-1],xtol=2e-13)
        y,h,n,_=sol.sol(ell); r,beta,a,ap=_static_collar_fields(sol.sol(ell),alpha,omega2)
        f=background.fields(x); expected_y=float(f['Y'])/(float(f['R'])+1)
        expected_beta=float(f['Rx'])/(alpha*math.sqrt(float(f['e2s'])))
        overlap.append(max(abs(y-expected_y)/(1+abs(expected_y)),
            abs(beta-expected_beta)/(1+abs(expected_beta)),abs(n-float(f['sigma'])),
            abs(a-float(f['sx'])/math.sqrt(float(f['e2s'])))))
    plot=[row(float(r)) for r in 1+np.geomspace(.01,99,161)]
    leading=lambda r:-5/(18*r**3)+1/(3*r**2)-1/18
    result=dict(alpha=alpha,integration_refinement=refinement,
        proper_distance_to_outer_radius=float(sol.t[-1]),checkpoints=checkpoints,plot=plot,
        maximum_overlap_difference=max(overlap),
        maximum_constraint_scaled=max(p['constraint_scaled'] for p in plot),
        maximum_fixed_radius_mass_coefficient_difference=max(abs(p['mass_excess_over_alpha']-leading(p['radius'])) for p in checkpoints),
        maximum_scaled_lapse_difference=max(abs(p['scaled_lapse_ratio']-12/(math.pi*(omega2-1))) for p in checkpoints))
    if refinement==2:
        from numpy.polynomial import chebyshev as cheb
        derivative_checks=[]
        for degree in (12,20,28):
            errors=[]
            for lo,hi in ((.2,2.2),(8.,10.),(18.,20.),(23.,25.),(32.,34.),(80.,82.)):
                nodes=np.cos(np.arange(degree+1)*math.pi/degree)
                jets=np.array([_static_collar_fields(sol.sol(lo+(z+1)*(hi-lo)/2),alpha,omega2) for z in nodes])
                coeff=[cheb.chebfit(nodes,jets[:,i]*(alpha if i==1 else 1),degree) for i in (1,2)]
                for z in np.linspace(-.8,.8,9):
                    state=sol.sol(lo+(z+1)*(hi-lo)/2)
                    r,beta,a,ap=_static_collar_fields(state,alpha,omega2);b=alpha*beta
                    bp,ad=[cheb.chebval(z,cheb.chebder(c))*2/(hi-lo) for c in coeff]
                    angular=ad+a*a+(bp+a*b)/r-1/r**4
                    radial=alpha*r/(r*r-alpha/3)*(2*state[0]/r**2+2*alpha*state[1]/r-1/(3*r**4))-a*b
                    errors.append((abs(angular),abs(ad-ap),abs(bp-radial)))
            maximum=np.max(errors,axis=0)
            derivative_checks.append(dict(degree=degree,angular_residual=float(maximum[0]),
                acceleration_derivative_error=float(maximum[1]),radial_derivative_error=float(maximum[2])))
        result['independent_local_polynomial_derivatives']=derivative_checks
    return result


def _collar_independent_proper_check(alpha,omega2):
    """Evolve b,a directly, without enforcing the mass/constraint reduction."""
    from scipy.optimize import brentq
    reduced,background,initial=_static_collar_solution(alpha,omega2,2)
    r,beta,a,ap=_static_collar_fields(initial,alpha,omega2)
    # Normalized y and beta retain the small initial radius derivatives.
    start=np.array([initial[0],beta,a,initial[2],initial[3],initial[1]])
    def rhs(distance,state):
        y,beta,a,n,x,h=state; r=1+alpha*y
        bprime=(r/(r*r-alpha/3)*(y*(2+alpha*y)/(r*r)-1/(3*r**4)-alpha*beta*beta)-a*beta)
        aprime=1/r**4-alpha*(bprime+a*beta)/r-a*a
        return beta,bprime,aprime,a,math.exp(-n),beta/6*(2*aprime+a*a-omega2*math.exp(-2*n))
    other=solve_ivp(rhs,(0,reduced.t[-1]),start,method='DOP853',
        rtol=2e-12,atol=2e-13,max_step=.04,dense_output=True)
    assert other.success
    errors=[];mass_identity=[]
    for ell in np.linspace(0,reduced.t[-1],201):
        z=reduced.sol(ell); expected=_static_collar_fields(z,alpha,omega2)
        q=other.sol(ell)
        reference=np.array([z[0],expected[1],expected[2],z[2],z[3],z[1]])
        errors.append(float(np.max(abs(q-reference)/(1+abs(reference)))))
        r=1+alpha*q[0]
        mass_identity.append(abs(q[5]-(q[0]**2/(2*r)-r*q[1]**2/2))/(1+abs(q[5])))
    return dict(alpha=alpha,maximum_state_scaled_difference=max(errors),
        independently_evolved_mass_geometric_identity=max(mass_identity),
        scope='Amplified alpha control resolves quantum corrections in an unconstrained independent ODE formulation; not a direct physical-alpha error bound.')


def polyakov_collar_checks():
    omega2=(2*math.pi/(math.pi+1.2))**2
    main=_static_collar_record(1e-9,omega2)
    tight=_static_collar_record(1e-9,omega2,2)
    controls=[_static_collar_record(alpha,omega2) for alpha in (.001,.0005,.00025)]
    independent=_collar_independent_proper_check(.001,omega2)
    keys=('proper_distance_from_start','optical_x','lapse_relative_to_center',
        'normal_radius_derivative','proper_log_lapse_derivative','mass_excess_over_alpha',
        'quantum_energy_over_alpha','quantum_radial_pressure_over_alpha','remaining_optical_path')
    refinement={key:max(abs(a[key]-b[key])/(1+abs(b[key])) for a,b in zip(main['checkpoints'],tight['checkpoints'])) for key in keys}
    result=dict(model='Finite static exterior collar with the SAME radius-independent Polyakov source and angular equation retained',
        parameters=dict(alpha=1e-9,omega_squared=omega2,start_x=.5,outer_radius=100.,
            total_assumed_optical_period=2*math.pi/math.sqrt(omega2)),
        physical=main,tighter_integration=tight,controls=controls,independent_proper_integration=independent,
        refinement_scaled_differences=refinement,
        leading_outer_expansion=dict(mass='M=1+alpha*(-5/(18R^3)+1/(3R^2)-1/18)+o(alpha) at fixed R>1',
            derivative='dM/dR=alpha*(5-4R)/(6R^4)+o(alpha)',
            scaled_lapse='alpha*N/(1-1/R) ->12/[pi*(Omega_squared-1)]',
            scaled_lapse_value=12/(math.pi*(omega2-1)),
            mass_coefficient_maximum=7/450,mass_coefficient_maximum_radius=1.25,
            nontrivial_mass_coefficient_zero_radius=(math.sqrt(21)-1)/2,
            scope='Matched small-alpha, fixed finite R expansion, not a uniform infinite-radius limit. The finite constant NS source still forbids standard finite-ADM spherical infinity.'),
        units='R in re, proper distance in re, x and t in center-normalized re/c units; M=Gmass/(c²re). Energy/pressure entries are 8piG re² T/(alpha*c⁴) in SI.',
        scope='Static continuation of the existing throat, not extension of the finite pulse or a stability test. A finite spherical collar is solved within the specified reduced action; full4D matter, nonspherical returning magnetic path, state self-consistency and formation remain unresolved. No physical source taper, shell or time offset introduced.',
        sources=['https://arxiv.org/html/1807.04726v3 (MMP matched regions and returning magnetic field lines)',
            'https://arxiv.org/abs/2102.03544 (different Boulware/regularized source; comparison only, not inserted)',
            'Our proper-distance continuation, constrained mass evolution and fixed-radius expansion'])
    records=[main,tight]+controls
    assert max(refinement.values())<2e-10
    assert max(r['maximum_overlap_difference'] for r in records)<2e-10
    assert max(r['maximum_constraint_scaled'] for r in records)<1e-12
    assert independent['maximum_state_scaled_difference']<2e-10
    assert independent['independently_evolved_mass_geometric_identity']<1e-7
    for key in ('maximum_fixed_radius_mass_coefficient_difference','maximum_scaled_lapse_difference'):
        assert all(b[key]<a[key] for a,b in zip(controls[:-1],controls[1:]))
        assert main[key]<controls[-1][key]
    derivatives=tight['independent_local_polynomial_derivatives'][-1]
    assert max(v for k,v in derivatives.items() if k!='degree')<2e-12
    for record in records:
        assert all(p['remaining_optical_path']>0 and p['normal_radius_derivative']>0
            and p['lapse_relative_to_center']>0 for p in record['plot'])
    result['verification']=dict(all_declared_checks_pass=True,
        maximum_refinement_scaled=max(refinement.values()),
        angular_residual_over_physical_alpha=derivatives['angular_residual']/1e-9,
        thresholds=dict(refinement=2e-10,overlap=2e-10,constraint=1e-12,
            independent_state=2e-10,independent_mass_identity=1e-7,polynomial_derivative=2e-12),
        scope='Refinement and independent derivative/ODE checks, not certified errors, stability, a full4D source or a closed exterior. The independent free evolution uses amplified alpha=.001.')
    result['matching_screens']=dict(required_geometry='re << R_boundary << physical mouth separation d; d is not solved here',
        optical='Remaining conformal path must be positive AND realized by an actual returning-field optical integral; positivity alone does not close the state.',
        localization='For a monopole B=q/(2R²), epsilon_B=2sqrt(2/q)*abs(dR/dproper_distance). Stored coefficient must be divided by sqrt(q); flux q is not fixed by this benchmark.',
        source='Spherical zero-mode separation is not a proof that this radius-independent action equals the complete four-dimensional renormalized source.')
    return result


def save_collar_benchmark():
    result=polyakov_collar_checks()
    project_path('polyakov_collar_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_STATIC_QUANTUM_COLLAR


# BEGIN_MAGNETIC_STATE_CLOSURE
def _magnetic_shapiro_log_integrand(z,nu):
    """Existing prolate path in eta=1-s; avoids loss of its tiny endpoint cap."""
    eta=math.exp(-z); s=1-eta
    if nu==2.: return 2/(1+s)
    h=eta*(2-eta); root=math.hypot(h,nu*s)
    u=(h+root)/nu
    um=(2-nu)*h/(nu*(u+1)-2*h); u=1+um
    up=s*(nu-2*u)/root
    difference=(um+eta)*(u+s)
    speed=.5*math.sqrt(difference*(up*up/(um*(u+1))+1/h))
    return speed*(2/(um+eta)+2/(u+s))*eta


def _magnetic_closure_paths(order,distance,cutoffs=(10.,30.,100.),quadrature=None):
    from scipy.integrate import quad
    from scipy.optimize import brentq
    from mmp_pulse import _magnetic_return_length,_magnetic_return_path
    if quadrature is None:
        nodes,weights=np.polynomial.legendre.leggauss(order)
        nu=nodes+1;weights=weights/2
    else:
        nu,weights=map(lambda v:np.asarray(v,float),quadrature)
        assert nu.ndim==weights.ndim==1 and nu.shape==weights.shape
        assert np.all((nu>0)&(nu<=2)) and np.all(weights>0) and abs(weights.sum()-1)<1e-13
    lengths=np.array([_magnetic_return_length(float(v)) for v in nu])
    pieces={};errors=[]
    for radius in cutoffs:
        isotropic=radius-1;xi=isotropic/distance
        correction=[]
        for v in nu:
            root_fn=lambda t:t*(2-xi*t)*(1-xi*t+2*xi)-2*v*(1-xi*t+xi)
            t=brentq(root_fn,0.,2.,xtol=1e-14)
            assert abs(root_fn(t))<2e-13
            value,error=quad(lambda z:_magnetic_shapiro_log_integrand(z,float(v)),
                0.,-math.log(xi*t),epsabs=2e-10,epsrel=2e-12,limit=150)
            correction.append(4*value);errors.append(4*error)
        e=(isotropic/(distance-isotropic))**2;slope=e/(1-e)
        cap_bound=2*isotropic*slope*slope/(math.sqrt(1+slope*slope)+1)
        pieces[str(radius)]=dict(shapiro=np.array(correction),isotropic=isotropic,
            flat_cap_excess_bound=cap_bound)
    # Independently check endpoint-safe expression on the old Cartesian path.
    direct_errors=[];axis_errors=[]
    for v in (.01,.5,1.,1.9):
        for z in (0.,.2,1.,3.):
            eta=math.exp(-z);rho,zz,speed=_magnetic_return_path(v,1-eta)
            direct=speed*(1/math.hypot(rho,zz+.5)+1/math.hypot(rho,zz-.5))*eta
            direct_errors.append(abs(direct-_magnetic_shapiro_log_integrand(z,v))/(1+abs(direct)))
    for radius in cutoffs:
        xi=(radius-1)/distance
        value=4*quad(lambda z:_magnetic_shapiro_log_integrand(z,2.),0.,-math.log(2*xi),epsabs=1e-11)[0]
        axis_errors.append(abs(value-4*math.log((1-xi)/xi)))
    assert max(direct_errors)<1e-12 and max(axis_errors)<1e-9
    return nu,weights,lengths,pieces,dict(maximum_quad_estimated_absolute_error=max(errors),
        stable_integrand_vs_existing_path=max(direct_errors),axis_analytic_error=max(axis_errors))


def _magnetic_closure_collar(omega2,alpha):
    from scipy.integrate import quad
    from scipy.optimize import brentq
    sol,_,_=_static_collar_solution(alpha,omega2,2)
    anchor=sol.y[:,-1];r,beta,a,ap=_static_collar_fields(anchor,alpha,omega2)
    clock=math.exp(anchor[2])/(alpha*beta)
    mass=1+alpha*alpha*anchor[1]
    rho_anchor=(r-mass+r*alpha*beta)/2
    data={}
    for radius in (10.,30.,100.):
        ell=float(sol.t[-1]) if radius==100. else brentq(lambda l:1+alpha*sol.sol(l)[0]-radius,0,sol.t[-1],xtol=2e-13)
        state=sol.sol(ell);r,beta,a,ap=_static_collar_fields(state,alpha,omega2)
        log_ratio=quad(lambda l:1/(1+alpha*sol.sol(l)[0]),ell,sol.t[-1],epsabs=2e-12,epsrel=2e-12)[0]
        data[str(radius)]=dict(optical_x=float(state[3]),lapse=math.exp(state[2]),normal_radius_derivative=alpha*beta,
            local_osculating_clock=math.exp(state[2])/(alpha*beta),
            isotropic_radius_from_collar=rho_anchor*math.exp(-log_ratio))
    return clock,data


def magnetic_state_closure_checks(orders=(64,128)):
    from functools import lru_cache
    from scipy.optimize import brentq
    alpha=1e-9;old_omega2=(2*math.pi/(math.pi+1.2))**2
    distance=1.2*12/(math.pi*alpha*(old_omega2-1))
    collar=lru_cache(maxsize=None)(lambda w:_magnetic_closure_collar(float(w),alpha))
    records=[];path_checks=[]
    for order in orders:
        nu,weights,f,pieces,checks=_magnetic_closure_paths(order,distance)
        path_checks.append(dict(order=order,**checks))
        def periods(w,cutoff,weak=True):
            clock,data=collar(w);piece=pieces[str(cutoff)];boundary=data[str(cutoff)]
            external=distance*f-2*piece['isotropic']+(piece['shapiro'] if weak else 0.)
            return 2*boundary['optical_x']+external/clock
        roots=[]
        for cutoff in (10.,30.,100.):
            residual=lambda w:w-float(weights@((2*math.pi/periods(w,cutoff))**2))
            lo,hi=1.1,3.5
            low,high=residual(lo),residual(hi)
            assert low<0<high
            root=brentq(residual,lo,hi,xtol=1e-13,rtol=1e-13)
            clock,data=collar(root);length=periods(root,cutoff);channel=(2*math.pi/length)**2
            mean=float(weights@channel);rms=math.sqrt(float(weights@((channel-mean)**2)))
            modes=[]
            for degree in range(1,7):
                polynomial=np.polynomial.legendre.legval(nu-1,[0]*degree+[1])
                coefficient=float((2*degree+1)*(weights@((channel-1)*polynomial))/3)
                modes.append(dict(degree=degree,source_coefficient=coefficient,
                    variance_fraction=coefficient**2/((2*degree+1)*(rms/3)**2)))
            boundary=data[str(cutoff)];piece=pieces[str(cutoff)]
            rho_error=abs(boundary['isotropic_radius_from_collar']-piece['isotropic'])
            # Radius-map error only, with a conservative weak optical factor;
            # does not bound nonspherical geometry, magnetic deformation or clock errors.
            radius_period_bound=4*rho_error*(1+4/(piece['isotropic']-rho_error))/clock
            axis_external=distance-2*piece['isotropic']+4*math.log((distance-piece['isotropic'])/piece['isotropic'])
            shortest=2*boundary['optical_x']+axis_external/clock
            roots.append(dict(cutoff=cutoff,omega_squared=root,root_residual=residual(root),
                bracket=[lo,hi],bracket_residuals=[low,high],clock_scale=clock,
                separation_over_clock=distance/clock,throat_optical_time=2*boundary['optical_x'],
                shortest_exterior_time=axis_external/clock,shortest_period=shortest,
                channel_omega_squared_mean=mean,channel_omega_squared_rms=rms,
                channel_source_mean_beta=(mean-1)/3,channel_source_rms_beta=rms/3,
                source_rms_over_mean=rms/(mean-1),source_angular_modes=modes,
                angular_dipole_coefficient=float(3*weights@(channel*(nu-1))),
                throat_to_shortest_exterior_time_ratio=2*boundary['optical_x']/(axis_external/clock),
                radial_localization_coefficient=2*math.sqrt(2)*boundary['normal_radius_derivative'],
                common_clock_local_calibration_relative_span=float(np.ptp([v['local_osculating_clock'] for v in data.values()])/clock),
                isotropic_radius_difference=rho_error,isotropic_radius_period_difference_bound=radius_period_bound,
                flat_cap_period_error_bound=piece['flat_cap_excess_bound']/clock,
                omitted_self_second_order_period_scale=2/(piece['isotropic']*clock),
                omitted_companion_clock_fraction_scale=1/distance,
                profile=dict(nu=nu.tolist(),weights=weights.tolist(),periods=length.tolist(),channel_omega_squared=channel.tolist())))
        leading_residual=lambda w:w-float(weights@((2*math.pi/(math.pi+distance*math.pi*alpha*(w-1)*f/12))**2))
        leading=brentq(leading_residual,1.1,3.5,xtol=1e-13,rtol=1e-13)
        raw=[];corrected=[]
        reference=roots[-1]['omega_squared']
        for cutoff in (10.,30.,100.):
            raw.append(float(weights@((2*math.pi/periods(reference,cutoff,False))**2)))
            corrected.append(float(weights@((2*math.pi/periods(reference,cutoff))**2)))
        old_mean=float(weights@((2*math.pi/periods(old_omega2,100.))**2))
        records.append(dict(order=order,roots=roots,leading_mean_field_omega_squared=leading,
            finest_cutoff_minus_leading=roots[-1]['omega_squared']-leading,
            cutoff_root_span=float(np.ptp([r['omega_squared'] for r in roots])),
            old_assumed_omega_squared=old_omega2,old_state_actual_path_mean=old_mean,
            old_state_mean_mismatch=old_mean-old_omega2,
            fixed_reference_omega_squared=reference,without_shapiro_mean_span=float(np.ptp(raw)),
            with_shapiro_mean_span=float(np.ptp(corrected))))
    roots=[root for record in records for root in record['roots']]
    maximum_residual=max(abs(row['root_residual']) for row in roots)
    root_refinement=max(abs(a['omega_squared']-b['omega_squared'])
        for a,b in zip(records[0]['roots'],records[-1]['roots']))
    rms_refinement=max(abs(a['channel_source_rms_beta']-b['channel_source_rms_beta'])
        for a,b in zip(records[0]['roots'],records[-1]['roots']))
    assert maximum_residual<2e-12 and root_refinement<2e-10 and rms_refinement<2e-10
    for row in roots:
        profile=row['profile'];weight=np.array(profile['weights']);period=np.array(profile['periods'])
        channel=np.array(profile['channel_omega_squared'])
        assert np.all(period>0) and np.all(np.isfinite(period)) and np.all(weight>0)
        assert abs(float(np.sum(weight))-1)<1e-14
        assert np.max(abs(channel-(2*math.pi/period)**2))<1e-13
        assert abs(float(weight@channel)-row['omega_squared'])<2e-12
        assert sum(mode['variance_fraction'] for mode in row['source_angular_modes'])<=1+1e-10
    result=dict(model='Finite-collar/weak-exterior magnetic-loop MONOPOLE state fixed point; nonspherical source not closed',
        parameters=dict(alpha=alpha,separation_over_re=distance,scaled_separation_alpha=alpha*distance,
            separation_convention='Fixed once so the prior single-delay1.2 corresponds to the shortest leading path in its old leading clock; not a measured or solved physical separation',
            cutoffs=[10.,30.,100.],clock_reference_cutoff=100.),
        quadrature_refinement=records,path_checks=path_checks,
        equations=dict(period='Lnu=2*x_b+[d*f(nu)-2*rho_b+Shapiro_nu]/A; common A=N100/b100; rho_b=Rb-1 with measured finite-collar coordinate discrepancy',
            state='Omega²=integral_0^2 dnu/2*(2pi/Lnu)²; collar recomputed at each trialOmega²',
            weak_exterior='Fixed point-magnetic paths, unit leading mass per mouth, dt=(1+2U)dl. Shapiro=4 integral_0^sc speed*(1/rhat_plus+1/rhat_minus) ds. Finite flat-cap excess is bounded separately.',
            endpoint='rho is isotropic, not areal; exact point-path cap root is reused. The finite collar isotropic relation d(logrho)/dl=1/R is integrated with one anchored normalization.',
            leading='Recovers z*(meanOmega²(delta/z)-1)=3 with ell_v=4/(pi*alpha),z=A/ell_v,delta=d/ell_v at leading order.'),
        scope='Only the spherical average and its finite-collar clock feedback are closed in this composite approximation. Angular state spread, omitted full4D polarization, companion metric, weak-metric magnetic-path deformation, exact junction, state preparation, pulse evolution, stability and formation remain unresolved. Cutoff agreement is not a physical accuracy bound. No manufactured independent channel clocks or renormalization of discarded channels.',
        sources=['https://arxiv.org/html/1807.04726v3#S5.SS5 (channel-dependent lengths and angular metric modes)',
            'https://arxiv.org/pdf/1409.3743 (weak optical metric; paths here are prescribed magnetic curves)',
            'Existing magnetic_return_path_checks/optical_matching_checks and our new same-source finite-collar feedback'])
    result['verification']=dict(all_declared_checks_pass=True,
        maximum_fixed_point_residual=maximum_residual,maximum_root_refinement_difference=root_refinement,
        maximum_source_rms_refinement_difference=rms_refinement,
        thresholds=dict(fixed_point_residual=2e-12,quadrature_root_and_rms=2e-10),
        scope='Positive periods, normalized full flux weights, bracketed roots, 64/128 quadrature refinement, Cartesian-path and analytic-axis cross-checks. Not certified errors or a stability test. Angular coefficients are source projections, not metric responses.')
    return result


def save_closure_benchmark():
    result=magnetic_state_closure_checks()
    project_path('polyakov_closure_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_MAGNETIC_STATE_CLOSURE


# BEGIN_ANGULAR_CASIMIR_RESPONSE
def angular_casimir_response_checks():
    """Sourced tensor particular on the leading product throat, not its collar."""
    import sympy as sp
    from mmp_pulse import _linear_diagonal_einstein_maxwell
    t,x,theta,phi=sp.symbols('t x theta phi',real=True)
    angular=sp.Function('h')(theta)
    metric=(-1/sp.cos(x)**2,1/sp.cos(x)**2,sp.Integer(1),sp.sin(theta)**2)
    perturbation=(angular,angular,sp.Integer(0),sp.Integer(0))
    f0=sp.zeros(4);f0[2,3]=2*sp.sin(theta);f0[3,2]=-f0[2,3]
    clean=lambda value:sp.simplify(sp.expand_trig(value).rewrite(sp.sin))
    einstein,maxwell,delta_f_squared=_linear_diagonal_einstein_maxwell(
        (t,x,theta,phi),metric,perturbation,f0,sp.zeros(4),sp.sin(theta)/sp.cos(x)**2,clean)
    expected=-(sp.diff(angular,theta,2)+sp.diff(angular,theta)/sp.tan(theta))/2
    residuals=[clean(value-(expected if pair in ((0,0),(1,1)) else 0)) for pair,value in einstein.items()]
    assert all(value==0 for value in residuals+maxwell) and delta_f_squared==0
    assert clean(sum(perturbation[i]/metric[i] for i in range(4)))==0
    record=json.loads(project_path('polyakov_closure_results.json').read_text(encoding='utf-8'))
    assert record['verification']['all_declared_checks_pass']
    alpha=record['parameters']['alpha'];cases=[]
    for quadrature in record['quadrature_refinement']:
        row=quadrature['roots'][-1];sigma=row['channel_source_rms_beta']
        modes=row['source_angular_modes'];coef=np.zeros(7)
        for mode in modes:
            j=mode['degree'];coef[j]=-2*mode['source_coefficient']/(j*(j+1))
        mu=np.linspace(-1,1,401);shape=np.polynomial.legendre.legval(mu,coef)
        # Verify the spherical operator via (1-mu²)h''-2mu*h', independently
        # of the coefficient division used to construct h.
        first=np.polynomial.legendre.legder(coef);second=np.polynomial.legendre.legder(first)
        laplace=(1-mu*mu)*np.polynomial.legendre.legval(mu,second)-2*mu*np.polynomial.legendre.legval(mu,first)
        source=np.r_[0.,[mode['source_coefficient'] for mode in modes]]
        operator_error=float(np.max(abs(-laplace/2+np.polynomial.legendre.legval(mu,source))))
        assert operator_error<1e-12
        tails=[]
        for cutoff in (2,6):
            remaining=sigma*sigma-sum(mode['source_coefficient']**2/(2*mode['degree']+1)
                for mode in modes if mode['degree']<=cutoff)
            assert remaining>0
            tails.append(dict(maximum_degree=cutoff,source_tail_rms=math.sqrt(remaining),
                uniform_metric_tail_bound=2*alpha*math.sqrt(remaining)/(cutoff+1)))
        cases.append(dict(order=quadrature['order'],alpha=alpha,
            metric_coefficients_over_alpha=coef.tolist(),angular_source_rms=sigma,
            uniform_all_modes_metric_bound=2*alpha*sigma,
            uniform_fixed_gauge_log_lapse_bound=alpha*sigma,
            maximum_six_mode_metric_magnitude=float(alpha*np.max(abs(shape))),
            sphere_operator_maximum_residual=operator_error,tail_bounds=tails,
            profile=dict(mu=mu.tolist(),metric_over_alpha=shape.tolist())))
    difference=float(np.max(abs(np.array(cases[0]['metric_coefficients_over_alpha'])-
        np.array(cases[-1]['metric_coefficients_over_alpha']))))
    assert difference<1e-11
    return dict(model='Leading product-throat angular Casimir tensor particular, with frozen source and no independently added homogeneous modes',
        conventions='Metric divided by re²: g0=sec²x(-dt²+dx²)+dOmega². Paper normalization F_theta_phi=2sin(theta), R-F²/4. Gravitationally normalized traceless source J_tt=J_xx=-alpha*(S-meanS); mixed/angular components and prescribed current vanish.',
        equations=dict(particular='h_tt=h_xx=H(mu)=-2alpha*sum_j>=1 s_j P_j(mu)/[j(j+1)]; all other metric and Maxwell perturbations zero',
            direct_check='delta[R_ab-.5*(F_ac F_b^c-.25*g_ab F²)]=-.5*Delta_S2 H for tt,xx; all other components zero; source is traceless so this also checks Einstein equations',
            fixed_gauge_lapse='delta ln N=-H(mu)*cos²x/2; fractional radial metric change=H(mu)*cos²x',
            bound='|H|<=2alpha*sigma_S; |H_tail_after_L|<=2alpha*sigma_tail/(L+1), by Cauchy-Schwarz, |P_j|<=1 and sum_j>L(2j+1)/[j²(j+1)²]=1/(L+1)²'),
        quadrature_refinement=cases,
        verification=dict(all_declared_checks_pass=True,ten_einstein_residuals=[str(v) for v in residuals],
            four_maxwell_residuals=[str(v) for v in maxwell],delta_f_squared=str(delta_f_squared),
            maximum_coefficient_refinement_difference=difference),
        scope='A sourced linear static particular on leading AdS2 times S2, distinct from the old homogeneous mass-six scalar mode. Not a solution on the finite collar or a choice of physical matching/boundary data. Bounds are analytic inequalities evaluated using numerical source moments, not certified continuum bounds. No global stability, formation, preparation, full quantum-current response or chronology change proved. Mean mode is excluded and handled by the existing spherical model.',
        sources=['https://arxiv.org/pdf/hep-th/9906056 Eqs3.28-3.30; our sourced particular and direct component check',
            'https://arxiv.org/html/1807.04726v3#S5.SS5; angular stress in the leading throat'])


def save_angular_casimir_benchmark():
    result=angular_casimir_response_checks()
    project_path('polyakov_angular_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_ANGULAR_CASIMIR_RESPONSE


# BEGIN_FLUX_CHANNEL_CURRENT
def _coupled_collar_seed(alpha,omega2,dipole,refinement=1):
    """Conditional outward seed test; stop at a declared coordinate size gate."""
    k=alpha/3;amplitude=alpha*dipole
    background,_,_=_static_collar_solution(alpha,omega2,2)
    def algebra(values):
        r,A,B,d,f,m,z,v=values
        Bp=r/(r*r-k)*(1-1/r**2-k/r**4-B*B)-A*B
        Ap=1/r**4-(Bp+A*B)/r-A*A
        p=-k*(A*A+d)/r**2;rho=k*(2*Ap+A*A-d)/r**2
        C=1-r**4*p/2;Cp=k*(r*B*(A*A+d)+r*r*A*(Ap-d))
        K=r*B+k*A;Dh=1+k*d+K*(A-B/r);Dv=-A*B*(r*r-k)/r
        h=(Dv*v+2*C/r**2*(m+K*z)+k*f)/Dh
        q=-A*(h-v)+B/r*(h+v)+2*C*z/r**2
        zp=-(A+Cp/C)*z+(2*m+h+v)/(C*r*r)+r*r*(rho*h-p*v)/(2*C)
        vp=r**3/(B*(r*r-k))*(-(2/r**2-1/r**4-rho)*v-rho*h+(2/r**4+rho)*m+
            k/r**2*(-(A+B/r)*q+(h+v)/r**2+(4*m+2*v)/r**4+(2*Ap+A*A)*(h-v)-f))
        direction=np.array([B,Ap,Bp,-2*A*d,-2*A*f,z,zp,vp])
        master=2*m+2*B*B*v/(B*B-r*Bp-1)
        return h,q,direction,(rho,p,C,Dh,master)
    def values(ell,state):
        bg=background.sol(ell);r,beta,A,Ap=_static_collar_fields(bg,alpha,omega2)
        invn2=math.exp(-2*bg[2])
        return np.r_[r,A,alpha*beta,omega2*invn2,invn2/k,state]
    def rhs(ell,state): return algebra(values(ell,state))[2][-3:]
    def gate(ell,state):
        h,_,_,_=algebra(values(ell,state))
        return .01-amplitude*max(abs(h),abs(state[0]),abs(state[2]))
    gate.terminal=True;gate.direction=-1
    initial=np.array([0.,0.,-math.exp(-2*background.y[2,0])])
    sol=solve_ivp(rhs,(0,background.t[-1]),initial,method='Radau',dense_output=True,events=gate,
        rtol=2e-9/refinement,atol=np.array([2e-11*alpha,2e-11*alpha,2e-11])/refinement,
        max_step=.15/refinement)
    assert sol.success and len(sol.t_events[0])==1
    rows=[];errors=[]
    for ell in np.linspace(0,sol.t[-1],121):
        z=values(ell,sol.sol(ell));r,A,B,d,f,m,mp,v=z
        h,q,direction,(rho,p,C,Dh,master)=algebra(z)
        # Differentiate the algebraic pressure constraint and mixed equation
        # independently along the background/state flow, with complex steps.
        derivatives=[]
        for index in range(8):
            probe=z.astype(complex);probe[index]+=1e-25j
            hp,qp,_,_=algebra(probe)
            derivatives.append([hp.imag/1e-25,qp.imag/1e-25])
        hprime,hsecond=np.array(derivatives).T@direction
        vp=direction[-1];Ap=direction[1]
        temporal=B/r*vp+(2/r**2-1/r**4-rho)*v+rho*h-2*m/r**4
        source_tt=k/r**2*(hsecond+A*(q-vp)+(2*Ap+A*A)*(h-v)-f)+m*rho
        radial=B/r*q-(h+v)/r**2+(2*m+v)/r**4
        source_ll=k/r**2*(-A*q+d*(h-v)-f)+m*p
        angular=hsecond+(2*A+B/r)*q-(A+B/r)*vp-(h+v)/r**2-(4*m+2*v)/r**4
        scale=1+sum(abs(item) for item in (h,q,hsecond,vp,m,v,f*k))
        errors.append(max(abs(hprime-q),abs(temporal-source_tt),abs(radial-source_ll),abs(angular))/scale)
        rows.append(dict(proper_distance=float(ell),radius=float(r),
            metric_time_fraction=float(amplitude*h),metric_radial_fraction=float(amplitude*v),
            magnetic_flux_fraction=float(amplitude*m),invariant_master=float(amplitude*master),
            pressure_constraint_denominator=float(Dh)))
    assert max(errors)<2e-8 and min(row['pressure_constraint_denominator'] for row in rows)>0
    return dict(alpha=alpha,refinement=refinement,solver_evaluations=sol.nfev,
        gate_radius=rows[-1]['radius'],gate_proper_distance=float(sol.t[-1]),
        maximum_independent_equation_scaled_residual=float(max(errors)),plot=rows,
        scope='Amplified-alpha numerical control with frozen dipole seed, m=m_l=0 at x=.5 and v=-alpha*s1/N²; h fixed by pressure constraint. Ends at1% coordinate perturbation, not a physical horizon or instability threshold. Physical end data are not chosen by this test.')


def _flux_collar_gravity_checks():
    """Check the reduced collar equations against all relevant 4D components."""
    import sympy as sp
    from mmp_pulse import _linear_diagonal_einstein_maxwell
    t,l,theta,phi=sp.symbols('t l theta phi',real=True)
    N,R,h,v,m=[sp.Function(name)(l) for name in ('N','R','h','v','m')]
    rho,p=sp.symbols('rho p',real=True)
    g=(-N*N,sp.Integer(1),R*R,R*R*sp.sin(theta)**2);Y=sp.cos(theta)
    perturb=(g[0]*h*Y,v*Y,sp.Integer(0),sp.Integer(0))
    f0=sp.zeros(4);f1=sp.zeros(4);f0[2,3]=2*sp.sin(theta);f0[3,2]=-f0[2,3]
    coordinates=(t,l,theta,phi);potential=(0,0,0,m*sp.sin(theta)**2)
    for a in range(4):
        for b in range(4):
            f1[a,b]=sp.diff(potential[b],coordinates[a])-sp.diff(potential[a],coordinates[b])
    clean=lambda value:sp.factor(sp.trigsimp(sp.expand_trig(value)))
    reversed_einstein,maxwell,_=_linear_diagonal_einstein_maxwell(
        coordinates,g,perturb,f0,f1,N*R*R*sp.sin(theta),clean)
    trace=sum(reversed_einstein[i,i]/g[i] for i in range(4))
    contraction=(v-h)*(rho+p)*Y/2
    einstein={pair:clean(value-(g[pair[0]]*trace/2+perturb[pair[0]]*(rho-p)/2-
        g[pair[0]]*contraction/2 if pair[0]==pair[1] else 0)) for pair,value in reversed_einstein.items()}
    A=sp.diff(N,l)/N;B=sp.diff(R,l)
    background={rho:(1-B*B)/R**2-2*sp.diff(R,l,2)/R-1/R**4,
        p:(B*B-1)/R**2+2*A*B/R+1/R**4}
    angular_background={sp.diff(N,l,2):N/R**4-(N*sp.diff(R,l,2)+sp.diff(N,l)*B)/R}
    expected={
        (0,0):N*N*((B/R)*sp.diff(v,l)+(2/R**2-1/R**4-rho)*v+rho*h-2*m/R**4)*Y,
        (1,1):((B/R)*sp.diff(h,l)-(h+v)/R**2+(2*m+v)/R**4)*Y,
        (1,2):(-A*(h-v)+(B/R)*(h+v)-sp.diff(h,l)+2*sp.diff(m,l)/R**2)*(-sp.sin(theta))/2,
        (2,2):R*R*(sp.diff(h,l,2)+(2*A+B/R)*sp.diff(h,l)-(A+B/R)*sp.diff(v,l)-
            (h+v)/R**2-(4*m+2*v)/R**4)/2*Y}
    expected[(3,3)]=expected[(2,2)]*sp.sin(theta)**2
    residuals=[clean((value-expected.get(pair,0)).subs(background).subs(angular_background))
        for pair,value in einstein.items()]
    assert all(value==0 for value in residuals),residuals
    expected_maxwell=N*(sp.diff(m,l,2)+A*sp.diff(m,l)-(2*m+h+v)/R**2)*sp.sin(theta)
    magnetic=[clean(value-(expected_maxwell if index==3 else 0)) for index,value in enumerate(maxwell)]
    assert all(value==0 for value in magnetic)
    return dict(ten_einstein_equation_reduction_residuals=[str(value) for value in residuals],
        four_maxwell_equation_reduction_residuals=[str(value) for value in magnetic],
        background_quantum_trace_terms_retained=True)


def flux_channel_current_checks():
    """Action-derived static dipole source/current on a spherical collar."""
    import sympy as sp
    x,theta,Q,k,w,dw=sp.symbols('x theta Q k w dw',real=True,nonzero=True)
    sigma,R,ht,hr,m=[sp.Function(name)(x) for name in ('sigma','R','Ht','Hr','m')]
    N=sp.exp(sigma);sx=sp.diff(sigma,x);sxx=sp.diff(sigma,x,2)
    rho=k*(2*sxx-sx*sx-w)/(R*R*N*N);p=-k*(sx*sx+w)/(R*R*N*N)
    st=k/R**2*(sp.diff(ht,x,2)-sx*sp.diff(hr,x)+(2*sxx-sx*sx)*(ht-hr)-dw)+m*N*N*rho/Q
    sr=k/R**2*(-sx*sp.diff(ht,x)+w*(ht-hr)-dw)+m*N*N*p/Q
    tau=R*R*p*sp.diff(m,x)/(2*Q)
    # Same channel action: label density and displacement work must vary together.
    lagrangian=N*N*R*R*m*(-rho*ht+p*hr)/(2*Q)+R**4*p*sp.diff(m,x)**2/(4*Q**2)
    current=sp.simplify((sp.diff(lagrangian,m)-sp.diff(sp.diff(lagrangian,sp.diff(m,x)),x))/(N*N*R*R))
    expected_current=-(rho*ht-p*hr)/(2*Q)-sp.diff(R**4*p*sp.diff(m,x),x)/(2*Q**2*N*N*R*R)
    action_error=sp.simplify(current-expected_current)
    assert action_error==0
    # Direct variation of divergence, using mixed tensors and volume density.
    coords=(sp.Symbol('t'),x,theta,sp.Symbol('phi'));Y=sp.cos(theta)
    metric=(-N*N,N*N,R*R,R*R*sp.sin(theta)**2)
    perturb=(metric[0]*ht*Y,metric[1]*hr*Y,sp.Integer(0),sp.Integer(0))
    base=sp.diag(N*N*rho,N*N*p,0,0);delta=sp.zeros(4)
    delta[0,0]=st*Y;delta[1,1]=sr*Y;delta[1,2]=delta[2,1]=-tau*sp.sin(theta)
    volume=N*N*R*R*sp.sin(theta);vf=(ht+hr)*Y/2
    wards=[]
    for b in (1,2):
        expression=0
        for a in range(4):
            mixed0=base[a,b]/metric[a]
            mixed1=delta[a,b]/metric[a]-perturb[a]*base[a,b]/metric[a]**2
            expression+=(sp.diff(volume*(mixed1+vf*mixed0),coords[a])-vf*sp.diff(volume*mixed0,coords[a]))/volume
            expression-=sp.diff(perturb[a],coords[b])*base[a,a]/(2*metric[a]**2)
            expression-=sp.diff(metric[a],coords[b])*(delta[a,a]/metric[a]**2-2*perturb[a]*base[a,a]/metric[a]**3)/2
        if b==2: expression-=Q*sp.sin(theta)*current
        wards.append(sp.simplify(sp.expand(expression)))
    assert wards==[0,0]
    C=1-R**4*p/(2*Q**2)
    direct=sp.diff(m,x,2)-N*N*(2*m+Q*(ht+hr))/R**2+N*N*R*R*current
    reduced=sp.diff(C*sp.diff(m,x),x)-N*N*(2*m+Q*(ht+hr))/R**2-N*N*R*R*(rho*ht-p*hr)/(2*Q)
    maxwell_error=sp.simplify(direct-reduced)
    assert maxwell_error==0
    record=json.loads(project_path('polyakov_closure_results.json').read_text(encoding='utf-8'))
    alpha=record['parameters']['alpha'];case=record['quadrature_refinement'][-1]['roots'][-1]
    omega2=case['omega_squared'];dipole=case['source_angular_modes'][0]['source_coefficient']
    sol,background,_=_static_collar_solution(alpha,omega2,2)
    def row(r,lapse,a,ap,b,ell):
        rr=alpha*(2*ap+a*a-omega2/lapse**2)/(3*r*r)
        pp=-alpha*(a*a+omega2/lapse**2)/(3*r*r)
        correction=-r**4*pp/2
        return dict(radius=r,proper_distance_from_collar_start=ell,lapse=lapse,
            lapse_log_derivative=a,radial_derivative=b,quantum_rho=rr,quantum_pressure=pp,
            magnetic_principal_correction=correction,principal_correction_over_alpha=correction/alpha,
            unbent_current_over_alpha_H=(rr+pp)/(2*alpha*lapse**2),
            unbent_current_for_leading_dipole=(rr+pp)*(-alpha*dipole)/(2*lapse**2))
    samples=[]
    for ell in np.linspace(0,sol.t[-1],301):
        state=sol.sol(ell);r,beta,a,ap=_static_collar_fields(state,alpha,omega2)
        samples.append(row(r,math.exp(state[2]),a,ap,alpha*beta,float(ell)))
    center=background.fields(0.);n2=float(center['e2s'])
    throat=row(float(center['R']),math.sqrt(n2),float(center['sx'])/math.sqrt(n2),
        (float(center['sxx'])-float(center['sx'])**2)/n2,float(center['Rx'])/math.sqrt(n2),None)
    assert all(item['quantum_pressure']<0 and item['magnetic_principal_correction']>0 for item in samples+[throat])
    gravity=_flux_collar_gravity_checks()
    controls=[_coupled_collar_seed(.001,omega2,dipole,level) for level in (1,4)]
    gate_difference=abs(controls[0]['gate_radius']-controls[1]['gate_radius'])
    assert gate_difference<2e-7
    return dict(model='Flux-following static Polyakov channel action: local dipole stress, tilt and induced current, evaluated on the same finite collar',
        parameters=dict(alpha=alpha,omega_squared=omega2,dipole_source_coefficient=dipole,Q=1.,collar_starts_at_x=.5,outer_radius=100.),
        equations=dict(metric='ds²=N²(-dt²+dx²)+R²dOmega²; δgtt=gtt*Ht*cos(theta), δgxx=gxx*Hr*cos(theta), δAphi=m*sin²(theta)',
            conventions='S=8piG*T, div Fcal=-Jcal, div S=Fcal.Jcal/2, Fcal_theta_phi=2Q sin(theta); k=alpha/3. Prime means conformal x derivative.',
            background='rho=k*(2sigma_second-sigma_prime²-w)/(R²N²), p=-k*(sigma_prime²+w)/(R²N²)',
            source_tt='δS_tt/cos(theta)=k/R²*[Ht_second-sigma_prime*Hr_prime+(2sigma_second-sigma_prime²)(Ht-Hr)-δw]+(m/Q)N²rho',
            source_xx='δS_xx/cos(theta)=k/R²*[-sigma_prime*Ht_prime+w(Ht-Hr)-δw]+(m/Q)N²p',
            tilt='δS_xtheta=-R²p*m_prime*sin(theta)/(2Q); δS_angular=0 at linear order',
            action='Γ2 density before sin³(theta) integration: N²R²*m*(-rho*Ht+p*Hr)/(2Q)+R⁴p*m_prime²/(4Q²), plus intrinsic2D metric response; convention δΓ=integral sqrt(g)*(Sδg+JδA)',
            current='Jphi=-(rho*Ht-p*Hr)/(2Q)-(R⁴p*m_prime)_prime/(2Q²N²R²)',
            maxwell='(C*m_prime)_prime-N²[2m+Q(Ht+Hr)]/R²=N²R²(rho*Ht-p*Hr)/(2Q), C=1-R⁴p/(2Q²)>1',
            principal='C-1=alpha*R²[(d_l lnN)²+w/N²]/(6Q²); l=proper distance',
            state='δw includes prescribed angular seed and complete-loop variation. For fixed spin/holonomy ground state, δw_geometry=-2w*δL/L. Local collar contributes δL=integral(Hr-Ht)dx/2; exterior/end terms unresolved.',
            nonvacuum_gravity='If E1=δ(R_ab-T_EM_ab/2), use δEinstein=E1-g tr(E1)/2-h R0/2+g*(h^ab E0_ab)/2, E0=R_ab-T_EM_ab/2. Electrovac-only trace conversion is insufficient here.'),
        coupled_proper_distance_equations=dict(convention='Q=1; h=Ht,v=Hr,m=magnetic potential; A=N_l/N,B=R_l; dots denote proper l derivatives. δS coefficients are divided by cos(theta).',
            mixed='h_dot=-A(h-v)+(B/R)(h+v)+(2C/R²)m_dot',
            radial='(B/R)h_dot-(h+v)/R²+(2m+v)/R⁴=δS_ll',
            temporal='(B/R)v_dot+(2/R²-1/R⁴-rho)v+rho*h-2m/R⁴=δS_tt/N²',
            angular='h_ddot+(2A+B/R)h_dot-(A+B/R)v_dot-(h+v)/R²-(4m+2v)/R⁴=0',
            maxwell='N^-1(C*N*m_dot)_dot-(2m+h+v)/R²=R²(rho*h-p*v)/2',
            source_tt='δS_tt/N²=k/R²*[h_ddot+A(h_dot-v_dot)+(2A_dot+A²)(h-v)-δw/N²]+m*rho',
            source_ll='δS_ll=k/R²*[-A*h_dot+(w/N²)(h-v)-δw/N²]+m*p',
            remaining='A consistent boundary-value problem still needs gauge regularization at B->0, complete-loop δw and coupled physical end data. These equations are derived and component-checked; no new full-collar angular solution is asserted.'),
        throat=throat,collar_samples=samples,conditional_coupled_seed_controls=controls,
        maximum_sampled_principal_correction=max(item['magnetic_principal_correction'] for item in samples+[throat]),
        verification=dict(all_declared_checks_pass=True,action_variation_residual=str(action_error),
            direct_radial_and_angular_ward_residuals=[str(v) for v in wards],maxwell_reduction_residual=str(maxwell_error),
            negative_pressure_and_positive_principal=True,amplified_control_gate_radius_refinement_difference=gate_difference,**gravity),
        scope='Derived local closure of the specified flux-channel approximation, not full4D quantum determinant. Main rows evaluate coefficients/current; separate amplified-alpha rows solve a conditional seed IVP, not physical mouth matching or a manufactured material. Positive principal coefficient excludes this static magnetic degeneracy only, not global/dynamic stability. Complete-loop state variations and physical boundary data remain unresolved. Physical-alpha raw gauge evolution was numerically ill-conditioned and is not certified by the amplified control.',
        sources=['https://arxiv.org/html/1807.04726v3 sections2.3,3.1 and5: lowest-Landau flux-channel framework',
            'QUANTUM_COMPLETION.md: existing pulled-back channel action and displacement work; present static collar source/current derivation is ours'])


def save_flux_channel_benchmark():
    result=flux_channel_current_checks()
    project_path('polyakov_channel_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_FLUX_CHANNEL_CURRENT


# BEGIN_INVARIANT_COLLAR
def _invariant_collar_algebra_checks():
    """Independent raw equations and gauge flow, evaluated at 70 digits."""
    import mpmath as mp
    with mp.workdps(70):
        def bg(r,A,B,d,k):
            Bp=r/(r*r-k)*(1-1/r**2-k/r**4-B*B)-A*B
            Ap=1/r**4-(Bp+A*B)/r-A*A
            rho,p=k/r**2*(2*Ap+A*A-d),-k/r**2*(A*A+d)
            C=1-r**4*p/2;Cp=k*(r*B*(A*A+d)+r*r*A*(Ap-d))
            K=r*B+k*A;H,D=1+k*d+K*(A-B/r),B*B-r*Bp-1
            return Bp,Ap,rho,p,C,Cp,K,H,D
        def raw(r,A,B,d,k,m,z,v,eta):
            Bp,Ap,rho,p,C,Cp,K,H,D=bg(r,A,B,d,k)
            h=(-A*B*(r*r-k)*v/r+2*C/r**2*(m+K*z)+eta)/H
            q=-A*(h-v)+B/r*(h+v)+2*C*z/r**2
            zp=-(A+Cp/C)*z+(2*m+h+v)/(C*r*r)+r*r*(rho*h-p*v)/(2*C)
            F=r**3/(r*r-k)*(-(2/r**2-1/r**4-rho)*v-rho*h
                +(2/r**4+rho)*m+k/r**2*(-(A+B/r)*q+(h+v)/r**2
                +(4*m+2*v)/r**4+(2*Ap+A*A)*(h-v))-eta/r**2)
            return mp.matrix([z,zp,F/B]),F
        def gauge(r,A,B,d,k):
            D=bg(r,A,B,d,k)[-1]
            return mp.matrix([-2,2/(r*B),2*D/B**2])
        reports=[]
        jets=[('1.2','.3','.001','1'),
            ('1.000000001','.55','3.3333333333333333333333333333333e-10','.5')]
        for jet in jets:
            r,A,k,d=map(mp.mpf,jet)
            B=-r*A+mp.sqrt(r*r*A*A+1-1/r**2-k*(A*A+d))
            point=[r,A,B,d,k];Bp,Ap,rho,p,C,Cp,K,H,D=bg(*point)
            flow=[B,Ap,Bp,-2*A*d,mp.mpf(0)]
            g=gauge(*point);expected=-g/(r*B)
            for i in range(5):
                for n in range(3):
                    def component(value):
                        args=point.copy();args[i]=value
                        return gauge(*args)[n]
                    expected[n]+=mp.diff(component,point[i])*flow[i]
            actual,_=raw(*point,*g,mp.mpf(0))
            gauge_error=max(abs(actual[n]-expected[n])/(1+abs(expected[n])) for n in range(3))
            matrix,forcing,_=_invariant_collar_matrix(r,A,B,d,Ap,Bp,k)
            matrix_error=mp.mpf(0)
            for i in range(2):
                P,X=mp.mpf(i==0),mp.mpf(i==1)
                velocity,F=raw(*point,P/2,X/2,mp.mpf(0),mp.mpf(0))
                projected=[X+2*B*F/D,2*velocity[1]-2*F/(r*D)]
                matrix_error=max(matrix_error,max(abs(projected[n]-matrix[n,i]) for n in range(2)))
            velocity,F=raw(*point,mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(1))
            projected=[2*B*F/D,2*velocity[1]-2*F/(r*D)]
            source_error=max(abs(projected[n]-forcing[n]) for n in range(2))
            errors=[gauge_error,matrix_error,source_error]
            assert max(errors)<mp.mpf('1e-50'),list(map(str,errors))
            reports.append(dict(radius=float(r),B=float(B),gauge_scaled_error=float(gauge_error),
                matrix_absolute_error=float(matrix_error),source_absolute_error=float(source_error),
                source_per_eta=[float(v) for v in forcing]))
        return dict(decimal_precision=70,all_declared_checks_pass=True,cases=reports,
            scope='Two independent constraint-compatible jets, including B~1.6e-9; local gauge and algebra checks, not a continuum stability theorem.')


def _invariant_collar_matrix(r,A,B,d,Ap,Bp,k):
    """Local gauge quotient; factored sources retain the small throat terms."""
    rho=k*(2*Ap+A*A-d)/r**2
    C=1+k*r*r*(A*A+d)/2
    Cp=k*(r*B*(A*A+d)+r*r*A*(Ap-d))
    K=r*B+k*A;H=1+k*d+K*(A-B/r);D=B*B-r*Bp-1;S=r*r-k
    E=A*A+d+(1-B*B)/r**2
    FP=r/S*(1/r**2+r*r*rho/2+2*k/r**4+k*E*C/(r*r*H))
    FX=r/S*k*C/r**2*(E*K/H-(A+B/r))
    hP=C/(r*r*H);hX=K*hP;J=(2/r**2+r*r*rho)/C
    matrix=np.array([[2*B*FP/D,1+2*B*FX/D],
        [2/(C*r*r)+J*hP-2*FP/(r*D),-A-Cp/C+J*hX-2*FX/(r*D)]])
    forcing=np.array([2*B*(B*B-1-r*A*B)/(r*D*H),
        -k*(Ap-d)*(2*r*(Bp+2*A*B)+k*(A*A+d))/(C*D*H)])
    return matrix,forcing,(D,H,C)


def _invariant_collar_setup(alpha,omega2,refinement):
    background,_,_=_static_collar_solution(alpha,omega2,refinement)
    k=alpha/3
    def field(ell):
        bg=background.sol(ell);r,beta,A,Ap=_static_collar_fields(bg,alpha,omega2)
        B=alpha*beta;eta=math.exp(-2*bg[2])
        # Reuse scaled mass/radius variables: 1-R^-2-B² loses digits here.
        Bp=alpha*r/(r*r-k)*(2*bg[0]/r**2+2*alpha*bg[1]/r-1/(3*r**4))-A*B
        matrix,forcing,denominators=_invariant_collar_matrix(r,A,B,omega2*eta,Ap,Bp,k)
        return r,A,B,eta,matrix,forcing,denominators
    return background,field


def _invariant_collar_seed(alpha,dipole,background,field,refinement):
    r,A,B,eta,_,_,(D,_,_)=field(0.)
    initial=np.array([-2*B*B*eta/D,2*B*eta/(r*D)])
    amplitude=alpha*dipole
    def rhs(ell,state):
        _,_,_,eta,matrix,forcing,_=field(ell)
        return matrix@state+forcing*eta
    def gate(ell,state): return .01-amplitude*abs(state[0])
    gate.terminal=True;gate.direction=-1
    sol=solve_ivp(rhs,(0,background.t[-1]),initial,method='DOP853',
        rtol=1e-10/refinement,atol=1e-15*alpha/refinement,
        max_step=.1/refinement,dense_output=True,events=gate)
    assert sol.success and len(sol.t_events[0])==1
    samples=[dict(proper_distance=float(ell),radius=float(field(ell)[0]),
        invariant_master=float(amplitude*sol.sol(ell)[0]),
        invariant_companion=float(amplitude*sol.sol(ell)[1]))
        for ell in np.linspace(0,sol.t[-1],161)]
    result=dict(alpha=alpha,refinement=refinement,solver_evaluations=sol.nfev,
        gate_radius=samples[-1]['radius'],gate_proper_distance=float(sol.t[-1]),
        initial_normalized_invariants=initial.tolist(),plot=samples,
        scope='Same arbitrary local seed as68, frozen delta_w=3*s1; stop at |P|=.01, an invariant-amplitude convention DIFFERENT from the old coordinate gate. Not a physical instability or transport threshold.')
    return result,sol


def _invariant_collar_boundary_map(alpha,dipole,background,field,refinement,_return_solutions=False):
    """Full linear response operator; neither endpoint value is selected."""
    from scipy.integrate import solve_bvp,quad
    end=float(background.t[-1])
    def homogeneous(ell,state): return field(ell)[4]@state
    solutions=[solve_ivp(homogeneous,interval,[0.,1.],method='DOP853',
        rtol=1e-11/refinement,atol=1e-12/refinement,max_step=.1/refinement,dense_output=True)
        for interval in ((0,end),(end,0))]
    forward,backward=solutions
    assert all(sol.success for sol in solutions)
    pf,xf=forward.y[:,-1];pb,xb=backward.y[:,-1]
    endpoint_matrix=np.array([[xb/pb,1/pf],[1/pb,xf/pf]])
    trace_integral=quad(lambda ell:float(np.trace(field(ell)[4])),0,end,
        epsabs=2e-10,epsrel=2e-11,limit=120)[0]
    wronskian_error=abs(math.log(-endpoint_matrix[1,0]/endpoint_matrix[0,1])-trace_integral)
    # Solve N*(P,X), with source divided by alpha, to resolve the weak outer
    # tail. Restore chi=alpha²*delta_w/3 and undo N afterwards.
    def rhs(positions,states):
        answer=np.empty_like(states)
        for index,ell in enumerate(positions):
            _,A,_,eta,matrix,forcing,_=field(ell)
            answer[:,index]=(matrix+A*np.eye(2))@states[:,index]+forcing*math.sqrt(eta)/alpha
        return answer
    def boundary(left,right): return np.array([left[0],right[0]])
    mesh=np.linspace(0,end,401)
    forced=solve_bvp(rhs,boundary,mesh,np.zeros((2,len(mesh))),
        tol=1e-8/refinement**2,max_nodes=6000)
    assert forced.success,forced.message
    ell=np.linspace(0,end,1201);weighted=forced.sol(ell)
    inverse_lapse=np.exp(-background.sol(ell)[2]);response=weighted*inverse_lapse
    chi=alpha*alpha*dipole
    rows=[dict(proper_distance=float(t),radius=float(field(t)[0]),
        invariant_master=float(chi*response[0,i]),invariant_companion=float(chi*response[1,i]),
        unit_source_master=float(response[0,i])) for i,t in enumerate(ell)]
    affine=forced.y[1,[0,-1]]*inverse_lapse[[0,-1]]
    # Independently compare collocation polynomials with the differential equations.
    weighted_rhs=rhs(ell,weighted)
    physical_rhs=np.column_stack([field(t)[4]@response[:,i]+field(t)[5]*field(t)[3]/alpha
        for i,t in enumerate(ell)])
    residual=(forced.sol(ell,1)-weighted_rhs)*inverse_lapse
    residual_max=float(np.max(abs(residual)/(1+abs(physical_rhs))))
    denominators=np.array([field(t)[-1] for t in ell])
    assert wronskian_error<2e-8 and residual_max<5e-7
    assert np.min(abs(denominators[:,0]))>.02 and np.min(denominators[:,1:])>0
    result=dict(refinement=refinement,outer_radius=float(field(end)[0]),
        proper_distance=end,endpoint_matrix=endpoint_matrix.tolist(),
        unit_source_endpoint_companion=affine.tolist(),
        source_multiplier_for_saved_dipole=chi,maximum_sampled_forced_master=float(max(abs(chi*response[0]))),
        homogeneous_shooting_end_values=dict(forward_P=float(pf),backward_P=float(pb)),
        wronskian_log_error=wronskian_error,collocation_nodes=len(forced.x),
        collocation_maximum_scaled_equation_residual=residual_max,
        minimum_abs_D=float(np.min(abs(denominators[:,0]))),
        minimum_H=float(np.min(denominators[:,1])),minimum_C=float(np.min(denominators[:,2])),plot=rows,
        convention='[X(0),X(L)] = endpoint_matrix@[P(0),P(L)] + chi*unit_source_endpoint_companion; chi=alpha²*delta_w/3. X uses increasing proper-distance orientation at BOTH ends. delta_w=3*s1 in plotted forced diagnostic.',
        collocation_scaling='Collocation solves N*(P,X) with forcing divided by alpha; saved rows undo N and apply chi. This resolves the weak intermediate tail that unweighted absolute-error control missed.',
        scope='Linear collar response operator only. The zero-P endpoints define an affine basis, not physical mouth conditions. Tiny basis response does not establish the actual solution, stability, global state or formation.')
    return (result,(forward,backward,forced)) if _return_solutions else result


def invariant_collar_checks():
    """Target-alpha gauge quotient, independent old-variable comparison and full-collar map."""
    old=json.loads(project_path('polyakov_channel_results.json').read_text(encoding='utf-8'))
    alpha=old['parameters']['alpha'];omega2=old['parameters']['omega_squared']
    dipole=old['parameters']['dipole_source_coefficient']
    algebra=_invariant_collar_algebra_checks()
    controls=[];target=[];maps=[];raw_errors=[]
    for a in (.001,alpha):
        for level in (1,2):
            bg,field=_invariant_collar_setup(a,omega2,2*level)
            result,sol=_invariant_collar_seed(a,dipole,bg,field,level)
            if a==alpha:
                target.append(result)
                maps.append(_invariant_collar_boundary_map(a,dipole,bg,field,level))
            else:
                controls.append(result)
                if level==2:
                    for row in old['conditional_coupled_seed_controls'][-1]['plot']:
                        ell=row['proper_distance']
                        if ell<=sol.t[-1]:
                            raw_errors.append(abs(a*dipole*sol.sol(ell)[0]-row['invariant_master']))
    gate_delta=abs(target[0]['gate_radius']-target[1]['gate_radius'])
    matrix_delta=float(np.max(abs(np.array(maps[0]['endpoint_matrix'])-maps[1]['endpoint_matrix'])))
    offdiag_relative=max(abs(maps[0]['endpoint_matrix'][i][j]/maps[1]['endpoint_matrix'][i][j]-1)
        for i,j in ((0,1),(1,0)))
    source_delta=float(np.max(abs(np.array(maps[0]['unit_source_endpoint_companion'])-
        maps[1]['unit_source_endpoint_companion'])))
    profile_delta=max(abs(a['unit_source_master']-b['unit_source_master'])
        for a,b in zip(maps[0]['plot'],maps[1]['plot']))
    regional=[]
    for lo,hi in ((0,1),(1,5),(5,10),(10,25),(25,100),(100,127)):
        pairs=[(a['unit_source_master'],b['unit_source_master'])
            for a,b in zip(maps[0]['plot'],maps[1]['plot']) if lo<=b['proper_distance']<hi]
        peak=max(abs(b) for a,b in pairs);difference=max(abs(a-b) for a,b in pairs)
        regional.append(dict(proper_distance_interval=[lo,hi],maximum_unit_source_master=peak,
            maximum_difference=difference,difference_relative_to_region_peak=difference/peak))
    assert gate_delta<2e-8 and max(raw_errors)<2e-8
    assert matrix_delta<2e-8 and offdiag_relative<2e-7 and source_delta<2e-8 and profile_delta<2e-8
    assert max(row['difference_relative_to_region_peak'] for row in regional)<.001
    return dict(model='Gauge-invariant static dipole response on the finite quantum collar',
        parameters=dict(alpha=alpha,omega_squared=omega2,dipole_source_coefficient=dipole,
            collar_start_x=.5,outer_radius=100.,Q=1.),
        equations=dict(invariants='P=2m+2B²*v/D; X=2m_l-2B*v/(R*D), D=B²-R*B_l-1',
            gauge='(m,m_l,v)=V*(-2,2/(RB),2D/B²), V_l=-V/(RB); both invariants vanish',
            quotient='Choose local representative v=0,m=P/2,m_l=X/2. F=B*v_l: P_l=X+2BF/D; X_l=2m_ll-2F/(RD). Coefficient derivatives multiply v=0, so no missing connection term.',
            forcing='eta=k*delta_w/N² in physical variables. bP=2B(B²-1-RAB)/(R*D*H); bX=-k(A_l-d)[2R(B_l+2AB)+k(A²+d)]/(C*D*H). Matrix implemented in _invariant_collar_matrix.',
            normalization='Seed IVPs divide fields by alpha*s1, giving eta=1/N² for delta_w=3*s1. Boundary affine basis instead divides the forcing by alpha; physical chi=alpha²*delta_w/3 restores it.',
            remaining='Only the collar portion of the complete-loop state/matching problem is determined. Fixed state susceptibility about uniform mean w; nonlinear angular products and changing loop state remain unresolved.'),
        algebra_checks=algebra,amplified_controls=controls,target_seed_refinement=target,
        boundary_map_refinement=maps,affine_profile_regional_refinement=regional,
        verification=dict(all_declared_checks_pass=True,target_gate_radius_difference=gate_delta,
            old_raw_variable_master_maximum_difference=max(raw_errors),endpoint_matrix_maximum_difference=matrix_delta,
            offdiagonal_maximum_relative_difference=offdiag_relative,affine_endpoint_maximum_difference=source_delta,
            affine_profile_maximum_difference=profile_delta),
        scope='Converged conditional target-alpha calculation and full finite-collar linear response operator. The assumed throat and state are retained. No end values are physically selected, no complete-loop feedback, dynamical/nonlinear stability, material generation, transport or 13Sep2027->12Sep2027 passage is proved.',
        sources=['https://arxiv.org/html/1807.04726v3: underlying lowest-Landau channel framework only',
            'Our quotient and cancellation-safe algebra derived from the component-checked increment68 action; not a published discovery'])


def save_invariant_collar_benchmark():
    result=invariant_collar_checks()
    project_path('polyakov_invariant_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_INVARIANT_COLLAR


# BEGIN_EMBEDDED_COLLAR
def _embedded_collar_geometry_checks():
    """Direct embedded metric/normal/curvature calculation, with arbitrary jets."""
    import sympy as s
    e,l,th,T,cap=s.symbols('e l th T cap',real=True)
    R,N,B=s.symbols('R N B',positive=True)
    A,Ap,Bp,h,q,v,m,z,ze,b,V=s.symbols('A Ap Bp h q v m z ze b V',real=True)
    coords=s.symbols('t',real=True),l,th,s.symbols('phi',real=True)
    trunc=lambda f:s.expand(f.subs(e,0)+e*s.diff(f,e).subs(e,0))
    simp=lambda f:s.factor(s.trigsimp(s.cancel(f)))
    lin=lambda f:simp(s.diff(f,e).subs(e,0))
    Rl=R+B*l+Bp*l*l/2;N2=N*N*(1+2*A*l+(Ap+2*A*A)*l*l)
    gd=[-N2*(1+e*(h+q*l)*s.cos(th)),1+e*v*s.cos(th),Rl*Rl,Rl*Rl*s.sin(th)**2]
    gi=[trunc(1/g) for g in gd]
    at=lambda f:f.subs({l:0,th:T},simultaneous=True)
    def pull(f):
        f0=f.subs(e,0)
        return simp(at(f0)+e*at(s.diff(f,e).subs(e,0)
            +ze*s.cos(th)*s.diff(f0,l)-b*s.sin(th)*s.diff(f0,th)))
    def gam(a,i,j):
        terms=((s.diff(gd[a],coords[i]) if a==j else 0)
            +(s.diff(gd[a],coords[j]) if a==i else 0)
            -(s.diff(gd[i],coords[a]) if i==j else 0))
        return pull(trunc(gi[a]*terms/2))
    tang=[s.Matrix([1,0,0,0]),s.Matrix([0,-e*ze*s.sin(T),1-e*b*s.cos(T),0]),s.Matrix([0,0,0,1])]
    accel=[s.zeros(4,1),s.Matrix([0,-e*ze*s.cos(T),e*b*s.sin(T),0]),s.zeros(4,1)]
    pg=[pull(g) for g in gd]
    ind=[trunc(sum(pg[i]*E[i]**2 for i in range(4))) for E in tang]
    raw=[0,1,e*ze*s.sin(th),0]
    norm2=trunc(sum(gi[i]*raw[i]**2 for i in range(4)))
    nc=[pull(trunc(s.sympify(raw[i])/s.sqrt(norm2))) for i in range(4)]
    def extr(j):
        E=tang[j]
        return trunc(-sum(nc[a]*(accel[j][a]+sum(gam(a,i,k)*E[i]*E[k]
            for i in (1,2) for k in (1,2))) for a in (1,2))) if j==1 else trunc(-nc[1]*gam(1,0,0))
    lapse=lin(s.log(-ind[0]/N**2)/2)/s.cos(T)
    area=R*R*s.sin(T)*s.sqrt(ind[1]/R**2*ind[2]/(R*R*s.sin(T)**2))
    angular=lin(ind[1]/R**2)/2/s.cos(T)
    kt=lin(extr(0)/ind[0])/s.cos(T);ko=lin(extr(1)/ind[1])/s.cos(T)
    pot=-2*s.cos(th)+e*(m+z*l)*s.sin(th)**2
    ffield=[s.diff(pot,c) for c in coords]
    flux=trunc(sum(pull(ffield[i])*tang[1][i] for i in (1,2)))
    contrast=s.integrate(lin(flux/2-area/R**2),(T,0,cap))/s.sin(cap)**2
    normalfield=lin(sum(pull(gi[i])*nc[i]*pull(ffield[i]) for i in (1,2)))/s.sin(T)**2
    got=list(map(simp,[lapse,angular,kt,ko,contrast,normalfield]));D=B*B-R*Bp-1
    expected=[h/2+A*ze,B*ze/R-b,q/2-A*v/2+Ap*ze,
        -B*v/(2*R)-D*ze/R**2,m/2-B*ze/R,z+2*ze/R**2]
    residuals=[simp(a-c) for a,c in zip(got,expected)]
    gauge={h:2*A*R*V/B,v:2*D*V/B**2,m:-2*V,z:2*V/(R*B),
        q:2*Ap*R*V/B+2*A*D*V/B**2,ze:-R*V/B,b:-V}
    gauge_residuals=[simp(a.subs(gauge,simultaneous=True)) for a in got]
    assert residuals==[0]*6 and gauge_residuals==[0]*6
    # A moving endpoint cancels the optical integral's pure coordinate term.
    radial=s.Function('radial')(l);lapse_fn=s.Function('N')(l)
    optical_gauge=s.simplify((2*s.diff(radial,l)-2*s.diff(lapse_fn,l)/lapse_fn*radial)/(2*lapse_fn)
        -s.diff(radial/lapse_fn,l))
    kk,dd,DD=s.symbols('k d D',nonzero=True)
    CC=1+kk*R**2*(A*A+dd)/2;HH=kk*(Ap+A*A-A*B/R)-DD
    Eb=A*A+dd+(1-B*B)/R**2;rho=kk*(2*Ap+A*A-dd)/R**2
    fp=R/(R**2-kk)*(1/R**2+R**2*rho/2+2*kk/R**4+kk*Eb*CC/(R**2*HH))
    safe=kk*(R**2*CC*(Ap+A*A-A*B/R+Eb)+HH*R**4*(Ap-dd)-CC*DD+2*HH)/(R**2*HH*(R**2-kk)*DD)
    factor_error=s.factor(CC/(R**2*HH)+R*fp/DD-safe)
    assert optical_gauge==0 and factor_error==0
    return dict(geometry_residuals=list(map(str,residuals)),gauge_residuals=list(map(str,gauge_residuals)),
        optical_gauge_residual=str(optical_gauge),optical_P_factorization_residual=str(factor_error),
        method='Arbitrary Taylor jets, pulled-back metric, normalized level-set normal and Christoffel extrinsic curvature; no equations of motion for geometry. Optical factorization uses background D+H=k(A_l+A²-AB/R).')


def _embedded_collar_coefficients(bg,alpha,omega2):
    """Six embedded observables and a regular optical integrand, Q=re=1."""
    k=alpha/3;r,beta,A,Ap=_static_collar_fields(bg,alpha,omega2)
    B=alpha*beta;inverse_lapse=math.exp(-bg[2]);d=omega2*inverse_lapse**2
    Bp=alpha*r/(r*r-k)*(2*bg[0]/r**2+2*alpha*bg[1]/r-1/(3*r**4))-A*B
    D=B*B-r*Bp-1;C=1+k*r*r*(A*A+d)/2;K=r*B+k*A
    H=1+k*d+K*(A-B/r);E=A*A+d+(1-B*B)/r**2
    hp=C/(r*r*H);hx=K*hp;qp=(-A+B/r)*hp;qx=(-A+B/r)*hx+C/r**2
    cut=np.array([[hp/2,hx/2,A,0],[0,0,B/r,-1],
        [qp/2,qx/2,Ap,0],[0,0,-D/r**2,0],[.25,0,-B/r,0],[0,.5,2/r**2,0]])
    eta_column=np.array([1/(2*H),0,(-A+B/r)/(2*H),0,0,0])
    # Explicit k factor retains the optical P response in the near-vacuum limit.
    optical_p=k*(r*r*C*(Ap+A*A-A*B/r+E)+H*r**4*(Ap-d)-C*D+2*H)/(r*r*H*(r*r-k)*D)
    FX=r/(r*r-k)*k*C/r**2*(E*K/H-(A+B/r))
    optical_x=hx+r*FX/D;optical_eta=(D+B*B-1-r*A*B)/(H*D)
    optical=-inverse_lapse*np.array([optical_p,optical_x,optical_eta])/2
    compatibility=np.zeros((2,6));compatibility[0,0]=1;compatibility[1,2]=1
    compatibility[0,4:6]=[-2*hp,-hx];compatibility[1,4:6]=[-2*qp,-qx]
    compatibility[0,3]=-A*(2*k*hp-r*r)/D
    compatibility[1,3]=-(2*C/r**2+2*k*A*(-A+B/r)*hp-Ap*r*r)/D
    error=float(np.max(abs(compatibility@cut)))
    assert error<2e-12
    return dict(radius=r,inverse_lapse=inverse_lapse,cut_matrix=cut,eta_column=eta_column,
        optical=optical,compatibility=compatibility,compatibility_error=error,
        background_quantum_pressure=-k*(A*A+d)/r**2)


def _embedded_collar_record(alpha,omega2,dipole,refinement):
    from scipy.integrate import quad_vec
    background,field=_invariant_collar_setup(alpha,omega2,2*refinement)
    response,(forward,backward,forced)=_invariant_collar_boundary_map(
        alpha,dipole,background,field,refinement,True)
    end=float(background.t[-1]);k=alpha/3
    endpoint_matrix=np.array(response['endpoint_matrix'])
    affine=np.array(response['unit_source_endpoint_companion'])
    blocks=[];cuts=[]
    for index,ell in enumerate((0.,end)):
        c=_embedded_collar_coefficients(background.sol(ell),alpha,omega2)
        inputs=np.zeros((4,7));inputs[0,index]=1
        inputs[1,:2]=endpoint_matrix[index];inputs[1,6]=alpha*k*affine[index]
        inputs[2,2+index]=1;inputs[3,4+index]=1
        assembled=c['cut_matrix']@inputs
        assembled[:,6]+=c['eta_column']*k*c['inverse_lapse']**2
        blocks.append(assembled)
        cuts.append(dict(endpoint=('inner','outer')[index],proper_distance=ell,radius=c['radius'],
            inverse_lapse=c['inverse_lapse'],cut_matrix=c['cut_matrix'].tolist(),
            eta_column=c['eta_column'].tolist(),compatibility_matrix=c['compatibility'].tolist(),
            compatibility_residual=c['compatibility_error'],
            background_quantum_pressure=c['background_quantum_pressure']))
    def integrand(ell):
        bg=background.sol(ell);c=_embedded_collar_coefficients(bg,alpha,omega2)
        q=c['optical'];invn=c['inverse_lapse']
        left=backward.sol(ell)/backward.y[0,-1]
        right=forward.sol(ell)/forward.y[0,-1]
        source=forced.sol(ell)*invn
        # Give the four quadratures comparable magnitudes before error control.
        return np.array([q[:2]@left/alpha,q[:2]@right/alpha,
            q[:2]@source/alpha,q[2]*invn**2])
    integrals,quadrature_error=quad_vec(integrand,0,end,
        epsabs=1e-12/refinement,epsrel=1e-10/refinement,limit=1000)
    integrals[:3]*=alpha
    optical_row=np.array([integrals[0],integrals[1],-cuts[0]['inverse_lapse'],
        cuts[1]['inverse_lapse'],0.,0.,alpha*k*integrals[2]+k*integrals[3]])
    cumulative=solve_ivp(lambda ell,_:integrand(ell),(0,end),np.zeros(4),method='DOP853',
        rtol=1e-10/refinement,atol=1e-12/refinement,max_step=.1/refinement,dense_output=True)
    assert cumulative.success
    cumulative_total=cumulative.y[:,-1].copy();cumulative_total[:3]*=alpha
    cumulative_error=float(np.max(abs(cumulative_total-integrals)))
    assert cumulative_error<2e-10
    samples=[]
    for ell in np.linspace(0,end,501):
        value=cumulative.sol(ell)
        per_dw=alpha*alpha*k*value[2]+k*value[3]
        samples.append(dict(proper_distance=float(ell),radius=float(field(ell)[0]),
            cumulative_optical_change_per_delta_w=float(per_dw),
            cumulative_optical_change_per_alpha_delta_w=float(per_dw/alpha)))
    matrix=np.vstack(blocks+[optical_row[None,:]])
    return dict(refinement=refinement,cuts=cuts,assembled_response_matrix=matrix.tolist(),
        optical_row=optical_row.tolist(),optical_quadrature_components=integrals.tolist(),
        scaled_quadrature_error_estimate=float(quadrature_error),
        independent_optical_integration_maximum_difference=cumulative_error,
        optical_change_for_saved_dipole_fixed_invariants_and_embeddings=float(optical_row[-1]*3*dipole),
        plot=samples,
        scope='One static collar, x=.5 toR100. All P and embedding endpoint values remain inputs. The last optical coefficient holds those inputs fixed; it is not the complete-loop susceptibility. Source scaling is linear about the uniform mean state.')


def embedded_collar_checks():
    record=json.loads(project_path('polyakov_invariant_results.json').read_text(encoding='utf-8'))
    alpha=record['parameters']['alpha'];omega2=record['parameters']['omega_squared']
    dipole=record['parameters']['dipole_source_coefficient']
    geometry=_embedded_collar_geometry_checks()
    cases=[_embedded_collar_record(alpha,omega2,dipole,level) for level in (1,2)]
    a,b=[np.array(case['assembled_response_matrix']) for case in cases]
    relative_columns=[float(np.max(abs(a[:,i]-b[:,i]))/np.max(abs(b[:,i]))) for i in range(7)]
    optical_relative=[abs(x/y-1) if y else abs(x-y) for x,y in zip(cases[0]['optical_row'],cases[1]['optical_row'])]
    assert max(relative_columns)<2e-7 and max(optical_relative)<2e-7
    return dict(model='Embedded finite-collar geometry and gauge-invariant optical state-response functional',
        parameters=dict(alpha=alpha,omega_squared=omega2,dipole_source_coefficient=dipole,
            Q=1.,collar_start_x=.5,outer_radius=100.),
        input_order=['P_inner','P_outer','Zeta_inner','Zeta_outer','beta_inner','beta_outer','delta_w'],
        observable_order=['log_lapse','angular_half_metric_fraction','mixed_time_extrinsic_curvature',
            'mixed_angular_extrinsic_curvature','cap_flux_area_contrast','normal_F_phi'],
        assembled_output_order='Six observables atinner, six atouter, then collar optical-length change. Units re=1 and common center-time normalization; increasing-l normal atbothends.',
        equations=dict(embedding='l=c+zeta*cosTheta,theta=Theta+beta*dTheta(cosTheta); Zeta=zeta+RBv/(2D), beta_invariant=beta+B²v/(2D).',
            raw_observables='h/2+A*zeta; B*zeta/R-beta; h_l/2-A*v/2+A_l*zeta; -B*v/(2R)-D*zeta/R²; m/2-B*zeta/R; m_l+2*zeta/R².',
            regular_observables='h0/2+A*Zeta; B*Zeta/R-beta_invariant; q0/2+A_l*Zeta; -D*Zeta/R²; P/4-B*Zeta/R; X/2+2*Zeta/R².',
            representative='h0=C(P+KX)/(R²H)+eta/H; q0=(-A+B/R)h0+C*X/R² is the actual h_l at the local v=0 representative, NOT the derivative of h0 along a sequence of representatives.',
            optical='deltaL=-integral[(h0+R*F0/D)/(2N)]dl+[Zeta/N]_inner^outer; F0=B*v_l in the local representative.',
            gauge='tau=RBv/(2D) transforms by radial displacement T. Original optical integrand changes by (T/N)_l; moving endpoints change by -[T/N], cancelling exactly.',
            state='delta_w_geometry=-2w*deltaL_full/L_full. This calculation supplies one collar contribution only. Outside paths, other side, spin/gauge transport and endpoint matching remain required.'),
        geometry_checks=geometry,refinement=cases,
        verification=dict(all_declared_checks_pass=True,maximum_assembled_column_relative_difference=max(relative_columns),
            maximum_optical_coefficient_relative_difference=max(optical_relative),
            maximum_cut_compatibility_residual=max(c['compatibility_residual'] for r in cases for c in r['cuts'])),
        scope='A boundary response and local optical functional within the assumed quantum-channel model, not a completed spacetime junction or source material. The retained nonzero quantum normal pressure prevents an abrupt exact source-free RN splice without further stress or a specified layer. No stability, formation or yesterday passage follows.',
        sources=['https://arxiv.org/html/1807.04726v3: channel/Casimir framework and global path dependence',
            'https://arxiv.org/abs/gr-qc/0201054: primary junction geometry framework',
            'Our proper-coordinate embedded-surface and optical-response derivations; reuses previous RN boundary meanings, not its vacuum equations'])


def save_embedded_collar_benchmark():
    result=embedded_collar_checks()
    project_path('polyakov_embedded_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_EMBEDDED_COLLAR


# BEGIN_CENTER_EXTERIOR
def _quantum_overlap_jet(state,alpha,omega2):
    """Quantum collar jet; caller supplies mp.workdps(70), b multiplies eta."""
    import mpmath as mp
    asmp=lambda value:value if isinstance(value,mp.mpf) else mp.mpf(str(value))
    aa,ww=map(asmp,(alpha,omega2));k=aa/3
    y,hm,n,_=map(asmp,state)
    r=1+aa*y;beta2=y*y/r**2-2*hm/r;beta=mp.sqrt(beta2);B=aa*beta
    d=ww*mp.exp(-2*n);bterm=y/r**2+aa*hm/r-d/6
    A=2*bterm/(r*beta+mp.sqrt(r*r*beta2+2*bterm/3))
    Ap=1/r**4-aa/(r*r-k)*(2*y/r**2+2*aa*hm/r-1/(3*r**4))-A*A
    Bp=aa*r/(r*r-k)*(2*y/r**2+2*aa*hm/r-1/(3*r**4))-A*B
    matrix,forcing,_=_invariant_collar_matrix(r,A,B,d,Ap,Bp,k)
    return dict(r=r,N=mp.exp(n),B=B,A=A,M=mp.matrix(matrix.tolist()),
                b=mp.matrix(forcing.tolist()))


def _rn_omission_diagnostic(alpha,omega2,clock,delta_w,cg):
    """Weight the discarded RN terms by a declared growing-only incident field."""
    import mpmath as mp
    from scipy.optimize import brentq
    assert alpha>0 and omega2>0 and clock>0 and delta_w!=0
    cases=[]
    with mp.workdps(70):
        aa,ww,cc,dw,gg=map(lambda value:mp.mpf(str(value)),(alpha,omega2,clock,delta_w,cg))
        k=aa/3;amplitude=-k*dw/cc**2
        for refinement in (2,4):
            background,_,_=_static_collar_solution(alpha,omega2,refinement)
            end=float(background.t[-1]);rows=[]
            for radius in (10.,30.,100.):
                ell=end if radius==100 else brentq(
                    lambda value:1+alpha*background.sol(value)[0]-radius,0,end,xtol=1e-12)
                state=background.sol(ell)
                # Re-evaluate the scaled background jets before subtracting matrices.
                jet=_quantum_overlap_jet(state,aa,ww)
                r,B,full=(jet[key] for key in ('r','B','M'))
                n=mp.mpf(str(state[2]))
                br=1-1/r;F=br**2;L=3*r-2
                rn=mp.matrix([[-2*br/(r*L),1],[6/(r*L),-1/r**2+2*br/(r*L)]])
                g=3*r*(r-1)**2/L
                gp=3*((r-1)*(3*r-1)*L-3*r*(r-1)**2)/L**2
                classical=gg*mp.matrix([g,br*(gp+2*g/(r*L))])
                omitted=(full-rn)*classical
                radial_correction=(1-B/br)*rn*classical
                pulled=omitted+radial_correction
                direct=(full-(B/br)*rn)*classical
                direct_scale=max(map(abs,direct))
                equality_error=max(map(abs,pulled-direct))/direct_scale if direct_scale else mp.mpf(0)
                assert equality_error<mp.mpf('1e-50')
                eta_rn=k*dw/(cc**2*F);retained=2*r*br/L*eta_rn
                ratios=[omitted[0]/retained,r*omitted[1]/retained]
                pulled_ratios=[pulled[0]/retained,r*pulled[1]/retained]
                rows.append(dict(radius=float(r),proper_distance=float(ell),
                    classical_invariants=list(map(float,classical)),
                    matrix_difference=[[float((full-rn)[i,j]) for j in range(2)] for i in range(2)],
                    omitted_rhs=list(map(float,omitted)),retained_state_rhs=[float(retained),0.],
                    signed_ratios_Pprime_and_R_Xprime=list(map(float,ratios)),
                    weighted_infinity_ratio=float(max(map(abs,ratios))),
                    radial_correction=list(map(float,radial_correction)),
                    pulled_classical_operator_mismatch=list(map(float,pulled)),
                    pulled_signed_ratios_Pprime_and_R_Xprime=list(map(float,pulled_ratios)),
                    pulled_direct_equality_relative_error=float(equality_error),
                    collar_eta_over_RN_eta_minus_one=float(cc**2*F*mp.exp(-2*n)-1)))
            cases.append(dict(background_refinement=refinement,rows=rows))
        amplitude=float(amplitude)
    errors=[];pulled_errors=[]
    for coarse,fine in zip(cases[0]['rows'],cases[1]['rows']):
        den=max(map(abs,fine['omitted_rhs']))
        errors.append(max(abs(a-b) for a,b in zip(coarse['omitted_rhs'],fine['omitted_rhs']))/den if den else 0.)
        den=max(map(abs,fine['pulled_classical_operator_mismatch']))
        pulled_errors.append(max(abs(a-b) for a,b in zip(coarse['pulled_classical_operator_mismatch'],fine['pulled_classical_operator_mismatch']))/den if den else 0.)
    assert max(errors)<5e-4 and max(pulled_errors)<5e-4
    return dict(parameters=dict(alpha=alpha,omega_squared=omega2,clock_scale=clock,
            delta_w_positive_mouth=delta_w,growing_coefficient=cg,source_amplitude=amplitude),
        refinement=cases,maximum_omitted_rhs_relative_refinement=max(errors),
        maximum_pulled_operator_relative_refinement=max(pulled_errors),
        all_declared_checks_pass=True,
        equations='O=(M_collar-M_RN)Y_class is matrix-only. Pullback correction=(1-B/sqrt(F))*M_RN*Y_class; their sum=(M_collar-(B/sqrt(F))*M_RN)Y_class. Y_class=c_g*(g,sqrt(F)*(g_prime+2g/[r(3r-2)])); g=3r(r-1)^2/(3r-2). S=b_RN*eta_RN, b_RN=(2r*sqrt(F)/(3r-2),0), eta_RN=k*delta_w/(clock²*F).',
        normalization='M=Q=re=1. Ratios are O_P/S_P and R*O_X/S_P; the factor R gives equal dimensions. S_X=0, so no componentwise division by it is attempted. delta_w uses the supplied positive-mouth Y=1-nu orientation.',
        precision='70-digit point algebra on two double-precision adaptive background histories; this is not a 70-digit background solution. Both matrix-only and pulled classical operator mismatches are refined, supporting three significant digits; the direct pullback identity is also checked.',
        scope='Growing-only classical test input, not a selected matched solution; the decaying amplitude remains free. Matrix-only differences omit the radial derivative pullback. The pulled mismatch includes that correction because d/dl=B*d/dr on the actual collar, but neither quantity is the full forced residual. Includes local mean-background and quantum-response differences relative to extremal RN; no global error, stability or failure theorem.')


def _center_even_response(alpha,omega2,collar,refinement):
    """Reflection-even static sector: P even and X odd, with X(0)=0."""
    bg=StaticBackground(alpha,omega2,.5);k=alpha/3
    def field(x):
        f=bg.fields(x);N=math.sqrt(f['e2s']);A=f['sx']/N;B=f['Rx']/N
        Ap=(f['sxx']-f['sx']**2)/N**2;Bp=(f['Rxx']-f['sx']*f['Rx'])/N**2
        matrix,forcing,_=_invariant_collar_matrix(f['R'],A,B,omega2/N**2,Ap,Bp,k)
        return N,matrix,forcing
    parity_error=0.;J=np.diag([1.,-1.])
    for x in np.linspace(0,.5,11):
        N,M,b=field(x);Nm,Mm,bm=field(-x)
        parity_error=max(parity_error,float(np.max(abs(Mm+J@M@J))),float(np.max(abs(bm+J@b))),abs(N-Nm))
    assert parity_error<1e-13
    def rhs(x,flat):
        N,M,b=field(x);out=N*M@flat.reshape(2,2)
        out[:,1]+=b/(3*alpha*N)  # Forced column divided by alpha²*delta_w.
        return out.ravel()
    sol=solve_ivp(rhs,(0,.5),np.array([[1.,0.],[0.,0.]]).ravel(),method='DOP853',
        rtol=1e-12/refinement,atol=1e-13/refinement,max_step=.01/refinement,dense_output=True)
    assert sol.success
    u=sol.y[:,-1].reshape(2,2);kt=u[1,0]/u[0,0];jt=u[1,1]-kt*u[0,1]
    K=np.asarray(collar['endpoint_matrix']);j=np.asarray(collar['unit_source_endpoint_companion'])
    den=kt-K[0,0];t=K[0,1]/den;s=(j[0]/3-jt)/den
    ke=K[1,1]+K[1,0]*t;je=j[1]/3+K[1,0]*s
    # Verify both independent inputs after eliminating the artificial inner cut.
    residuals=[kt*t-K[0,0]*t-K[0,1],kt*s+jt-K[0,0]*s-j[0]/3,
        ke-K[1,0]*t-K[1,1],je-K[1,0]*s-j[1]/3]
    assert max(abs(v) for v in residuals)<1e-15
    return dict(refinement=refinement,center_to_inner=u.tolist(),inner_slope=kt,inner_source_scaled=jt,
        inner_from_outer=t,inner_from_source_scaled=s,outer_slope=ke,outer_source_scaled=je,
        center_from_outer=t/u[0,0],center_from_source_scaled=(s-u[0,1])/u[0,0],
        parity_maximum_residual=parity_error,composition_residuals=residuals,
        convention='Xinner=kt*Pinner+alpha²*delta_w*jt; Pinner=t*Pouter+alpha²*delta_w*s; Xouter=ke*Pouter+alpha²*delta_w*je. P even, X odd under proper-distance reflection with fixed internal flux labels.')


def _center_rn_algebra():
    """Reuse the quotient at k=0 and verify the source-affine RN map exactly."""
    import sympy as s
    r=s.symbols('r',positive=True);B=(r-1)/r;F=B**2;L=3*r-2
    M,b,_=_invariant_collar_matrix(r,1/r**2,B,s.Integer(0),-2*(r-1)/r**4,(r-1)/r**3,s.Integer(0))
    M=s.Matrix(M);b=s.Matrix(b);g=3*r*(r-1)**2/L
    decay=r*(30*r**2-45*r+17)/(10*(r-1)**3*L)
    hp=-4*r/(3*L);xp=2*r/(L*B)
    residuals=[]
    for h in (g,decay):
        state=s.Matrix([h,B*s.diff(h,r)+2*B*h/(r*L)])
        residuals.extend(s.factor(v) for v in B*state.diff(r)-M*state)
    yp=s.Matrix([hp,xp]);residuals.extend(s.factor(v) for v in B*yp.diff(r)-M*yp+b/F)
    momentum=s.factor(B*xp-2*F*hp/(r*L)-2*r/L)
    residuals.append(s.factor(momentum-F*s.diff(hp,r)))
    kappa=-F*s.diff(decay,r)/decay
    residuals.append(s.factor(F*s.diff(g,r)+kappa*g-3/decay))
    assert residuals==[0]*8
    return dict(residuals=list(map(str,residuals)),growing=str(g),decaying=str(decay),
        particular_per_amplitude=str(hp),particular_X_per_amplitude=str(xp),
        momentum='Pi=sqrt(F)*X-2F*P/[r(3r-2)]-2r*A_dip/(3r-2)',
        outer_relation='Pi+kappa*P=3*c_g/d_mode+A_dip*(Pi_p+kappa*h_p)+Xi_ext; kappa=-F*d_mode_prime/d_mode. Xi_ext is not selected by this RN calculation.')


def center_exterior_checks():
    """Close the center condition, orient the full-flux source, retain exterior uncertainty."""
    import sympy as s
    closure=json.loads(project_path('polyakov_closure_results.json').read_text(encoding='utf-8'))
    invariant=json.loads(project_path('polyakov_invariant_results.json').read_text(encoding='utf-8'))
    alpha=closure['parameters']['alpha'];separation=closure['parameters']['separation_over_re']
    ref=closure['quadrature_refinement'][-1]['roots'][-1];omega2=ref['omega_squared'];clock=ref['clock_scale']
    assert invariant['parameters']['alpha']==alpha and invariant['parameters']['omega_squared']==omega2
    moments=[]
    for case in closure['quadrature_refinement']:
        root=case['roots'][-1];p=root['profile'];weights=np.asarray(p['weights']);nu=np.asarray(p['nu'])
        w=np.asarray(p['channel_omega_squared']);period=np.asarray(p['periods'])
        assert abs(sum(weights)-1)<1e-14
        dw=3*np.sum(weights*(1-nu)*w)
        moment=3*np.sum(weights*(1-nu)/(root['clock_scale']*period/separation)**2)
        amplitude=-alpha*dw/(3*root['clock_scale']**2)
        cross=-4*math.pi**2*alpha*moment/(3*separation**2)
        assert abs(amplitude/cross-1)<1e-12 and abs(dw+root['angular_dipole_coefficient'])<1e-12
        moments.append(dict(order=case['order'],delta_w_positive_mouth=float(dw),full_flux_moment=float(moment),
            source_amplitude=float(amplitude),independent_amplitude_relative_difference=float(abs(amplitude/cross-1))))
    dw=moments[-1]['delta_w_positive_mouth'];amplitude=moments[-1]['source_amplitude']
    cases=[_center_even_response(alpha,omega2,collar,level) for level,collar in
        zip((1,2),invariant['boundary_map_refinement'])]
    names=('inner_slope','inner_source_scaled','inner_from_outer','inner_from_source_scaled',
        'outer_slope','outer_source_scaled','center_from_outer','center_from_source_scaled')
    differences={name:abs(cases[0][name]/cases[1][name]-1) for name in names}
    assert max(differences.values())<2e-7
    algebra=_center_rn_algebra();r=s.symbols('r',positive=True)
    decay=s.sympify(algebra['decaying'],locals={'r':r})
    R=100.;F=(1-1/R)**2;B=math.sqrt(F);L=3*R-2
    dm=float(decay.subs(r,R));kappa=float((-((r-1)/r)**2*s.diff(decay,r)/decay).subs(r,R))
    hp=-4*R/(3*L);pip=8*F/(3*L**2);even=cases[-1];cg=-2/separation**2
    den=B*even['outer_slope']-2*F/(R*L)+kappa
    c_class=3/dm/den;c_source=(-alpha/(3*clock**2)*(pip+kappa*hp+2*R/L)
        -B*alpha**2*even['outer_source_scaled'])/den
    outer=dict(radius=R,incident_growing_amplitude=cg,kappa=kappa,denominator=den,
        P_outer_per_cg=c_class,P_outer_per_delta_w=c_source,P_outer_per_Xi=1/den,
        conditional_classical_P_outer=c_class*cg,conditional_state_P_outer=c_source*dw,
        conditional_classical_P_center=even['center_from_outer']*c_class*cg,
        direct_center_state_P=alpha**2*dw*even['center_from_source_scaled'],
        relation='Pouter=c_class*c_g+c_source*delta_w+Xi_ext/den. The quantum source coefficient is conditional on an RN approximation; unknown Xi_ext can contribute at that same order. No physical exterior or moving interface is selected.')
    omissions=_rn_omission_diagnostic(alpha,omega2,clock,dw,cg)
    assert abs(moments[0]['source_amplitude']/amplitude-1)<1e-10
    return dict(model='Reflection-even center condition and full-flux source-affine RN overlap',
        parameters=dict(alpha=alpha,omega_squared=omega2,separation_over_re=separation,clock_scale=clock,
            delta_w_positive_mouth=dw,source_amplitude=amplitude),source_projection=moments,refinement=cases,
        rn_algebra=algebra,conditional_outer_map=outer,rn_omission_diagnostic=omissions,
        verification=dict(all_declared_checks_pass=True,coefficient_relative_refinement=differences,
            amplitude_relative_refinement=abs(moments[0]['source_amplitude']/amplitude-1)),
        scope='Smooth reflection-even static sector only, fixed internal flux labels. The source orientation is Y=1-nu, opposite the saved mu=nu-1 convention. RN outer map is conditional and retains Xi_ext. Quantum polarization, complete-loop feedback, moving-cut matching, formation, transport and dynamical stability remain unsolved. No calendar-time travel or physical new discovery is established.',
        sources=['https://arxiv.org/html/1807.04726v3: assumed magnetic/Casimir channel framework; published solution itself does not violate ambient causality',
            'https://arxiv.org/abs/gr-qc/9507012: RN polar dipole framework',
            'Our reflection-sector composition, full reference-state projection and source-affine map, checked against existing component-verified equations'])


def save_center_exterior_benchmark():
    result=center_exterior_checks()
    project_path('polyakov_center_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_CENTER_EXTERIOR


# BEGIN_OVERLAP_TRANSPORT
def _overlap_asymptotic_checks():
    """First-order outer expansion and its projected RN-basis cancellation."""
    import sympy as s
    r,k=s.symbols('r k',positive=True)
    br=(r-1)/r;L=3*r-2
    mu=-s.Rational(1,6)+1/r**2-s.Rational(5,6)/r**3
    b1=-mu/(r-1)
    a1=mu/(r-1)**2-1/(2*r**4*(r-1))
    B=br+k*b1;A=1/r**2+k*a1;Ap=B*s.diff(A,r);Bp=B*s.diff(B,r)
    equations=[B*B+2*r*A*B-1+1/r**2+k*A*A,
        Ap-(1/r**4-(Bp+A*B)/r-A*A),
        Bp-(r/(r*r-k)*(1-1/r**2-k/r**4-B*B)-A*B)]
    residuals=[s.factor(s.diff(value,k).subs(k,0)) for value in equations]
    assert residuals==[0]*3
    raw,_,_=_invariant_collar_matrix(r,A,B,s.Integer(0),Ap,Bp,k)
    matrix=s.Matrix(raw)
    zeroth=matrix.subs(k,0).applyfunc(s.factor)
    first=matrix.diff(k).subs(k,0).applyfunc(s.factor)
    pulled=(first-b1/br*zeroth).applyfunc(s.factor)
    g=3*r*(r-1)**2/L
    decay=r*(30*r*r-45*r+17)/(10*(r-1)**3*L)
    x=lambda h:br*(s.diff(h,r)+2*h/(r*L))
    basis=s.Matrix([[g,decay],[x(g),x(decay)]])
    inverse=br*s.Matrix([[-x(decay),decay],[x(g),-g]])/3
    inverse_residual=(inverse*basis-s.eye(2)).applyfunc(s.factor)
    assert inverse_residual==s.zeros(2)
    omitted=(pulled*basis[:,0]).applyfunc(s.factor)
    radial=(inverse*omitted/br).applyfunc(s.factor)
    expected=-(465*r**5-1722*r**4+2066*r**3-944*r**2+100*r+20)/(15*r**3*(r-1)**2*L**3)
    projection_residual=s.factor(radial[0]-expected)
    leading_omitted=[s.limit(omitted[0],r,s.oo),s.limit(r*omitted[1],r,s.oo)]
    leading_parts=[s.limit(r*r*inverse[0,i]*omitted[i]/br,r,s.oo) for i in range(2)]
    growing_limit=s.limit(r**3*radial[0],r,s.oo)
    decaying_limit=s.limit(radial[1]/r,r,s.oo)
    assert projection_residual==0 and leading_omitted==[-s.Rational(4,9),s.Rational(4,9)]
    assert leading_parts==[-s.Rational(4,27),s.Rational(4,27)]
    assert growing_limit==-s.Rational(31,27) and decaying_limit==-s.Rational(4,9)
    return dict(background=dict(mu=str(mu),mass='1+k*mu+O(k²)',
            B_first_order=str(B),A_first_order=str(A),k='alpha/3'),
        background_equation_residuals=list(map(str,residuals)),
        inverse_residual=[[str(v) for v in row] for row in inverse_residual.tolist()],
        growing_projection_residual=str(projection_residual),
        pulled_proper_mismatch_per_k_cg=list(map(str,omitted)),
        radial_growing_projection_per_k_cg=str(expected),
        radial_decaying_projection_per_k_cg=str(radial[1]),
        large_radius=dict(proper_components_scaled_1_and_r=list(map(str,leading_omitted)),
            growing_adjoint_terms_times_r_squared=list(map(str,leading_parts)),
            leading_cancellation=str(sum(leading_parts)),
            growing_projection_times_r_cubed=str(growing_limit),
            decaying_projection_divided_by_r=str(decaying_limit)),
        comparison_points=[dict(radius=value,growing_projection_per_k_cg=float(expected.subs(r,value)))
            for value in (10,30,100)],all_declared_checks_pass=True,
        convention='Growing-only classical test input; projection coefficients are derivatives with respect to r, divided by k*c_g. The proper mismatch includes the B/sqrt(F) radial pullback. This checks the homogeneous operator correction, not the full forced transport.',
        scope='First order in k at fixed r>1, with d=w/N²=O(alpha²) omitted at this order and mean mass calibrated to the throat. The displayed large-r limits describe these expansion coefficients only; the approximation is neither near-throat nor infinity uniform. Running decaying-basis coefficients may grow while their physical profile correction remains subleading. No exterior matching or stability claim.')


def _rn_overlap_basis(r):
    """RN coordinates on the existing (P,X) state space, not junction data."""
    import mpmath as mp
    B=1-1/r;F=B*B;L=3*r-2
    g=3*r*(r-1)**2/L;d=r*(30*r*r-45*r+17)/(10*(r-1)**3*L)
    gp=g*(1/r+2/(r-1)-3/L)
    dp=d*(1/r+(60*r-45)/(30*r*r-45*r+17)-3/(r-1)-3/L)
    xg=B*(gp+2*g/(r*L));xd=B*(dp+2*d/(r*L))
    basis=mp.matrix([[g,d],[xg,xd]])
    inverse=B/3*mp.matrix([[-xd,d],[xg,-g]])
    particular=mp.matrix([-4*r/(3*L),2*r/(L*B)])
    matrix=mp.matrix([[-2*B/(r*L),1],[6/(r*L),-1/r**2+2*B/(r*L)]])
    forcing=mp.matrix([2*r*B/L,0])
    return basis,inverse,particular,matrix,forcing,B


def _overlap_transport_record(alpha,omega2,clock,cg,amplitude,refinement):
    """Integrate small corrections directly, avoiding subtraction of unit transfers."""
    import mpmath as mp
    from scipy.interpolate import CubicSpline
    from scipy.optimize import brentq
    background,_,_=_static_collar_solution(alpha,omega2,2*refinement)
    end=float(background.t[-1])
    at=lambda radius:brentq(lambda ell:1+alpha*background.sol(ell)[0]-radius,0,end,xtol=1e-12)
    start,mid=at(10.),at(30.)
    maximum_identity_error=0.
    with mp.workdps(70):
        aa,cc=mp.mpf(str(alpha)),mp.mpf(str(clock))
        grid=np.linspace(start,end,512*refinement+1)
        jets=[_quantum_overlap_jet(background.sol(ell),alpha,omega2) for ell in grid]
        # Recover delta=log[N/(clock*B_RN)] from its small derivative rather
        # than subtracting O(20) logarithms. The existing67 clock is N100/B100.
        slopes=np.array([float((j['A']-j['B']/(j['r']*(j['r']-1)))/aa) for j in jets])
        primitive=CubicSpline(grid,slopes).antiderivative()
        last=jets[-1];anchor=float(mp.log(last['B']/(1-1/last['r']))/aa)
        scaled_delta=anchor+primitive(grid)-primitive(end)
        clock_offset=float(mp.log(last['N']/(cc*last['B'])))
        clock_error=max(float(abs(mp.log(j['N']/(cc*(1-1/j['r'])))-aa*mp.mpf(str(value))))
            for j,value in zip(jets,scaled_delta))
        assert clock_error<5e-12
        def coefficients(jet,delta):
            nonlocal maximum_identity_error
            phi,inv,yp,mrn,brn,br=_rn_overlap_basis(jet['r'])
            Q=jet['M']-jet['B']/br*mrn
            E=inv*Q*phi
            clock_over_N_squared=mp.exp(-2*aa*mp.mpf(str(delta)))/br**2
            z=inv*(Q*yp-jet['b']*clock_over_N_squared+jet['B']/br*brn/br**2)
            zclock=inv*jet['b']*clock_over_N_squared
            error=max(abs(v) for v in inv*phi-mp.eye(2))
            maximum_identity_error=max(maximum_identity_error,float(error))
            return np.array(E.tolist(),float),np.array(list(z),float),np.array(list(zclock),float)
        sampled=np.array([np.r_[E.ravel(),z,zc] for E,z,zc in map(coefficients,jets,scaled_delta)])
        interpolant=CubicSpline(grid,sampled,axis=0)
        def rhs(ell,values):
            local=interpolant(ell);E=local[:4].reshape(2,2);z=local[4:6];zc=local[6:]
            D=values[:4].reshape(2,2);a=values[4:6]
            return np.r_[(E@D+E/alpha).ravel(),E@a+z/alpha,
                (E/alpha).ravel(),z/alpha,np.trace(E),E@values[13:15]+zc,zc]
        def evolve(lo,hi):
            sol=solve_ivp(rhs,(lo,hi),np.zeros(17),method='DOP853',
                rtol=1e-10/refinement,atol=1e-13/refinement,max_step=2/refinement,dense_output=True)
            assert sol.success
            return sol
        whole=evolve(start,end);right=evolve(mid,end)
        def unpack(values):return alpha*values[:4].reshape(2,2),alpha*values[4:6]
        D,a=unpack(whole.y[:,-1]);Dl,al=unpack(whole.sol(mid));Dr,ar=unpack(right.y[:,-1])
        composed=Dr+Dl+Dr@Dl;affine_composed=ar+al+Dr@al
        def rel(error,reference):return float(np.max(abs(error))/max(np.max(abs(reference)),1e-100))
        composition_error=rel(composed-D,D);source_error=rel(affine_composed-a,a)
        aclock=whole.y[13:15,-1]
        clock_composition_error=rel(right.y[13:15,-1]+(np.eye(2)+Dr)@whole.sol(mid)[13:15]-aclock,aclock)
        assert composition_error<2e-7 and source_error<2e-7 and maximum_identity_error<1e-50
        assert clock_composition_error<2e-7
        # Compare the small determinant correction without adding it to one.
        det_minus_one=np.trace(D)+np.linalg.det(D)
        determinant_error=abs(det_minus_one-math.expm1(whole.y[12,-1]))/max(abs(det_minus_one),1e-100)
        assert determinant_error<2e-7
        bornD=alpha*whole.y[6:10,-1].reshape(2,2);borna=alpha*whole.y[10:12,-1]
        rows=[]
        for ell in np.linspace(start,end,121):
            delta,affine=unpack(whole.sol(ell));r=1+alpha*background.sol(ell)[0]
            phi,inv,yp,*_=_rn_overlap_basis(mp.mpf(str(r)))
            rows.append(dict(radius=float(r),proper_distance=float(ell),
                growing_shift_from_incident=float(delta[0,0]*cg),
                growing_shift_from_incident_over_A=float(delta[0,0]*cg/amplitude),
                growing_shift_from_source=float(amplitude*affine[0]),
                RN_particular_growing_projection_per_A=float((inv*yp)[0])))
    return dict(refinement=refinement,radius_interval=[10.,100.],proper_interval=[start,end],
        transfer_correction=D.tolist(),affine_correction_per_A=a.tolist(),
        first_order_transfer_correction=bornD.tolist(),first_order_affine_correction_per_A=borna.tolist(),
        growing_shift_from_incident=D[0,0]*cg,growing_shift_from_source=amplitude*a[0],
        source_order_incident_growing_shift=D[0,0]*cg/amplitude,
        growing_shift_per_unknown_decaying_input=D[0,1],
        composition_relative_error=composition_error,affine_composition_relative_error=source_error,
        determinant_relative_error=determinant_error,basis_inverse_maximum_error=maximum_identity_error,
        beyond_first_order_transfer_relative_size=rel(D-bornD,D),
        beyond_first_order_affine_relative_size=rel(a-borna,a),solver_evaluations=whole.nfev+right.nfev,
        coefficient_sample_count=len(grid),
        reconstructed_clock_log_maximum_difference=clock_error,
        background_clock_offset_estimate=clock_offset,
        affine_clock_offset_response_per_A=aclock.tolist(),
        clock_composition_relative_error=clock_composition_error,
        plot=rows,
        convention='Y=Phi_RN*c+A_dip*y_particular. c100=(I+D)c10+A_dip*[a+(1-exp(-2epsilon))*a_clock], epsilon=log[N100/(clock*Bfull100)]. The zero-offset column follows the existing clock calibration; epsilon is retained explicitly, not silently selected from rounded logs. D and a are integrated after division byalpha. Growing/source columns are diagnostics, not selected physical inner data.')


def overlap_transport_checks():
    import sympy as s
    from scipy.integrate import quad
    ref=json.loads(project_path('polyakov_center_results.json').read_text(encoding='utf-8'))
    p=ref['parameters'];alpha=p['alpha'];omega2=p['omega_squared'];clock=p['clock_scale']
    cg=ref['conditional_outer_map']['incident_growing_amplitude'];amplitude=p['source_amplitude']
    cases=[_overlap_transport_record(alpha,omega2,clock,cg,amplitude,level) for level in (1,2)]
    errors={}
    for name in ('transfer_correction','affine_correction_per_A','affine_clock_offset_response_per_A'):
        a,b=[np.asarray(c[name]) for c in cases]
        errors[name]=float(np.max(abs((a-b)/b)))
    assert max(errors.values())<2e-6
    asymptotic=_overlap_asymptotic_checks();radius=s.symbols('r',positive=True)
    projection=s.sympify(asymptotic['radial_growing_projection_per_k_cg'],locals={'r':radius})
    predicted=alpha/3*quad(s.lambdify(radius,projection,'numpy'),10,100,epsabs=1e-14,epsrel=1e-12)[0]
    prediction_error=abs(cases[-1]['transfer_correction'][0][0]/predicted-1)
    assert prediction_error<2e-6
    return dict(model='Full quantum-collar transport in an affine RN basis across the outer overlap',
        parameters=p,incident_growing_amplitude=cg,refinement=cases,asymptotic_checks=asymptotic,
        equations=dict(basis='Y=(P,X)=Phi_RN*c+A_dip*yp; Phi columns are (g,Xg),(d,Xd).',
            flow='dc/dl=E*c+A_dip*z; Q=Mfull-(Bfull/B_RN)*M_RN; E=Phi_inverse*Q*Phi.',
            source='z=Phi_inverse*(Q*yp-bfull*clock²/N²+(Bfull/B_RN)*bRN/F). Includes source-affine momentum and radial pullback.',
            clock='delta0_l=A-Bfull/[r(r-1)], delta0_100=log(Bfull100/B_RN100); actual delta=delta0+epsilon, epsilon=log[N100/(clock*Bfull100)]. Integrate delta0/alpha. Source correction is z0+(1-exp(-2epsilon))*Phi_inverse*bfull*exp(-2delta0)/F.',
            inverse='det(Phi)=-3/B_RN; growing adjoint row=B_RN/3*(-Xd,d).',
            boundary='c_g100=(1+Dgg)*c_g10+Dgd*c_d10+A_dip*[a_g+(1-exp(-2epsilon))*a_clock,g]. Outer c_g/c_d, clock offset and nonradial return state are not selected.'),
        verification=dict(all_declared_checks_pass=True,relative_refinement=errors,
            independent_first_order_growing_transfer=predicted,
            full_growing_transfer_relative_difference_from_first_order=prediction_error),
        scope='R10 toR100 only, inside the existing spherical mean-background collar and linear dipole model. This transports all retained local quantum terms, not the nonspherical far exterior or its complete state. RN functions are a coordinate basis on (P,X), not a physical electrovac junction. No global matching, dynamical stability, manufacture or13Sep2027->12Sep2027 passage is established.',
        sources=['https://arxiv.org/html/1807.04726v3: assumed magnetic-channel framework',
            'https://arxiv.org/abs/gr-qc/9507012: RN polar dipole framework',
            'Our affine variation-of-constants calculation, reusing the previously component-checked collar equations'])


def save_overlap_transport_benchmark():
    result=overlap_transport_checks()
    project_path('polyakov_overlap_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_OVERLAP_TRANSPORT


# BEGIN_RETURN_CURRENT_FIELD
def _return_current_kernel(nu, s, circumference):
    """Paired half-path Biot--Savart integrand; lengths in units of d."""
    h=1-s*s;root=np.hypot(h,nu*s);u=(h+root)/nu
    um=(2-nu)*h/(nu*(u+1)-2*h)
    C=u*u*(1+3*s*s)-s*s+s**4
    return -um*(u+1)*u*u*(1+s*s)*(u*u+3*s*s)/(4*circumference**2*C**1.5*(u*u+s*s))


def return_current_field_checks():
    from scipy.integrate import quad
    from mmp_pulse import _magnetic_loop_source,_magnetic_return_path,_magnetic_return_length
    path=project_path('polyakov_closure_results.json');closure=json.loads(path.read_text(encoding='utf-8'))
    center=json.loads(project_path('polyakov_center_results.json').read_text(encoding='utf-8'))
    alpha=closure['parameters']['alpha'];distance=closure['parameters']['separation_over_re']
    records=[]
    for saved in closure['quadrature_refinement']:
        state=saved['roots'][-1];profile=state['profile'];clock=state['clock_scale']
        nu=np.asarray(profile['nu']);weights=np.asarray(profile['weights'])
        length=clock*np.asarray(profile['periods'])/distance
        assert abs(weights.sum()-1)<1e-14 and np.all((nu>0)&(nu<2))
        assert np.allclose((2*math.pi/length)**2,
            np.asarray(profile['channel_omega_squared'])*(distance/clock)**2,rtol=1e-14)
        # f(nu)>=1, positive Shapiro delay; this does not bound unknown actual channels.
        lower_length=clock*state['throat_optical_time']/distance+1-2*(state['cutoff']-1)/distance
        assert lower_length>0 and np.all(length>=lower_length)
        for order in (saved['order'],2*saved['order']):
            nodes,w=np.polynomial.legendre.leggauss(order);s=(nodes+1)/2;w=w/2
            kernel=_return_current_kernel(nu[:,None],s[None,:],length[:,None])
            flux_integrand=kernel@w;terms=weights*flux_integrand;integral=float(terms.sum())
            assert np.all(kernel<0)
            records.append(dict(flux_order=saved['order'],path_order=order,integral=integral,
                lower_circumference=lower_length,nu=nu.tolist(),weights=weights.tolist(),
                circumference=length.tolist(),flux_integrand=flux_integrand.tolist(),
                weighted_contribution=terms.tolist(),cumulative_integral=np.cumsum(terms).tolist()))
    fine=records[-1];integral=fine['integral']
    refinement=max(abs(row['integral']/integral-1) for row in records)
    assert refinement<1e-9

    # Independent Cartesian field, derivative, force/current and volume Jacobian.
    direct=records[0];nodes,w=np.polynomial.legendre.leggauss(direct['path_order']);nodes=(nodes+1)/2;w=w/2
    flux_direct=[];kernel_errors=[];jacobian_errors=[];current_errors=[]
    for nu,length in zip(direct['nu'],direct['circumference']):
        total=0.
        for s,weight in zip(nodes,w):
            rho,z,speed=_magnetic_return_path(nu,float(s))
            row=_magnetic_loop_source([rho,0.,z],0.,circumference=length)
            rp=math.hypot(rho,z-.5);rm=math.hypot(rho,z+.5)
            value=speed/row['magnitude']*row['current_physical'][1]*rho*(rp**-3+rm**-3)
            analytic=float(_return_current_kernel(nu,s,length))
            kernel_errors.append(abs(value/analytic-1));total+=weight*value
            h=1-s*s;root=math.hypot(h,nu*s);u=(h+root)/nu
            volume=(u*u-s*s)**2/(16*root)
            jacobian_errors.append(abs(speed/row['magnitude']/volume-1))
            expected=-np.cross(row['tangent'],row['curvature'])/(12*length*length)
            current_errors.append(float(np.linalg.norm(row['current_physical']-expected)/np.linalg.norm(expected)))
        flux_direct.append(total)
    direct_integral=float(np.dot(direct['weights'],flux_direct))
    direct_error=abs(direct_integral/direct['integral']-1)
    # Cartesian cross products lose digits near the axis. Resolve the worst
    # point independently at high precision; retain its float error as data.
    import mpmath as mp
    i,j=divmod(int(np.argmax(kernel_errors)),len(nodes))
    nu=direct['nu'][i];s=float(nodes[j]);length=direct['circumference'][i]
    with mp.workdps(50):
        v,t,L=map(lambda x:mp.mpf(str(x)),(nu,s,length))
        h=1-t*t;root=mp.sqrt(h*h+v*v*t*t);u=(h+root)/v
        rho=mp.sqrt((u*u-1)*h)/2;z=u*t/2
        plus=mp.matrix([rho,0,z-mp.mpf('.5')]);minus=mp.matrix([rho,0,z+mp.mpf('.5')])
        rp=mp.norm(plus);rm=mp.norm(minus);b=plus/rp**3-minus/rm**3;mag=mp.norm(b)
        db=(mp.eye(3)/rp**3-3*plus*plus.T/rp**5)-(mp.eye(3)/rm**3-3*minus*minus.T/rm**5)
        tangent=b/mag;curvature=(mp.eye(3)-tangent*tangent.T)*db*tangent/mag
        current=-(tangent[2]*curvature[0]-tangent[0]*curvature[2])/(12*L*L)
        du=t*(v-2*u)/root
        drho=(u*du*h-t*(u*u-1))/(4*rho);dz=(u+t*du)/2
        speed=mp.sqrt(drho*drho+dz*dz)
        high_value=speed/mag*current*rho*(rp**-3+rm**-3)
        high_error=float(abs(high_value/mp.mpf(str(float(_return_current_kernel(nu,s,length))))-1))
    assert high_error<1e-11 and max(jacobian_errors)<1e-8
    assert max(current_errors)<1e-12 and direct_error<1e-10
    compatibility=[]
    for nu in (.1,1.,1.9):
        rho,z,_=_magnetic_return_path(nu,.3)
        old=_magnetic_loop_source([rho,0.,z],3.)
        explicit=_magnetic_loop_source([rho,0.,z],3.,circumference=old['circumference'])
        compatibility.append(bool(np.array_equal(old['current_physical'],explicit['current_physical'])))
    assert all(compatibility)

    # Endpoint limits independently obtained from the prolate expression.
    far_length=2+math.asinh(math.sqrt(3))/math.sqrt(3)
    low_slope=-quad(lambda s:(1-s*s)*(1+s*s)/(2*far_length**2*(1+3*s*s)**1.5),0,1,epsabs=1e-14)[0]
    state=closure['quadrature_refinement'][-1]['roots'][-1];clock=state['clock_scale']
    axis_length=state['shortest_period']*clock/distance
    high_slope=-(math.pi/32+1/6)/(2*axis_length**2)
    endpoint=[]
    for epsilon in (1e-3,1e-4,1e-5):
        low=quad(lambda s:float(_return_current_kernel(epsilon,s,far_length/epsilon)),0,1,epsabs=1e-13)[0]/epsilon
        high=quad(lambda s:float(_return_current_kernel(2-epsilon,s,axis_length)),0,1,epsabs=1e-13)[0]/epsilon
        endpoint.append(dict(epsilon=epsilon,low_flux_slope=low,high_flux_slope=high,
            low_relative_error=abs(low/low_slope-1),high_relative_error=abs(high/high_slope-1)))
    assert max(endpoint[-1]['low_relative_error'],endpoint[-1]['high_relative_error'])<2e-5
    # Strict bounds ONLY on the extrapolated flat point-current in excised balls.
    epsilon=(state['cutoff']-1)/distance
    wmax=(2*math.pi/fine['lower_circumference'])**2
    factor=(1+2*epsilon+2*epsilon**2)/(1-2*epsilon)**3
    upper_bound=wmax*factor*epsilon**2/(48*math.pi**2)
    lower_bound=wmax*factor*epsilon**5/(120*math.pi**2*(1-epsilon)**3)
    nu=np.asarray(fine['nu']);weights=np.asarray(fine['weights']);W=(2*math.pi/np.asarray(fine['circumference']))**2
    w0=float(weights@W);w2=float(5*weights@(W*(3*(1-nu)**2-1)/2))
    near_coefficient=w0-w2/5
    magnetic=16*math.pi**2*alpha*integral/distance**2
    incident=-2/distance**2;amplitude=center['parameters']['source_amplitude']
    return dict(model='Free-space Biot-Savart contribution from prescribed returning-channel current',
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),parameters=closure['parameters'],
        refinement=records,integral=integral,magnetic_growing_contribution=magnetic,
        classical_companion_growing_input=incident,relative_to_classical_input=magnetic/incident,
        source_amplitude=amplitude,magnetic_contribution_over_A=magnetic/amplitude,
        endpoint=dict(low_flux_slope=low_slope,high_flux_slope=high_slope,checks=endpoint,
            small_nu_integral_coefficient=low_slope/4,
            limitation='These are formal continuum-channel endpoint limits. At finite q, nu*q of order one loses localization; its unmodeled QFT contribution is not bounded by the formal tail.'),
        mouth_excision=dict(isotropic_cutoff=state['cutoff']-1,cutoff_over_separation=epsilon,
            Wbar_upper_bound=wmax,upper_ball_absolute_bound=upper_bound,lower_ball_absolute_bound=lower_bound,
            relative_bound_to_integral=(upper_bound+lower_bound)/abs(integral),
            Wbar_monopole=w0,Wbar_quadrupole=w2,local_dipole_current_coefficient=near_coefficient,
            upper_ball_leading_estimate=-epsilon**2*near_coefficient/(48*math.pi**2),
            scope='Bounds remove only the assumed flat point-current. Unknown throat/interface and curved-region currents have no bound here.'),
        equations=dict(current='jhat_phi=-(t cross curvature)_phi/(12 Lbar(nu)^2), unit tangent t; jhat=j*d^3.',
            volume='dV/d^3=2*pi*(speed/|b|)*dnu*ds. Flux weights integrate dnu/2.',
            kernel='G=-integral_0^1 [(u^2-1)u^2(1+s^2)(u^2+3s^2)/(4 Lbar^2 C^(3/2)(u^2+s^2))] ds; C=u^2(1+3s^2)-s^2+s^4.',
            normalization='I=integral_0^2 G(nu)dnu/2; delta Bphysical_z=g^2 I/d^2; delta cg_magnetic=16*pi^2*alpha*I/d^2.',
            additive_data='cg_magnetic=classical_companion+formal_volume_current+unknown_correction_where_the_model_fails/interface. Replacing an invalid modeled region requires actual-minus-modeled current, not double counting. The coupled metric/master boundary data are still unselected.'),
        verification=dict(all_declared_checks_pass=True,relative_quadrature_change=refinement,
            independent_cartesian_integral=direct_integral,cartesian_integral_relative_error=direct_error,
            maximum_kernel_relative_error=max(kernel_errors),maximum_volume_jacobian_relative_error=max(jacobian_errors),
            worst_cartesian_point_high_precision_relative_error=high_error,
            maximum_current_identity_relative_error=max(current_errors),default_source_compatibility=all(compatibility)),
        scope='Leading static flat-region current on prescribed reference geometry and current67 periods. Source-generated particular only, without artificial patch clamping. This does not solve the coupled nonspherical metric, channel-length susceptibility, unknown finite-q endpoint state, interface fields, formation, dynamical stability or time travel. The source extends to infinity; a vacuum dipole far boundary cannot be assumed.',
        sources=['https://arxiv.org/html/1807.04726v3: magnetic channels and finite-q localization, Appendix D/E. The paper describes a causality-preserving long wormhole.',
            'Our current from the assumed channel action, compact Green functional and numerical checks; not a published or experimentally verified discovery.'])


def save_return_current_benchmark():
    result=return_current_field_checks()
    project_path('polyakov_return_current_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_RETURN_CURRENT_FIELD


# BEGIN_RETURN_PATH_FEEDBACK
def _return_channel_coefficients(order):
    saved=json.loads(project_path('polyakov_closure_results.json').read_text(encoding='utf-8'))
    row=next(r for r in saved['quadrature_refinement'] if r['order']==order)['roots'][-1]
    nu=np.asarray(row['profile']['nu']);distance=saved['parameters']['separation_over_re']
    W=np.asarray(row['profile']['channel_omega_squared'])*(distance/row['clock_scale'])**2
    # Divide the known quadratic endpoint zero before interpolating.
    coefficients=np.polynomial.legendre.legfit(nu-1,W/nu**2,order-1)
    return coefficients


def _return_ring_kernel(rho,z,source_rho,source_z):
    """Azimuthal Green integral divided by 4; elliptic arguments are parameters."""
    from scipy.special import ellipkm1,ellipe,hyp2f1
    D=(rho+source_rho)**2+(z-source_z)**2
    m=4*rho*source_rho/D
    complement=((rho-source_rho)**2+(z-source_z)**2)/D
    small=m<1e-3
    H=np.empty_like(m)
    H[small]=math.pi*m[small]/16*hyp2f1(1.5,1.5,3.,m[small])
    H[~small]=((1+complement[~small])*ellipkm1(complement[~small])-2*ellipe(m[~small]))/m[~small]
    assert np.all(complement>0) and np.all(H>0)
    return H/np.sqrt(D)


def _return_current_stream(nu,s,coefficients,order=32):
    """Free-space w at a reference path point; delta nu=8*pi²*alpha*w."""
    from mmp_pulse import _magnetic_return_path
    rho,z,_=_magnetic_return_path(nu,s)
    nodes,weights=np.polynomial.legendre.leggauss(order);p=(nodes[:,None]+1)/2;t=(nodes[None,:]+1)/2
    weights=weights[:,None]*weights[None,:]/4
    value=0.
    # Eight Duffy triangles cover the full source domain 0<nu'<2,-1<s'<1.
    # The self-ring logarithm lies at p=0, never at a quadrature node.
    for side in (-1,1):
        for b in (-1-s,1-s):
            for swap in (False,True):
                x=p*(t if swap else np.ones_like(t))
                if side<0:
                    v=nu*(1-x);dv=nu
                else:
                    c=1-nu/2;v=nu/(1-c*x);dv=nu*c/(1-c*x)**2
                ss=s+b*p*(np.ones_like(t) if swap else t)
                h=1-ss*ss;root=np.hypot(h,v*ss);u=(h+root)/v
                um=(2-v)*h/(v*(u+1)-2*h)
                rr=np.sqrt(um*(u+1)*h)/2;zz=u*ss/2
                C=u*u*(1+3*ss*ss)-ss*ss+ss**4
                W=v*v*np.polynomial.legendre.legval(v-1,coefficients)
                assert np.all(W>0)
                jvol=-W*rr*u*(1+ss*ss)*(u*u-ss*ss)**2/(64*math.pi**2*C**1.5*root)
                green=_return_ring_kernel(rho,z,rr,zz)
                value+=float(np.sum(weights*abs(b)*dv*p*jvol*green))
    return rho*value/math.pi
def _return_far_limit_checks():
    """Formal long-loop tail of the selected flat free-space current response."""
    import sympy as s
    from scipy.integrate import quad
    mu,z=s.symbols('mu z',real=True);f=s.Function('f')
    C=2+math.asinh(math.sqrt(3))/math.sqrt(3)
    length_quad=quad(lambda value:2*math.sqrt(1+3*value*value),0,1,epsabs=1e-14)[0]
    forcing=lambda value:(1-value*value)**2*(1+value*value)/(4*C*C*(1+3*value*value)**1.5)
    # The even Dirichlet solution; each integration weight is nonnegative.
    green=-(1-mu)*s.Integral(f(z),(z,0,mu))-s.Integral((1-z)*f(z),(z,mu,1))
    ode=s.simplify(s.diff(green,mu,2)-f(mu))
    symmetry=s.simplify(s.diff(green,mu).subs(mu,0).doit())
    boundary=s.simplify(green.subs(mu,1).doit())
    c=s.symbols('C',positive=True)
    jtail=-(1-mu*mu)**s.Rational(5,2)*(1+mu*mu)/(4*c*c*(1+3*mu*mu)**s.Rational(3,2))
    exact_forcing=(1-mu*mu)**2*(1+mu*mu)/(4*c*c*(1+3*mu*mu)**s.Rational(3,2))
    maxwell=s.simplify((1-mu*mu)*exact_forcing+s.sqrt(1-mu*mu)*jtail)
    curvature=3*s.sqrt(1-mu*mu)*(1+mu*mu)/(1+3*mu*mu)**s.Rational(3,2)
    current=s.simplify(jtail+4*s.pi**2*(1-mu*mu)**2/c**2*curvature/(48*s.pi**2))
    assert [ode,symmetry,boundary,maxwell,current]==[0]*5
    def stream(value):
        value=abs(float(value))
        first=quad(forcing,0,value,epsabs=1e-14,epsrel=1e-12)[0]
        second=quad(lambda t:(1-t)*forcing(t),value,1,epsabs=1e-14,epsrel=1e-12)[0]
        return -(1-value)*first-second
    midplane=stream(0.)
    weight=lambda value:(1-value**4)/(1+3*value*value)**1.5
    coefficient=6*quad(lambda value:stream(value)*weight(value),0,1,epsabs=1e-14,epsrel=1e-12)[0]
    # Independent IVP evaluation of the same boundary-value solution and length.
    sol=solve_ivp(lambda value,y:[y[1],forcing(value)],(0,1),[midplane,0.],
        method='DOP853',rtol=1e-12,atol=1e-14,max_step=.025,dense_output=True)
    assert sol.success
    ivp_coefficient=6*quad(lambda value:sol.sol(value)[0]*weight(value),0,1,epsabs=1e-14)[0]
    grid=np.linspace(-1,1,61);values=[stream(value) for value in grid]
    profile_error=max(abs(value-sol.sol(abs(position))[0]) for position,value in zip(grid,values))
    coefficient_error=abs(coefficient-ivp_coefficient)
    assert abs(C-length_quad)<2e-14 and profile_error<2e-13 and coefficient_error<2e-13
    assert values[0]==values[-1]==0 and max(values[1:-1])<0 and coefficient<0
    return dict(far_length_coefficient=C,far_midplane_stream=midplane,
        length_inverse_nu_squared_coefficient=coefficient,
        relative_length_inverse_nu_coefficient=coefficient/C,
        equations=dict(length='Lbar~C/nu, C=2+asinh(sqrt(3))/sqrt(3)',
            current='r³*jhat_phi=-(1-mu²)^(5/2)*(1+mu²)/[4C²(1+3mu²)^(3/2)], mu=cos(theta)',
            source_ode="w_infinity''(mu)=(1-mu²)²*(1+mu²)/[4C²(1+3mu²)^(3/2)]; w_infinity(-1)=w_infinity(1)=0",
            variation='lim_nu->0 nu²*deltaLbar/lambda = 6*integral_0^1 w_infinity(mu)*(1-mu^4)/(1+3mu²)^(3/2) dmu',
            normalization='Delta_*w=-rho*jhat_phi; delta_nu=lambda*w, lambda=8pi²*alpha. Unit tangent in jhat_phi=-(t cross curvature)_phi/(12Lbar²).'),
        symbolic_residuals=dict(green_ode=str(ode),even_derivative=str(symmetry),dirichlet=str(boundary),
            angular_maxwell=str(maxwell),current_normalization=str(current)),
        verification=dict(all_declared_checks_pass=True,length_quadrature_difference=abs(C-length_quad),
            independent_IVP_profile_difference=float(profile_error),independent_length_coefficient_difference=float(coefficient_error),
            independent_IVP_boundary_residual=float(abs(sol.y[0,-1]))),
        plot=[dict(mu=float(position),stream=value) for position,value in zip(grid,values)],
        scope='Formal continuum-channel far limit of the flat point-pair free-space response, not a finite-q completed exterior. Stream negativity gives this exterior shape contribution only. deltaLbar/Lbar grows as lambda/nu, so linearization is not uniform at nu=0; finite-q localization also fails there. First frequency-moment response remains integrable, but an unweighted mean deltaLbar does not. Throat, interface and metric feedback are not supplied.')


def _return_length_feedback(nu,coefficients,source_order,path_order):
    nodes,weights=np.polynomial.legendre.leggauss(path_order);s=(nodes+1)/2
    length_weight=_return_path_length_weight(nu,s)
    stream=np.array([_return_current_stream(nu,float(t),coefficients,source_order) for t in s])
    return float(weights@(stream*length_weight))  # Reflection doubles the half path.


def _return_path_length_weight(nu,s):
    h=1-s*s;root=np.hypot(h,nu*s);u=(h+root)/nu
    C=u*u*(1+3*s*s)-s*s+s**4
    return 3*u*(1+s*s)*(u*u-s*s)**2/(4*root*C**1.5)


def _return_feedback_record(flux_order,path_order,source_order):
    closure=json.loads(project_path('polyakov_closure_results.json').read_text(encoding='utf-8'))
    state=next(r for r in closure['quadrature_refinement'] if r['order']==flux_order)['roots'][-1]
    profile=state['profile'];nu=np.asarray(profile['nu']);weights=np.asarray(profile['weights'])
    distance=closure['parameters']['separation_over_re'];clock=state['clock_scale']
    length=clock*np.asarray(profile['periods'])/distance
    coefficients=_return_channel_coefficients(flux_order)
    variation=np.array([_return_length_feedback(float(v),coefficients,source_order,path_order) for v in nu])
    coupling=8*math.pi**2*closure['parameters']['alpha']
    fractional_length=coupling*variation/length
    change=-2*np.asarray(profile['channel_omega_squared'])*fractional_length
    assert np.all(variation<0) and np.all(change>0)
    return dict(flux_order=flux_order,path_order=path_order,source_order=source_order,
        nu=nu.tolist(),weights=weights.tolist(),length_response_per_coupling=variation.tolist(),
        fractional_period_change=fractional_length.tolist(),channel_omega_squared_change=change.tolist(),
        fixed_collar_mean_omega_squared_change=float(weights@change),
        largest_sampled_fractional_period_change=float(np.max(abs(fractional_length))),
        positive_mouth_dipole_change=float(3*weights@(change*(1-nu))))
def return_path_feedback_checks(records=None):
    from scipy.integrate import quad
    from mmp_pulse import _magnetic_loop_source,_magnetic_return_path
    if records is None:
        records=[_return_feedback_record(64,24,32),_return_feedback_record(128,48,64)]
    assert [(r['flux_order'],r['path_order'],r['source_order']) for r in records]==[(64,24,32),(128,48,64)]
    closure_path=project_path('polyakov_closure_results.json');closure=json.loads(closure_path.read_text(encoding='utf-8'))
    state=closure['quadrature_refinement'][-1]['roots'][-1];profile=state['profile']
    weights=np.asarray(profile['weights']);nu=np.asarray(profile['nu']);period=np.asarray(profile['periods'])
    alpha=closure['parameters']['alpha'];distance=closure['parameters']['separation_over_re']
    clock=state['clock_scale'];omega2=state['omega_squared'];fine=records[-1]
    for record in records:
        assert np.all(np.asarray(record['length_response_per_coupling'])<0)
        assert np.all(np.asarray(record['fractional_period_change'])<0)
        assert math.isclose(np.dot(record['weights'],record['channel_omega_squared_change']),
            record['fixed_collar_mean_omega_squared_change'],rel_tol=1e-14)
    changes={key:abs(records[0][key]/fine[key]-1) for key in
        ('fixed_collar_mean_omega_squared_change','positive_mouth_dipole_change')}
    assert max(changes.values())<2e-6
    coefficients=_return_channel_coefficients(128);coarse_coefficients=_return_channel_coefficients(64)
    testnu=np.r_[np.geomspace(1e-10,.1,100),np.linspace(.1,2,100)]
    interpolation_error=float(np.max(abs(np.polynomial.legendre.legval(testnu-1,coarse_coefficients)/
        np.polynomial.legendre.legval(testnu-1,coefficients)-1)))
    assert interpolation_error<3e-7
    angular_errors=[]
    for rho,z,rr,zz in ((1.,.1,.6,-.2),(.01,.2,2.,.7),(1.,0.,1.01,.01)):
        angular=quad(lambda phi:math.cos(phi)/math.sqrt(rho*rho+rr*rr-2*rho*rr*math.cos(phi)+(z-zz)**2),
            0,2*math.pi,epsabs=1e-12,epsrel=1e-12)[0]
        ring=4*float(_return_ring_kernel(rho,z,np.array(rr),np.array(zz)))
        angular_errors.append(abs(ring/angular-1))
    assert max(angular_errors)<1e-9
    length_errors=[]
    for v,s in ((.1,.3),(.5,.7),(1.,0.),(1.9,.8)):
        rho,z,speed=_magnetic_return_path(v,s)
        row=_magnetic_loop_source([rho,0,z],0.,circumference=1.)
        direct=np.cross(row['tangent'],row['curvature'])[1]*speed/(rho*row['magnitude'])
        h=1-s*s;root=math.hypot(h,v*s);u=(h+root)/v;C=u*u*(1+3*s*s)-s*s+s**4
        compact=3*u*(1+s*s)*(u*u-s*s)**2/(4*root*C**1.5)
        length_errors.append(abs(direct/compact-1))
    assert max(length_errors)<1e-10
    def cartesian_stream(rho,z):
        rp=math.hypot(rho,z-.5);rm=math.hypot(rho,z+.5)
        return _return_current_stream((z+.5)/rm-(z-.5)/rp,rm-rp,coefficients,128)
    rho,z=1.,.1;value=cartesian_stream(rho,z)
    v=(z+.5)/math.hypot(rho,z+.5)-(z-.5)/math.hypot(rho,z-.5)
    W=v*v*np.polynomial.legendre.legval(v-1,coefficients)
    rhs=-rho*_magnetic_loop_source([rho,0,z],0.,circumference=2*math.pi/math.sqrt(W))['current_physical'][1]
    lhs=[]
    for h in (.08,.04,.02):
        plus=cartesian_stream(rho+h,z);minus=cartesian_stream(rho-h,z)
        lhs.append((plus+minus+cartesian_stream(rho,z+h)+cartesian_stream(rho,z-h)-4*value)/h**2-(plus-minus)/(2*h*rho))
    pde_error=abs(((4*lhs[-1]-lhs[-2])/3)/rhs-1)
    assert pde_error<2e-6 and 3.8<abs((lhs[0]-rhs)/(lhs[1]-rhs))<4.2
    far=_return_far_limit_checks()
    far_sample=_return_current_stream(1e-5,0.,coefficients,128)
    assert abs(far_sample/far['far_midplane_stream']-1)<2e-4

    # Reuse the actual mean collar solution. The exterior reference lengths
    # remain fixed while the throat optical time and common clock respond.
    external=clock*(period-state['throat_optical_time']);derivatives=[]
    for step in (1e-4,5e-5):
        samples=[]
        for sign in (-1,1):
            newclock,data=_magnetic_closure_collar(omega2+sign*step,alpha)
            T=2*data['100.0']['optical_x']+external/newclock
            samples.append(dict(frequency=(2*math.pi/T)**2,clock=newclock,
                throat_time=2*data['100.0']['optical_x']))
        channel=(samples[1]['frequency']-samples[0]['frequency'])/(2*step)
        derivatives.append(dict(step=step,mean_map_derivative=float(weights@channel),
            channel_derivative=channel.tolist(),clock_log_derivative=(samples[1]['clock']-samples[0]['clock'])/(2*step*clock),
            throat_optical_time_derivative=(samples[1]['throat_time']-samples[0]['throat_time'])/(2*step)))
    derivative=derivatives[-1];slope=derivative['mean_map_derivative'];denominator=1-slope
    slope_error=abs(derivatives[0]['mean_map_derivative']/slope-1)
    assert slope_error<1e-7 and denominator>0
    forced=fine['fixed_collar_mean_omega_squared_change'];shift=forced/denominator
    total=np.asarray(fine['channel_omega_squared_change'])+shift*np.asarray(derivative['channel_derivative'])
    residual=abs(float(weights@total)/shift-1)
    assert residual<1e-12
    return dict(model='First magnetic path-length feedback with a separate mean-throat response',
        input_sha256=hashlib.sha256(closure_path.read_bytes()).hexdigest(),parameters=closure['parameters'],
        refinement=records,far_limit=far,
        mean_closure=dict(reference_omega_squared=omega2,derivative_refinement=derivatives,
            fixed_collar_forcing=forced,implicit_denominator=denominator,mean_omega_squared_shift=shift,
            retained_fraction_of_fixed_collar_forcing=1/denominator,
            clock_fractional_shift=derivative['clock_log_derivative']*shift,
            throat_optical_time_shift=derivative['throat_optical_time_derivative']*shift,
            channel_omega_squared_total_change=total.tolist(),
            positive_mouth_dipole_total_change=float(3*weights@(total*(1-nu))),
            scope='Only the existing radial mean collar adjusts. The nonspherical metric, interface and finite-q state do not. Fprime is a static closure-map derivative; its magnitude is not a dynamical stability criterion. This is a linear correction, not a completed nonlinear solution.'),
        verification=dict(all_declared_checks_pass=True,relative_moment_refinement=changes,
            channel_interpolation_relative_difference=interpolation_error,angular_Green_relative_error=max(angular_errors),
            path_weight_cartesian_relative_error=max(length_errors),maxwell_rhs=rhs,maxwell_finite_difference_lhs=lhs,
            maxwell_extrapolated_relative_error=pde_error,mean_map_derivative_relative_change=slope_error,
            linear_closure_relative_residual=residual,small_nu_midplane_stream=far_sample),
        equations=dict(stream='Delta_*w=-rho*jhat_phi; delta nu=lambda*w, lambda=8*pi²*alpha. w=rho/pi integral jhat*V*H(m)/sqrt(D) dnu ds.',
            length='deltaLbar/lambda=integral w*[3u(1+s²)(u²-s²)²/(4root*C^(3/2))] ds; point-mouth endpoint work vanishes.',
            state='delta omega_nu²=-2omega_nu²*deltaLbar/Lbar at fixed common clock and collar.',
            mean='delta w=delta F_magnetic/[1-Fprime(w)]; delta omega_nu²_total=delta omega_nu²_magnetic+(partial_w omega_nu²)*delta w.'),
        scope='Formal leading static flat-source Green response and its first path/state correction; same continuum channel assumptions as73. The longest-channel perturbation is nonuniform: deltaL/L~−0.00798909*lambda/nu. Finite-q channels near nu*q~1 require a different quantum description. Unknown actual-minus-modeled currents, curved-region and nonspherical metric/clock response, physical boundaries, formation and dynamics remain. No stable medium or backward travel is established.',
        sources=['https://arxiv.org/html/1807.04726v3: magnetic-channel model and its finite-q limit.',
            'https://arxiv.org/html/2510.11663v1: azimuthal Green representation; our noncompact source needs the separately checked convergent limit.',
            'https://dlmf.nist.gov/19.2: elliptic-integral definitions. Path feedback, tail and conditional mean update are our calculations.'])


def save_return_feedback_benchmark(records=None):
    result=return_path_feedback_checks(records)
    project_path('polyakov_feedback_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
# END_RETURN_PATH_FEEDBACK


# BEGIN_GRAVITY_PATH_FEEDBACK
def _gravity_path_length_change(nu):
    from scipy.integrate import quad
    def integrand(s):
        h=1-s*s;root=math.hypot(h,nu*s);u=(h+root)/nu
        return 2*nu*(1/u-1)*float(_return_path_length_weight(nu,s))
    return 2*quad(integrand,0,1,epsabs=1e-11,epsrel=2e-12,limit=150)[0]


def gravity_path_feedback_checks():
    from scipy.integrate import quad
    from mmp_pulse import _weak_magnetic_return_path
    closure_path=project_path('polyakov_closure_results.json');feedback_path=project_path('polyakov_feedback_results.json')
    closure=json.loads(closure_path.read_text(encoding='utf-8'))
    feedback=json.loads(feedback_path.read_text(encoding='utf-8'))
    distance=closure['parameters']['separation_over_re'];epsilon=1/distance
    records=[]
    for saved in closure['quadrature_refinement']:
        state=saved['roots'][-1];p=state['profile'];nu=np.asarray(p['nu']);weights=np.asarray(p['weights'])
        periods=np.asarray(p['periods']);omega2=np.asarray(p['channel_omega_squared'])
        length=periods*state['clock_scale']/distance
        bending=np.array([_gravity_path_length_change(float(v)) for v in nu])
        fractional=epsilon*bending/length;change=-2*omega2*fractional
        assert np.all(bending<0) and np.all(change>0)
        records.append(dict(flux_order=saved['order'],nu=nu.tolist(),weights=weights.tolist(),
            path_change_per_mass=bending.tolist(),fractional_period_change=fractional.tolist(),
            channel_omega_squared_change=change.tolist(),mean_forcing=float(weights@change),
            positive_mouth_dipole_forcing=float(3*weights@(change*(1-nu)))))
    fine=records[-1];state=closure['quadrature_refinement'][-1]['roots'][-1]
    p=state['profile'];nu=np.asarray(p['nu']);weights=np.asarray(p['weights']);T=np.asarray(p['periods'])
    om=np.asarray(p['channel_omega_squared']);clock=state['clock_scale']
    refinement={key:abs(records[0][key]/fine[key]-1) for key in ('mean_forcing','positive_mouth_dipole_forcing')}
    assert max(refinement.values())<1e-9
    independent=[]
    for v in (.01,.1,.5,1.,1.5,1.9):
        value=_gravity_path_length_change(v)
        direct=2*quad(lambda s:float(np.imag(_weak_magnetic_return_path(v,s,1j*1e-25)[2])/1e-25),
            0,1,epsabs=1e-10,epsrel=2e-12,limit=150)[0]
        error=abs(direct/value-1)
        assert error<1e-10
        independent.append(dict(nu=v,shape_integral=value,coordinate_variation=direct,relative_difference=error))
    far_length=2+math.asinh(math.sqrt(3))/math.sqrt(3);endpoints=[]
    for v in (1e-3,1e-4,1e-5):
        low=_gravity_path_length_change(v);high=_gravity_path_length_change(2-v)
        endpoints.append(dict(epsilon=v,low_nu_times_response=v*low,high_response_over_gap=high/v,
            response_over_leading_far_length=v*low/far_length))
    assert abs(endpoints[-1]['low_nu_times_response']/(-2*far_length)-1)<1e-5
    assert abs(endpoints[-1]['high_response_over_gap']+1)<1e-5

    mean=feedback['mean_closure'];denominator=mean['implicit_denominator']
    slope=np.asarray(mean['derivative_refinement'][-1]['channel_derivative'])
    assert np.array_equal(nu,np.asarray(feedback['refinement'][-1]['nu']))
    magnetic=np.asarray(feedback['refinement'][-1]['channel_omega_squared_change'])
    gravitational=np.asarray(fine['channel_omega_squared_change'])
    forcing=magnetic+gravitational;mean_forcing=float(weights@forcing);shift=mean_forcing/denominator
    corrected=forcing+slope*shift
    assert abs(float(weights@corrected)/shift-1)<1e-12
    # Independent reference clock and throat optical time remain boundary data.
    # Their mean-state derivatives are already in denominator, so do not count twice.
    clock_channel=2*om*(T-state['throat_optical_time'])/T
    throat_channel=-2*om/T
    clock_forcing=float(weights@clock_channel);throat_forcing=float(weights@throat_channel)
    perturbation=1e-5;ext=T-state['throat_optical_time']
    clock_fd=float(weights@(((2*math.pi/(state['throat_optical_time']+ext*math.exp(-perturbation)))**2-
        (2*math.pi/(state['throat_optical_time']+ext*math.exp(perturbation)))**2)/(2*perturbation)))
    throat_fd=float(weights@(((2*math.pi/(T+perturbation))**2-(2*math.pi/(T-perturbation))**2)/(2*perturbation)))
    assert abs(clock_fd/clock_forcing-1)<1e-8 and abs(throat_fd/throat_forcing-1)<1e-8
    return dict(model='Weak-gravity Maxwell path bending added to magnetic-current path feedback',
        inputs=dict(closure_sha256=hashlib.sha256(closure_path.read_bytes()).hexdigest(),
            magnetic_feedback_sha256=hashlib.sha256(feedback_path.read_bytes()).hexdigest()),
        parameters=dict(**closure['parameters'],gravitational_mass_over_re=1.,epsilon_m_over_d=epsilon),
        refinement=records,independent_path_variation=independent,endpoint_checks=endpoints,
        finite_cut_checks=_gravity_cut_checks(),
        causality_scope=dict(throat_optical_time_over_separation=state['throat_optical_time']*clock/distance,
            assumed_time='Static lapse N>0 and a single globally shared time t, with no temporal mouth identification offset.',
            time_function='g_inverse(dt,dt)=−1/N²<0. Future causal curves have increasing t, so this assumed static construction cannot yield backward travel.',
            implication='Small static source/path corrections alone do not introduce the relative mouth-time offset required by a chronology-changing proposal. Existence and stability of a physical matched metric remain unproved.'),
        conditional_mean=dict(magnetic_forcing=mean['fixed_collar_forcing'],gravitational_forcing=fine['mean_forcing'],
            combined_forcing=mean_forcing,implicit_denominator=denominator,
            magnetic_mean_shift=mean['mean_omega_squared_shift'],gravitational_mean_shift=fine['mean_forcing']/denominator,
            combined_mean_shift=shift,channel_omega_squared_total_change=corrected.tolist(),
            positive_mouth_dipole_total_change=float(3*weights@(corrected*(1-nu))),
            source_driven_clock_fractional_shift=mean['derivative_refinement'][-1]['clock_log_derivative']*shift),
        matching_sensitivity=dict(clock_forcing_per_log_offset=clock_forcing,
            throat_forcing_per_optical_time_offset=throat_forcing,
            mean_shift_per_independent_log_clock_offset=clock_forcing/denominator,
            mean_shift_per_independent_throat_time_offset=throat_forcing/denominator,
            mean_shift_for_test_clock_offset_epsilon=epsilon*clock_forcing/denominator,
            convention='eta is an independent log shift of the clock converting exterior lengths into throat time; theta is an independent added throat optical time, both at fixed mean state. delta w=combined_selected_shift+C_eta*eta+C_theta*theta. They are not selected physical data. The mean-state clock response is already included in 1-Fprime.'),
        equations=dict(metric='ds²=−(1−2U)dt²+(1+2U)dx², U=m/rplus+m/rminus; curl[(1−2U)Bflux]=0 to first order.',
            flux='delta nu=epsilon*Sg, Sg=2nu(1/u−1), epsilon=m/d. Fixed charges, no extra point multipoles or imposed field at infinity.',
            path='delta f/epsilon=integral Sg*K_L ds; point-mouth endpoint work vanishes. b(nu)~−2*C_L/nu at0; b(nu)~−(2−nu) at2.',
            accounting='67 already includes the first gravitational optical/Shapiro delay on the reference curves. This adds only the first field-path bending. Varying that O(m/d) optical factor along the O(m/d) displaced path is beyond this first variation, with near-interface matching still separate.',
            closure='delta w=[Fmag+Fgrav+Cclock*eta+Cthroat*theta]/[1-Fprime]; gravity/current cross terms are not included.'),
        verification=dict(all_declared_checks_pass=True,relative_flux_refinement=refinement,
            maximum_independent_path_relative_difference=max(row['relative_difference'] for row in independent),
            clock_sensitivity_finite_difference_relative_error=abs(clock_fd/clock_forcing-1),
            throat_sensitivity_finite_difference_relative_error=abs(throat_fd/throat_forcing-1)),
        scope='Selected static weak-metric and point-mouth harmonic boundary choice, not a solved two-body Einstein–Maxwell–quantum exterior. Full nonspherical throat response, physical finite-mouth interfaces, independent clock matching, orbital motion, source preparation and finite-q endpoint state remain unresolved. Combined figures are conditional first-order components, not complete physical error bars or dynamical stability. No backward time travel is demonstrated.',
        sources=['https://arxiv.org/html/1903.00223: curved-spacetime Maxwell constitutive equations; not a material recipe.',
            'https://arxiv.org/html/1807.04726v3: magnetic channels and need to support or orbit the two mouths. Bending/state update and sensitivity calculation are ours.'])


def save_gravity_feedback_benchmark():
    result=gravity_path_feedback_checks()
    project_path('polyakov_gravity_feedback_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
def _gravity_cut_checks():
    """Moving-endpoint check at fixed isotropic spheres in the point-field model."""
    from scipy.integrate import quad
    from scipy.optimize import brentq
    from mmp_pulse import _weak_magnetic_return_path
    rows=[];step=1e-25
    for nu in (.5,1.,1.5):
        point=_gravity_path_length_change(nu);limit=-nu*(2-nu);limit_errors=[]
        def integrand(s):
            h=1-s*s;u=(h+math.hypot(h,nu*s))/nu
            return 2*nu*(1/u-1)*float(_return_path_length_weight(nu,s))
        for a in (.02,.01,.005):
            # t=(1-s_cut)/a resolves a small sphere without subtracting distances.
            t=brentq(lambda t:t*(2-a*t)*(1+a*(2-t))-2*nu*(1+a*(1-t)),
                0,2,xtol=1e-14)
            sc=1-a*t;_,_,speed,u,up=_weak_magnetic_return_path(nu,sc)
            C=u*u*(1+3*sc*sc)-sc*sc+sc**4
            ue=-2*(u-1)*(u*u-sc*sc)/(u*u+sc*sc)
            endpoint=ue*(u-sc)**2/(2*math.sqrt(C)*(1-up))
            truncated=2*quad(integrand,0,sc,epsabs=1e-12,epsrel=2e-12)[0]
            cut=truncated+2*endpoint
            # Integrate the small omitted cap directly, avoiding cut-point cancellation.
            caps=2*quad(integrand,sc,1,epsabs=1e-18,epsrel=2e-12)[0]
            correction=2*endpoint-caps
            ft=(2-2*a*t)*(1+a*(2-t))-a*t*(2-a*t)+2*nu*a
            sc_epsilon=-2*a*a*t*(2-a*t)*(2-t)/ft
            direct_integral=2*quad(lambda s:float(np.imag(
                _weak_magnetic_return_path(nu,s,1j*step)[2])/step),
                0,sc,epsabs=1e-12,epsrel=2e-12)[0]
            direct=direct_integral+2*speed*sc_epsilon
            error=abs(cut-direct);decomposition=abs((cut-point)-correction)
            endpoint_error=abs(sc_epsilon-ue/(1-up))
            scaled=correction/a**4;limit_errors.append(abs(scaled/limit-1))
            assert error<2e-11 and decomposition<2e-11 and endpoint_error<2e-14
            assert endpoint<0 and correction<0 and cut<0
            rows.append(dict(nu=nu,isotropic_cutoff_over_separation=a,s_cut=sc,
                point_path_change_per_mass=point,truncated_curvature_integral=truncated,
                two_endpoint_work=2*endpoint,finite_cut_path_change_per_mass=cut,
                direct_curve_variation=direct,direct_integral=direct_integral,
                direct_moving_limit_term=2*speed*sc_epsilon,
                finite_cut_minus_point=correction,correction_over_cutoff_fourth_power=scaled,
                point_limit_coefficient=limit,point_limit_prediction=limit*a**4,
                point_limit_relative_difference=limit_errors[-1],
                direct_variation_difference=error,point_decomposition_difference=decomposition,
                endpoint_derivative_difference=endpoint_error))
        assert limit_errors[2]<limit_errors[1]<limit_errors[0] and limit_errors[-1]<.02
    return dict(records=rows,
        equations=dict(units='All coordinates use separation d. epsilon=m/d; b=delta f/epsilon=delta physical_length/m. a=R_isotropic/d is fixed during the epsilon derivative.',
            cutoff='u(s_cut)-s_cut=2a; t=(1-s_cut)/a obeys t(2-at)[1+a(2-t)]=2nu[1+a(1-t)].',
            endpoint='u_epsilon=-2(u-1)(u^2-s^2)/(u^2+s^2); E=u_epsilon*(u-s)^2/[2sqrt(C)*(1-u_s)]. b_cut=integral_-sc^sc Sg*K_L ds+2E.',
            direct='b_cut=2[integral_0^sc partial_epsilon(speed) ds+speed(sc)*partial_epsilon(sc)]. The last term is the moving integration limit; 2E instead includes the entire endpoint work in the curvature formula.',
            limit='b_cut-b_point=-nu*(2-nu)*a^4+O_nu(a^5). Each removed curvature cap contributes -(3/2)nu*(2-nu)*a^4; each endpoint work is -2nu*(2-nu)*a^4.'),
        verification=dict(all_declared_checks_pass=True,
            maximum_direct_variation_difference=max(row['direct_variation_difference'] for row in rows),
            maximum_point_decomposition_difference=max(row['point_decomposition_difference'] for row in rows),
            maximum_endpoint_derivative_difference=max(row['endpoint_derivative_difference'] for row in rows),
            fixed_nu_cutoff_limit_converges=True),
        scope='Only fixed global isotropic spheres in the selected weak-gravity point-source field model. The small-cut expansion is checked at the listed fixed flux labels, not a uniform physical-interface error bound. No actual throat continuation, mouth displacement, clock matching or finite-q state is supplied; the existing optical/Shapiro contribution is not added again.')


# END_GRAVITY_PATH_FEEDBACK


def distributed_clock_spectrum_checks():
    """Exact spectral signs for the zero-offset, fixed-delay linear JT/CFT model."""
    import mpmath as mp
    import sympy as s
    z,x,o=s.symbols('z x o');aa=2-o-z
    bb=2*z-1-aa*(z-2);cc=-z-aa*(1-2*z);dd=-aa*z
    remainder=x*(x-o)/((x-1)**2*(z+x))-1/x-aa/x**2
    assert s.factor(remainder-(bb*x*x+cc*x+dd)/(x*x*(x-1)**2*(x+z)))==0
    q=s.symbols('q');resonant_R=(1+q)*(1-q+q*q)/((1-q)*(1+q*q))
    resonance_derivative=s.simplify((-s.pi*q*s.diff(resonant_R,q)).subs(q,-1))
    assert resonance_derivative==3*s.pi/4
    # At L=4pi, n=2: the apparent 0/0 spectral weight is pi/32.
    resonance_weight=s.simplify(s.pi**2*(1-s.Rational(1,4))/(6*4*s.pi))
    assert resonance_weight==s.pi/32
    records=[];coefficient_min=math.inf
    with mp.workdps(50):
        cache={}
        def data(length,N):
            key=(length,N)
            if key not in cache:
                omega=2*mp.pi/length;oo=omega**2;angle=mp.pi**2/length;terms=[]
                for n in range(1,N+1):
                    frequency=n*omega;xx=frequency**2;cos2=mp.cos(n*angle)**2
                    # Entire continuation through frequency=1; mp.sinc(v)=sin(v)/v.
                    ratio=-mp.pi*mp.sinc(mp.pi*(frequency-1)/2)/(2*(frequency+1))
                    hh=xx*(xx-oo)*ratio**2
                    assert hh>=0
                    terms.append((xx,cos2,mp.sin(n*angle)**2,hh))
                cache[key]=(oo,angle,terms)
            return cache[key]
        def spectral(length,zz,N):
            oo,angle,terms=data(length,N);a1=2-oo-zz
            s2=mp.pi**2/6-mp.pi*angle/2+angle**2/2-mp.cos(angle)**2
            s4=mp.pi**4/90-mp.pi**2*angle**2/6+mp.pi*angle**3/3-angle**4/6-mp.cos(angle)**2
            anti=8/(3*length)*(s2/oo+a1*s4/oo**2+mp.fsum(
                hh/(zz+xx)-cos2/xx-a1*cos2/xx**2 for xx,cos2,_,hh in terms[1:]))
            sine2=mp.pi*angle/2-angle**2/2
            sine4=mp.pi**2*angle**2/6-mp.pi*angle**3/3+angle**4/6
            common=8*(zz+oo)/(3*length)*(sine2/oo-zz*sine4/oo**2+mp.fsum(
                sin2*zz**2/(xx**2*(zz+xx)) for xx,_,sin2,_ in terms))
            X=(N+1)**2*oo;assert X>max(1,abs(zz))
            b1=2*zz-1-a1*(zz-2);c1=-zz-a1*(1-2*zz);d1=-a1*zz
            K=(abs(b1)+abs(c1)/X+abs(d1)/X**2)/((1-1/X)**2*(1-abs(zz)/X))
            anti_tail=8*K/(15*length*oo**3*N**5)
            common_tail=8*abs(zz**2*(zz+oo))/(15*length*oo**3*(1-abs(zz)/X)*N**5)
            return anti,common,anti_tail,common_tail
        cases=[('D0',[mp.mpf(1)],[mp.pi]),
            ('distributed_long_tail',[mp.mpf('.7'),mp.mpf('.2'),mp.mpf('.1')],[mp.pi+mp.mpf('.1'),mp.mpf(10),mp.mpf(40)]),
            ('exact_resonant_channel',[mp.mpf('.75'),mp.mpf('.25')],[mp.pi,4*mp.pi])]
        D0_errors=[]
        for name,weights,lengths in cases:
            mean=mp.fsum(w*(2*mp.pi/L)**2 for w,L in zip(weights,lengths))
            b0=mp.pi*(mean-1)/6;assert b0>0
            for lam in (mp.mpc('.4','.7'),mp.mpc('.2','2.3'),mp.mpf('.6')):
                zz=lam**2;plus=[];minus=[];acc=[]
                for L in lengths:
                    oo=(2*mp.pi/L)**2;ex=mp.exp(-mp.pi*lam);ed=mp.exp(-(L-mp.pi)*lam)
                    den=1-mp.exp(-L*lam)
                    plus.append((zz+oo)*(1+ex)*(1+ed)/den)
                    minus.append((zz+oo)*(1-ex)*(1-ed)/den)
                    acc.append(spectral(L,zz,1024))
                    coefficient_min=min(coefficient_min,float(min(t[3] for t in data(L,1024)[2][1:])))
                direct_a=(b0*zz*(zz+1)+lam*mp.fsum(w*v for w,v in zip(weights,plus))/3)/(zz+1)**2
                direct_s=b0*(zz+1)+mp.fsum(w*v for w,v in zip(weights,minus))/(3*lam)
                series_a=b0+mp.fsum(w*v[0] for w,v in zip(weights,acc))
                series_s=b0*(zz+1)+mp.fsum(w*v[1] for w,v in zip(weights,acc))
                bound_a=mp.fsum(w*v[2] for w,v in zip(weights,acc));bound_s=mp.fsum(w*v[3] for w,v in zip(weights,acc))
                error_a=abs(direct_a-series_a);error_s=abs(direct_s-series_s)
                assert error_a<=bound_a+mp.mpf('1e-39') and error_s<=bound_s+mp.mpf('1e-39')
                if mp.im(zz):
                    assert mp.im(direct_a)*mp.im(zz)<0 and mp.im(direct_s)*mp.im(zz)>0
                else: assert direct_a>0 and direct_s>0
                if name=='D0':
                    H=lam*(zz+1)*mp.sinh(mp.pi*lam/2)+4/(3*mp.pi)*(zz+4)*mp.cosh(mp.pi*lam/2)
                    D0_errors.append(float(abs(direct_a-b0*lam*H/(mp.sinh(mp.pi*lam/2)*(zz+1)**2))))
                records.append(dict(case=name,weights=list(map(float,weights)),loop_periods=list(map(float,lengths)),
                    mean_omega_squared=float(mean),b0=float(b0),lambda_value=[float(mp.re(lam)),float(mp.im(lam))],
                    antisymmetric_value=[float(mp.re(direct_a)),float(mp.im(direct_a))],
                    symmetric_value=[float(mp.re(direct_s)),float(mp.im(direct_s))],
                    series_terms=1024,antisymmetric_error=float(error_a),antisymmetric_tail_bound=float(bound_a),
                    symmetric_error=float(error_s),symmetric_tail_bound=float(bound_s)))
        assert max(D0_errors)<1e-45
    return dict(model='Linear zero-offset NS-vacuum clock feedback with positive weighted, fixed channel delays; no finite-offset dynamics inferred',
        conventions=dict(delays='D_nu>=0, L_nu=pi+D_nu, Omega_nu=2pi/L_nu; weights positive and normalized.',
            opening='b0=pi*(mean(Omega_nu^2)-1)/6>0; b0 and bulk A,B use the existing source/alpha convention.',
            coordinates='u=tau+chi,v=tau-chi; clocks r(s)=s+xi_R(s),l(s)=s+xi_L(s). C=A-B, S=A+B; P=d^3/dx^3+d/dx. xi_antisymmetric=(xi_R-xi_L)/2 and j_antisymmetric=(delta_J_R-delta_J_L)/2.',
            clock='xi_antisymmetric_prime(s)=[C(s+pi/2)+C(s-pi/2)]/(2*b0); differentiated clock equations alone must retain the original bulk/mouth integration constraints.'),
        equations=dict(characteristics='G_diff=P(lambda)+(1/(3b0))*mean[(lambda^2+Omega^2)*(1+exp(-pi*lambda))*(1+exp(-D*lambda))/(1-exp(-L*lambda))]; G_sum changes both plus signs to minus.',
            channel_recurrence='(1-E_L)j_antisymmetric=(1/3)*partial*(partial^2+Omega^2)*(1+E_D)*xi_antisymmetric. b0*partial^2*(partial^2+1)*xi_antisymmetric=-(1+E_pi)*mean(j_antisymmetric).',
            transfer='With incoming per-pass source i_nu, D_A*xi=F_boundary-(1+exp(-pi*lambda))*mean[i_nu/(1-exp(-L_nu*lambda))], D_A=b0*lambda*G_diff. Initial-state and one-pass sources have different numerators.',
            spectral='z=lambda^2; F_A=b0*lambda*G_diff/(z+1)^2=b0+mean[sum_(n>=2)c_n/(z+n^2*Omega^2)]. c_n=8*n^2*Omega^4*(n^2-1)*cos^2(n*pi^2/L)/[3L*(n^2*Omega^2-1)^2]>=0.',
            symmetric='F_S=b0*G_sum/lambda=b0*(z+1)+(8/3)*mean[L^-1*sum_(n>=1)(z+Omega^2)*sin^2(n*pi^2/L)/(z+n^2*Omega^2)].'),
        proof=dict(kernel='The periodic resolvent identity gives lambda*R_plus=4/L+(8/L)*sum_(n>=1) z*cos^2(n*pi^2/L)/(z+n^2*Omega^2). R_minus has no constant and replaces cos^2 by sin^2.',
            subtraction='Multiplication by z+Omega^2 removes the n=1 vacuum pole. Dividing by (z+1)^2, subtracting the value and derivative at -1 and using b_L=pi*(Omega^2-1)/6 yields the displayed spectral identity. Its constant at positive infinity is b_L. Resonances use the continuous positive residue below, rather than substituting 0/0.',
            antisymmetric_sign='Im(F_A)=-Im(z)*mean[sum c_n/abs(z+n^2*Omega^2)^2]; F_A>0 for real z>=0 since b0>0. Therefore G_diff has no zero with Re(lambda)>0.',
            symmetric_sign='Im(F_S)=Im(z)*[b0+(8/3)*mean[L^-1*sum (n^2-1)*Omega^2*sin^2(n*pi^2/L)/abs(z+n^2*Omega^2)^2]]. This is strictly of the sign of Im(z), and F_S>0 for real z>=0. G_sum likewise has no right-half-plane zero.',
            continuum='At any positive z0, each positive spectral sum equals its bounded closed form minus bounded b_L, uniformly for L>=pi. Tonelli therefore applies to any normalized positive delay distribution. The bounded ratio (z0+a)/(z+a), a>=0, on compact subsets off the negative real axis gives local uniform absolute convergence. Unbounded delays may produce a spectral cut; no isolated neutral-mode counting is asserted there.',
            truncation='Numerics subtract cos^2/x+A*cos^2/x^2 and restore their exact Fourier sums, with A=2-Omega^2-z. B=2z-1-A(z-2), C=-z-A(1-2z), D=-Az. With x=n^2*Omega^2 and X=(N+1)^2*Omega^2, the remaining rational numerator is B*x^2+C*x+D over x^2*(x-1)^2*(x+z). Its absolute tail is at most8K/(15L*Omega^6*N^5), K=(abs(B)+abs(C)/X+abs(D)/X^2)/[(1-1/X)^2*(1-abs(z)/X)], requiring X>max(1,abs(z)). This follows from sum_(n>N)n^-6<=1/(5N^5). The sine tail is at most8*abs(z^2*(z+Omega^2))/[15L*Omega^6*(1-abs(z)/X)*N^5]. These are analytic truncation bounds, not interval bounds on floating-point arithmetic.'),
        resonance=dict(example_loop_period='4pi',resolvent_derivative_at_lambda_i=str(resonance_derivative),
            spectral_weight_at_n2=str(resonance_weight),general_weight='At n*Omega=1: c_n=pi^2*(1-Omega^2)/(6L). n=1 contributes zero.',
            gauge='For a finite nonresonant ensemble, the eliminated C characteristic contains (lambda^2+1)^2, generically giving double +/-i roots. C=exp(ix) reconstructs zero; C=x*exp(ix) reconstructs i*exp(i*tau)*(sin(chi)+chi*sec(chi)), a bounded AdS-frame perturbation. At L=2pi*n,n>=2, F_A has a positive pole at z=-1 and one factor cancels. Transport poles, accidental degeneracies and homogeneous state modes require the uneliminated equations for complete neutral-mode counting; none changes the right-half-plane exclusion.'),
        cases=records,verification=dict(all_declared_checks_pass=True,precision_decimal_digits=50,
            rational_tail_identity='0',exact_resonance_derivative='3*pi/4',exact_resonance_weight='pi/32',
            minimum_sampled_positive_spectral_numerator=coefficient_min,maximum_D0_identity_error=max(D0_errors),
            maximum_antisymmetric_identity_error=max(r['antisymmetric_error'] for r in records),
            maximum_symmetric_identity_error=max(r['symmetric_error'] for r in records),
            all_identity_errors_within_analytic_truncation_bounds=True),
        scope='Analytic absence of exponentially growing characteristic modes in this specified linear, zero-offset, free-NS-vacuum, common-mouth, fixed-delay leading-JT closure. This is not damping, an all-time response bound, nonlinear stability, robustness to a shifted/excited state, finite-collar dynamics, a solved nonspherical4D exterior or chronology protection. Multiplying transport denominators can add or hide neutral modes; no frequency sweep substitutes for the proof.',
        sources=['https://arxiv.org/html/1807.04726v3: NS Casimir normalization and leading throat framework.',
            'https://arxiv.org/pdf/math-ph/0412028: chiral stress transformation and positive-energy CFT input. The distributed-delay spectral identity and its application here are our derivation.'])


def stationary_offset_gate_checks():
    """Calibrate the existing stationary timing bound with the saved full profile."""
    from functools import lru_cache
    from scipy.optimize import brentq
    source = project_path('polyakov_closure_results.json')
    baseline = json.loads(source.read_text(encoding='utf-8'))
    alpha = baseline['parameters']['alpha']
    distance = baseline['parameters']['separation_over_re']

    @lru_cache(None)
    def collar(w):
        clock, cuts = _magnetic_closure_collar(w, alpha)
        return clock, 2*cuts['100.0']['optical_x']

    rows = []
    for profile in baseline['quadrature_refinement']:
        old = profile['roots'][-1]
        weights = np.array(old['profile']['weights'])
        external = old['clock_scale']*(np.array(old['profile']['periods'])-old['throat_optical_time'])
        # Fixed weak exterior paths; the mean collar/clock responds to the offset.
        emin = old['clock_scale']*old['shortest_exterior_time']
        assert weights.min()>0 and abs(weights.sum()-1)<1e-14 and external.min()>=emin
        for ratio in (0., 1., 3., 10.):
            shift = ratio*distance
            def residual(w):
                clock, throat = collar(w)
                return w-float(weights@((2*math.pi/(throat+(external+shift)/clock))**2))
            w = brentq(residual, 1.0001, 2., xtol=2e-13, rtol=2e-14)
            clock, throat = collar(w)
            long = throat+(external+shift)/clock
            short = throat+(external-shift)/clock
            pairs = (1-(short/long)**2)/24
            # Positive neutral L0=1 pair mixtures exactly balance each channel.
            balanced_short = 2*math.pi*(pairs-1/24)/short**2
            vacuum_long = -math.pi/(12*long**2)
            balance_error = float(np.max(abs(balanced_short-vacuum_long)))
            gate = w*throat**2/math.pi**2-1
            exit_delay = clock*throat-shift
            assert short.min()>0 and pairs.min()>=0 and pairs.max()<=1/24
            assert balance_error<1e-14 and abs(residual(w))<1e-11
            assert gate>=0 and exit_delay>emin
            if ratio==0: assert abs(w-old['omega_squared'])<2e-12
            rows.append(dict(order=profile['order'],offset_over_d=ratio,omega_squared=w,
                clock_over_d=clock/distance,throat_time_over_d=clock*throat/distance,
                earlier_direction_exit_delay_over_d=exit_delay/distance,
                fastest_complete_return_over_d=(exit_delay+emin)/distance,
                minimum_short_period=float(short.min()),maximum_pair_probability=float(pairs.max()),
                finite_collar_timing_gate=gate,root_residual=residual(w),balance_error=balance_error))
    differences = [abs(rows[i][key]-rows[i+4][key]) for i in range(4)
                   for key in ('omega_squared','earlier_direction_exit_delay_over_d')]
    assert max(differences)<1e-9
    result = dict(model='Stationary balanced mean-source offset branch; fixed weak exterior paths',
        input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        parameters=dict(alpha=alpha,separation_over_re=distance,cutoff=100.),cases=rows,
        verification=dict(all_declared_checks_pass=True,maximum_64_128_difference=max(differences)),
        equations=dict(periods='P_plus/minus=Tth+(E_nu+/-S)/A; require P_minus>0',
            closure='w=<[2pi/P_plus]^2>; p_nu=[1-(P_minus/P_plus)^2]/24',
            timing='Positive NS excitation energies and total flux balance imply w<=[2pi/(Tth+(Emin+S)/A)]^2. If sqrt(w)*Tth>=pi, A*Tth>=S+Emin.'),
        scope='Calibration of the existing stationary bound, not a new no-go theorem. The four offsets are separate assumed states, not a preparation trajectory. No driven 4D matching, full angular response, finite-q control, material production or past travel is established. Times use d/c; no physical d or calendar conversion is selected.')
    project_path('stationary_offset_gate_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__ == '__main__':
    result = polyakov_response_checks()
    project_path('polyakov_response_results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('plot','conventions','initial_data')}, indent=2))
