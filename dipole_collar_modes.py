"""Exact angular zero modes inside the partner focus; radial coupling remains.

Flat point-pair model only. Reads frozen increment88 results, never reruns them.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import eigh_tridiagonal

from incoming_modes import collar_mode


def correction(q, b, d, t):
    """Stable real F and its radial derivative g, using mpmath arithmetic."""
    q, b, d, t = map(mp.mpf, (q, b, d, t))
    if not (q > 0 and q == int(q) and 0 < b < d and 0 <= t <= 1):
        raise ValueError('Require integer q>0, 0<b<d, and 0<=t<=1')
    e, c = b/d, 1-2*t
    D = mp.sqrt(1+2*e*c+e*e)
    F = q/2*mp.log1p(4*e*e*t/((1+e+D)*(D+1-e)*(1+e)))
    g = 2*q*e*t/(d*D*(1+e)*(1+e+D))
    return F, g


def beta_nodes(a, b, order):
    """Small Gaussian quadrature for normalized Beta(a,b); no huge Gamma weights."""
    reflected = a > b
    if reflected:
        a, b = b, a
    S = a+b
    n = np.arange(1, order, dtype=float)
    diagonal = np.r_[a/S, ((S-2)*a+2*n*(n+S-1))/((2*n+S-2)*(2*n+S))]
    off = np.sqrt(n*(n+a-1)*(n+b-1)*(n+S-2)/((2*n+S-3)*(2*n+S-1)))/(2*n+S-2)
    nodes, vectors = eigh_tridiagonal(diagonal, off)
    # Preserve tiny distance from the south pole before converting to high precision.
    ts = [1-mp.mpf(float(x)) if reflected else mp.mpf(float(x)) for x in nodes]
    weights = [mp.mpf(float(x*x)) for x in vectors[0]]
    total = mp.fsum(weights)
    return ts, [w/total for w in weights]


def angular_moments(q, k, b, d, order=40):
    """Centered normalized moments, including basis change too small for subtraction."""
    if not (isinstance(q, int) and isinstance(k, int) and 0 <= k < q and order >= 4):
        raise ValueError('Require integer 0<=k<q and quadrature order>=4')
    with mp.workdps(50):
        b, d = mp.mpf(b), mp.mpf(d)
        ts, weights = beta_nodes(k+1, q-k, order)
        mean_t = mp.mpf(k+1)/(q+1)
        F0, g0 = correction(q, b, d, mean_t)
        pairs = [correction(q, b, d, t) for t in ts]
        offsets = [F-F0 for F, g in pairs]
        logZ = mp.log1p(mp.fsum(w*mp.expm1(2*x) for w, x in zip(weights, offsets)))
        corrected = [w*mp.exp(2*x-logZ) for w, x in zip(weights, offsets)]
        dg = [g-g0 for F, g in pairs]
        mean_dg = mp.fsum(w*x for w, x in zip(corrected, dg))
        centered = [x-mean_dg for x in dg]
        variance = mp.fsum(w*x*x for w, x in zip(corrected, centered))
        # Direct difference of unit vectors, not 2-2*overlap at tiny separation.
        basis_distance = mp.sqrt(mp.fsum(w*mp.expm1(x-logZ/2)**2 for w, x in zip(weights, offsets)))
        beta_variance = mp.mpf((k+1)*(q-k))/((q+1)**2*(q+2))
        sampled_variance = mp.fsum(w*(t-mean_t)**2 for w, t in zip(weights, ts))
        leading = q*b/(d*d)*mp.sqrt(beta_variance)
        e = b/d
        logM = -q*mp.log1p(-e*e)
        lower_factor = mp.exp(-3*mp.log1p(e)-logM/2)
        upper_factor = mp.exp(-3*mp.log1p(-e)+logM/2)
        profiles = []
        for t in (max(mp.mpf('0'), mean_t-2*mp.sqrt(beta_variance)), mean_t,
                  min(mp.mpf('1'), mean_t+2*mp.sqrt(beta_variance))):
            F, g = correction(q, b, d, t)
            profiles.append(dict(t=float(t), relative_amplitude_minus_one=float(mp.expm1(F-F0-logZ/2)),
                                 radial_log_derivative=float(g-g0-mean_dg)))
        return dict(q=q, k=k, b=float(b), d=float(d), order=order,
                    log_normalization=float(2*F0+logZ), mean_g=float(g0+mean_dg),
                    radial_residual_norm=float(mp.sqrt(variance)), leading_radial_residual=float(leading),
                    relative_to_leading_minus_one=float(mp.sqrt(variance)/leading-1),
                    normalized_basis_distance=float(basis_distance),
                    analytic_radial_residual_lower=float(leading*lower_factor),
                    analytic_radial_residual_upper=float(leading*upper_factor),
                    analytic_enclosure_factors=[float(lower_factor),float(upper_factor)],
                    berry_diagonal_residual=float(abs(mp.fsum(w*x for w, x in zip(corrected, centered)))),
                    beta_variance_relative_error=float(abs(sampled_variance/beta_variance-1)),
                    profiles=profiles)


def corrected_mode(q, k, rho, z, energy, d, order=40):
    """Reduced f about the upper focus; z is relative to that focus."""
    b = math.hypot(rho, z)
    t = (1-z/b)/2
    data = angular_moments(q, k, b, d, order)
    with mp.workdps(40):
        factor = float(mp.exp(correction(q, b, d, t)[0]-mp.mpf(data['log_normalization'])/2))
    return collar_mode(q, k, rho, z, energy)*factor


def moderate_checks():
    """Differentiate normalization integral independently, then the full spinor."""
    rows = []
    with mp.workdps(38):
        sx, sy, sz = mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]), mp.diag([1, -1])
        q, d, E = 5, mp.mpf(1), mp.mpf('1.3')
        for k in (0, 2, 4):
            beta = mp.beta(k+1, q-k)
            density = lambda t: t**k*(1-t)**(q-1-k)/beta
            def norm(radius):
                return mp.quad(lambda t: density(t)*mp.exp(2*correction(q, radius, d, t)[0]), [0, .25, .5, .75, 1])
            for b in map(mp.mpf, ('.05', '.3', '.7')):
                Z = norm(b)
                # Numerical derivative of the integral itself; does not use supplied g.
                logZprime = mp.diff(norm, b)/Z
                mean_g = mp.quad(lambda t: density(t)*mp.exp(2*correction(q,b,d,t)[0])*correction(q,b,d,t)[1], [0,.25,.5,.75,1])/Z
                data = angular_moments(q, k, b, d)
                # First jet of the independently differentiated normalizer is exact
                # for this first-order operator at radius b, without nested quadrature.
                def f(rho, z):
                    r = mp.sqrt(rho*rho+z*z)
                    theta = mp.atan2(rho, z)
                    t = mp.sin(theta/2)**2
                    angular = mp.sqrt(rho/(2*r*r*beta))*mp.sin(theta/2)**k*mp.cos(theta/2)**(q-1-k)
                    angular *= mp.exp(correction(q,r,d,t)[0]-(mp.log(Z)+logZprime*(r-b))/2+1j*E*r)
                    return angular*mp.matrix([mp.cos(theta/2),mp.sin(theta/2)])
                residuals, angular_errors, derivative_errors, tau_errors, implementation_errors = [], [], [], [], []
                for t in map(mp.mpf, ('.17', '.53', '.83')):
                    theta = 2*mp.asin(mp.sqrt(t))
                    rho, z = b*mp.sin(theta), b*mp.cos(theta)
                    v = f(rho,z)
                    dr = mp.matrix([mp.diff(lambda x:f(x,z)[i],rho) for i in range(2)])
                    dz = mp.matrix([mp.diff(lambda x:f(rho,x)[i],z) for i in range(2)])
                    A = q/2*((z+d)/mp.sqrt(rho*rho+(z+d)**2)-z/b)
                    Hf = -1j*sx*dr-1j*sz*dz+sy*(mp.mpf(k)+mp.mpf('.5')-A)/rho*v
                    F, g = correction(q,b,d,t)
                    target = -1j*(g-mean_g)*v
                    residuals.append(float(mp.norm(Hf-E*v-target)/mp.norm(v)))
                    angular_errors.append(float(abs(mp.diff(lambda th: correction(q,b,d,mp.sin(th/2)**2)[0],theta)+(A-q*t)/mp.sin(theta))))
                    derivative_errors.append(float(abs(mp.diff(lambda r: correction(q,r,d,t)[0],b)-g)))
                    D = mp.sqrt(1+2*(b/d)*(1-2*t)+(b/d)**2)
                    tau_errors.append(float(abs(mp.diff(lambda tt: correction(q,b,d,tt)[1],t)-q*b/(d*d*D**3))))
                    reference = np.array([complex(component) for component in v])
                    numerical = corrected_mode(q,k,float(rho),float(z),float(E),float(d))
                    implementation_errors.append(float(np.linalg.norm(numerical-reference)/np.linalg.norm(reference)))
                independent_flux = Z/mp.exp(data['log_normalization'])
                independent_variance = mp.quad(lambda t: density(t)*mp.exp(2*correction(q,b,d,t)[0])*(correction(q,b,d,t)[1]-mean_g)**2, [0,.25,.5,.75,1])/Z
                rows.append(dict(k=k,b=float(b),max_full_dirac_residual=max(residuals),
                                 max_angular_ode_residual=max(angular_errors),max_F_derivative_error=max(derivative_errors),
                                 max_g_tau_derivative_error=max(tau_errors),
                                 normalization_derivative_error=float(abs(logZprime/2-mean_g)),
                                 unit_flux_error=float(abs(independent_flux-1)),
                                 implementation_error=max(implementation_errors),
                                 residual_norm_quadrature_relative_error=float(abs(data['radial_residual_norm']/mp.sqrt(independent_variance)-1))))
        # F is analytic in t at both axes for b<d, so original patchwise powers survive.
        b = mp.mpf('.7')
        north, south = correction(q,b,d,0)[0], correction(q,b,d,1)[0]
        axis_error = max(abs(north), abs(south+q/2*mp.log(1-b*b)))
        gauge_error = abs(((mp.mpf('.5')-q)-(mp.mpf('.3')-q))-(mp.mpf('.5')-mp.mpf('.3')))
    return rows, float(axis_error), float(gauge_error)


def checks():
    frozen_path = Path(__file__).with_name('incoming_modes_results.json')
    frozen = json.loads(frozen_path.read_text())
    saved = frozen['saved_collar']
    if not (frozen['passed'] and saved['passed']):
        raise ValueError('Frozen incoming-mode controls failed')
    baseline_hash = hashlib.sha256(Path(__file__).with_name('incoming_modes.py').read_bytes()).hexdigest()
    if baseline_hash != frozen['source_sha256']:
        raise ValueError('Monopole helper source differs from the frozen input')
    p = saved['parameters']
    q,b,d = p['q'],p['isotropic_b'],p['separation_d']
    actual, refinements = [], []
    for k in (0,16,q//2,q-1):
        coarse = angular_moments(q,k,b,d,24)
        fine = angular_moments(q,k,b,d,40)
        refinements.append(abs(coarse['radial_residual_norm']/fine['radial_residual_norm']-1))
        old = next(row for row in saved['channels'] if row['k']==k)['flat_rms_defect']
        fine['previous_potential_residual']=old
        fine['new_over_previous_residual']=fine['radial_residual_norm']/old
        fine['leading_new_over_previous']=2/math.sqrt(q+1)
        actual.append(fine)
    moderate, axes, gauge = moderate_checks()
    tests = {
        'full_reduced_dirac_residual_identity': max(r['max_full_dirac_residual'] for r in moderate)<1e-30,
        'exact_angular_equation_and_radial_F_derivative': max(max(r['max_angular_ode_residual'],r['max_F_derivative_error']) for r in moderate)<1e-30,
        'independent_normalizer_derivative': max(r['normalization_derivative_error'] for r in moderate)<1e-30,
        'radial_log_derivative_slope_identity': max(r['max_g_tau_derivative_error'] for r in moderate)<1e-30,
        'unit_flux_and_implemented_mode': max(max(r['unit_flux_error'],r['implementation_error']) for r in moderate)<1e-10,
        'independent_residual_norm_integral': max(r['residual_norm_quadrature_relative_error'] for r in moderate)<1e-10,
        'regular_axis_factors_and_gauge_covariance': max(axes,gauge)<1e-30,
        'actual_collar_quadrature_refinement': max(refinements)<1e-10,
        'actual_collar_beta_variance': max(r['beta_variance_relative_error'] for r in actual)<1e-10,
        'actual_collar_small_b_leading': max(abs(r['relative_to_leading_minus_one']) for r in actual)<2e-7,
        'actual_collar_exact_analytic_enclosure': all(r['analytic_radial_residual_lower']<=r['radial_residual_norm']<=r['analytic_radial_residual_upper'] for r in actual),
        'zero_diagonal_radial_connection': max(r['berry_diagonal_residual'] for r in actual)<1e-50,
    }
    tests={key:bool(value) for key,value in tests.items()}
    assert all(tests.values()), tests
    return dict(passed=True, tests=tests, moderate_checks=moderate,axis_error=axes,gauge_error=gauge,
                actual_collar=actual,actual_refinement_relative_errors=refinements,
                saved_collar_parameters=p, frozen_input_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest(),
                upstream_closure_input_sha256=saved['input_sha256'],
                monopole_helper_source_sha256=baseline_hash,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                units='Saved-collar lengths in re; radial residual in inverse re; normalized basis distance dimensionless.',
                equations=dict(F='(q/2)*log1p[2*eps^2*(1-cos(theta))/((1+eps+D)*(D+1-eps)*(1+eps))]',
                               g='q*eps*(1-cos(theta))/[d*D*(1+eps)*(1+eps+D)]',
                               D='sqrt(1+2*eps*cos(theta)+eps^2), eps=b/d',
                               angular='h_corrected=h_monopole*exp(F)/sqrt(<exp(2F)>_Beta)',
                               radial='partial_b Y=(g-<g>_corrected)Y; <Y,partial_bY>=0',
                               full_residual='(H-E)Psi=-i exp(iEb)partial_bY/b',
                               norm='sphere-weighted residual norm=sqrt(Var_corrected(g))',
                               slope='partial_t g=(q*b/d^2)*D^(-3), t=sin^2(theta/2)',
                               enclosure='M=(1-eps^2)^(-q); sigma_lead*(1+eps)^(-3)/sqrt(M) <= sigma <= sigma_lead*(1-eps)^(-3)*sqrt(M)',
                               leading='q*b/d^2*sqrt(Var_Beta(t)); new/old leading=2/sqrt(q+1)'),
                scope=['Exact angular basis on upper-centered spheres 0<b<d in flat point-pair field.',
                       'Radial coupling remains; diagonal zero is not absence of coupling.',
                       'Radial derivative is orthogonal to all angular zero modes: same k by normalization, different k by azimuthal phase.',
                       'Analytic enclosures are exact formulas; reported decimal endpoints are not interval-arithmetic certificates.',
                       'Fixed-order numerical quadrature is validated at the reported moderate-q tests and saved collar, not uniformly for all q and b/d.',
                       'No global return amplitude, continuum current, curved finite-mouth match or time travel result.',
                       'No q-mode continuation asserted across b=d, where the second focus changes enclosed flux.'])


if __name__=='__main__':
    result=checks()
    Path(__file__).with_name('dipole_collar_modes_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
