"""Flux operators and a small stationary Weyl absorber calibration.

The frozen planar problem is an exact local reduction of H_mu, not the
two-mouth exterior. Absorbed flux below is numerical loss, never escape.
Run .venv/Scripts/python.exe -B -W error open_boundary_spinor.py.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm


SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.diag([1., -1.])


def boundary_flux(f, normal):
    """Return density, outward and inward CHARACTERISTIC flux densities.

    f has shape (2, ...), with reduced measure d_rho d_z. These positive
    characteristic parts are not asymptotic scattering-channel probabilities.
    """
    f = np.asarray(f, dtype=complex)
    normal = np.asarray(normal, dtype=float)
    if f.ndim < 1 or f.shape[0] != 2 or not np.all(np.isfinite(f)):
        raise ValueError('Require a finite two-component spinor')
    if normal.shape != (2,) or not np.all(np.isfinite(normal)) or not np.isclose(normal @ normal, 1., atol=1e-12, rtol=0):
        raise ValueError('Require a unit (n_rho,n_z) normal')
    action = np.einsum('ab,b...->a...', normal[0]*SX+normal[1]*SZ, f)
    outgoing = np.sum(abs((f+action)/2)**2, axis=0)
    incoming = np.sum(abs((f-action)/2)**2, axis=0)
    return {'density': np.sum(abs(f)**2, axis=0), 'outgoing': outgoing,
            'incoming': incoming, 'net': outgoing-incoming}


def outgoing_ratio(energy, tangent, W):
    """b/a at a +z face for constant real W and exp(i*k*rho-i*E*t).

    Positive E only. Includes the decaying evanescent branch, but threshold
    flux normalization is singular. The full varying exterior needs its own map.
    """
    if not np.all(np.isfinite([energy, tangent, W])) or energy <= 0:
        raise ValueError('Require finite real coefficients and positive energy')
    p = np.sqrt(complex(energy*energy-tangent*tangent-W*W))
    return (tangent+1j*W)/(energy+p)


def flux_modes(energy, tangent, W):
    """Columns have +1/-1 z flux and fixed component phases at x=0."""
    outgoing_ratio(energy, tangent, W)
    p2 = energy*energy-tangent*tangent-W*W
    if p2 <= 0:
        raise ValueError('Unit-flux modes require E^2 > k^2 + W^2')
    p = np.sqrt(p2)
    scale = np.sqrt(2*p*(energy+p))
    return p, np.array([[energy+p, tangent-1j*W],
                        [tangent+1j*W, energy+p]], dtype=complex)/scale


def modal_amplitudes(f, energy, tangent, W):
    """Complex +z/-z amplitudes at the same reference plane as f."""
    _, modes = flux_modes(energy, tangent, W)
    return np.linalg.solve(modes, np.asarray(f, dtype=complex))


def cap_scattering(energy, tangent, W, width, strength, power=3, wall_phase=0., rtol=2e-11):
    """Solve a finite CAP eta=strength*(x/width)^power backed by j_z=0.

    Integrate two complex components backwards, normalize to incident flux 1,
    and independently quadrature the loss density. No dense spatial matrices.
    """
    _, modes = flux_modes(energy, tangent, W)
    if not np.all(np.isfinite([width, strength, power, wall_phase, rtol])) or width <= 0 or strength < 0 or power < 0 or not 0 < rtol < 1:
        raise ValueError('Invalid absorber parameters')
    eta = lambda x: strength*(x/width)**power
    generator = lambda x: 1j*(energy+1j*eta(x))*SZ+tangent*SY-W*SX
    wall = np.array([1., np.exp(1j*wall_phase)])
    # ponytail: bounded-layer shooting; use stabilized impedance for large eta*L.
    solution = solve_ivp(lambda x, f: generator(x) @ f, (width, 0.), wall,
                         method='DOP853', rtol=rtol, atol=rtol/100, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    incident, reflected = np.linalg.solve(modes, solution.y[:, -1])
    reflection = reflected/incident
    density = lambda x: float(np.sum(abs(solution.sol(x)/incident)**2))
    absorbed, quadrature_error = quad(lambda x: 2*eta(x)*density(x), 0., width,
                                      epsabs=2e-10, epsrel=2e-10, limit=200)
    final = wall/incident
    wall_flux = float(np.real(final.conj() @ SZ @ final))
    return {'energy': energy, 'tangent': tangent, 'W': W, 'width': width,
            'strength': strength, 'power': power, 'wall_phase': wall_phase,
            'r': [float(reflection.real), float(reflection.imag)],
            'reflection_probability': float(abs(reflection)**2),
            'reflection_phase': float(np.angle(reflection)),
            'numerically_absorbed_flux': absorbed, 'wall_flux': wall_flux,
            'flux_balance_error': abs(1-abs(reflection)**2-absorbed-wall_flux),
            'loss_quadrature_error_estimate': quadrature_error, 'ode_evaluations': solution.nfev}


def checks():
    # Reuse the existing symbolic full-variable-W operator check, not its runs.
    from spinful_return import operator_identity
    import sympy as s
    operator_residual = operator_identity()
    # Exact evanescent solution: E=2, k=1, W=2, p=i; it carries no z flux.
    x = s.symbols('x', nonnegative=True)
    beta_exact = (4+3*s.I)/5
    evanescent = s.Matrix([1, beta_exact])*s.exp(-x)
    sx = s.Matrix([[0, 1], [1, 0]])
    sy = s.Matrix([[0, -s.I], [s.I, 0]])
    sz = s.diag(1, -1)
    evanescent_residual = [s.simplify(v) for v in
        -s.I*sz*evanescent.diff(x)+(sx+2*sy)*evanescent-2*evanescent]
    evanescent_current = s.simplify((evanescent.conjugate().T*sz*evanescent)[0])
    evanescent_decay = s.simplify(s.exp(s.I*s.sqrt(s.Integer(-1))*x)-s.exp(-x))
    evanescent_ratio_error = abs(outgoing_ratio(2., 1., 2.)-complex(beta_exact))
    rejected = {}
    for name, tangent, W in [('threshold', 0., 2.), ('evanescent', 1., 2.)]:
        try:
            flux_modes(2., tangent, W)
        except ValueError:
            rejected[name] = True
        else:
            rejected[name] = False
    rng = np.random.default_rng(71087)
    f = rng.normal(size=(2, 13))+1j*rng.normal(size=(2, 13))
    theta = .73
    n = np.array([np.sin(theta), np.cos(theta)])
    flux = boundary_flux(f, n)
    direct = np.real(np.einsum('ai,ab,bi->i', f.conj(), n[0]*SX+n[1]*SZ, f))
    projector_error = float(np.max(abs(flux['net']-direct)))
    characteristic = []
    mode_error = 0.
    for velocity in (1., .8, .2, .05):
        energy = 2.
        transverse = energy*np.sqrt(1-velocity**2)
        tangent, W = .6*transverse, .8*transverse
        p, modes = flux_modes(energy, tangent, W)
        mode_error = max(mode_error, float(np.max(abs(modes.conj().T @ SZ @ modes-SZ))))
        for col, sign in ((0, 1), (1, -1)):
            mode_error = max(mode_error, float(np.max(abs((sign*p*SZ+tangent*SX+W*SY) @ modes[:, col]-energy*modes[:, col]))))
        beta = outgoing_ratio(energy, tangent, W)
        r_exact = -(modes[1, 0]-beta*modes[0, 0])/(modes[1, 1]-beta*modes[0, 1])
        r_local = -modes[1, 0]/modes[1, 1]
        characteristic.append({'normal_group_velocity': velocity,
            'local_characteristic_reflection': float(abs(r_local)**2),
            'analytic_reflection': (1-velocity)/(1+velocity),
            'exact_outgoing_reflection': float(abs(r_exact)**2)})
    normal = cap_scattering(2., 0., 0., 4., 1., wall_phase=.37)
    normal_exact = np.exp(2j*2.*4.+.37j-2*4./4.)
    normal_error = abs(complex(*normal['r'])-normal_exact)
    constant = cap_scattering(2., .6, .8, 2., .7, power=0, wall_phase=.31)
    _, modes = flux_modes(2., .6, .8)
    wall = np.array([1., np.exp(.31j)])
    exact_start = expm(-2*(1j*(2.+.7j)*SZ+.6*SY-.8*SX)) @ wall
    a, b = np.linalg.solve(modes, exact_start)
    constant_error = abs(complex(*constant['r'])-b/a)
    smooth = [cap_scattering(2., .6, .8, length, 1.) for length in (4., 8., 16.)]
    refined = cap_scattering(2., .6, .8, 16., 1., rtol=2e-12)
    refinement = abs(complex(*smooth[-1]['r'])-complex(*refined['r']))
    grazing = cap_scattering(2., 1.2*np.sqrt(.96), 1.6*np.sqrt(.96), 16., 1.)
    # Exact free propagation independently tests preserved transmission phase.
    p, modes = flux_modes(2., .6, .8)
    length = 1.7
    evolution = expm(length*(2j*SZ+.6*SY-.8*SX)) @ modes[:, 0]
    transmitted, reflected = modal_amplitudes(evolution, 2., .6, .8)
    phase_error = max(abs(transmitted-np.exp(1j*p*length)), abs(reflected))
    cases = [normal, constant, *smooth, refined, grazing]
    tests = {'full_variable_W_operator_square': operator_residual == ['0', '0'],
        'boundary_characteristic_identity': projector_error < 1e-13,
        'spinor_eigenvectors_and_unit_flux': mode_error < 1e-12,
        'local_boundary_reflection_formula': max(abs(r['local_characteristic_reflection']-r['analytic_reflection']) for r in characteristic) < 1e-13,
        'exact_outgoing_map': max(r['exact_outgoing_reflection'] for r in characteristic) < 1e-25,
        'exact_evanescent_branch': evanescent_residual == [0, 0] and evanescent_current == 0 and evanescent_decay == 0 and evanescent_ratio_error < 1e-15,
        'nonpropagating_unit_flux_rejected': all(rejected.values()),
        'normal_incidence_complex_CAP_solution': normal_error < 2e-10,
        'constant_CAP_independent_matrix_exponential': constant_error < 2e-10,
        'independently_integrated_loss_balance': max(r['flux_balance_error'] for r in cases) < 2e-8,
        'complex_amplitude_refinement': refinement < 2e-9,
        'free_transmission_complex_phase': phase_error < 1e-13,
        'smooth_CAP_calibration': smooth[-1]['reflection_probability'] < 1e-6 and smooth[-1]['reflection_probability'] < smooth[0]['reflection_probability']}
    tests = {key: bool(value) for key, value in tests.items()}
    return {'passed': all(tests.values()), 'tests': tests,
        'errors': {'projector': projector_error, 'flux_modes': mode_error,
                   'normal_complex_amplitude': float(normal_error), 'constant_complex_amplitude': float(constant_error),
                   'refinement_complex_amplitude': float(refinement), 'free_transmission_phase': float(phase_error)},
        'characteristic_boundary': characteristic, 'normal': normal, 'constant': constant,
        'evanescent': {'energy': 2, 'tangent': 1, 'W': 2, 'p': 'I',
            'beta': str(beta_exact), 'operator_residual': [str(v) for v in evanescent_residual],
            'normal_current': str(evanescent_current), 'spatial_factor': 'exp(-x)',
            'decay_identity_residual': str(evanescent_decay), 'ratio_error': evanescent_ratio_error},
        'unit_flux_rejected': rejected,
        'smooth': smooth, 'refined': refined, 'grazing': grazing,
        'scope': 'Positive-energy, frozen planar spinor boundary calibration only. No physical return matrix, infinity flux, finite-mouth matching, vacuum state or escape probability is computed.',
        'phase_convention': 'Modes from flux_modes at CAP entrance x=0; wall b=exp(i*wall_phase)*a. Complex r stored as [real,imag].',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    result = checks()
    Path(__file__).with_name('open_boundary_spinor_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))
    if not result['passed']:
        raise SystemExit('Boundary benchmark failed; diagnostics retained.')
