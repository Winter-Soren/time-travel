"""Globally conforming compact trial fields on the actual sewn geometry.

Exterior tensor-polynomial bumps are lifted through both exact mouth traces.
Their physical Euclidean forms reuse137. No absorbing boundary, polygonal cut,
or projected first-order square is used. This is a finite variational pilot;
neither a continuum source nor the local-subtracted regulator limit is proved.
"""
import gc
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.integrate import simpson
from scipy.linalg import eigh

from balanced_dihole import symbolic_fields
from dihole_curved_exterior import coefficient_values
from dihole_euclidean_form import assemble, mask_derivative, relative_error
from dihole_bridge_scattering import collar_mask
from dihole_stretched_collar import exact_target_parameters, field_values, gauss, radius_profile
from dihole_sewn_quadrature import exterior_nodes, self_check as quadrature_check
from project_paths import project_path


RHO_BOUNDS = (.005, .5)
Z_BOUNDS = (-.85, .85)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.])
SX4, SY4, SZ4 = [np.kron(np.diag([1., -1.]), sigma) for sigma in (SX, SY, SZ)]


def scalar_trials(points, degrees=(2, 3)):
    """Double-zero tensor Legendre bumps and their exact spatial derivatives."""
    pairs = []
    for coordinate, bounds, number in zip(np.asarray(points).T, (RHO_BOUNDS, Z_BOUNDS), degrees):
        lo, hi = bounds
        x = 2*(coordinate-lo)/(hi-lo)-1
        inside = (x > -1) & (x < 1)
        bump = np.where(inside, (1-x*x)**2, 0.)
        bump_derivative = np.where(inside, -4*x*(1-x*x), 0.)
        value, derivative = [], []
        for n in range(number):
            poly = Legendre.basis(n)
            value.append(bump*poly(x))
            derivative.append((bump_derivative*poly(x)+bump*poly.deriv()(x))*2/(hi-lo))
        pairs.append((np.array(value).T, np.array(derivative).T))
    (r, dr), (z, dz) = pairs
    return tuple(np.einsum('pi,pj->pij', a, b).reshape(len(points), -1)
                 for a, b in ((r, z), (dr, z), (r, dz)))


def spin_values(scalars):
    return np.einsum('pj,ab->pabj', scalars, np.eye(4)).reshape(len(scalars), 4, -1).astype(complex)


def matrix_action(matrix, value):
    return np.einsum('ab,pbj->paj', matrix, value)


def exterior_data(mu, order=24, degrees=(2, 3)):
    points, weights = exterior_nodes(order, rho_bounds=RHO_BOUNDS, z_bounds=Z_BOUNDS)
    rho = points[:, 0, None, None]
    fields = coefficient_values(points)
    value, dr, dz = [spin_values(x) for x in scalar_trials(points, degrees)]
    radial = -1j*(matrix_action(SX4, dr+value/(2*rho))+matrix_action(SZ4, dz))
    angular = matrix_action(SY4, value)*(fields['NP']*(mu-fields['Aphysical']))[:, None, None]/rho
    return dict(value=value, radial=radial, angular=angular, P=fields['P'], W=fields['W'],
                f=np.zeros(len(points)), fp=np.zeros(len(points)), weight=weights*points[:, 0])


def trace_values(theta, side, degrees=(2, 3)):
    """b=v_b/sqrt(sin): exact cut trace and theta derivative in bridge frame."""
    sign = 1 if side == 'right' else -1
    points = np.column_stack((.05*np.sin(theta), sign*(.5+.05*np.cos(theta))))
    value, dr, dz = [spin_values(x) for x in scalar_trials(points, degrees)]
    c, s = np.cos(theta/2), np.sin(theta/2)
    rotation = np.zeros((len(theta), 4, 4), complex)
    small = np.array([[c, sign*s], [-sign*s, c]]).transpose(2, 0, 1)
    rotation[:, :2, :2] = rotation[:, 2:, 2:] = small
    derivative = .05*(np.cos(theta)[:, None, None]*dr-sign*np.sin(theta)[:, None, None]*dz)
    derivative += .5j*sign*matrix_action(np.kron(np.eye(2), SY), value)
    return tuple(.05*np.einsum('pab,pbj->paj', rotation, v) for v in (value, derivative))


def bridge_data(mu, phase, order=32, degrees=(2, 3)):
    function, parameters = symbolic_fields()[0], exact_target_parameters()
    nodes, weights = gauss(order)
    maximum_cosine = np.sqrt(1-(RHO_BOUNDS[0]/.05)**2)
    cosine, theta_weights = maximum_cosine*nodes, maximum_cosine*weights
    theta = np.arccos(cosine)
    left, left_derivative = trace_values(theta, 'left', degrees)
    right, right_derivative = [np.exp(-1j*phase)*v for v in trace_values(theta, 'right', degrees)]
    rows = {key: [] for key in ('value', 'radial', 'angular', 'P', 'W', 'f', 'fp', 'weight')}
    for lo, hi in ((-.1, -.05), (-.05, .05), (.05, .1)):
        for s, weight in zip((hi+lo)/2+(hi-lo)*nodes/2, (hi-lo)*weights/2):
            t = (s+.1)/.2
            value = (1-t)*left+t*right
            derivative_theta = (1-t)*left_derivative+t*right_derivative
            derivative_s = (right-left)/.2
            radius = float(radius_profile(min(.05, max(0., abs(s)-.05)))[0])
            fields = field_values(function, parameters, radius, cosine)
            lapse, meridian, potential = [np.asarray(fields[k], float) for k in
                                         ('lapse', 'meridian_factor', 'unit_A_phi')]
            angular = (-1j*matrix_action(SX4, derivative_theta
                        +cosine[:, None, None]/(2*np.sin(theta)[:, None, None])*value)
                       +matrix_action(SY4, value)
                       *(lapse*meridian*(mu-potential)/np.sin(theta))[:, None, None])/radius
            rows['value'].append(value)
            rows['radial'].append(-1j*matrix_action(SZ4, derivative_s))
            rows['angular'].append(angular)
            for key, data in (('P', meridian), ('W', meridian/lapse),
                              ('f', np.full(len(theta), collar_mask(s))),
                              ('fp', np.full(len(theta), mask_derivative(s))),
                              ('weight', weight*theta_weights)):
                rows[key].append(data)
    return {key: np.concatenate(value) for key, value in rows.items()}


def combine(exterior, bridge):
    return {key: np.concatenate((exterior[key], bridge[key])) for key in exterior}


def frequency_polynomial(data):
    """The fibre form and variation are exactly quadratic in real kappa."""
    zero, plus, minus = [assemble(data, k) for k in (0., 1., -1.)]
    inverse = np.linalg.inv(np.linalg.cholesky(zero['norm']))
    result = {}
    for key in ('q', 'derivative'):
        coefficients = (zero[key], (plus[key]-minus[key])/2,
                        (plus[key]+minus[key])/2-zero[key])
        result[key] = [inverse@matrix@inverse.conj().T for matrix in coefficients]
    return result


def fibre_response(polynomial, kappa, mass=1.):
    p, dp = [sum(kappa**j*coefficient for j, coefficient in enumerate(polynomial[key]))
             for key in ('q', 'derivative')]
    eigenvalues, vectors = eigh(p)
    if eigenvalues[0] <= 0:
        raise ArithmeticError('This finite pencil is not positive definite; inspect its zero modes.')
    inserted = np.einsum('ij,ij->j', vectors.conj(), dp@vectors).real
    rational = 6*mass**6/np.prod(eigenvalues[:, None]+mass**2*np.arange(4), axis=1)
    return -.5*float(inserted@rational)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), layout='constrained')
    rho, z = np.meshgrid(np.linspace(*RHO_BOUNDS, 121), np.linspace(*Z_BOUNDS, 241))
    value = scalar_trials(np.column_stack((rho.ravel(), z.ravel())))[0][:, 0].reshape(rho.shape)
    value = np.ma.masked_where((rho*rho+(z-.5)**2 < .05**2)
                              |(rho*rho+(z+.5)**2 < .05**2), value)
    axes[0].contourf(rho, z, value, levels=12, cmap='Blues')
    for center in (-.5, .5):
        axes[0].add_patch(Circle((0, center), .05, color='white', ec='#ba6b17', lw=2))
    axes[0].annotate('Exact mouth cuts', xy=(.04, .52), xytext=(.16, .67),
                     arrowprops=dict(arrowstyle='->'), fontsize=10)
    axes[0].text(.2, -.03, 'Field also extends\nthrough the bridge', ha='center', fontsize=10)
    axes[0].set(xlabel='Exterior rho / D', ylabel='Exterior z / D',
                title='One global exterior trial amplitude', xlim=(0, .5))
    phase_names = ['0', 'pi']
    base = [next(r['integral_33'] for r in result['records'] if r['mu'] == .5 and r['phase'] == phase)
            for phase in (0., np.pi)]
    enriched = [r['integral_33'] for r in result['basis_refinement']]
    for values, label, marker in ((base, '24 global trials', 'o'), (enriched, '48 global trials', 's')):
        axes[1].plot(phase_names, values, marker=marker, label=label)
    axes[1].axhline(0, color='grey', lw=.7)
    axes[1].set(xlabel='Handle phase', ylabel='Finite-window response (hbar c / D)',
                title='Basis sensitivity: mu = 1/2, MD = 1')
    axes[1].legend()
    fig.suptitle('Fields now cross both openings; supporting stress remains undetermined', fontsize=12)
    fig.supxlabel('Conforming finite trial spaces; omitted fields and frequency tails are not bounded.', fontsize=9)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_sewn_form_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    start = perf_counter()
    quadrature = quadrature_check()
    records, diagnostics, central = [], [], {}
    for mu in (-.5, .5, 1.5):
        exterior = exterior_data(mu, order=24)
        for phase in (0., np.pi):
            bridge = bridge_data(mu, phase, order=48)
            data = combine(exterior, bridge)
            base = assemble(data, .7, mass=1.3)
            plus, minus = [assemble(data, .7, epsilon=e) for e in (1e-5, -1e-5)]
            polynomial = frequency_polynomial(data)
            if mu == .5:
                central[phase] = dict(form=base, polynomial=polynomial)
            check_q = sum(.7**j*c for j, c in enumerate(polynomial['q']))
            inverse = np.linalg.inv(np.linalg.cholesky(base['norm']))
            varied_mass = assemble(data, .7, epsilon=.001, mass=1.3)
            diagnostics.append(dict(mu=mu, phase=phase, dimension=len(base['norm']),
                global_skew_error=relative_error(base['weak_dirac'], -base['weak_dirac'].conj().T),
                mass_shift_error=relative_error(base['direct_massive'], base['q']+1.3**2*base['norm']),
                varied_mass_shift_error=relative_error(varied_mass['direct_massive'], varied_mass['q']+1.3**2*varied_mass['norm']),
                derivative_fd_error=relative_error(base['derivative'], (plus['q']-minus['q'])/2e-5),
                frequency_polynomial_error=relative_error(check_q, inverse@base['q']@inverse.conj().T),
                min_eigenvalue=float(eigh(base['q'], base['norm'], eigvals_only=True)[0]),
                norm_condition=float(np.linalg.cond(base['norm'])), exterior_nodes=len(exterior['P']),
                bridge_nodes=len(bridge['P'])))
            kappas = np.geomspace(.01, 10., 33)
            responses = [fibre_response(polynomial, k) for k in kappas]
            records.append(dict(mu=mu, phase=phase, kappas=kappas.tolist(), responses=responses,
                                integral_33=float(simpson(responses, x=kappas)/np.pi),
                                integral_17=float(simpson(responses[::2], x=kappas[::2])/np.pi)))
            del data, bridge
            gc.collect()
        del exterior
    # One shared geometry/phase control resolves integration error separately
    # from basis truncation; it is not an exterior or spectral remainder.
    refined = combine(exterior_data(.5, order=32), bridge_data(.5, 0., order=64))
    refined_form = assemble(refined, .7)
    refined_polynomial = frequency_polynomial(refined)
    original_form = central[0.]['form']
    refinement = dict(form_change=relative_error(refined_form['q'], original_form['q']),
                      derivative_change=relative_error(refined_form['derivative'], original_form['derivative']),
                      norm_change=relative_error(refined_form['norm'], original_form['norm']),
                      refined_skew_error=relative_error(refined_form['weak_dirac'], -refined_form['weak_dirac'].conj().T),
                      fibre_response_change=abs(fibre_response(refined_polynomial, .7)-fibre_response(
                          central[0.]['polynomial'], .7)))
    del refined
    gc.collect()
    basis_refinement = []
    for phase in (0., np.pi):
        enlarged = combine(exterior_data(.5, order=24, degrees=(3, 4)),
                           bridge_data(.5, phase, order=48, degrees=(3, 4)))
        form = assemble(enlarged, .7)
        polynomial = frequency_polynomial(enlarged)
        old_indices = np.array([12*component+4*i+j for component in range(4)
                                for i in range(2) for j in range(3)])
        old = central[phase]['form']
        kappas = np.geomspace(.01, 10., 33)
        values = [fibre_response(polynomial, k) for k in kappas]
        basis_refinement.append(dict(phase=phase, dimension=len(form['norm']),
            nested_form_error=relative_error(form['q'][np.ix_(old_indices, old_indices)], old['q']),
            nested_norm_error=relative_error(form['norm'][np.ix_(old_indices, old_indices)], old['norm']),
            min_eigenvalue=float(eigh(form['q'], form['norm'], eigvals_only=True)[0]),
            global_skew_error=relative_error(form['weak_dirac'], -form['weak_dirac'].conj().T),
            integral_33=float(simpson(values, x=kappas)/np.pi),
            integral_17=float(simpson(values[::2], x=kappas[::2])/np.pi),
            kappas=kappas.tolist(), responses=values))
        del enlarged
        gc.collect()
    tests = dict(exact_domain_quadrature=quadrature['passed'],
        sewn_flux=all(r['global_skew_error'] < 1e-9 for r in diagnostics+basis_refinement),
        common_mass_shift=all(max(r['mass_shift_error'], r['varied_mass_shift_error']) < 1e-10 for r in diagnostics),
        metric_variation=all(r['derivative_fd_error'] < 2e-8 for r in diagnostics),
        frequency_polynomial=all(r['frequency_polynomial_error'] < 1e-12 for r in diagnostics),
        positive_finite_form=all(r['min_eigenvalue'] > 0 for r in diagnostics),
        quadrature_refinement=max(refinement[k] for k in ('form_change', 'derivative_change', 'norm_change')) < 2e-5,
        nested_trial_space=all(max(r['nested_form_error'], r['nested_norm_error']) < 1e-12 for r in basis_refinement),
        finite_responses=all(np.all(np.isfinite(r['responses'])) for r in records))
    names = ('dihole_sewn_quadrature.py', 'dihole_euclidean_form.py', 'dihole_curved_exterior.py',
             'dihole_stretched_collar.py', 'dihole_bridge_scattering.py', 'balanced_dihole.py')
    result = dict(complete=True, passed=all(tests.values()), tests=tests, diagnostics=diagnostics,
        refinement=refinement, basis_refinement=basis_refinement, quadrature=quadrature,
        records=records, seconds=perf_counter()-start,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
        scope=['Global H1 trial fields with nonzero exactly sewn traces and connected exterior support.',
               'Exact circular cuts; compact rectangle support, no absorbing numerical boundary.',
               'Six exterior scalar tensor bumps times four components, linear bridge trace lift.',
               'Finite physical-form Galerkin pilot at MD=1; three sectors and finite frequency window.',
               'No basis/exterior/axis, omitted angular, frequency-tail, IR, or regulator-removal error bound.',
               'Neither displayed sign nor magnitude is the renormalized supporting source.'])
    project_path('dihole_sewn_form_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise RuntimeError('Sewn physical form checks failed; evidence retained.')
    plot(result)
    return result


if __name__ == '__main__':
    run()
