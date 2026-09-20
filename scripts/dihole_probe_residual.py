"""Compact saddle trial tested in the actual scaled curved dihole operator.

The source/trial are localized away from the axis and horizons. Norms below
are ordinary quadrature estimates, not interval certificates for the curved
resolvent. No physical exterior domain, quantum stress or throat is supplied.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from balanced_dihole import symbolic_fields
from dihole_null_scaling import physical_coefficients
from null_saddle_leads import AU, AV, SX, SZ, SY, _hermite_functions
from null_saddle_resolvent import real_generator
from project_paths import project_path, source_matches


def cutoff(x, inner=5., outer=7.):
    """C2 compact quintic cutoff, plus its first derivative."""
    if not 0 < inner < outer:
        raise ValueError('Require 0<inner<outer.')
    x = np.asarray(x, float)
    t = np.clip((abs(x)-inner)/(outer-inner), 0., 1.)
    value = 1-10*t**3+15*t**4-6*t**5
    derivative = (-30*t*t+60*t**3-30*t**4)*np.sign(x)/(outer-inner)
    return value, derivative


def trial_grid(a, b, order, inner=5., outer=7.):
    if outer < 7.:
        raise ValueError('Trial integration must contain the fixed source support |u|,|v|<=7.')
    nodes, weights = np.polynomial.legendre.leggauss(order)
    cuts = sorted(set((-outer, -inner, -7., -5., 5., 7., inner, outer)))
    x = np.concatenate([(lo+hi)/2+(hi-lo)*nodes/2 for lo, hi in zip(cuts, cuts[1:])])
    w = np.concatenate([(hi-lo)*weights/2 for lo, hi in zip(cuts, cuts[1:])])
    count = len(a)
    extended = _hermite_functions(x, count)
    basis = extended[:, :count]
    derivative = -extended[:, 1:count+1]*np.sqrt(np.arange(1, count+1)/2)
    derivative[:, 1:] += extended[:, :count-1]*np.sqrt(np.arange(1, count)/2)
    chi1, dchi1 = cutoff(x, inner, outer)
    chi, chiu, chiv = np.outer(chi1, chi1), np.outer(dchi1, chi1), np.outer(chi1, dchi1)
    source_chi1 = cutoff(x)[0]
    source_chi = np.outer(source_chi1, source_chi1)
    fields = []
    for vector in (a, b):
        value = np.einsum('pi,ijs,qj->pqs', basis, vector, basis, optimize=True)
        du = np.einsum('pi,ijs,qj->pqs', derivative, vector, basis, optimize=True)
        dv = np.einsum('pi,ijs,qj->pqs', basis, vector, derivative, optimize=True)
        fields.append((value, du, dv))
    source = np.zeros_like(fields[0][0])
    source[:, :, 0] = np.outer(basis[:, 0], basis[:, 0])
    u, v = np.meshgrid(x, x, indexing='ij')
    return dict(u=u, v=v, weights=np.outer(w, w), chi=chi, chi_u=chiu, chi_v=chiv,
                fields=fields, source=source, source_chi=source_chi,
                order=order, inner=inner, outer=outer, basis=basis)


def evaluate(grid, function, parameters, coefficients, q, kappa=1., return_residuals=False):
    if isinstance(q, bool) or int(q) != q or q <= 0:
        raise ValueError('Magnetic flux q must be a positive integer.')
    if parameters['q'] != 1. or parameters['physical_A_factor'] != -1/(2*parameters['Q']):
        raise ValueError('physical_coefficients must return unit-q a; do not count q twice.')
    c, p = coefficients, parameters
    length = c['width_prefactor']*q**(-1/3)
    energy = c['v_star']/length
    mu = math.floor(q*c['a_star']+energy/c['b_star'])+.5
    eta = c['b_star']*(mu-q*c['a_star'])/energy
    rho_min = c['rho_star']-math.sqrt(2)*length*grid['outer']
    if rho_min <= 0:
        raise ValueError('The compact rotated square reaches the physical axis.')
    u, v, chi = grid['u'], grid['v'], grid['chi']
    rho = c['rho_star']+length*(u+v)/math.sqrt(2)
    z = length*(u-v)/math.sqrt(2)
    actual = physical_coefficients(function, p, rho, z)
    speed = actual['v']/c['v_star']
    connection_u = length*(actual['v_rho']+actual['v_z'])/(2*math.sqrt(2)*c['v_star'])
    connection_v = length*(actual['v_rho']-actual['v_z'])/(2*math.sqrt(2)*c['v_star'])
    mass = actual['b']*(mu-q*actual['a'])/energy
    norm = lambda field: float(np.sqrt(np.sum(grid['weights']*np.sum(field**2, axis=2))))
    residuals, components, split_error = [], [], 0.
    # Real residuals for psi=i*a and adjoint chi=i*b. The imaginary energy
    # changes sign; the spatial generator is the same skew expression.
    for number, (value, du, dv) in enumerate(grid['fields']):
        sign = -1 if number == 0 else 1
        kinetic = du@AU.real.T+dv@AV.real.T
        local = grid['source']-kinetic+(1-2*u*v)[:, :, None]*(value@SY.imag.T)+sign*kappa*value
        geometry = (-(speed-1)[:, :, None]*kinetic
                    -connection_u[:, :, None]*(value@AU.real.T)
                    -connection_v[:, :, None]*(value@AV.real.T)
                    +(mass-(1-2*u*v))[:, :, None]*(value@SY.imag.T))
        edge = -speed[:, :, None]*(grid['chi_u'][:, :, None]*(value@AU.real.T)
                                  +grid['chi_v'][:, :, None]*(value@AV.real.T))
        source_difference = (grid['source_chi']-chi)[:, :, None]*grid['source']
        parts = [chi[:, :, None]*local, chi[:, :, None]*geometry, edge, source_difference]
        residual = sum(parts)
        # Independent direct action in the unrotated physical coordinates.
        cut_value = chi[:, :, None]*value
        cut_du = chi[:, :, None]*du+grid['chi_u'][:, :, None]*value
        cut_dv = chi[:, :, None]*dv+grid['chi_v'][:, :, None]*value
        dr = (cut_du+cut_dv)/(math.sqrt(2)*length)
        dz = (cut_du-cut_dv)/(math.sqrt(2)*length)
        direct = (grid['source_chi'][:, :, None]*grid['source']
                  -(actual['v'][:, :, None]*dr+actual['v_rho'][:, :, None]*cut_value/2)@SX.real.T/energy
                  -(actual['v'][:, :, None]*dz+actual['v_z'][:, :, None]*cut_value/2)@SZ.real.T/energy
                  +mass[:, :, None]*(cut_value@SY.imag.T)+sign*kappa*cut_value)
        split_error = max(split_error, float(np.max(abs(residual-direct))))
        residuals.append(residual)
        components.append(dict(saddle_residual=norm(parts[0]), curved_coefficient_mismatch=norm(parts[1]),
                               cutoff_derivative=norm(parts[2]), source_cutoff_difference=norm(parts[3]),
                               full_residual=norm(residual)))
    a, b = grid['fields'][0][0], grid['fields'][1][0]
    source = grid['source_chi'][:, :, None]*grid['source']
    Gtrial = float(np.sum(grid['weights']*np.sum(source*chi[:, :, None]*a, axis=2)))
    correction = -float(np.sum(grid['weights']*np.sum(chi[:, :, None]*b*residuals[0], axis=2)))
    radius = components[0]['full_residual']*components[1]['full_residual']/kappa
    row = dict(q=int(q), mu=mu, eta_actual=eta, eta_trial=1., length_scale=length,
        energy_scale=energy, kappa=kappa, inner=grid['inner'], outer=grid['outer'],
        minimum_support_rho=rho_min, quadrature_order_per_segment=grid['order'],
        source_norm=norm(source), trial_norm=norm(chi[:, :, None]*a),
        primal=components[0], adjoint=components[1], imaginary_G_trial=Gtrial,
        imaginary_G_correction=correction, corrected_imaginary_G=Gtrial+correction,
        numerical_dual_radius=radius, numerical_solution_error_ratio=components[0]['full_residual']/kappa,
        direct_physical_vs_split_error=split_error,
        independent_metric_speed_error=actual['independent_v_error'])
    return (row, residuals) if return_residuals else row


def correct_trial(a, b, function, calibration, q, steps=2, order=96, inner=5., outer=7.):
    """Bounded defect correction using the saddle inverse as a preconditioner.

    No convergence theorem is assumed: every iterate is tested again in the
    actual curved differential expression with the SAME compact source.
    """
    count = len(a)
    factor = splu(real_generator(count)-sparse.eye(2*count*count))
    a, b = a.copy(), b.copy()
    rows = []
    for step in range(steps+1):
        grid = trial_grid(a, b, order, inner, outer)
        row, residuals = evaluate(grid, function, calibration['parameters'], calibration['coefficients'], q, return_residuals=True)
        row['correction_step'] = step
        rows.append(row)
        print(json.dumps(dict(q=q, step=step, residual=row['primal']['full_residual'],
            corrected_G=row['corrected_imaginary_G'], radius=row['numerical_dual_radius'])), flush=True)
        if step == steps:
            break
        projections = [np.stack([grid['basis'].T@(grid['weights']*r[:, :, s])@grid['basis']
                                  for s in range(2)], axis=2) for r in residuals]
        a += factor.solve(-projections[0].ravel()).reshape(a.shape)
        b += factor.solve(projections[1].ravel(), trans='T').reshape(b.shape)
    return rows, a, b


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.9), layout='constrained')
    calibration = json.loads(project_path('dihole_null_scaling_results.json').read_text())
    center = calibration['coefficients']['rho_star']
    ax = axes[0]
    for row, color, label in ((result['rows'][8], '#64a7a0', 'Initial patch: q=8,000'),
                              (result['refined_corrected']['wide_q512000'], '#9769b3', 'Corrected patch: q=512,000')):
        radius = row['outer']*row['length_scale']*math.sqrt(2)
        vertices = [(center-radius, 0), (center, radius), (center+radius, 0), (center, -radius)]
        ax.add_patch(Polygon(vertices, facecolor=color, edgecolor=color, alpha=.32, label=label))
    ax.plot([center], [0], 'ko', ms=4)
    ax.axvline(0, color='#555555', lw=1)
    ax.set(xlim=(0, 1.2), ylim=(-.55, .55), xlabel=r'$\rho/D$', ylabel='z/D',
           title='Trial support inside the curved exterior', aspect='equal')
    ax.legend(fontsize=8, loc='lower right')
    ax = axes[1]
    for key, color, label in (('q8000', '#558e87', 'q=8,000; cutoff 5/7'),
                              ('q512000', '#536cb2', 'q=512,000; cutoff 5/7'),
                              ('wide_q512000', '#9862b0', 'q=512,000; cutoff 7/9')):
        rows = result['corrections'][key]
        ax.semilogy([x['correction_step'] for x in rows], [x['primal']['full_residual'] for x in rows],
                     'o-', color=color, label=label)
    ax.set(xlabel='Corrections using the curved-equation residual', ylabel='Residual norm (quadrature estimate)',
           title='Check every correction in the actual equations', xticks=[0, 1, 2])
    ax.legend(fontsize=8)
    ax = axes[2]; ax.axis('off')
    ax.text(0, .96, 'What improved', va='top', fontsize=13, weight='bold')
    ax.text(0, .84, 'Curved metric, gauge field and quantized mode\nare included in the differential equation.\nTwo corrections reduce the selected residual\nfrom 0.03978 to 0.001669 (about 24-fold).', va='top', fontsize=10.5, linespacing=1.6)
    ax.text(0, .53, 'What the check does not prove', va='top', fontsize=13, weight='bold')
    ax.text(0, .40, 'Quadrature agreement is not an interval certificate.\nThe exterior operator domain still needs completion.\nNo supporting quantum stress or stable throat.\nNo physical material or backward passage.', va='top', fontsize=10.5, linespacing=1.6)
    ax.text(0, .10, 'Source cutoff stays fixed when the trial is widened.', va='top', fontsize=10, color='#555555')
    fig.suptitle('From the ideal saddle to the actual curved coefficients | A numerical matching step', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path('dihole_probe_residual_progress.'+ext), dpi=175)
    plt.close(fig)


def run():
    saved = json.loads(project_path('null_saddle_resolvent_results.json').read_text())
    calibration = json.loads(project_path('dihole_null_scaling_results.json').read_text())
    assert saved['passed'] and calibration['passed']
    assert source_matches('null_saddle_resolvent.py', saved['source_sha256'])
    assert source_matches('dihole_null_scaling.py', calibration['source_sha256'])
    fields_path = project_path('null_saddle_resolvent_fields.npz')
    assert hashlib.sha256(fields_path.read_bytes()).hexdigest() == saved['field_sha256']
    with np.load(fields_path) as fields:
        a, b = fields['k1_n128_a'], fields['k1_n128_b']
    function = symbolic_fields()[0]
    rows = []
    for order in (64, 96, 144):
        grid = trial_grid(a, b, order)
        for q in (8000, 64000, 512000, 4096000):
            row = evaluate(grid, function, calibration['parameters'], calibration['coefficients'], q)
            rows.append(row)
            print(json.dumps({key: row[key] for key in ('q', 'quadrature_order_per_segment',
                'eta_actual', 'primal', 'corrected_imaginary_G', 'numerical_dual_radius')}), flush=True)
    corrections, refined, saved_fields = {}, {}, {}
    for key, q, inner, outer in (('q8000', 8000, 5., 7.), ('q512000', 512000, 5., 7.),
                                 ('wide_q512000', 512000, 7., 9.)):
        corrections[key], aa, bb = correct_trial(a, b, function, calibration, q, inner=inner, outer=outer)
        saved_fields[key+'_a'], saved_fields[key+'_b'] = aa, bb
        refined[key] = evaluate(trial_grid(aa, bb, 144, inner, outer), function,
                                calibration['parameters'], calibration['coefficients'], q)
    convergence = []
    for key in corrections:
        coarse, fine = corrections[key][-1], refined[key]
        convergence.append(dict(case=key, residual_norm_change=abs(coarse['primal']['full_residual']-fine['primal']['full_residual']),
            corrected_G_change=abs(coarse['corrected_imaginary_G']-fine['corrected_imaginary_G']),
            numerical_radius_change=abs(coarse['numerical_dual_radius']-fine['numerical_dual_radius'])))
    all_rows = rows+[r for group in corrections.values() for r in group]+list(refined.values())
    # Independent exact endpoint checks on the C2 polynomial.
    import sympy as sp
    t = sp.symbols('t')
    polynomial = 1-10*t**3+15*t**4-6*t**5
    cutoff_smooth = (polynomial.subs(t, 0) == 1 and polynomial.subs(t, 1) == 0
                     and all(sp.diff(polynomial, t, j).subs(t, endpoint) == 0 for j in (1, 2) for endpoint in (0, 1)))
    axis_guard = False
    try:
        evaluate(grid, function, calibration['parameters'], calibration['coefficients'], 1)
    except ValueError as error:
        axis_guard = 'axis' in str(error)
    tests = dict(compact_C2_cutoff=bool(cutoff_smooth), inadmissible_q1_rejected=axis_guard,
        physical_support=all(r['minimum_support_rho'] > 0 for r in all_rows),
        quantized_labels=all((2*r['mu']) % 2 == 1 and .5 <= r['mu'] <= r['q']-.5 for r in all_rows),
        detuning_retained=all(abs(r['eta_actual']-1) <= calibration['coefficients']['b_star']/(2*r['energy_scale'])+1e-12 for r in all_rows),
        direct_physical_action=all(r['direct_physical_vs_split_error'] < 2e-13 for r in all_rows),
        fixed_source_with_wider_trial=abs(refined['wide_q512000']['source_norm']-refined['q512000']['source_norm']) < 2e-13,
        corrected_quadrature_refinement=max(r['residual_norm_change'] for r in convergence) < 2e-8,
        correction_reduces_actual_residual=all(group[-1]['primal']['full_residual'] < group[0]['primal']['full_residual'] for group in corrections.values()),
        wider_trial_reduces_cutoff_residual=refined['wide_q512000']['primal']['cutoff_derivative'] < refined['q512000']['primal']['cutoff_derivative'])
    assert all(tests.values()), tests
    fields_output = project_path('dihole_probe_residual_fields.npz')
    np.savez_compressed(fields_output, **saved_fields)
    result = dict(passed=True, tests=tests, rows=rows, corrections=corrections,
        refined_corrected=refined, corrected_quadrature_comparisons=convergence,
        definitions=dict(source='h=chi_source*f0, unnormalized; source cutoff fixed at inner5/outer7 in every case.',
            trial='i*chi_trial*a; primal and adjoint corrections keep the same h.',
            scaled_operator='Hq=T^-1 Hphysical T/Eq; Tphi=Lq^-1 phi, flat canonical measure.',
            conditional_bound='For an agreeing self-adjoint realization: |G-Gtrial-<dual_trial,rho>|<=||rho||||sigma||/kappa.',
            physical_units='At physical frequency i*Eq*kappa, Green responses/radii for unitary-transformed h divide by Eq. If normalizing h, scalar values/radii also divide by ||h||^2.',
            scope='All curved integral norms and radii are quadrature estimates. This is not a certified curved-resolvent disk.'),
        limitations=['The natural full exterior spin/gauge/operator domain has not been completed here; no throat-cut domain is supplied.',
            'The compact trial reaches neither horizon nor a second mouth. This is a local response test, not demonstrated transport through the geometry.',
            'Q/D is fixed. Changing q requires recalibration or physical rescaling to preserve 4piG/g_R^2=(2Q/q)^2; no material preparation is inferred.',
            'No renormalized quantum stress, backreaction, stable throat, obtainable medium or backward travel.'],
        field_sha256=hashlib.sha256(fields_output.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name: hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                      for name in ('balanced_dihole.py', 'dihole_return_geometry.py', 'dihole_null_scaling.py', 'null_saddle_leads.py',
                                   'null_saddle_resolvent.py',
                                   'null_saddle_resolvent_results.json', 'null_saddle_resolvent_fields.npz',
                                   'dihole_null_scaling_results.json', 'project_paths.py')})
    project_path('dihole_probe_residual_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    draw(result)
    print(json.dumps(dict(passed=True, tests=tests, quadrature=convergence), indent=2))
    return result


if __name__ == '__main__':
    run()
