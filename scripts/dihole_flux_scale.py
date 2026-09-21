"""Physical charge/size consistency for the prescribed Q/D=1/5 family.

The integer q is the charge-included magnetic flux, not a species count.
This converts an existing normalization constraint into dimensional tests;
it neither computes quantum stress nor establishes an attainable medium.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from project_paths import project_path, source_matches


NAME = 'dihole_flux_scale'
ALPHA_EXAMPLE = 7.2973525643e-3
PLANCK_LENGTH_M = 1.616255e-35
NIST = 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scale(q, alpha):
    """D/lP for nonzero integer flux and a declared positive coupling."""
    if isinstance(q, bool) or not isinstance(q, (int, np.integer)) or q == 0:
        raise ValueError('Require a nonzero integer flux q, not a species count.')
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError('Require a finite positive physical coupling alpha_R.')
    value = 2.5*abs(q)/math.sqrt(alpha)
    if not math.isfinite(value):
        raise ValueError('Scale exceeds the finite numerical range.')
    return value


def implied_alpha(q, d_over_planck):
    scale(q, 1.)  # Use the same flux validation.
    if not math.isfinite(d_over_planck) or d_over_planck <= 0:
        raise ValueError('Require a finite positive D/lP.')
    value = (2.5*abs(q)/d_over_planck)**2
    if not math.isfinite(value) or value == 0:
        raise ValueError('Coupling exceeds the finite positive numerical range.')
    return value


def load_inputs():
    records, hashes = {}, {}
    for name in ('balanced_dihole', 'dihole_local_source'):
        path = project_path(name+'_results.json')
        record = json.loads(path.read_text(encoding='utf-8'))
        if not record['passed'] or not source_matches(name+'.py', record['source_sha256']):
            raise ArithmeticError('Historical source/results provenance failed: '+name)
        for filename, digest in record.get('input_sha256', {}).items():
            if not source_matches(filename, digest):
                raise ArithmeticError('Historical input provenance failed: '+filename)
            hashes[filename] = sha(project_path(filename))
        hashes[path.name] = sha(path)
        hashes[name+'.py'] = sha(project_path(name+'.py'))
        records[name] = record
    hashes['project_paths.py'] = sha(project_path('project_paths.py'))
    return records, hashes


def controls(required):
    q, charge, gravitational, coupling = sp.symbols('q Q G alpha', positive=True)
    e_squared = 4*sp.pi*coupling
    # Maxwell coefficient after converting Fhat=(2Q/q)F_charge.
    physical_relation = {charge**2: gravitational*q**2/(4*coupling)}
    kinetic = ((2*charge/q)**2/(16*sp.pi*gravitational)).subs(physical_relation)
    d2 = 25*gravitational*q**2/(4*coupling)
    errors = []
    for flux in (-7, 1, 9, 512000):
        for alpha in (1e-8, ALPHA_EXAMPLE, .1):
            d = scale(flux, alpha)
            errors.append(abs(implied_alpha(flux, d)/alpha-1))
    rejected = 0
    for args in ((True, .1), (0, .1), (1.5, .1), (1, 0.), (1, float('nan'))):
        try:
            scale(*args)
        except ValueError:
            rejected += 1
    return {k: bool(v) for k, v in dict(
        canonical_Maxwell_coefficient=sp.simplify(kinetic-1/(4*e_squared)) == 0,
        required_source_substitution=sp.simplify(d2/gravitational-25*q**2/(4*coupling)) == 0,
        scale_coupling_roundtrip=max(errors) < 8e-16,
        orientation_independence=scale(-1, ALPHA_EXAMPLE)==scale(1, ALPHA_EXAMPLE),
        fixed_q_rescaling_changes_coupling=math.isclose(implied_alpha(1, 2*scale(1, ALPHA_EXAMPLE)), ALPHA_EXAMPLE/4, rel_tol=1e-15),
        species_not_implicitly_multiplied=math.isclose(scale(2, ALPHA_EXAMPLE)/scale(1, ALPHA_EXAMPLE), 2.),
        invalid_parameters_rejected=rejected == 5,
        required_target_negative=required < 0).items()}


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), layout='constrained')
    flux = np.geomspace(1., 1e34, 250)
    unit = result['illustration']['D_m']
    axes[0].loglog(flux, unit*flux, color='#267f97')
    for row in result['target_size_examples']:
        axes[0].plot(row['continuous_flux_requirement'], row['D_m'], 'o', color='#c07c27')
        axes[0].annotate(row['label'], (row['continuous_flux_requirement'], row['D_m']),
                         xytext=(-8, 8), textcoords='offset points', ha='right', fontsize=9)
    axes[0].plot(1, unit, 'o', color='#267f97')
    axes[0].annotate('Current q = 1', (1, unit), xytext=(7, 7), textcoords='offset points', fontsize=9)
    axes[0].set(xlabel='Integer flux q (log scale; line guides the eye)',
                ylabel='Weyl-coordinate scale D (metres)', title='Fixed shape and illustrative fixed coupling')
    axes[0].grid(alpha=.2)
    axes[1].axis('off')
    example = result['illustration']
    text = ('What the normalization fixes\n\n'
        'D / Planck length = 5 |q| / (2 sqrt(alpha_R))\n\n'
        f"For q = 1 and alpha_R = {ALPHA_EXAMPLE:.8f}:\n"
        f"D = {example['D_over_planck']:.2f} Planck lengths\n"
        f"Central areal radius = {example['central_areal_radius_over_planck']:.2f} Planck lengths\n\n"
        'To change the size:\n'
        'change q, change the coupling, or change the shape.\n'
        'Each changes the physical model being tested.\n\n'
        'No quantum source value or support bound is supplied.')
    axes[1].text(.02, .96, text, va='top', fontsize=10, linespacing=1.45)
    fig.suptitle('The present field calculation cannot simply be enlarged to human size')
    fig.supxlabel('Low-energy alpha is used only as a fixed illustrative value; no running-coupling or electron model is asserted.', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext), dpi=160)
    plt.close(fig)


def run():
    records, inputs = load_inputs()
    geometry, local = records['balanced_dihole'], records['dihole_local_source']
    p = geometry['parameters']
    if p['q'] != 1 or p['Q']/p['D'] != .2:
        raise ArithmeticError('This conversion is for the saved q=1, Q/D=1/5 family.')
    required = local['required_extra_source']
    cuts = {r['epsilon']: r for r in geometry['cut_summaries'] if r['angular_order']==96}
    ratio = scale(1, ALPHA_EXAMPLE)
    d_m = ratio*PLANCK_LENGTH_M
    result = dict(passed=False, source_sha256=sha(__file__), input_sha256=inputs,
        parameters=dict(q=1, Q_over_D=.2, alpha_R=ALPHA_EXAMPLE,
                        alpha_status='Illustrative fixed coupling, numerically equal to low-energy CODATA alpha; not a high-energy prediction.'),
        constants=dict(alpha=ALPHA_EXAMPLE, alpha_standard_uncertainty=1.1e-12,
            planck_length_m=PLANCK_LENGTH_M, planck_length_standard_uncertainty_m=1.8e-40,
            source=NIST, retrieved_on='2026-09-22', metrology='Central values only; no uncertainty enclosure claimed.'),
        normalization=dict(action='(R-Fhat^2)/(16 pi G)', charge_included_connection='A=-q Fhat_potential/(2Q)',
            alpha_R='e_R^2/(4 pi) in hbar=c=1 canonical Maxwell units',
            scale='D/lP=5 |q|/(2 sqrt(alpha_R)) for Q/D=1/5',
            magnetic_flux='integral F_charge=2 pi q; q is integer flux, not species count',
            required_source='C_required=-(25 a/4) q^2/alpha_R, a=0.00680448466178458'),
        illustration=dict(D_over_planck=ratio, D_m=d_m, Q_over_planck=.2*ratio,
            central_areal_radius_over_planck=cuts[.025]['area_radius']*ratio,
            cut_areal_radius_over_planck=cuts[.05]['area_radius']*ratio,
            required_dimensionless_source=required*ratio**2),
        target_size_examples=[dict(label=label, D_m=size, continuous_flux_requirement=size/d_m,
            alpha_if_q_stays_one=implied_alpha(1, size/PLANCK_LENGTH_M),
            scope='Continuous requirement locates the scale; allowed q is integer and the source must be recomputed.')
            for label, size in (('D = 1 fm', 1e-15), ('D = 1 nm', 1e-9), ('D = 1 micrometre', 1e-6), ('D = 1 m', 1.))],
        required_q_squared_over_alpha_coefficient=-25*required/4,
        actual_physical_source=None, physical_support_bound=None,
        sources=['https://arxiv.org/pdf/hep-th/9906160', NIST],
        scope=[
            'This applies an existing normalization constraint; it is not a new physical law.',
            'No physical fermion species, mass or high-energy running coupling is identified.',
            'D is a Weyl-coordinate scale, not proper distance or the opening radius.',
            'Changing q changes the Dirac operator and mode population; q=1 source data cannot be scaled to large q.',
            'At fixed physical coupling, changing separation with fixed charge changes Q/D and requires a new geometry/source calculation.',
            'The source remains unknown. Even support would leave current, backreaction, stability and chronology unresolved.'])
    result['tests'] = controls(required)
    result['tests'].update(exact_current_provenance=all(sha(project_path(n))==v for n,v in inputs.items()),
        actual_radius_distinguished=cuts[.025]['area_radius']>.025 and cuts[.05]['area_radius']>.05,
        no_fabricated_source=result['actual_physical_source'] is None and result['physical_support_bound'] is None)
    result['passed'] = bool(all(result['tests'].values()))
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Physical flux-scale checks failed; evidence saved.')
    plot(result)
    print(json.dumps({k:result[k] for k in ('passed','tests','illustration','target_size_examples')}, indent=2))


if __name__ == '__main__':
    run()
