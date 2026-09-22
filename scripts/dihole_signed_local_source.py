"""Signed enclosure of the local magnetic reference, not the full quantum source."""
import json

import numpy as np
from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_slice_channel import sha
from dihole_strong_field_response import local_response_per_q2
from project_paths import project_path


NAME = 'dihole_signed_local_source'


def enclosure(iv, q, field_lower, field_upper, weight_upper):
    if q <= 0 or field_lower <= 0 or field_upper < field_lower or weight_upper < 0:
        raise ValueError('Positive flux/field bounds and nonnegative weight required.')
    constant = 12*iv.log(iv.glaisher)-iv.mpf(1)/2-2*iv.euler
    lo = constant+iv.log(iv.mpf(q)/iv.mpf(field_upper))
    hi = constant+iv.log(iv.mpf(q)/iv.mpf(field_lower))
    coefficient = iv.mpf(weight_upper)/(12*iv.pi**2)
    # Only an upper bound on the positive weight is available.
    left = coefficient*(lo if lo.a < 0 else iv.mpf(0))
    right = coefficient*(hi if hi.b > 0 else iv.mpf(0))
    return dict(q=q, bracket_interval=[lower(lo), upper(hi)],
        local_response_per_q2_interval=[lower(left) if lo.a < 0 else 0.,
                                        upper(right) if hi.b > 0 else 0.],
        uniformly_nonnegative=bool(lo.a >= 0),
        uniformly_nonpositive=bool(hi.b <= 0),
        full_source=None, global_remainder=None, atom_transit=None)


def controls():
    mp.dps = 100
    # Functional equation at s=-1 supplies an independent constant evaluation.
    direct = 12*mp.log(mp.glaisher)-mp.mpf(1)/2-2*mp.euler
    zeta = mp.log(2*mp.pi)-mp.euler-mp.mpf(1)/2-mp.diff(mp.zeta, 2)/mp.zeta(2)
    q, field, weight = mp.mpf(1000), mp.mpf(8), mp.mpf('.002')
    expected = weight*(zeta+mp.log(q/field))/(12*mp.pi**2)
    reused = local_response_per_q2(float(q), float(q*mp.sqrt(2)*mp.exp(-mp.euler)),
        np.array([float(field)]), np.array([float(weight/field)]),
        np.array([1.]), float(mp.log(mp.glaisher)))
    return dict(zeta_constant_error=float(abs(direct-zeta)),
        existing_reference_formula_relative_error=float(abs(mp.mpf(reused)/expected-1)))


def run():
    names = ('dihole_heat_source_graph', 'dihole_magnetic_remainder',
             'dihole_strong_field_geometry', 'dihole_strong_field_response')
    saved, inputs = {}, {}
    for name in names:
        path = project_path(name+'_results.json')
        row = json.loads(path.read_text())
        if not row['passed']:
            raise ArithmeticError('Passing input required: '+name)
        current = {**row['input_sha256'], path.name:sha(path), name+'.py':row['source_sha256']}
        for filename, digest in current.items():
            if filename in inputs and inputs[filename] != digest:
                raise ArithmeticError('Conflicting provenance: '+filename)
            inputs[filename] = digest
        saved[name] = row
    if not all(sha(project_path(n)) == h for n,h in inputs.items()):
        raise ArithmeticError('Exact current inputs required.')
    iv = MPIntervalContext(); iv.dps = 100
    geometry = saved['dihole_strong_field_geometry']['interval_zero_exclusion']['records'][-1]
    if not geometry['passed'] or geometry['failed_tiles']:
        raise ArithmeticError('Closed-collar field enclosure required.')
    field_lower = geometry['q1_B_D0_squared_lower']
    constants = saved['dihole_magnetic_remainder']['geometry_constants']
    field_upper = constants['magnetic_unit_flux_upper_D2']
    weight = constants['weighted_perpendicular_field_upper_per_q2_D']
    qs = sorted({r['q'] for r in saved['dihole_heat_source_graph']['records']})
    records = [enclosure(iv,q,field_lower,field_upper,weight) for q in qs]
    sign_controls = [enclosure(iv,q,field_lower,field_upper,weight) for q in (1,2,3,4)]
    constant = 12*iv.log(iv.glaisher)-iv.mpf(1)/2-2*iv.euler
    positive_threshold = iv.mpf(field_upper)*iv.exp(-constant)
    negative_threshold = iv.mpf(field_lower)*iv.exp(-constant)
    ref = saved['dihole_strong_field_response']
    numerical = ref['samples'][-1]['atomic_fixed_physical_reference_response_per_q2']
    target = ref['required_source_per_q2']
    atomic = records[-1]
    # The target is an older quadrature estimate; do not certify its uncertainty.
    needed = iv.mpf(target)-iv.mpf(atomic['local_response_per_q2_interval'])
    diagnostic = dict(numerical_local_source_per_q2=numerical,
        numerical_target_per_q2=target,
        necessary_omitted_contribution_if_target_exact=[lower(needed),upper(needed)],
        target_geometry_uncertainty_certified=False,
        actual_omitted_contribution=None)
    control = controls()
    invalid = 0
    for args in ((0,1,2,1),(1,0,2,1),(1,2,1,1),(1,1,2,-1)):
        try: enclosure(iv,*args)
        except ValueError: invalid += 1
    checks = dict(exact_inputs=all(sha(project_path(n)) == h for n,h in inputs.items()),
        zeta_constant_identity=control['zeta_constant_error'] < 1e-90,
        existing_formula_agreement=control['existing_reference_formula_relative_error'] < 1e-14,
        invalid_inputs_rejected=invalid == 4,
        entire_collar_field_bounds=field_lower > 0 and field_upper >= field_lower,
        saved_flux_references_nonnegative=all(r['uniformly_nonnegative'] for r in records),
        low_flux_sign_transition=sign_controls[0]['uniformly_nonpositive']
            and sign_controls[-1]['uniformly_nonnegative'],
        numerical_reference_inside=atomic['local_response_per_q2_interval'][0] <= numerical
            <= atomic['local_response_per_q2_interval'][1],
        missing_physics_explicit=all(r['full_source'] is None and r['global_remainder'] is None
            and r['atom_transit'] is None for r in records))
    result = dict(passed=all(checks.values()), checks=checks, controls=control,
        field_unit_interval=[field_lower,field_upper], positive_weight_upper=weight,
        constant_interval=[lower(constant),upper(constant)],
        all_point_nonnegative_q_threshold_interval=[lower(positive_threshold),upper(positive_threshold)],
        all_point_nonpositive_q_threshold_interval=[lower(negative_threshold),upper(negative_threshold)],
        records=records, sign_transition_controls=sign_controls, diagnostic=diagnostic,
        input_sha256=inputs, source_sha256=sha(__file__),
        units='Local reference energy variation divided by q^2 in hbar*c/D; q is exact integer flux.',
        fixed_reference='mu_H*D=q*sqrt(2)*exp(-Euler); same physical reference and illustrative coupling as162.',
        scope=['Entire-collar enclosure of the defined constant-magnetic-field reference only.',
            'The actual curved, varying-field global heat response is not this reference.',
            'Old geometry quadratures remain diagnostics, not certified source/target intervals.',
            'No signed global remainder, self-consistent support, stability, atom passage or clock gain.'],
        sources=['https://arxiv.org/abs/0906.2430'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:
        raise ArithmeticError('Signed local reference checks failed; result retained.')
    plot(result)
    print(json.dumps({k:result[k] for k in ('passed','checks','records','diagnostic')},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row, diagnostic = result['records'][-1], result['diagnostic']
    lo, hi = row['local_response_per_q2_interval']
    fig, axes = plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    axes[0].hlines(0,lo,hi,color='#267c84',lw=7)
    axes[0].plot(diagnostic['numerical_local_source_per_q2'],0,'o',color='#222222')
    axes[0].set(xlim=(-.04*hi,1.06*hi),ylim=(-.8,.8),yticks=[],
        title='Local magnetic reference: certified nonnegative',
        xlabel='Reference response / q squared')
    target = diagnostic['numerical_target_per_q2']
    axes[1].axvline(0,color='#777777',lw=1)
    axes[1].plot(target,0,'x',color='#b14343',ms=10,label='Required response (numerical)')
    axes[1].plot(0,0,'o',color='#267c84',label='Local reference interval near zero')
    axes[1].set(xlim=(1.2*target,abs(target)*.12),ylim=(-.8,.8),yticks=[],
        title='The missing contribution must be negative',
        xlabel='Response / q squared')
    axes[1].legend(fontsize=8,loc='upper center')
    fig.supxlabel('Units: hbar c / D. Defined local reference only; omitted global response unknown.\nNo supported opening, atom crossing or time shift is established.',fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()

