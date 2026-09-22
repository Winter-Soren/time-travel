"""All-order local magnetic remainder, separated from the unknown full heat.

The reference uses the actual collar field pointwise. It is not asserted to
equal the heat kernel of a varying field on the full curved manifold.
"""
import json

import numpy as np
from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_local_source import mask_jets
from dihole_regulator_heat_tail import matching_weights
from dihole_slice_channel import sha
from dihole_stretched_collar import gauss, smooth_step
from project_paths import project_path


NAME='dihole_magnetic_remainder'


def geometric_ceiling(iv, geometry):
    """Both collars, complete sphere, exact integral of f*(R')^2."""
    width=iv.mpf(1)/20
    angular=4*iv.log(2)-2
    N=iv.mpf(geometry['N_upper']);arx=iv.mpf(geometry['polar_flux_Rx_upper_D'])
    B=iv.mpf(geometry['magnetic_unit_flux_upper_D2'])
    moment=4*iv.pi*(width/6)*N*N*arx*arx*angular
    return dict(radial_mask_rate_integral_over_D=upper(width/6),
        angular_polar_factor_upper=upper(angular),weighted_perpendicular_field_upper_per_q2_D=upper(moment),
        magnetic_a4_upper_per_q2_D=upper(iv.mpf(8)/3*moment),
        magnetic_C_upper_per_q4_D5=upper(iv.mpf(16)/45*B*B*moment),
        magnetic_unit_flux_upper_D2=geometry['magnetic_unit_flux_upper_D2'])


def scalar_controls():
    rows=[];terms=200
    for x in map(mp.mpf,('1e-8','1e-4','.01','.1','1','3','10','100','1000')):
        g=mp.coth(x)-x/mp.sinh(x)**2
        deficit=2*x/3-g;ceiling=4*x**3/45
        partial=mp.fsum(4*x*(mp.pi*n)**2/(x*x+(mp.pi*n)**2)**2 for n in range(1,terms+1))
        tail=4*x/(mp.pi**2*terms)
        rows.append(dict(x=float(x),deficit_nonnegative=bool(deficit>=0),
            deficit_below_cubic=bool(deficit<=ceiling),deficit_ceiling_ratio=float(deficit/ceiling),
            partial_fraction_encloses=bool(partial<=g<=partial+tail)))
    radial=[]
    for order in (64,128):
        nodes,weights=gauss(order);u=.025*(nodes+1)
        rate=smooth_step(u/.05)[0];mask=mask_jets(u)[0]
        value=.025*np.dot(weights,mask*rate*rate)
        radial.append(dict(order=order,value=float(value),absolute_error=abs(float(value)-1/120)))
    angular=mp.quad(lambda x:(1-abs(x))/(1+abs(x)),[-1,0,1])
    return dict(scalar=rows,mask_rate=radial,angular_integral_error=float(abs(angular-(4*mp.log(2)-2))))


def run():
    mp.dps=100;iv=MPIntervalContext();iv.dps=65
    path=project_path('dihole_regulator_heat_tail_results.json');prior=json.loads(path.read_text())
    inputs={**prior['input_sha256'],path.name:sha(path),'dihole_regulator_heat_tail.py':prior['source_sha256']}
    if not prior['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing full-field heat inputs required.')
    geometry=json.loads(project_path('dihole_buffer_smoothing_results.json').read_text())['geometry']
    constants=geometric_ceiling(iv,geometry);controls=scalar_controls();records=[]
    C=iv.mpf(constants['magnetic_C_upper_per_q4_D5'])
    for previous in prior['records']:
        q=previous['q'];rows=[]
        for target in previous['matching_rows']:
            _,_,J1=matching_weights(iv,target['target_regulator_M_D'],previous['cutoff_s0_over_D_squared'])
            magnitude=C*q**4*J1/(32*iv.pi**2)
            heat=iv.mpf(target['long_heat_correction_upper_D'])
            rows.append(dict(nominal_M_squared_s0=target['nominal_M_squared_s0'],
                target_regulator_M_D=target['target_regulator_M_D'],
                pure_magnetic_short_interval_D=[lower(-magnitude),0.],
                pure_magnetic_short_magnitude_upper_per_q2=upper(magnitude/(q*q)),
                long_heat_upper_per_q2=target['long_heat_correction_upper_per_q2'],
                known_error_interval_D=[lower(-magnitude-heat),upper(heat)],
                known_error_lower_per_q2=lower((-magnitude-heat)/(q*q)),
                known_error_upper_per_q2=upper(heat/(q*q)),
                C_geometry_error_coefficient_upper_D4=target['Cstar_error_coefficient_upper_D4']))
        records.append(dict(q=q,cutoff_s0_over_D_squared=previous['cutoff_s0_over_D_squared'],
            magnetic_C_upper_D5=upper(C*q**4),matching_rows=rows,
            actual_curved_remainder_C_geometry=None,finite_target_regulator_source=None,
            renormalized_source=None,stability=None,atom_transit=None))
    local=json.loads(project_path('dihole_local_source_results.json').read_text())
    a4_estimate=-local['records'][-1]['charge_included_F2_energy_derivative']*2/3
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        all_order_magnetic_sign=all(r['deficit_nonnegative'] for r in controls['scalar']),
        all_order_cubic_bound=all(r['deficit_below_cubic'] for r in controls['scalar']),
        independent_partial_fractions=all(r['partial_fraction_encloses'] for r in controls['scalar']),
        actual_mask_rate_integral=max(r['absolute_error'] for r in controls['mask_rate'])<2e-12,
        polar_integral=controls['angular_integral_error']<1e-85,
        positive_geometry_constants=all(v>0 for v in constants.values()),
        saved_a4_estimate_within_ceiling=0<a4_estimate<constants['magnetic_a4_upper_per_q2_D'],
        signed_magnetic_intervals=all(row['pure_magnetic_short_interval_D'][0]<0==row['pure_magnetic_short_interval_D'][1] for r in records for row in r['matching_rows']),
        magnetic_bound_decreases=all(all(a['pure_magnetic_short_magnitude_upper_per_q2']>b['pure_magnetic_short_magnitude_upper_per_q2'] for a,b in zip(r['matching_rows'],r['matching_rows'][1:])) for r in records),
        known_error_contains_both_signs=all(row['known_error_interval_D'][0]<0<row['known_error_interval_D'][1] for r in records for row in r['matching_rows']),
        missing_physics_remains_explicit=all(r['actual_curved_remainder_C_geometry'] is None and r['finite_target_regulator_source'] is None and r['renormalized_source'] is None for r in records))
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,geometry_constants=constants,controls=controls,
        records=records,local_a4_gauge_numerical_estimate=a4_estimate,input_sha256=inputs,source_sha256=sha(__file__),
        reference='Pointwise constant-magnetic-field response deltaTheta0 from162, evaluated with the actual collar B and metric variation.',
        remaining_remainder='deltaK-deltaTheta0-H_gravity-a6*s/(16*pi^2); its bound C_geometry is unknown.',
        scope=['The pure magnetic reference correction is nonpositive and bounded to all orders in s*B. No truncation or small-s*B assumption enters that inequality.',
               'Actual directed collar bounds, polar regularity and an exact mask/radial-rate integral supply the prefactor; no geometry quadrature is promoted to a certificate.',
               'The known error interval combines only the magnetic reference short correction and177 long heat term; it is not a full source interval.',
               'Curvature, field gradients and nonlocal effects remain in C_geometry. Removing pure F^4 does not imply a uniform q^2 bound on this remainder.',
               'The target-M source and same physical counterterms from177 are still required. No support, material recipe, stability, atom passage or time gain is established.'],
        sources=['https://dlmf.nist.gov/4.36.E3','https://arxiv.org/abs/0810.4889','https://arxiv.org/pdf/0906.2430'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Magnetic remainder controls failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,constants=constants,atomic=records[-1]['matching_rows']),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=result['records'][-1]['matching_rows'];z=[r['nominal_M_squared_s0'] for r in rows]
    fig,ax=plt.subplots(figsize=(8.5,4.6),layout='constrained')
    ax.semilogy(z,[r['long_heat_upper_per_q2'] for r in rows],'o-',label='Bounded long-time heat correction')
    ax.semilogy(z,[r['pure_magnetic_short_magnitude_upper_per_q2'] for r in rows],'s-',label='Bounded magnetic short-time reference')
    ax.set(xlabel='Target regulator squared times proper-time cutoff',ylabel='Correction magnitude ceiling times D / q²',
           title='Two controlled pieces; curved-field remainder still unknown')
    ax.legend(fontsize=9);ax.grid(alpha=.2)
    fig.supxlabel('Actual atomic-flux reference. The curves are error ceilings, not supporting stress.\nThe total source, stable opening and atom passage remain unproved.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
