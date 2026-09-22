"""Recompute the same full defect in a source-centred local U(1) gauge.

The phase is extended globally, is time independent, and is held fixed during
the metric variation. Neither the curvature nor the physical model changes.
"""
import copy
import json
import math

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import assemble, unpack
from dihole_local_transport import packed, point_bound, transport_bounds
from dihole_regular_jets import operator_bounds
from dihole_scaled_graph import record_bound
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_centered_gauge'


def centred_jets(iv, prior, local, magnetic):
    radius=iv.mpf(2)**-32
    old=list(map(iv.mpf,prior['jets']['total_order_derivative_upper']['gauge_per_abs_q']))
    values=[4*radius*old[n+1]/(n+2)+2*n*old[n]/(n+1) for n in range(20)]
    curvature=radius*iv.mpf(local['real_geometry']['metric_norm_upper'])*iv.mpf(magnetic)/2
    values[0]=iv.mpf(min(upper(values[0]),upper(curvature)))
    updated=copy.deepcopy(prior)
    updated['jets']['total_order_derivative_upper']['gauge_per_abs_q']=[upper(x) for x in values]
    return updated,dict(coordinate_inner_radius=lower(radius),coordinate_outer_radius=lower(2*radius),
        gauge_derivatives_through=19,gauge_per_abs_q_derivative_upper=[upper(x) for x in values],
        curvature_zero_order_upper=upper(curvature),original_zero_order_upper=upper(old[0]),
        zero_order_ceiling_ratio_lower=lower(old[0]/values[0]))


def saved_normal(iv, local):
    record=local['normal_geometry']
    return dict(radius=iv.mpf(record['coordinate_velocity_radius_lower']),
        rho=iv.mpf(record['normal_metric_radius_lower']),
        Z=[iv.mpf(0)]+[unpack(iv,x) for x in record['scaled_phase_normalized_frechet_upper']],
        exp=[iv.mpf(0)]+[unpack(iv,x) for x in record['exp_normalized_frechet_upper']],
        inverse=[iv.mpf(0)]+[unpack(iv,x) for x in record['inverse_normalized_frechet_upper']])


def controls():
    """Independent polynomial pure-gauge and constant-curvature examples."""
    mp.dps=75
    source=tuple(map(mp.mpf,('.1','-.2','.3')))
    points=[tuple(map(mp.mpf,p)) for p in [('.4','.2','-.1'),('-.2','.3','.7')]]
    def chi(x):return 7*x[0]+11*x[1]+13*x[2]+17*x[0]*x[1]+19*x[2]**2
    def A(x,B):return (-B*x[1]/2+7+17*x[1],B*x[0]/2+11+17*x[0],13+38*x[2])
    def phase(x,B):
        h=[x[i]-source[i] for i in range(3)]
        return mp.quad(lambda t:sum(h[i]*A([source[j]+t*h[j] for j in range(3)],B)[i] for i in range(3)),[0,1])
    def change(x,i,z):return [z if j==i else x[j] for j in range(3)]
    def new_A(x,B):return tuple(A(x,B)[i]-mp.diff(lambda z:phase(change(x,i,z),B),x[i]) for i in range(3))
    errors=[];covariance=[];curvature=[]
    for B in (mp.mpf(0),mp.mpf(5)):
        for x in points:
            h=[x[i]-source[i] for i in range(3)]
            exact=(-B*h[1]/2,B*h[0]/2,mp.mpf(0))
            observed=new_A(x,B)
            errors.extend(abs(a-b) for a,b in zip(observed,exact))
            exact_phase=chi(x)-chi(source)-B*source[1]*h[0]/2+B*source[0]*h[1]/2
            errors.append(abs(phase(x,B)-exact_phase))
            curvature.append(abs(mp.diff(lambda z:new_A(change(x,0,z),B)[1],x[0])
                -mp.diff(lambda z:new_A(change(x,1,z),B)[0],x[1])-B))
            psi=lambda v:mp.exp(-sum(z*z for z in v)/3)
            transformed=lambda v:mp.exp(-1j*phase(v,B))*psi(v)
            for i in range(3):
                lhs=mp.diff(lambda z:transformed(change(x,i,z)),x[i])-1j*observed[i]*transformed(x)
                rhs=mp.exp(-1j*phase(x,B))*(mp.diff(lambda z:psi(change(x,i,z)),x[i])-1j*A(x,B)[i]*psi(x))
                covariance.append(abs(lhs-rhs))
    return dict(radial_gauge_identity_error=float(max(errors)),curvature_preservation_error=float(max(curvature)),
                covariant_derivative_identity_error=float(max(covariance)))


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_scaled_graph_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_scaled_graph.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    local=json.loads(project_path('dihole_local_transport_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    heat=json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    updated,gauge=centred_jets(iv,prior,local,magnetic['geometry_constants']['magnetic_unit_flux_upper_D2'])
    normal=saved_normal(iv,local)
    radius=iv.mpf(gauge['coordinate_inner_radius'])
    domain_checks=dict(normal_support_inside=bool((2*normal['radius']).b<radius.a),
        evaluation_cube_inside=bool((iv.mpf(2)**-35).b<radius.a),
        phase_extension_inside_buffer=bool((2*radius).b<iv.mpf(local['normal_geometry']['coordinate_clearance']).a))
    if not all(domain_checks.values()):raise ArithmeticError('Local gauge domain not covered.')
    records=[]
    for old,mr,hr in zip(previous['records'],magnetic['records'],heat['records']):
        if not old['q']==mr['q']==hr['q']:raise ArithmeticError('Flux records do not align.')
        operator=operator_bounds(iv,updated['jets'],old['q'])
        point=point_bound(iv,updated,operator,local['real_geometry'])['exported']
        transport=transport_bounds(iv,updated,operator,normal,local['real_geometry'])['exported']
        defect=assemble(iv,updated,local,operator,transport,point,hr,mr)
        scaled=record_bound(iv,updated,operator,defect,mr)
        ratio=iv.log(unpack(iv,old['heat_defect_B4_uniform_upper'])/unpack(iv,scaled['heat_defect_B4_uniform_upper']),10)
        records.append(dict(q=old['q'],operator=operator,point_evaluation=point,heat_transport=transport,
            heat_defect=defect,scaled_bound=scaled,
            gauge_ceiling_log10_change_interval=[lower(ratio),upper(ratio)],
            previous_B4_upper=old['heat_defect_B4_uniform_upper']))
    control=controls()
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),**domain_checks,
        gauge_order_explicit=gauge['gauge_derivatives_through']==19 and len(gauge['gauge_per_abs_q_derivative_upper'])==20,
        operator_orders_retained=all(len(v)==19 for r in records for v in r['operator']['coefficient_operator_norm_upper'].values()),
        transport_orders_retained=all(r['heat_transport']['derivative_orders']==[16,14,12,10,8] for r in records),
        zero_order_gauge_bound_improves=gauge['zero_order_ceiling_ratio_lower']>1,
        radial_gauge_control=control['radial_gauge_identity_error']<1e-60,
        curvature_control=control['curvature_preservation_error']<1e-60,
        gauge_covariance_control=control['covariant_derivative_identity_error']<1e-60,
        physical_source_not_claimed=all(r['scaled_bound']['physical_source'] is None and r['scaled_bound']['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,gauge=gauge,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,bound_encoding=previous['bound_encoding'],
        scope=['Source-specific ordinary-coordinate radial gauge on a fixed spatial ball, extended by a time-independent global phase held fixed under metric variation.',
               'All source derivatives act on the first kernel slot with the source fixed. No derivative of the source-dependent gauge is taken along the diagonal.',
               'The metric, curvature, magnetic flux, physical regulator, payload and Q4 time interval are unchanged.',
               'Electromagnetic transport along curved geodesics is retained. Only gauge derivatives through19 are claimed.',
               'These are conservative ceilings, not signed source values or physical stability and atom-transport results.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Centred-gauge checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,gauge=gauge,
        records=[dict(q=r['q'],new_B4_log10=r['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper'],
            ceiling_change=r['gauge_ceiling_log10_change_interval'],
            Cstar_log10=r['scaled_bound']['exact_low_order_remainder_Cstar_upper']['log10_upper']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-1]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['Previous local gauge','Source-centred gauge'],[row['previous_B4_upper']['log10_upper'],
        row['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper']],color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of B4 ceiling, D = 1',title='The same physical field in a better local gauge')
    fig.supxlabel('Atomic-flux reference. Conservative error ceilings remain too loose for support.\nNo stable opening, atom passage or clock gain has been established.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
