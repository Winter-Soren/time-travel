"""Sharpen the unchanged positive symmetric frame using real metric bounds."""
import json
import math

from mpmath import mp

from dihole_centered_gauge import centred_jets, saved_normal
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import assemble, unpack
from dihole_local_transport import packed, point_bound, transport_bounds
from dihole_regular_jets import operator_bounds
from dihole_scaled_graph import record_bound
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_positive_frame'


def frame_bounds(iv, fields, real):
    oldE,oldI,g,h=(list(map(iv.mpf,fields[k])) for k in ('coframe','inverse_coframe','metric','inverse_metric'))
    M,H,m=(iv.mpf(real[k]) for k in ('metric_norm_upper','inverse_metric_norm_upper','metric_lower'))
    gop=[M,4*iv.mpf(real['first_component_upper'])]+[4*x for x in g[2:]]
    hop=[H];E=[iv.sqrt(M)];I=[iv.sqrt(H)]
    for n in range(1,len(g)):
        inverse_metric=H*sum(math.comb(n,j)*gop[j]*hop[n-j] for j in range(1,n+1))
        hop.append(iv.mpf(min(upper(4*h[n]),upper(inverse_metric))))
        direct=(gop[n]+sum(math.comb(n,j)*E[j]*E[n-j] for j in range(1,n)))/(2*iv.sqrt(m))
        E.append(iv.mpf(min(upper(4*oldE[n]),upper(direct))))
        inverse_root=iv.sqrt(M)/2*(hop[n]+sum(math.comb(n,j)*I[j]*I[n-j] for j in range(1,n)))
        inverse_product=I[0]*sum(math.comb(n,j)*E[j]*I[n-j] for j in range(1,n+1))
        I.append(iv.mpf(min(upper(4*oldI[n]),upper(inverse_root),upper(inverse_product))))
    # Internal entries above are operator norms. Only these exported minima
    # become component bounds; feeding them back without factor4 is unsafe.
    exportedE=[min(upper(a),upper(b)) for a,b in zip(oldE,E)]
    exportedI=[min(upper(a),upper(b)) for a,b in zip(oldI,I)]
    return dict(coframe=exportedE,inverse_coframe=exportedI,
        internal_metric_operator_upper=[packed(iv,x) for x in gop],
        internal_inverse_metric_operator_upper=[packed(iv,x) for x in hop],
        internal_frame_operator_upper=[packed(iv,x) for x in E],
        internal_inverse_frame_operator_upper=[packed(iv,x) for x in I])


def controls(iv):
    mp.dps=80
    E0=mp.diag([2,3,4,5])
    E1=mp.matrix([[0,mp.mpf('.1'),0,0],[mp.mpf('.1'),0,mp.mpf('.2'),0],[0,mp.mpf('.2'),0,mp.mpf('.1')],[0,0,mp.mpf('.1'),0]])
    E2=mp.matrix([[mp.mpf('.01')*(i+j+1) for j in range(4)] for i in range(4)])
    frame=lambda x:E0+x*E1+x*x*E2
    metric=lambda x:frame(x)*frame(x)
    def derivatives(fun):
        return [mp.matrix([[mp.diff(lambda x:fun(x)[i,j],mp.mpf(0),n) for j in range(4)] for i in range(4)]) for n in range(7)]
    g,h=derivatives(metric),derivatives(lambda x:metric(x)**-1)
    actualE,actualI=derivatives(frame),derivatives(lambda x:frame(x)**-1)
    component=lambda a:float(max(abs(x) for x in a))*(1+1e-14)
    fields=dict(metric=[component(a) for a in g],inverse_metric=[component(a) for a in h],
        coframe=[10*math.factorial(n)*100**n for n in range(7)],
        inverse_coframe=[10*math.factorial(n)*100**n for n in range(7)])
    real=dict(metric_norm_upper=36,inverse_metric_norm_upper=1,metric_lower=1,first_component_upper=component(g[1]))
    bounds=frame_bounds(iv,fields,real)
    component_ratios=[];operator_ratios=[]
    for actual,key,opkey in ((actualE,'coframe','internal_frame_operator_upper'),(actualI,'inverse_coframe','internal_inverse_frame_operator_upper')):
        for n,a in enumerate(actual):
            component_ratios.append(float(max(abs(x) for x in a)/bounds[key][n]))
            eigenvalues=mp.eigsy(a,eigvals_only=True)
            operator_ratios.append(float(max(abs(x) for x in eigenvalues)/mp.mpf(unpack(iv,bounds[opkey][n])._mpi_[0])))
    # Independent first-derivative Sylvester solve in the diagonal E0 basis.
    recovered=mp.matrix([[g[1][i,j]/(E0[i,i]+E0[j,j]) for j in range(4)] for i in range(4)])
    return dict(noncommuting_control=float(mp.norm(E0*E1-E1*E0)),
        sylvester_first_derivative_error=float(mp.norm(recovered-E1)),
        frame_component_max_ratio=max(component_ratios),frame_operator_max_ratio=max(operator_ratios))


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_centered_gauge_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_centered_gauge.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    local=json.loads(project_path('dihole_local_transport_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    heat=json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    updated,gauge=centred_jets(iv,prior,local,magnetic['geometry_constants']['magnetic_unit_flux_upper_D2'])
    frames=frame_bounds(iv,prior['jets']['total_order_derivative_upper'],local['real_geometry'])
    for key in ('coframe','inverse_coframe'):
        updated['jets']['total_order_derivative_upper'][key]=frames[key]
    normal=saved_normal(iv,local)
    records=[]
    for old,mr,hr in zip(previous['records'],magnetic['records'],heat['records']):
        if not old['q']==mr['q']==hr['q']:raise ArithmeticError('Flux records do not align.')
        operator=operator_bounds(iv,updated['jets'],old['q'])
        point=point_bound(iv,updated,operator,local['real_geometry'])['exported']
        transport=transport_bounds(iv,updated,operator,normal,local['real_geometry'])['exported']
        defect=assemble(iv,updated,local,operator,transport,point,hr,mr)
        scaled=record_bound(iv,updated,operator,defect,mr)
        baseline=old['scaled_bound']['heat_defect_B4_uniform_upper']
        ratio=iv.log(unpack(iv,baseline)/unpack(iv,scaled['heat_defect_B4_uniform_upper']),10)
        records.append(dict(q=old['q'],operator=operator,point_evaluation=point,heat_transport=transport,
            heat_defect=defect,scaled_bound=scaled,previous_B4_upper=baseline,
            frame_ceiling_log10_change_interval=[lower(ratio),upper(ratio)]))
    control=controls(iv)
    original=prior['jets']['total_order_derivative_upper']
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        full_frame_orders=all(len(frames[k])==21 for k in ('coframe','inverse_coframe')),
        component_bounds_do_not_increase=all(a<=b for k in ('coframe','inverse_coframe') for a,b in zip(frames[k],original[k])),
        zero_frame_bounds_improve=all(frames[k][0]<original[k][0] for k in ('coframe','inverse_coframe')),
        noncommuting_control=control['noncommuting_control']>0,
        sylvester_control=control['sylvester_first_derivative_error']<1e-65,
        component_control=control['frame_component_max_ratio']<1,
        operator_control=control['frame_operator_max_ratio']<1,
        gauge_unchanged=gauge==previous['gauge'],
        derivative_schedule=all(r['heat_transport']['derivative_orders']==[16,14,12,10,8] for r in records),
        finite_positive_bounds=all(int(r['scaled_bound']['heat_defect_B4_uniform_upper']['mantissa'])>0 for r in records),
        physical_source_not_claimed=all(r['scaled_bound']['physical_source'] is None and r['scaled_bound']['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,frame_bounds=frames,gauge=gauge,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,bound_encoding=previous['bound_encoding'],
        scope=['The same positive symmetric coframe is bounded through order20 using Sylvester equations and inverse identities; it is not rotated or changed.',
               'Internal recurrences use matrix operator norms. Component-only ceilings are converted with factor4 before internal use.',
               'The source-centred gauge, physical geometry, field, regulator and atomic target remain unchanged.',
               'All downstream operator, transport and full-defect bounds are recomputed, with the original normal-neighbourhood estimates retained.',
               'Conservative ceilings do not establish supporting stress, stability, atom passage or clock gain.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Positive-frame checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,
        frame_first={k:frames[k][:2] for k in ('coframe','inverse_coframe')},
        records=[dict(q=r['q'],B4_log10=r['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper'],
            improvement=r['frame_ceiling_log10_change_interval'],Cstar_log10=r['scaled_bound']['exact_low_order_remainder_Cstar_upper']['log10_upper']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-1]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['Previous frame bounds','Real-metric frame bounds'],[row['previous_B4_upper']['log10_upper'],
        row['scaled_bound']['heat_defect_B4_uniform_upper']['log10_upper']],color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of B4 ceiling, D = 1',title='Sharper geometry bounds for the same atomic target')
    fig.supxlabel('Conservative error ceilings. Physical support remains unproved.\nNo field, opening or payload parameter was changed.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
