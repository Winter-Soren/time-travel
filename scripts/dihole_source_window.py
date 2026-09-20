"""Finite frequency-window interaction response of the actual curved pilot.

Reuse each exterior assembly while shifting only its physical mass matrix.
The displayed Dirac prefactor is conditional on the common-domain scattering
formula; the finite absorbing DG exterior does not establish that premise.
No omitted frequency/angular tail or absolute local source is supplied.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson
from scipy.linalg import block_diag
from scipy.sparse.linalg import splu

from balanced_dihole import symbolic_fields
from dihole_bridge_scattering import BridgeMatrices, complex_record, embed
from dihole_curved_exterior import curved_model, graded_mesh, attenuation_ledger
from dihole_return_pilot import cut_basis, cut_factors, sew
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


def prepare_exterior(radius, size=4):
    model, mass, geometry = curved_model(kappa=.1, triangulation=graded_mesh(Rtarget=radius))
    ports = [model.modal_port(kind, lambda points,kind=kind:cut_basis(kind,points,size))
             for kind in ('upper','lower')]
    factors = [cut_factors(p,size) for p in ports]
    return dict(model=model, mass=mass, geometry=geometry, ports=ports,
                zero_operator=(model.matrix_full-.1*mass).tocsc(),
                lift=np.column_stack([p['G'] for p in ports]),
                receiver=np.vstack([p['C'] for p in ports]),
                fi=block_diag(factors[0][1],factors[1][1]),
                fo=block_diag(factors[0][0],factors[1][0]))


def return_map(prepared, kappa, zero=None):
    model, mass = prepared['model'], prepared['mass']
    model.matrix_full=(prepared['zero_operator']+kappa*mass).tocsc()
    model.matrix=model.matrix_full[model.active][:,model.active].tocsc()
    model.lu=splu(model.matrix)
    model.kappa, model.energy=float(kappa),1j*kappa
    fields=model.solve(prepared['lift'])
    value=prepared['receiver']@fields
    residual=model.matrix@fields[model.active]-prepared['lift'][model.active]
    diagnostics=dict(maximum_singular=float(np.linalg.svd(value,compute_uv=False)[0]),
                     relative_linear_residual=float(np.linalg.norm(residual)/np.linalg.norm(prepared['lift'][model.active])))
    if zero is not None:
        size=len(zero)-1
        source=prepared['fi']@np.r_[zero,np.zeros(size)]
        source/=np.linalg.norm(source)
        incoming=model.modal_incoming(prepared['ports'],[source[:size+1],source[size+1:]])
        diagnostics['zero_input_ledger']=attenuation_ledger(
            model,mass,kappa,fields@source,incoming,prepared['lift']@source)
    return value,diagnostics


def connected_response(prepared, returned, bridge, minus, plus, epsilon=1e-4):
    rows=[]
    for phase in (0.,np.pi):
        mapped=[sew(s,phase,4,prepared['fi'],prepared['fo']) for s in (bridge,minus,plus)]
        b,bm,bp=mapped
        denominator=np.eye(len(b))-b@returned
        derivative=-float(np.real(np.trace(np.linalg.solve(denominator,((bp-bm)/(2*epsilon))@returned))))
        fd=float((np.linalg.slogdet(np.eye(len(b))-bp@returned)[1]
                  -np.linalg.slogdet(np.eye(len(b))-bm@returned)[1])/(2*epsilon))
        rows.append(dict(phase=phase,real_logdet=float(np.linalg.slogdet(denominator)[1]),
                         dlog=derivative,determinant_derivative_error=abs(derivative-fd),
                         denominator_condition=float(np.linalg.cond(denominator))))
    return rows


def integrate_records(records, radius, stride=1):
    rows=sorted(records,key=lambda r:r['grid_index'])[::stride]
    x=np.array([r['log_kappa'] for r in rows]); k=np.exp(x)
    values=np.array([[p['dlog'] for p in r['domains'][str(radius)]['phases']] for r in rows])
    answer=-2/np.pi*simpson(k[:,None]*values,x=x,axis=0)
    logs=np.array([[p['real_logdet'] for p in r['domains'][str(radius)]['phases']] for r in rows])
    energy=-2/np.pi*simpson(k[:,None]*logs,x=x,axis=0)
    return dict(points=len(rows),conditional_interaction_energy=energy.tolist(),
                conditional_metric_response=answer.tolist(),
                phase0_minus_phasepi_response=float(answer[0]-answer[1]))


def run():
    started=perf_counter()
    names=('dihole_return_pilot.py','dihole_curved_exterior.py','dihole_bridge_scattering.py',
           'dihole_cylinder_transfer.py','dihole_stretched_collar.py','balanced_dihole.py',
           'exterior_boundary_map.py','two_mouth_mesh.py','spherical_weyl_leads.py','project_paths.py')
    hashes={str(project_path(n)):hashlib.sha256(project_path(n).read_bytes()).hexdigest() for n in names}
    prior_path=project_path('dihole_return_pilot_results.json')
    hashes[str(prior_path)]=hashlib.sha256(prior_path.read_bytes()).hexdigest()
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output=project_path('dihole_source_window_results.json')
    result=dict(complete=False,source_sha256=source_hash,input_sha256=hashes,
                window=[.01,10.],sector=dict(q=1,mu=.5,angular_pairs=4),records=[])
    if output.exists():
        saved=json.loads(output.read_text())
        if saved.get('source_sha256')!=source_hash or saved.get('input_sha256')!=hashes:
            raise RuntimeError('Existing checkpoint has different source/inputs; preserve it before a changed run.')
        if saved.get('complete'):
            print('This exact source/input run is already complete; no recomputation.');return saved
        result=saved
    function,_=symbolic_fields()
    bridge_model=BridgeMatrices(function,exact_target_parameters(),4,48)
    zero=bridge_model.cut_zero()[0]
    systems={str(radius):prepare_exterior(radius) for radius in (8.,16.)}
    result['geometry']={key:p['geometry'] for key,p in systems.items()}
    # Check shifted assembly against frozen131 at a frequency outside this grid.
    prior=json.loads(prior_path.read_text())
    anchors=[]
    for key,system in systems.items():
        actual,diag=return_map(system,.1,zero)
        old=next(r for r in prior['records'] if r['graded'] and r['outer_radius']==float(key))['exterior_return']
        expected=np.array(old['real'])+1j*np.array(old['imag'])
        anchors.append(dict(radius=float(key),map_change=float(np.max(abs(actual-expected))),diagnostics=diag))
    result['assembly_anchors']=anchors
    grid=np.linspace(np.log(.01),np.log(10.),33)
    order=list(range(0,33,2))+list(range(1,33,2))
    done={r['grid_index'] for r in result['records']}
    for index in order:
        if index in done:continue
        kappa=float(np.exp(grid[index])); t=perf_counter()
        bridge,_=embed(bridge_model.coefficients,5,4,1j*kappa)
        varied=[embed(lambda s,e=e:bridge_model.coefficients(s,e),5,4,1j*kappa)[0] for e in (-1e-4,1e-4)]
        row=dict(grid_index=index,log_kappa=float(grid[index]),kappa=kappa,domains={})
        for key,system in systems.items():
            returned,diag=return_map(system,kappa,zero if index in (0,32) else None)
            phases=connected_response(system,returned,bridge,*varied)
            row['domains'][key]=dict(phases=phases,diagnostics=diag)
        row['runtime_seconds']=perf_counter()-t
        result['records'].append(row)
        result['completed_points']=len(result['records'])
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(dict(done=len(result['records']),kappa=kappa,
                              dlog_R16=[p['dlog'] for p in row['domains']['16.0']['phases']],
                              seconds=row['runtime_seconds'])),flush=True)
    result['integrals']={key:dict(coarse=integrate_records(result['records'],float(key),2),
                                 fine=integrate_records(result['records'],float(key))) for key in systems}
    domains=[d for r in result['records'] for d in r['domains'].values()]
    phases=[p for d in domains for p in d['phases']]
    ledgers=[d['diagnostics']['zero_input_ledger'] for d in domains if 'zero_input_ledger' in d['diagnostics']]
    tests=dict(shifted_assembly_matches131=max(a['map_change'] for a in anchors)<2e-11,
               passive_finite_return=max(d['diagnostics']['maximum_singular'] for d in domains)<1+1e-8,
               linear_solves=max(d['diagnostics']['relative_linear_residual'] for d in domains)<2e-10,
               extreme_frequency_flux=max(l['attenuation_identity_error'] for l in ledgers)<2e-9,
               finite_metric_response=all(np.isfinite(p['dlog']) for p in phases),
               finite_determinant_derivative=max(p['determinant_derivative_error'] for p in phases)<3e-6,
               complete_nested_grid={r['grid_index'] for r in result['records']}==set(range(33)),
               frozen_inputs=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()))
    result.update(complete=True,passed=all(tests.values()),tests=tests,runtime_seconds=perf_counter()-started,
                  normalization='Conditional Dirac interaction: -(2/pi) integral Re logdet dkappa; response uses its compact metric derivative. Units hbar=c=1, energy D^-1.',
                  scope=['Actual coupled finite bridge/exterior geometry, one mu sector, finite angular basis.',
                         'Finite window .01<=kappa D<=10; neither frequency tail is bounded.',
                         'Coarse/fine quadrature and R8/R16 changes are sensitivities, not certified errors.',
                         'Additional spatial/angular limits and finite-difference error remain; no forced sign or convergence test.',
                         'Incoming-zero outer boundary is not transparent and does not supply a physical sewn self-adjoint spectral vacuum.',
                         'Interaction contributions omit absolute local renormalized source terms; phase difference can cancel common local terms only in a valid common continuum state construction.',
                         'No supporting material, current, stability, formation or chronology-changing solution supplied.'])
    render(result)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,integrals=result['integrals']),indent=2))
    if not result['passed']:raise RuntimeError('Frequency-window algebraic controls failed; results retained.')
    return result


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=sorted(result['records'],key=lambda r:r['kappa'])
    fig,axes=plt.subplots(1,2,figsize=(9.2,4.1),layout='constrained')
    for radius,style in (('8.0','--'),('16.0','-')):
        for j,label in enumerate(('phase 0','phase pi')):
            axes[0].semilogx([r['kappa'] for r in rows],[-2/np.pi*r['domains'][radius]['phases'][j]['dlog'] for r in rows],
                             ls=style,color=('tab:blue','tab:orange')[j],label=f'{label}, R={radius}')
    axes[0].set(xlabel='Imaginary frequency kappa D',ylabel='Conditional interaction-response integrand')
    for j,label in enumerate(('phase 0','phase pi')):
        vals=[result['integrals'][r]['fine']['conditional_metric_response'][j] for r in ('8.0','16.0')]
        axes[1].plot([8,16],vals,'o-',label=label)
    axes[1].set(xlabel='Outer radius R / D',ylabel='Finite-window metric response [1/D]')
    for ax in axes:ax.axhline(0,color='gray',lw=.7);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Finite connected interaction contribution; absolute supporting stress remains unknown',fontsize=10)
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_source_window_progress.'+ext),dpi=180)
    plt.close(fig)


if __name__=='__main__':
    run()
