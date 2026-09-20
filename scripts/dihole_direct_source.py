"""Missing direct determinant response of the actual finite-angular bridge.

Accumulate the symmetric first-order secular determinant without inverting
an exponentially small transmission matrix. The power subtraction below is
a finite-angular pencil convention, NOT four-dimensional Hadamard stress.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import expm, solve_sylvester

from balanced_dihole import symbolic_fields
from curved_throat_modes import _constant_exact
from dihole_bridge_scattering import BridgeMatrices, collar_mask, embed
from dihole_stretched_collar import exact_target_parameters
from project_paths import project_path


def direct_response(coefficients, plus, minus, kappa, mask=collar_mask,
                    interval=(-.1,.1),rtol=2e-11,atol=2e-13):
    if not np.isfinite(kappa) or kappa<=0:
        raise ValueError('Use a positive finite imaginary frequency.')
    size=plus*minus
    initial=np.zeros(2*size+3,complex)

    def rhs(s,state):
        r=state[:size].reshape(plus,minus)
        z=state[size:2*size].reshape(plus,minus)
        wp,wm,x=coefficients(s); f=mask(s)
        a,b,c,d=-kappa*wp,-1j*x,1j*x.conj().T,kappa*wm
        ar=a-r@c; cd=c@r+d
        dr=ar@r+b-r@d
        dz=ar@z-z@cd-2*f*a@r-f*b+f*r@c@r+2*f*r@d
        trace_w=np.trace(wp)+np.trace(wm)
        direct=np.trace(c@r)+kappa*trace_w/2
        subtracted=np.trace(c@(z-f*r))
        volume=f*trace_w
        return np.r_[dr.ravel(),dz.ravel(),direct,subtracted,volume]

    boundaries=[interval[0]]+[s for s in (-.05,.05) if interval[0]<s<interval[1]]+[interval[1]]
    state=initial; evaluations=0
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        solution=solve_ivp(rhs,(left,right),state,method='DOP853',rtol=rtol,atol=atol)
        if not solution.success:raise RuntimeError(solution.message)
        state=solution.y[:,-1]; evaluations+=solution.nfev
    if not np.all(np.isfinite(state)):raise ArithmeticError('Nonfinite direct response.')
    g,sub,volume=state[-3:]
    result=dict(kappa=float(kappa),direct_logdet=float(g.real),
                direct_dlog=float((sub-kappa*volume).real),
                power_subtracted_dlog=float(sub.real),volume_factor=float(volume.real),
                scalar_imaginary_error=float(max(abs(state[-3:].imag))),nfev=evaluations)
    return result,state[:size].reshape(plus,minus)


def controls(model):
    scalar=[]; periodic=[]; chiral=[]
    for kappa in (.1,.7,3.):
        w,m,length=2.3,1.7,.4
        coefficients=lambda s:(np.array([[w]]),np.array([[w]]),np.array([[m]]))
        answer,r=direct_response(coefficients,1,1,kappa,lambda s:1.,(0.,length))
        expected=_constant_exact(1j*kappa*w,m,length)
        values=[]
        for e in (-1e-5,1e-5):
            t=_constant_exact(1j*kappa*w*np.exp(-2*e),m*np.exp(-e),length)[0,1]
            values.append(-np.log(abs(t)))
        scalar.append(dict(kappa=kappa,log_error=abs(answer['direct_logdet']+np.log(abs(expected[0,1]))),
                           derivative_error=abs(answer['direct_dlog']-(values[1]-values[0])/2e-5),
                           reflection_error=float(abs(r[0,0]-expected[1,1]))))
        k=np.array([[-kappa*w,-1j*m],[1j*m,kappa*w]])
        transfer=expm(k*length)
        for phase in (0.,.9,np.pi):
            env=np.array([[0,np.exp(-1j*phase)],[np.exp(1j*phase),0]])
            joined=answer['direct_logdet']+np.linalg.slogdet(np.eye(2)-expected@env)[1]
            exact=np.linalg.slogdet(transfer-np.exp(1j*phase)*np.eye(2))[1]
            periodic.append(float(abs(joined-exact)))
        one,_=direct_response(lambda s:(np.array([[w]]),np.empty((0,0)),np.empty((1,0))),
                              1,0,kappa,lambda s:1.,(0.,length))
        chiral.append(abs(one['direct_logdet']-kappa*w*length/2)+abs(one['power_subtracted_dlog']))
    answer,right=direct_response(model.coefficients,5,4,.1)
    old,_=embed(model.coefficients,5,4,.1j)
    reflection_error=float(np.max(abs(right-old[4:,5:])))
    varied=[direct_response(lambda s,e=e:model.coefficients(s,e),5,4,.1)[0]['direct_logdet'] for e in (-1e-4,1e-4)]
    derivative_error=abs(answer['direct_dlog']-(varied[1]-varied[0])/2e-4)
    uv_errors=[]
    for s in (-.08,0.,.08):
        base=[]
        for e in (0.,-.1,.1):
            wp,wm,x=model.coefficients(s,e)
            y=solve_sylvester(wp,wm,x)
            base.append(np.trace(x.conj().T@y).real)
        uv_errors.append(float(max(abs(np.array(base)-base[0]))))
    return dict(scalar=scalar,periodic_determinant_errors=periodic,
                unequal_port_chiral_errors=chiral,actual_reflection_error=reflection_error,
                actual_derivative_error=derivative_error,uv_mass_invariance_errors=uv_errors)


def integrals(rows,prior,stride=1):
    chosen=sorted(rows,key=lambda x:x['grid_index'])[::stride]
    old={r['grid_index']:r for r in prior['records']}
    k=np.array([r['kappa'] for r in chosen]); x=np.log(k)
    direct=np.array([r['direct_dlog'] for r in chosen])
    sub=np.array([r['power_subtracted_dlog'] for r in chosen])
    loop=np.array([[p['dlog'] for p in old[r['grid_index']]['domains']['16.0']['phases']] for r in chosen])
    integrate=lambda values:(-2/np.pi*simpson(k[:,None]*values,x=x,axis=0)).tolist()
    return dict(points=len(chosen),direct=integrate(direct[:,None])[0],
                direct_power_subtracted=integrate(sub[:,None])[0],
                full_unrenormalized=integrate(direct[:,None]+loop),
                full_power_subtracted=integrate(sub[:,None]+loop),
                loop_only=integrate(loop))


def render(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=sorted(result['records'],key=lambda x:x['kappa'])
    fig,axes=plt.subplots(1,2,figsize=(9,4.1),layout='constrained')
    for key,label in (('direct_dlog','Direct, before subtraction'),('power_subtracted_dlog','Direct, leading power removed')):
        axes[0].semilogx([r['kappa'] for r in rows],[-2/np.pi*r[key] for r in rows],label=label)
    fine=result['integrals']['fine']
    axes[1].bar(['loop only','direct + loop'],[fine['loop_only'][0],fine['full_power_subtracted'][0]],color=['tab:blue','tab:orange'])
    axes[0].set(xlabel='Imaginary frequency kappa D',ylabel='Determinant-response integrand')
    axes[1].set(ylabel='Finite-window response [1/D]',title='Phase 0, leading power removed')
    for ax in axes:ax.axhline(0,color='gray',lw=.7);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Finite-angular determinant convention; not covariant renormalized supporting stress',fontsize=10)
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_direct_source_progress.'+ext),dpi=180)
    plt.close(fig)


def run():
    start=perf_counter()
    names=('dihole_bridge_scattering.py','dihole_cylinder_transfer.py','dihole_stretched_collar.py',
           'balanced_dihole.py','spherical_weyl_leads.py','curved_throat_modes.py','project_paths.py',
           'dihole_source_window_results.json')
    hashes={str(project_path(n)):hashlib.sha256(project_path(n).read_bytes()).hexdigest() for n in names}
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    prior=json.loads(project_path('dihole_source_window_results.json').read_text())
    assert prior['complete'] and prior['passed']
    function,_=symbolic_fields(); model=BridgeMatrices(function,exact_target_parameters(),4,48)
    checked=controls(model)
    result=dict(complete=False,source_sha256=source_hash,input_sha256=hashes,controls=checked,records=[])
    output=project_path('dihole_direct_source_results.json')
    if output.exists():
        saved=json.loads(output.read_text())
        if saved['source_sha256']!=source_hash or saved['input_sha256']!=hashes:
            raise RuntimeError('Existing evidence uses different source/inputs; preserve before a changed run.')
        if saved.get('complete'):
            if not saved.get('passed'):raise RuntimeError('Saved exact run failed controls.')
            print('Exact source/input calculation already complete.');return saved
        result['records']=saved['records']
    done={r['grid_index'] for r in result['records']}
    for old in sorted(prior['records'],key=lambda r:r['grid_index']):
        if old['grid_index'] in done:continue
        row,_=direct_response(model.coefficients,5,4,old['kappa'])
        row['grid_index']=old['grid_index']; result['records'].append(row)
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(row),flush=True)
    result['integrals']=dict(coarse=integrals(result['records'],prior,2),fine=integrals(result['records'],prior))
    tests=dict(constant_scalar=max(max(r['log_error'],r['derivative_error'],r['reflection_error']) for r in checked['scalar'])<2e-8,
               exact_periodic_determinant=max(checked['periodic_determinant_errors'])<2e-10,
               unequal_port_centering=max(checked['unequal_port_chiral_errors'])<2e-11,
               actual_bridge_reflection=checked['actual_reflection_error']<2e-9,
               actual_variation=checked['actual_derivative_error']<2e-6,
               invariant_uv_mass_term=max(checked['uv_mass_invariance_errors'])<2e-10,
               real_accumulators=max(r['scalar_imaginary_error'] for r in result['records'])<2e-10,
               complete_grid={r['grid_index'] for r in result['records']}==set(range(33)),
               frozen_inputs=all(hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in hashes.items()))
    tests={name:bool(value) for name,value in tests.items()}
    result.update(complete=True,passed=all(tests.values()),tests=tests,runtime_seconds=perf_counter()-start,
                  sector=dict(q=1,mu=.5,angular_pairs=4),window=prior['window'],
                  convention='Symmetric first-order pencil logdet: -logdet Trl - .5 integral Tr K + logdet(I-B Rext). Conditional paired-Weyl prefactor -2/pi.',
                  scope=['The direct factor was absent from132. R16 loop responses are reused exactly; no exterior rerun.',
                         'Leading-power subtraction removes -kappa integral f Tr W from the direct metric derivative only.',
                         'This is a fixed-angular secular determinant convention; canonical measure/contact and covariant four-dimensional Hadamard/anomaly terms are not supplied.',
                         'No physical all-mode source, current, backreaction, stability or chronology change is proved.',
                         'The finite outer/near geometry, angular/frequency tails and quadrature changes remain uncontrolled errors.'])
    render(result)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,integrals=result['integrals']),indent=2))
    if not result['passed']:raise RuntimeError('Direct determinant controls failed; evidence retained.')
    return result


if __name__=='__main__':run()
