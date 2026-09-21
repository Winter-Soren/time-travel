"""Paired-frequency tails and finite selected massless source bounds.

The infrared smoothing is spectral, not a physical Pauli-Villars mass.
No omitted sectors, regulator matching, or full supporting stress is supplied.
"""
import json
import math

import numpy as np
from scipy.integrate import quad

from dihole_channel_bounds import MPIntervalContext, upper
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_paired_frequency'


def frequency_integrals(iv,delta,omega):
    if delta<=0 or omega<=0:raise ValueError('Positive spectral smoothing and split required.')
    d,k=iv.mpf(delta),iv.mpf(omega);root=iv.sqrt(k*k+d*d)
    angle=iv.atan2(k,d)
    return (k-d*angle,iv.ln((k+root)/d)-k/root,
        angle/(2*d)-k/(2*(k*k+d*d)),
        iv.atan2(d,k)/(2*d)+k/(2*(k*k+d*d)))


def bounds(iv,base,metric,constants):
    q=base['q'];w=iv.mpf(constants['W_lower'][0]);W=iv.mpf(constants['W_upper'][1])
    lam_value=base['frequency_rows'][1]['filter_Lambda_D'];lam=iv.mpf(lam_value)
    a=iv.mpf(base['graph_h_coefficient_upper']);b=iv.mpf(base['graph_zero_coefficient_upper_D'])
    alpha=iv.mpf(metric['relative_velocity_derivative_upper']);beta=iv.mpf(metric['reduced_graph_zero_upper_D'])
    M=iv.mpf(metric['profile_metric_norm_upper']);ad=iv.mpf(metric['defect_graph_h_upper']);bd=iv.mpf(metric['defect_graph_zero_upper_D'])
    optical=iv.mpf(base['inner_optical_volume_upper_over_D'])
    A2=iv.mpf(upper(lam*optical/(2*iv.sqrt(iv.pi))))
    B2=iv.mpf(upper(lam**3*optical/(4*iv.sqrt(iv.pi))+lam*25600/(2*iv.sqrt(iv.pi)*w)))
    A,B=iv.sqrt(A2),iv.sqrt(B2)
    aH=4*(1+iv.sqrt(W/w));bH=(1200+4*iv.sqrt(400+31*iv.mpf(q)))/w
    BH=(1+a)*B+b*A;DH=aH*BH+bH*A;Dh=alpha*B+beta*A
    C=BH*DH+B*Dh
    K1=(M+aH+alpha)*a+ad
    K2=(M+aH+alpha)*b+(bH+beta)*a+bd
    K3=(bH+beta)*b
    cases=[]
    for delta in (1.,.1*lam_value):
        splits=[]
        for factor in (1.,4.,16.,64.):
            omega=factor*lam_value;I1,I2,I3,tail=frequency_integrals(iv,delta,omega)
            middle=2*A2/iv.pi*(K1*I1+K2*I2+K3*I3)
            high=4*C/iv.pi*tail
            splits.append(dict(split_omega_D=omega,finite_window_error_upper_D=upper(middle),
                paired_tail_error_upper_D=upper(high),combined_error_upper_D=upper(middle+high)))
        whole=upper(C/iv.mpf(delta))
        best=min([whole]+[x['combined_error_upper_D'] for x in splits])
        cases.append(dict(spectral_delta_D=delta,splits=splits,all_frequency_pair_bound_D=whole,
            best_selected_energy_derivative_error_upper_D=best))
    return dict(q=q,band=base['band'],filter_Lambda_D=lam_value,
        filter_HS_squared_upper=upper(A2),h_filter_HS_squared_upper_D2=upper(B2),
        full_graph_h_coefficient_upper=upper(aH),full_graph_zero_upper_D=upper(bH),
        HY_HS_upper_D=upper(BH),dotH_Y_HS_upper_D=upper(DH),doth_X_HS_upper_D=upper(Dh),
        paired_energy_tail_constant_upper_D3=upper(C),
        finite_massless_symmetrized_insertion_upper_D=upper(A*DH),
        finite_reduced_symmetrized_insertion_upper_D=upper(A*Dh),
        band_massless_symmetrized_insertion_upper_D=upper(A*DH*base['band']['count']),
        cases=cases)


def integral_controls():
    errors=[];iv=MPIntervalContext();iv.dps=60
    for d,k in ((.2,.1),(1.,1.),(.1,10.),(2.,100.)):
        exact=frequency_integrals(iv,d,k)
        numerical=[quad(lambda x:x*x/(x*x+d*d)**power,0,k,epsabs=1e-11)[0] for power in (1.,1.5,2.)]
        numerical.append(quad(lambda x:x*x/(x*x+d*d)**2,k,np.inf,epsabs=1e-11)[0])
        errors.extend(abs(value-float(exact[j].mid))/max(1,abs(value)) for j,value in enumerate(numerical))
    moments=[]
    for lam in (.3,2.,8.):
        spectral=quad(lambda p:p*p*math.exp(-p*p/lam**2)/(2*math.pi),-np.inf,np.inf,epsabs=1e-11)[0]
        spatial=quad(lambda x:lam**6*x*x*math.exp(-lam*lam*x*x)/(2*math.pi),-np.inf,np.inf,epsabs=1e-11)[0]
        exact=lam**3/(4*math.sqrt(math.pi))
        moments.extend((abs(spectral/exact-1),abs(spatial/exact-1)))
    rejected=[]
    for delta,omega in ((0,1),(-1,1),(1,0)):
        try:frequency_integrals(iv,delta,omega)
        except ValueError:rejected.append(True)
        else:rejected.append(False)
    return dict(frequency_integral_relative_errors=errors,gaussian_derivative_moment_errors=moments,invalid_parameters_rejected=rejected)


def matrix_controls():
    # Opposite-sign eigenvalues and an off-diagonal variation expose ordering errors.
    H=np.diag([1.,-2.]);V=np.array([[0.,1.],[1.,0.]]);Y=np.ones((2,1))/np.sqrt(2)
    U=np.linalg.norm(Y,'fro');B=np.linalg.norm(H@Y,'fro');D=np.linalg.norm(V,2)*U
    def insertion(k):
        G=np.linalg.inv(H-1j*k*np.eye(2))
        return np.trace(Y.conj().T@G@V@G@Y)
    paired=[]
    for y in (.2,1.,8.):
        S=np.linalg.inv(H@H+y*y*np.eye(2))
        exact=4j*y*np.real(np.trace((H@S@Y).conj().T@V@S@Y))
        pair=insertion(y)-np.conj(insertion(y))
        paired.append(dict(y=y,identity_error=float(abs(pair-exact)),imaginary_bound_holds=bool(abs(insertion(y).imag)<=2*B*D/y**3)))
    energy=[]
    for delta in (.2,1.,2.):
        integrand=lambda omega:2/math.pi*omega**2/math.hypot(omega,delta)*insertion(math.hypot(omega,delta)).imag
        integrated=quad(integrand,0,np.inf,epsabs=1e-11)[0]
        exact=(math.sqrt(1+delta*delta)-math.sqrt(4+delta*delta))/3
        tails=[]
        for omega in (1.,4.,16.):
            actual=quad(lambda x:abs(integrand(x)),omega,np.inf,epsabs=1e-11)[0]
            limit=4*B*D/math.pi*(math.atan(delta/omega)/(2*delta)+omega/(2*(omega*omega+delta*delta)))
            tails.append(dict(omega=omega,absolute_tail=actual,bound=limit,enclosed=bool(actual<=limit)))
        energy.append(dict(delta=delta,integrated_derivative=integrated,frechet_derivative=exact,error=abs(integrated-exact),tails=tails))
    sign=np.diag([1.,-1.]);selected_sign=float(np.real(np.trace(Y.conj().T@V@sign@Y)))
    unsmoothed=quad(lambda y:2/math.pi*y*insertion(y).imag,0,np.inf,epsabs=1e-11)[0]
    sign_derivative=quad(lambda y:-2/math.pi*insertion(y).real,0,np.inf,epsabs=1e-11)[0]
    # Physical mass has a lapse-dependent commutator, absent in scalar smoothing.
    hw=np.array([[2.,1.],[1.,-1.]]);zero=np.zeros((2,2));eye=np.eye(2)
    H0=np.block([[hw,zero],[zero,-hw]]);beta=np.block([[zero,eye],[eye,zero]])
    N=np.diag([1.,3.,1.,3.]);m=.7;Hm=H0+m*N@beta
    proper=H0@H0+m*(H0@N-N@H0)@beta+m*m*N@N
    return dict(paired=paired,energy=energy,selected_symmetrized_sign_insertion=selected_sign,
        unsmoothed_selected_energy_derivative=unsmoothed,selected_sign_derivative=sign_derivative,
        selected_sign_bound_holds=bool(abs(selected_sign)<=U*D),
        proper_mass_square_error=float(np.linalg.norm(Hm@Hm-proper)),
        scalar_mass_shortcut_mismatch=float(np.linalg.norm(Hm@Hm-(H0@H0+m*m*np.eye(4)))))


def run():
    path=project_path('dihole_metric_resolvent_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256']);inputs.update({path.name:sha(path),'dihole_metric_resolvent.py':previous['source_sha256']})
    collar=json.loads(project_path('dihole_collar_resolvent_results.json').read_text())
    constants=json.loads(project_path('dihole_channel_bounds_results.json').read_text())['uniform_bound']['constants']
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing input provenance required.')
    iv=MPIntervalContext();iv.dps=60
    rows=[bounds(iv,b,m,constants) for b,m in zip(collar['records'],previous['records'],strict=True)]
    integrals=integral_controls();matrix=matrix_controls()
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        frequency_closed_integrals=max(integrals['frequency_integral_relative_errors'])<1e-10,
        gaussian_first_graph_moments=max(integrals['gaussian_derivative_moment_errors'])<1e-10,
        rejects_zero_smoothing=all(integrals['invalid_parameters_rejected']),
        exact_paired_identity=max(x['identity_error'] for x in matrix['paired'])<1e-13,
        imaginary_resolvent_bound=all(x['imaginary_bound_holds'] for x in matrix['paired']),
        weighted_energy_frechet_identity=max(x['error'] for x in matrix['energy'])<1e-11,
        all_matrix_tails_enclosed=all(y['enclosed'] for x in matrix['energy'] for y in x['tails']),
        energy_is_not_filtered_sign_source=abs(matrix['unsmoothed_selected_energy_derivative']+1/3)<1e-12 and abs(matrix['selected_symmetrized_sign_insertion'])<1e-12,
        sign_derivative_is_different=abs(matrix['selected_sign_derivative']-2/3)<1e-12,
        selected_sign_bound=matrix['selected_sign_bound_holds'],
        proper_mass_identity=matrix['proper_mass_square_error']<1e-13 and matrix['scalar_mass_shortcut_mismatch']>1,
        finite_frequency_integrated_bounds=all(math.isfinite(c['best_selected_energy_derivative_error_upper_D']) and c['best_selected_energy_derivative_error_upper_D']>0 for x in rows for c in x['cases']),
        paired_tails_decrease=all(all(a['paired_tail_error_upper_D']>b['paired_tail_error_upper_D'] for a,b in zip(c['splits'],c['splits'][1:])) for x in rows for c in x['cases']))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,records=rows,integral_controls=integrals,matrix_controls=matrix,
        source_sha256=sha(__file__),input_sha256=inputs,total_renormalized_source=None,source_sign=None,atom_transit=None,
        scope=['One Weyl field, one collar selection, one declared axial sector; band values multiply only the exact saved sector count.',
            'Frozen Y=JX and X in the selected energy derivative; no moving embedding or filter energy terms are silently included.',
            'Finite delta>0 spectral sqrt(H^2+delta^2) smoothing only; no delta-removal or physical proper-mass identification.',
            'The symmetrized selected sign insertion uses sign(0)=0 and has its own direct bound, without spectral smoothing. It is not the selected energy derivative for a noncommuting spatial filter.',
            'No fermionic -1/2, Weyl doubling, second collar, physical hbar*c/D, PV signed masses or counterterms are folded into the exported finite-selection numbers.',
            'Proper PV uses H_m=H_0+m N beta, its full metric derivative, common signed-mass sum and fixed counterterms. The remaining modes and exterior are unbounded by this calculation.',
            'Atomic entries evaluate bounding formulas only; no supported opening, hydrogen evolution, clock advantage or new physical law.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Paired-frequency checks failed; evidence retained.')
    plot(result)
    print(json.dumps({'passed':result['passed'],'tests':tests,'atomic_bound':rows[-1]},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    atomic=result['records'][-1];lam=atomic['filter_Lambda_D'];case=atomic['cases'][0]
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    axes[0].loglog([x['split_omega_D']/lam for x in case['splits']],
        [x['paired_tail_error_upper_D'] for x in case['splits']],'o-',color='#187c80')
    axes[0].set(title='Atomic reference: paired frequency tail',xlabel='Starting frequency / Gaussian scale',ylabel='Per-sector energy-response tail bound times D')
    qs=[x['q'] for x in result['records']]
    axes[1].loglog(qs,[x['finite_massless_symmetrized_insertion_upper_D']/x['q'] for x in result['records']],'o-',color='#7840aa')
    axes[1].set(title='Separate finite massless source selection',xlabel='Integer flux q (bounding formulas)',ylabel='Absolute insertion bound times D / q')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Frequency pairing gives a finite high-frequency remainder')
    fig.supxlabel('Fixed spatial filter. Left: spectral delta D = 1. Right: no IR smoothing. Neither is the full renormalized stress.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
