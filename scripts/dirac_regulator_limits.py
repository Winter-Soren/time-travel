"""Controls for the common source's frequency and regulator-limit proof.

The matrix is the existing noncommuting-lapse algebra benchmark. The momentum
integrals are flat symbol controls, not a stress calculation on the sewn bridge.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import expit

from dirac_green_mass_limit import benchmark
from project_paths import project_path


def high_weight(eigenvalues,threshold):
    if not np.isfinite(threshold) or threshold<=0:
        raise ValueError('A finite positive Hamiltonian energy threshold is required.')
    t=(abs(eigenvalues)-threshold)/threshold
    q=np.zeros_like(t)
    q[t>=1]=1
    mask=(t>0)&(t<1)
    q[mask]=expit(1/(1-t[mask])-1/t[mask])
    return q


def split_source(largest_lapse,kappa,eta,threshold=.5,check_split=True):
    """High-H split after variation, retaining the exact lapse contact term."""
    H,beta,N,V,dN,A,P,_=benchmark(largest_lapse,kappa)
    sqrtN=np.diag(np.sqrt(np.diag(N)))
    invsqrtN=np.diag(1/np.sqrt(np.diag(N)))
    f=np.linalg.solve(N,dN)
    dotA=beta@invsqrtN@V@invsqrtN-(f@A+A@f)/2
    high=np.zeros_like(A);low=np.zeros_like(A);reference=np.zeros_like(A)
    contacts=[];norm_ratios=[]
    eye=np.eye(len(H))
    for j,cj in enumerate((1,-3,3,-1)):
        mass=np.sqrt(eta+j)
        for sign in (-1,1):
            m=sign*mass
            Hm=H+m*N@beta
            eigenvalues,vectors=np.linalg.eigh(Hm)
            weights=high_weight(eigenvalues,threshold)
            inverse=np.zeros_like(eigenvalues,dtype=complex)
            np.divide(weights,eigenvalues+1j*kappa,out=inverse,where=weights!=0)
            rhigh=(vectors*inverse)@vectors.conj().T
            norm_ratios.append(float(np.linalg.norm(rhigh,2)*threshold))
            high+=cj/2*sqrtN@rhigh@sqrtN@beta
            if check_split and (eta>0 or kappa!=0):
                rlow=(vectors*((1-weights)/(eigenvalues+1j*kappa)))@vectors.conj().T
                low+=cj/2*sqrtN@rlow@sqrtN@beta
                band=(vectors*(1-weights))@vectors.conj().T
                lhs=np.trace(dotA@sqrtN@rlow@sqrtN@beta)
                rhs=np.trace(rlow@(V+m*dN@beta))-np.trace(band@f)
                contacts.append(float(abs(lhs-rhs)/max(1,abs(lhs),abs(rhs))))
        if check_split and (eta>0 or kappa!=0):
            reference-=cj*A@np.linalg.inv(P+(eta+j)*eye)
    source=dotA@high
    upper=len(H)*8*np.linalg.norm(N,2)*np.linalg.norm(dotA,2)/threshold
    return dict(source=source,norm=float(np.linalg.norm(source,'nuc')),
                norm_bound=float(upper),high_norm_ratio=max(norm_ratios),
                contact_error=max(contacts,default=0.),
                split_error=float(np.linalg.norm(high+low-reference)/max(1,np.linalg.norm(reference)))
                            if check_split and (eta>0 or kappa!=0) else None)


def symbol_frequency_bound(kappa,eta,mass_scale=1.):
    """One scalar multiplicity, flat d^3p/(2pi)^3, source symbol q^2 R_PV."""
    if kappa<=0 or eta<0 or mass_scale<=0:
        raise ValueError('Positive frequency and mass scale, nonnegative eta required.')
    shifts=eta+np.arange(4)*mass_scale**2
    def integrand(angle):
        cosine=np.cos(angle)**2
        return np.sin(angle)**2*cosine/np.prod(1+shifts*cosine/kappa**2)
    integral,error=quad(integrand,0,np.pi/2,epsabs=2e-14,epsrel=2e-13)
    value=3*mass_scale**6/(np.pi**2*kappa**3)*integral
    bound=3*mass_scale**6/(16*np.pi*kappa**3)
    return value,bound,error


def run():
    x,eta,h=sp.symbols('x eta h',positive=True)
    coefficients=(1,-3,3,-1)
    rational=sum(c/(x+eta+j*h) for j,c in enumerate(coefficients))
    product=6*h**3/sp.prod(x+eta+j*h for j in range(4))
    moments=[sp.expand(sum(c*(eta+j*h)**r for j,c in enumerate(coefficients))) for r in range(4)]
    matrices=[];corner=[];errors=[];ratios=[];bounds=[];convergence=[]
    for largest in (2.,10.,100.):
        for kappa in (0.,1e-6,.03,1.):
            baseline=split_source(largest,kappa,0.)
            # The high part is regular at eta=0; do not evaluate the singular
            # low inverse merely to test its continuity at smaller shifts.
            tiny=split_source(largest,kappa,1e-12,check_split=False)
            tiny_change=float(np.linalg.norm(tiny['source']-baseline['source'],'nuc'))
            corner.append(tiny_change)
            convergence.append(dict(largest_lapse=largest,kappa=kappa,eta=1e-12,
                                    distance_from_eta_zero=tiny_change))
            for shift in (1e-8,1e-5,.01):
                row=split_source(largest,kappa,shift)
                change=float(np.linalg.norm(row['source']-baseline['source'],'nuc'))
                matrices.append(dict(largest_lapse=largest,kappa=kappa,eta=shift,
                    high_source_trace_norm=row['norm'],finite_matrix_trace_bound=row['norm_bound'],
                    distance_from_eta_zero=change,contact_error=row['contact_error'],split_error=row['split_error']))
                errors.append(row['contact_error']);ratios.append(row['high_norm_ratio'])
                bounds.append(row['norm']/row['norm_bound'])
    rows=[];quadrature=[]
    for kappa in np.geomspace(1,100,31):
        values=[]
        for shift in (0.,.1,1.):
            value,bound,error=symbol_frequency_bound(kappa,shift)
            values.append(dict(eta=shift,value=value,analytic_bound=bound,ratio=value/bound))
            quadrature.append(error)
        rows.append(dict(kappa=float(kappa),values=values,
                         integrated_two_sided_tail_bound=3/(32*np.pi**2*kappa**2)))
    # Regulator-scale derivative heat moments; the t=0 cancellation is evaluated stably.
    exact=[3*np.log(3)-4*np.log(2),2*np.log(2)-np.log(3),1/3,11/18]
    heat=[]
    for n in range(4):
        def integrand(t):
            if t==0: return 1. if n==0 else 0.
            return t**(n-2)*np.exp(-t)*(-np.expm1(-t))**2
        value,error=quad(integrand,0,np.inf,epsabs=2e-13,epsrel=2e-13)
        heat.append(dict(n=n,value=value,exact=float(exact[n]),absolute_error=abs(value-exact[n]),quadrature_error=error))
    beta0=(9*np.log(3)-12*np.log(2))/2
    beta1=6*np.log(2)-3*np.log(3)
    checks=dict(
        exact_common_pv_resolvent=sp.cancel(rational-product)==0,
        first_three_mass_moments_vanish=all(m==0 for m in moments[:3]),
        first_surviving_mass_moment=moments[3]==-6*h**3,
        high_resolvent_uniform_bound=max(ratios)<=1+1e-12,
        retained_lapse_contact_identity=max(errors)<2e-8,
        high_low_split_after_variation=max(r['split_error'] for r in matrices)<2e-6,
        finite_matrix_trace_bound=max(bounds)<=1,
        finite_matrix_eta_zero_continuity=max(corner)<2e-5,
        flat_symbol_frequency_bound=all(v['ratio']<=1+1e-12 for r in rows for v in r['values']),
        flat_symbol_asymptotic_order=abs(rows[-1]['values'][0]['ratio']-1)<.001,
        flat_symbol_quadrature=max(quadrature)<1e-11,
        heat_moment_integrals=max(r['absolute_error'] for r in heat)<1e-11,
        declared_counterterm_derivatives=max(abs(3*exact[0]-2*beta0),abs(3*exact[1]-beta1),abs(3*exact[2]-1))<1e-13)
    checks={key:bool(value) for key,value in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,matrix_controls=matrices,
        small_shift_high_only_controls=convergence,
        frequency_symbol_controls=rows,heat_scale_moments=heat,
        maximum_low_contact_relative_error=max(errors),maximum_small_eta_matrix_change=max(corner),
        first_unsubtracted_scale_derivative_coefficient='11 deltaA6 / (96 pi^2 M^3)',
        leading_integrated_scale_tail='11 deltaA6 / (192 pi^2 M^2), when the next heat coefficient is retained',
        scope='Finite noncommuting-lapse algebra and flat principal-symbol/heat-moment controls. The actual geometry proof is separate; its error constants and supporting stress are not evaluated.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                       for name in ('dirac_green_mass_limit.py','project_paths.py')})
    project_path('dirac_regulator_limits_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']: raise ArithmeticError('Regulator controls failed; results retained.')
    plot(result)
    print(json.dumps({key:result[key] for key in ('passed','checks','maximum_low_contact_relative_error','maximum_small_eta_matrix_change')},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    rows=result['frequency_symbol_controls'];kappa=[r['kappa'] for r in rows]
    for j,shift in enumerate((0.,.1,1.)):
        axes[0].loglog(kappa,[r['values'][j]['value'] for r in rows],label=f'Flat symbol: eta = {shift:g}')
    axes[0].loglog(kappa,[r['values'][0]['analytic_bound'] for r in rows],'k--',label='Analytic frequency bound proportional to kappa^-3')
    axes[0].set(title='Large-frequency contribution is integrable',xlabel='Frequency kappa',ylabel='Flat source-symbol integral (one multiplicity)')
    for largest in (2.,10.,100.):
        selected=[r for r in result['matrix_controls'] if r['largest_lapse']==largest and r['kappa']==0]
        axes[1].loglog([r['eta'] for r in selected],[r['distance_from_eta_zero'] for r in selected],'o-',label=f'Matrix max lapse {largest:g}')
    axes[1].set(title='High-energy part remains regular at zero frequency',xlabel='Proper squared-mass shift eta',ylabel='Trace-norm difference from eta = 0')
    for ax in axes: ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Controls for the source-limit proof\nThese curves do not calculate the bridge stress',fontsize=12)
    fig.savefig(project_path('dirac_regulator_limits_progress.png'),dpi=160)
    fig.savefig(project_path('dirac_regulator_limits_progress.pdf'))
    plt.close(fig)


if __name__=='__main__': run()
