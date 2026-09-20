"""Positivity-compatible Schur pilot for the Dirichlet continuum spectrum.

The analytic comparison includes the infinite spatial complement. Numerical
Gram integration and matrix operations are not enclosed, so the displayed
eigenvalue/low-heat intervals remain conditional numerical diagnostics.
"""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, solve
from mpmath.ctx_iv import MPIntervalContext

from dihole_continuum_heat_tail import target
from dihole_spatial_complement import evaluators, coefficients, basis_grid, square_action
from dihole_quadrature_remainder import dyadic, outward
from dihole_heat_angular import gershgorin_lower
import dihole_interval_arithmetic as ia
from project_paths import project_path


def hermitian(matrix):
    return (matrix+matrix.conj().T)/2


def parallel_sum(C,B):
    """C:B; C must be positive. No commutation or inverse of B is assumed."""
    return hermitian(C-C@solve(C+B,C,assume_a='pos',check_finite=False))


def low_values(matrix,count=8):
    return eigh(hermitian(matrix),subset_by_index=(0,count-1),check_finite=False)[0]


def schur_lower(M,G,D,gamma,spectral):
    if not spectral<gamma:raise ValueError('Comparison denominator is not positive.')
    C=M-np.diag(D)
    B=G/(gamma-spectral)
    return hermitian(M-spectral*np.eye(len(M))-parallel_sum(C,B))


def controls():
    rng=np.random.default_rng(126);worst=0.;cases=0
    for _ in range(6):
        p,q=4,7
        random=lambda a,b:rng.normal(size=(a,b))+1j*rng.normal(size=(a,b))
        W=random(p,p);C=W@W.conj().T+2*np.eye(p)
        R=random(q,p);W=random(q,q);extra=W@W.conj().T
        Dp=np.diag(np.linspace(.5,2,p));Dq=np.diag(np.linspace(4,8,q))
        M=Dp+C;AQ=Dq+R@solve(C,R.conj().T,assume_a='pos')+extra
        A=np.block([[M,R.conj().T],[R,AQ]])
        G=R.conj().T@R
        for spectral in (-1.,0.,1.,2.,3.):
            T=R.conj().T@solve(Dq-spectral*np.eye(q),R,assume_a='pos')
            exact=R.conj().T@solve(AQ-spectral*np.eye(q),R,assume_a='pos')
            woodbury=R.conj().T@solve(Dq-spectral*np.eye(q)+R@solve(C,R.conj().T,assume_a='pos'),R,assume_a='pos')
            U=parallel_sum(C,T)
            error=float(np.max(abs(U-woodbury)));worst=max(worst,error)
            B=G/(4-spectral);upper=parallel_sum(C,B)
            assert min(low_values(U-exact,1))>-2e-12
            assert min(low_values(T-U,1))>-2e-12
            assert min(low_values(upper-U,1))>-2e-12
            actual_count=int(np.count_nonzero(eigh(A,check_finite=False)[0]<spectral))
            retained_count=int(np.count_nonzero(eigh(M,check_finite=False)[0]<spectral))
            lower_count=int(np.count_nonzero(eigh(M-spectral*np.eye(p)-upper,check_finite=False)[0]<0))
            assert retained_count<=actual_count<=lower_count
            cases+=1
    assert worst<2e-12
    return dict(passed=True,noncommuting_complex_cases=cases,woodbury_max_error=worst)


def comparison(iv):
    path=project_path('dihole_spatial_complement_results.json')
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='4f3694b383a435607b88056cf131f30d010f7e096ac4cb886f0881a790781c34'
    old=json.loads(path.read_text());row=next(g for g in old['gaps'] if g['count']==18)
    cover=old['covers'][-1];eps=dyadic(iv,row['epsilon'])
    kappa=(1-eps)*dyadic(iv,cover['w_squared_lower'])
    V=-dyadic(iv,old['endomorphism_lower'])+(1/eps-1)*dyadic(iv,cover['w_squared_upper'])*dyadic(iv,cover['connection_squared_upper'])
    # One common global comparison D=kappa*L-V. All values record outward intervals.
    diagonal=[outward(kappa*iv.pi**2*(i*i+j*j)/64-V) for i in range(1,19) for j in range(1,19)]
    gamma=outward(kappa*iv.pi**2*362/64-V)
    D=np.repeat([(x['lower']+x['upper'])/2 for x in diagonal],2)
    return D,gamma['lower'],dict(kappa=outward(kappa),V=outward(V),gamma=gamma,diagonal=diagonal,
        epsilon=row['epsilon'],scope='Enclosed global comparison constants. The numerical pilot uses midpoint retained diagonals; it is not a certificate.')


def residual_gram(function,constants,M,offset,order):
    grid,basis=basis_grid(18,order);u,v,weight=grid
    action=square_action(*basis,coefficients(function,constants,u,v,offset))
    residual=action-basis[0]@M
    weighted=np.repeat(np.sqrt(weight),2)[:,None]*residual.reshape(-1,len(M))
    # Direct positive Gram; do not subtract two large Gram matrices or clip eigenvalues.
    return hermitian(weighted.T@weighted)


def spectrum_pilot(M,G,D,gamma):
    C=M-np.diag(D);cmin=float(low_values(C,1)[0]);assert cmin>0
    upper=low_values(M);assert upper[-1]<gamma
    zero=low_values(schur_lower(M,G,D,gamma,0.));assert zero[0]>0
    brackets=[]
    for index,ritz in enumerate(upper):
        lo,hi=0.,float(ritz)
        # Positive kth lower-Schur eigenvalue implies kth full eigenvalue >= trial.
        for _ in range(20):
            mid=(lo+hi)/2
            value=float(eigh(schur_lower(M,G,D,gamma,mid),subset_by_index=(index,index),check_finite=False)[0][0])
            if value>0:lo=mid
            else:hi=mid
        margin=float(eigh(schur_lower(M,G,D,gamma,lo),subset_by_index=(index,index),check_finite=False)[0][0])
        brackets.append(dict(index=index+1,conditional_lower=lo,ritz_upper=float(ritz),search_upper=hi,lower_schur_margin=margin))
    shifts=[]
    for z in (1.,4.):
        L=schur_lower(M,G,D,gamma,-z)
        assert low_values(L,1)[0]>0
        # Response for each normalized retained Ritz state; finite-dimensional brackets.
        vectors=eigh(M,subset_by_index=(0,7),check_finite=False)[1]
        high=np.diag(vectors.T@solve(L,vectors,assume_a='pos',check_finite=False))
        low=1/(upper+z)
        assert np.min(high-low)>-1e-10
        shifts.append(dict(positive_shift=z,conditional_lower=low.tolist(),conditional_upper=high.tolist()))
    heat=[]
    lower=np.array([b['conditional_lower'] for b in brackets])
    for tau in (.5,1.,2.):
        lo=float(np.exp(-tau*upper).sum());hi=float(np.exp(-tau*lower).sum())
        heat.append(dict(heat_parameter=tau,conditional_low_eight_trace_lower=lo,
            conditional_low_eight_trace_upper=hi,relative_width=(hi-lo)/lo))
    raw=low_values(M-G/gamma,1)[0]
    return dict(compatibility_minimum=cmin,raw_schur_minimum_at_zero=float(raw),
        strengthened_schur_minima_at_zero=zero.tolist(),eigenvalue_brackets=brackets,
        retained_resolvent_brackets=shifts,low_eight_heat_trace=heat)


def positive_certificate(interval,label,saved):
    midpoint=(interval[0]+interval[1])/2
    _,Q=eigh(hermitian(midpoint),check_finite=False)
    gram=ia.matmul(ia.point(Q.T),ia.point(Q))
    gram_lower=gershgorin_lower(gram)
    transformed=ia.matmul(ia.point(Q.T),ia.matmul(interval,ia.point(Q)))
    lower=gershgorin_lower(transformed)
    assert gram_lower>0 and lower>0,(label,gram_lower,lower)
    saved[label+'_basis']=Q
    print(json.dumps(dict(stage='conditional_matrix_certificate',part=label,gram_lower=gram_lower,congruence_lower=lower)),flush=True)
    return dict(label=label,gram_gershgorin_lower=gram_lower,congruence_gershgorin_lower=lower)


def conditional_certificate(M,G,comp,saved):
    """Certify only the matrix implication under explicit, still-unproved input budgets."""
    em,eg,lam=1e-4,.01,.18
    eye=np.eye(len(M));m=ia.point(M)
    dlo=np.diag(np.repeat([d['lower'] for d in comp['diagonal']],2))
    dhi=np.diag(np.repeat([d['upper'] for d in comp['diagonal']],2))
    minus=ia.sub(ia.sub(m,ia.point(em*eye)),ia.point(dhi))
    plus=ia.sub(ia.add(m,ia.point(em*eye)),ia.point(dlo))
    Gplus=ia.add(ia.point(G),ia.point(eg*eye))
    checks=[positive_certificate(minus,'compatibility',saved),positive_certificate(Gplus,'gram_upper',saved)]
    denominator=ia.sub(ia.point(comp['gamma']['lower']),ia.point(lam))
    B=ia.div(Gplus,denominator);S=ia.add(plus,B)
    midpoint=lambda a:(a[0]+a[1])/2
    X=solve(hermitian(midpoint(S)),midpoint(plus),assume_a='pos',check_finite=False)
    saved['conditional_trial_inverse_action']=X
    CX=ia.matmul(plus,ia.point(X))
    XSC=ia.matmul(ia.point(X.T),ia.matmul(S,ia.point(X)))
    shift=ia.add(ia.point(lam),ia.mul(ia.point(2.),ia.point(em)))
    diagonal=ia.sub(ia.point(dlo),ia.mul(shift,ia.point(eye)))
    trial=ia.sub(ia.add(ia.add(diagonal,CX),(CX[0].T,CX[1].T)),XSC)
    checks.append(positive_certificate(trial,'schur_lower',saved))
    return dict(passed=True,angular_offset=0,candidate_lower=lam,
        required_M_operator_error=em,required_G_upper_error=eg,checks=checks,
        interpretation='IF ||M_exact-M_saved||<=eM and G_exact<=G_saved+eG I, THEN the exact central-sector Dirichlet eigenvalue is >=0.18. Only this implication and its matrix arithmetic are certified here; the two integration-error premises remain unproved.',
        arithmetic='Outward IEEE754 interval contractions; saved X and trial bases are exact dyadics. No interval inversion or assumed exact orthogonality.')


def run():
    started=time.perf_counter();checks=controls();iv=MPIntervalContext();iv.dps=60
    constants={k:float(v.a) for k,v in target(iv).items()}
    _,function,_=evaluators(iv);D,gamma,comp=comparison(iv)
    old=project_path('dihole_heat_angular_state.npz')
    assert hashlib.sha256(old.read_bytes()).hexdigest()=='218545a22a2f1d9d757f7d7535c5534c82a2233e1cdac0b938dd885ceed01bdd'
    with np.load(old) as state:matrices=[state[f'M{i}'].copy() for i in range(3)]
    rows=[];grams={}
    for offset in (-128,0,128):
        M=matrices[0]+offset*matrices[1]+offset**2*matrices[2]
        G=residual_gram(function,constants,M,offset,60)
        refined=residual_gram(function,constants,M,offset,72)
        change=float(np.linalg.norm(refined-G,ord='fro'))
        pilot=spectrum_pilot(M,refined,D,gamma)
        # A fixed, larger-margin spectral check is more useful for later certification than a root at numerical zero.
        candidate=np.floor(1000*pilot['eigenvalue_brackets'][0]['conditional_lower'])/1000
        margin=float(low_values(schur_lower(M,refined,D,gamma,candidate),1)[0])
        pilot.update(offset=offset,quadrature_orders=[60,72],gram_refinement_frobenius_change=change,
            fixed_ground_candidate=candidate,fixed_ground_schur_margin=margin)
        rows.append(pilot);grams[f'j{offset}_gram']=refined
        print(json.dumps(dict(stage='compatible_schur',offset=offset,
            ground=pilot['eigenvalue_brackets'][0],gram_change=change,fixed_candidate=candidate,margin=margin)),flush=True)
    certificate=conditional_certificate(matrices[0],grams['j0_gram'],comp,grams)
    state_path=project_path('dihole_schur_spectrum_state.npz')
    np.savez_compressed(state_path,comparison_diagonal=D,**grams)
    names=['dihole_spatial_complement.py','dihole_spatial_complement_results.json',
        'dihole_heat_angular_state.npz','dihole_continuum_heat_tail.py','dihole_finite_heat.py',
        'dihole_quadrature_remainder.py','dihole_optical_heat.py','project_paths.py',
        'dihole_interval_arithmetic.py','dihole_heat_angular.py']
    result=dict(passed=True,controls=checks,comparison=comp,sectors=rows,conditional_matrix_certificate=certificate,
        target='Exact Dirichlet A4 Friedrichs square, N18 retained sines with both spins, q512000, mu0=398855.5.',
        theorem='A>=D, C=M-D_P>0, gamma=inf(D_Q)>lambda imply Sigma(lambda)<=C-C(C+G/(gamma-lambda))^-1C, with G=R*R. R=QAP.',
        qualification='The theorem includes the infinite spatial complement, but the saved eigenvalue/resolvent/low-heat intervals are numerical diagnostics conditional on enclosing the matrix and integral inputs. They are not certified continuum intervals.',
        residual_gram_note='For exact A and a fixed approximate Mhat, the full residual Gram equals R*R+(M-Mhat)^2, so a certified full Gram is a valid upper bound even with retained matrix defects.',
        limitations=['No new physical law or spectral novelty is claimed; the sharpening uses standard positive block forms and matrix parallel addition.',
            'Only three angular sectors and eight low levels per sector are tested; the heat values are partial traces, not full traces or pointwise diagonals.',
            'Quadrature refinement changes do not certify integration, coefficient evaluation, or matrix arithmetic.',
            'Natural exterior boundary comparison, renormalized stress, material support, backreaction, stability and backward passage remain open.'],
        source='https://doi.org/10.1184/R1/6479477',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        state_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_schur_spectrum_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,seconds=result['runtime_seconds'],controls=checks)),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11.4,4.8))
    center=next(s for s in result['sectors'] if s['offset']==0)
    rows=center['eigenvalue_brackets'];x=np.array([r['index'] for r in rows])
    low=np.array([r['conditional_lower'] for r in rows]);high=np.array([r['ritz_upper'] for r in rows])
    left.errorbar(x,(low+high)/2,yerr=(high-low)/2,fmt='o',capsize=5,label='Conditional interval')
    left.set(xlabel='Ordered level (central angular sector)',ylabel='Scaled squared energy',title='Low continuum spectrum: numerical brackets')
    left.legend(fontsize=8)
    for row in result['sectors']:
        heat=row['low_eight_heat_trace']
        right.plot([h['heat_parameter'] for h in heat],[100*h['relative_width'] for h in heat],'o-',label=f"Angular offset {row['offset']:+d}")
    right.set(xlabel='Heat parameter',ylabel='Conditional interval width / lower value (%)',title='Uncertainty of the eight-level partial trace')
    right.legend(fontsize=8)
    for ax in (left,right):ax.grid(alpha=.2)
    fig.suptitle('Using positivity to constrain omitted-mode coupling',fontsize=14)
    fig.text(.5,.02,'Analytic comparison verified; numerical integration is not enclosed. These are not full heat or stability results.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.94))
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_schur_spectrum_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
