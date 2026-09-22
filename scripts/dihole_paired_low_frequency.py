"""Pair the selected proper insertion at low frequency, with zero states retained.

The fixed-regulator signed integral is bounded; no complete physical source,
regulator limit, atom propagation or clock gain is inferred.
"""
import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from dihole_channel_bounds import MPIntervalContext, calibration, lower, upper
from dihole_local_zero_shift import components
from dihole_slice_channel import sha
from dirac_green_mass_limit import benchmark
from project_paths import project_path


NAME = 'dihole_paired_low_frequency'


def global_lapse(iv):
    k, mass, ag, _, _ = calibration(iv)
    epsilon = iv.mpf(1) / 40
    h = epsilon * (2*k-epsilon)
    floor = (ag/k)**2 * h / (h+2*mass*(mass+k))
    return dict(lapse_lower=lower(floor),inverse_lapse_upper=upper(1/floor),
                geometric_mass_over_D=[lower(mass),upper(mass)],
                focal_half_distance_over_D=.5,minimum_focal_distance_over_D=.025,
                ratio_identity_enclosed=bool((h/(h+2*mass*(mass+k))-iv.mpf(39)/359).a<=0
                                            <=(h/(h+2*mass*(mass+k))-iv.mpf(39)/359).b))


def low_constants(iv, M, ca, da, local_lower, local_upper, nu, hs2, eta_cap=0):
    if M.a<=0 or local_lower.a<=0 or local_upper.a<local_lower.b or nu.a<=0 or hs2.a<0:
        raise ValueError('Positive regulator and lapse bounds, nonnegative HS square required.')
    parts=components(iv,iv.mpf(0),iv.mpf(0),iv.mpf(1),iv.mpf(eta_cap),M)
    regular=sum(abs(parts['residues'][j])/(iv.sqrt(j)*M) for j in (1,2,3))
    singular=parts['residues'][0]
    t=1+1/iv.sqrt(2)+1/iv.sqrt(3)
    P=ca*M+da
    lipschitz=nu*t*(ca+da/M)+4/local_lower
    frozen=hs2*singular*local_upper*P/2
    slope=hs2*(singular*local_upper*lipschitz+P*regular)/iv.pi
    return dict(regular_remainder_upper_over_D=upper(regular),
                compact_left_leg_upper_D=upper(P),compact_left_leg_lipschitz_upper=upper(lipschitz),
                paired_low_constant_upper_D=upper(frozen),paired_low_slope_upper=upper(slope))


def combined(iv, constants, coefficients, cutoff):
    K=iv.mpf(cutoff)
    if K.a<=0:raise ValueError('Require a positive frequency split.')
    low=iv.mpf(constants['paired_low_constant_upper_D'])+K*constants['paired_low_slope_upper']
    high=sum(iv.mpf(v)/((int(p)-1)*K**(int(p)-1)) for p,v in coefficients.items())/iv.pi
    return dict(split_K_D=cutoff,paired_low_integral_absolute_upper_D=upper(low),
                absolute_high_integral_upper_D=upper(high),
                complete_selected_signed_integral_absolute_upper_D=upper(low+high))


def bounds(iv, previous, proper, lapse, global_bound):
    M=previous['fixed_regulator_M_D']
    args=list(map(iv.mpf,(M,proper['proper_metric_graph_a_upper'],
                         proper['proper_metric_graph_zero_upper_D'],lapse['lapse_lower'],
                         lapse['lapse_upper'],global_bound['inverse_lapse_upper'],
                         proper['filter_HS_squared_upper'])))
    constants=low_constants(iv,*args)
    coefficients=previous['zero_shift_frequency_tail']['absolute_coefficient_upper_D_power']
    # Optimize only the proved upper envelope. The resulting binary cutoff
    # is re-evaluated with interval arithmetic; it is not a physical optimum.
    scaled={int(p):v/M**(int(p)-1) for p,v in coefficients.items()}
    target=constants['paired_low_slope_upper']*M*math.pi
    balance=lambda z:sum(a*math.exp(-p*z) for p,a in scaled.items())-target
    left,right=-1.,1.
    while balance(left)<0:left-=1
    while balance(right)>0:right+=1
    ratio=math.exp(brentq(balance,left,right,xtol=1e-13))
    rows=[combined(iv,constants,coefficients,f*M) for f in (.1,1.,10.)]
    selected=combined(iv,constants,coefficients,ratio*M)
    cap=previous['uniform_eta_frequency_tail']['eta_cap_D2']
    uniform=low_constants(iv,*args,eta_cap=cap)
    uniform['eta_cap_D2']=cap
    uniform['same_split']=combined(iv,uniform,previous['uniform_eta_frequency_tail']['absolute_coefficient_upper_D_power'],ratio*M)
    return dict(q=previous['q'],fixed_regulator_M_D=M,
                fixed_filter_Lambda_D=previous['fixed_filter_Lambda_D'],
                **constants,cutoff_controls=rows,selected_split=selected,
                selected_K_over_M=ratio,uniform_eta_interval=uniform)


def matrix_controls(eta_ratio=0.):
    rows=[];integrals=[];spectral=[];rng=np.random.default_rng(173)
    f=(rng.normal(size=(6,2))+1j*rng.normal(size=(6,2)))/4
    norm=lambda a:float(np.linalg.norm(a,2))
    hs2=float(np.linalg.norm(f,'fro')**2)
    t=1+1/math.sqrt(2)+1/math.sqrt(3)
    for largest in (2.,10.,100.):
        H,beta,N,V,dN,A0,_,_=benchmark(largest,0.)
        root=np.diag(np.sqrt(np.diag(N)));B=np.diag(1/np.sqrt(np.diag(N)))
        F=np.linalg.solve(N,dN);eye=np.eye(6);nu=1/min(np.diag(N));Nmax=max(np.diag(N))
        D=beta@B@V@B-(A0@F-F@A0)/2
        ca=norm(F);da=norm(D);slope=norm(F)*nu
        Y=root@beta@f
        for M in (.7,2.):
            all_m=np.sqrt(eta_ratio+np.arange(4))*M
            m0=all_m[0];masses=all_m[1:];alpha=math.sqrt(6)*M**3
            Hm=H+m0*N@beta
            values,vectors=np.linalg.eigh(Hm)
            # The existing massless benchmark has two exact zero eigenvalues.
            zero=abs(values)<1e-12
            if eta_ratio==0:values[zero]=0.
            spectral.append(dict(largest_lapse=largest,M=M,eta_over_M_squared=eta_ratio,
                                 zero_states=int(sum(zero)),eigen_residual=norm(Hm@vectors-vectors*values)))
            residues=[alpha*math.prod((all_m[j]+all_m[k])/((k-j)*M*M) for k in range(4) if j!=k) for j in range(4)]
            regular=sum(abs(w)/m for w,m in zip(residues[1:],masses))
            singular=residues[0]
            P=ca*M+da
            lip=nu*t*(ca+da/M)+slope
            def data(kappa):
                A=A0+1j*kappa*beta@np.linalg.inv(N)
                dot=D-F@A
                R=[np.linalg.inv(A+m*eye) for m in masses]
                K=alpha*R[0]@R[1]@R[2]
                L=sum(w*r for w,r in zip(residues[1:],R))
                p=dot@K@f
                return A,dot,K,L,p
            _,_,K0,_,p0=data(0.)
            ell0=root@p0
            def resolvent(kappa):
                return (vectors*(1/(values+1j*kappa)))@vectors.conj().T
            def source(kappa):
                A,dot,K,L,p=data(kappa)
                u=singular*root@resolvent(kappa)@Y+L@f
                direct=-float(np.vdot(K@f-m0*u,dot@u).real)
                moved=float(singular*np.vdot(root@p,resolvent(kappa)@Y).real+np.vdot(p,L@f).real)
                return direct,moved
            def paired(kappa):
                paired_inverse=(vectors*(2*values/(values*values+kappa*kappa)))@vectors.conj().T
                value=singular*float(np.vdot(ell0,paired_inverse@Y).real)
                for sign in (-1,1):
                    _,_,_,L,p=data(sign*kappa)
                    value+=float(singular*np.vdot(root@(p-p0),resolvent(sign*kappa)@Y).real+np.vdot(p,L@f).real)
                return value
            for kappa in (-3.,-.1,-1e-4,1e-4,.1,3.):
                A,dot,K,L,p=data(kappa);delta=K-K0
                direct,moved=source(kappa)
                factor=alpha*np.linalg.inv(A+m0*eye)@np.linalg.inv(A+masses[0]*eye)@np.linalg.inv(A+masses[1]*eye)@np.linalg.inv(A+masses[2]*eye)
                rows.append(dict(largest_lapse=largest,M=M,kappa=kappa,eta_over_M_squared=eta_ratio,
                    K_norm=norm(K),action_norm_over_M=norm(A@K)/M,
                    split_relative_error=norm(factor-(singular*np.linalg.inv(A+m0*eye)+L))/max(1,norm(factor)),
                    K_difference_ratio=norm(delta)/(abs(kappa)*nu*t/M),
                    K_graph_difference_ratio=norm(A@delta)/(abs(kappa)*nu*t),
                    left_leg_difference_ratio=float(np.linalg.norm(p-p0,'fro'))/(abs(kappa)*lip*math.sqrt(hs2)),
                    local_integration_by_parts_relative_error=abs(direct-moved)/max(1,abs(direct),abs(moved)),
                    pair_relative_error=abs(source(abs(kappa))[0]+source(-abs(kappa))[0]-paired(abs(kappa)))/max(1,abs(paired(abs(kappa))))))
            for cutoff in (.2,1.):
                value,error=quad(paired,0,cutoff,epsabs=2e-9,epsrel=2e-9)
                value/=2*math.pi
                absolute,_=quad(lambda k:abs(paired(k)),0,cutoff,epsabs=2e-9,epsrel=2e-9)
                absolute/=2*math.pi
                limit=hs2*(singular*Nmax*P/2+cutoff*(singular*Nmax*lip+P*regular)/math.pi)
                signs=np.sign(values);weights=np.zeros_like(values)
                np.divide(cutoff,abs(values),out=weights,where=values!=0)
                weights=2*signs*np.arctan(weights)
                analytic=float(np.vdot(ell0,(vectors*weights)@vectors.conj().T@Y).real)
                frozen,_=quad(lambda k:float(np.vdot(ell0,(vectors*(2*values/(values*values+k*k)))@vectors.conj().T@Y).real),0,cutoff)
                integrals.append(dict(largest_lapse=largest,M=M,cutoff=cutoff,eta_over_M_squared=eta_ratio,
                    paired_integral=value,absolute_paired_integral=absolute,absolute_bound=limit,
                    enclosed=bool(abs(value)<=limit),absolute_pair_enclosed=bool(absolute<=limit),
                    quadrature_error=error/(2*math.pi),frozen_spectral_integral_error=abs(analytic-frozen),
                    frozen_bound_holds=bool(abs(analytic)<=math.pi*Nmax*P*hs2)))
    return rows,integrals,spectral


def run():
    path=project_path('dihole_local_zero_shift_results.json');previous=json.loads(path.read_text())
    inputs=dict(previous['input_sha256'])
    inputs.update({path.name:sha(path),'dihole_local_zero_shift.py':previous['source_sha256'],
                   'dirac_green_mass_limit.py':sha(project_path('dirac_green_mass_limit.py'))})
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing prior inputs required.')
    proper=json.loads(project_path('dihole_proper_channel_results.json').read_text())
    massive=json.loads(project_path('dihole_massive_channel_results.json').read_text())
    iv=MPIntervalContext();iv.dps=60
    global_bound=global_lapse(iv)
    records=[bounds(iv,p,a,m['lapse'],global_bound) for p,a,m in zip(previous['records'],proper['records'],massive['records'],strict=True)]
    matrices,integrals,spectral=matrix_controls()
    positive,positive_integrals,positive_spectral=matrix_controls(.01)
    matrices+=positive;integrals+=positive_integrals;spectral+=positive_spectral
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        positive_global_lapse=0<global_bound['lapse_lower']<min(m['lapse']['lapse_lower'] for m in massive['records']),
        exact_geometric_ratio=global_bound['ratio_identity_enclosed'],
        benchmark_zero_states_retained=all(r['zero_states']==2 and r['eigen_residual']<1e-12 for r in spectral if r['eta_over_M_squared']==0),
        positive_shift_benchmark=all(r['zero_states']==0 and r['eigen_residual']<1e-12 for r in spectral if r['eta_over_M_squared']>0),
        bounded_action_filter=max(r['K_norm'] for r in matrices)<=1+1e-12,
        first_action_bound=max(r['action_norm_over_M'] for r in matrices)<=1+1e-12,
        positive_mass_remainder_identity=max(r['split_relative_error'] for r in matrices)<1e-9,
        resolvent_product_lipschitz=max(r['K_difference_ratio'] for r in matrices)<=1+1e-10,
        graph_product_lipschitz=max(r['K_graph_difference_ratio'] for r in matrices)<=1+1e-10,
        compact_leg_lipschitz=max(r['left_leg_difference_ratio'] for r in matrices)<=1+1e-10,
        local_integration_by_parts=max(r['local_integration_by_parts_relative_error'] for r in matrices)<1e-7,
        paired_source_identity=max(r['pair_relative_error'] for r in matrices)<1e-7,
        paired_low_frequency_enclosed=all(r['enclosed'] for r in integrals),
        absolute_paired_integral_enclosed=all(r['absolute_pair_enclosed'] for r in integrals),
        frozen_spectral_integral=max(r['frozen_spectral_integral_error'] for r in integrals)<1e-9,
        zero_state_safe_spectral_bound=all(r['frozen_bound_holds'] for r in integrals),
        full_selected_frequency_bound=all(math.isfinite(r['selected_split']['complete_selected_signed_integral_absolute_upper_D']) for r in records),
        uniform_shift_covers_endpoint=all(r['uniform_eta_interval']['same_split']['complete_selected_signed_integral_absolute_upper_D']>=r['selected_split']['complete_selected_signed_integral_absolute_upper_D'] for r in records),
        chosen_split_improves_samples=all(r['selected_split']['complete_selected_signed_integral_absolute_upper_D']<=min(x['complete_selected_signed_integral_absolute_upper_D'] for x in r['cutoff_controls']) for r in records))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,global_lapse=global_bound,records=records,
        matrix_controls=matrices,integral_controls=integrals,spectral_controls=spectral,
        source_sha256=sha(__file__),input_sha256=inputs,physical_source=None,source_sign=None,atom_transit=None,
        scope=['The paired selected insertion is bounded over the complete frequency axis at fixed M and Gaussian selection, uniformly over the saved eta interval.',
               'Low-frequency pairing is essential; no assertion that the two unsigned integrals or the absolute pointwise integral exist.',
               'Zero states are retained with the symmetric frequency prescription. No gap or vanishing zero projection is assumed.',
               'Uniform integral boundedness does not establish uniform integrability near zero or interchange of the eta and frequency limits. Regulator removal, omitted-selection sum, counterterm matching and full physical stress remain unresolved.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Paired low-frequency checks failed; evidence retained.')
    plot(result,previous['records'][-1])
    print(json.dumps(dict(passed=result['passed'],tests=tests,global_lapse=global_bound,atomic=records[-1]),indent=2))


def plot(result,old):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    iv=MPIntervalContext();iv.dps=40
    atomic=result['records'][-1];M=atomic['fixed_regulator_M_D']
    ratios=np.geomspace(.03,10,120)
    rows=[combined(iv,atomic,old['zero_shift_frequency_tail']['absolute_coefficient_upper_D_power'],float(r*M)) for r in ratios]
    fig,ax=plt.subplots(figsize=(8.8,4.7),layout='constrained')
    for key,label in (('paired_low_integral_absolute_upper_D','Paired low-frequency bound'),
                      ('absolute_high_integral_upper_D','Absolute high-frequency tail'),
                      ('complete_selected_signed_integral_absolute_upper_D','Combined selected bound')):
        ax.loglog(ratios,[r[key] for r in rows],label=label)
    ax.axvline(atomic['selected_K_over_M'],color='gray',linestyle=':',label='Chosen frequency split')
    ax.set(xlabel='Frequency split K / fixed regulator M',ylabel='Bound multiplied by D',
           title='One selected contribution now has a full-frequency bound')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.supxlabel('Atomic reference. Fixed regulator and spatial filter; total support, stability and atom passage remain unresolved.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
