"""Exact continuum spectral-envelope control for the conditional IR theorem.

The uniform density on [0,1] is an illustrative measure, not the unknown
spectral measure of the sewn bridge. No supporting stress is calculated.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from project_paths import project_path


def envelope(mass,epsilon):
    if not np.isfinite(mass+epsilon) or mass<0 or epsilon<=0:
        raise ValueError('Finite nonnegative mass and positive frequency cutoff required.')
    points=sorted({v for v in (mass,epsilon) if 0<v<1})
    return quad(lambda p:2/np.pi*np.arctan2(epsilon,np.hypot(p,mass)),
                0.,1.,points=points,epsabs=2e-13,epsrel=2e-11,limit=200)


def run():
    rows=[];analytic_errors=[];quadrature_errors=[]
    for epsilon in np.logspace(-6,-1,11):
        exact=2/np.pi*(np.arctan(epsilon)+epsilon/2*np.log1p(1/epsilon**2))
        numerical,error=envelope(0.,epsilon)
        analytic_errors.append(abs(numerical-exact));quadrature_errors.append(error)
        row=dict(epsilon=float(epsilon),gapless_uniform_bound=float(exact),paths={},zero_atom={})
        for power in (.5,1.,2.):
            mass=epsilon**power;value,error=envelope(mass,epsilon)
            quadrature_errors.append(error)
            row['paths'][str(power)]=dict(mass=float(mass),envelope=float(value))
            row['zero_atom'][str(power)]=float(2/np.pi*np.arctan2(epsilon,mass))
        rows.append(row)
    counter=[];target=4/np.pi*np.exp(-1j*np.pi/4)
    c8=np.pi/(8*np.sin(np.pi/8));cD=3**(3/8)/2
    for mass in np.logspace(-12,-2,6):
        z=np.sqrt(1j*mass);denominator=z*np.arctan(1/z)
        f=lambda x:2j*np.sqrt(mass)*x/(x**4+1j*mass)
        value=quad(lambda x:f(x).real,0.,1.,points=[mass**.25],epsabs=2e-12,epsrel=2e-12)[0]
        value+=1j*quad(lambda x:f(x).imag,0.,1.,points=[mass**.25],epsabs=2e-12,epsrel=2e-12)[0]
        v_error_sq=quad(lambda x:mass**2/(x**8+mass**2),0.,1.,points=[mass**.25],epsabs=2e-13)[0]
        coefficient=2*np.sqrt(mass)/denominator
        counter.append(dict(mass=float(mass),normalized_denominator_error=float(abs(value-denominator/np.sqrt(mass))),
            coefficient_real=float(coefficient.real),coefficient_imag=float(coefficient.imag),
            coefficient_error=float(abs(coefficient-target)),v_error_squared=float(v_error_sq),
            v_error_squared_bound=float(c8*mass**.25),
            scaled_inverse_operator_error_bound=float(cD*mass**.25+abs(coefficient-target)+2*abs(target)*np.sqrt(c8)*mass**.125)))
    checks=dict(
        gapless_integral_matches_exact=max(analytic_errors)<1e-11,
        all_paths_obey_mass_independent_bound=all(p['envelope']<=r['gapless_uniform_bound']+1e-11 for r in rows for p in r['paths'].values()),
        all_gapless_paths_decrease_to_zero=all(rows[i]['paths'][str(p)]['envelope']<rows[i+1]['paths'][str(p)]['envelope'] for p in (.5,1.,2.) for i in range(len(rows)-1)),
        smallest_cutoff_envelope_below_1e_5=rows[0]['gapless_uniform_bound']<1e-5,
        atom_equal_scales_half=all(abs(r['zero_atom']['1.0']-.5)<1e-15 for r in rows),
        atom_joint_limit_is_path_dependent=rows[0]['zero_atom']['2.0']-rows[0]['zero_atom']['0.5']>.99,
        numerical_quadrature_control=max(quadrature_errors)<1e-10,
        counterexample_denominator_identity=max(r['normalized_denominator_error'] for r in counter)<1e-10,
        counterexample_weighted_vector_bound=all(r['v_error_squared']<=r['v_error_squared_bound']+1e-12 for r in counter),
        counterexample_non_beta_phase=counter[0]['coefficient_error']<1e-5,
        counterexample_operator_bound_decreases=all(counter[i]['scaled_inverse_operator_error_bound']<counter[i+1]['scaled_inverse_operator_error_bound'] for i in range(len(counter)-1)))
    checks={key:bool(value) for key,value in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,rows=rows,
        abstract_unbounded_weight_counterexample=counter,
        abstract_counterexample_limit=dict(scaled_inverse_coefficient_real=float(target.real),scaled_inverse_coefficient_imag=float(target.imag),
            signed_mass_averaged_zero_probe_response=float(-1/np.sqrt(2)),scope='Bounded nonlocal K on L2(0,1), not the differential Dirac operator of the sewn geometry.'),
        maximum_exact_integral_error=max(analytic_errors),maximum_reported_quadrature_error=max(quadrature_errors),
        measure='Illustrative uniform continuous density on 0<=p<=1, plus a separate unit zero-atom control. No measure fitted to the sewn geometry.',
        scope='Controls a bounded-energy spectral envelope only. Atom path dependence is not a proof that the beta-odd physical stress fails; cancellation needs separate analysis.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs_sha256={'project_paths.py':hashlib.sha256(project_path('project_paths.py').read_bytes()).hexdigest()})
    project_path('dirac_infrared_corner_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']: raise ArithmeticError('Infrared controls failed; results retained.')
    plot(result)
    print(json.dumps({k:result[k] for k in ('passed','checks','maximum_exact_integral_error','maximum_reported_quadrature_error','abstract_counterexample_limit')},indent=2))
    print('Smallest-cutoff gapless envelope:',rows[0]['gapless_uniform_bound'])


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=result['rows'];eps=[r['epsilon'] for r in rows]
    fig,axes=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
    axes[0].loglog(eps,[r['gapless_uniform_bound'] for r in rows],'k--',label='Exact mass-independent bound')
    for power in (.5,1.,2.):
        label='m = epsilon'+('' if power==1 else f'^{power:g}')
        axes[0].loglog(eps,[r['paths'][str(power)]['envelope'] for r in rows],label=label)
        axes[1].semilogx(eps,[r['zero_atom'][str(power)] for r in rows],label=label)
    axes[0].set(title='Gapless continuum: all paths vanish',ylabel='Low-frequency envelope')
    axes[1].set(title='Exact zero atom: individual limits differ',ylabel='Normalized atom factor')
    for ax in axes[:2]:
        ax.set_xlabel('Frequency cutoff epsilon (decreases to the left)')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    counter=result['abstract_unbounded_weight_counterexample'];masses=[r['mass'] for r in counter]
    axes[2].loglog(masses,[r['scaled_inverse_operator_error_bound'] for r in counter],label='Analytic operator bound')
    axes[2].loglog(masses,[r['coefficient_error'] for r in counter],label='Rank-one coefficient error')
    axes[2].set(title='Abstract counterexample: wrong zero state',xlabel='Proper mass m (decreases to the left)',ylabel='Distance bound from non-beta limit')
    axes[2].grid(alpha=.2);axes[2].legend(fontsize=8)
    fig.suptitle('Low-energy limit: no spectral gap required; exact zero states need separate control',fontsize=11)
    fig.supxlabel('Illustrative continuum/atom measures and an abstract nonlocal counterexample. None is a bridge stress calculation.',fontsize=8)
    for ext in ('png','pdf'): fig.savefig(project_path('dirac_infrared_corner_progress.'+ext),dpi=160)
    plt.close(fig)


if __name__=='__main__': run()
