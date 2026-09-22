"""Full-field heat variation and a finite-regulator matching tail.

The comparison filter Lambda need not equal the physical regulator M.
The short-time remainder and the finite-M physical source remain unknown.
"""
import json
import math

import numpy as np
from mpmath import mp
from scipy.linalg import expm

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_pv_counterterms import MU_H_D0
from dihole_slice_channel import sha
from dirac_green_mass_limit import benchmark
from project_paths import project_path


NAME='dihole_regulator_heat_tail'


def insertion_integral(iv, record):
    """Use176's direct B trace-norm envelope, not the converse of its I bound."""
    low=record['selected_low_detail'];K=iv.mpf(record['selected']['cutoff_K_D'])
    x0=iv.mpf(low['scalar_HS_upper']);xs=iv.mpf(low['radial_HS_upper_D'])
    fN=iv.mpf(record['tail']['first_order_coefficient_N'])
    fP=iv.mpf(record['tail']['first_order_coefficient_P'])
    below=(fN*x0*x0*K*K/2+fP*x0*xs*K)/iv.pi
    above=sum(iv.mpf(a)/(int(p)-2)/K**(int(p)-2)
              for p,a in record['tail']['coefficient_upper_D_power'].items())/iv.pi
    return dict(low_upper_D2=upper(below),high_upper_D2=upper(above),
                complete_upper_D2=upper(below+above))


def heat_multiplier(iv, comparison_mass, time):
    """Positive-term spectral ceiling for ||A J(t) J(t)*||, units 1/D."""
    mass=iv.mpf(comparison_mass);s=iv.mpf(time)
    if mass.a<=0 or s.a<=0:
        raise ValueError('Positive comparison mass and proper time required.')
    y=mass*mass*s
    return mass*sum(p*((iv.mpf(k)+iv.mpf('.5'))/(iv.e*y))**(iv.mpf(k)+iv.mpf('.5'))
                    for k,p in enumerate((iv.mpf(1),iv.mpf(11)/6,iv.mpf(1),iv.mpf(1)/6)))


def matching_weights(iv, mass, cutoff):
    """Elementary positive-integral weights; no infrared shift is present."""
    M=iv.mpf(mass);s=iv.mpf(cutoff)
    if M.a<=0 or s.a<=0:
        raise ValueError('Positive target mass and cutoff required.')
    z=M*M*s
    tail=sum(iv.mpf(c)*iv.exp(-j*z)/j for j,c in enumerate((3,-3,1),1))/(M*M)
    J0=iv.mpf(11)/(6*M*M)-tail
    J1=sum(iv.mpf(c)*(1-(1+j*z)*iv.exp(-j*z))/(j*j)
           for j,c in enumerate((3,-3,1),1))/M**4
    return tail,J0,J1


def local_coefficients(q, local, uv):
    """Numerical saved geometry coefficients, not certified enclosures."""
    row=next(r for r in local['records'] if r['order']==128)
    action=list(map(mp.mpf,row['action_derivative']))
    a2=-action[0]/3
    a4=action[1]/30-action[2]/10-mp.mpf(2)/3*mp.mpf(row['charge_included_F2_energy_derivative'])*q*q
    poly=uv['fixed_geometry_charge_polynomial']
    a6=mp.mpf(poly['gravity'])+mp.mpf(poly['gauge'])*q*q
    return a2,a4,a6


def matching_center(q, mass, cutoff, local, uv):
    """Correction relative to Q_M-C_M; the required Q_M is not evaluated."""
    M=mp.mpf(mass);s=mp.mpf(cutoff);z=M*M*s
    a2,a4,a6=local_coefficients(q,local,uv)
    I1=mp.fsum(c*mp.e1(j*z) for j,c in enumerate((3,-3,1),1))
    I2=mp.fsum(c*(mp.exp(-j*z)/s-j*M*M*mp.e1(j*z)) for j,c in enumerate((3,-3,1),1))
    J0=mp.fsum(mp.mpf(c)/j*(-mp.expm1(-j*z)) for j,c in enumerate((3,-3,1),1))/(M*M)
    leading=a6*J0/(32*mp.pi**2);local_tail=(a2*I2+a4*I1)/(32*mp.pi**2)
    mu=mp.mpf(q)*mp.mpf(MU_H_D0)
    counterterm=((6*mp.log(2)-3*mp.log(3))*a2*M*M
                 +a4*(mp.log(M*M/(mu*mu))+mp.log(mp.mpf(3)/8)))/(32*mp.pi**2)
    return dict(a2_estimate_over_D=float(a2),a4_estimate_D=float(a4),a6_estimate_D3=float(a6),
        leading_short_correction_estimate_D=float(leading),local_tail_estimate_D=float(local_tail),
        correction_center_estimate_D=float(leading-local_tail),counterterm_estimate_D=float(counterterm),
        fixed_reference_mu_D=float(mu),I1_estimate=float(I1),I2_estimate_D2=float(I2))


def matrix_controls():
    rows=[];eye=np.eye(6)
    for largest in (10.,1000.):
        for kappa in (0.,.3,2.):
            _,beta,N,V,dN,A,P,dP=benchmark(largest,kappa)
            C=np.diag(np.sqrt(np.diag(N)));inverse=np.linalg.inv(C);F=np.linalg.solve(N,dN)
            dA=beta@inverse@V@inverse-(F@A+A@F)/2
            for mass in (.7,2.):
                K=np.sqrt(6)*mass**3*np.linalg.multi_dot([np.linalg.inv(A+mass*np.sqrt(j)*eye) for j in (1,2,3)])
                D=K.conj().T@dA@K;B=beta@C@D@C;nu=np.linalg.norm(np.linalg.inv(N),2)
                for time in (.03,1.,10.):
                    E=expm(-time*P/2);J=np.linalg.solve(K,E)
                    direct=-time*np.trace(dP@E@E).real
                    cyclic=2*time*np.trace(D@A@J@J.conj().T).real
                    iv=MPIntervalContext();iv.dps=50
                    multiplier=upper(heat_multiplier(iv,mass,time))
                    bound=2*time*nu*np.linalg.norm(B,'nuc')*multiplier
                    rows.append(dict(largest_lapse=largest,kappa=kappa,mass=mass,time=time,
                        cyclic_scaled_error=abs(direct-cyclic)/max(1,abs(direct)),
                        weight_identity_error=float(np.linalg.norm(D-inverse@beta@B@inverse)),
                        multiplier_ratio=float(np.linalg.norm(A@J@J.conj().T,2))/multiplier,
                        heat_bound_ratio=float(abs(direct)/bound)))
    return rows


def scalar_controls():
    rows=[]
    for y in map(mp.mpf,('0.0001','.01','1','10','100')):
        p=lambda z:(z+1)*(z+2)*(z+3)/6
        root=mp.findroot(lambda z:1/(2*z)+sum(1/(z+j) for j in (1,2,3))-y,
                         (1/(2*y),mp.mpf(7)/(2*y)),solver='bisect',tol=mp.mpf('1e-70'),maxsteps=300)
        exact=mp.sqrt(root)*p(root)*mp.exp(-y*root)
        iv=MPIntervalContext();iv.dps=60;bound=upper(heat_multiplier(iv,1.,float(y)))
        rows.append(dict(y=float(y),maximum_ratio=float(exact/bound)))
    weights=[]
    for z in map(mp.mpf,('.1','1','10','100','200')):
        # Scale time by s0=1; the integration interval stays well conditioned.
        r=lambda s:3*mp.exp(-z*s)-3*mp.exp(-2*z*s)+mp.exp(-3*z*s)
        exact=[mp.fsum(mp.mpf(c)*mp.exp(-j*z)/j for j,c in enumerate((3,-3,1),1))/z,
               mp.fsum(mp.mpf(c)*(-mp.expm1(-j*z))/j for j,c in enumerate((3,-3,1),1))/z,
               mp.fsum(mp.mpf(c)*(1-(1+j*z)*mp.exp(-j*z))/(j*j) for j,c in enumerate((3,-3,1),1))/z**2]
        quadrature=[mp.exp(-z)*mp.quad(lambda u:mp.exp(z)*r(1+u/z),[0,1,10,mp.inf])/z,
                    mp.quad(r,[0,1]),mp.quad(lambda s:s*r(s),[0,1])]
        weights.append(dict(z=float(z),maximum_relative_error=float(max(abs(a-b)/abs(a) for a,b in zip(exact,quadrature)))))
    return rows,weights


def run():
    mp.dps=85;iv=MPIntervalContext();iv.dps=60
    path=project_path('dihole_full_frequency_results.json');prior=json.loads(path.read_text())
    inputs={**prior['input_sha256'],path.name:sha(path),'dihole_full_frequency.py':prior['source_sha256']}
    if not prior['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing full-frequency inputs required.')
    geometry=json.loads(project_path('dihole_buffer_smoothing_results.json').read_text())['geometry']
    nu=json.loads(project_path('dihole_paired_low_frequency_results.json').read_text())['global_lapse']['inverse_lapse_upper']
    local=json.loads(project_path('dihole_local_source_results.json').read_text())
    uv=json.loads(project_path('dihole_uv_bound_results.json').read_text())
    records=[]
    for source in (prior['records'][0],prior['records'][2]):
        q=source['q'];comparison=source['regulator_M_D'];W=insertion_integral(iv,source)
        magnetic=iv.mpf(q)*iv.mpf(geometry['magnetic_unit_flux_upper_D2'])
        cutoff=lower(iv.mpf('.01')/magnetic)
        L=upper(heat_multiplier(iv,comparison,cutoff));rows=[]
        for nominal_z in (1,10,30,50,100,200):
            mass=upper(iv.sqrt(iv.mpf(nominal_z)/iv.mpf(cutoff)))
            tail,J0,J1=matching_weights(iv,mass,cutoff)
            E=iv.mpf(nu)*iv.mpf(W['complete_upper_D2'])*iv.mpf(L)*tail
            coefficient=J1/(32*iv.pi**2)
            rows.append(dict(nominal_M_squared_s0=nominal_z,target_regulator_M_D=mass,
                actual_M_squared_s0_interval=[lower(iv.mpf(mass)**2*iv.mpf(cutoff)),upper(iv.mpf(mass)**2*iv.mpf(cutoff))],
                long_heat_correction_upper_D=upper(E),long_heat_correction_upper_per_q2=upper(E/(q*q)),
                long_weight_upper_D_squared=upper(tail),short_constant_weight_upper_D_squared=upper(J0),
                short_linear_weight_upper_D4=upper(J1),Cstar_error_coefficient_upper_D4=upper(coefficient),
                local=matching_center(q,mass,cutoff,local,uv)))
        records.append(dict(q=q,comparison_filter_Lambda_D=comparison,W=W,inverse_lapse_upper=nu,
            cutoff_s0_over_D_squared=cutoff,magnetic_product_upper=upper(magnetic*iv.mpf(cutoff)),
            magnetic_scale_choice_valid=bool((magnetic*iv.mpf(cutoff)).b<=iv.mpf('.01').a),
            spectral_multiplier_upper_D=L,heat_variation_at_cutoff_upper_D=upper(2*iv.mpf(cutoff)*iv.mpf(nu)*iv.mpf(W['complete_upper_D2'])*iv.mpf(L)),
            matching_rows=rows,actual_short_time_Cstar=None,finite_target_regulator_source=None,renormalized_source=None))
    matrices=matrix_controls();multipliers,weights=scalar_controls()
    rejected=[]
    for mass,time in ((0.,1.),(1.,0.),(-1.,1.)):
        try:heat_multiplier(iv,mass,time)
        except ValueError:rejected.append(True)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        positive_inputs_required=len(rejected)==3,
        finite_integrated_insertion=all(math.isfinite(r['W']['complete_upper_D2']) and r['W']['complete_upper_D2']>0 for r in records),
        zero_shift_only=all(prior['records'][i]['eta_ratio_interval']==[0.,0.] for i in (0,2)),
        heat_cyclic_identity=max(r['cyclic_scaled_error'] for r in matrices)<1e-10,
        lapse_weight_identity=max(r['weight_identity_error'] for r in matrices)<1e-12,
        finite_matrix_multiplier=max(r['multiplier_ratio'] for r in matrices)<=1+1e-12,
        finite_matrix_heat_bound=max(r['heat_bound_ratio'] for r in matrices)<=1+1e-12,
        continuous_spectral_maximum=max(r['maximum_ratio'] for r in multipliers)<=1,
        matching_weight_integrals=max(r['maximum_relative_error'] for r in weights)<1e-55,
        magnetic_cutoff_scale=all(r['magnetic_scale_choice_valid'] for r in records),
        distinct_comparison_and_target=all(row['target_regulator_M_D']!=r['comparison_filter_Lambda_D'] for r in records for row in r['matching_rows']),
        positive_error_weights=all(row['Cstar_error_coefficient_upper_D4']>0 and row['long_heat_correction_upper_D']>0 for r in records for row in r['matching_rows']),
        long_heat_bounds_decrease=all(all(a['long_heat_correction_upper_D']>b['long_heat_correction_upper_D'] for a,b in zip(r['matching_rows'],r['matching_rows'][1:])) for r in records))
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,records=records,matrix_controls=matrices,
        spectral_controls=multipliers,weight_controls=weights,input_sha256=inputs,source_sha256=sha(__file__),
        scope=['Actual full-field heat-variation and long-time correction bounds, conditional on the previously stated operator model; no angular or exterior truncation.',
               'eta=0 after the established ordered infrared prescription. A fixed positive eta would change the matching kernel.',
               'Lambda is an auxiliary comparison filter. The matching identity needs the physical source at Mtarget, which is not supplied by the saved Lambda source bound.',
               'Numerical local coefficients and correction centers are not interval-certified geometry data.',
               'The magnetic cutoff scale is a declared choice, not a proof of a local heat remainder bound.',
               'Short-time Cstar, finite target-regulator source, renormalized stress/sign, stability and atom transport remain uncomputed.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Heat matching controls failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,atomic_W=records[-1]['W'],atomic_rows=[{k:v for k,v in row.items() if k!='local'} for row in records[-1]['matching_rows']]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=result['records'][-1];rows=r['matching_rows'];z=[row['nominal_M_squared_s0'] for row in rows]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.semilogy(z,[row['long_heat_correction_upper_per_q2'] for row in rows],'o-',color='#28718e')
    ax.set(xlabel='Target regulator squared times proper-time cutoff',ylabel='Long-time correction ceiling times D / q²',
           title='One previously unknown correction now has a finite bound')
    ax.grid(alpha=.2)
    fig.supxlabel('Atomic-flux supporting field. Short-time remainder and total physical source remain unknown.\nThis is not a passage probability or an energy prediction.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
