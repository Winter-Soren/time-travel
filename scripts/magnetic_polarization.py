"""Local magnetic log-action benchmark and a separate classical-RN trace audit.

This is not a full fermion determinant, a combined LLL source, or mouth matching.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sy


def metric_from_horizon_offset(offset, epsilon):
    """Evaluate f at R=1+offset without subtracting order-one RN terms.

    Passing the offset directly retains values smaller than spacing near R=1.
    """
    u=np.asarray(offset,dtype=float)
    if np.any(u<=-1):raise ValueError('Require positive areal radius')
    flat=u.ravel();phi=np.empty_like(flat);small=abs(flat)<.01
    z=flat[small]
    # (log(1+u)-u)/u^2 = -1/2+u/3-u^2/4+... .
    value=np.zeros_like(z)
    for power in range(12,-1,-1):value=value*z+(-1.)**(power+1)/(power+2)
    phi[small]=value
    z=flat[~small];phi[~small]=(np.log1p(z)-z)/(z*z)
    result=(flat/(1+flat))**2*(1+4*epsilon*phi)
    return result.reshape(u.shape)


def magnetic_polarization_checks():
    source=project_path('mmp_pulse_results.json')
    saved=json.loads(source.read_text(encoding='utf-8'))['parameters']
    q_value=int(saved['magnetic_charge_q']);g_value=float(saved['gauge_coupling_g'])
    epsilon_value=g_value*g_value/(24*math.pi**2)
    if not 0<epsilon_value<.5:raise ValueError('Benchmark requires the positive-curvature perturbative branch')
    R,X,Xstar,g,q,epsilon=sy.symbols('R X Xstar g q epsilon',positive=True)
    pi=sy.pi;k=1/(96*pi**2)
    lagrangian=-X/(4*g**2)+k*X*sy.log(X/Xstar)
    lx=sy.diff(lagrangian,X)
    rho=-lagrangian;pr=lagrangian;pt=lagrangian-2*X*lx
    trace=sy.simplify(-rho+pr+2*pt)
    x_field=q**2/(2*R**4);x_reference=q**2/2
    # Set re=1 and impose re^2=pi*G*q^2/g^2 while deriving normalized sources.
    G=g**2/(pi*q**2)
    def normalize(stress):
        value=sy.simplify(8*pi*G*stress.subs({X:x_field,Xstar:x_reference}))
        return sy.simplify(value.subs(g**2,24*pi**2*epsilon))
    density,radial,angular=map(normalize,(rho,pr,pt))
    mass=1+2*epsilon-(sy.Rational(1,2)+2*epsilon*(sy.log(R)+1))/R
    f=sy.simplify(1-2*mass/R)
    f1,f2=sy.diff(f,R),sy.diff(f,R,2)
    einstein_tt=f1/R+(f-1)/R**2
    einstein_rr=einstein_tt
    einstein_angular=f2/2+f1/R
    residuals={
        'magnetic_radial_nec':sy.simplify(rho+pr),
        'gauge_trace':sy.simplify(trace+X/(24*pi**2)),
        'density_normalization':sy.simplify(density-(1+4*epsilon*sy.log(R))/R**4),
        'radial_stress_conservation':sy.simplify(sy.diff(radial,R)+
            (density+radial)*f1/(2*f)+2*(radial-angular)/R),
        'mass_equation':sy.simplify(sy.diff(mass,R)-R**2*density/2),
        'einstein_tt':sy.simplify(einstein_tt+density),
        'einstein_rr':sy.simplify(einstein_rr-radial),
        'einstein_angular':sy.simplify(einstein_angular-angular),
        'horizon_value':sy.simplify(f.subs(R,1)),
        'horizon_first_derivative':sy.simplify(f1.subs(R,1)),
        'horizon_second_derivative':sy.simplify(f2.subs(R,1)-2*(1-2*epsilon)),
        'product_sphere_equation':sy.simplify(density.subs(R,1)-1),
        'product_ads_equation':sy.simplify(angular.subs(R,1)-(1-2*epsilon)),
        'weak_coupling_metric':sy.simplify(f.subs(epsilon,0)-(1-1/R)**2),
        'linear_mass_correction':sy.simplify(sy.diff(mass,epsilon,2))}
    # Proper-distance Einstein components give the flare-out identity directly.
    radius,N,Rl,Rll,Nl=sy.symbols('radius N Rl Rll Nl',positive=True)
    rho_g=(1-Rl**2-2*radius*Rll)/radius**2
    pr_g=(-1+Rl**2+2*radius*Rl*Nl/N)/radius**2
    residuals['proper_distance_flare_out']=sy.simplify(Rll/N-Rl*Nl/N**2+
                                                     radius*(rho_g+pr_g)/(2*N))

    # This second audit uses CLASSICAL extremal RN, not the corrected f above.
    f_rn=(1-1/R)**2;fp=sy.diff(f_rn,R);fpp=sy.diff(f_rn,R,2)
    ricci_scalar=sy.simplify(-fpp-4*fp/R+2*(1-f_rn)/R**2)
    ricci_squared=sy.simplify(2*(fpp/2+fp/R)**2+2*((1-f_rn-R*fp)/R**2)**2)
    kretschmann=sy.simplify(fpp**2+4*(fp/R)**2+4*((1-f_rn)/R**2)**2)
    weyl=sy.simplify(kretschmann-2*ricci_squared+ricci_scalar**2/3)
    euler=sy.simplify(kretschmann-4*ricci_squared+ricci_scalar**2)
    gravity_trace=sy.simplify((weyl/20-sy.Rational(11,360)*euler)/(16*pi**2))
    gravity_polynomial=(42*R**2-84*R+53)/(720*pi**2*R**8)
    gauge_trace=-q**2/(48*pi**2*R**4)
    lll_trace=q*(4*R-6)/(96*pi**2*R**6)
    required_trace=sy.simplify(gravity_trace-lll_trace)
    residuals.update(classical_rn_ricci_scalar=ricci_scalar,
        classical_rn_ricci_squared=sy.simplify(ricci_squared-4/R**8),
        classical_rn_weyl_squared=sy.simplify(weyl-48*(R-1)**2/R**8),
        classical_rn_euler=sy.simplify(euler-(48*R**2-96*R+40)/R**8),
        gravity_trace_polynomial=sy.simplify(gravity_trace-gravity_polynomial),
        gravity_trace_second_primary=sy.simplify(gravity_trace-(7*kretschmann+8*ricci_squared)/(360*(4*pi)**2)),
        required_trace_horizon=sy.simplify(required_trace.subs(R,1)-q/(48*pi**2)-11/(720*pi**2)))
    assert all(value==0 for value in residuals.values()),residuals

    radii=(1.,1.000001,1.01,1.1,1.5,2.,10.,100.)
    table=[];mass_errors=[];metric_errors=[];near=[]
    with mp.workdps(70):
        gm=mp.mpf(g_value);qm=mp.mpf(q_value);ep=gm**2/(24*mp.pi**2)
        gravity=gm**2/(mp.pi*qm**2);km=1/(96*mp.pi**2);xs=qm**2/2
        def mass_rhs(r):
            invariant=qm**2/(2*r**4)
            density_m=invariant/(4*gm**2)-km*invariant*mp.log(invariant/xs)
            return 4*mp.pi*gravity*r*r*density_m
        for radius in radii:
            rm=mp.mpf(radius)
            quadrature_mass=mp.mpf('.5')+mp.quad(mass_rhs,[1,rm])
            analytic_mass=1+2*ep-(mp.mpf('.5')+2*ep*(mp.log(rm)+1))/rm
            integral_metric=1-2*quadrature_mass/rm
            value=float(metric_from_horizon_offset(radius-1,epsilon_value))
            mass_error=float(abs(quadrature_mass-analytic_mass));mass_errors.append(mass_error)
            metric_error=float(abs(mp.mpf(value)-integral_metric));metric_errors.append(metric_error)
            norm=(1+4*epsilon_value*math.log(radius))/radius**4
            table.append(dict(radius_over_re=radius,metric_f=value,rn_metric=(1-1/radius)**2,
                geometric_mass_over_re=float(analytic_mass),mass_correction_over_rn=float(analytic_mass-(1-1/(2*rm))),
                rho_times_8piG_re2=norm,pr_times_8piG_re2=-norm,
                pperp_times_8piG_re2=(1+4*epsilon_value*math.log(radius)-2*epsilon_value)/radius**4,
                radial_nec_times_8piG_re2=0.,trace_times_8piG_re2=-4*epsilon_value/radius**4,
                mass_quadrature_absolute_error=mass_error,metric_vs_mass_quadrature_absolute_error=metric_error))
        for offset in (1e-24,1e-18,1e-15,1e-12,1e-9,1e-6,.001,.01,.1):
            u=mp.mpf(offset)
            reference=(u*u+4*ep*(mp.log1p(u)-u))/(1+u)**2
            safe=float(metric_from_horizon_offset(offset,epsilon_value))
            relative=float(abs(mp.mpf(safe)/reference-1))
            radius=1+offset
            naive=1-2*(1+2*epsilon_value)/radius+(1+4*epsilon_value*(math.log(radius)+1))/radius**2
            near.append(dict(horizon_offset=offset,safe_metric=safe,high_precision_metric=float(reference),
                safe_relative_error=relative,naive_expanded_metric=naive,
                radius_addition_rounds_to_horizon=radius==1.))
    assert max(mass_errors)<1e-60 and max(metric_errors)<5e-14
    assert max(row['safe_relative_error'] for row in near)<5e-13
    assert float(metric_from_horizon_offset(0.,epsilon_value))==0
    weak_offsets=np.array([0.,1e-12,.01,.5,1.,9.,99.])
    weak_error=float(np.max(abs(metric_from_horizon_offset(weak_offsets,0.)-
                                (weak_offsets/(1+weak_offsets))**2)))
    assert weak_error==0.
    trace_table=[]
    for radius in (1.,1.5,2.,10.,100.):
        substitution={R:sy.Rational(str(radius)),q:q_value}
        def number(expression):return float(expression.subs(substitution).evalf(40))
        trace_table.append(dict(radius_over_re=radius,
            gauge_trace_times_re4=number(gauge_trace),gravitational_trace_times_re4=number(gravity_trace),
            lowest_landau_trace_times_re4=number(lll_trace),full_required_trace_times_re4=number(gauge_trace+gravity_trace),
            local_log_plus_lll_trace_times_re4=number(gauge_trace+lll_trace),
            required_additional_trace_times_re4=number(required_trace)))
    return dict(model='Local one-loop magnetic log action; separate classical-RN four-dimensional trace audit',
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        parameters=dict(magnetic_charge_q=q_value,gauge_coupling_at_Xstar=g_value,
            epsilon=epsilon_value,k=1/(96*math.pi**2),finite_collar_radius_interval=[1.,100.],
            re_over_sqrtG=math.sqrt(math.pi)*q_value/g_value,
            subtraction='Xstar=q^2/(2re^4), tied to the stated magnetic reference scale; g is the renormalized coupling at that reference.'),
        equations=dict(conventions='Signature -+++; X=F_ab F^ab=q^2/(2r^4); re^2=pi G q^2/g^2; one massless unit-charge Dirac species.',
            local_action='L=-X/(4g^2)+X*log(X/Xstar)/(96pi^2); stress=-4L_X F_ac F_b^c+g_ab L.',
            stresses='rho=-L; pr=L; pperp=L-2X L_X; rho+pr=0; trace=4(L-X L_X)=-X/(24pi^2).',
            mass='R=r/re; mu=Gm/re; mu prime=(1/2+2epsilon log R)/R^2; mu(1)=1/2.',
            metric='f=1-2(1+2epsilon)/R+[1+4epsilon(log R+1)]/R^2; f(1)=f_prime(1)=0, f_second(1)=2(1-2epsilon).',
            stable_metric='u=R-1; f=u^2/(1+u)^2*[1+4epsilon*(log1p(u)-u)/u^2], with a local series for the ratio.',
            product='On the branch connected to classical RN, RS=re and LA^2/RS^2=1/(1-2epsilon). This is a product solution of the stated local action.',
            near_horizon_clock='f approximately (r-RS)^2/LA^2; tau=t/ell and rho_th=ell*(r-RS)/LA^2. The ideal global product lapse is N=LA/(ell*cos chi), so its full null optical crossing remains pi*ell. This is a local calibration statement, not a recomputed matched collar.',
            flare_out='For ds^2=-N(l)^2dt^2+dl^2+R(l)^2dOmega^2, (R_prime/N)_prime=-4piG*R*(rho+pr)/N. A regular N>0 center with R_prime=0 and this sector alone forces constant R, not strict opening.'),
        product_solution=dict(sphere_radius_over_re=1.,ads_radius_squared_over_sphere_squared=1/(1-2*epsilon_value),
            ads_radius_over_re=1/math.sqrt(1-2*epsilon_value),horizon_second_derivative=2*(1-2*epsilon_value),
            formal_extrapolated_mass_coefficient=1+2*epsilon_value,
            formal_extrapolated_mass_shift=2*epsilon_value,
            finite_collar_R100_mass_correction=table[-1]['mass_correction_over_rn'],
            maximum_loop_density_fraction_on_collar=4*epsilon_value*math.log(100)),
        finite_collar_table=table,near_horizon_checks=near,
        classical_rn_trace_audit=dict(background='CLASSICAL extremal RN f=(1-1/R)^2; traces in re^-4 units. This is not the backreacted log-action metric.',
            invariants={name:str(expression) for name,expression in dict(Ricci_scalar=ricci_scalar,Ricci_squared=ricci_squared,Kretschmann=kretschmann,Weyl_squared=weyl,Euler_density=euler).items()},
            equations=dict(gravity='(C^2/20-11E4/360)/(16pi^2)=(42R^2-84R+53)/(720pi^2 R^8); box R4=0.',
                gauge='-q^2/(48pi^2 R^4)',lll='q(4R-6)/(96pi^2 R^6)',
                remainder='Full required trace minus [local gauge-log trace+LLL trace] = gravitational trace-LLL trace.',
                horizon_remainder='q/(48pi^2)+11/(720pi^2)'),
            table=trace_table,
            scope='The nonzero remainder is a necessary TRACE constraint only. It does not supply a conserved source tensor, quantum state or nonlocal determinant. The local log action and existing LLL contribution cannot simply be added and called the complete four-dimensional source; renormalization and double counting remain unresolved.'),
        verification=dict(all_declared_checks_pass=True,symbolic_residuals={name:str(value) for name,value in residuals.items()},
            quadrature_precision_digits=70,maximum_mass_quadrature_error=max(mass_errors),
            maximum_metric_vs_quadrature_error=max(metric_errors),
            maximum_near_horizon_relative_error=max(row['safe_relative_error'] for row in near),
            weak_coupling_metric_error=weak_error,
            thresholds=dict(mass_quadrature_absolute=1e-60,metric_quadrature_absolute=5e-14,near_horizon_relative=5e-13)),
        scope='Conditional local leading-log benchmark, in a specified finite-subtraction convention. The radial NEC vanishes, so this magnetic sector alone cannot open a regular static throat; this does not rule out combinations with a nonlocal Casimir source. The finite collar is R=1..100. The algebraic mass coefficient at infinity is a formal extrapolation of this local action, not a matched physical ADM mass or verified asymptotic quantum state. Higher derivative, nonlocal, gravitational-anomaly and state-dependent terms are not solved by this action. No loop term should be counted again in a separately renormalized source. No stability, source completion, preparation, material or human passage is established.',
        sources=[dict(url='https://arxiv.org/pdf/1808.03489',equation='4.20 and Appendix A',supports='Mostly-plus Dirac gauge trace -F^2/(24pi^2), fixing the local log coefficient within the ansatz.'),
            dict(url='https://arxiv.org/pdf/0707.3751',equation='2.1-2.3',supports='Four-dimensional curvature anomaly coefficients and the distinction between the trace and the full state-dependent stress tensor.'),
            dict(url='https://arxiv.org/pdf/1911.02271',equation='35',supports='Independent Dirac gravitational trace written with Riemann/Ricci invariants.'),
            dict(url='https://arxiv.org/html/1807.04726v3',equation='2.1-2.3',supports='Magnetic charge and gauge coupling normalization; this code derives the stated local-action solution, not a new MMP solution.')])


if __name__=='__main__':
    result=magnetic_polarization_checks()
    project_path('magnetic_polarization_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(product=result['product_solution'],verification=result['verification'],
                         trace_table=result['classical_rn_trace_audit']['table']),indent=2))
