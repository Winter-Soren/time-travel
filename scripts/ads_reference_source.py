"""AdS-cap reference for a compact product-throat source difference.

The cap/reference pilot reuses saved increment98 exterior matrices. The
frequency integration computes fresh physical exterior Green responses on
the same declared geometry.
Both Weyl components are retained.  The massless reflecting reference couples
their chiral zero modes; a separate reflecting zero for each Weyl is invalid.
The output is a retained-sector, selected-frequency reference-relative source
diagnostic, not an angular/frequency-complete absolute stress or a collar.

Primary radial domain/reference: https://arxiv.org/pdf/2208.08252, Eqs21--25,
SecIII.C--D.  The invariant massless vacuum has half-integer frequencies;
an arbitrary self-adjoint boundary phase need not define that vacuum.

Run: .venv\\Scripts\\python.exe -B -W error ads_reference_source.py
"""
from project_paths import project_path
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import block_diag
from numpy.polynomial.legendre import leggauss

from exterior_boundary_map import _complex_json
from throat_null_response import probe_geometry


def _array(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


@lru_cache(maxsize=4096)
def cap_reflection(energy,lam,ell,rho_c,start=1e-5,rtol=3e-12,atol=3e-14):
    """Right normalizable cap f_minus/f_plus in the sigma1 mass convention.

    delta=pi/2-sigma and r=-i*y.  The regular Riccati equation is
    y_delta=lambda*csc(delta)*(1-y*y)-2*k*y, k=-i*E*ell, y(0)=1.
    For E=i*kappa>0 it preserves 0<y<1.  The left cap has the same
    reflection in the right-Weyl pair convention.  The left-Weyl cap is
    its imaginary-frequency adjoint after swapping the channel spaces.
    """
    energy=complex(energy);lam=float(lam);ell=float(ell);rho_c=float(rho_c)
    if (not all(math.isfinite(v) for v in (energy.real,energy.imag,lam,ell,rho_c,start,rtol,atol))
        or lam<0 or ell<=0 or rho_c<=0 or energy.imag<0 or not 0<rtol<1 or not 0<atol<1):
        raise ValueError('Finite energy in the closed upper half-plane and positive geometry/tolerances required')
    end=math.pi/2-math.atan(rho_c)
    if not 0<start<end:raise ValueError('The Frobenius start must lie inside the cap')
    if lam==0:
        value=-1j*np.exp(2j*energy*ell*end)
        return value,dict(nfev=0,delta_start=0.,delta_end=end,exact_massless=True)
    if lam<.5:raise ValueError('The unique normalizable massive domain is used only for lambda>=1/2')
    k=-1j*energy*ell;d=1+2*lam
    a1=-2*k/d;a2=2*k*k/d**2
    a3=(2*lam*k/(3*d)-4*k**3/d**3)/(3+2*lam)
    initial=1+a1*start+a2*start**2+a3*start**3
    def rhs(delta,y):return lam/math.sin(delta)*(1-y*y)-2*k*y
    solved=solve_ivp(rhs,(start,end),np.array([initial],complex),method='DOP853',rtol=rtol,atol=atol)
    if not solved.success:raise RuntimeError('AdS cap integration failed: '+solved.message)
    value=-1j*solved.y[0,-1]
    if not np.isfinite(value):raise ArithmeticError('Nonfinite cap reflection')
    return value,dict(nfev=solved.nfev,delta_start=start,delta_end=end,exact_massless=False,
        canonical_falloff='f~delta^lambda*(1,-i), physical spinor~delta^(lambda+1/2)',
        boundary_series_order=3,reflection_modulus=float(abs(value)))


def reference_map(kappa,nmax,mu,rho_c=2.):
    """Combined R/L cap map from core-exiting to core-entering channels.

    Each chirality uses the saved exterior ordering.  For mu=1/2, R columns
    are [upper0,upper n,lower n], R rows [upper n,lower0,lower n]; L swaps
    these spaces.  The two -i zero reflections connect different Weyls.
    """
    if kappa<=0 or not math.isfinite(float(kappa)):raise ValueError('Positive finite kappa required')
    if mu not in (-.5,.5,1.5) or int(nmax)!=nmax or nmax<1:raise ValueError('Unsupported retained sector')
    ell=.2*math.hypot(1.,rho_c);zero=mu==.5;count=2*nmax+int(zero)
    right=np.zeros((count,count),complex);caps=[]
    for n in range(1,nmax+1):
        r,diag=cap_reflection(1j*kappa,math.sqrt(n*(n+1)),ell,rho_c)
        if zero:upper_row,upper_col,lower=n-1,n,nmax+n
        else:upper_row=upper_col=n-1;lower=nmax+n-1
        right[upper_row,upper_col]=r;right[lower,lower]=r
        caps.append(dict(n=n,reflection=_complex_json(r),**diag))
    full=block_diag(right,right.conj().T)
    if zero:
        r,_=cap_reflection(1j*kappa,0.,ell,rho_c)
        full[count+0,0]=r
        full[nmax,count+nmax]=r
    return full,dict(ell=ell,half_cap_optical_length=ell*(math.pi/2-math.atan(rho_c)),
        global_ads_optical_length=math.pi*ell,has_cross_chirality_zero=zero,caps=caps)


def _raw_environment(matrices):
    values={name:_array(value) for name,value in matrices.items()}
    right=np.linalg.solve(values['F_mouth_right_out'],values['S_right']@values['F_mouth_right_in'])
    left=np.linalg.solve(values['F_mouth_left_out'],values['S_left_physical']@values['F_mouth_left_in'])
    return block_diag(right,left),values


def reference_relative(Sconnected,Sreference,B,pressure,shift):
    """Cancellation-safe compact metric derivative of a common-core ratio.

    (I-Sc B)^-1 Sc-(I-Sr B)^-1 Sr
       =(I-Sc B)^-1 (Sc-Sr) (I-B Sr)^-1.
    This is ordinary complex linear algebra, not an adjoint identity.
    """
    eye=np.eye(len(B));Mc=eye-Sconnected@B;Mr=eye-Sreference@B
    difference=np.linalg.solve((eye-B@Sreference).T,
        np.linalg.solve(Mc,Sconnected-Sreference).T).T
    direct=np.linalg.solve(Mc,Sconnected)-np.linalg.solve(Mr,Sreference)
    p=np.array([np.trace(difference@d).real/math.pi for d in pressure])
    s=np.array([np.trace(difference@d).real/math.pi for d in shift])
    pc=np.array([np.trace(np.linalg.solve(Mc,Sconnected@d)).real/math.pi for d in pressure])
    pr=np.array([np.trace(np.linalg.solve(Mr,Sreference@d)).real/math.pi for d in pressure])
    sc=np.array([np.trace(np.linalg.solve(Mc,Sconnected@d)).real/math.pi for d in shift])
    sr=np.array([np.trace(np.linalg.solve(Mr,Sreference@d)).real/math.pi for d in shift])
    scv=np.linalg.svd(Mc,compute_uv=False);srv=np.linalg.svd(Mr,compute_uv=False)
    return dict(pressure_integrands=p.tolist(),shift_integrands=s.tolist(),
        null_plus_integrands=(p-2*s).tolist(),null_minus_integrands=(p+2*s).tolist(),
        connected_pressure_trace=pc.tolist(),reference_pressure_trace=pr.tolist(),
        connected_shift_trace=sc.tolist(),reference_shift_trace=sr.tolist(),
        minimum_singular_values=[float(scv[-1]),float(srv[-1])],
        condition_numbers=[float(scv[0]/scv[-1]),float(srv[0]/srv[-1])],
        resolvent_difference_identity_error=float(np.max(abs(difference-direct))),
        connected_logabs=float(np.linalg.slogdet(Mc)[1]),reference_logabs=float(np.linalg.slogdet(Mr)[1]))


def saved_case(config,row,save_matrices=False):
    """No DG rebuild: transform frozen98 actual physical R/L maps to raw bases."""
    connected,values=_raw_environment(row['matrices'])
    reference,info=reference_map(row['kappa'],config['nmax'],config['mu'],config['rho_c'])
    branches=[];matrices=[]
    for index,phase in enumerate((0.,math.pi)):
        raw=values['B_raw_right'][index];dp=values['pressure_raw_right'][index];ds=values['shift_raw_right'][index]
        B=block_diag(raw,raw.conj().T)
        P=[block_diag(d,d.conj().T) for d in dp]
        Q=[block_diag(d,-d.conj().T) for d in ds]
        result=reference_relative(connected,reference,B,P,Q)
        branches.append(dict(phase=phase,**result))
        if save_matrices:matrices.append(dict(B=_complex_json(B),pressure=_complex_json(P),shift=_complex_json(Q)))
    output=dict(kappa=row['kappa'],phases=branches,cap_diagnostics=info,
        reference_phase_pressure_error=float(np.max(abs(np.array(branches[0]['reference_pressure_trace'])-branches[1]['reference_pressure_trace']))),
        reference_phase_shift_error=float(np.max(abs(np.array(branches[0]['reference_shift_trace'])-branches[1]['reference_shift_trace']))))
    if save_matrices:output['matrices']=dict(S_connected=_complex_json(connected),S_reference=_complex_json(reference),phases=matrices)
    return output


def cap_checks():
    ell=.2*math.sqrt(5);rho_c=2.;cut=math.atan(rho_c);rows=[]
    sigma3=np.diag([1.,-1.]);sigma2=np.array([[0.,-1j],[1j,0.]])
    for n in (1,2,8):
        lam=math.sqrt(n*(n+1))
        energy=(lam+.5)/ell
        ground,_=cap_reflection(energy,lam,ell,rho_c)
        for kappa in (.01,.4,2.):
            value,info=cap_reflection(1j*kappa,lam,ell,rho_c)
            fine,_=cap_reflection(1j*kappa,lam,ell,rho_c,start=5e-6,rtol=5e-13,atol=5e-15)
            # Independent normalizable-column equation after removing delta^lambda.
            start=1e-5;end=info['delta_end'];k=kappa*ell;d=1+2*lam
            y=1-2*k/d*start+2*k*k/d**2*start**2+(2*lam*k/(3*d)-4*k**3/d**3)/(3+2*lam)*start**3
            def rhs(delta,f):
                return (-1j*(1j*kappa)*ell*sigma3-lam/math.sin(delta)*sigma2-lam/delta*np.eye(2))@f
            solved=solve_ivp(rhs,(start,end),np.array([1.,-1j*y]),method='DOP853',rtol=5e-13,atol=5e-15)
            if not solved.success:raise RuntimeError(solved.message)
            direct=solved.y[1,-1]/solved.y[0,-1]
            rows.append(dict(n=n,kappa=kappa,reflection=_complex_json(value),
                cutoff_tolerance_change=float(abs(value-fine)),direct_column_error=float(abs(value-direct)),
                ground_mode_error=float(abs(ground-np.exp(-1j*cut))),
                positive_real_y=bool(0<(-value.imag)<1),real_part_abs=float(abs(value.real))))
    sigma1=np.array([[0.,1.],[1.,0.]])
    U=(np.eye(2)+1j*sigma1)/math.sqrt(2)
    representation=dict(
        kinetic_rotation_error=float(np.max(abs(U@sigma3@U.conj().T-sigma2))),
        mass_rotation_error=float(np.max(abs(U@sigma1@U.conj().T-sigma1))),
        right_normalizable_error=float(np.max(abs(U@np.array([1.,-1j])-np.array([math.sqrt(2),0.])))),
        left_normalizable_error=float(np.max(abs(U@np.array([1.,1j])-np.array([0.,1j*math.sqrt(2)])))),
        formula='U=(I+i*sigma1)/sqrt2: U H U_dagger=-i*sigma2*partial_sigma+lambda*sec(sigma)*sigma1, the primary paper Eq25',
        domain='Right normalizable vector maps to (sqrt2,0), left to (0,i*sqrt2): type-I endpoints Phi2(right)=0, Phi1(left)=0, including the selected massless limit.')
    return dict(rows=rows,representation_domain=representation,exact_ground='(cos sigma)^lambda*(exp(i sigma/2),exp(-i sigma/2)); E*ell=lambda+1/2',
        start_scope='Observed Frobenius-cutoff/tolerance refinement, not interval-certified truncation error.')


def massless_checks():
    ell=.2*math.sqrt(5);rho_c=2.;cut=math.atan(rho_c);core=2*ell*cut
    total=math.pi*ell;half=ell*(math.pi/2-cut);If=probe_geometry(ell,rho_c)[0]['I_f']
    rows=[]
    for kappa in (.01,.4,2.,12.):
        r=-1j*math.exp(-2*kappa*half);t=math.exp(-kappa*core)
        S=np.array([[0.,r],[r,0.]],complex);B=t*np.eye(2);dB=2*kappa*If*B
        M=np.eye(2)-S@B;det=np.linalg.det(M)
        expected=1+math.exp(-2*kappa*total)
        derivative=np.trace(np.linalg.solve(M,S@dB)).real/math.pi
        exact=-4*kappa*If/(math.pi*(math.exp(2*kappa*total)+1))
        rows.append(dict(kappa=kappa,determinant_error=float(abs(det-expected)),
            pressure_derivative_error=abs(derivative-exact),pressure_integrand=derivative))
    def pressure(k):
        e=math.exp(-2*k*total);return -4*k*If*e/(math.pi*(1+e))
    energy,energy_error=quad(lambda k:-math.log1p(math.exp(-2*k*total))/math.pi,0,np.inf,epsabs=2e-13)
    response,response_error=quad(pressure,0,np.inf,epsabs=2e-13)
    loop=2*total
    return dict(global_ads_optical_length=total,I_f=If,rows=rows,
        reference_energy=energy,reference_energy_exact=-math.pi/(24*total),energy_quadrature_estimate=energy_error,
        reference_pressure=response,reference_pressure_exact=-math.pi*If/(12*total**2),pressure_quadrature_estimate=response_error,
        connected_NS_Dirac_length=loop,connected_NS_Dirac_pressure=-math.pi*If/(3*loop**2),
        matching_local_pressure_error=abs(response+math.pi*If/(3*loop**2)),
        scope='Full 1+1 Dirac zero sector with two reflecting endpoints, equivalent to one chiral NS circle of length 2*L_AdS; no extra doubled chirality.')


def high_frequency_checks():
    """Three authorized new exterior samples to test the OLD cutoff, not a tail bound."""
    from azimuthal_source import AzimuthalSector
    started=time.perf_counter();case=AzimuthalSector(.5,h=.05,R=2.,nmax=8,L=40.,rho_c=2.)
    rows=[]
    for kappa in (12.,24.,48.):
        row=case.evaluate(kappa,save_matrices=True,adjoint_control=False)
        result=saved_case(case.config,row,save_matrices=kappa==24.)
        result['exterior_seconds']=row['seconds']
        result['exterior_equation_relative_residuals']=row['exterior_equation_relative_residuals']
        rows.append(result)
        print(json.dumps(dict(high_frequency=kappa,primary_pi=result['phases'][1]['pressure_integrands'],
            shift=result['phases'][1]['shift_integrands'],exterior_seconds=row['seconds'])),flush=True)
    return dict(config=case.config,samples=rows,seconds=time.perf_counter()-started,
        scope='Selected-frequency cutoff exploration only; three decaying samples do not certify the omitted integral or angular tail.')


FIELDS=('pressure','shift','null_plus','null_minus')
PHASES=('0','pi')


def integrate_reference(case,order,edges):
    """One bounded quadrature, both phase settings sharing each exterior solve."""
    started=time.perf_counter();nodes,weights=leggauss(order);panels=[]
    minimum=np.full((2,2),np.inf);maximum_condition=np.zeros((2,2))
    identity_error=equation_error=reference_flux=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        ks=(lo+hi)/2+(hi-lo)*nodes/2;ws=(hi-lo)*weights/2
        values={phase:{field:[] for field in FIELDS} for phase in PHASES}
        singular=[];conditions=[]
        for kappa in ks:
            raw=case.evaluate(float(kappa),save_matrices=True,adjoint_control=False)
            row=saved_case(case.config,raw)
            case.cache.clear() # keep the current finite-volume model, not every large matrix record
            equation_error=max(equation_error,max(raw['exterior_equation_relative_residuals']))
            singular.append([p['minimum_singular_values'] for p in row['phases']])
            conditions.append([p['condition_numbers'] for p in row['phases']])
            minimum=np.minimum(minimum,singular[-1]);maximum_condition=np.maximum(maximum_condition,conditions[-1])
            for label,phase in zip(PHASES,row['phases']):
                for field in FIELDS:values[label][field].append(phase[field+'_integrands'])
                identity_error=max(identity_error,phase['resolvent_difference_identity_error'])
                reference_flux=max(reference_flux,max(abs(v) for v in phase['reference_shift_trace']))
        responses={phase:{field:np.einsum('i,ij->j',ws,values[phase][field]).tolist() for field in FIELDS} for phase in PHASES}
        panels.append(dict(bounds=[lo,hi],kappa=ks.tolist(),weights=ws.tolist(),integrands=values,
            responses=responses,minimum_singular_values=singular,condition_numbers=conditions))
        print(json.dumps(dict(integration_h=case.config['h'],order=order,panel=[lo,hi],
            phase_pi_pressure=responses['pi']['pressure'],elapsed=time.perf_counter()-started)),flush=True)
    totals={phase:{field:np.sum([p['responses'][phase][field] for p in panels],axis=0).tolist() for field in FIELDS} for phase in PHASES}
    return dict(config=case.config,order=order,edges=list(edges),panels=panels,totals=totals,
        diagnostics=dict(minimum_singular_values=minimum.tolist(),maximum_condition_numbers=maximum_condition.tolist(),
            maximum_equation_relative_residual=equation_error,resolvent_identity_error=identity_error,
            reference_shift_trace_error=reference_flux),seconds=time.perf_counter()-started)


def integration_checks(saved,high):
    from azimuthal_source import AzimuthalSector
    edges=(1e-6,.02,.2,1.,4.,12.,24.,48.,96.)
    started=time.perf_counter();case=AzimuthalSector(.5,h=.05,R=2.,nmax=8,L=40.,rho_c=2.)
    coarse=integrate_reference(case,8,edges);fine=integrate_reference(case,12,edges)
    small=AzimuthalSector(.5,h=.035,R=2.,nmax=8,L=40.,rho_c=2.)
    finer=integrate_reference(small,8,edges)
    def difference(a,b):
        return {phase:{field:(np.asarray(b['totals'][phase][field])-a['totals'][phase][field]).tolist()
                       for field in FIELDS} for phase in PHASES}
    envelope={phase:{field:np.sum([abs(np.asarray(b['responses'][phase][field])-a['responses'][phase][field])
                  for a,b in zip(coarse['panels'],fine['panels'])],axis=0).tolist()
                  for field in FIELDS} for phase in PHASES}
    larger=AzimuthalSector(.5,h=.05,R=2.,nmax=12,L=40.,rho_c=2.);angular=[]
    old2=next(row for c in saved['mesh_refinement'] if c['config']['mu']==.5 and c['config']['h']==.05
              for row in c['samples'] if row['kappa']==2.)
    old2_config=next(c['config'] for c in saved['mesh_refinement'] if c['config']['mu']==.5 and c['config']['h']==.05)
    bases={2.:saved_case(old2_config,old2),12.:next(r for r in high['samples'] if r['kappa']==12.)}
    for kappa in (2.,12.):
        row=saved_case(larger.config,larger.evaluate(kappa,save_matrices=True,adjoint_control=False),save_matrices=False)
        angular.append(dict(kappa=kappa,nmax_from=8,nmax_to=12,baseline=bases[kappa],larger=row,
            differences={phase:{field:(np.asarray(row['phases'][j][field+'_integrands'])-
                bases[kappa]['phases'][j][field+'_integrands']).tolist() for field in FIELDS} for j,phase in enumerate(PHASES)}))
        print(json.dumps(dict(angular_kappa=kappa,phase_pi_pressure=row['phases'][1]['pressure_integrands'])),flush=True)
    return dict(primary='fine',quadrature_absolute_target=2e-6,coarse=coarse,fine=fine,finer_mesh=finer,
        comparisons=dict(quadrature_difference=difference(coarse,fine),quadrature_panel_envelope=envelope,
            matched_order_mesh_difference=difference(coarse,finer),
            final_panel_contribution=fine['panels'][-1]['responses'],
            final_two_panel_contribution={phase:{field:(np.asarray(fine['panels'][-2]['responses'][phase][field])+
                fine['panels'][-1]['responses'][phase][field]).tolist() for field in FIELDS} for phase in PHASES}),
        selected_angular=angular,seconds=time.perf_counter()-started,
        scope='One retained mu=1/2 sector including its actual two Weyl responses. Both physical phases use the same local geometry. Fixed finite frequency window, observed quadrature/mesh differences, and two angular-cutoff samples; no certified omitted tail or full four-dimensional source.')


def checks():
    started=time.perf_counter();caps=cap_checks();zero=massless_checks()
    path=project_path('azimuthal_source_results.json');saved=json.loads(path.read_text());records=[]
    for case in saved['mesh_refinement']:
        config=case['config']
        rows=[saved_case(config,row,save_matrices=(config['h']==.05 and row['kappa']==.4)) for row in case['samples']]
        records.append(dict(config=config,samples=rows))
        print(json.dumps(dict(mu=config['mu'],h=config['h'],samples=[dict(kappa=r['kappa'],primary_pi=r['phases'][1]['pressure_integrands'],
            reference=r['phases'][1]['reference_pressure_trace']) for r in rows])),flush=True)
    high=high_frequency_checks()
    integration=integration_checks(saved,high)
    allrows=[r for c in records for r in c['samples']]+high['samples']
    tests=dict(exact_ground_mode=max(r['ground_mode_error'] for r in caps['rows'])<2e-10,
        primary_representation_domain=max(v for k,v in caps['representation_domain'].items() if k.endswith('_error'))<1e-14,
        cap_cutoff_tolerance_refinement=max(r['cutoff_tolerance_change'] for r in caps['rows'])<2e-10,
        independent_normalizable_column=max(r['direct_column_error'] for r in caps['rows'])<2e-10,
        imaginary_cap_contractivity=all(r['positive_real_y'] and r['real_part_abs']<1e-12 for r in caps['rows']),
        massless_reflecting_determinant=max(r['determinant_error'] for r in zero['rows'])<1e-13,
        massless_compact_derivative=max(r['pressure_derivative_error'] for r in zero['rows'])<1e-13,
        massless_normalization=abs(zero['reference_energy']-zero['reference_energy_exact'])<1e-12,
        massless_source_integral=abs(zero['reference_pressure']-zero['reference_pressure_exact'])<1e-12,
        matched_zero_local_source=zero['matching_local_pressure_error']<1e-12,
        common_core_resolvent_identity=max(p['resolvent_difference_identity_error'] for r in allrows for p in r['phases'])<2e-9,
        reference_holonomy_independent=max(r['reference_phase_pressure_error'] for r in allrows)<2e-10,
        reference_shift_zero=max(abs(v) for r in allrows for p in r['phases'] for v in p['reference_shift_trace'])<2e-10,
        sampled_invertibility=min(p['minimum_singular_values'][1] for r in allrows for p in r['phases'])>1e-7,
        high_frequency_equation_residual=max(max(r['exterior_equation_relative_residuals']) for r in high['samples'])<1e-8)
    integrated=[integration[name] for name in ('coarse','fine','finer_mesh')]
    tests.update(integrated_equation_residual=max(r['diagnostics']['maximum_equation_relative_residual'] for r in integrated)<1e-8,
        integrated_sewing_invertibility=min(v for r in integrated for pair in r['diagnostics']['minimum_singular_values'] for v in pair)>1e-8,
        integrated_resolvent_identity=max(r['diagnostics']['resolvent_identity_error'] for r in integrated)<2e-9,
        integrated_reference_shift_zero=max(r['diagnostics']['reference_shift_trace_error'] for r in integrated)<2e-9,
        observed_quadrature_panel_envelope=max(v for ph in integration['comparisons']['quadrature_panel_envelope'].values()
            for values in ph.values() for v in values)<integration['quadrature_absolute_target'],
        integrated_null_identities=max(abs(np.asarray(r['totals'][phase]['null_plus'])-
            np.asarray(r['totals'][phase]['pressure'])+2*np.asarray(r['totals'][phase]['shift'])).max()
            for r in integrated for phase in PHASES)<1e-12)
    tests={name:bool(value) for name,value in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,stage='AdS-reference compact source difference: single-sector finite-window integration with pilot and observed controls',
        parameters=dict(q=1,mu=[-.5,.5,1.5],b=.2,d=1.,rho_c=2.,primary_phase=math.pi,connection='explicit reflected spin/gauge endpoint branch'),
        cap_checks=caps,massless_checks=zero,mesh_refinement=records,high_frequency=high,integration=integration,probes=saved['probes'],
        observable_definitions=dict(pressure='+(1/pi) Re Tr[((I-Sc B)^-1 Sc-(I-Sref B)^-1 Sref) dB_pressure]',
            shift='Same trace with dB_shift; equals minus the weighted lower orthonormal radial-flux difference',
            null_plus='pressure-2*shift',null_minus='pressure+2*shift',
            reference_zero='The two -i cap reflections couple physical opposite Weyl zero modes; D0=1+exp(-2*kappa*pi*ell)',
            local_source='Integral of N*A*r^2*f*(T_connected-T_invariant_reference), not divided by probe volume'),
        scope=['The pilot reuses frozen98 matrices; the integrated calculation computes fresh actual right/left physical exterior responses on the same declared geometry, retaining finite-facet mismatch and applying no silent doubling.',
            'The standard invariant AdS2 x S2 reference has radial null contraction zero by its AdS2 symmetry; using the difference as an absolute local source is conditional on the correctly matched Hadamard states, full angular/frequency completion and a common covariant renormalization.',
            'The integrated calculation retains only mu=1/2 and n<=8, with selected n<=12 comparisons. It is not a complete absolute stress, an NEC conclusion, a matched exterior, a smooth collar or an Einstein solution.',
            'The compact probe vanishes before both cuts, so the cap/exterior maps and canonical boundary bases are fixed; common core determinant factors cancel between environments. Individual connected/reference traces are environmental determinant contributions, not separately renormalized local stresses.',
            'The finite product throat still requires its unprovided physical junction and source. The reference does not create that junction.',
            'The integral covers kappa in [1e-6,96]; the last panels diagnose the inadequacy of the old cutoff12, but do not certify either omitted frequency tail or the missing angular modes.'],
        sources=[dict(url='https://arxiv.org/pdf/2208.08252',use='Canonical Dirac radial operator, unique normalizable massive domain, and the half-integer invariant massless vacuum; Eqs21--25 and SecIII.C--D.')],
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in
            set(saved['input_source_sha256'])|{'azimuthal_source_results.json','azimuthal_source.py','azimuthal_modes.py',
                'curved_throat_modes.py','throat_null_response.py','dirac_sewing_audit.py'}},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.perf_counter()-started)
    return result


if __name__=='__main__':
    result=checks();path=project_path('ads_reference_source_results.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    assert json.loads(path.read_text())==result
    print(json.dumps(dict(passed=result['passed'],tests=result['tests'],seconds=result['seconds'])),flush=True)
    if not result['passed']:raise SystemExit('A recorded numerical gate failed; evidence remains in the results JSON.')
