"""Compact throat variation probing an energy-plus-radial-pressure contrast.

At fixed sigma, N_eta=N exp(eta*f), A_eta=A exp(-eta*f), and r=b is fixed.
The radial Weyl frequency/mass coefficients change by exp(-2eta*f) and
exp(-eta*f). Compact support leaves the endpoint bases and exterior map
unchanged. The result is a smeared difference of rho+p_r, conditional on
the same continuum/state assumptions as increments94/95, not absolute NEC.
A second compact stationary-shift variation retains the endpoint time
offset and probes minus the lower orthonormal radial energy flux. The two
variations give opposite radial-null PHASE DIFFERENCES, not absolute NEC.

Run: .venv\\Scripts\\python.exe -B -W error throat_null_response.py
"""
from project_paths import project_path
from functools import lru_cache
import hashlib
import json
import math
import operator
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad,solve_ivp

from curved_throat_modes import throat_scattering,_energy,_positive
from curved_exterior_closure import raw_mouth_map
from exterior_boundary_map import _complex_json
from exterior_tail_transport import facet_map
from imaginary_weyl_leads import tail_map
from local_lapse_response import LocalLapse
from phase_vacuum_contrast import PhaseContrast


PROBES=(dict(name='central',center_fraction=0.,width_fraction=.65),
        dict(name='upper_side',center_fraction=.60,width_fraction=.25))


def bump(sigma,cut,index):
    probe=PROBES[index]
    u=((sigma/cut-probe['center_fraction'])/probe['width_fraction'])**2
    return math.exp(1-1/(1-u)) if u<1 else 0.


def probe_geometry(ell,rho_c,b=.2):
    cut=math.atan(rho_c)
    result=[]
    for index,probe in enumerate(PROBES):
        center=cut*probe['center_fraction'];width=cut*probe['width_fraction']
        integral=quad(lambda s:bump(s,cut,index),center-width,center+width,
                      epsabs=2e-14,epsrel=2e-13)[0]
        weighted=quad(lambda s:bump(s,cut,index)/math.cos(s)**2,center-width,center+width,
                      epsabs=2e-14,epsrel=2e-13)[0]
        points=np.linspace(-cut,cut,161)
        result.append(dict(**probe,center=center,width=width,support=[center-width,center+width],
            endpoint_clearance=min(center-width+cut,cut-center-width),
            I_f=ell*integral,full_geometric_weight=4*math.pi*b**4/ell*weighted,
            sigma=points.tolist(),profile=[bump(float(s),cut,index) for s in points],
            formula='u=((sigma-center)/width)^2; f=exp(1-1/(1-u)) for u<1, otherwise0.',
            normalization='Peak-one weight; no extra angular degeneracy or division by geometric volume is applied to the sector trace.'))
    return result


@lru_cache(maxsize=2048)
def mode_variation(energy,n,ell=.2*math.sqrt(5),rho_c=2.,q=1,eta=0.,selected=0,
                   rtol=3e-11,atol=3e-13):
    """Riccati scattering with shared first variations for the two probes.

    A finite eta solve uses only selected probe; tangents are computed at that
    background. Production derivatives use eta=0, so both share one baseline.
    S maps [plus_left,minus_right] to [minus_left,plus_right].
    """
    energy=_energy(energy);ell=_positive(ell,'ell');rho_c=_positive(rho_c,'rho_c')
    n=operator.index(n);q=operator.index(q);selected=operator.index(selected)
    eta=float(eta)
    if n<0 or q<=0 or not 0<=selected<len(PROBES) or not math.isfinite(eta) or abs(eta)>.1:
        raise ValueError('Require n>=0,q>0, a known probe and a small finite |eta|<=.1')
    cut=math.atan(rho_c);lam=math.sqrt(n*(n+q))
    if not 0<rtol<1 or not 0<atol<1:raise ValueError('Positive ODE tolerances below1 required')

    def rhs(sigma,state):
        right,logt,left=state[:3]
        tangents=state[3:].reshape(len(PROBES),3)
        f=np.array([bump(sigma,cut,j) for j in range(len(PROBES))])
        e=energy*ell*math.exp(-2*eta*f[selected])
        mass=lam/math.cos(sigma)*math.exp(-eta*f[selected])
        tt=np.exp(2*logt)
        base=np.array([2j*e*right-1j*mass*(1+right*right),
                       1j*e-1j*mass*right,-1j*mass*tt])
        dr,dt,dl=tangents.T
        derivatives=np.column_stack(((2j*e-2j*mass*right)*dr-4j*f*e*right+1j*f*mass*(1+right*right),
            -2j*f*e+1j*f*mass*right-1j*mass*dr,
            1j*f*mass*tt-2j*mass*tt*dt))
        return np.r_[base,derivatives.reshape(-1)]
    solved=solve_ivp(rhs,(-cut,cut),np.zeros(3*(len(PROBES)+1),complex),method='DOP853',
        rtol=rtol,atol=atol,max_step=2*cut/24)
    if not solved.success:raise RuntimeError('Metric-variation throat solve failed: '+solved.message)
    right,logt,left=solved.y[:3,-1]
    tangents=solved.y[3:,-1].reshape(len(PROBES),3)
    transmission=np.exp(logt)
    S=np.array([[left,transmission],[transmission,right]],complex)
    derivatives=np.array([[[dl,transmission*dt],[transmission*dt,dr]] for dr,dt,dl in tangents])
    if not np.all(np.isfinite(S)) or not np.all(np.isfinite(derivatives)):
        raise ArithmeticError('Nonfinite throat scattering or metric derivative')
    return S,derivatives,dict(n=n,energy=_complex_json(energy),ell=ell,rho_c=rho_c,eta=eta,
        selected=selected,log_transmission=_complex_json(logt),
        log_transmission_derivatives=_complex_json(tangents[:,1]),nfev=solved.nfev,
        optical_length=2*ell*cut,rtol=rtol,atol=atol,
        transmission_underflowed=bool(transmission==0),
        normalization='Canonical radial current per common coordinate time; the n=0 paired matrix is only a control, one allowed chirality is retained.')


def raw_variation(kappa,nmax,rho_c=2.,eta=0.,selected=0,rtol=3e-11,atol=3e-13):
    """Phase-free canonical B and its derivatives in frozen93 mouth ordering."""
    ell=.2*math.hypot(1.,rho_c);count=2*nmax+1
    B=np.zeros((count,count),complex);derivatives=np.zeros((len(PROBES),count,count),complex)
    records=[]
    for n in range(nmax+1):
        S,dS,diag=mode_variation(1j*kappa,n,ell,rho_c,1,eta,selected,rtol,atol)
        if n==0:
            B[0,nmax]=S[1,0];derivatives[:,0,nmax]=dS[:,1,0]
        else:
            j=n-1;upper=n;lower=nmax+n
            B[upper,j]=S[1,1];B[upper,lower]=S[1,0]
            B[lower,j]=S[0,1];B[lower,lower]=S[0,0]
            derivatives[:,upper,j]=dS[:,1,1];derivatives[:,upper,lower]=dS[:,1,0]
            derivatives[:,lower,j]=dS[:,0,1];derivatives[:,lower,lower]=dS[:,0,0]
        records.append(dict(n=n,S=_complex_json(S),derivatives=_complex_json(dS),diagnostics=diag))
    return B,derivatives,records


def phase_map(matrix,phase,nmax):
    incoming=np.r_[np.full(nmax+1,np.exp(1j*phase)),np.ones(nmax)]
    outgoing=np.r_[np.full(nmax,np.exp(1j*phase)),np.ones(nmax+1)]
    return incoming[:,None]*matrix*outgoing.conj()[None,:]


def shift_derivatives(raw,kappa,nmax,probes):
    """Fixed-clock shift v=eta*(N/A)*f: forward +kappa*I, reverse -kappa*I."""
    derivatives=np.zeros((len(PROBES),len(raw),len(raw)),complex)
    for j,probe in enumerate(probes):
        factor=kappa*probe['I_f']
        derivatives[j,:nmax+1,nmax:]=factor*raw[:nmax+1,nmax:]
        derivatives[j,nmax+1:,:nmax]=-factor*raw[nmax+1:,:nmax]
    return derivatives


def finite_shift_mode(energy,n,ell,rho_c,eta,selected,metric_variation=False,shift_multiplier=1.):
    """Independent finite-shift solve after time diagonalization, with endpoint offset."""
    cut=math.atan(rho_c);lam=math.sqrt(n*(n+1))
    if abs(eta*shift_multiplier)>=1:raise ValueError('The stationary shift must be subluminal')
    def rhs(sigma,state):
        right,logt,left,offset=state
        f=bump(sigma,cut,selected);speed=eta*shift_multiplier*f;denominator=1-speed*speed
        deformation=eta*f if metric_variation else 0.
        e=energy*ell*math.exp(-2*deformation)/denominator
        mass=lam/math.cos(sigma)*math.exp(-deformation)/math.sqrt(denominator)
        return np.array([2j*e*right-1j*mass*(1+right*right),1j*e-1j*mass*right,
                         -1j*mass*np.exp(2*logt),ell*math.exp(-2*deformation)*speed/denominator],complex)
    solved=solve_ivp(rhs,(-cut,cut),np.zeros(4,complex),method='DOP853',
                     rtol=3e-13,atol=3e-15,max_step=2*cut/48)
    if not solved.success:raise RuntimeError(solved.message)
    right,logt,left,offset=solved.y[:,-1]
    # The common endpoint clock makes the two off-diagonal factors opposite.
    return np.array([[left,np.exp(logt+1j*energy*offset)],
                     [np.exp(logt-1j*energy*offset),right]],complex)


class ThroatNull(LocalLapse):
    """Reuse only frozen LocalLapse._green; no exterior lapse is changed."""

    def __init__(self,h=.05,R=2.,nmax=8,L=40.,rho_c=2.):
        self.base=PhaseContrast(h,R,nmax,L,rho_c)
        self.active=self.base.model.active
        self.config=self.base.config.copy()
        self.probes=probe_geometry(.2*math.hypot(1.,rho_c),rho_c)
        self.cache={}

    def evaluate(self,kappa,audit=False):
        kappa=_positive(kappa,'kappa');key=(kappa,bool(audit))
        if key in self.cache:return self.cache[key]
        started=time.perf_counter();base=self.base;nmax=self.config['nmax']
        matrix=(base.unshifted+kappa*base.mass)[self.active][:,self.active].tocsc()
        rawZ,tail=tail_map(kappa,self.config['R'],self.config['L'],nmax=nmax)
        Z,_,geometry=facet_map(base.outer,base.outer_basis,rawZ)
        green=self._green(matrix,Z);S=green['S']
        raw,d_raw,mode_records=raw_variation(kappa,nmax,self.config['rho_c'])
        d_shift=shift_derivatives(raw,kappa,nmax,self.probes)
        Bs=[];dBs=[];shiftBs=[];matrices=[];logs=[];conditions=[];minimum=[]
        for phase in (0.,math.pi):
            b=phase_map(raw,phase,nmax);db=np.array([phase_map(d,phase,nmax) for d in d_raw])
            B=np.linalg.solve(base.Fout.T,(base.Fin@b).T).T
            dB=np.array([np.linalg.solve(base.Fout.T,(base.Fin@d).T).T for d in db])
            shiftB=np.array([np.linalg.solve(base.Fout.T,(base.Fin@phase_map(d,phase,nmax)).T).T for d in d_shift])
            M=np.eye(len(B))-S@B
            sign,logabs=np.linalg.slogdet(M);sv=np.linalg.svd(M,compute_uv=False)
            if sign==0 or not math.isfinite(logabs):raise ArithmeticError('Singular sampled sewing')
            Bs.append(B);dBs.append(dB);shiftBs.append(shiftB);matrices.append(M);logs.append(float(logabs))
            minimum.append(float(sv[-1]));conditions.append(float(sv[0]/sv[-1]))
        traces=[];shift_traces=[]
        for j in range(len(PROBES)):
            traces.append(np.trace(np.linalg.solve(matrices[1],S@dBs[1][j])
                                   -np.linalg.solve(matrices[0],S@dBs[0][j])))
            shift_traces.append(np.trace(np.linalg.solve(matrices[1],S@shiftBs[1][j])
                                   -np.linalg.solve(matrices[0],S@shiftBs[0][j])))
        values=[float(t.real/math.pi) for t in traces]
        shift_values=[float(t.real/math.pi) for t in shift_traces]
        row=dict(kappa=kappa,probes=[dict(name=p['name'],trace=_complex_json(t),response_integrand=v,
            shift_trace=_complex_json(st),shift_response_integrand=sv,
            null_plus_integrand=v-2*sv,null_minus_integrand=v+2*sv)
            for p,t,v,st,sv in zip(PROBES,traces,values,shift_traces,shift_values)],log_determinants=logs,
            minimum_sewing_singular_values=minimum,sewing_condition_numbers=conditions,
            outer_geometry=geometry,tail_accepted_steps=tail['accepted_steps'])
        if audit:
            baseline_errors=[];base.problem['kappa']=kappa
            for phase,B in zip((0.,math.pi),Bs):
                old,_,_=raw_mouth_map(base.problem,phase,self.config['rho_c'],'curved')
                oldfacet=np.linalg.solve(base.Fout.T,(base.Fin@old).T).T
                baseline_errors.append(float(np.max(abs(oldfacet-B))))
            fd=[]
            for j,p in enumerate(PROBES):
                estimates=[]
                for eta in (.02,.01,.005):
                    contrasts=[];bs=[]
                    for signed in (-eta,eta):
                        trial,_,_=raw_variation(kappa,nmax,self.config['rho_c'],signed,j)
                        triallogs=[];trialbs=[]
                        for phase in (0.,math.pi):
                            B=np.linalg.solve(base.Fout.T,(base.Fin@phase_map(trial,phase,nmax)).T).T
                            sign,logabs=np.linalg.slogdet(np.eye(len(B))-S@B)
                            if sign==0:raise ArithmeticError('Singular finite-metric sewing')
                            triallogs.append(float(logabs));trialbs.append(B)
                        contrasts.append(triallogs[1]-triallogs[0]);bs.append(trialbs)
                    derivative=-(contrasts[1]-contrasts[0])/(2*eta*math.pi)
                    derivative_maps=[(bs[1][phase]-bs[0][phase])/(2*eta) for phase in range(2)]
                    estimates.append(dict(eta=eta,response_integrand_derivative=derivative,
                        absolute_error=abs(derivative-values[j]),relative_error=abs(derivative-values[j])/max(abs(values[j]),1e-15),
                        matrix_derivative_relative_error=max(float(np.linalg.norm(derivative_maps[phase]-dBs[phase][j])
                            /max(np.linalg.norm(dBs[phase][j]),1e-15)) for phase in range(2))))
                richardson=(4*estimates[-1]['response_integrand_derivative']-estimates[-2]['response_integrand_derivative'])/3
                fd.append(dict(name=p['name'],analytic_response=values[j],finite_differences=estimates,
                    richardson_response=richardson,richardson_relative_error=abs(richardson-values[j])/max(abs(values[j]),1e-15),
                    successive_absolute_error_ratios=[estimates[k]['absolute_error']/max(estimates[k+1]['absolute_error'],1e-30) for k in (0,1)]))
            row['audit']=dict(baseline_frozen93_matrix_errors=baseline_errors,finite_differences=fd,
                matrices=dict(S_exterior=_complex_json(S),B_zero=_complex_json(Bs[0]),B_pi=_complex_json(Bs[1]),
                    derivatives_zero=_complex_json(dBs[0]),derivatives_pi=_complex_json(dBs[1]),
                    shift_derivatives_zero=_complex_json(shiftBs[0]),shift_derivatives_pi=_complex_json(shiftBs[1]),
                    F_in=_complex_json(base.Fin),F_out=_complex_json(base.Fout),B_raw=_complex_json(raw),
                    derivatives_raw=_complex_json(d_raw),shift_derivatives_raw=_complex_json(d_shift)),modes=mode_records)
        row['seconds']=time.perf_counter()-started
        self.cache[key]=row;return row


def independent_transfer_tangent_control():
    """Low-action 2x2 transfer derivative, independently audited by formation_updates."""
    ell,rho_c,n=.7,.8,1
    cut=np.arctan(rho_c);lam=np.sqrt(n*(n+1))
    s3=np.diag([1.,-1.]);s2=np.array([[0.,-1j],[1j,0.]])
    rows=[]
    for energy in (.4,.3j):
        def rhs(s,y):
            T,U,V=y.reshape(3,2,2)
            e=energy*ell;m=lam/np.cos(s);Q=1j*e*s3+m*s2
            tangent=-2j*e*s3-m*s2
            return np.array([Q@T,Q@U+bump(s,cut,0)*tangent@T,
                             Q@V+bump(s,cut,1)*tangent@T]).ravel()
        y0=np.zeros((3,2,2),complex);y0[0]=np.eye(2)
        sol=solve_ivp(rhs,(-cut,cut),y0.ravel(),method='DOP853',
                      rtol=3e-13,atol=3e-15,max_step=.02)
        if not sol.success:raise RuntimeError(sol.message)
        T,U,V=sol.y[:,-1].reshape(3,2,2)
        a,b,c,d=T.ravel();determinant=a*d-b*c
        direct=np.array([[-c/d,1/d],[determinant/d,b/d]])
        ds=[];det_tangents=[]
        for Z in (U,V):
            aa,bb,cc,dd=Z.ravel();ddeter=aa*d+a*dd-bb*c-b*cc
            ds.append(np.array([[-(cc*d-c*dd)/d**2,-dd/d**2],
                [(ddeter*d-determinant*dd)/d**2,(bb*d-b*dd)/d**2]]))
            det_tangents.append(abs(ddeter))
        S,dS,_=mode_variation(energy,n,ell,rho_c,rtol=3e-13,atol=3e-15)
        row=dict(energy=_complex_json(energy),scattering_error=float(np.max(abs(S-direct))),
            tangent_error=float(np.max(abs(dS-np.array(ds)))),
            determinant_error=float(abs(determinant-1)),
            determinant_tangent_error=float(max(det_tangents)))
        if complex(energy).imag==0:
            row['unitary_tangent_error']=float(max(np.max(abs(Z.conj().T@S+S.conj().T@Z)) for Z in dS))
        rows.append(row)
    return dict(ell=ell,rho_c=rho_c,n=n,controls=rows)


def analytic_controls():
    """Actual zero-mode ODE versus optical length, plus an independent circle integral."""
    ell=.2*math.sqrt(5);rho_c=2.;cut=math.atan(rho_c)
    probes=probe_geometry(ell,rho_c);zero=[];baseline=[]
    for energy in (.4,.4j,2j):
        S,dS,diag=mode_variation(energy,0)
        expected=np.exp(1j*energy*2*ell*cut)
        for j,probe in enumerate(probes):
            tangent=-2j*energy*probe['I_f']*expected
            eta=.02
            direct,_,_=mode_variation(energy,0,eta=eta,selected=j)
            length=ell*quad(lambda s:math.exp(-2*eta*bump(s,cut,j)),-cut,cut,
                            points=probe['support'],epsabs=2e-13,epsrel=2e-13)[0]
            zero.append(dict(energy=_complex_json(energy),probe=probe['name'],
                transmission_error=float(abs(S[1,0]-expected)),
                tangent_error=float(abs(dS[j,1,0]-tangent)),
                finite_metric_transmission_error=float(abs(direct[1,0]-np.exp(1j*energy*length))),
                I_f=probe['I_f'],log_transmission_derivative=_complex_json(-2j*energy*probe['I_f'])))
    for energy in (.4,.4j,2j):
        for n in (0,1,3,8):
            S,dS,_=mode_variation(energy,n)
            old,diag=throat_scattering(energy,n,ell,rho_c)
            baseline.append(dict(energy=_complex_json(energy),n=n,
                scattering_error=float(np.max(abs(S-old))),
                smallest_imaginary_passivity_eigenvalue=(float(np.min(np.linalg.eigvalsh(np.eye(2)-S.conj().T@S)))
                    if complex(energy).imag>0 else None),
                real_unitarity_tangent_error=(float(max(np.max(abs(Z.conj().T@S+S.conj().T@Z)) for Z in dS))
                    if complex(energy).imag==0 else None)))
    shift=[]
    for energy in (.4,.3j):
        for n in (0,1):
            S,_,_=mode_variation(energy,n,.7,.8,rtol=3e-13,atol=3e-15)
            for j,probe in enumerate(probe_geometry(.7,.8)):
                I=probe['I_f'];expected=np.array([[0.,1j*energy*I*S[0,1]],[-1j*energy*I*S[1,0],0.]])
                estimates=[]
                for eta in (.02,.01,.005):
                    plus=finite_shift_mode(energy,n,.7,.8,eta,j)
                    minus=finite_shift_mode(energy,n,.7,.8,-eta,j)
                    derivative=(plus-minus)/(2*eta)
                    estimates.append(derivative)
                richardson=(4*estimates[-1]-estimates[-2])/3
                shift.append(dict(energy=_complex_json(energy),n=n,probe=probe['name'],
                    finite_difference_errors=[float(np.max(abs(d-expected))) for d in estimates],
                    richardson_error=float(np.max(abs(richardson-expected))),
                    derivative=_complex_json(expected),endpoint_clock_offset_derivative=I,
                    minimum_stationary_lapse_factor=math.sqrt(1-.02**2)))
    combined=[]
    energy=.3j;n=1;ell_control=.7;rho_control=.8;index=0
    S,dS,_=mode_variation(energy,n,ell_control,rho_control,rtol=3e-13,atol=3e-15)
    I=probe_geometry(ell_control,rho_control)[index]['I_f']
    shift_matrix=np.array([[0.,1j*energy*I*S[0,1]],[-1j*energy*I*S[1,0],0.]])
    for multiplier in (-2.,2.):
        expected=dS[index]+multiplier*shift_matrix;estimates=[]
        for eta in (.02,.01,.005):
            plus=finite_shift_mode(energy,n,ell_control,rho_control,eta,index,True,multiplier)
            minus=finite_shift_mode(energy,n,ell_control,rho_control,-eta,index,True,multiplier)
            estimates.append((plus-minus)/(2*eta))
        richardson=(4*estimates[-1]-estimates[-2])/3
        combined.append(dict(shift_multiplier=multiplier,energy=_complex_json(energy),n=n,probe='central',
            finite_difference_errors=[float(np.max(abs(d-expected))) for d in estimates],
            richardson_error=float(np.max(abs(richardson-expected))),
            scope='One actual combined metric path: N_eta=N exp(eta f), A_eta=A exp(-eta f), v_eta=multiplier*eta*(N_eta/A_eta)*f.'))
    circle=[]
    with mp.workdps(45):
        length=mp.mpf('1.7')
        for p in probes:
            I=mp.mpf(str(p['I_f']))
            numeric=mp.quad(lambda k:-2*k*I/(mp.pi*mp.sinh(k*length)) if k else -2*I/(mp.pi*length),
                            [0,1,4,mp.inf])
            exact=-mp.pi*I/(2*length**2)
            circle.append(dict(probe=p['name'],total_length=float(length),I_f=float(I),
                imaginary_axis_integral=str(numeric),exact=str(exact),absolute_error=float(abs(numeric-exact)),
                exact_shift=str(exact/2),exact_null_plus='0',exact_null_minus=str(2*exact),
                energy_contrast='E_pi-E_0=-pi/(4L), one complex chiral sector with symmetric sea.',
                integrand='-2*kappa*I_f/[pi*sinh(kappa*L)]',
                length_derivative='dL/deta=-2I_f; dDeltaE/deta=-pi*I_f/(2L^2).'))
    return dict(zero_mode=zero,frozen93_baseline=baseline,closed_circle=circle,finite_shift=shift,combined_metric=combined,
                independent_transfer=independent_transfer_tangent_control())


def integrate_response(evaluator,order=8,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
    nodes,weights=np.polynomial.legendre.leggauss(order);edges=np.asarray(edges,float)
    if len(edges)<2 or edges[0]<=0 or not np.all(np.diff(edges)>0):
        raise ValueError('Positive increasing frequency edges required')
    panels=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        frequencies=(lo+hi)/2+(hi-lo)/2*nodes
        rows=[evaluator.evaluate(float(k)) for k in frequencies]
        values=np.array([[p['response_integrand'] for p in r['probes']] for r in rows])
        shifts=np.array([[p['shift_response_integrand'] for p in r['probes']] for r in rows])
        panels.append(dict(lower=float(lo),upper=float(hi),kappa=frequencies.tolist(),
            response_integrands=values.tolist(),responses=((hi-lo)/2*(weights@values)).tolist(),
            shift_response_integrands=shifts.tolist(),shift_responses=((hi-lo)/2*(weights@shifts)).tolist(),
            null_plus_integrands=(values-2*shifts).tolist(),null_minus_integrands=(values+2*shifts).tolist(),
            null_plus_responses=((hi-lo)/2*(weights@(values-2*shifts))).tolist(),
            null_minus_responses=((hi-lo)/2*(weights@(values+2*shifts))).tolist(),
            minimum_sewing_singular_value=min(min(r['minimum_sewing_singular_values']) for r in rows)))
    return dict(config=evaluator.config.copy(),order=int(order),edges=edges.tolist(),
        probe_order=[p['name'] for p in PROBES],panels=panels,
        total_responses=np.sum([p['responses'] for p in panels],axis=0).tolist(),
        total_shift_responses=np.sum([p['shift_responses'] for p in panels],axis=0).tolist(),
        total_null_plus_responses=np.sum([p['null_plus_responses'] for p in panels],axis=0).tolist(),
        total_null_minus_responses=np.sum([p['null_minus_responses'] for p in panels],axis=0).tolist(),
        scope='Smeared energy-plus-radial-pressure PHASE DIFFERENCE over the displayed frequency window only.')


def _decode(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def checks():
    started=time.perf_counter();root=Path(__file__).resolve().parent
    inputs=('local_lapse_response.py','local_lapse_response_results.json',
        'phase_vacuum_contrast.py','phase_vacuum_contrast_results.json',
        'curved_exterior_closure.py','curved_throat_modes.py','imaginary_exterior.py',
        'imaginary_weyl_leads.py','exterior_boundary_map.py','exterior_tail_transport.py',
        'magnetic_tail_map.py','spherical_weyl_leads.py','two_mouth_mesh.py')
    hashes={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs}
    evaluators={};integrations={}
    def get(h=.05,R=2.,N=8,L=40.):
        key=(h,R,N,L)
        if key not in evaluators:evaluators[key]=ThroatNull(h,R,N,L)
        return evaluators[key]
    def integral(h=.05,order=8,edges=(1e-6,.002,.02,.2,1.,4.,8.,12.)):
        key=(h,order,tuple(edges))
        if key not in integrations:
            stamp=time.perf_counter();record=integrate_response(get(h=h),order,edges)
            record['seconds']=time.perf_counter()-stamp;integrations[key]=record
            print(json.dumps(dict(h=h,order=order,edges=edges,responses=record['total_responses'],
                shift=record['total_shift_responses'],null_plus=record['total_null_plus_responses'],
                null_minus=record['total_null_minus_responses'],seconds=record['seconds'])),flush=True)
        return integrations[key]

    analytic=analytic_controls();primary=get()
    coarse=integral();fine=integral(order=16)
    mesh=[integral(h=h) for h in (.08,.05,.035)]
    fd_cases=[primary.evaluate(k,audit=True) for k in (.01,.4,2.)]
    prior=json.loads((project_path('local_lapse_response_results.json')).read_text(encoding='utf-8'))
    exterior_preservation=[]
    for now,old in zip(fd_cases,prior['finite_difference_checks']):
        if now['kappa']!=old['kappa']:raise AssertionError('Prior comparison frequencies changed')
        exterior_preservation.append(dict(kappa=now['kappa'],maximum_matrix_difference=float(np.max(abs(
            _decode(now['audit']['matrices']['S_exterior'])-_decode(old['audit']['matrices']['S_exterior']))))))
    frequencies=(1e-4,.01,.4,2.,8.)
    reference=[primary.evaluate(k) for k in frequencies];point_controls=[]
    for name,config in (('angular',dict(N=12)),('matching_radius',dict(R=4.)),('far_tail',dict(L=80.))):
        evaluator=get(**config);rows=[evaluator.evaluate(k) for k in frequencies]
        differences=[[a['response_integrand']-b['response_integrand'] for a,b in zip(row['probes'],ref['probes'])]
                     for row,ref in zip(rows,reference)]
        extra={kind:[[a[field]-b[field] for a,b in zip(row['probes'],ref['probes'])] for row,ref in zip(rows,reference)]
               for kind,field in (('shift','shift_response_integrand'),('null_plus','null_plus_integrand'),('null_minus','null_minus_integrand'))}
        point_controls.append(dict(name=name,config=evaluator.config,frequencies=list(frequencies),samples=rows,
            response_differences=differences,maximum_absolute_difference=float(np.max(abs(np.array(differences)))),
            other_response_differences=extra,other_maximum_absolute_differences={k:float(np.max(abs(np.array(v)))) for k,v in extra.items()},
            scope='Selected-frequency sensitivity only; not an integrated error bound.'))
    ode_controls=[]
    for row in fd_cases:
        k=row['kappa'];S=_decode(row['audit']['matrices']['S_exterior'])
        B,dB,_=raw_variation(k,8,rtol=3e-12,atol=3e-14)
        ts=[]
        for phase in (0.,math.pi):
            bf=np.linalg.solve(primary.base.Fout.T,(primary.base.Fin@phase_map(B,phase,8)).T).T
            dfs=[np.linalg.solve(primary.base.Fout.T,(primary.base.Fin@phase_map(d,phase,8)).T).T for d in dB]
            M=np.eye(len(B))-S@bf
            ts.append([np.trace(np.linalg.solve(M,S@df)) for df in dfs])
        response=(np.asarray(ts[1])-np.asarray(ts[0])).real/math.pi
        old=np.array([r['response_integrand'] for r in row['probes']])
        ode_controls.append(dict(kappa=k,refined_responses=response.tolist(),absolute_differences=abs(response-old).tolist(),
            rtol=3e-12,atol=3e-14))
    ir=integral(edges=(1e-8,1e-7,1e-6));uv=integral(edges=(12.,16.))
    ir_mesh=[dict(h=h,samples=[get(h=h).evaluate(k) for k in (1e-8,1e-6,1e-4,.001,.01)]) for h in (.08,.05,.035)]
    all_rows=[r for e in evaluators.values() for r in e.cache.values()]
    fd=[r for row in fd_cases for r in row['audit']['finite_differences']]
    quad_error=np.array(fine['total_responses'])-np.array(coarse['total_responses'])
    mesh_changes=[(np.array(b['total_responses'])-np.array(a['total_responses'])).tolist() for a,b in zip(mesh[:-1],mesh[1:])]
    extra_quad={kind:(np.array(fine[key])-np.array(coarse[key])).tolist() for kind,key in (
        ('shift','total_shift_responses'),('null_plus','total_null_plus_responses'),('null_minus','total_null_minus_responses'))}
    extra_mesh={kind:[(np.array(b[key])-np.array(a[key])).tolist() for a,b in zip(mesh[:-1],mesh[1:])] for kind,key in (
        ('shift','total_shift_responses'),('null_plus','total_null_plus_responses'),('null_minus','total_null_minus_responses'))}
    transfer=analytic['independent_transfer']['controls']
    tests=dict(
        compact_positive_metric=min(p['endpoint_clearance'] for p in primary.probes)>0 and math.exp(-.02)>0,
        actual_zero_mode_optical_identity=max(max(r['transmission_error'],r['tangent_error'],r['finite_metric_transmission_error']) for r in analytic['zero_mode'])<2e-10,
        independent_closed_circle=max(r['absolute_error'] for r in analytic['closed_circle'])<1e-35,
        independent_full_transfer=max(max(r['scattering_error'],r['tangent_error'],r['determinant_error'],r['determinant_tangent_error']) for r in transfer)<1e-10,
        real_energy_unitary_tangent=max(r.get('unitary_tangent_error',0) for r in transfer)<1e-10,
        independent_finite_shift=max(r['richardson_error'] for r in analytic['finite_shift'])<2e-8,
        independent_combined_metric=max(r['richardson_error'] for r in analytic['combined_metric'])<2e-8,
        frozen93_baseline=max(r['scattering_error'] for r in analytic['frozen93_baseline'])<2e-8
            and max(max(r['audit']['baseline_frozen93_matrix_errors']) for r in fd_cases)<2e-8,
        fixed_exterior_matches95=max(r['maximum_matrix_difference'] for r in exterior_preservation)<1e-12,
        direct_metric_determinant_Richardson=max(r['richardson_relative_error'] for r in fd)<1e-6,
        direct_metric_second_order=min(min(r['successive_absolute_error_ratios']) for r in fd)>3.7
            and max(max(r['successive_absolute_error_ratios']) for r in fd)<4.3,
        refined_throat_ODE=max(max(r['absolute_differences']) for r in ode_controls)<2e-8,
        finite_frequency_quadrature=max(max(abs(quad_error)),max(abs(v) for a in extra_quad.values() for v in a))<2e-7,
        sampled_sewings_nonsingular=min(min(r['minimum_sewing_singular_values']) for r in all_rows)>1e-7,
        frozen_dependencies_preserved=hashes=={name:hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs})
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,probes=primary.probes,
        parameters=dict(q=1,mu=.5,b=.2,d=1.,rho_c=2.,ell=.2*math.sqrt(5),phase_difference='pi minus0',
            units='hbar=c=d=1; eta and f dimensionless; integrated responses restore hbar*c/d.',
            source='N_eta=N*exp(eta*f), A_eta=A*exp(-eta*f), fixed sphere radius and magnetic connection.',
            shift_source='Separate v_eta=eta*(N/A)*f in ds^2=-N^2 dt^2+A^2(d_sigma-v_eta dt)^2+r^2 dOmega^2; endpoint clocks are fixed.',
            baseline='N=b*sec(sigma)/ell, A=b*sec(sigma), r=b, |sigma|<=atan(rho_c).'),
        analytic_controls=analytic,quadrature_refinement=[coarse,fine],mesh_refinement=mesh,
        representative=fd_cases[1],finite_difference_checks=fd_cases,exterior_preservation=exterior_preservation,
        selected_frequency_reference=reference,selected_frequency_sensitivities=point_controls,
        ODE_refinement=ode_controls,infrared_extension=ir,ultraviolet_extension=uv,infrared_mesh_samples=ir_mesh,
        diagnostics=dict(finite_window_responses=fine['total_responses'],frequency_quadrature_differences=quad_error.tolist(),
            finite_window_shift_responses=fine['total_shift_responses'],finite_window_null_plus_responses=fine['total_null_plus_responses'],
            finite_window_null_minus_responses=fine['total_null_minus_responses'],other_frequency_quadrature_differences=extra_quad,
            other_mesh_successive_changes=extra_mesh,
            mesh_successive_changes=mesh_changes,maximum_Richardson_relative_error=max(r['richardson_relative_error'] for r in fd),
            minimum_sampled_sewing_singular_value=min(min(r['minimum_sewing_singular_values']) for r in all_rows),
            maximum_sampled_sewing_condition_number=max(max(r['sewing_condition_numbers']) for r in all_rows),
            frequency_geometry_evaluations=len(all_rows),combined_continuum_bound=None,omitted_IR_bound=None,
            omitted_UV_bound=None,integrated_N_R_L_sensitivity_bound=None),
        equations=dict(radial='d_sigma f=(i E ell exp(-2eta*f_probe) sigma3 + lambda sec(sigma) exp(-eta*f_probe) sigma2) f.',
            tangent='Derivative includes BOTH frequency -2f and mass -f variations, with fixed endpoint frames.',
            sewing='M_phi=I-S_ext B_phi; B_facet=F_in B_raw F_out^(-1). S_ext and F_in/F_out stay fixed.',
            response='dDeltaE/deta=(1/pi) integral ReTr[M_pi^(-1)S_ext dB_pi - M_0^(-1)S_ext dB_0] dkappa.',
            continuum='Conditionally: integral d_sigma dOmega N A r^2 f_probe Delta(rho+p_r). This is half the sum of the two normalized opposite-null contractions.',
            shift='D_shift=-integral N A r^2 f_probe Delta(T_lower_hat0hatr); forward cross dB=+kappa I_f B, reverse=-kappa I_f B, same-end reflections=0.',
            endpoint_clock='T=t+F, DeltaF=integral ell*eta*f/(1-eta^2*f^2) d_sigma. Diagonal Riccati coefficients e=Eell/(1-eta^2*f^2), m=lambda sec(sigma)/sqrt(1-eta^2*f^2) have no linear shift change; opposite transmission factors exp(-/+iE DeltaF) must remain.',
            null_contrasts='Q_plus=Q_rho_plus_p-2D_shift; Q_minus=Q_rho_plus_p+2D_shift, for orthonormal k_plus/minus=(1,+/-1). Absolute errors matter near cancellation.',
            zero_mode='d(log transmission)/deta=2*kappa*I_f, I_f=ell*integral f_probe d_sigma.',
            unchanged_geometry='For a fixed-radius spherical cylinder, geometric rho+p_r=0 for every positive N and A. This phase contrast is not its absolute supporting stress.'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_source_sha256=hashes,
        seconds=time.perf_counter()-started,
        scope=['Smeared energy-plus-radial-pressure phase-contrast diagnostic in a retained q=1, mu=.5 dissipative DG sector on a prescribed synthetic geometry.',
            'Both phases have the same local metric deformation, action and endpoint data. Only the assumed cross-throat holonomy differs; local determinant factors cancel in the difference.',
            'Only one allowed zero-mode chirality is retained. No extra azimuthal, angular, or Dirac degeneracy is multiplied into the trace.',
            'The additional fixed-clock compact-shift derivative supplies minus the lower orthonormal radial energy-flux PHASE DIFFERENCE; its combination yields both radial-null contrasts, still not absolute NEC.',
            'A nonzero phase contrast cannot be called support for the unchanged constant-radius cylinder, whose geometric rho+p_r remains zero.',
            'The synthetic cylinder-to-flat interface and endpoint spin/gauge identification are assumed; the necessary junction source and physical collar remain unsolved.',
            'The local stress interpretation requires continuum self-adjoint sewing, a common regulator/state prescription and convergence; the current dissipative finite-element diagnostic does not prove these.',
            'Responses are not divided by a geometric volume. Profile weights are supplied for reproducibility, not as a density estimator.',
            'Sampled N/R/L and ODE sensitivities, frequency quadrature and mesh changes are distinct. They do not furnish a combined continuum or infinite-frequency bound.',
            'The finite-mesh IR plateau is not proof of a physical spectral gap. Uncomputed tails, absolute stress and the full field/sector sum remain unresolved.',
            'No physical material, prepared state, self-consistent Einstein solution, stable passage or time-travel result is established.'],
        sources=['https://arxiv.org/pdf/1106.0173','https://arxiv.org/pdf/1606.07834',
                 'https://arxiv.org/pdf/0711.2395','https://arxiv.org/html/1807.04726v3'])
    return result


if __name__=='__main__':
    result=checks()
    target=project_path('throat_null_response_results.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(target.read_text(encoding='utf-8'))==result
    print(json.dumps({k:result[k] for k in ('passed','tests','diagnostics','seconds','source_sha256')},indent=2))
    if not result['passed']:raise SystemExit('A compact throat metric-variation control failed; diagnostics retained.')
