"""Numerical and exact-special-function checks of scalar fields on an Ori core.

These are test-field equations restricted to an EXACT invariant mode family in
four dimensions, not nonlinear Einstein-scalar simulations. No integration is
continued through the chronology boundary.
"""
from __future__ import annotations
import json, math, platform, time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import mpmath as mp

OUT=Path(__file__).resolve().parent
mp.mp.dps=60
S_END=12*math.log(10)
SAMPLE=np.linspace(0,S_END,601)
A_VALUES=(0.,.001,.01,.02,.04,.06,.08,.10,.12)
K_VALUES=(1,2,4,8)
M_VALUES=(0.,.1,1.,10.)
RTOL=5e-11
ATOL=5e-13

def exponent(a: float)->float:
    if not 0 <= a < .125:
        raise ValueError('This audit uses 0 <= a < 1/8.')
    return (math.sqrt(1-8*a)+math.sqrt(1+8*a)-2)/4

def rhs(s: float, Y: np.ndarray, p: complex, m: float)->np.ndarray:
    u,w=Y
    return np.array([w,-p*w-m*m*math.exp(-s)*u],dtype=complex)

def energy(s: np.ndarray, u: np.ndarray,w: np.ndarray,k:int,m:float)->np.ndarray:
    """T_ab u^a u^b for two real canonical scalars along a regular observer."""
    T=-np.exp(-s)
    dT=-w/T
    dz=1j*k*u
    along_u=(dT+dz)/np.sqrt(2+T)
    along_n=((T+1)*dT-dz)/np.sqrt(2+T)
    return .5*(abs(along_u)**2+abs(along_n)**2+m*m*abs(u)**2)

def exact_massive(a:float,k:int,m:float,s_values):
    """Independent solution using 0F1 and a two-by-two initial-data match."""
    aa=mp.mpf(str(a)); mm=mp.mpf(str(m)); kk=mp.mpf(k)
    real_p=(mp.sqrt(1-8*aa)+mp.sqrt(1+8*aa)-2)/4
    pp=real_p+2j*kk
    def basis(S):
        z=mm*mm*S
        f1=mp.hyp0f1(1-pp,z)
        f1p=mm*mm/(1-pp)*mp.hyp0f1(2-pp,z)
        pref=mp.exp(pp*mp.log(-S))
        f2=pref*mp.hyp0f1(1+pp,z)
        f2p=pp/S*f2+pref*mm*mm/(1+pp)*mp.hyp0f1(2+pp,z)
        return f1,f2,-S*f1p,-S*f2p
    f1,f2,w1,w2=basis(mp.mpf(-1))
    coeff=mp.lu_solve(mp.matrix([[f1,f2],[w1,w2]]),mp.matrix([1,-pp]))
    vals=[]
    for sv in s_values:
        S=-mp.exp(-mp.mpf(str(float(sv))))
        f1,f2,w1,w2=basis(S)
        vals.append([complex(coeff[0]*f1+coeff[1]*f2),
                     complex(coeff[0]*w1+coeff[1]*w2)])
    return np.array(vals).T, float(abs(coeff[1]))

start=time.perf_counter()
cases=[]
trajectories={}
for a in A_VALUES:
    for k in K_VALUES:
        for m in M_VALUES:
            p=exponent(a)+2j*k
            sol=solve_ivp(rhs,(0,S_END),np.array([1,-p],complex),args=(p,m),
                          method='DOP853',rtol=RTOL,atol=ATOL,t_eval=SAMPLE)
            if not sol.success: raise RuntimeError(sol.message)
            u,w=sol.y
            rho=energy(SAMPLE,u,w,k,m)
            assert np.all(np.isfinite(rho)) and np.all(rho>0)
            analytic,coef=exact_massive(a,k,m,SAMPLE[[0,100,300,600]])
            numerical=sol.y[:,[0,100,300,600]]
            mismatch=float(np.max(abs(numerical-analytic)/(1+abs(analytic))))
            assert mismatch<1e-6, (a,k,m,mismatch)
            slope=float(np.polyfit(SAMPLE[-121:],np.log(rho[-121:]),1)[0])
            expected=2-2*exponent(a)
            assert abs(slope-expected)<1e-3, (a,k,m,slope,expected)
            item=dict(a=a,k=k,mass=m,exponent=exponent(a),success=bool(sol.success),
                      nfev=sol.nfev,energy_amplification=float(rho[-1]/rho[0]),
                      final_log_energy_slope=slope,expected_asymptotic_slope=expected,
                      exact_relative_error=mismatch,singular_branch_coefficient=coef)
            if m==0:
                uexact=np.exp(-p*SAMPLE)
                item['massless_analytic_error']=float(np.max(abs(u-uexact)/(1+abs(uexact))))
            cases.append(item)
            if (a,k,m) in ((.02,1,0.),(.02,1,1.),(.02,1,10.),(.12,8,10.)):
                key=f'a{a}_k{k}_m{m}'
                trajectories[key]=dict(s=SAMPLE.tolist(),u_real=u.real.tolist(),
                    u_imag=u.imag.tolist(),energy=rho.tolist())

# Independent implicit solver on a selected massive curved-background run.
a,k,m=.02,2,1.
p=exponent(a)+2j*k
s_ind=np.linspace(0,16,321)
def rhs_real(s,y):
    val=rhs(s,np.array([y[0]+1j*y[1],y[2]+1j*y[3]]),p,m)
    return [val[0].real,val[0].imag,val[1].real,val[1].imag]
init=[1.,0.,(-p).real,(-p).imag]
rad=solve_ivp(rhs_real,(0,16),init,method='Radau',rtol=1e-10,atol=1e-12,t_eval=s_ind)
assert rad.success
explicit=solve_ivp(rhs,(0,16),np.array([1,-p],complex),args=(p,m),method='DOP853',
                   rtol=1e-11,atol=1e-13,t_eval=s_ind)
rad_complex=np.vstack([rad.y[0]+1j*rad.y[1],rad.y[2]+1j*rad.y[3]])
solver_error=float(np.max(abs(rad_complex-explicit.y)/(1+abs(explicit.y))))
assert solver_error<1e-7

# A fine-tuned regular massless branch is a necessary negative/control test.
control=solve_ivp(rhs,(0,S_END),np.array([1,0],complex),args=(exponent(.02)+2j,0.),
                   method='DOP853',rtol=1e-11,atol=1e-13,t_eval=SAMPLE)
assert control.success
assert np.max(abs(control.y[0]-1))<1e-12 and np.max(abs(control.y[1]))<1e-12
rho_control=energy(SAMPLE,control.y[0],control.y[1],1,0.)
assert np.max(rho_control)<=1+1e-12

# Reducing initial amplitude cannot change the linear amplification factor.
selected=next(c for c in cases if (c['a'],c['k'],c['mass'])==(.02,1,0.))
amplitudes=[1.,1e-4,1e-8,1e-12]
a_rows=[dict(initial_amplitude=eps,final_energy_relative_to_unit_initial_energy=
             eps*eps*selected['energy_amplification']) for eps in amplitudes]
summary=dict(
  scope='Exact-mode test-field evolution on a fixed curved vacuum core; not full numerical relativity.',
  count=len(cases),nonzero_curvature_cases=sum(c['a']>0 for c in cases),
  massless_cases=sum(c['mass']==0 for c in cases),massive_cases=sum(c['mass']>0 for c in cases),
  interval=dict(initial_S=-1.,final_S=-1e-12,log_time_end=S_END),
  tolerances=dict(rtol=RTOL,atol=ATOL),
  min_amplification=min(c['energy_amplification'] for c in cases),
  max_amplification=max(c['energy_amplification'] for c in cases),
  max_independent_exact_error=max(c['exact_relative_error'] for c in cases),
  max_asymptotic_slope_error=max(abs(c['final_log_energy_slope']-c['expected_asymptotic_slope']) for c in cases),
  min_singular_branch_coefficient=min(c['singular_branch_coefficient'] for c in cases),
  independent_solver_relative_error=solver_error,
  regular_branch_control=dict(bounded=True,initial_energy=float(rho_control[0]),final_energy=float(rho_control[-1])),
  selected_case=selected,
  amplitude_scaling=a_rows,
  elapsed_seconds=time.perf_counter()-start,
  versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
  interpretation='Instability witness for this background/mode family. Final coupled geometry and all alternative media remain undetermined.'
)
(OUT/'mode_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
(OUT/'sample_trajectories.json').write_text(json.dumps(trajectories,indent=2)+'\n')
(OUT/'numerical_results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
