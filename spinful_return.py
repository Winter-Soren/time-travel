"""Spinful exterior operator, exact transverse benchmark, and local 2D evolution.

The periodic propagation test is not a mouth-to-mouth scattering matrix.
No outgoing-boundary return probability or Casimir update is inferred here.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import betaincc, gammaincc


def mode_log_norm(q, k):
    """Norm of the raw density in r=rho/(q*d); k=0 has no finite norm."""
    if not (isinstance(q,int) and q>0 and isinstance(k,int) and 1<=k<q):
        raise ValueError('Normalizable physical slice modes require integer 1<=k<q')
    with mp.workdps(40):
        a = mp.mpf(q-k); b = mp.mpf(2*k)
        return float(b*mp.log(q)-mp.log(2)+mp.log(mp.beta(a,b))+mp.log1p(a/(a+b)))


def mode_density(q, k, r):
    r = np.asarray(r,dtype=float)
    if np.any(r<=0) or np.any(~np.isfinite(r)): raise ValueError('Require finite r>0')
    return np.exp(-(2*k+1)*np.log(r)-2*q*np.arcsinh(1/(2*q*r))-mode_log_norm(q,k))


def mode_cdf(q, k, r):
    mode_log_norm(q,k)  # Shared domain validation; excludes the marginal channel.
    r = np.asarray(r,dtype=float)
    if np.any(r<=0) or np.any(~np.isfinite(r)): raise ValueError('Require finite r>0')
    y_complement = -np.expm1(-2*np.arcsinh(1/(2*q*r)))
    a,b = q-k,2*k
    return ((a+b)*betaincc(b,a,y_complement)+a*betaincc(b,a+1,y_complement))/(2*a+b)


def operator_identity():
    import sympy as s
    r,z = s.symbols('r z',positive=True); W = s.Function('W')(r,z)
    f = s.Matrix([s.Function('f')(r,z),s.Function('g')(r,z)])
    sx=s.Matrix([[0,1],[1,0]]); sy=s.Matrix([[0,-s.I],[s.I,0]]); sz=s.diag(1,-1)
    H=lambda v: -s.I*sx*v.diff(r)-s.I*sz*v.diff(z)+sy*W*v
    expected=-f.diff(r,2)-f.diff(z,2)+W*W*f+sz*W.diff(r)*f-sx*W.diff(z)*f
    residual=[s.simplify(v) for v in H(H(f))-expected]
    assert residual==[0,0]
    return [str(v) for v in residual]


def local_evolution(n=128,dt=.002,mu=16.5,inner=.03,outer=4.,zouter=1.5,finish=.4):
    """Unitary Fourier-Strang bulk check in the scaled dipole tail, on a torus."""
    if not (n>=32 and dt>0 and mu>.5 and 0<inner<.5 and outer>2 and zouter>1 and finish>0):
        raise ValueError('Invalid local evolution resolution/domain')
    rho=np.linspace(inner,outer,n,endpoint=False); z=np.linspace(-zouter,zouter,n,endpoint=False)
    dr=rho[1]-rho[0]; dz=z[1]-z[0]; R,Z=np.meshgrid(rho,z,indexing='ij')
    W=mu*(1-R*R/(R*R+Z*Z)**1.5)/R
    initial=np.exp(-mu*np.log(R)-mu/R+mu-Z*Z/(4*.12**2))
    f=np.stack((np.zeros_like(initial),initial)).astype(complex)
    normalization=math.sqrt(float(np.sum(abs(f)**2)*dr*dz)); f/=normalization
    count=math.ceil(finish/dt); step=finish/count
    cv,sv=np.cos(step*W/2),np.sin(step*W/2)
    kr=2*math.pi*np.fft.fftfreq(n,dr)[:,None]; kz=2*math.pi*np.fft.fftfreq(n,dz)[None,:]
    momentum=np.hypot(kr,kz); ck=np.cos(step*momentum); sk=step*np.sinc(step*momentum/math.pi)
    edge=(R<inner+.08)|(R>outer-.08)|(abs(Z)>zouter-.08)
    def reference(time, residual=False):
        g=np.stack((np.zeros_like(R),np.exp(-mu*np.log(R)-mu/R+mu-(Z+time)**2/(4*.12**2))/normalization)).astype(complex)
        if not residual: return g
        a,b=np.fft.fft2(g,axes=(-2,-1))
        Hg=np.fft.ifft2(np.stack((kz*a+kr*b,kr*a-kz*b)),axes=(-2,-1))
        Hg+=np.stack((-1j*W*g[1],1j*W*g[0]))
        derivative=-(Z+time)/(2*.12**2)*g
        return math.sqrt(float(np.sum(abs(1j*derivative-Hg)**2)*dr*dz))
    def observables():
        density=np.sum(abs(f)**2,axis=0); norm=float(density.sum()*dr*dz)
        return {'norm':norm,'mean_rho':float(np.sum(R*density)*dr*dz),
            'mean_z':float(np.sum(Z*density)*dr*dz),
            'spin_z':float(np.sum(abs(f[0])**2-abs(f[1])**2)*dr*dz),
            'edge_probability':float(density[edge].sum()*dr*dz)}
    history=[dict(t=0.,**observables())]
    for j in range(count):
        a,b=f; f=np.stack((cv*a-sv*b,sv*a+cv*b))
        a,b=np.fft.fft2(f,axes=(-2,-1))
        f=np.fft.ifft2(np.stack((ck*a-1j*sk*(kz*a+kr*b),ck*b-1j*sk*(kr*a-kz*b))),axes=(-2,-1))
        a,b=f; f=np.stack((cv*a-sv*b,sv*a+cv*b))
        if (j+1)%max(1,count//20)==0 or j==count-1:
            history.append(dict(t=(j+1)*step,**observables()))
    residual_integrals=[]
    for order in (24,48):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        residual_integrals.append(float(finish/2*np.dot(weights,[reference((v+1)*finish/2,True) for v in nodes])))
    reference_error=math.sqrt(float(np.sum(abs(f-reference(finish))**2)*dr*dz))
    return {'grid':n,'step':step,'mu':mu,'inner':inner,'outer':outer,'zouter':zouter,'finish':finish,
        'maximum_norm_error':max(abs(row['norm']-1) for row in history),
        'maximum_edge_probability':max(row['edge_probability'] for row in history),
        'translated_reference_error':reference_error,'duhamel_residual_integrals':residual_integrals,
        'duhamel_scope':'The exact residual integral bounds exact semidiscrete Hamiltonian evolution. The saved integral is numerical quadrature; the Strang trajectory is only observed below it. Splitting and quadrature errors are not rigorously enclosed.',
        'history':history,'rho':rho.tolist(),'z':z.tolist(),
        'final_density':np.sum(abs(f)**2,axis=0).tolist(),
        'scope':'Specified packet evolved briefly on a periodic computational box. Length unit q*d/(2mu), time unit length/c; neither asymptotic scattering nor a throat injection.'}


def spinful_checks():
    q=100000; r=np.geomspace(.003,20.,500); profiles=[]; references=[]
    for k in (1,2,8,16):
        density=mode_density(q,k,r)
        quantiles=[brentq(lambda x: float(mode_cdf(q,k,x))-p,1e-4,100.,xtol=2e-13) for p in (.05,.5,.95)]
        profiles.append({'k':k,'mu':k+.5,'r_over_qd':r.tolist(),'density':density.tolist(),
            'radial_quantiles_05_50_95':quantiles,'mean_flux_sphere':2*(k+1)/(q+1),
            'mean_flux_slice':2*k/q,'theta_density_mode_flux':(2*k+1)/q,
            'flux_density_mode_sphere':2*k/(q-1),'second_radial_moment_finite':k>1})
        # Compactified quadrature is independent of the incomplete-Beta CDF.
        integrand=lambda t: float(mode_density(q,k,t/(1-t)))/(1-t)**2
        norm=quad(integrand,0.,1.,epsabs=2e-10,epsrel=2e-10)[0]
        cdf=quad(integrand,0.,.5,epsabs=2e-10,epsrel=2e-10)[0]
        meanA=quad(lambda t: integrand(t)/(2*math.hypot(t/(1-t),1/(2*q))),
                   0.,1.,epsabs=2e-9,epsrel=2e-10)[0]
        references.append({'k':k,'norm_error':abs(norm-1),
            'cdf_error':abs(cdf-float(mode_cdf(q,k,1.))), 'mean_A_error':abs(meanA-k)})
    shell=[]
    for L in (10.,100.,1000.):
        # Marginal k0 raw norm gains ln2 in each doubled far-radius shell.
        value=quad(lambda x: math.exp(-2*q*math.asinh(1/(2*q*x)))/x,L,2*L,
                   epsabs=1e-12,epsrel=1e-12)[0]
        shell.append({'L_over_qd':L,'raw_norm_in_doubled_shell':value})
    asymptotic=[]
    for charge in (100,1000,100000):
        sample=np.geomspace(.1,10.,31)
        errors=abs(mode_cdf(charge,1,sample)-gammaincc(2,1/sample))
        asymptotic.append({'q':charge,'maximum_cdf_difference':float(errors.max())})
    with mp.workdps(50):
        residuals=[]
        for k in (0,1,16):
            mu=mp.mpf(k)+mp.mpf('.5')
            logu=lambda x:-mu*mp.log(x)-q*mp.asinh(1/(2*x))
            for point in (mp.mpf('.2'),mp.mpf(q)/(2*mu),mp.mpf(q)*10):
                W=(mu-mp.mpf(q)/(2*mp.sqrt(point*point+mp.mpf('.25'))))/point
                residuals.append(float(abs(mp.diff(logu,point)+W)/(1+abs(W))))
    operator=operator_identity()
    configurations=[dict(n=128,dt=.004),dict(n=128,dt=.002),dict(n=256,dt=.002),
                    dict(n=256,dt=.001),dict(n=320,dt=.001,inner=.015,outer=4.8,zouter=1.8)]
    evolution=[local_evolution(**cfg) for cfg in configurations]
    def difference(a,b):
        return max(abs(a['history'][-1][key]-b['history'][-1][key]) for key in ('mean_rho','mean_z','spin_z'))
    differences={'time_coarse':difference(evolution[0],evolution[1]),
        'spatial':difference(evolution[1],evolution[2]),'time_fine':difference(evolution[2],evolution[3]),
        'outer_domain':difference(evolution[3],evolution[4])}
    tests={
        'exact_spinful_operator_square':operator==['0','0'],
        'exact_transverse_zero_mode':max(residuals)<1e-35,
        'normalization_cdf_and_flux_mean':max(max(row['norm_error'],row['cdf_error'],row['mean_A_error']) for row in references)<3e-8,
        'marginal_channel_not_normalized':abs(shell[-1]['raw_norm_in_doubled_shell']-math.log(2))<.001 and shell[-1]['raw_norm_in_doubled_shell']>shell[0]['raw_norm_in_doubled_shell'],
        'far_dipole_limit':asymptotic[-1]['maximum_cdf_difference']<1e-8 and asymptotic[1]['maximum_cdf_difference']<asymptotic[0]['maximum_cdf_difference']/90,
        'unitary_bulk_evolution':max(e['maximum_norm_error'] for e in evolution)<2e-11,
        'bulk_evolution_refinement':max(differences.values())<2e-4,
        'finite_box_packet_control':max(e['maximum_edge_probability'] for e in evolution)<2e-6,
        'semidiscrete_reference_diagnostic':all(e['translated_reference_error']<e['duhamel_residual_integrals'][-1] and abs(e['duhamel_residual_integrals'][0]-e['duhamel_residual_integrals'][1])<2e-7 for e in evolution),
    }
    return {'passed':all(tests.values()),'tests':tests,'q':q,'profiles':profiles,'references':references,
        'planar_normalizable_zero_mode_count':q-1,'marginal_k':0,
        'operator_square_residuals':operator,'max_log_zero_mode_residual':max(residuals),
        'marginal_shells':shell,'far_limit':asymptotic,'evolution':evolution,'evolution_differences':differences,
        'equations':{'operator':'Hmu=-i sigma_x d_rho-i sigma_z d_z+sigma_y (mu-Aphi)/rho; measure d_rho d_z; mu=k+1/2.',
            'gauge':'Aphi->Aphi+n, mu->mu+n; regular midplane patch n=-q. Standard Cartesian curl, opposite the printed MMP Eq3.14 convention.',
            'current':'d_t(f†f)+d_rho(f†sigma_x f)+d_z(f†sigma_z f)=0.',
            'zero_mode':'u=rho^(-mu)*exp[-q asinh(d/(2rho))]; integral|u|²d_rho finite iff1/2<mu<q+1/2.',
            'moments':'Moment p finite iff -2(q-mu+1/2)<p<2mu-1.',
            'far_density':'rbar=rho/(q*d); p(rbar)=rbar^(-2mu)*exp(-1/rbar)/Gamma(2mu-1).'},
        'limitations':['The planar zero modes reproduce the known Aharonov-Casher integer-flux result, not a new law.',
            'The planar count is not the number of 3D returning channels; marginality is not an escape probability.',
            'No global transverse gap separates these zero modes from the continuum.',
            'Bulk unitary evolution uses a short specified packet in a periodic box, not incoming throat modes or outgoing scattering conditions.',
            'Full return amplitude/phase, continuum vacuum stress, physical finite-mouth matching and quantum current remain unresolved.',
            'No stable medium or 13 Sep2027 to12 Sep2027 passage is established.'],
        'sources':['https://arxiv.org/html/1807.04726v3','https://arxiv.org/pdf/math-ph/0109015','https://arxiv.org/pdf/1801.01346'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result=spinful_checks()
    Path(__file__).with_name('spinful_return_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','tests','references','marginal_shells','far_limit','evolution_differences')},indent=2))
    print(json.dumps([{'grid':e['grid'],'dt':e['step'],'norm_error':e['maximum_norm_error'],
        'edge_probability':e['maximum_edge_probability'],'last':e['history'][-1]} for e in result['evolution']],indent=2))
    if not result['passed']: raise SystemExit('A spinful benchmark control failed; diagnostics retained.')
