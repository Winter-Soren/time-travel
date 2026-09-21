"""Uniform exterior coercivity for the actual sewn geometry's zero states.

The all-position/all-sector claim follows from the analytic inequalities in
PHYSICAL_EUCLIDEAN_OPERATOR.md and the interval constants below. Sampled
curves only check the formulas; no zero eigenstate or supporting stress is fitted.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from balanced_dihole import parameters, symbolic_fields
from project_paths import project_path


def certificate():
    """Outward interval arithmetic for the linked exact Q=1/5, D=1, q=1."""
    iv, old = mp.iv, mp.iv.dps
    iv.dps = 70
    try:
        k, Q, R0 = iv.mpf(1)/2, iv.mpf(1)/5, iv.mpf(100)
        M = (iv.sqrt(k*k+4*Q*k)-k)/2
        a = iv.sqrt(k*k-M*M)
        B = M/((k+a)*(M+k))
        b, F, e = B*B/4, 1/(Q*B), 1-k*k
        assert F.a > iv.mpf('20.5') and F.b < iv.mpf(21)
        g = F-iv.mpf('20.5')
        assert (B*a).a > 0 and (B*a).b < 1
        U0 = 2*M*(1+M)/e
        U1 = 2*M*(1+2*M+k*k)/e**2
        S1 = 1/iv.sqrt(e)
        T1 = a*a*(1/(e*iv.sqrt(e))+2/e**2)
        A = 1/R0**2+b*(1+U0/R0)
        C = (2*b*(1+U0/R0)+b*U1/R0+b*U0*T1/R0**3
             +(U1+B*a*U0*S1)/R0**3+U0*T1/R0**5)
        K1 = U1+B*a*U0*S1/2+U0*T1/R0**2
        A1, A0 = 2*A*C+A*A, F*(C+A+A*K1/R0)/R0**2
        amplitude_margin = b*R0**2-2*F/g
        derivative_margin = b**4*R0**4/8-A1/g-A0/g**2
        c0 = g*g*b**4/8
        assert amplitude_margin.a > 0 and derivative_margin.a > 0 and c0.a > 0
        values = dict(B=B, beta=b, flux=F, half_integer_gap=g, U0=U0, U1=U1,
                      S1=S1, T1=T1, N_upper=A, N_gradient=C, K_gradient=K1,
                      angular_gradient_linear=A1, angular_gradient_constant=A0,
                      amplitude_margin=amplitude_margin,
                      derivative_margin=derivative_margin, coercivity=c0)
        return dict(passed=True, digits=70, radius=100,
                    calibration='Q=1/5, D=1, q=1; all lengths in D units',
                    intervals={key:str(value) for key,value in values.items()},
                    display_values={key:float(value.mid) for key,value in values.items()},
                    assertion='For every rho>=100, z real and mu in Z+1/2: W_mu^2-v*|grad W_mu| >= c0*rho^6.')
    finally:
        iv.dps = old


def normalized_fields():
    """Stable exact formulas: no subtraction of nearly equal focal distances."""
    R, s, M, a, B, rho = sp.symbols('R s M a B rho', real=True)
    d = R*R-M*M-a*a
    u, t = 2*M*(R+M)/d, a*a*s/d
    den = 1+u+t
    G = ((1+u)**2+t)/den
    L = (1+t+B*a*u*s)/den
    K = (1+t+B*a*u*s/2)/den
    N = L+B*B*rho*rho*G/4
    # Compare against the original published-coordinate equations symbolically.
    r, sigma = R+M, d+2*M*(R+M)+a*a*s
    w = (r*r-a*a)**2+d*a*a*s
    original_N = (d+a*a*s+2*B*M*r*a*s+B*B*w*s/4)/sigma
    gauge_ratio = (B*M*r*a+B*B*w/4)*s/(original_N*sigma)
    identities = dict(
        exact_lapse=sp.cancel((N-original_N).subs(rho*rho,d*s)) == 0,
        exact_gauge=sp.cancel((1-K/N-gauge_ratio).subs(rho*rho,d*s)) == 0)
    values = [N,K]+[sp.diff(f,x) for f in (N,K) for x in (rho,R,s)]
    fn = sp.lambdify((rho,R,s,M,a,B),values,'numpy',cse=True)
    return fn, identities


def coefficients(fn,p,rho,z,mu):
    """Canonical flat-measure angular potential and its exact Weyl gradient."""
    rho,z = np.broadcast_arrays(np.asarray(rho,float),np.asarray(z,float))
    if np.any(rho<=0) or not np.all(np.isfinite(rho+z)):
        raise ValueError('Finite positions with rho>0 required.')
    M,a,k,B = [p[key] for key in ('M','a','k','B')]
    R = (np.hypot(rho,z-k)+np.hypot(rho,z+k))/2
    d = (R-k)*(R+k)
    s, c = rho*rho/d, z/R
    h, sigma = d+k*k*s, d+2*M*(R+M)+a*a*s
    gradR = np.array([rho*R/h,c*d/h])
    grads = np.array([2*rho*c*c/h,-2*R*c*s/h])
    N,K,Nrho,NR,Ns,Krho,KR,Ks = [np.broadcast_to(v,rho.shape) for v in fn(rho,R,s,M,a,B)]
    gradN = NR*gradR+Ns*grads+np.array([Nrho,np.zeros_like(Nrho)])
    gradK = KR*gradR+Ks*grads+np.array([Krho,np.zeros_like(Krho)])
    F, delta = p['q']/(p['Q']*B), mu-p['q']/(p['Q']*B)
    W = (delta*N*N+F*N*K)/rho
    gradW = (2*delta*N*gradN+F*(K*gradN+N*gradK))/rho
    gradW[0] -= W/rho
    v = (h/sigma)**2
    return dict(N=N,K=K,v=v,W=W,gradN=gradN,gradK=gradK,gradW=gradW,
                coercivity=W*W-v*np.linalg.norm(gradW,axis=0),R=R,c=c)


def run():
    cert = certificate()
    fn, identities = normalized_fields()
    p, cut = parameters(), symbolic_fields()[0]
    val, rows = cert['display_values'], []
    errors = []
    # Independent original-coordinate field and derivative check, before tail sampling.
    for rho,z,mu in ((1.,0.,.5),(3.,-8.,21.5),(100.,0.,20.5),(200.,1600.,20.5)):
        f = coefficients(fn,p,rho,z,mu)
        r,c = f['R']+p['M'],f['c']
        n,P,ahat,nr,_,ar,nc,_,ac = cut(r,c,p['M'],p['a'],p['B'])
        h = f['R']**2-p['k']**2*c*c
        gr = np.array([rho*f['R']/h,c*(f['R']**2-p['k']**2)/h])
        gc = np.array([-rho*c/h,f['R']*(1-c*c)/h])
        aa = ahat*p['physical_A_factor']
        gn,ga = nr*gr+nc*gc,(ar*gr+ac*gc)*p['physical_A_factor']
        w = n*n*(mu-p['q']*aa)/rho
        gw = (2*n*gn*(mu-p['q']*aa)-n*n*p['q']*ga)/rho
        gw[0] -= w/rho
        errors.append(max(abs(f['N']-n)/max(1,abs(n)),abs(f['v']-n/P),
                          abs(f['W']-w)/max(1,abs(w)),
                          np.linalg.norm(f['gradW']-gw)/max(1,np.linalg.norm(gw))))
    ratios = np.array([0.,.1,1.,10.,1e3,1e6,1e9])
    minima, derivative_ratios = [], []
    for rho in np.geomspace(100,1e5,31):
        sector_minima = []
        for mu in (-1000.5,-.5,.5,20.5,21.5,1000.5):
            f = coefficients(fn,p,rho,rho*ratios,mu)
            delta = abs(mu-val['flux'])
            bound = (delta*val['angular_gradient_linear']+val['angular_gradient_constant'])*rho*rho
            derivative_ratios.append(float(np.max(np.linalg.norm(f['gradW'],axis=0)/bound)))
            sector_minima.append(float(np.min(f['coercivity']/(val['coercivity']*rho**6))))
        minima.append(min(sector_minima))
        rows.append(dict(rho=float(rho),minimum_sampled_coercivity_ratio=minima[-1],
                         sampled_lapse_over_rho_squared=(f['N']/rho**2).tolist()))
    sx,sy,sz = [sp.Matrix(x) for x in ([[0,1],[1,0]],[[0,-sp.I],[sp.I,0]],[[1,0],[0,-1]])]
    checks = dict(**identities,
        interval_all_positions_all_sectors=cert['passed'],
        pauli_cross_term=(-sp.I*sx*sy==sz and -sp.I*sz*sy==-sx),
        original_field_and_gradient_agreement=max(errors)<2e-11,
        sampled_coercivity_obeys_analytic_bound=min(minima)>=1,
        sampled_gradient_obeys_analytic_bound=max(derivative_ratios)<=1,
        actual_parameter_matches_exact_interval=abs(p['B']-val['B'])<1e-14)
    checks = {key:bool(value) for key,value in checks.items()}
    result = dict(passed=all(checks.values()),checks=checks,certificate=cert,
                  original_formula_max_relative_error=float(max(errors)),
                  maximum_sampled_gradient_bound_ratio=max(derivative_ratios),
                  minimum_sampled_coercivity_ratio=min(minima),rows=rows,
                  sampled_z_over_rho=ratios.tolist(),
                  scope='Analytic exterior bound plus interval-certified constants; samples are controls. Applies to zero solutions of the existing compact sewing; does not assert existence or absence of zero states, calculate physical stress, or establish a supported bridge.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  inputs_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                                 for name in ('balanced_dihole.py','project_paths.py')},
                  sources=['https://arxiv.org/pdf/hep-th/9906160'])
    project_path('dihole_zero_tail_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']:
        raise ArithmeticError('Zero-tail controls failed; results retained.')
    plot(result)
    print(json.dumps({key:result[key] for key in ('passed','checks','original_formula_max_relative_error',
          'maximum_sampled_gradient_bound_ratio','minimum_sampled_coercivity_ratio')},indent=2))
    print(json.dumps(cert['display_values'],indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows, val = result['rows'], result['certificate']['display_values']
    rho = [r['rho'] for r in rows]
    fig,axes = plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for j,ratio in enumerate(result['sampled_z_over_rho']):
        if j in (0,2,4,6):
            axes[0].semilogx(rho,[r['sampled_lapse_over_rho_squared'][j] for r in rows],label=f'z/rho = {ratio:g}')
    axes[0].axhline(val['beta'],color='k',ls='--',label='Analytic lower bound B^2/4')
    axes[0].set(title='Actual exterior lapse',ylabel='N / rho^2')
    axes[1].semilogx(rho,[r['minimum_sampled_coercivity_ratio'] for r in rows],label='Minimum of sampled positions and sectors')
    axes[1].axhline(1,color='k',ls='--',label='Certified bound: all z, all sectors')
    axes[1].set(title='Positive confinement beyond rho = 100 D',ylabel='(W^2 - v |grad W|) / (c0 rho^6)',ylim=(0,9))
    for ax in axes:
        ax.set_xlabel('Weyl radius rho / D');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Zero-state tail condition: proved for the specified geometry\nSupporting quantum stress and a stable passage remain unestablished',fontsize=12)
    fig.savefig(project_path('dihole_zero_tail_progress.png'),dpi=160)
    fig.savefig(project_path('dihole_zero_tail_progress.pdf'))
    plt.close(fig)


if __name__=='__main__':
    run()
