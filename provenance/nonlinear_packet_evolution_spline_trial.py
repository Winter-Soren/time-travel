"""Nonlinear spherical canonical evolution of the action144 candidate.

Finite radial domain with even center and asymptotic r^-6 scalar tails.
This does not supply a material, an EFT validity bound or a time machine.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from nonlinear_packet_constraints import background, evaluate, solve
from project_paths import project_path


def derivatives(x, values):
    """Derivatives of the same even-center cubic interpolant used in constraints."""
    spline = CubicSpline(x,values,bc_type=((1,0.),(1,-6*values[-1]/x[-1])))
    return spline(x,1),spline(x,2)


def constraints(tau, zeta, mu, x, base):
    """m=(P+3c/N^2)trK=m0+P*H*mu; sigma=(beta_r-beta/r)/N=H*s."""
    P,c,X,H,R = [base[k] for k in ('P','c','X','H','radius')]
    zx,zxx = derivatives(x,zeta)
    mux,_ = derivatives(x,mu)
    radial = np.divide(zx,x,out=np.array(zxx),where=x!=0)
    R3 = np.exp(-2*(tau+zeta))/R**2*(-4*(zxx+2*radial)-2*zx*zx)
    # Integrate the regular momentum equation using Gaussian quadrature of
    # spline-interpolated fields; direct division of a Simpson integral by
    # x^3 would lose center accuracy in the first few cells.
    nodes,weights = np.polynomial.legendre.leggauss(4)
    step = x[1]-x[0]
    quadrature = x[:-1,None]+step*(nodes+1)/2
    zs = CubicSpline(x,zeta,bc_type=((1,0.),(1,-6*zeta[-1]/x[-1])))
    ms = CubicSpline(x,mu,bc_type=((1,0.),(1,-6*mu[-1]/x[-1])))
    integrand = quadrature**3*np.exp(3*zs(quadrature))*ms(quadrature,1)
    integral = np.r_[0.,np.cumsum(step/2*(integrand@weights))]
    shear = np.divide(-np.exp(-3*zeta)*integral,x**3,out=np.zeros_like(x),where=x!=0)
    dy = np.zeros_like(x)
    m0 = 3*H*(P+3*c)
    dm = P*H*mu
    m = m0+dm
    for _ in range(6):
        y = 1+dy
        dk = (dm-9*H*c*dy)/(P+3*c*y)
        gravity = (P+9*c)*(6*H*dk+dk*dk)/3+3*c*dy*(3*H+dk)**2
        residual = (.5*P*R3+gravity-P*H*H*shear*shear/3
                    +dy*(.5*X-1.5*X*X-.75*X*X*dy))
        slope = c*m*m*(P-9*c*y)/(P+3*c*y)**3-1.5*X*X*y+.5*X
        if np.any(slope >= -.5) or np.any(y <= 0):
            raise ArithmeticError(f'Canonical lapse branch failed at Ht={tau:.8g}: '
                                  f'min_y={np.min(y):.6g}, max_slope={np.max(slope):.6g}, '
                                  f'max_abs_zeta={np.max(abs(zeta)):.6g}, max_abs_mu={np.max(abs(mu)):.6g}.')
        dy -= residual/slope
    if np.max(abs(residual)) > 2e-12:
        raise ArithmeticError('Local nonlinear lapse constraint did not converge.')
    y = 1+dy
    sqrt_y = np.sqrt(y)
    n = -dy/(sqrt_y*(1+sqrt_y))
    lapse = 1+n
    dk = (dm-9*H*c*dy)/(P+3*c*y)
    integrand = np.divide(lapse*shear,x,out=np.zeros_like(x),where=x!=0)
    spline = CubicSpline(x,integrand,bc_type=((2,0.),(1,-4*integrand[-1]/x[-1]))).antiderivative()
    w = spline(x)-spline(x[-1])-shear[-1]/3
    return dict(zx=zx,zxx=zxx,mux=mux,curvature=R3,shear=shear,
                lapse=lapse,lapse_deviation=n,inverse_lapse_squared_deviation=dy,
                dk=dk,w=w,lapse_residual=residual)


def rhs(tau, values, x, base):
    zeta,mu = values.reshape(2,-1)
    f = constraints(tau,zeta,mu,x,base)
    P,c,X,H,R = [base[k] for k in ('P','c','X','H','radius')]
    N,n,dy,dk,s,w,zx = [f[k] for k in ('lapse','lapse_deviation',
        'inverse_lapse_squared_deviation','dk','shear','w','zx')]
    nx,nxx = derivatives(x,n)
    radial = np.divide(nx,x,out=np.array(nxx),where=x!=0)
    lap_n = np.exp(-2*(tau+zeta))/R**2*(nxx+2*radial+zx*nx)
    dm = P*H*mu
    m0 = 3*H*(P+3*c)
    trace_bracket = (3*H*dm+m0*dk+dm*dk+P*H*H*s*s
        +.75*(2*(X-1)*X*dy+X*X*dy*dy)+.5*P*f['curvature'])
    # Preserve the exact momentum cancellation before discretization. Directly
    # adding w+N*s/3 to dk/(3H) leaves a spurious grid-dependent kinetic term.
    correction = nx*(mu+s)-3*N*zx*s
    integral = CubicSpline(x,correction,bc_type=((2,0.),(1,-10*correction[-1]/x[-1]))).antiderivative()
    tail = correction[-1]*x[-1]/9
    I = integral(x[-1])-integral(x)+tail
    zeta_tau = n-N*c*((1+dy)*mu+3*dy)/(P+3*c*(1+dy))+x*w*zx-I/3
    mu_tau = x*w*f['mux']+lap_n/(H*H)-N*trace_bracket/(2*P*H*H)
    return np.r_[zeta_tau,mu_tau]


def initial(x, epsilon, base):
    solution,_ = solve(epsilon,xmax=max(12.,float(x[-1])+1.),tolerance=1e-11)
    initial_fields = evaluate(solution,base,epsilon,np.maximum(x,solution.x[0]))
    N,e = initial_fields['lapse'],initial_fields['e']
    P,c = base['P'],base['c']
    mu = ((P+3*c/N**2)*e/N-3*(P+3*c))/P
    zeta = -epsilon*np.exp(-x*x/2)
    return np.r_[zeta,mu],initial_fields


def evolve(epsilon=.001, *, step=.05, xmax=12., end=.5, rtol=2e-9, atol=2e-12, max_step=.005):
    if not (0 <= epsilon <= .001 and 0 < step <= .1 and xmax >= 10 and 0 < end <= 1):
        raise ValueError('Require tested amplitude/grid/domain/time ranges.')
    base = background()
    x = np.linspace(0.,xmax,round(xmax/step)+1)
    values,old = initial(x,epsilon,base)
    solution = solve_ivp(lambda t,y:rhs(t,y,x,base),(0.,end),values,method='DOP853',
        t_eval=np.linspace(0.,end,11),rtol=rtol,atol=atol,max_step=max_step)
    if not solution.success:
        raise ArithmeticError(solution.message)
    return dict(x=x,tau=solution.t,zeta=solution.y[:len(x)].T,
                mu=solution.y[len(x):].T,nfev=solution.nfev,base=base)


if __name__ == '__main__':
    for step in (.05,.025):
        state=evolve(step=step)
        print(json.dumps(dict(step=step,nfev=state['nfev'],
            center_zeta=state['zeta'][:,0].tolist(),center_mu=state['mu'][:,0].tolist())),flush=True)
