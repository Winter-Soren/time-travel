"""Finite initial energy on a compact part of Ori's spacelike initial slice."""
from pathlib import Path
import json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
OUT=Path(__file__).resolve().parent

def integral(n,a=.02,k=1):
    z,w=leggauss(n)
    r=(z+1)/2; wr=w/2
    theta=2*np.pi*np.arange(2*n)/(2*n)
    rr,th=np.meshgrid(r,theta,indexing='ij')
    x=rr*np.cos(th); y=rr*np.sin(th)
    e=1/16; P=2*e-a; Q=2*e+a
    b=(1-math.sqrt(1+8*a))/8
    d=(1-math.sqrt(1-8*a))/8
    s=-2*(b+d); p=s+2j*k
    f=.5*a*(x*x-y*y)
    T=-1+f-e*(x*x+y*y)
    S=T+b*x*x+d*y*y
    D=1+(e-P*P)*x*x+(e-Q*Q)*y*y
    assert np.all(S<0) and np.all(D>0)
    normal_dot=(-(T-f)-2*P*b*x*x-2*Q*d*y*y)*p/(S*np.sqrt(D))+1j*k/np.sqrt(D)
    rho=(-S)**(2*s)*(abs(normal_dot)**2+s*s/(2*S))
    assert np.all(rho>0)
    # z angle contributes 2 pi; sqrt(det h)=sqrt(D).
    val=2*np.pi*np.sum(wr[:,None]*rr*np.sqrt(D)*rho)*(2*np.pi/(2*n))
    return float(val)
coarse=integral(40); fine=integral(80)
assert abs(coarse-fine)<1e-9
out={
 'initial_slice':'t=-1, x^2+y^2<=1, 0<=z<2*pi',
 'a':.02,'k':1,'field_amplitude':1,
 'energy_in_dimensionless_model_units':fine,
 'quadrature_refinement_difference':abs(coarse-fine),
 'arbitrary_smallness':'Multiplying the two scalar components by epsilon multiplies this energy by epsilon^2.',
 'scope':'Compact core only; no matching to a complete external universe or gravitational constraint solve.'}
(OUT/'initial_energy_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
