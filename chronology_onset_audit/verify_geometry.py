"""Exact symbolic checks for an Ori vacuum core and a scalar instability witness.

Coordinates and all parameters are dimensionless; an overall length scale may be
restored. This does not construct a laboratory exterior or solve backreaction.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parent
T, x, y, z = sp.symbols('T x y z', real=True)
a, e, b, d, k = sp.symbols('a e b d k', real=True)
coords = (T, x, y, z)
f = a * (x*x-y*y) / 2
g = sp.Matrix([[0,0,0,-1],[0,1,0,0],[0,0,1,0],[-1,0,0,f-T]])
gi = g.inv()
Gamma = [[[sp.simplify(sum(gi[i,l]*(sp.diff(g[l,j],coords[h])+
             sp.diff(g[l,h],coords[j])-sp.diff(g[j,h],coords[l]))/2
             for l in range(4))) for h in range(4)] for j in range(4)] for i in range(4)]
Ricci = sp.Matrix(4,4, lambda i,j: sp.simplify(sum(
    sp.diff(Gamma[h][i][j],coords[h])-sp.diff(Gamma[h][i][h],coords[j])+
    sum(Gamma[h][h][l]*Gamma[l][i][j]-Gamma[h][j][l]*Gamma[l][i][h]
        for l in range(4)) for h in range(4))))
assert Ricci == sp.zeros(4)
assert g.det() == -1

def gradnorm(F):
    return sp.expand(sum(gi[i,j]*sp.diff(F,coords[i])*sp.diff(F,coords[j])
                         for i in range(4) for j in range(4)))
def box(F):
    return sp.simplify(sum(sp.diff(gi[i,j]*sp.diff(F,coords[j]),coords[i])
                          for i in range(4) for j in range(4)))

# Ori's chronology-adapted coordinate t, not the mode coordinate S below.
t = sp.symbols('t',real=True)
clock = T-f+e*(x*x+y*y)
clock_norm = sp.expand(gradnorm(clock).subs(T,t+f-e*(x*x+y*y)))
expected = t+((2*e-a)**2-e)*x*x+((2*e+a)**2-e)*y*y
assert sp.simplify(clock_norm-expected)==0

# A quadratic scalar coordinate closes the exact wave operator on functions of S.
S = T+b*x*x+d*y*y
conditions = (4*b*b-b-a/2, 4*d*d-d+a/2)
Grobner = sp.groebner(conditions,b,d,order='lex',domain=sp.QQ_I.frac_field(a))
norm_residual = sp.Poly(gradnorm(S)-S,x,y,T)
assert all(sp.simplify(Grobner.reduce(c)[1])==0 for c in norm_residual.coeffs())
kappa = 1+2*b+2*d
assert sp.simplify(box(S)-kappa)==0
assert sp.simplify(sum(gi[i,3]*sp.diff(S,coords[i]) for i in range(4))+1)==0
assert box(z)==0 and gradnorm(z)==0
p = -2*(b+d)+2*sp.I*k
assert sp.simplify(p*(p+kappa-1-2*sp.I*k))==0

# Check the full four-dimensional wave equation, not only an asserted ODE.
loggrad = [p*sp.diff(S,c)/S+(sp.I*k if i==3 else 0) for i,c in enumerate(coords)]
kg_over_phi = sp.factor(sum(
    sp.diff(gi[i,j],coords[i])*loggrad[j] + gi[i,j]*(
    sp.diff(loggrad[j],coords[i])+loggrad[i]*loggrad[j])
    for i in range(4) for j in range(4)))
kg_num = sp.together(kg_over_phi).as_numer_denom()[0]
coefficients = sp.Poly(kg_num,x,y,T,k).coeffs()
kg_checks = [sp.simplify(Grobner.reduce(c)[1]) for c in coefficients]
assert all(c==0 for c in kg_checks)

broot=(1-sp.sqrt(1+8*a))/8
droot=(1-sp.sqrt(1-8*a))/8
exponent=sp.simplify(-2*(broot+droot))
assert sp.simplify(conditions[0].subs(b,broot))==0
assert sp.simplify(conditions[1].subs(d,droot))==0

# A smooth, finite-acceleration observer reaching the central null circle.
ga = g.subs({x:0,y:0})
u=sp.Matrix([1,0,0,1])/sp.sqrt(2+T)
n=sp.Matrix([T+1,0,0,-1])/sp.sqrt(2+T)
assert sp.simplify((u.T*ga*u)[0])==-1
assert sp.simplify((n.T*ga*n)[0])==1
assert sp.simplify((u.T*ga*n)[0])==0
acc=sp.Matrix([sp.simplify(u[0]*sp.diff(u[i],T)+sum(
    Gamma[i][j][h].subs({x:0,y:0})*u[j]*u[h]
    for j in range(4) for h in range(4))) for i in range(4)])
acc2=sp.simplify((acc.T*ga*acc)[0])
acc2_limit=sp.limit(acc2,T,0,dir='-')
assert acc2_limit.is_finite

# Positive measured energy of two real canonical scalars, phi1+i phi2.
s=sp.symbols('s',real=True)
derivative=sp.Matrix([(s+2*sp.I*k)/T,0,0,sp.I*k])
P=(s+2*sp.I*k)/T
Z=sp.I*k
rho_factor=sp.simplify(((P+Z)*sp.conjugate(P+Z)+
                       ((T+1)*P-Z)*sp.conjugate((T+1)*P-Z))/(2*(2+T)))
expected_rho=k*k*(T+2)/(T*T)+s*s*(T*T+2*T+2)/(2*T*T*(T+2))
assert sp.simplify(rho_factor-expected_rho)==0
limit_coeff=sp.limit(T*T*rho_factor,T,0,dir='-')
assert sp.simplify(limit_coeff-(2*k*k+s*s/2))==0

res={
  'metric_determinant':str(g.det()),
  'all_16_Ricci_components_zero':True,
  'chronology_time_gradient_norm':str(clock_norm),
  'massless_wave_equation_all_terms_zero':True,
  'quadratic_coordinate_conditions':[str(c) for c in conditions],
  'b_small_root':str(broot),'d_small_root':str(droot),
  'amplitude_exponent_s':str(exponent),
  'mode':'phi = epsilon * exp(i k z) * (-S)^(s+2 i k); S=T+b*x^2+d*y^2 < 0',
  'kappa':str(kappa),
  'observer_squared_norm':-1,
  'observer_acceleration_squared':str(acc2),
  'observer_acceleration_squared_limit_at_T0':str(acc2_limit),
  'energy_over_epsilon_squared_minusT_to_2s':str(expected_rho),
  'energy_asymptotic_coefficient':str(limit_coeff),
  'energy_power':'rho ~ epsilon^2 * (2*k^2+s^2/2) * (-T)^(2*s-2)',
  'scope':'Exact background and exact test-field checks, not nonlinear Einstein-scalar evolution.'
}
(OUT/'symbolic_results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
