"""Try a positive power-law kinetic coupling to suppress the scalar instability.

Z(S)=(-S)^(-q) is treated as a prescribed positive scalar background. We do NOT
claim to have supplied its dynamics or gravitational source. This generous
fixed-background repair already fails to keep both mode branches bounded.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp

OUT=Path(__file__).resolve().parent
s,k,q,lam=sp.symbols('s k q lam',real=True)
p=s+2*sp.I*k
poly=lam**2-(p+q)*lam+sp.I*q*k
# If roots are lambda_1,lambda_2, their real parts sum to s+q.
energy_exponent_sum=sp.expand(2*(s+q)-2*q-4)
assert energy_exponent_sum==2*s-4
rows=[]
for a in [0.,.001,.01,.02,.04,.06,.08,.10,.12]:
    ss=(np.sqrt(1-8*a)+np.sqrt(1+8*a)-2)/4
    for kk in [1,2,4,8]:
        for qq in np.linspace(-20,20,81):
            pp=ss+2j*kk
            roots=np.roots([1.,-(pp+qq),1j*qq*kk])
            exps=2*roots.real-qq-2
            nonconstant=np.abs(roots)>1e-10
            failed=bool(np.any(nonconstant & (exps<0)))
            assert failed
            assert abs(sum(exps)-(2*ss-4))<1e-10
            rows.append(dict(a=a,k=kk,q=float(qq),
                 roots=[[float(z.real),float(z.imag)] for z in roots],
                 leading_energy_exponents=exps.tolist(),
                 has_nonconstant_divergent_branch=failed))
summary={
 'candidate':'Prescribed positive scalar kinetic coupling Z(S)=(-S)^(-q)',
 'field_equation':'S*u_pp + (kappa-q-2*i*k)*u_p + i*q*k*u/S = 0',
 'indicial_polynomial':str(poly),
 'root_sum':'lambda_1+lambda_2 = s+q+2*i*k',
 'sum_of_leading_energy_exponents':str(energy_exponent_sum),
 'analytic_result':'For s<=0 and nonzero k, at least one nonconstant mode has divergent measured energy for every real q. q=0 has a regular constant branch but also the divergent original branch.',
 'algebraic_grid_cases':len(rows),
 'grid_cases_with_divergent_branch':sum(r['has_nonconstant_divergent_branch'] for r in rows),
 'time_domain_simulations':0,
 'scope':'Fixed-background diagnostic; Z dynamics and its own gravitational contribution have not been solved.'
}
(OUT/'kinetic_repair_cases.json').write_text(json.dumps(rows,indent=2)+'\n')
(OUT/'kinetic_repair_results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
