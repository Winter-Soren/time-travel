"""Exact classical balanced dihole exterior and finite throat-cut data.

Emparan, hep-th/9906160, equations8--16. This is an Einstein--Maxwell
black-hole exterior with Melvin asymptotics, not a traversable wormhole.
The new numerical work evaluates its actual nonspherical finite-cut data.
No quantum source, horizon replacement or stability is established.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp


def parameters(charge_radius=.2, weyl_separation=1.):
    if not 0 < charge_radius < weyl_separation:
        raise ValueError('Require 0<Q<D for this balanced branch.')
    k = weyl_separation/2
    mass = (math.sqrt(k*k+4*charge_radius*k)-k)/2
    a = math.sqrt(k*k-mass*mass)
    rp = mass+k
    field = mass/((k+a)*rp)  # (k-a)/(M*rplus), without cancellation.
    return dict(Q=charge_radius, D=weyl_separation, k=k, M=mass, a=a,
                rplus=rp, B=field, q=1., physical_A_factor=-1/(2*charge_radius))


def symbolic_fields():
    r, c, mass, a, field = sp.symbols('r c M a B', real=True)
    delta = r*r-2*mass*r-a*a
    sigma = r*r-a*a*c*c
    h = delta+(mass*mass+a*a)*(1-c*c)
    w = (r*r-a*a)**2+delta*a*a*(1-c*c)
    lapse = (delta+a*a*(1-c*c)+2*field*mass*r*a*(1-c*c)
             +field**2*w*(1-c*c)/4)/sigma
    meridian = lapse*sigma**2/h**2
    potential = -(2*mass*r*a+field*w/2)*(1-c*c)/(lapse*sigma)
    # The M=0 member is exactly Melvin, including its gauge orientation.
    rho2 = (r*r-a*a)*(1-c*c)
    melvin = 1+field**2*rho2/4
    assert sp.cancel(lapse.subs(mass, 0)-melvin) == 0
    assert sp.cancel(meridian.subs(mass, 0)-melvin) == 0
    assert sp.cancel(potential.subs(mass, 0)+field*rho2/(2*melvin)) == 0
    expr = [lapse, meridian, potential]
    values = expr+[sp.diff(f,r) for f in expr]+[sp.diff(f,c) for f in expr]
    cut_function = sp.lambdify((r,c,mass,a,field), values, 'numpy', cse=True)
    theta = sp.symbols('theta', real=True)
    metric = [-lapse**2, lapse**2*sigma**4/(h**3*delta),
              lapse**2*sigma**4/h**3, delta*(1-c*c)/lapse**2]
    # Actual derivatives of the published metric, not the field equations.
    bulk = [f.subs(c,sp.cos(theta)) for f in metric+[potential]]
    deriv = bulk+[sp.diff(f,r) for f in bulk]+[sp.diff(f,theta) for f in bulk]
    deriv += [sp.diff(f,r,2) for f in bulk]
    deriv += [sp.diff(f,r,theta) for f in bulk]
    deriv += [sp.diff(f,theta,2) for f in bulk]
    bulk_function = sp.lambdify((r,theta,mass,a,field), deriv, 'numpy', cse=True)
    return cut_function, bulk_function


def curvature_check(function, p, r, theta, mass=None, field=None):
    m = p['M'] if mass is None else mass
    b = p['B'] if field is None else field
    values = np.asarray(function(r,theta,m,p['a'],b),float).reshape(6,5)
    metric = values[0,:4]
    inverse = 1/metric
    dg = np.zeros((4,4,4)); ddg = np.zeros((4,4,4,4))
    for i,j in enumerate((1,2)):
        dg[j] = np.diag(values[1+i,:4])
    ddg[1,1] = np.diag(values[3,:4])
    ddg[1,2] = ddg[2,1] = np.diag(values[4,:4])
    ddg[2,2] = np.diag(values[5,:4])
    connection = np.zeros((4,4,4)); derivative = np.zeros((4,4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                bracket = dg[b,a,c]+dg[c,a,b]-dg[a,b,c]
                connection[a,b,c] = inverse[a]*bracket/2
                for d in (1,2):
                    derivative[d,a,b,c] = (inverse[a]*(ddg[d,b,a,c]+ddg[d,c,a,b]-ddg[d,a,b,c])
                        -inverse[a]**2*dg[d,a,a]*bracket)/2
    ricci = np.zeros((4,4))
    for a in range(4):
        for b in range(4):
            ricci[a,b] = sum(derivative[c,c,a,b]-derivative[b,c,a,c]
                +sum(connection[c,c,d]*connection[d,a,b]-connection[c,b,d]*connection[d,a,c]
                     for d in range(4)) for c in range(4))
    tensor = np.zeros((4,4))
    tensor[1,3] = values[1,4]; tensor[2,3] = values[2,4]
    tensor[3,1] = -tensor[1,3]; tensor[3,2] = -tensor[2,3]
    contracted = tensor @ np.diag(inverse) @ tensor.T
    f2 = np.sum(tensor*tensor*inverse[:,None]*inverse[None,:])
    source = 2*(contracted-np.diag(metric)*f2/4)
    frame = np.sqrt(np.abs(metric[:,None]*metric[None,:]))
    residual = (ricci-source)/frame
    reference = max(float(np.max(np.abs(source/frame))),1.)
    maxwell_terms = []
    for i,second in ((1,3),(2,5)):
        log_derivative = .5*np.sum(np.diag(dg[i])*inverse)-dg[i,i,i]*inverse[i]-dg[i,3,3]*inverse[3]
        first = values[i,4]
        maxwell_terms.extend([inverse[i]*inverse[3]*values[second,4],
                              inverse[i]*inverse[3]*first*log_derivative])
    return dict(r=r,theta=theta,M=m,B=b,ricci_covariant=ricci.tolist(),
                Maxwell_Ricci_source_covariant=source.tolist(),
                max_orthonormal_Einstein_residual=float(np.max(abs(residual))),
                relative_Einstein_residual=float(np.max(abs(residual)))/reference,
                scalar_curvature=float(np.sum(np.diag(ricci)*inverse)),
                Maxwell_phi_divergence=sum(maxwell_terms),
                relative_Maxwell_divergence=abs(sum(maxwell_terms))/max(sum(abs(v) for v in maxwell_terms),1.),
                metric_diagonal=metric.tolist())


def cut_profile(function,p,epsilon,mu):
    """Outward proper normal points away from the upper horizon endpoint."""
    if not 0 < epsilon < p['k']:
        raise ValueError('Cuts must be disjoint and surround only one endpoint.')
    mu = np.asarray(mu,float)
    if np.any(abs(mu)>1): raise ValueError('Local cosine must be in [-1,1].')
    k=p['k']; other=np.sqrt(epsilon**2+4*k*k+4*k*epsilon*mu)
    r=p['M']+(epsilon+other)/2; c=(other-epsilon)/(2*k)
    re=(1+(epsilon+2*k*mu)/other)/2
    ce=((epsilon+2*k*mu)/other-1)/(2*k)
    rmu=k*epsilon/other; cmu=epsilon/other
    raw=function(r,c,p['M'],p['a'],p['B'])
    n,h,ahat,nr,hr,ar,nc,hc,ac=[np.broadcast_to(np.asarray(v,float),mu.shape) for v in raw]
    ne=nr*re+nc*ce; he=hr*re+hc*ce
    ae=(ar*re+ac*ce)*p['physical_A_factor']
    amu=(ar*rmu+ac*cmu)*p['physical_A_factor']
    area=epsilon**2*h/n
    return dict(mu=mu,r=r,cos_theta=c,lapse=n,meridian_factor=h,
        Aphysical=ahat*p['physical_A_factor'],A_epsilon=ae,A_mu=amu,
        area_density_per_unit_solid_angle=area,
        meridional_radius=epsilon*h,azimuthal_radius=epsilon/n,
        K_time=ne/(h*n),K_meridional=(1/epsilon+he/h)/h,
        K_azimuthal=(1/epsilon-ne/n)/h,
        F_normal_phi=ae/h,B_normal=-amu/area,
        current_trace_normalization=np.sqrt(n*area))


def cut_summary(function,p,epsilon,order):
    mu,weights=np.polynomial.legendre.leggauss(order)
    z=cut_profile(function,p,epsilon,mu)
    area=2*math.pi*float(weights @ z['area_density_per_unit_solid_angle'])
    flux=-2*math.pi*float(weights @ z['A_mu'])
    poles=cut_profile(function,p,epsilon,np.array([-1.,1.]))
    pole_flux=2*math.pi*float(poles['Aphysical'][0]-poles['Aphysical'][1])
    return dict(epsilon=epsilon,angular_order=order,area=area,
        area_radius=math.sqrt(area/(4*math.pi)),flux=flux,pole_flux=pole_flux,
        lapse_min=float(np.min(z['lapse'])),lapse_max=float(np.max(z['lapse'])),
        meridional_radius_min=float(np.min(z['meridional_radius'])),
        meridional_radius_max=float(np.max(z['meridional_radius'])),
        azimuthal_radius_min=float(np.min(z['azimuthal_radius'])),
        azimuthal_radius_max=float(np.max(z['azimuthal_radius'])),
        max_extrinsic_angular_shear=float(np.max(abs(z['K_meridional']-z['K_azimuthal']))),
        K_time_min=float(np.min(z['K_time'])),K_time_max=float(np.max(z['K_time'])),
        pole_conical_factors=(poles['lapse']*poles['meridian_factor']).tolist())


def run():
    p=parameters();cut,bulk=symbolic_fields()
    bulk_rows=[curvature_check(bulk,p,r,t) for r,t in ((.72,.25),(.8,.9),(1.1,1.4),(2.,.6),(3.,1.2))]
    # Exact flat and Melvin limits independently exercise curvature construction.
    limits=[curvature_check(bulk,p,1.3,.8,mass=0.,field=f) for f in (0.,p['B'])]
    rows=[cut_summary(cut,p,e,n) for e in (.025,.05,.1) for n in (48,96)]
    profile=cut_profile(cut,p,.05,np.linspace(-1,1,257))
    compare=[]
    for e in (.025,.05,.1):
        a,b=[z for z in rows if z['epsilon']==e]
        compare.append(dict(epsilon=e,area_change=abs(a['area']-b['area']),flux_change=abs(a['flux']-b['flux'])))
    # Check normal derivatives without differentiating the symbolic formulas again.
    mu=np.array([-.8,-.1,.6]);step=2e-6
    z=cut_profile(cut,p,.05,mu)
    around=[cut_profile(cut,p,.05+j*step,mu) for j in (-2,-1,1,2)]
    fd=lambda key:(around[0][key]-8*around[1][key]+8*around[2][key]-around[3][key])/(12*step)
    derivative_error=max(float(np.max(abs(fd('Aphysical')-z['A_epsilon']))),
        float(np.max(abs(fd('lapse')/(z['meridian_factor']*z['lapse'])-z['K_time']))))
    balance=(1+p['M']**2/p['a']**2)**2/(1+p['B']*p['M']*p['rplus']/p['a'])**4
    charge=p['M']*p['rplus']/p['k']
    tests=dict(exact_symbolic_Melvin_limit=True,
        balanced_axis=abs(balance-1)<2e-13,
        charge_radius_calibration=abs(charge-p['Q'])<2e-14,
        bulk_Einstein_Maxwell=max(z['relative_Einstein_residual'] for z in bulk_rows+limits)<2e-9,
        source_free_Maxwell=max(z['relative_Maxwell_divergence'] for z in bulk_rows+limits)<2e-9,
        flat_curvature=max(abs(v) for row in limits[0]['ricci_covariant'] for v in row)<2e-10,
        trace_free_curvature=max(abs(z['scalar_curvature']) for z in bulk_rows+limits)<2e-9,
        metric_signature=all(z['metric_diagonal'][0]<0 and min(z['metric_diagonal'][1:])>0 for z in bulk_rows+limits),
        cut_flux=all(abs(z['flux']-2*math.pi)<2e-10 and abs(z['pole_flux']-2*math.pi)<2e-10 for z in rows),
        regular_cut_poles=max(abs(v-1) for z in rows for v in z['pole_conical_factors'])<2e-10,
        angular_quadrature=max(max(z['area_change'],z['flux_change']) for z in compare)<2e-10,
        normal_derivative_control=derivative_error<2e-7,
        positive_finite_cut_lapse=float(np.min(profile['lapse']))>0,
        distinct_angular_geometry=float(np.ptp(profile['lapse']))>1e-4)
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,parameters=p,
        bulk_checks=bulk_rows,flat_and_Melvin_controls=limits,cut_summaries=rows,
        selected_cut_profile={k:v.tolist() for k,v in profile.items()},
        angular_quadrature_sensitivity=compare,normal_derivative_error=derivative_error,
        axis_balance_factor=balance,charge_radius_reconstructed=charge,
        normalization=dict(action='(R-Fhat_ab Fhat^ab)/(16*pi*G)',
            project_gauge='Aphysical=-q*Ahat/(2Q), upper flux+2*pi*q, lower flux-2*pi*q.',
            required_coupling='4*pi*G/g_R^2=(2Q/q)^2; gauge rescaling is accompanied by this relation.',
            time='Killing time normalized to the asymptotic Melvin axis. There is no asymptotically flat infinity.',
            length='D=2k is Weyl coordinate separation, not proper horizon separation; the exact extremal horizons are infinitely far down their throats.',
            cut='Euclidean auxiliary Weyl radius epsilon about upper endpoint; epsilon is neither proper distance nor areal radius.',
            normal='n points away from the horizon toward exterior. Reverse sign for the exterior-domain inner normal.',
            spinor_trace='chi_Sigma=sqrt(N*dA/dOmega)*U_normal*Psi. Fixed round reference dOmega; actual angular frame/gauge transitions are still required.'),
        scope=['A known exact classical Einstein-Maxwell exterior is evaluated; this is not a new solution of physics.',
            'Axis balance needs a Melvin external field. It removes a conical strut, not black-hole horizons or dynamical instabilities.',
            'Finite cuts have angular lapse, metric, extrinsic curvature and gauge data. The spherical108 collar and old flat infinity map cannot be silently reused.',
            'Flux and bulk equations are checked numerically, with independent limiting cases; this is not a rigorous numerical error enclosure.',
            'No horizon excision and throat sewing, quantum state/source/current, self-consistent backreaction, matter preparation or time travel is established.',
            'The RN gate is a separate isolated comparator; these exact pair parameters are not the old flat model parameters under a coordinate rename.'],
        sources=['https://arxiv.org/pdf/hep-th/9906160 (metric8-10, balance12, charge16, near-horizon14-15)'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if not result['passed']: raise AssertionError(tests)
    Path(__file__).with_name('balanced_dihole_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,parameters=p,
        max_Einstein_error=max(z['relative_Einstein_residual'] for z in bulk_rows),
        max_Maxwell_error=max(z['relative_Maxwell_divergence'] for z in bulk_rows),
        selected_cut=next(z for z in rows if z['epsilon']==.05 and z['angular_order']==96),
        derivative_error=derivative_error),indent=2))
    return result


if __name__=='__main__':
    run()
