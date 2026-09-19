"""Untwisted Ori scalar image kernel: geometry and fixed-background stress.

Uses a=1/50, mu=-11251/11250, hbar=c=1. This is a source calculation on
the fixed curved geometry, not a solved semiclassical Einstein evolution.
"""
from pathlib import Path
import json
import mpmath as mp

OUT = Path(__file__).resolve().parent


def rosen(point, nu):
    T, x, y, z = map(mp.mpf, point)
    u = mp.exp(-z/2)
    r = [(1-v)/2 for v in nu]
    e = [u**v for v in r]
    V = -2*T/u-(r[0]*x*x+r[1]*y*y)/(2*u)
    return u, V, (x/e[0], y/e[1]), e[0]*e[1], [u**v/v for v in nu]


def image_geometry(point, other, winding, nu):
    shifted = list(other)
    shifted[3] += 2*mp.pi*winding
    u, V, X, E, J = rosen(point, nu)
    up, Vp, Xp, Ep, Jp = rosen(shifted, nu)
    if u == up:
        raise ValueError('This Jacobi formula requires distinct u; use its continuous limit otherwise')
    I = [J[i]-Jp[i] for i in range(2)]
    sigma = (u-up)*(-(V-Vp)+sum((X[i]-Xp[i])**2/(2*I[i]) for i in range(2)))
    delta_sqrt = abs(u-up)/mp.sqrt(E*Ep*I[0]*I[1])
    return sigma, delta_sqrt


def image_kernel(point, other, winding, nu):
    """Real nonsingular value; distributional future i0 is a separate condition."""
    sigma, delta_sqrt = image_geometry(point, other, winding, nu)
    return delta_sqrt/(8*mp.pi**2*sigma)


def mixed_image_derivatives(point, winding, nu):
    """Differentiate the full kernel before restricting to a chosen point."""
    point = list(map(mp.mpf, point))
    def derivative(i, j):
        def varied(a, b):
            p, q = point.copy(), point.copy()
            p[i], q[j] = a, b
            return image_kernel(p, q, winding, nu)
        return mp.diff(varied, (point[i], point[j]), (1, 1))
    return mp.matrix([[derivative(i, j) for j in range(4)] for i in range(4)])


def canonical_stress(point, mixed):
    T, x, y, z = map(mp.mpf, point)
    a, mu = mp.mpf(1)/50, -mp.mpf(11251)/11250
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = mp.matrix([[0, 0, 0, -1], [0, 1, 0, 0], [0, 0, 1, 0], [-1, 0, 0, F-T]])
    contracted = (T-F)*mixed[0, 0]-mixed[0, 3]-mixed[3, 0]+mixed[1, 1]+mixed[2, 2]
    return mixed-g*contracted/2


def axis_coefficients(n, nu):
    h = mp.pi*n/2
    s = mp.sinh(h)
    f = mp.sqrt(nu[0]*nu[1])/(64*mp.pi**2*s*mp.sqrt(mp.sinh(nu[0]*h)*mp.sinh(nu[1]*h)))
    L = mp.coth(h)+sum(v*mp.coth(v*h) for v in nu)/2
    Q = mp.csch(h)**2+sum(v*v*mp.csch(v*h)**2 for v in nu)/2
    a = f/mp.sinh(h)**2
    b = f*(L*mp.coth(h)+mp.csch(h)**2)/4
    cx, cy = [f*v/(4*s*mp.sinh(v*h)) for v in nu]
    d = f*(L*L+Q)/8
    j = a-2*b+cx+cy
    return dict(TT=a,Tz=b+j/2,xx=cx-j/2,yy=cy-j/2,zz=d+j/2,trace=-j,phi2=2*f)


def axis_stress(T, coefficients):
    a, b, cx, cy, d = [coefficients[k] for k in ('TT','Tz','xx','yy','zz')]
    return mp.matrix([[a/T**3,0,0,b/T**2],[0,cx/T**2,0,0],
                      [0,0,cy/T**2,0],[b/T**2,0,0,d/T]])


def response_checks(minimal, conformal, mu):
    """Necessary optical/observer diagnostics, not a coupled Einstein solution."""
    import sympy as sp
    from scipy.integrate import solve_ivp
    T,x,y,z = sp.symbols('T x y z', real=True)
    coordinates = (T,x,y,z)
    g = sp.Matrix([[0,0,0,-1],[0,1,0,0],[0,0,1,0],[-1,0,0,sp.Function('F')(x,y)-T]])
    inverse = g.inv()
    def connection(a,b,c):
        return sp.simplify(sum(inverse[a,d]*(sp.diff(g[d,c],coordinates[b])
            +sp.diff(g[d,b],coordinates[c])-sp.diff(g[b,c],coordinates[d])) for d in range(4))/2)
    assert g[0,0] == 0 and all(connection(a,0,0)==0 for a in range(4))
    assert all(connection(3,i,j)==0 for i in (1,2) for j in (1,2))
    s,s0,beta = sp.symbols('s s0 beta', positive=True)
    theta = beta*(s**-2-s0**-2)/2
    area_change = beta*(s0-s)**2/(2*s*s0**2)
    assert sp.simplify(-sp.diff(theta,s)-beta/s**3)==0
    assert sp.simplify(-sp.diff(area_change,s)-theta)==0
    assert theta.subs(s,s0)==0 and area_change.subs(s,s0)==0
    delta,q,A,B,D = sp.symbols('delta q A B D', positive=True)
    P = A*(2*q+delta)**2-2*B*delta*(2*q+delta)+D*delta**2
    t = -q*delta-delta**2/4
    uT,uz = q+delta/2,2/(4*q+delta)
    rho = A*uT**2/t**3+2*B*uT*uz/t**2+D*uz**2/t
    assert sp.simplify(rho+16*P/(delta**3*(4*q+delta)**3))==0
    assert sp.limit(delta**3*rho,delta,0)==-A/q
    denominator = delta**3*(4*q+delta)
    derivative_polynomial = (48*A*q**3+16*q**2*(3*A-2*B)*delta
        +4*q*(4*A-5*B+D)*delta**2+2*(A-2*B+D)*delta**3)
    assert sp.expand((sp.diff(denominator,delta)*P-denominator*sp.diff(P,delta))
                     -delta**2*derivative_polynomial)==0
    optical = []
    for epsilon in map(mp.mpf,('.001','.002','.01')):
        x0 = mp.sqrt(2*epsilon)
        def radius(r):
            return x0*mp.sqrt(r)*(mp.besselk(2,x0)*mp.besseli(1,x0/mp.sqrt(r))
                +mp.besseli(2,x0)*mp.besselk(1,x0/mp.sqrt(r)))
        assert abs(radius(1)-1)<mp.mpf('1e-50')
        assert abs(mp.diff(radius,mp.mpf(1)))<mp.mpf('1e-50')
        r = mp.mpf('.1')
        assert abs(mp.diff(radius,r,2)-epsilon*radius(r)/(2*r**3))<mp.mpf('1e-48')
        solved = solve_ivp(lambda r,b:(b[1],float(epsilon)*b[0]/(2*r**3)),
            (1.,.1),(1.,0.),method='DOP853',rtol=2e-12,atol=2e-13,max_step=.01)
        assert solved.success and solved.t[-1]==.1
        error = abs(solved.y[0,-1]-float(radius(r)))
        assert error<1e-10
        linear = 1+epsilon*(1-r)**2/(2*r)
        optical.append(dict(epsilon=float(epsilon),r=float(r),linear_area=float(linear),
            prescribed_area=float(radius(r)**2),radius_ODE_error=error,
            area_remainder=float(radius(r)**2-linear)))
    assert 3.9<optical[1]['area_remainder']/optical[0]['area_remainder']<4.1
    records = {}
    for name,c in (('minimal',minimal),('conformal',conformal)):
        a,b,d = [c[k] for k in ('TT','Tz','zz')]
        assert a>0 and a*d-b*b>0
        assert min(3*a-2*b,4*a-5*b+d,a-2*b+d)>0
        q = 1/mp.sqrt(2)
        def ratio(delta,eta):
            p = a*(2*q+delta)**2-2*b*delta*(2*q+delta)+d*delta**2
            return 4*eta*p/(abs(mu)*delta**3*(4*q+delta))
        events = []
        for eta in map(mp.mpf,('1','1e-6','1e-12','1e-60')):
            leading = (4*eta*a*q/abs(mu))**(mp.mpf(1)/3)
            log_delta = mp.findroot(lambda v:mp.log(ratio(mp.exp(v),eta)),mp.log(leading))
            delta = mp.exp(log_delta)
            assert abs(ratio(delta,eta)-1)<mp.mpf('1e-48')
            s = q*delta+delta**2/4
            velocity = mp.matrix([q+delta/2,0,0,2/(4*q+delta)])
            measured = (velocity.T*axis_stress(-s,c)*velocity)[0]
            target = mu/eta*velocity[3]**2
            assert abs(measured/target-1)<mp.mpf('1e-48')
            change = eta*a*(1-s)**2/(2*s)
            events.append(dict(eta=mp.nstr(eta,8),delta=mp.nstr(delta,35),s=mp.nstr(s,35),
                exact_ratio=mp.nstr(ratio(delta,eta),8),leading_delta=mp.nstr(leading,35),
                leading_relative_error=mp.nstr(abs(leading/delta-1),12),
                null_area_first_order_change=mp.nstr(change,20),optical_chi=mp.nstr(eta*a/s,20)))
        records[name]=dict(observer_q=mp.nstr(q,35),events=events,
            positivity_determinant=mp.nstr(a*d-b*b,25),
            proper_length_asymptotic_coefficient=mp.nstr((32*mp.pi*a*q/abs(mu))**(mp.mpf(1)/3),25))
    return dict(
        scope='Leading Einstein-term optical response conditional on regular aligned rank-one classical support, plus exact quantum stress measurements on the prescribed background. Neither is a self-consistent backreacted solution',
        null_geometry='k=partial_T is null and affine, k_flat=-dz is closed, and transverse optical matrix is zero. All checked symbolically',
        support_response='delta[rho*(l.k)^2]=0 when (l.k)_0=0, for regular first-order variations of aligned rank-one null support; other support models need their own response',
        first_order_optics='beta=eta*A_xi, s=-T; initial theta1(-s0)=0 and area ratio1. Then theta1=beta*(1/s^2-1/s0^2)/2; area/area0=1+beta*(s0-s)^2/(2*s*s0^2)+O(beta^2) on fixed interior intervals',
        optical_control='For 0<s<=s0, chi=abs(beta)/s bounds the displayed first-order area change by chi/2. theta1^2/(2*abs(deltaRkk))=chi*(1-s^2/s0^2)^2/8. This is not a bound on full metric, unknown shear or state-response errors',
        prescribed_optical_ODE='For positive epsilon=beta/s0, radius b obeys b_rr=epsilon*b/(2*r^3), b(1)=1,b_r(1)=0. Bessel solution and numerical integration compared; imposed zero shear and fixed Ricci profile are not a solved gravity system',
        prescribed_optical_checks=optical,
        observer_formula='delta=-lambda>0; P=A*(2q+delta)^2-2B*delta*(2q+delta)+D*delta^2; rho_quantum_model=-16*P/[delta^3*(4q+delta)^3]',
        observer_physical_units='rho_quantum=(hbar*c/ell^4)*rho_model; rho_classical=c^4/(8*pi*G*ell^2)*4mu/(4q+delta)^2; eta=8*pi*l_P^2/ell^2',
        source_ratio='R=4*eta*P/[abs(mu)*delta^3*(4q+delta)]; R=1 has a unique positive root for both tested couplings. It measures source comparability for this observer, not a universal validity threshold',
        observer_asymptotic='rho_quantum_model~-A/(q*delta^3), R~4*eta*A*q/(abs(mu)*delta^3). At R=1, ell*delta~(32*pi*A*q/abs(mu))^(1/3)*l_P^(2/3)*ell^(1/3)',
        observer_cases=records,
        remaining='Full metric/constraint evolution, support conservation, metric-dependent state response and stress correlations, shear, and higher-curvature gravitational terms are uncomputed. Small optical area alone cannot establish stability; fixed-background divergence is not an injury or horizon-destruction prediction',
        sources=['https://arxiv.org/pdf/math/9909158 (null optical equation 2.4)',
                 'https://arxiv.org/pdf/gr-qc/0209075 (semiclassical response and validity, section III)'])


def abiotic_probe_checks():
    """Unbound tracer-pair size test on the existing classical geometry."""
    from scipy.integrate import solve_ivp
    q = 1/mp.sqrt(2)
    directions = {}
    for name,nu in (('x',mp.mpf(227)/75),('y',mp.mpf(223)/75)):
        plus,minus = (1+nu)/2,(1-nu)/2
        strength = (nu**2-1)/4
        def separation(lam):
            r = (4*q-lam)/(4*q+1)
            return (plus*r**minus-minus*r**plus)/nu
        assert abs(separation(-1)-1)<mp.mpf('1e-50')
        assert abs(mp.diff(separation,mp.mpf(-1)))<mp.mpf('1e-50')
        errors = [abs(mp.diff(separation,lam,2)
            -strength*separation(lam)/(4*q-lam)**2)
            for lam in map(mp.mpf,('-.7','-.1','0'))]
        assert max(errors)<mp.mpf('1e-48')
        integration = solve_ivp(lambda lam,state:(state[1],
            float(strength)*state[0]/(4*float(q)-lam)**2),
            (-1.,0.),(1.,0.),method='DOP853',rtol=2e-12,atol=2e-13,max_step=.02)
        assert integration.success and integration.t[-1]==0
        ode_error = abs(integration.y[0,-1]-float(separation(0)))
        assert ode_error<1e-11
        growth = separation(0)-1
        assert growth>0
        directions[name]=dict(nu=mp.nstr(nu,25),
            classical_boundary_separation_factor=mp.nstr(separation(0),35),
            pre_boundary_factor_at_lambda_minus001=mp.nstr(separation(mp.mpf('-.01')),35),
            change_for_initial_1_micrometre_in_nanometres=mp.nstr(1000*growth,20),
            tidal_eigenvalue_at_classical_boundary=mp.nstr(-strength/(4*q)**2,25),
            Jacobi_residual_max=mp.nstr(max(errors),8),independent_ODE_error=ode_error)
    return dict(
        payload='Neutral, nonspinning passive classical tracer. Linear Jacobi proxy for initially comoving unbound tracers; the1 micrometre numbers are infinitesimal-separation rescalings until ell is specified. Require separation small against local curvature/variation scales, e.g. d0/ell much less than1 here. Not a bonded grain or organic tissue model',
        geometry='Same critical untwisted prescribed classical metric; q=1/sqrt2 is the positive magnitude of the specific Killing momentum, with p_z/m=-q in model units. Same initial four-velocity is held fixed when test mass changes',
        trajectory='T=q*lambda-lambda^2/4, z=-2log(1-lambda/(4q)), x=y=0. Point-particle mass cancels from the geodesic equation when self-force and nongravitational forces are neglected',
        parallel_frame='The transverse unit vectors partial_x,partial_y are parallel along the central trajectory; E_x=2*(mu-a)/(4q-lambda)^2, E_y=2*(mu+a)/(4q-lambda)^2, longitudinal eigenvalue0',
        separation_formula='Start lambda=-1 with separation ratio1 and derivative0. For nu_i, p_plus/minus=(1+/-nu_i)/2, r=(4q-lambda)/(4q+1): f_i=(p_plus*r^p_minus-p_minus*r^p_plus)/nu_i',
        Jacobi_equation='f_i_second=-(E_i)*f_i=(nu_i^2-1)*f_i/[4*(4q-lambda)^2]',
        directions=directions,
        size_scaling='Absolute differential displacement and tidal acceleration scale linearly with initial separation; fractional unbound-tracer separation change is independent of that size. Physical delta_acceleration_i=-c^2*E_i*separation_i/ell^2',
        time_scale='Physical proper time is ell*lambda/c; ell is the model length normalization, not automatically an aperture or grain size',
        boundary_scope='lambda0 is only the smooth classical boundary limit. The scalar state/backreaction already obstructs treating this as a validated physical crossing; no calendar-date map is supplied',
        remaining='Internal bonding, elastic stress, heating, electromagnetic interactions, quantum wavepacket spread and self-force are omitted. Scaling a tracer result to a grain, molecule, animal or human requires separate tests. Decreasing payload mass or size does not change the source-matching or horizon-regularity conditions',
        sources=['https://arxiv.org/abs/1102.0529 (point-particle/small-body approximation and self-force)',
                 'https://www.preposterousuniverse.com/wp-content/uploads/grnotes-three.pdf (geodesic deviation, equation3.113)'])


def run_checks():
    with mp.workdps(60):
        nu = (mp.mpf(227)/75, mp.mpf(223)/75)
        pairs = []
        for n in (1, 2):
            for eps in map(mp.mpf, ('.1', '.01', '.001')):
                scale = mp.exp(-mp.pi*n)
                p = [-eps/2, 0, 0, 0]
                q = [-(1+eps)*scale*eps/2, 0, 0, -2*mp.log(1+eps)]
                local, _ = image_geometry(p, q, 0, nu)
                winding, amplitude = image_geometry(p, q, n, nu)
                assert local > 0 and abs(local-eps**2*(1-scale)) < mp.mpf('1e-55')
                assert abs(winding) < mp.mpf('1e-55') and amplitude > 0
                pairs.append(dict(winding=n,epsilon=float(eps),
                    point=[float(v) for v in p],other=[float(v) for v in q],
                    local_sigma=float(local),winding_sigma_residual=float(winding),
                    singularity_prefactor=float(amplitude/(8*mp.pi**2))))
        coefficients = [axis_coefficients(n,nu) for n in range(1,31)]
        trace_identity_errors = []
        for n, row in enumerate(coefficients, 1):
            h = mp.pi*n/2
            positive_trace = row['phi2']/2*(mp.mpf('.5')+sum(
                v*(mp.cosh(h)*mp.cosh(v*h)-1)/(4*mp.sinh(h)*mp.sinh(v*h)) for v in nu))
            trace_identity_errors.append(abs(row['trace']-positive_trace))
        assert max(trace_identity_errors) < mp.mpf('1e-50')
        summed = {key:mp.fsum(row[key] for row in coefficients) for key in coefficients[0]}
        phi_coefficient, tt_coefficient = summed['phi2'], summed['TT']
        point = [mp.mpf(-1), 0, 0, 0]
        mixed = mixed_image_derivatives(point, 1, nu)+mixed_image_derivatives(point, -1, nu)
        stress = canonical_stress(point, mixed)
        assert mp.norm(mixed-mixed.T) < mp.mpf('1e-50')
        tt_error = abs(stress[0, 0]+coefficients[0]['TT'])
        assert tt_error < mp.mpf('1e-50') and stress[0, 0] < 0
        probe = list(map(mp.mpf,('-0.7','0','0','0.3')))
        direct = sum((mixed_image_derivatives(probe,n,nu) for n in (-4,-3,-2,-1,1,2,3,4)),mp.zeros(4))
        four = {key:mp.fsum(row[key] for row in coefficients[:4]) for key in coefficients[0]}
        full_error = mp.norm(canonical_stress(probe,direct)-axis_stress(probe[0],four))
        assert full_error < mp.mpf('1e-48')
        kg_residuals = []
        p = list(map(mp.mpf,('-.7','.12','-.17','.2')))
        q = list(map(mp.mpf,('-.9','-.07','.09','-.11')))
        for n in (-1,1):
            for first in (True,False):
                point = p if first else q
                fn = lambda *args: image_kernel(args,q,n,nu) if first else image_kernel(p,args,n,nu)
                def derivative(orders):
                    return mp.diff(fn,tuple(point),orders)
                T,x,y,z = point
                a,mu = mp.mpf(1)/50,-mp.mpf(11251)/11250
                F = ((a-mu)*x*x-(a+mu)*y*y)/2
                residual = ((T-F)*derivative((2,0,0,0))+derivative((1,0,0,0))
                    -2*derivative((1,0,0,1))+derivative((0,2,0,0))+derivative((0,0,2,0)))
                kg_residuals.append(abs(residual))
        assert max(kg_residuals) < mp.mpf('1e-48')
        h = mp.pi/2
        L = mp.coth(h)+sum(v*mp.coth(v*h) for v in nu)/2
        Q = mp.csch(h)**2+sum(v*v*mp.csch(v*h)**2 for v in nu)/2
        aa,bb = mp.csch(h)**2,(L*mp.coth(h)+mp.csch(h)**2)/4
        cc = [v/(4*mp.sinh(h)*mp.sinh(v*h)) for v in nu]
        dd = (L*L+Q)/8
        jj = aa+2*bb+sum(cc)
        multiplier = max(aa,bb+jj/2,cc[0]+jj/2,cc[1]+jj/2,dd+jj/2,jj,2)
        N = 30; alpha = (sum(nu)+2)/4
        correction = ((1-mp.exp(-mp.pi*(N+1)))
            *mp.sqrt((1-mp.exp(-nu[0]*mp.pi*(N+1)))*(1-mp.exp(-nu[1]*mp.pi*(N+1)))))
        f_tail = (mp.sqrt(nu[0]*nu[1])/(16*mp.pi**2)
            *mp.exp(-mp.pi*alpha*(N+1))/(correction*(1-mp.exp(-mp.pi*alpha))))
        coefficient_tail = multiplier*f_tail
        assert summed['trace'] > coefficient_tail and summed['TT'] > coefficient_tail
        conformal_tt = tt_coefficient+phi_coefficient/3
        assert conformal_tt > 0
        transverse_Q = []
        q_hessian_errors = []
        for i,v in enumerate(nu):
            weights = []
            for n,row in enumerate(coefficients,1):
                h = mp.pi*n/2
                weights.append(row['phi2']*(v*(mp.cosh(h)*mp.cosh(v*h)-1)
                    /(mp.sinh(h)*mp.sinh(v*h))-1)/4)
            transverse_Q.append(mp.fsum(weights))
            def diagonal_Q(x):
                p = [mp.mpf(-1),mp.mpf(0),mp.mpf(0),mp.mpf(0)]
                p[i+1] = x
                return image_kernel(p,p,1,nu)+image_kernel(p,p,-1,nu)
            q_hessian_errors.append(abs(mp.diff(diagonal_Q,mp.mpf(0),2)+weights[0]))
        assert max(q_hessian_errors) < mp.mpf('1e-50')
        Hx,Hy = transverse_Q
        mu = -mp.mpf(11251)/11250
        improvement = dict(TT=2*phi_coefficient,Tz=3*phi_coefficient/2+Hx+Hy,
            xx=-phi_coefficient-Hy,yy=-phi_coefficient-Hx,
            zz=(mp.mpf('1.5')-mu)*phi_coefficient+Hx+Hy)
        conformal = {k:summed[k]+v/6 for k,v in improvement.items()}
        def null_polynomial(c,q,direction):
            return -c['TT']*(1+q)**2/4+c['Tz']*(1+q)+c[direction]*q-c['zz']
        minimal_witness = null_polynomial(summed,mp.mpf(4),'yy')
        assert minimal_witness > 0
        direction = max(('xx','yy'),key=lambda k:conformal[k])
        linear = conformal['Tz']+conformal[direction]-conformal['TT']/2
        maximizing_q = max(mp.mpf(0),2*linear/conformal['TT'])
        conformal_max = null_polynomial(conformal,maximizing_q,direction)
        assert conformal['TT'] > 0 and conformal_max < 0
        assert abs(-improvement['TT']/4+improvement['Tz']-improvement['zz']
                   -(mu-mp.mpf('.5'))*phi_coefficient) < mp.mpf('1e-50')
        planck_factor = mp.sqrt(8*mp.pi*(-conformal_max)/abs(mu))
        return dict(
            model='Massless scalar periodic image state on the fixed untwisted negative-mu Ori geometry; original coordinate order T,x,y,z; hbar=c=1',
            parameters=dict(a='1/50',mu='-11251/11250',nu=['227/75','223/75']),
            state_scope='Analytically constructed quasifree Hadamard state only on the pre-horizon domain D/Gamma. Correct future-oriented cover modes, infrared distribution estimate, full regulated kernel/CCR, compact C-infinity image tails and Fejer positivity independently reviewed; proof in RESEARCH.md. This script checks finite nonsingular geometry and stress formulas, not that distributional proof',
            cover_prescription='Cover W=sqrtDelta/(8*pi^2*sigma_+) with the original future orientation and no added smooth bisolution. Its massless Hadamard tail and subtracted cover stress vanish. Curvature-squared ambiguity tensors vanish by null alignment; gravitational couplings are fixed. Under this prescription the nonidentity stress is the full renormalized stress on D/Gamma',
            winding_null_pairs=pairs,
            winding_interpretation='Locally spacelike axis pairs approach the same base point, but one nonzero deck image is null related. Its connecting axis null geodesic stays inT<0. Nonzero kernel amplitude cannot be removed by the local zero-winding subtraction',
            scalar_square_coefficient=mp.nstr(phi_coefficient,45),
            scalar_square_formula='phi2=coefficient/(-T) on the axis in D/Gamma under the stated cover subtraction',
            TT_coefficient=mp.nstr(tt_coefficient,45),
            TT_formula='T_TT=C_TT/T^3<0 forT<0; C_TT=sum_n>=1 sqrtDelta_n/[64*pi^2*sinh(n*pi/2)^4]>0',
            full_axis_coefficients={key:mp.nstr(value,45) for key,value in summed.items()},
            full_axis_tensor='Covariant(T,x,y,z): [[a/T^3,0,0,b/T^2],[0,cx/T^2,0,0],[0,0,cy/T^2,0],[b/T^2,0,0,d/T]]. Coefficients keysTT,Tz,xx,yy,zz; trace=trace_coefficient/T^2',
            coefficient_formula='Forh=n*pi/2,s=sinhh: f=sqrt(nux*nuy)/(64*pi^2*s*sqrt(sinh(nux*h)*sinh(nuy*h))); L=cothh+sum(nui*coth(nui*h))/2; Q=cschh^2+sum(nui^2*csch(nui*h)^2)/2. Raw mixed coefficients a=f*cschh^2,b=f*(L*cothh+cschh^2)/4,ci=f*nui/(4*s*sinh(nui*h)),d=f*(L^2+Q)/8; j=a-2b+cx+cy. Stress coefficients a,b+j/2,cx-j/2,cy-j/2,d+j/2; trace=-j',
            coefficient_tail_bound=mp.nstr(coefficient_tail,12),
            tail_scope='All seven coefficient sums after30 terms, bounded using monotone hyperbolic factors and a geometric f tail. Does not include arithmetic rounding; each tensor component also requires its displayed power ofT',
            nonminimal_coupling=dict(convention='T_ab(xi)=T_ab(0)+xi*(G_ab*Q+g_ab*BoxQ-nabla_a nabla_b Q), Q=<phi^2>; same field equation becauseR=0',
                TT_formula='[a+2*xi*Cphi]/T^3',trace_formula='(1-6*xi)*minimal_trace',
                minimal_trace_positive_identity='Each paired coefficient -j=f*[1/2+sum_i nui*(cosh(h)*cosh(nui*h)-1)/(4*sinh(h)*sinh(nui*h))]>0 for h,nui>0',
                xi_required_for_zero_trace='1/6',conformal_TT_coefficient=mp.nstr(conformal_tt,45),
                conclusion='For this image state, zero trace forcesxi1/6 and leavesT_TT<0. No realxi gives the required negative-null tensor, assuming the stated cover/renormalization prescription'),
            positive_trace_identity_max_error=mp.nstr(max(trace_identity_errors),8),
            transverse_scalar_square=dict(Hx=mp.nstr(Hx,45),Hy=mp.nstr(Hy,45),
                formula='On axis Q_ii=-H_i/T^2; H_i=sum_n phi2_n*[nui*(coshh*cosh(nui*h)-1)/(sinhh*sinh(nui*h))-1]/4',
                independent_first_pair_hessian_error=mp.nstr(max(q_hessian_errors),8)),
            curvature_coupling_coefficient_increments={k:mp.nstr(v,45) for k,v in improvement.items()},
            NEC_compensation=dict(
                null_family='s=-T>0; k=(s*(1+Rx^2+Ry^2)/2,sqrt(s)*Rx,sqrt(s)*Ry,1), so g(k,k)=0; also test kz=0 separately',
                contraction='s*Tscalar(k,k)=-A*(1+q)^2/4+B*(1+q)+Cx*Rx^2+Cy*Ry^2-D, q=Rx^2+Ry^2',
                minimal_witness_Rx_Ry=[0,2],minimal_positive_contraction=mp.nstr(minimal_witness,45),
                minimal_result='At every s>0 this direction has positive scalar null stress but negative target mu/eta, so the compensator must violate the NEC there for every finite positive eta',
                conformal_axis_coefficients={k:mp.nstr(v,45) for k,v in conformal.items()},
                conformal_maximum_scaled_null_contraction=mp.nstr(conformal_max,45),
                conformal_maximum_direction=direction,conformal_maximum_R_squared=mp.nstr(maximizing_q,45),
                conformal_local_condition='All axis null contractions of the required compensator are nonnegative iff abs(mu)/eta <= -conformal_maximum_scaled_null_contraction/s. Necessary local matter test only; does not construct a source',
                physical_normalization='For ds_phys^2=ell^2 ds_model^2, eta=8*pi*l_P^2/ell^2 with l_P^2=hbar*G/c^3',
                conformal_ell_bound_in_Planck_lengths_at_s1=mp.nstr(planck_factor,30),
                length_bound_scope='Local NEC condition implies ell/l_P <= this coefficient/sqrt(s). This formal scale bound does not establish semiclassical control or human passage',
                complete_domain_witness='With Rx=Ry=0, C0(xi)=-A/4+B-D+xi*(mu-1/2)*Cphi. Compensator contraction mu/eta-C0/s is negative whenever s>eta*max(0,-C0)/abs(mu)',
                complete_domain_exclusion='For every fixed finite xi and eta>0 finite, take Rx=Ry=0. Scalar contraction is finite C0(xi)/s while target is mu/eta<0. As s grows unbounded in D, required compensator becomes negative. Thus no NEC-obeying compensator repairs this state on the entire D; finite intervals, altered early states, other geometries and backreaction are not excluded'),
            first_order_response=response_checks(summed,conformal,mu),
            abiotic_test_probe=abiotic_probe_checks(),
            independent_four_image_pair_full_tensor_error=mp.nstr(full_error,8),
            generic_full_kernel_KG_residuals=[mp.nstr(v,8) for v in kg_residuals],
            independent_first_image_pair_TT_error=mp.nstr(tt_error,8),
            first_image_pair_stress_at_T_minus1=[[mp.nstr(stress[i,j],35) for j in range(4)] for i in range(4)],
            target_comparison='Required Ori source mu/kappa*dz^2 has T_TT=0. This scalar state alone cannot equal that target for any real curvature coupling xi under the stated prescription. Other states, additional sources and changed geometry are not excluded',
            limitations='Analytic state construction holds strictly before the horizon on the prescribed geometry; it is not experimental preparation or a solved source/geometry pair. Gravitational backreaction, other states, additional matter, and stable traversal have not been solved',
            sources=['https://arxiv.org/pdf/1202.0540 (Jacobi/worldfunction framework)',
                     'https://arxiv.org/pdf/2309.02252 (plane-wave mode/Wightman method)',
                     'https://arxiv.org/pdf/1401.2026 (free-field state requirements)',
                     'https://arxiv.org/pdf/gr-qc/0404074 (scalar stress convention)',
                     'https://arxiv.org/pdf/gr-qc/9603012 (horizon quantum-regularity obstruction)'])


if __name__ == '__main__':
    result = run_checks()
    (OUT/'ori_quantum_source_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
