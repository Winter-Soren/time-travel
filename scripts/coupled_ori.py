"""Exact curved Ori with null scalars or Maxwell radiation, plus a spectator test.

Local verification/specialization of Ori's null-fluid family; no novelty claim.
Dimensionless coordinates,
signature -+++, positive coupling kappa. The support fields AND geometry solve
Einstein's equations. The additional spectator is linear on that solution;
its nonlinear gravitational response is NOT evolved here.
Run: .venv/Scripts/python.exe coupled_ori.py
"""
from project_paths import project_path
from pathlib import Path
import json
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parent


def geometry(g, coords):
    gi = g.inv()
    gamma = [[[sp.simplify(sum(gi[i, l]*(sp.diff(g[l, j], coords[h]) +
        sp.diff(g[l, h], coords[j])-sp.diff(g[j, h], coords[l]))/2
        for l in range(4))) for h in range(4)] for j in range(4)] for i in range(4)]
    ricci = sp.Matrix(4, 4, lambda i, j: sp.simplify(sum(
        sp.diff(gamma[h][i][j], coords[h])-sp.diff(gamma[h][i][h], coords[j]) +
        sum(gamma[h][h][l]*gamma[l][i][j]-gamma[h][j][l]*gamma[l][i][h]
            for l in range(4)) for h in range(4))))
    return gi, gamma, ricci


def symbolic_checks():
    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    a, mu, e, b, d, n = sp.symbols('a mu e b d n', real=True)
    eps, k, kappa = sp.symbols('eps k kappa', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, 0],
                   [0, 0, 1, 0], [-1, 0, 0, F-T]])
    gi, gamma, ricci = geometry(g, coords)
    fields = [eps*sp.cos(k*z), eps*sp.sin(k*z)]
    source = sp.Matrix(4, 4, lambda i, j: sp.trigsimp(kappa*sum(
        sp.diff(field, coords[i])*sp.diff(field, coords[j]) for field in fields)))
    residual = (ricci-source).subs(mu, kappa*eps**2*k**2).applyfunc(sp.simplify)
    assert residual == sp.zeros(4)
    assert g.det() == -1

    def box(field):
        return sp.simplify(sum(sp.diff(gi[i, j]*sp.diff(field, coords[j]), coords[i])
                              for i in range(4) for j in range(4)))

    def norm(field):
        return sp.expand(sum(gi[i, j]*sp.diff(field, coords[i])*sp.diff(field, coords[j])
                             for i in range(4) for j in range(4)))

    assert all(box(field) == 0 and norm(field) == 0 for field in fields)
    # Scalar invariants can vanish even when tidal curvature does not.
    riemann = {}
    for i in range(4):
        for j in range(4):
            for h in range(4):
                for l in range(4):
                    value = sp.simplify(sum(g[i, m]*(
                        sp.diff(gamma[m][j][l], coords[h])-
                        sp.diff(gamma[m][j][h], coords[l]) + sum(
                            gamma[m][h][v]*gamma[v][j][l]-
                            gamma[m][l][v]*gamma[v][j][h] for v in range(4)))
                        for m in range(4)))
                    if value != 0:
                        riemann[i, j, h, l] = value
    kretschmann = sp.simplify(sum(v*w*math.prod(gi[i, j] for i, j in zip(p, q))
        for p, v in riemann.items() for q, w in riemann.items()))
    assert kretschmann == 0
    assert sp.simplify(riemann[1, 3, 1, 3]+(a-mu)/2) == 0
    assert sp.simplify(riemann[2, 3, 2, 3]-(a+mu)/2) == 0

    t = sp.symbols('t', real=True)
    clock = T-F+e*(x*x+y*y)
    clock_norm = sp.expand(norm(clock).subs(T, t+F-e*(x*x+y*y)))
    expected_clock = t+((2*e-a+mu)**2-e)*x*x+((2*e+a+mu)**2-e)*y*y
    assert sp.simplify(clock_norm-expected_clock) == 0

    # Exact transverse closure for the independent spectator, not an axis ODE.
    S = T+b*x*x+d*y*y
    conditions = [4*b*b-b-(a-mu)/2, 4*d*d-d+(a+mu)/2]
    basis = sp.groebner(conditions, b, d, domain=sp.QQ_I.frac_field(a, mu))
    assert all(basis.reduce(c)[1] == 0 for c in sp.Poly(norm(S)-S, x, y, T).coeffs())
    assert box(S) == 1+2*b+2*d
    p = -2*(b+d)+2*sp.I*n
    loggrad = [p*sp.diff(S, c)/S+(sp.I*n if i == 3 else 0)
               for i, c in enumerate(coords)]
    residual_kg = sp.together(sum(sp.diff(gi[i, j], coords[i])*loggrad[j] +
        gi[i, j]*(sp.diff(loggrad[j], coords[i])+loggrad[i]*loggrad[j])
        for i in range(4) for j in range(4))).as_numer_denom()[0]
    assert all(sp.simplify(basis.reduce(c)[1]) == 0
               for c in sp.Poly(residual_kg, x, y, T, n).coeffs())
    return dict(all_16_Einstein_scalar_residuals_zero=True,
        both_support_wave_equations_zero=True, determinant=-1,
        Ricci_nonzero={'zz': str(ricci[3, 3])}, Kretschmann=str(kretschmann),
        nonzero_tidal_components={'xzxz': str(riemann[1, 3, 1, 3]),
                                 'yzyz': str(riemann[2, 3, 2, 3])},
        time_function_gradient_squared=str(clock_norm),
        spectator_full_4D_wave_residual_zero=True,
        support_intensity='mu = kappa * epsilon^2 * k^2, k integer for period 2*pi',
        regular_initial_region='e=1/16; a>=0, mu>=0, a+mu<1/8; t<0',
        central_loop_squared_norm='-t; smooth null crossing in this exact solution',
        support_observer_energy='epsilon^2*k^2/(2+T), x=y=0, -1<=T<0')


def maxwell_checks():
    """Verify an alternative ordinary electromagnetic source, not added matter.

    Signature -+++, c=8*pi*G=1, rationalized Maxwell normalization. This is a
    local specialization of the Einstein-Maxwell plane-wave family; see
    Balakin et al. (2000), doi:10.1088/0264-9381/17/24/306, section 3.1.
    No preparation mechanism, generic stability or originality is implied.
    """
    t, x, y, z = coords = sp.symbols('t x y z', real=True)
    a, e, amplitude = sp.symbols('a e A', real=True)
    k = sp.symbols('k', integer=True, positive=True)
    P, Q = 2*e-a, 2*e+a
    mu = amplitude**2
    H = (e-mu/2)*(x*x+y*y)-t
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, P*x],
                   [0, 0, 1, Q*y], [-1, P*x, Q*y, H]])
    gi, _, ricci = geometry(g, coords)
    assert g.det() == -1
    potential = sp.Matrix([0, amplitude*sp.sin(k*z)/k,
                           -amplitude*sp.cos(k*z)/k, 0])
    field = sp.Matrix(4, 4, lambda i, j:
        sp.diff(potential[j], coords[i])-sp.diff(potential[i], coords[j]))
    raised = (gi*field*gi).applyfunc(sp.trigsimp)
    expected_raised = sp.zeros(4)
    expected_raised[0, 1] = -amplitude*sp.cos(k*z)
    expected_raised[0, 2] = -amplitude*sp.sin(k*z)
    expected_raised -= expected_raised.T
    assert raised == expected_raised
    assert (potential.subs(z, z+2*sp.pi)-potential).applyfunc(sp.simplify) == sp.zeros(4, 1)
    # sqrt(-g)=1; check both source-free Maxwell equations with the actual metric.
    closure = [sp.simplify(sp.diff(field[j, h], coords[i]) +
        sp.diff(field[h, i], coords[j])+sp.diff(field[i, j], coords[h]))
        for i in range(4) for j in range(i+1, 4) for h in range(j+1, 4)]
    divergence = [sp.simplify(sum(sp.diff(raised[i, j], coords[i])
        for i in range(4))) for j in range(4)]
    invariant = sp.trigsimp(sum(field[i, j]*raised[i, j]
        for i in range(4) for j in range(4)))
    wedge = sp.simplify(field[0, 1]*field[2, 3]-field[0, 2]*field[1, 3]
                        + field[0, 3]*field[1, 2])
    assert all(v == 0 for v in closure+divergence) and invariant == wedge == 0
    stress = (field*gi*field.T-g*invariant/4).applyfunc(sp.trigsimp)
    expected_stress = sp.zeros(4)
    expected_stress[3, 3] = mu
    curvature_scalar = sp.simplify(sp.trace(gi*ricci))
    residual = (ricci-g*curvature_scalar/2-stress).applyfunc(sp.simplify)
    assert stress == expected_stress and curvature_scalar == 0 and residual == sp.zeros(4)

    # Test transverse localization within this same source-free aligned ansatz.
    Ex, Ey = [sp.Function(name, real=True)(x, y, z) for name in ('Ex', 'Ey')]
    local_field = sp.zeros(4)
    local_field[3, 1], local_field[3, 2] = Ex, Ey
    local_field -= local_field.T
    local_raised = (gi*local_field*gi).applyfunc(sp.simplify)
    local_divergence = [sp.simplify(sum(sp.diff(local_raised[i, j], coords[i])
        for i in range(4))) for j in range(4)]
    local_closure = [sp.simplify(sp.diff(local_field[j, h], coords[i]) +
        sp.diff(local_field[h, i], coords[j])+sp.diff(local_field[i, j], coords[h]))
        for i in range(4) for j in range(i+1, 4) for h in range(j+1, 4)]
    div, curl = sp.diff(Ex, x)+sp.diff(Ey, y), sp.diff(Ey, x)-sp.diff(Ex, y)
    assert local_divergence == [div, 0, 0, 0] and local_closure == [0, 0, 0, -curl]
    laplacians = [sp.diff(v, x, 2)+sp.diff(v, y, 2) for v in (Ex, Ey)]
    assert sp.simplify(laplacians[0]-sp.diff(div, x)+sp.diff(curl, y)) == 0
    assert sp.simplify(laplacians[1]-sp.diff(div, y)-sp.diff(curl, x)) == 0
    intensity = Ex**2+Ey**2
    gradient_norm2 = sum(sp.diff(v, c)**2 for v in (Ex, Ey) for c in (x, y))
    intensity_laplacian = sp.diff(intensity, x, 2)+sp.diff(intensity, y, 2)
    assert sp.simplify(intensity_laplacian -
        2*(gradient_norm2+Ex*laplacians[0]+Ey*laplacians[1])) == 0
    local_invariant = sp.simplify(sum(local_field[i, j]*local_raised[i, j]
        for i in range(4) for j in range(4)))
    local_stress = (local_field*gi*local_field.T).applyfunc(sp.simplify)
    expected_local_stress = sp.zeros(4)
    expected_local_stress[3, 3] = intensity
    assert local_invariant == 0 and local_stress == expected_local_stress
    expected_local_residual = sp.zeros(4)
    expected_local_residual[3, 3] = mu-intensity
    assert (ricci-local_stress).applyfunc(sp.simplify) == expected_local_residual
    return dict(all_16_Einstein_Maxwell_residuals_zero=True,
        all_4_Maxwell_divergences_zero=True, all_4_dF_components_zero=True,
        potential_periodic=True, determinant=-1, electromagnetic_invariants_zero=True,
        coordinates='t,x,y,z; P=2e-a, Q=2e+a; mu=A^2',
        metric='g_tz=-1, g_xx=g_yy=1, g_xz=P*x, g_yz=Q*y, g_zz=(e-mu/2)*(x^2+y^2)-t',
        potential='A_x=A*sin(k*z)/k, A_y=-A*cos(k*z)/k; k positive integer, z period 2*pi',
        field_nonzero={'zx': 'A*cos(k*z)', 'zy': 'A*sin(k*z)'},
        raised_field_nonzero={'tx': '-A*cos(k*z)', 'ty': '-A*sin(k*z)'},
        stress_nonzero={'zz': 'A^2'}, Ricci_nonzero={'zz': str(ricci[3, 3])},
        convention='Signature -+++; c=8*pi*G=1; T_ab=F_ac*F_b^c-g_ab*F_cd*F^cd/4; antisymmetry implied',
        energy='For any unit timelike U, rho=A^2*(U^z)^2>=0; ordinary null Maxwell radiation',
        symmetry='Metric and stress are z independent; the electromagnetic field itself rotates with z',
        slicing='t is a time function only where D=(e-mu/2)*(x^2+y^2)-t-P^2*x^2-Q^2*y^2>0',
        scope='Alternative exact source for this support geometry; neither added to the scalars nor a stability proof',
        limitations=['Presupposes the Ori geometry and periodic spatial identification',
                    'Uniform transverse idealization; no finite apparatus or exterior matching supplied',
                    'Does not establish attainable initial data, quantum viability or safe traversal'],
        precedent='Balakin, van Holten and Kerner, CQG 17 (2000) 5009-5024, section 3.1; conventions converted locally',
        source='https://doi.org/10.1088/0264-9381/17/24/306',
        source_preprint='https://arxiv.org/abs/gr-qc/0009016',
        originality='Local verification of a specialization of a published solution family; no novelty claim',
        localization=dict(divergence_and_curl_form_verified=True,
            harmonic_component_identities_verified=True, intensity_laplacian_identity_verified=True,
            unchanged_geometry_Einstein_residual_verified=True,
            ansatz='Real F=dz wedge(Ex(x,y,z)*dx+Ey(x,y,z)*dy), no t dependence, same metric structure',
            Maxwell_divergence='(d_x Ex+d_y Ey,0,0,0) in coordinate order t,x,y,z',
            exterior_derivative='dF=(d_y Ex-d_x Ey)*dx wedge dy wedge dz; other independent components vanish',
            harmonicity='Source-free Maxwell equations imply Delta_perp Ex=Delta_perp Ey=0 at each fixed z',
            compact_support=dict(method='Analytic integration-by-parts deduction, not a numerical proof',
                assumptions=['Real C2 fields on the complete transverse R2 at each fixed z',
                             'Smooth compact transverse support with no singularities or surface sources',
                             'Source-free Maxwell equations and the stated pure aligned ansatz'],
                proof='Integral |grad Ei|^2 = -Integral Ei*Delta_perp Ei = 0; each Ei is constant, and compact support forces zero',
                conclusion='The only smooth compactly supported transverse field in this ansatz is Ex=Ey=0'),
            unchanged_geometry='Einstein equations require Ex^2+Ey^2=mu=A^2; harmonicity gives Delta_perp(intensity)=2*(|grad Ex|^2+|grad Ey|^2), so both fields must be transversely uniform',
            scope='Does not exclude general finite Maxwell pulses or beams with additional field components, sources, boundaries or a broader geometry'))


def shifted_null_orbit_checks():
    """Check a displaced null circle and the isometry behind that displacement.

    Constant cross terms shift the circle relative to a chosen finite patch;
    they are not a physical perturbation of the complete vacuum core.
    """
    t, x, y, z = coords = sp.symbols('t x y z', real=True)
    a = sp.symbols('a', real=True, nonzero=True)
    e, cx, cy = sp.symbols('e cx cy', real=True)
    P, Q, H = (2*e-a)*x+cx, (2*e+a)*y+cy, e*(x*x+y*y)-t
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, P],
                   [0, 0, 1, Q], [-1, P, Q, H]])
    _, connection, ricci = geometry(g, coords)
    assert ricci == sp.zeros(4) and g.det() == -1
    xn, yn = cx/a, -cy/a
    tn = e*(xn*xn+yn*yn)
    orbit = {x: xn, y: yn, t: tn}
    acceleration = sp.Matrix([connection[i][3][3] for i in range(4)])
    expected = sp.Matrix([(P*P+Q*Q-H)/2-e*(x*P+y*Q),
                          (-a*x+cx)/2, (a*y+cy)/2, -sp.Rational(1, 2)])
    assert (acceleration-expected).applyfunc(sp.simplify) == sp.zeros(4, 1)
    assert sp.simplify(H.subs(orbit)) == 0
    assert acceleration.subs(orbit).applyfunc(sp.simplify) == sp.Matrix([0, 0, 0, -sp.Rational(1, 2)])
    # This affine change preserves the periodic z identification exactly.
    T, X, Y, Z = new_coords = sp.symbols('T X Y Z', real=True)
    old_coords = sp.Matrix([T+2*e*xn*X+2*e*yn*Y+tn, X+xn, Y+yn, Z])
    jacobian = old_coords.jacobian(new_coords)
    transformed = jacobian.T*g.subs(dict(zip(coords, old_coords)))*jacobian
    original = g.subs({cx: 0, cy: 0}).subs(dict(zip(coords, new_coords)))
    assert (transformed-original).applyfunc(sp.simplify) == sp.zeros(4)
    return dict(all_16_vacuum_Ricci_components_zero=True, determinant=-1,
        null_orbit_and_nonaffine_geodesic_verified=True, affine_isometry_verified=True,
        assumptions='a != 0; real e,cx,cy; z period 2*pi; r0 > 0',
        metric='g_tz=-1, g_xx=g_yy=1, g_xz=(2e-a)*x+cx, g_yz=(2e+a)*y+cy, g_zz=e*(x^2+y^2)-t',
        orbit='xN=cx/a, yN=-cy/a, tN=e*(xN^2+yN^2)',
        geodesic_acceleration='Gamma^mu_zz=(0,0,0,-1/2) on the null circle; z is nonaffine',
        patch_condition='Circle lies in x^2+y^2 < r0^2 iff sqrt(cx^2+cy^2) < abs(a)*r0; equality is the boundary',
        isometry='X=x-xN, Y=y-yN, T=t-2e*xN*x-2e*yN*y+e*(xN^2+yN^2), Z=z removes cx,cy',
        scope='Exact algebra in the prescribed vacuum family. The bound retains this circle inside a fixed patch; it is not physical perturbation tolerance, generic stability, a chronology-horizon proof, or removal of all closed causal curves',
        provenance='Local symbolic deduction from the Ori vacuum metric; no novelty or construction claim')


def exponent(a, mu):
    if not (a >= 0 and mu >= 0 and a+mu < 1/8):
        raise ValueError('Require a>=0, mu>=0, a+mu<1/8 for this real-root sector.')
    # Retain the a^2 term even when the two O(a) contributions would cancel.
    C = 1-8*mu
    D = math.sqrt((1-8*(a+mu))*(1+8*(a-mu)))
    delta = 8*mu+32*a*a/(C+D)
    return -delta/(2*(1+math.sqrt(1-delta)))


def rk4_modes(cases, steps, end=12*math.log(10)):
    p = np.array([exponent(a, mu)+2j*n for a, mu, n in cases])
    u = np.ones(len(cases), dtype=complex)
    w = -p.copy()
    h = end/steps
    def rhs(u, w):
        return w, -p*w
    for _ in range(steps):
        u1, w1 = rhs(u, w)
        u2, w2 = rhs(u+h*u1/2, w+h*w1/2)
        u3, w3 = rhs(u+h*u2/2, w+h*w2/2)
        u4, w4 = rhs(u+h*u3, w+h*w3)
        u += h*(u1+2*u2+2*u3+u4)/6
        w += h*(w1+2*w2+2*w3+w4)/6
    exact_u = np.exp(-p*end)
    exact_w = -p*exact_u
    errors = np.maximum(abs(u-exact_u)/(1+abs(exact_u)), abs(w-exact_w)/(1+abs(exact_w)))
    assert np.all(np.isfinite(u)) and np.all(np.isfinite(w))
    return u, w, errors


def numerical_checks():
    assert math.isclose(exponent(1e-20, 0), -4e-40, rel_tol=1e-14)
    # Genuine independent RK4 integration, compared to the exact nonconstant branch.
    cases = [(a, mu, n) for a in (0, .02, .06, .10)
             for mu in (0, .001, .01, .02) for n in (1, 2, 4)]
    end = 12*math.log(10)
    coarse = rk4_modes(cases, 8192)
    fine = rk4_modes(cases, 16384)
    finer = rk4_modes(cases, 32768)
    errors = [float(np.max(run[2])) for run in (coarse, fine, finer)]
    ratios = [errors[i]/errors[i+1] for i in range(2)]
    assert errors[2] < 1e-7 and all(14 < r < 18 for r in ratios), (errors, ratios)
    T = -math.exp(-end)
    rows = []
    for idx, (a, mu, n) in enumerate(cases):
        s = exponent(a, mu)
        u, w = finer[0][idx], finer[1][idx]
        dT, dz = -w/T, 1j*n*u
        rho = (abs(dT+dz)**2+abs((T+1)*dT-dz)**2)/(2*(2+T))
        rho_initial = n*n+s*s/2
        rho_exact = (-T)**(2*s)*(n*n*(T+2)/(T*T) +
                    s*s*(T*T+2*T+2)/(2*T*T*(T+2)))
        assert abs(rho/rho_exact-1) < 1e-7
        assert s <= 1e-15
        rows.append(dict(a=a, mu=mu, spectator_wave_number=n, exponent_s=s,
            energy_power=2*s-2, energy_amplification=float(rho/rho_initial),
            scaled_exact_state_error=float(finer[2][idx])))
    # Analytic support-matter response: ds/dmu < 0 throughout this domain.
    derivatives = [-1/math.sqrt(1+8*a-8*mu)-1/math.sqrt(1-8*a-8*mu)
                   for a, mu, _ in cases]
    assert max(derivatives) < 0
    return dict(count=len(cases), nonzero_support_cases=sum(mu > 0 for _, mu, _ in cases),
        method='classical RK4; three step refinements; independent exact-mode reference',
        steps=[8192, 16384, 32768], max_exact_errors=errors, refinement_ratios=ratios,
        end_S=-1e-12, minimum_a_plus_mu_margin=float(min(.125-a-mu for a, mu, _ in cases)),
        analytic_obstruction='s<=0 and ds/dmu<0: added positive null support does not remove this spectator instability',
        cases=rows)


def twist_checks():
    """Test Ori's published spatial twist, including the exact scalar support."""
    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    a, mu, w, b, c, d, r = sp.symbols('a mu w b c d r', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, w*y],
                   [0, 0, 1, -w*x], [-1, w*y, -w*x, F-T+w*w*(x*x+y*y)]])
    gi, _, ricci = geometry(g, coords)
    expected = sp.zeros(4)
    expected[3, 3] = mu
    assert ricci == expected and g.det() == -1
    assert gi[3, 3] == 0 and gi[0, 3] == -1
    assert all(sp.simplify(sp.diff(gi[i, 3], coords[i])) == 0 for i in range(4))
    S = T+b*x*x+2*c*x*y+d*y*y
    norm = sp.expand(sum(gi[i, j]*sp.diff(S, coords[i])*sp.diff(S, coords[j])
                         for i in range(4) for j in range(4)))
    closure = [4*(b*b+c*c)-4*w*c-b-(a-mu)/2,
               4*(d*d+c*c)+4*w*c-d+(a+mu)/2,
               8*c*(b+d)+4*w*(b-d)-2*c]
    assert sp.expand(norm-S-closure[0]*x*x-closure[1]*y*y-closure[2]*x*y) == 0
    box_S = sp.simplify(sum(sp.diff(gi[i, j]*sp.diff(S, coords[j]), coords[i])
                            for i in range(4) for j in range(4)))
    assert box_S == 1+2*b+2*d
    assert sp.expand(sum(closure[:2])-(4*(b*b+d*d+2*c*c)-(b+d)+mu)) == 0
    # Explicit real branch: r=4*tr(B)-1; its squared value solves this polynomial.
    den = r*r+16*w*w
    h, v, cc = (r+1)/4, a*r/den, -2*w*a/den
    substitution = {b: (h+v)/2, d: (h-v)/2, c: cc}
    polynomial = (r*r-1+8*mu)*(r*r+16*w*w)+16*a*a
    assert sp.factor((closure[0]-closure[1]).subs(substitution)) == 0
    assert sp.factor(closure[2].subs(substitution)) == 0
    assert sp.factor(sum(closure[:2]).subs(substitution)-polynomial/(8*den)) == 0

    rows, max_closure_error, max_mode_error = [], 0., 0.
    end = 12*math.log(10)
    sample = np.linspace(0, end, 101)
    for aa in (0., .02, .06, .10):
        for mm in (0., .01, .02):
            for ww in (0., .1, .5, 1., 3., 10.):
                C, W = 1-8*mm, 16*ww*ww
                discriminant = (C+W)**2-64*aa*aa
                assert discriminant > 0
                # Rationalized to avoid subtracting two large terms for large twist.
                delta = 8*mm+32*aa*aa/(C+W+math.sqrt(discriminant))
                u = 1-delta
                assert 0 < u <= 1
                rr = -math.sqrt(u)
                hh = delta/(4*(1+math.sqrt(u)))
                vv, cc = aa*rr/(u+W), -2*ww*aa/(u+W)
                bb, dd = (hh+vv)/2, (hh-vv)/2
                residuals = [4*(bb*bb+cc*cc)-4*ww*cc-bb-(aa-mm)/2,
                             4*(dd*dd+cc*cc)+4*ww*cc-dd+(aa+mm)/2,
                             8*cc*(bb+dd)+4*ww*(bb-dd)-2*cc]
                max_closure_error = max(max_closure_error, max(map(abs, residuals)))
                assert abs(hh-4*(bb*bb+dd*dd+2*cc*cc)-mm) < 1e-13
                s = -2*hh
                for n in (1, 2, 4):
                    p = s+2j*n
                    sol = solve_ivp(lambda t, Y: (Y[1], -p*Y[1]), (0, end),
                        np.array([1, -p], dtype=complex), t_eval=sample,
                        method='DOP853', rtol=1e-11, atol=1e-13)
                    assert sol.success
                    ref = np.vstack([np.exp(-p*sample), -p*np.exp(-p*sample)])
                    error = float(np.max(abs(sol.y-ref)/(1+abs(ref))))
                    max_mode_error = max(max_mode_error, error)
                    rows.append(dict(a=aa, mu=mm, twist=ww, n=n, exponent_s=s,
                        energy_power=2*s-2, scaled_exact_state_error=error))
    assert max_closure_error < 1e-13 and max_mode_error < 1e-8
    return dict(all_16_Einstein_scalar_residuals_zero=True,
        support_and_periodicity='Same two periodic positive-kinetic scalars; Rzz=mu',
        exact_closure_verified=True, real_branch_domain='a>=0, mu>=0, a+mu<1/8, finite real twist',
        analytic_trace_identity='tr(B) = 4*tr(B^2) + mu >= 0; s=-2*tr(B)<=0',
        initial_time_function='(grad S)^2=S; S<0 admits no closed causal curves',
        conclusion='Spatial twist cannot remove this spectator instability in the real quadratic sector',
        geometry_cases=72, mode_evolutions=len(rows), nonzero_twist_modes=sum(row['twist'] != 0 for row in rows),
        max_closure_residual=max_closure_error, max_exact_mode_error=max_mode_error,
        end_S=-1e-12, cases=rows,
        source='Ori 2005, footnote 2; local extension and diagnostic are not claimed original')


def null_screen_jacobi_checks():
    """Verify a conditional null-orbit focusing test, not horizon formation.

    At a regular central null orbit, independent x/y reflections make the
    transverse coordinate frame parallel along xi=partial_z. In that frame
    the screen curvature K is constant. This is stronger than invariance of
    the axis alone; numerical applications must check the reflection residuals.
    """
    h, lam = sp.symbols('h lam', real=True)
    b = sp.symbols('b', real=True, nonzero=True)
    X, Y = sp.symbols('X Y', positive=True)
    ht, bt, lt, Xt, Yt, sx, sy, u, v = sp.symbols('ht bt lt Xt Yt sx sy u v', real=True)
    axis_g = sp.Matrix([[h, 0, 0, b], [0, X, 0, 0], [0, 0, Y, 0], [b, 0, 0, lam]])
    inverse = axis_g.inv()
    dg = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for i, j, value in [(0, 0, ht), (0, 3, bt), (1, 1, Xt), (2, 2, Yt), (3, 3, lt)]:
        dg[0, i, j] = dg[0, j, i] = value
    for derivative, i, j, value in [(1, 0, 1, sx), (1, 1, 3, u),
                                     (2, 0, 2, sy), (2, 2, 3, v)]:
        dg[derivative, i, j] = dg[derivative, j, i] = value

    def connection(a, c, d):
        return sp.simplify(sum(inverse[a, m]*(dg[c, m, d]+dg[d, m, c]-dg[m, c, d])/2
                               for m in range(4)))

    assert all(connection(a, 3, i) == 0 for a in range(4) for i in (1, 2))
    acceleration = [sp.simplify(connection(a, 3, 3).subs(lam, 0)) for a in range(4)]
    assert acceleration == [0, 0, 0, -lt/(2*b)]

    # Same root-coordinate exact family as ori_evolution.background(), allowing
    # symbolic a, e and mu. Differentiate before restricting to the null orbit.
    t, x, y, z = coords = sp.symbols('t x y z', real=True)
    a, mu, e = sp.symbols('a mu e', real=True)
    P, Q = 2*e-a, 2*e+a
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, P*x], [0, 0, 1, Q*y],
                   [-1, P*x, Q*y, (e-mu/2)*(x*x+y*y)-t]])
    _, gamma, ricci = geometry(g, coords)

    def screen_component(i, j):
        value = sum(g[i, m]*(sp.diff(gamma[m][3][3], coords[j])-
            sp.diff(gamma[m][3][j], z)+sum(gamma[m][j][n]*gamma[n][3][3]-
            gamma[m][3][n]*gamma[n][3][j] for n in range(4))) for m in range(4))
        return sp.simplify(value.subs({t: 0, x: 0, y: 0}))

    K = sp.Matrix(2, 2, lambda i, j: screen_component(i+1, j+1))
    nonaffinity = sp.simplify(gamma[3][3][3].subs({t: 0, x: 0, y: 0}))
    assert K == sp.diag((mu-a)/2, (mu+a)/2)
    assert nonaffinity == -sp.Rational(1, 2) and sp.trace(K) == ricci[3, 3] == mu

    Z, r, kappa_orbit = sp.symbols('Z r kappa_orbit', real=True)
    f = sp.Function('f')(Z)
    J = sp.exp(kappa_orbit*Z/2)*f
    transformed = sp.simplify((sp.diff(J, Z, 2)-kappa_orbit*sp.diff(J, Z)+r*J)/sp.exp(kappa_orbit*Z/2))
    assert sp.simplify(transformed-sp.diff(f, Z, 2)-(r-kappa_orbit**2/4)*f) == 0

    # Integrate the Jacobi initial-value problem independently of its reference
    # formula. Absence of an event on this finite interval is not a proof for
    # all windings; that conclusion comes from the exact sinh/linear formula.
    rows = []
    for mm in (sp.Integer(0), sp.Rational(1, 25), sp.Rational(21, 200), sp.Rational(3, 25)):
        eigenvalues = np.array([float(value) for value in K.subs({a: sp.Rational(1, 50), mu: mm}).diagonal()])
        kap = float(nonaffinity)
        gaps = eigenvalues-kap*kap/4
        expected_first = math.pi/math.sqrt(max(gaps)) if max(gaps) > 0 else None
        end = max(40., 1.1*expected_first) if expected_first is not None else 40.
        sample = np.linspace(0., end, 401)

        def ode(_, state):
            return np.r_[state[2:], kap*state[2:]-eigenvalues*state[:2]]

        def event_x(_, state):
            return state[0]

        def event_y(_, state):
            return state[1]

        solution = solve_ivp(ode, (0., end), (0., 0., 1., 1.), method='DOP853',
            t_eval=sample, events=(event_x, event_y), rtol=2e-11, atol=1e-13, max_step=.2)
        assert solution.success
        positions, velocities = [], []
        for gap in gaps:
            if gap > 0:
                frequency = math.sqrt(gap)
                amplitude, derivative = np.sin(frequency*sample)/frequency, np.cos(frequency*sample)
            elif gap < 0:
                rate = math.sqrt(-gap)
                amplitude, derivative = np.sinh(rate*sample)/rate, np.cosh(rate*sample)
            else:
                amplitude, derivative = sample, np.ones_like(sample)
            envelope = np.exp(kap*sample/2)
            positions.append(envelope*amplitude)
            velocities.append(envelope*(derivative+kap*amplitude/2))
        reference = np.vstack(positions+velocities)
        error = float(np.max(abs(solution.y-reference)/(1+abs(reference))))
        assert error < 2e-10
        events = sorted(float(value) for values in solution.t_events for value in values if value > 1e-7)
        event_error = None
        if expected_first is None:
            assert not events
        else:
            assert events
            event_error = abs(events[0]-expected_first)
            assert event_error < 1e-7
        rows.append(dict(a=.02, mu=float(mm), screen_eigenvalues=eigenvalues.tolist(),
            kappa_orbit=kap, eigenvalue_bound_satisfied=bool(max(gaps) <= 0),
            weaker_trace_bound_satisfied=bool(float(mm) <= kap*kap/2),
            analytic_conjugate_pairs_exist=expected_first is not None,
            analytic_first_conjugate_delta_z=expected_first,
            analytic_first_conjugate_windings=None if expected_first is None else expected_first/(2*math.pi),
            numerical_positive_event_times=events, numerical_end_z=end,
            scaled_exact_state_error=error, first_event_absolute_error=event_error))
    assert rows[2]['eigenvalue_bound_satisfied']  # Critical linear solution, not oscillation.
    assert rows[3]['weaker_trace_bound_satisfied'] and not rows[3]['eigenvalue_bound_satisfied']

    return dict(general_reflection_jet_transport_coefficients_zero=8,
        general_null_orbit_nonaffinity=str(-lt/(2*b)),
        exact_screen_curvature=[str(K[0, 0]), str(K[1, 1])], exact_kappa_orbit=str(nonaffinity),
        exact_screen_trace_equals_Rzz=True, jacobi_transform_verified=True,
        convention='Signature -+++; gravitational coupling 8*pi*G=1 and scalar winding number=1; kappa_orbit denotes nonaffinity, not gravitational coupling; R_izjz uses the existing code convention; xi=partial_z has period 2*pi',
        assumptions=['A regular C2 metric at an actual lambda=gzz=0 orbit, with nonzero xi and b=gtz',
                     'z independence and independent x and y reflection isometries throughout the orbit',
                     'Positive transverse metric; e_x=partial_x/sqrt(gxx), e_y=partial_y/sqrt(gyy) are parallel along xi',
                     'Numerical applications must check these symmetries; invariance of the axis alone is insufficient'],
        jacobi_equation='J_zz-kappa_orbit*J_z+K*J=0; J=exp(kappa_orbit*z/2)*Y gives Y_zz+(K-kappa_orbit^2*I/4)*Y=0',
        necessary_achronality_condition='largest_eigenvalue(K)<=kappa_orbit^2/4; exactly equivalent to absence of conjugate pairs along this constant-screen orbit, not to full achronality',
        weaker_trace_condition='At the null orbit tr(K)=Rzz=gravitational_coupling*wave_number^2*|Phi_H|^2<=kappa_orbit^2/2 is necessary but not sufficient, even for absence of conjugate pairs; coupling and wave_number both equal 1 here',
        affine_scope='ds/dz=C*exp(kappa_orbit*z), C>0; nonzero kappa_orbit gives finite affine length in one direction. Any finite conjugacy interval occurs before that endpoint, even after multiple windings',
        killing_horizon_scope='A full null hypersurface generated by xi has zero expansion and shear; Raychaudhuri requires Rzz=0 and Phi_H=0. An isolated closed null orbit does not meet that additional hypothesis',
        numerical_scope='Only the stated finite intervals are integrated. Analytic conjugacy classification follows from the constant-coefficient Jacobi solutions',
        example_scope='mu=.105 and .12 are exact local geometric benchmarks, not assertions that the original compact-onset construction or the current evolved initial data reach those cases',
        limits='Not applicable as a null-generator theorem at lambda>0; does not prove formation, full achronality, generic stability, curvature singularity, or impossibility of other chronology onsets',
        sources=['https://arxiv.org/pdf/gr-qc/0609119 (Theorem 2.30 and Lemma 2.35: conjugacy, causal deformation and screen Jacobi fields)',
                 'https://arxiv.org/pdf/gr-qc/0111067 (Eq.2.6: R_ll=0 on a non-expanding horizon)',
                 'https://arxiv.org/pdf/math/9909158 (Theorem 4.1: explicit completeness hypotheses for null splitting)'],
        cases=rows)


def twisted_null_screen_jacobi_checks():
    """Check removal of one conjugacy obstruction by Ori's existing rotation.

    Exact certificates supplement finite ODE checks. This does not establish
    full achronality, formation, or nonlinear stability; the spectator grows.
    """
    from scipy.linalg import expm

    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    a, mu, w = sp.symbols('a mu w', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, w*y], [0, 0, 1, -w*x],
                   [-1, w*y, -w*x, F-T+w*w*(x*x+y*y)]])
    inverse, connection, ricci = geometry(g, coords)
    axis = {T: 0, x: 0, y: 0}
    Omega = sp.Matrix(2, 2, lambda i, j: sp.simplify(connection[i+1][3][j+1].subs(axis)))

    def screen_component(i, j):
        value = sum(g[i, m]*(sp.diff(connection[m][3][3], coords[j])-
            sp.diff(connection[m][3][j], z)+sum(connection[m][j][n]*connection[n][3][3]-
            connection[m][3][n]*connection[n][3][j] for n in range(4))) for m in range(4))
        return sp.simplify(value.subs(axis))

    K = sp.Matrix(2, 2, lambda i, j: screen_component(i+1, j+1))
    kappa_orbit = sp.simplify(connection[3][3][3].subs(axis))
    assert Omega == sp.Matrix([[0, w], [-w, 0]])
    assert K == sp.diag((mu-a)/2, (mu+a)/2) and kappa_orbit == -sp.Rational(1, 2)
    assert ricci == sp.diag(0, 0, 0, mu) and g.det() == -1

    # Check every Omega term and the scalar envelope, including its sign.
    f = sp.Matrix([sp.Function('f1')(z), sp.Function('f2')(z)])
    j = sp.exp(kappa_orbit*z/2)*f
    original = j.diff(z, 2)+(2*Omega-kappa_orbit*sp.eye(2))*j.diff(z)+(
        K+Omega*Omega-kappa_orbit*Omega)*j
    transformed = f.diff(z, 2)+2*Omega*f.diff(z)+(K+Omega*Omega-kappa_orbit**2*sp.eye(2)/4)*f
    assert (original/sp.exp(kappa_orbit*z/2)-transformed).applyfunc(sp.simplify) == sp.zeros(2, 1)

    aa, mm = sp.Rational(1, 50), sp.Rational(3, 25)
    fixed_K = K.subs({a: aa, mu: mm})
    om = Omega.subs(w, sp.Rational(1, 10))
    frequency = sp.sqrt(3)/20
    sn, co = sp.sin(frequency*z)/frequency, sp.cos(frequency*z)
    Y = sp.Matrix([[sn/3+2*z*co/3, -z*sn/10], [z*sn/10, -sn+2*z*co]])
    normalized_potential = fixed_K+om*om-kappa_orbit**2*sp.eye(2)/4
    assert (Y.diff(z, 2)+2*om*Y.diff(z)+normalized_potential*Y).applyfunc(sp.simplify) == sp.zeros(2)
    assert Y.subs(z, 0) == sp.zeros(2) and Y.diff(z).subs(z, 0) == sp.eye(2)
    determinant = (4*z*z-sn*sn)/3
    assert sp.trigsimp(Y.det()-determinant) == 0

    # The pointwise identity below integrates to a positive index form for
    # endpoint-zero vectors whenever the real symmetric Riccati witness exists.
    l11, l12, l22, q1, q2, v1, v2 = sp.symbols('l11 l12 l22 q1 q2 v1 v2', real=True)
    L = sp.Matrix([[l11, l12], [l12, l22]])
    q, velocity = sp.Matrix([q1, q2]), sp.Matrix([v1, v2])
    dq = velocity+Omega*q
    riccati_K = kappa_orbit*L-L*L-(Omega*L-L*Omega)
    identity = (dq.dot(dq)-(q.T*riccati_K*q)[0]-(dq-L*q).dot(dq-L*q)-
                (2*(velocity.T*L*q)[0]-kappa_orbit*(q.T*L*q)[0]))
    assert sp.expand(identity) == 0
    assert sp.expand(sp.trace(L*L)-sp.trace(L)**2/2-(l11-l22)**2/2-2*l12*l12) == 0
    assert sp.trace(Omega*L-L*Omega) == 0

    # Validate these real branches directly. The old sufficient sweep domain
    # a+mu<1/8 does not include this benchmark; U=0 at w=.1 is admissible here.
    witnesses = []
    for ww in (sp.Rational(1, 10), sp.Rational(1, 5)):
        C, W = 1-8*mm, 16*ww*ww
        discriminant = (C+W)**2-64*aa*aa
        U = sp.simplify((C-W+sp.sqrt(discriminant))/2)
        assert discriminant >= 0 and U >= 0 and U+W > 0
        rr, hh = -sp.sqrt(U), (1-sp.sqrt(U))/4
        vv, cc = aa*rr/(U+W), -2*ww*aa/(U+W)
        B = sp.Matrix([[(hh+vv)/2, cc], [cc, (hh-vv)/2]])
        L = -2*B
        om = Omega.subs(w, ww)
        assert (om*L-L*om+L*L-kappa_orbit*L+fixed_K).applyfunc(sp.simplify) == sp.zeros(2)
        if ww == sp.Rational(1, 10):
            assert L == sp.Matrix([[-sp.Rational(1, 4), sp.Rational(1, 20)],
                                   [sp.Rational(1, 20), -sp.Rational(1, 4)]])
        S = T+(sp.Matrix([x, y]).T*B*sp.Matrix([x, y]))[0]
        inv = inverse.subs({a: aa, mu: mm, w: ww})
        gradient = sp.Matrix([sp.diff(S, coordinate) for coordinate in coords])
        assert sp.simplify((gradient.T*inv*gradient)[0]-S) == 0
        box_S = sp.simplify(sum(sp.diff(inv[i, j]*gradient[j], coords[i])
                                 for i in range(4) for j in range(4)))
        assert sp.simplify(box_S-1-2*sp.trace(B)) == 0
        assert (gradient.T*inv[:, 3])[0] == -1 and inv[3, 3] == 0
        assert sum(sp.diff(inv[i, 3], coords[i]) for i in range(4)) == 0
        n = sp.symbols('n', integer=True)
        exponent_p = -2*sp.trace(B)+2*sp.I*n
        # Full wave residual divided by exp(i*n*z)*(-S)^p, using the verified
        # norm, Box(S), mixed-gradient, inverse-zz and divergence identities.
        assert sp.simplify(exponent_p*(box_S+exponent_p-1-2*sp.I*n)/S) == 0
        assert float(-2*sp.trace(B)) < 0
        witnesses.append(dict(w=float(ww), branch_U_exact=str(U),
            branch_domain_checked=True, symmetric_Riccati_residual_exactly_zero=True,
            witness_L_exact=[[str(L[i, j]) for j in range(2)] for i in range(2)],
            witness_L_numeric=np.array(L, dtype=float).tolist(),
            norm_S_equals_S=True, spectator_wave_residual_zero=True,
            spectator_exponent_s=float(-2*sp.trace(B))))

    # Independent ODE and matrix-exponential checks of the normalized Jacobi
    # matrix. The analytic certificates, not an event scan, prove no conjugacy.
    rows = []
    for ww in (0., .1, .2):
        om = np.array([[0., ww], [-ww, 0.]])
        kval = np.array(fixed_K, dtype=float)
        kap = float(kappa_orbit)
        potential = kval+om@om-kap*kap*np.eye(2)/4
        generator = np.block([[np.zeros((2, 2)), np.eye(2)], [-potential, -2*om]])

        def ode(_, state):
            yy, dy = state[:4].reshape(2, 2), state[4:].reshape(2, 2)
            return np.r_[dy.ravel(), (-2*om@dy-potential@yy).ravel()]

        def event(_, state):
            return np.linalg.det(state[:4].reshape(2, 2))

        sample = np.linspace(0., 80., 801)
        solution = solve_ivp(ode, (0., 80.), np.r_[np.zeros(4), np.eye(2).ravel()],
            method='DOP853', rtol=1e-11, atol=1e-13, max_step=.1, t_eval=sample, events=event)
        assert solution.success
        reference = []
        for value in sample:
            ref = expm(generator*value)[:, 2:]
            reference.append(np.r_[ref[:2].ravel(), ref[2:].ravel()])
        reference = np.array(reference).T
        error = float(np.max(abs(solution.y-reference)/(1+abs(reference))))
        assert error < 1e-8
        events = [float(value) for value in solution.t_events[0] if value > 1e-6]
        if ww == 0:
            assert events and abs(events[0]-math.pi/math.sqrt(.07-.0625)) < 1e-8
        else:
            assert not events
        rows.append(dict(w=ww, numerical_end_z=80., scaled_expm_state_error=error,
            positive_determinant_events=events,
            analytic_no_conjugate_pairs=ww != 0,
            analytic_certificate='untwisted oscillatory mode' if ww == 0 else
                'exact positive determinant and symmetric Riccati witness' if ww == .1 else
                'exact symmetric Riccati witness'))

    return dict(a=.02, mu=.12, exact_Omega='[[0,w],[-w,0]]',
        exact_K=['(mu-a)/2', '(mu+a)/2'], exact_kappa_orbit='-1/2',
        geometric_identities_verified=True, normalized_equation_verified=True,
        full_equation='j_zz+(2*Omega-kappa_orbit*I)*j_z+(K+Omega^2-kappa_orbit*Omega)*j=0',
        normalization='j=exp(kappa_orbit*z/2)*Y; Y(0)=0, Y_z(0)=I; Y_zz+2*Omega*Y_z+(K+Omega^2-kappa_orbit^2*I/4)*Y=0',
        exact_w_point1=dict(frequency='sqrt(3)/20', fundamental_matrix_equation_and_initial_conditions_verified=True,
            determinant='det(Y)=(4*z^2-(sin(omega*z)/omega)^2)/3',
            analytic_bound='|sin(omega*z)|<=|omega*z| gives det(Y)>=z^2; det(j)=exp(kappa_orbit*z)*det(Y)>0 for every finite z!=0'),
        index_form=dict(pointwise_identity_verified=True,
            witness_equation='[Omega,L]+L^2-kappa_orbit*L+K=0, L real symmetric',
            factorization='For q=0 at both endpoints: integral exp(-kappa_orbit*z)*(|D_z q|^2-q^T K q) = integral exp(-kappa_orbit*z)*|D_z q-Lq|^2, D_z=partial_z+Omega',
            implication='A Dirichlet Jacobi field makes the left side zero; the right side forces D_z q=Lq, whose zero initial value implies q=0. Thus no conjugate pairs occur on any finite interval'),
        trace_obstruction=dict(commutator_trace_and_square_inequality_verified=True,
            equation='B_z+[Omega,B]+B^2-kappa_orbit*B+K=0; theta_z=kappa_orbit*theta-tr(B^2)-mu',
            inequality='theta_z<=-(theta-kappa_orbit)^2/2-(mu-kappa_orbit^2/2)',
            implication='mu>kappa_orbit^2/2 forces conjugacy regardless of screen rotation; from a point, first delta_z<=pi*sqrt(2/(mu-kappa_orbit^2/2))',
            qualification='Necessary trace bound only; equality and smaller mu require separate analysis. The Jacobi bundle is optically twist-free; Omega is screen-frame transport, not optical vorticity'),
        branch_witnesses=witnesses, numerical_cases=rows,
        scope='Exact central null-support family, positive gravitational coupling=1 and winding number=1. These branches were validated outside the old sufficient spectator sweep sector',
        limits='No-conjugacy does not prove full achronality, chronology formation, generic or nonlinear stability. Finite event scans alone do not prove no-conjugacy. Negative spectator exponents remain; no claim about the live evolved initial data',
        sources=['https://arxiv.org/pdf/gr-qc/0503077 (footnote 2: the metric extension)',
                 'https://arxiv.org/pdf/gr-qc/0609119 (Theorem 2.30 and Lemma 2.35: conjugacy and screen Jacobi fields)',
                 'https://arxiv.org/pdf/math/9909158 (Eqs.2.2-2.4: matrix Riccati and Raychaudhuri equations)'])


def massive_spectator_checks():
    """Finite scalar mass leaves the singular branch on this fixed geometry.

    This is a separated linear spectator test, not a massive support solution
    or a test of the nonminimal field's improved stress-energy tensor.
    """
    import mpmath as mp

    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    a, mu, w, b, c, d, n, mass2, xi = sp.symbols('a mu w b c d n mass2 xi', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, w*y], [0, 0, 1, -w*x],
                   [-1, w*y, -w*x, F-T+w*w*(x*x+y*y)]])
    inverse, _, ricci = geometry(g, coords)
    scalar_R = sp.simplify(sp.trace(inverse*ricci))
    assert scalar_R == 0 and ricci == sp.diag(0, 0, 0, mu) and g.det() == -1
    S = T+b*x*x+2*c*x*y+d*y*y
    closure = [4*(b*b+c*c)-4*w*c-b-(a-mu)/2,
               4*(d*d+c*c)+4*w*c-d+(a+mu)/2,
               8*c*(b+d)+4*w*(b-d)-2*c]
    q, r = sp.symbols('q r')
    f = sp.Function('f')
    field = sp.exp(sp.I*n*z)*f(S)
    residual = sp.simplify(sum(sp.diff(inverse[i, j]*sp.diff(field, coords[j]), coords[i])
        for i in range(4) for j in range(4))/sp.exp(sp.I*n*z)-(mass2+xi*scalar_R)*f(S))
    first = sp.Subs(sp.diff(f(q), q), q, S)
    second = sp.Subs(sp.diff(f(q), q, 2), q, S)
    p = -2*(b+d)+2*sp.I*n
    expected = (S+closure[0]*x*x+closure[1]*y*y+closure[2]*x*y)*second+(1-p)*first-mass2*f(S)
    assert sp.simplify(residual-expected) == 0
    # Multiplying the x=-S equation by x leaves mass2*x below indicial order.
    X = sp.symbols('X', positive=True)
    trial = X**r
    indicial = sp.simplify((X*sp.diff(trial, X, 2)+(1-p)*sp.diff(trial, X)+mass2*trial)*X/trial)
    assert sp.simplify(indicial-r*(r-p)-mass2*X) == 0

    end = 12*math.log(10)
    sample = np.linspace(0., end, 61)

    def compare(exponent_p, mass, reference):
        exact = np.array([reference(t) for t in sample], dtype=complex).T
        solution = solve_ivp(lambda t, v: (v[1], -exponent_p*v[1]-mass*mass*math.exp(-t)*v[0]),
            (0., end), exact[:, 0], t_eval=sample, method='DOP853', rtol=1e-11, atol=1e-13)
        assert solution.success and np.all(np.isfinite(solution.y))
        error = float(np.max(abs(solution.y-exact)/(1+abs(exact))))
        assert error < 2e-8
        return dict(max_scaled_state_error=error,
            end_log_amplitude_slope=float(np.real(-solution.y[1, -1]/solution.y[0, -1])),
            end_amplitude=float(abs(solution.y[0, -1])),
            end_log_time_derivative=float(abs(solution.y[1, -1])))

    rows, resonances = [], []
    modes = [('untwisted_a_0.02_mu_0.04', complex(exponent(.02, .04), 2)),
             ('rotated_a_0.02_mu_0.12_w_0.1', complex(-.5, 2)),
             ('flat_nonzero_winding', 2j), ('rotated_zero_winding', complex(-.5, 0))]
    with mp.workdps(60):
        for label, exponent_p in modes:
            pp = mp.mpc(exponent_p)
            for mass in (0., 1., 10.):
                mm = mp.mpf(mass)

                def reference(t):
                    xx = mp.exp(-mp.mpf(t))
                    hyper = mp.hyp0f1(1+pp, -mm*mm*xx)
                    derivative = -mm*mm/(1+pp)*mp.hyp0f1(2+pp, -mm*mm*xx)
                    return xx**pp*hyper, -xx**pp*(pp*hyper+xx*derivative)

                check = compare(exponent_p, mass, reference)
                assert abs(check['end_log_amplitude_slope']-exponent_p.real) < 1e-8
                rows.append(dict(geometry=label, p_real=exponent_p.real, p_imag=exponent_p.imag,
                    dimensionless_mass=mass, predicted_log_amplitude_slope=exponent_p.real, **check))
        for exponent_p, order in ((0., 0), (-1., 1)):

            def bessel_reference(t):
                xx = mp.exp(-mp.mpf(t))
                return xx**(mp.mpf(exponent_p)/2)*mp.bessely(order, 2*mp.sqrt(xx))

            check = compare(exponent_p, 1., lambda t: (bessel_reference(t), mp.diff(bessel_reference, t)))
            resonances.append(dict(p=exponent_p, dimensionless_mass=1.,
                reference='x^(p/2)*Y_(-p)(2*sqrt(x))', **check))

    return dict(full_4D_reduction_verified=True, exact_R_scalar=0,
        equation='(Box-m^2-xi*R)chi=0; chi=exp(i*n*z)*F(S), S=T+x^T*B*x, z has period 2*pi',
        hypotheses='Fixed exact null-supported Ori geometry, real symmetric B satisfying the quadratic closure, mu>=0, integer n, finite constant m>=0 and xi; S<0',
        reduced_equation='S*F_SS+(1+2*tr(B)-2*i*n)*F_S-m^2*F=0',
        x_equation='x=-S>0; p=-2*tr(B)+2*i*n; x*G_xx+(1-p)*G_x+m^2*G=0',
        log_time_equation='tau=-log(x); u_tautau+p*u_tau+m^2*exp(-tau)*u=0',
        indicial_polynomial='r*(r-p); roots 0 and p are independent of finite mass',
        noninteger_basis=['0F1(;1-p;-m^2*x)', 'x^p*0F1(;1+p;-m^2*x)'],
        exceptions=dict(repeated_root='p=0: second solution log(x) for m=0, Y0(2*m*sqrt(x)) for m>0; logarithmic divergence remains',
            negative_integer='p=-N<0: m>0 basis x^(-N/2)*{J_N,Y_N}(2*m*sqrt(x)) retains an x^(-N) singularity and logarithmic terms; m=0 basis is {1,x^(-N)}',
            real_sector='tr(B)=4*tr(B^2)+mu and mu>=0 give 0<=tr(B)<=1/2, so the only negative-integer p in this real sector is -1',
            bounded_amplitude='tr(B)=0 and n!=0 has bounded oscillatory amplitude but a derivative growing as 1/x'),
        conclusion='Finite positive mass does not remove the singular branch; xi*R is identically zero. The regular branch requires setting the independent singular coefficient to zero',
        scope='Generic means nonzero singular coefficient within this separated two-solution sector only. No generic nonlinear instability theorem, matter production recipe, or human-traversal claim',
        support_limit='Giving the supporting scalars mass changes their field equations and stress tensor; this calculation does not supply a new massive support solution',
        energy_limit='No claim about the improved nonminimal stress-energy tensor or cancellations in its observer energy',
        numerical_method='DOP853 versus independent 60-digit 0F1/Bessel references; rtol=1e-11, atol=1e-13, 61 samples',
        numerical_end_S=-1e-12, dimensionless_mass_values=[0., 1., 10.],
        cases=rows, resonance_cases=resonances,
        max_scaled_state_error=max(row['max_scaled_state_error'] for row in rows+resonances),
        sources=['https://arxiv.org/pdf/1401.2026 (Eq.1 and Lagrangian 113: scalar mass and curvature coupling)',
                 'https://arxiv.org/pdf/gr-qc/0503077 (geometry and stability discussion; this massive reduction is our calculation)'])


def negative_null_support_checks():
    """A conditional branch bound with formal negative null stress, not a medium."""
    import mpmath as mp

    a = sp.Rational(1, 50)
    w, mu, b, c, d = sp.symbols('w mu b c d', real=True)
    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, w*y], [0, 0, 1, -w*x],
                   [-1, w*y, -w*x, F-T+w*w*(x*x+y*y)]])
    inverse, connection, ricci = geometry(g, coords)
    assert ricci == sp.diag(0, 0, 0, mu) and sp.simplify(sp.trace(inverse*ricci)) == 0
    closure = [4*(b*b+c*c)-4*w*c-b-(a-mu)/2,
               4*(d*d+c*c)+4*w*c-d+(a+mu)/2,
               8*c*(b+d)+4*w*(b-d)-2*c]
    h = b+d
    assert sp.expand(sum(closure[:2])-(4*(b*b+d*d+2*c*c)-h+mu)) == 0
    assert sp.expand(b*b+d*d+2*c*c-h*h/2-(b-d)**2/2-2*c*c) == 0
    q = sp.symbols('q', nonnegative=True)
    assert sp.expand((-sp.Rational(1, 2)-q)-2*(-sp.Rational(1, 2)-q)**2+1+3*q+2*q*q) == 0

    # r=4*tr(B)-1=-3 gives s=1; check the full quadratic closure.
    critical_mu = -1-2*a*a/(9+16*w*w)
    critical_v, critical_c = -3*a/(9+16*w*w), -2*w*a/(9+16*w*w)
    critical = {b: (-sp.Rational(1, 2)+critical_v)/2,
                d: (-sp.Rational(1, 2)-critical_v)/2, c: critical_c, mu: critical_mu}
    assert all(sp.simplify(value.subs(critical, simultaneous=True)) == 0 for value in closure)
    assert sp.simplify((-2*h).subs(critical, simultaneous=True)) == 1
    assert sp.simplify((9-(1-8*mu))*(9+16*w*w)+16*a*a-
                       8*(9+16*w*w)*(mu-critical_mu)) == 0
    mu0 = critical_mu.subs(w, 0)
    plus_root, minus_root = sp.sqrt(1-8*mu0+8*a), sp.sqrt(1-8*mu0-8*a)
    assert sp.simplify(plus_root+minus_root) == 6
    branches = []
    for sign_b in (-1, 1):
        for sign_d in (-1, 1):
            bb, dd = (1+sign_b*plus_root)/8, (1+sign_d*minus_root)/8
            assert all(sp.simplify(value.subs({b: bb, d: dd, c: 0, w: 0, mu: mu0},
                                             simultaneous=True)) == 0 for value in closure)
            s = sp.simplify(-2*(bb+dd))
            branches.append(dict(signs=[sign_b, sign_d], b_exact=str(bb), d_exact=str(dd),
                s_exact=str(s), s=float(s), kinetic_energy_power=float(2*s-2)))
    # These radii expose the bad branch's pre-existing singular initial data
    # on the complete selected slice, not a global instability counterexample.
    radius_x = sp.sqrt(1/(plus_root/4))
    radius_y = sp.sqrt(1/(minus_root/4))
    assert sp.simplify(-1+(plus_root/4)*radius_x**2) == 0
    assert sp.simplify(-1+(minus_root/4)*radius_y**2) == 0

    # Audit a smooth cutoff on the common -- slice using exact null rays.
    # Losing its domain of dependence is not a claim that growth disappears.
    bm, dm = (1-plus_root)/8, (1-minus_root)/8
    bp, dp = (1+plus_root)/8, (1+minus_root)/8
    selected = {w: 0, mu: mu0}
    metric0, inverse0 = g.subs(selected), inverse.subs(selected)
    Sm, Sp = T+bm*x*x+dm*y*y, T+bp*x*x+dp*y*y
    assert metric0.det() == -1
    for time_function in (Sm, Sp):
        gradient = sp.Matrix([sp.diff(time_function, coordinate) for coordinate in coords])
        assert sp.simplify((gradient.T*inverse0*gradient)[0]-time_function) == 0
    initial_Sp = sp.expand(Sp.subs(T, -1-bm*x*x-dm*y*y))
    assert initial_Sp == -1+sp.Rational(227, 300)*x*x+sp.Rational(223, 300)*y*y
    tau = sp.symbols('tau', real=True)
    old_coordinates = sp.Matrix([tau-bm*x*x-dm*y*y, x, y, z])
    jacobian = old_coordinates.jacobian((tau, x, y, z))
    pulled_back = sp.simplify(jacobian.T*metric0.subs(T, old_coordinates[0])*jacobian)
    expected_metric = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, 2*bm*x],
        [0, 0, 1, 2*dm*y], [-1, 2*bm*x, 2*dm*y, 4*bm*bm*x*x+4*dm*dm*y*y-tau]])
    assert pulled_back == expected_metric and sp.simplify(pulled_back.inv()[0, 0]-tau) == 0
    ray_a, ray_b, ray_z = sp.symbols('ray_a ray_b ray_z', real=True)
    v = sp.symbols('v', positive=True)  # v=-A/B; the inward family has 0<v<1.
    arrivals, ray_checks = [], []
    for axis, index, beta_minus, beta_plus in (('x', 1, bm, bp), ('y', 2, dm, dp)):
        ray_c, ray_d = -2*beta_minus, 2*beta_plus
        ray_h = ray_c+ray_d
        assert ray_d-ray_c == sp.Rational(1, 2) and ray_c > 0
        radial = ray_a*sp.exp(ray_c*ray_z)+ray_b*sp.exp(-ray_d*ray_z)
        ray_T = ray_c*ray_d*radial**2-sp.diff(radial, ray_z)**2
        path = sp.Matrix([ray_T, radial if index == 1 else 0, radial if index == 2 else 0, ray_z])
        path_sub = dict(zip(coords, path))
        velocity = path.diff(ray_z)
        ray_metric = metric0.subs(path_sub, simultaneous=True)
        assert sp.simplify((velocity.T*ray_metric*velocity)[0]) == 0
        assert sp.simplify((ray_metric*velocity)[3]) == 0
        for component in range(4):
            acceleration = path[component].diff(ray_z, 2)+velocity[component]/2
            acceleration += sum(connection[component][j][k].subs(selected).subs(path_sub, simultaneous=True)
                                *velocity[j]*velocity[k] for j in range(4) for k in range(4))
            assert sp.simplify(sp.expand(acceleration)) == 0
        ray_tau = Sm.subs(path_sub, simultaneous=True)
        expected_tau = -ray_h*ray_b**2*sp.exp(-2*ray_d*ray_z)/2+2*ray_c*ray_h*ray_a*ray_b*sp.exp(-ray_z/2)
        assert sp.simplify(ray_tau-expected_tau) == 0
        normalized_B2 = 2/(ray_h*(1+4*ray_c*v))
        assert sp.simplify(expected_tau.subs({ray_z: 0, ray_a: -v*ray_b})
                           .subs(ray_b**2, normalized_B2)+1) == 0
        arrival_z = -sp.log(v)/ray_h
        assert sp.simplify(radial.subs({ray_a: -v*ray_b, ray_z: arrival_z})) == 0
        assert sp.simplify(ray_T.subs({ray_a: -v*ray_b, ray_z: arrival_z})
                           +ray_h**2*ray_b**2*v**(2*ray_d/ray_h)) == 0
        ray_checks.append(dict(axis=axis, c_exact=str(ray_c), d_exact=str(ray_d), H_exact=str(ray_h),
            null_norm_zero=True, Killing_momentum_zero=True, four_geodesic_residuals_zero=True,
            initial_slice_and_central_arrival_exact=True))
        cc, dd, hh = float(ray_c), float(ray_d), float(ray_h)
        for radius in (1., 1.1, 1.14):
            constant = 1-hh*radius*radius/2
            linear = 2+2*cc*hh*radius*radius
            u = -2*constant/(linear+math.sqrt(linear*linear-4*constant))
            B = math.sqrt(2/(hh*(1-4*cc*u)))
            zf = math.log(-1/u)/hh
            Tf = -2*hh*(-u)**(2*dd/hh)/(1-4*cc*u)
            assert -1 < u < 0 and Tf < 0
            assert abs(B*(1+u)-radius) < 1e-14
            endpoint_error = abs(B*(u*math.exp(cc*zf)+math.exp(-dd*zf)))
            assert endpoint_error < 1e-13
            arrivals.append(dict(axis=axis, initial_radius=radius, central_T=Tf,
                unwrapped_z_travel=zf, windings=zf/(2*math.pi), endpoint_error=endpoint_error))

    def compare(exponent_p, decades):
        sample = np.linspace(0., decades*math.log(10), 10*decades+1)
        solution = solve_ivp(lambda t, v: (v[1], -exponent_p*v[1]), (0., sample[-1]),
            np.array([1, -exponent_p]), t_eval=sample, method='DOP853', rtol=1e-11, atol=1e-24)
        exact = np.vstack([np.exp(-exponent_p*sample), -exponent_p*np.exp(-exponent_p*sample)])
        assert solution.success and np.all(np.isfinite(solution.y))
        relative = float(np.max(abs(solution.y-exact)/abs(exact)))
        scaled = float(np.max(abs(solution.y-exact)/(1+abs(exact))))
        defect = float(np.max(abs(solution.y[0]+solution.y[1]/exponent_p)))
        assert scaled < 1e-8 and defect < 1e-12
        return dict(max_relative_state_error=relative, max_scaled_state_error=scaled,
                    max_constant_mode_invariant_defect=defect)

    rows = []
    for ww in (0., .1, 1.):
        threshold = float(critical_mu.subs(w, ww))
        for offset in (.01, 0., -.01):
            mm = threshold+offset
            C, W = 1-8*mm, 16*ww*ww
            discriminant = (C+W)**2-64*float(a)**2
            assert discriminant > 0
            Uplus = (C-W+math.sqrt(discriminant))/2
            assert Uplus > 0
            s = (math.sqrt(Uplus)-1)/2
            check = compare(complex(s, 2), 6)
            assert check['max_relative_state_error'] < 1e-8
            assert (s < 1 if offset > 0 else abs(s-1) < 1e-14 if offset == 0 else s > 1)
            rows.append(dict(w=ww, mu=mm, s=s, kinetic_energy_power=2*s-2, **check))
    precision_probe = compare(complex(1, 2), 8)

    # A massive zero-winding resonant solution at s=1 has a logarithmic derivative.
    mass_exception = []
    with mp.workdps(60):
        for xx in (mp.mpf('1e-4'), mp.mpf('1e-8'), mp.mpf('1e-12')):
            field = lambda value: mp.sqrt(value)*mp.bessely(1, 2*mp.sqrt(value))
            residual = xx*mp.diff(field, xx, 2)+field(xx)
            assert abs(residual) < mp.mpf('1e-45')
            mass_exception.append(dict(x=float(xx), field=float(field(xx)),
                derivative=float(mp.diff(field, xx)), equation_residual=float(abs(residual))))

    # Normalize the null direction and absorb its scale into A>0. The k+/-
    # sum proves the compensation obstruction for every orientation of ell.
    eta = sp.diag(-1, 1, 1, 1)
    A, Ccas = sp.symbols('A Ccas', positive=True)
    nx, ny, nz = sp.symbols('nx ny nz', real=True)
    ell_cov = sp.Matrix([-1, nx, ny, nz])
    target, casimir = -A*ell_cov*ell_cov.T, Ccas*sp.diag(-1, 1, 1, -3)
    assert target.rank() == 1 and casimir.rank() == 4
    assert sp.expand(sp.trace(eta*target*eta*target)-A*A*(-1+nx*nx+ny*ny+nz*nz)**2) == 0
    assert sp.trace(eta*casimir*eta*casimir) == 12*Ccas*Ccas
    target_sum = 0
    for sign in (-1, 1):
        ray = sp.Matrix([1, sign, 0, 0])
        assert (ray.T*eta*ray)[0] == 0 and (ray.T*casimir*ray)[0] == 0
        target_sum += (ray.T*target*ray)[0]
    assert sp.expand(target_sum+2*A*(1+nx*nx)) == 0

    # Additivity of independently supplied ideal tensors is an assumption;
    # this does not compute the vacuum state of intersecting physical plates.
    orthogonal_sources = []
    for axis in range(3):
        plate_tensor = Ccas*sp.diag(-1, 1, 1, 1)
        plate_tensor[axis+1, axis+1] = -3*Ccas
        orthogonal_sources.append(plate_tensor)
    multiple_casimir = sum(orthogonal_sources, sp.zeros(4))
    assert multiple_casimir == Ccas*sp.diag(-3, -1, -1, -1)
    forward_ray = sp.Matrix([1, 0, 0, 1])
    forward_tensor = (eta*forward_ray)*(eta*forward_ray).T
    fixed_target = -A*forward_tensor
    five_rays = [sp.Matrix(ray) for ray in ((1, 1, 0, 0), (1, -1, 0, 0),
                 (1, 0, 1, 0), (1, 0, -1, 0), (1, 0, 0, -1))]
    five_beams = sum((A*(eta*ray)*(eta*ray).T for ray in five_rays), sp.zeros(4))
    assert multiple_casimir.subs(Ccas, 2*A)+five_beams == fixed_target
    cone_speed = 2*sp.sqrt(2)/3
    cone_rays = [sp.Matrix(ray) for ray in ((1, cone_speed, 0, -sp.Rational(1, 3)),
                 (1, -cone_speed, 0, -sp.Rational(1, 3)),
                 (1, 0, cone_speed, -sp.Rational(1, 3)),
                 (1, 0, -cone_speed, -sp.Rational(1, 3)))]
    four_beams = sum((3*A*(eta*ray)*(eta*ray).T/4 for ray in cone_rays), sp.zeros(4))
    assert sp.simplify(multiple_casimir.subs(Ccas, 4*A/3)+four_beams-fixed_target) == sp.zeros(4)
    required_radiation = fixed_target-multiple_casimir
    assert sp.simplify(required_radiation.extract([0, 3], [0, 3]).det()-Ccas*(3*Ccas-4*A)) == 0
    # The cone construction plus a forward beam proves sufficiency throughout
    # Ccas>=4*A/3. This minimizes Ccas in this ansatz, not the number of beams.
    general_cone = four_beams.subs(A, 3*Ccas/4)
    assert sp.simplify(multiple_casimir+general_cone+(3*Ccas/4-A)*forward_tensor-fixed_target) == sp.zeros(4)
    observer = sp.Matrix(sp.symbols('u0 u1 u2 u3', real=True))
    weight = sp.symbols('weight', positive=True)
    for ray in five_rays+cone_rays:
        assert (ray.T*eta*ray)[0] == 0 and ray[0] == 1
        beam = weight*(eta*ray)*(eta*ray).T
        contraction = (eta*ray).dot(observer)
        assert sp.expand((observer.T*beam*observer)[0]-weight*contraction**2) == 0
        assert sp.simplify(-eta*beam*observer+weight*contraction*ray) == sp.zeros(4, 1)

    return dict(a=.02, exact_R_scalar=0, exact_Rzz='mu',
        hypotheses='Fixed Ori geometry, a=.02, real symmetric quadratic B, finite real w, formal mu<0; positive-kinetic massless spectator unless explicitly discussing the massive exception',
        spectator='chi=exp(i*n*z)*(-S)^p, S=T+x^T*B*x<0, p=s+2*i*n, integer n, s=-2*tr(B)',
        bounded_kinetic_criterion='For nonzero p and nonzero branch coefficient, regular-observer kinetic energy scales as (-S)^(2*s-2), so it is bounded iff s>=1',
        untwisted=dict(reality='mu<=1/8-|a|; all mu<0 allowed for a=.02',
            favorable_exponent='s=(sqrt(1-8*mu+8*a)+sqrt(1-8*mu-8*a)-2)/4',
            critical_mu_exact=str(mu0), critical_mu=float(mu0), all_four_branches=branches),
        rotating=dict(favorable_threshold='mu<=-1-2*a^2/(9+16*w^2)', critical_closure_exactly_zero=True,
            branches='C=1-8*mu, W=16*w^2; U=(C-W +/- sqrt((C+W)^2-64*a^2))/2; retain U>=0 and r=+/-sqrt(U), r^2+W>0; s=-(r+1)/2',
            selection='For a=.02 and mu<0, the larger U is positive with positive discriminant; the favorable branch uses r=-sqrt(Uplus). The r=+sqrt(Uplus) counterpart has s_other=-1-s',
            threshold_certificate='Q(9)=8*(9+W)*(mu-mu_critical). Since C>1 and a=.02, Uminus<=C-Uplus=16*a^2/(Uplus+W)<32*a^2/(C+W)<0.0128<9; thus Q(9)<=0 iff Uplus>=9 iff s>=1',
            trace_bound='mu=h-4*tr(B^2)<=h-2*h^2; h=-1/2-q, q>=0 gives mu<=-1-3*q-2*q^2<=-1',
            rotation_limit='For a!=0 no finite w reaches mu=-1 at s=1; the threshold approaches -1 as |w| grows'),
        initial_slice_caveat=dict(slice='At the untwisted critical parameter on the complete S_--=-1 slice, S_++=-1+(sqrt(1-8*mu+8*a)/4)*x^2+(sqrt(1-8*mu-8*a)/4)*y^2',
            S_plus_zero_axis_radii=[float(radius_x), float(radius_y)],
            implication='The ++ mode is already singular on this complete initial slice; on a small patch its singular surface may enter through the boundary. Other formal B branches are not a proof of global instability from common smooth initial data'),
        cutoff_causal_audit=dict(mu_exact=str(mu0), b_minus_exact=str(bm), d_minus_exact=str(dm),
            b_plus_exact=str(bp), d_plus_exact=str(dp),
            time_function_identities='grad(S_minus)^2=S_minus and grad(S_plus)^2=S_plus, checked with the exact inverse metric',
            coordinate_change='tau=S_minus; ds^2=(dx+2*b_minus*x*dz)^2+(dy+2*d_minus*y*dz)^2-2*dtau*dz-tau*dz^2; g^tautau=tau<0 on tau<0',
            initial_intersection='S_plus on tau=-1 is -1+(227/300)*x^2+(223/300)*y^2; its zero ellipse contains singular ++ data',
            ray_equations='On either axis, c=-2*beta_minus, d=2*beta_plus=c+1/2, H=c+d; r=A*exp(c*Z)+B*exp(-d*Z), T=c*d*r^2-(dr/dZ)^2, Z=z-z0',
            parameterization='A<0, B>0, u=A/B in (-1,0); B^2=2/[H*(1-4*c*u)] fixes tau(0)=-1; R^2=2*(1+u)^2/[H*(1-4*c*u)]',
            nonaffinity='q_second+Gamma(q_first,q_first)=-(1/2)*q_first; p_z=0. Z is unwrapped even though z is periodic',
            time_along_ray='tau=-H*B^2*exp(-2*d*Z)/2+2*c*H*A*B*exp(-Z/2)<0 and increases for A<0,B>0',
            central_arrival='Zf=log(-1/u)/H; Tf=tau_f=-2*H*(-u)^(2*d/H)/(1-4*c*u)<0',
            limit='As u approaches 0 from below, R^2 approaches 2/H=1/(beta_plus-beta_minus) and Tf approaches 0 from below',
            exact_axis_checks=ray_checks, arrival_table=arrivals,
            table_scope='Explicit axis null arrivals, not a minimization over all causal paths or a numerical PDE evolution',
            domain_statement='For fixed-background normally hyperbolic scalar evolution, matching both the field and its normal derivative on an initial set U guarantees equality only on their common domain of dependence',
            cutoff_conclusion='A fixed smooth transverse cutoff that changes the ++ data before its compact singular ellipse has mismatching points strictly inside that ellipse. The explicit rays put the center outside the guaranteed matching domain at a negative time; later central events retain such points in their causal past. Moving the cutoff toward the ellipse delays this loss toward zero, but no single such cutoff imports the exact ++ asymptotic arbitrarily near N by finite propagation alone',
            assumptions='Fixed untwisted critical negative-null-stress geometry; smooth data on tau=-1; cutoff matching the exact field and normal derivative near the central circle and regularizing the pole before the ellipse. For a z-independent cutoff, translation of z0 covers every central phase',
            limit_of_conclusion='Loss of the equality guarantee does not prove the cutoff solution differs or that its growth disappears. This is no stability verdict, no nonlinear constraint solution, and no exclusion of other smooth-data instabilities',
            source='https://arxiv.org/pdf/0806.1036 (Theorem 3.2.11: smooth Cauchy evolution and causal support for normally hyperbolic operators on globally hyperbolic spacetimes; the rays and cutoff application are our derivation)'),
        regularity=dict(critical_nonzero_winding='At s=1 and n!=0 the derivative is bounded but oscillates without a limit; bounded kinetic energy is not smooth crossing',
            repeated_massless_root='n=0,p=0 has independent solutions 1 and log(-S); the logarithmic mode remains',
            massless_critical_zero_winding='n=0,s=1,m=0 has regular solutions 1 and x',
            massive_critical_zero_winding='For n=0,s=1,m>0, sqrt(x)*Y1(2*m*sqrt(x)) contains x*log(x) and a logarithmically divergent derivative; equality does not regularize generic massive zero-winding data',
            massive_check=dict(mass=1., equation='x*G_xx+G=0', samples=mass_exception)),
        source=dict(required='T_ab=(mu/kappa)*(dz)_a*(dz)_b, kappa>0; negative null dust with T_kk=(mu/kappa)*(k^z)^2',
            canonical_sources='Positive-kinetic scalar and ordinary Maxwell null contractions are sums of nonnegative squares. Positive-distribution particle dark matter cannot supply this NEC-violating tensor',
            ghost_limit='Wrong-sign supporting scalars formally provide the sign but introduce ghost degrees of freedom; quantum vacuum stability is a separate unresolved viability condition',
            quantum_limit='A quantum source requires a stress tensor and state on the actual curved periodic background. Boundary-free Minkowski quantum inequalities cannot be applied for arbitrary sampling durations here',
            amount_limit='Dimensionless mu is not a mass in kilograms; no fabrication recipe, attainable quantity, or physical source has been established'),
        casimir=dict(convention='Orthonormal signature -+++; standard zero-temperature ideal flat-space parallel-plate bulk contribution Ccas*diag(-1,1,1,-3), Ccas=pi^2*hbar*c/(720*plate_gap^4)>0',
            ranks={'negative_null': 1, 'ideal_Casimir': 4},
            contracted_squares={'negative_null': '0 for null ell', 'ideal_Casimir': '12*Ccas^2'},
            normalization='ell^a=(1,nx,ny,nz), nx^2+ny^2+nz^2=1, target=-A*ell_a*ell_b, A>0',
            exact_sum_proof='k_plus=(1,1,0,0), k_minus=(1,-1,0,0); Casimir contraction vanishes on each. Target contractions sum to -2*A*(1+nx^2)<0, while any NEC-respecting added stress has a nonnegative sum',
            conclusion='A single ideal Casimir plate pair plus NEC-respecting added stress cannot equal this nonzero negative-null target locally, for any null-dust orientation',
            limit='Assumes the supplied standard flat-space bulk tensor and common local frame, without curvature, finite-size, edge, or plate-material backreaction corrections. This is not a theorem about every actual plate pair in curved geometry; it does not exclude other negative-energy contributions, different plate configurations, or different target metrics'),
        multi_source_casimir=dict(scope='Local algebra of three equal independent ideal Casimir contributions with orthogonal normals, plus positive null-beam radiation, all in one orthonormal -+++ frame; independence and additivity are assumed',
            ideal_input='Each contribution uses the zero-temperature ideal flat-space parallel-plate bulk tensor; their sum is Ccas*diag(-3,-1,-1,-1)',
            target='-A*ell_a*ell_b, A>0, ell^a=(1,0,0,1); a spatial rotation gives any target direction',
            five_beam_example=dict(each_Ccas='2*A', directions=['+x', '-x', '+y', '-y', '-z'], each_positive_weight='A',
                radiation_tensor_over_A=[[str(value/A) for value in five_beams.row(i)] for i in range(4)],
                target_tensor_over_A=[[str(value/A) for value in fixed_target.row(i)] for i in range(4)],
                exact_equality=True),
            four_beam_at_minimum_C=dict(each_Ccas='4*A/3', each_positive_weight='3*A/4',
                contravariant_null_rays=[[str(value) for value in ray] for ray in cone_rays],
                exact_equality=True, clarification='Minimum refers to the coefficient Ccas within this ansatz, not a minimum number of beams'),
            sharp_coefficient_bound=dict(criterion='Ccas>=4*A/3, equivalently 0<A<=3*Ccas/4',
                necessity='A nonnegative sum of null-beam outer products is a positive semidefinite component matrix. The required radiation has its (0,z) principal minor equal to Ccas*(3*Ccas-4*A), which must be nonnegative',
                sufficiency='Use the four listed cone rays with weight 9*Ccas/16 each, plus a forward +z ray of weight 3*Ccas/4-A. All weights are nonnegative exactly in the stated range and the tensor equality is exact',
                limit='Bound is only for three equal orthogonal ideal contributions and nonnegative null-beam radiation; it is not a bound for every canonical matter tensor or plate arrangement'),
            ordinary_matter_checks=dict(all_nine_listed_rays_null_and_future=True,
                NEC='For each positive beam, T_ab*q^a*q^b=weight*(k dot q)^2>=0; verified algebraically',
                DEC='For future timelike u and future null k, k dot u<0, so -T^a_b*u^b=-weight*(k dot u)*k^a is future causal. The flux identity is exact; positive sums preserve the future causal cone'),
            interpretation='The one-pair exclusion does not extend to these multiple independently supplied tensors. This is an algebraic stress budget, not a realizable material or stable-medium result',
            physical_limitations=['No joint quantum state or boundary-value problem for actual intersecting plate pairs has been constructed; their real vacuum stress is not supplied by this algebra',
                'Three plate pairs acting on the same electromagnetic field do not create three independent quantum fields; the matrix sum does not establish additivity of the actual joint boundary problem',
                'No proof that the assumed independent contributions can coexist in a common curved periodic background',
                'No finite-size, edge, curvature, plate-material stress or gravitational backreaction corrections',
                'No regional conservation, state fluctuations, applicable quantum energy restrictions, preparation, nonlinear stability or human-traversal result'],
            source='https://arxiv.org/pdf/hep-th/0702091 (Eq.1, citing Brown-Maclay: each ideal tensor input; the multiple-source decompositions and sharp restricted bound are our algebra)'),
        numerical=dict(method='DOP853 versus exact separated mode; rtol=1e-11, atol=1e-24, 61 samples, winding n=1',
            end_S=-1e-6, relative_error_acceptance=1e-8, cases=rows,
            max_relative_state_error=max(row['max_relative_state_error'] for row in rows),
            deeper_precision_probe=dict(end_S=-1e-8, s=1., n=1, **precision_probe,
                limitation='Small absolute constant-mode contamination is amplified in relative error as the selected mode decays. The reported nine-case relative-accuracy validation stops at S=-1e-6; the threshold itself is exact algebra')),
        conclusion='Formal negative null support makes one selected massless branch satisfy a bounded-kinetic-energy criterion. It establishes neither a stable material nor global or nonlinear stability or human traversal',
        missing=['A physically allowed source realizing the complete stress tensor', 'Quantum state, fluctuations and applicable energy inequalities',
                 'Finite-region matching, preparation and common admissible initial data', 'Nonlinear stability and controlled human passage'],
        sources=['https://arxiv.org/pdf/gr-qc/0503077 (geometry; negative-support branch criterion is our calculation)',
                 'https://arxiv.org/pdf/hep-th/0702091 (Eq.1: ideal Casimir tensor and Brown-Maclay reference)',
                 'https://doi.org/10.1103/PhysRev.184.1272 (Brown-Maclay: conducting-plate vacuum stress)',
                 'https://arxiv.org/pdf/1106.1367 (Vlasov matter energy conditions)',
                 'https://arxiv.org/pdf/hep-ph/0311312 (wrong-sign kinetic ghosts and gravitational vacuum instability)',
                 'https://arxiv.org/pdf/gr-qc/9607003 (Eq.1 and applicability: quantum negative-energy restrictions)'])


def quantum_horizon_geometry_checks():
    """Critical formal source: observer measurements and horizon hypotheses.

    Global statements assume the complete analytic R x R^2 x S^1 extension;
    they do not describe an unspecified finite core or the live PDE solutions.
    """
    from scipy.linalg import expm
    from ori_evolution import axis_geodesic_velocity

    T, x, y, z = coords = sp.symbols('T x y z', real=True)
    mu, w = sp.symbols('mu w', real=True)
    a = sp.Rational(1, 50)
    q = sp.symbols('q', positive=True)
    lam, tau = sp.symbols('lambda tau', real=True)
    F = ((a-mu)*x*x-(a+mu)*y*y)/2
    g = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, w*y], [0, 0, 1, -w*x],
                   [-1, w*y, -w*x, F-T+w*w*(x*x+y*y)]])
    inverse, connection, ricci = geometry(g, coords)
    path = sp.Matrix([q*lam-lam**2/4, 0, 0, -2*sp.log(1-lam/(4*q))])
    velocity = path.diff(lam)
    sub = dict(zip(coords, path))
    metric = g.subs(sub)
    assert sp.simplify((velocity.T*metric*velocity)[0]+1) == 0
    assert sp.simplify((metric*velocity)[3]+q) == 0
    for i in range(4):
        residual = path[i].diff(lam, 2)+sum(connection[i][j][k].subs(sub)
            *velocity[j]*velocity[k] for j in range(4) for k in range(4))
        assert sp.simplify(residual) == 0
    electric = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            for h in (0, 3):
                for k in (0, 3):
                    curvature = sum(g[i, m]*(sp.diff(connection[m][h][k], coords[j])
                        -sp.diff(connection[m][h][j], coords[k])+sum(
                        connection[m][j][v]*connection[v][h][k]
                        -connection[m][k][v]*connection[v][h][j] for v in range(4)))
                        for m in range(4))
                    electric[i, j] += sp.simplify(curvature.subs(sub)*velocity[h]*velocity[k])
    electric = electric.applyfunc(sp.simplify)
    expected = sp.diag(0, (mu-a)*velocity[3]**2/2, (mu+a)*velocity[3]**2/2, 0)
    assert (electric-expected).applyfunc(sp.simplify) == sp.zeros(4)
    rho_kappa = sp.factor((velocity.T*ricci.subs(sub)*velocity)[0])
    tidal_squared = sp.factor(sp.trace(electric*inverse.subs(sub)*electric*inverse.subs(sub)))
    assert sp.simplify(rho_kappa-4*mu/(4*q-lam)**2) == 0
    assert sp.simplify(tidal_squared-8*(mu*mu+a*a)/(4*q-lam)**4) == 0
    observer_error = 0.
    for qq in (.5, 2**-.5, 1.):
        for ll in (-2., -1., -.1, 0.):
            tt = qq*ll-ll*ll/4
            gm = np.array([[0, 0, 0, -1], [0, 1, 0, 0], [0, 0, 1, 0], [-1, 0, 0, -tt]], float)
            dg = np.zeros((4, 4)); dg[3, 3] = -1.
            actual, _ = axis_geodesic_velocity(gm, dg, momentum=-qq)
            exact = np.array([qq-ll/2, 0, 0, 2/(4*qq-ll)])
            observer_error = max(observer_error, float(np.max(abs(actual-exact))))
    assert observer_error < 1e-14

    denominator = 9+16*w*w
    critical_mu = -1-2*a*a/denominator
    v, c = -3*a/denominator, -2*w*a/denominator
    B = sp.Matrix([[-sp.Rational(1, 4)+v/2, c], [c, -sp.Rational(1, 4)-v/2]])
    J = sp.Matrix([[0, -1], [1, 0]])
    Q = 2*B-w*J
    r = sp.Matrix([x, y])
    old = sp.Matrix([tau-(r.T*B*r)[0], x, y, z])
    jacobian = old.jacobian((tau, x, y, z))
    pulled = (jacobian.T*g.subs({T: old[0], mu: critical_mu})*jacobian).applyfunc(sp.simplify)
    qr = Q*r
    expected_metric = sp.Matrix([[0, 0, 0, -1], [0, 1, 0, qr[0]],
        [0, 0, 1, qr[1]], [-1, qr[0], qr[1], (qr.T*qr)[0]-tau]])
    assert (pulled-expected_metric).applyfunc(sp.simplify) == sp.zeros(4)
    assert sp.simplify(pulled.inv()[0, 0]-tau) == 0
    generator = sp.Matrix([0, -qr[0], -qr[1], 1])
    assert sp.simplify((generator.T*pulled*generator)[0]+tau) == 0
    assert (pulled.subs(tau, 0)*generator+sp.Matrix([1, 0, 0, 0])).applyfunc(sp.simplify) == sp.zeros(4, 1)
    assert ((B+sp.eye(2)/4)**2-sp.eye(2)*a*a/(4*denominator)).applyfunc(sp.simplify) == sp.zeros(2)
    assert sp.Rational(1, 2)-a/3 == sp.Rational(37, 75)
    branch_r = sp.symbols('branch_r', negative=True)
    W = 16*w*w
    branch_mu = (1-branch_r**2-16*a*a/(branch_r**2+W))/8
    branch_h = (branch_r+1)/4
    branch_v = a*branch_r/(branch_r**2+W)
    branch_c = -2*a*w/(branch_r**2+W)
    branch_b, branch_d = (branch_h+branch_v)/2, (branch_h-branch_v)/2
    closures = [4*(branch_b**2+branch_c**2)-4*w*branch_c-branch_b-(a-branch_mu)/2,
                4*(branch_d**2+branch_c**2)+4*w*branch_c-branch_d+(a+branch_mu)/2,
                8*branch_c*branch_h+4*w*branch_v-2*branch_c]
    assert all(sp.simplify(value) == 0 for value in closures)
    C = sp.symbols('C', real=True)
    P8 = (8-C)*(8+W)+16*a*a
    assert sp.simplify(P8+(C-9)*(8+W)+8+W-16*a*a) == 0
    assert -8+16*a*a < 0  # C>=9 and W>=0 imply P(8)<0: the larger root U>8.
    branch_eigen_bound = (1-sp.sqrt(8))/8+a/(2*sp.sqrt(8))
    assert branch_eigen_bound < 0
    flow_samples = []
    for ww in (0., .1, 1.):
        flow = -np.array(Q.subs(w, ww), float)
        rate = .5-float(a)/math.sqrt(9+16*ww*ww)
        start = np.array([2., 1.])
        for turns in range(5):
            zz = -2*math.pi*turns
            point = expm(flow*zz)@start
            bound = np.linalg.norm(start)*math.exp(rate*zz)
            assert np.linalg.norm(point) <= bound*(1+1e-12)
            flow_samples.append(dict(w=ww, past_windings=turns,
                transverse_position=point.tolist(), radius=float(np.linalg.norm(point)), radius_bound=bound))
    return dict(
        model='Complete analytic extension R_tau x R^2_transverse x S^1_z, period2*pi; a=1/50, mu=-1-2*a^2/(9+16*w^2), any finite real w',
        observer=dict(path='T=q*lambda-lambda^2/4, x=y=0, z=-2*log(1-lambda/(4*q)); q=-p_z>0, lambda<=0, horizon atlambda0',
            unit_norm_and_four_geodesic_equations=True, conserved_momentum='p_z=-q',
            rho_times_kappa=str(rho_kappa), tidal_squared=str(tidal_squared),
            existing_observer_helper_error=observer_error,
            scope='Physical proper seconds=ell*lambda/c; multiply rho_times_kappa by c^4/(8*pi*G*ell^2). Multiply stored tidal_squared byell^-4, or its square root byell^-2. Both finite atlambda0 for fixedq>0, not human-safe bounds or a quantum-state result'),
        closure=dict(B=[[str(value) for value in B.row(i)] for i in range(2)],
            pulled_back_metric='|dr+Q*r*dz|^2-2*d_tau*dz-tau*dz^2, Q=2B-wJ',
            gradient_tau_squared='tau', generator='L=partial_z-Q*r.partial_r, norm(L)=-tau',
            symmetric_flow_eigenvalue_lower_bound='m=1/2-a/sqrt(9+16*w^2)>=37/75>0',
            backward_bound='|r(z)|<=exp(m*(z-z0))*|r(z0)| forz<=z0', numerical_flow_samples=flow_samples),
        extended_negative_branch=dict(scope='Same complete-extension proof for a=1/50, allmu<=-1, finite realw, choosing r=-sqrt(U) for the larger real root of (U-C)(U+W)+16*a^2=0, C=1-8mu,W=16w^2',
            closure='h=(r+1)/4, v=a*r/(r^2+W), c=-2*a*w/(r^2+W), b=(h+v)/2,d=(h-v)/2; all three closure equations checked exactly',
            eigenvalues_B='(r+1)/8 +/- a/(2*sqrt(r^2+W))',
            bound_argument='C>=9,W>=0 give P(8)<0, so largerrootU>8 and r<-sqrt8',
            largest_eigenvalue_upper_bound=str(branch_eigen_bound),
            largest_eigenvalue_upper_bound_decimal=float(branch_eigen_bound),
            conclusion='B is negative definite throughout this branch. Making mu more negative or changing finite w does not remove its compact generation; this does not classify other branches or extensions'),
        global_hypotheses=dict(initial_surface='Sigma={tau=-1}, complete R^2 x S^1',
            time_orientation='partial_tau is future null; for future causal curves withdz>0, d_tau/dz>=-tau/2, so they cannot reenter tau<0 fromtau>=0',
            negative_slab_estimate='For future causal curves parameterized bytau<0: 0<=dz/dtau<=2/(-tau), |dr/dtau+Q*r*dz/dtau|<=1/sqrt(-tau). Gronwall and finite winding prevent escape in any closed negative-tau slab',
            domain='D(Sigma)={tau<0}, D+(Sigma)={-1<=tau<0}; tau>=0 has past-inextendible constant-tau L curves avoidingSigma',
            horizon='H+(Sigma)={tau=0}; L is its null generator with nonaffinity-1/2',
            compact_generation='All past generators enter and remain in the same K={tau0,|r|<=R0,z modulo2*pi}, anyfixedR0>0, and approach the axis nullcircle',
            limits='These are analytic global-extension arguments aided by exact identities, not a numerical theorem about finite domains, other exterior gluings, or perturbed PDE evolutions'),
        scalar_quantum_obstruction=dict(
            input='Free Hermitian Klein-Gordon field with a distributional two-point bisolution that is initially Hadamard onD(Sigma). R=0 also removes a finite xi*R coupling term from this equation',
            theorem='Kay-Radzikowski-Wald Theorem2/2prime: subtracting a local Hadamard parametrix fails to give a locally L2 function near every horizon base-point diagonal',
            application='The checked complete critical-branch geometry has a compactly generated horizon with base set the central periodic null circle. Renormalized scalar square and stress are singular or ill-defined there',
            limits='Not a theorem of stress divergence along every approaching observer, not automatically a Maxwell/Dirac/interacting-field theorem, and not a prediction of nonlinear backreaction, injury, or all possible time-travel models',
            sources=['https://arxiv.org/pdf/gr-qc/9603012 (Theorems2 and2prime, pp28-29)',
                     'https://arxiv.org/pdf/gr-qc/9606027 (finite limiting stress need not define regular horizon stress)']),
        interpretation='Finite classical observer measurements coexist with a compactly generated Cauchy horizon and the stated scalar quantum obstruction. No source, stable crossing, or calendar-date map has been constructed')


if __name__ == '__main__':
    result = dict(symbolic=symbolic_checks(), numerical=numerical_checks(), twist=twist_checks(),
        maxwell=maxwell_checks(), shifted_null_orbit=shifted_null_orbit_checks(),
        null_screen_jacobi=null_screen_jacobi_checks(),
        twisted_null_screen_jacobi=twisted_null_screen_jacobi_checks(),
        massive_spectator=massive_spectator_checks(),
        negative_null_support=negative_null_support_checks(),
        quantum_horizon_geometry=quantum_horizon_geometry_checks(),
        provenance='Local derivation from Ori 2005 and canonical Einstein-scalar theory; no originality claim',
        scope='Exact coupled support family plus linear spectator test, NOT generic nonlinear evolution',
        missing=['Exterior matching and formation from accessible initial matter',
                 'Nonlinear gravity of the growing spectator',
                 'Generic perturbations, quantum effects and human traversal'])
    (project_path('coupled_ori_results.json')).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: ({key: value for key, value in v.items() if key != 'cases'}
                         if k in ('numerical', 'twist') else v)
                      for k, v in result.items()}, indent=2))
