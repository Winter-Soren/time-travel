"""Third-jet local Dirac A6 calculation for the fixed smooth paired collars.

The coefficient is a numerical estimate, not a short-time heat remainder bound.
No historical geometry, state, source, or counterterm record is changed.
"""
import hashlib
import json
import math
import os
from pathlib import Path
from time import perf_counter

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_key, '1')

import numpy as np
import sympy as sp

from balanced_dihole import symbolic_fields
from dihole_local_source import local_first, mask_jets
from dihole_stretched_collar import exact_target_parameters, gauss, radius_profile
from dihole_heat_cutoff import authenticated_coefficients
from project_paths import project_path


NAME = 'dihole_uv_bound'
ACTIVE = (1, 2)
JET2 = {(1, 1): 3, (1, 2): 4, (2, 1): 4, (2, 2): 5}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collar_logs(function, p, u, theta, radius_order=128):
    """Reuse analytic second jets; all arrays are flattened point arrays."""
    u, theta = np.broadcast_arrays(u, theta)
    radius, first, second = radius_profile(u, order=radius_order)
    x, xr, xt = local_first(function, p, radius, theta)
    tiny = 1e-25
    xrr = np.imag(local_first(function, p, radius+1j*tiny, theta)[1])/tiny
    _, cross, xtt = local_first(function, p, radius, theta+1j*tiny)
    cross, xtt = np.imag(cross)/tiny, np.imag(xtt)/tiny
    jets = np.array([x, xr*first, xt, xrr*first**2+xr*second, cross*first, xtt])
    zero = np.zeros_like(radius)
    log_r = np.array([np.log(radius), first/radius, zero,
                      second/radius-first**2/radius**2, zero, zero])
    log_sin = np.array([np.log(np.sin(theta)), zero, 1/np.tan(theta),
                        zero, zero, -1/np.sin(theta)**2])
    logs = np.stack((jets[:, 0], jets[:, 1], jets[:, 1]+log_r,
                     log_r+log_sin-jets[:, 0]), axis=1)
    return logs, jets[:, 2]*p['physical_A_factor'], mask_jets(u), radius


def prepare_third(function, p, order, step_fraction=2e-4):
    """Fourth-order centered coordinate differences of analytic second jets.

    The mixed third jets are symmetrized only after their independent routes
    have been compared. Difference errors are exposed by step refinement.
    """
    nodes, weights = gauss(order)
    u, theta = np.broadcast_arrays(.025*(nodes[:, None]+1), np.arccos(nodes)[None, :])
    u, theta = u.ravel(), theta.ravel()
    logs, potential, mask, radius = collar_logs(function, p, u, theta)
    du = min(.05*step_fraction, float(min(np.min(u), .05-np.max(u)))/3)
    dt = step_fraction
    derivatives = []
    for direction, step in ((0, du), (1, dt)):
        values = []
        for multiple in (-2, -1, 1, 2):
            row = collar_logs(function, p, u+multiple*step if direction == 0 else u,
                              theta+multiple*step if direction == 1 else theta)
            values.append(np.concatenate((row[0][3:].reshape(12, -1), row[2][3:]), axis=0))
        derivatives.append((values[0]-8*values[1]+8*values[2]-values[3])/(12*step))
    radial, angular = [d[:12].reshape(3, 4, -1) for d in derivatives]
    mismatch = max(float(np.max(abs(radial[1]-angular[0])/np.maximum(1., abs(radial[1])))),
                   float(np.max(abs(radial[2]-angular[1])/np.maximum(1., abs(radial[2])))))
    third = np.stack((radial[0], (2*radial[1]+angular[0])/3,
                      (radial[2]+2*angular[1])/3, angular[2]))
    mask3 = derivatives[0][12]
    measure = (4*np.pi*(.025*weights[:, None]*weights[None, :])).ravel()
    measure *= np.exp(2*logs[0, 1])*radius**2
    return dict(logs=logs, third=third, potential=potential, mask=mask, mask3=mask3,
                measure=measure, order=order, step_fraction=step_fraction,
                coordinate_steps=[du, dt], mixed_third_scaled_mismatch=mismatch)


def diagonal_metric_jets(logs, third):
    n = logs.shape[-1]
    g = np.exp(2*logs[0])
    dg = np.zeros((4, 4, 4, n), dtype=g.dtype)
    ddg = np.zeros((4, 4, 4, 4, n), dtype=g.dtype)
    dddg = np.zeros((4, 4, 4, 4, 4, n), dtype=g.dtype)
    for a in range(4):
        for i in ACTIVE:
            dg[i, a, a] = 2*g[a]*logs[i, a]
            for j in ACTIVE:
                lij = logs[JET2[i, j], a]
                ddg[i, j, a, a] = g[a]*(2*lij+4*logs[i, a]*logs[j, a])
                for k in ACTIVE:
                    index = (i == 2)+(j == 2)+(k == 2)
                    dddg[i, j, k, a, a] = g[a]*(2*third[index, a]
                        +4*(logs[i, a]*logs[JET2[j, k], a]+logs[j, a]*logs[JET2[i, k], a]
                            +logs[k, a]*lij)+8*logs[i, a]*logs[j, a]*logs[k, a])
    return g, dg, ddg, dddg


def curvature(g, dg, ddg, dddg, potential):
    """Coordinate curvature and covariant derivatives, then orthonormalize.

    R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+Gamma^a_ce Gamma^e_db
            -Gamma^a_de Gamma^e_cb. Unit spheres have positive Ricci curvature.
    """
    n, inv = g.shape[-1], 1/g
    connection = np.zeros((4, 4, 4, n), dtype=g.dtype)
    dc = np.zeros((4, 4, 4, 4, n), dtype=g.dtype)
    ddc = np.zeros((4, 4, 4, 4, 4, n), dtype=g.dtype)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                bracket = dg[b, a, c]+dg[c, a, b]-dg[a, b, c]
                connection[a, b, c] = inv[a]*bracket/2
                for h in ACTIVE:
                    bh = ddg[h, b, a, c]+ddg[h, c, a, b]-ddg[h, a, b, c]
                    ih = -inv[a]**2*dg[h, a, a]
                    dc[h, a, b, c] = (ih*bracket+inv[a]*bh)/2
                    for j in ACTIVE:
                        bj = ddg[j, b, a, c]+ddg[j, c, a, b]-ddg[j, a, b, c]
                        bhj = dddg[h, j, b, a, c]+dddg[h, j, c, a, b]-dddg[h, j, a, b, c]
                        ij = -inv[a]**2*dg[j, a, a]
                        ihj = 2*inv[a]**3*dg[h, a, a]*dg[j, a, a]-inv[a]**2*ddg[h, j, a, a]
                        ddc[h, j, a, b, c] = (ihj*bracket+ih*bj+ij*bh+inv[a]*bhj)/2
    product = np.einsum('acep,edbp->abcdp', connection, connection, optimize=True)
    rup = dc.transpose(1, 3, 0, 2, 4)-dc.transpose(1, 3, 2, 0, 4)+product-product.swapaxes(2, 3)
    dproduct = (np.einsum('hacep,edbp->habcdp', dc, connection, optimize=True)
                +np.einsum('acep,hedbp->habcdp', connection, dc, optimize=True))
    drup = (ddc.transpose(0, 2, 4, 1, 3, 5)-ddc.transpose(0, 2, 4, 3, 1, 5)
             +dproduct-dproduct.swapaxes(3, 4))
    rc = g[:, None, None, None]*rup
    grad = g[None, :, None, None, None]*drup
    for h in ACTIVE:
        grad[h] += dg[h, np.arange(4), np.arange(4)][:, None, None, None]*rup
    grad -= np.einsum('ehap,ebcdp->habcdp', connection, rc, optimize=True)
    grad -= np.einsum('ehbp,aecdp->habcdp', connection, rc, optimize=True)
    grad -= np.einsum('ehcp,abedp->habcdp', connection, rc, optimize=True)
    grad -= np.einsum('ehdp,abcep->habcdp', connection, rc, optimize=True)
    frame = 1/np.sqrt(g)
    rf = np.einsum('abcdp,ap,bp,cp,dp->abcdp', rc, frame, frame, frame, frame, optimize=True)
    drf = np.einsum('habcdp,hp,ap,bp,cp,dp->habcdp', grad, frame, frame, frame, frame, frame, optimize=True)
    ric = np.einsum('abadp->bdp', rf)
    dric = np.einsum('habadp->hbdp', drf)
    scalar, dscalar = np.einsum('aap->p', ric), np.einsum('haap->hp', dric)
    magnetic = np.zeros((4, 4, n), dtype=g.dtype)
    dm = np.zeros((4, 4, 4, n), dtype=g.dtype)
    for a in ACTIVE:
        magnetic[a, 3], magnetic[3, a] = potential[a], -potential[a]
        for h in ACTIVE:
            dm[h, a, 3], dm[h, 3, a] = potential[JET2[h, a]], -potential[JET2[h, a]]
    dm -= np.einsum('ehap,ebp->habp', connection, magnetic, optimize=True)
    dm -= np.einsum('ehbp,aep->habp', connection, magnetic, optimize=True)
    magnetic = np.einsum('abp,ap,bp->abp', magnetic, frame, frame, optimize=True)
    dm = np.einsum('habp,hp,ap,bp->habp', dm, frame, frame, frame, optimize=True)
    return dict(riemann=rf, grad_riemann=drf, ricci=ric, grad_ricci=dric,
                scalar=scalar, grad_scalar=dscalar, magnetic=magnetic, grad_magnetic=dm)


def gamma_matrices():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.diag([1., -1.])
    return np.array([np.kron(sx, np.eye(2)), np.kron(sy, sx), np.kron(sy, sy), np.kron(sy, sz)])


def bundle_tensors(t):
    gamma = gamma_matrices()
    bivector = np.array([[(a@b-b@a)/2 for b in gamma] for a in gamma])
    ident = np.eye(4)
    E = (-t['scalar'][None, None]*ident[:, :, None]/4
         -.5j*np.einsum('abij,abp->ijp', bivector, t['magnetic']))
    dE = (-t['grad_scalar'][:, None, None]*ident[None, :, :, None]/4
          -.5j*np.einsum('abij,habp->hijp', bivector, t['grad_magnetic']))
    omega = (.25*np.einsum('abij,abcdp->cdijp', bivector, t['riemann'])
             -1j*t['magnetic'][:, :, None, None]*ident[None, None, :, :, None])
    domega = (.25*np.einsum('abij,habcdp->hcdijp', bivector, t['grad_riemann'])
              -1j*t['grad_magnetic'][:, :, :, None, None]*ident[None, None, None, :, :, None])
    return E, dE, omega, domega


def heat_densities(t):
    """Integrated-by-parts A6 density, before its compact metric variation.

    Vassilevich Eq4.29, translated to our positive-sphere R_abcd convention:
    every odd explicit Riemann occurrence changes sign. Total divergences
    have zero COMPACT VARIATION; no heat smearing f is inserted in this formula.
    """
    e = lambda recipe, *values: np.einsum(recipe, *values, optimize=True)
    r, dr = t['scalar'], t['grad_scalar']
    ric, dri = t['ricci'], t['grad_ricci']
    riem, drm = t['riemann'], t['grad_riemann']
    ric2, riem2 = e('ijx,ijx->x', ric, ric), e('ijklx,ijklx->x', riem, riem)
    scalar_part = (-11*e('ix,ix->x', dr, dr)+6*e('ijkx,ijkx->x', dri, dri)
        -28*e('njkx,kjnx->x', dri, dri)-3*e('hijklx,hijklx->x', drm, drm)
        +35*r**3/9-14*r*ric2/3+14*r*riem2/3
        -208*e('jkx,jnx,knx->x', ric, ric, ric)/9
        +64*e('ijx,klx,ikjlx->x', ric, ric, riem)/3
        -16*e('jkx,jnlix,knlix->x', ric, riem, riem)/3
        +44*e('ijknx,ijabx,knabx->x', riem, riem, riem)/9
        +80*e('ijknx,ilkmx,jlnmx->x', riem, riem, riem)/9)
    E, dE, omega, domega = bundle_tensors(t)
    divergence = e('jijabx->iabx', domega)
    om2 = e('ijabx,ijbax->x', omega, omega)
    Etrace = e('aax->x', E)
    bundle_part = (-4*e('hijabx,hijbax->x', domega, domega)
        +2*e('iabx,ibax->x', divergence, divergence)
        -12*e('ijabx,jkbcx,kicax->x', omega, omega, omega)
        +6*e('ijknx,ijabx,knbax->x', riem, omega, omega)
        -4*e('jkx,jnabx,knbax->x', ric, omega, omega)+5*r*om2
        -30*e('iabx,ibax->x', dE, dE)+60*e('abx,bcx,cax->x', E, E, E)
        +30*e('abx,ijbcx,ijcax->x', E, omega, omega)
        -12*e('hx,haax->x', dr, dE)+30*r*e('abx,bax->x', E, E)
        +Etrace*(5*r*r-2*ric2+2*riem2))
    f2 = e('ijx,ijx->x', t['magnetic'], t['magnetic'])
    return dict(a2=-r/3,
        a4=(5*r*r-8*ric2-7*riem2)/360+2*f2/3,
        a6_identity=4*scalar_part/5040, a6_bundle=bundle_part/360,
        a6=4*scalar_part/5040+bundle_part/360)


def scalar_heat_densities(t):
    """Independent traced reduction: gravity plus the explicitly q^2 part.

    This polynomial in real tensor components also permits a complex-step
    metric variation without imaginary Clifford roundoff contamination.
    """
    e = lambda recipe, *values: np.einsum(recipe, *values, optimize=True)
    r, dr = t['scalar'], t['grad_scalar']
    ric, dri = t['ricci'], t['grad_ricci']
    riem, drm = t['riemann'], t['grad_riemann']
    F, dF = t['magnetic'], t['grad_magnetic']
    ric2, riem2 = e('ijx,ijx->x', ric, ric), e('ijklx,ijklx->x', riem, riem)
    derivative = (19*e('ix,ix->x', dr, dr)-4*e('ijkx,ijkx->x', dri, dri)
                  -84*e('njkx,kjnx->x', dri, dri)+16*e('hijklx,hijklx->x', drm, drm))/5040
    cubic = (-35*r**3+168*r*ric2+147*r*riem2
        -1664*e('jkx,jnx,knx->x', ric, ric, ric)
        +1536*e('ijx,klx,ikjlx->x', ric, ric, riem)
        +120*e('jkx,jnlix,knlix->x', ric, riem, riem)
        -404*e('ijknx,ijabx,knabx->x', riem, riem, riem)
        -872*e('ijknx,ilkmx,jlnmx->x', riem, riem, riem))/90720
    f2 = e('ijx,ijx->x', F, F)
    divF = e('aabx->bx', dF)
    gauge = (-44*e('hijx,hijx->x', dF, dF)-8*e('ix,ix->x', divF, divF)
             -20*r*f2+16*e('abx,acx,bcx->x', ric, F, F)
             +36*e('abcdx,abx,cdx->x', riem, F, F))/360
    gravity = derivative+cubic
    return dict(a2=-r/3, a4=(5*r*r-8*ric2-7*riem2)/360+2*f2/3,
                a6_gravity=gravity, a6_gauge=gauge, a6=gravity+gauge)


def synthetic_curvature(sectional, magnetic=None):
    eye = np.eye(4)
    riem = sectional*(np.einsum('ac,bd->abcd', eye, eye)-np.einsum('ad,bc->abcd', eye, eye))
    riem = riem[..., None]
    ric = np.einsum('abadp->bdp', riem)
    return dict(riemann=riem, grad_riemann=np.zeros((4, 4, 4, 4, 4, 1)),
        ricci=ric, grad_ricci=np.zeros((4, 4, 4, 1)), scalar=np.einsum('aap->p', ric),
        grad_scalar=np.zeros((4, 1)), magnetic=np.zeros((4, 4, 1)) if magnetic is None else magnetic,
        grad_magnetic=np.zeros((4, 4, 4, 1)))


def tensor_controls():
    gamma = gamma_matrices()
    clifford = max(np.max(abs(a@b+b@a-(2*np.eye(4) if i == j else 0)))
                   for i, a in enumerate(gamma) for j, b in enumerate(gamma))
    sphere4 = heat_densities(synthetic_curvature(1.))
    magnetic = np.zeros((4, 4, 1)); magnetic[0, 1, 0] = 1; magnetic[1, 0, 0] = -1
    flat = heat_densities(synthetic_curvature(0., magnetic))
    # Product S2(a) x S2(b): an independent metric-jet/connection control.
    a, b = 1.3, .9
    theta, phi = .71, 1.12
    logs = np.zeros((6, 4, 1)); third = np.zeros((4, 4, 1))
    logs[0, :, 0] = np.log([a*np.sin(theta), a, b, b*np.sin(phi)])
    logs[1, 0, 0], logs[2, 3, 0] = 1/np.tan(theta), 1/np.tan(phi)
    logs[3, 0, 0], logs[5, 3, 0] = -1/np.sin(theta)**2, -1/np.sin(phi)**2
    third[0, 0, 0] = 2*np.cos(theta)/np.sin(theta)**3
    third[3, 3, 0] = 2*np.cos(phi)/np.sin(phi)**3
    t = curvature(*diagonal_metric_jets(logs, third), np.zeros((6, 1)))
    product = heat_densities(t)
    expected = 4*(-(a**-6+b**-6)/252+(a**-2*b**-4+a**-4*b**-2)/360)
    return dict(clifford_error=float(clifford), S4_a6=float(sphere4['a6'].real[0]),
        S4_expected=31/315, S4_error=float(abs(sphere4['a6'][0]-31/315)),
        flat_magnetic_a6=float(flat['a6'].real[0]), flat_magnetic_a4=float(flat['a4'][0]),
        product_a6=float(product['a6'].real[0]), product_expected=expected,
        product_a6_error=float(abs(product['a6'][0]-expected)),
        product_scalar_error=float(abs(t['scalar'][0]-2/a**2-2/b**2)),
        product_parallel_curvature_error=float(np.max(abs(t['grad_riemann']))))


def integrated_coefficients(prepared, epsilon=0., chunk=256):
    values = {key: 0.j for key in ('a2', 'a4', 'a6_gravity', 'a6_gauge', 'a6')}
    sign = np.array([1., -1., 0., 0.])
    for lo in range(0, len(prepared['measure']), chunk):
        at = slice(lo, lo+chunk)
        logs = prepared['logs'][:, :, at]+epsilon*prepared['mask'][:, None, at]*sign[None, :, None]
        third = prepared['third'][:, :, at].astype(np.result_type(logs.dtype, float), copy=True)
        third[0] += epsilon*sign[:, None]*prepared['mask3'][None, at]
        t = curvature(*diagonal_metric_jets(logs, third), prepared['potential'][:, at])
        row = scalar_heat_densities(t)
        for key in values:
            values[key] += row[key]@prepared['measure'][at]
    return values


def paired_symbolic_control():
    epsilon, f = sp.symbols('epsilon f', real=True)
    diagonal = sp.symbols('g0:4', positive=True)
    metric = sp.diag(*(diagonal[i]*sp.exp((2 if i == 0 else -2 if i == 1 else 0)*epsilon*f)
                       for i in range(4)))
    determinant = sp.simplify(metric.det()/sp.prod(diagonal)-1)
    principal = [-sp.diff(1/metric[i, i], epsilon).subs(epsilon, 0) for i in range(4)]
    contracted = sp.simplify(sum(diagonal[i]*principal[i] for i in range(4)))
    return dict(determinant_residual=str(determinant), gaussian_Hessian_residual=str(contracted),
                passed=bool(determinant == 0 and contracted == 0))


def actual_tensor_controls(prepared):
    """Independent rank-four Clifford contraction and charge tests on the collar."""
    tensor = curvature(*diagonal_metric_jets(prepared['logs'], prepared['third']),
                       prepared['potential'])
    scalar = scalar_heat_densities(tensor)
    errors, rows = [], []
    for charge in (0, 1, -1, 2):
        charged = dict(tensor, magnetic=charge*tensor['magnetic'],
                       grad_magnetic=charge*tensor['grad_magnetic'])
        direct = heat_densities(charged)['a6']
        predicted = scalar['a6_gravity']+charge**2*scalar['a6_gauge']
        error = float(np.max(abs(direct-predicted)/np.maximum(1., abs(predicted))))
        errors.append(error)
        rows.append(dict(q=charge, matrix_vs_scalar_polynomial_scaled_error=error))
    dr = tensor['grad_riemann']
    di = tensor['grad_ricci']
    divergence = np.einsum('aabcdp->bcdp', dr)
    left = np.einsum('bcdp,bcdp->p', divergence, divergence)
    right = 2*(np.einsum('abcp,abcp->p', di, di)-np.einsum('njkp,kjnp->p', di, di))
    bianchi = float(np.max(abs(left-right)/np.maximum(1., abs(right))))
    return dict(charge_cases=rows, max_matrix_vs_scalar_scaled_error=max(errors),
        contracted_Bianchi_scaled_error=bianchi,
        covariant_time_curvature_derivative_max=float(np.max(abs(dr[0]))),
        covariant_time_field_derivative_max=float(np.max(abs(tensor['grad_magnetic'][0]))),
        note='All four covariant derivative directions are retained, including static time and azimuth.')


def run():
    start = perf_counter()
    saved, inputs, history = authenticated_coefficients()
    inputs['dihole_heat_cutoff.py'] = sha(project_path('dihole_heat_cutoff.py'))
    function, _ = symbolic_fields()
    parameters = exact_target_parameters()
    symbolic = paired_symbolic_control()
    controls = tensor_controls()
    actual = actual_tensor_controls(prepare_third(function, parameters, 8))
    records = []
    for order, step in ((64, 2e-4), (64, 1e-4), (64, 4e-4), (96, 2e-4), (128, 2e-4)):
        prepared = prepare_third(function, parameters, order, step)
        integral = integrated_coefficients(prepared, 1e-25j)
        records.append(dict(order=order, step_fraction=step,
            coordinate_steps=prepared['coordinate_steps'],
            mixed_third_scaled_mismatch=prepared['mixed_third_scaled_mismatch'],
            derivatives={key: float(value.imag/1e-25) for key, value in integral.items()}))
    final = records[-1]['derivatives']
    scale = max(1., abs(final['a6']))
    quadrature_error = abs(records[-2]['derivatives']['a6']-final['a6'])/scale
    step_error = max(abs(row['derivatives']['a6']-records[0]['derivatives']['a6'])
                     for row in records[1:3])/scale
    recovery = dict(a2_absolute_error=abs(final['a2']-saved['a2_action_derivative']),
                    a4_absolute_error=abs(final['a4']-saved['a4_action_derivative']))
    finite_prepared = prepare_third(function, parameters, 32)
    reference = integrated_coefficients(finite_prepared, 1e-25j)['a6'].imag/1e-25
    differences = []
    for epsilon in (1e-4, 5e-5):
        plus = integrated_coefficients(finite_prepared, epsilon)['a6'].real
        minus = integrated_coefficients(finite_prepared, -epsilon)['a6'].real
        value = (plus-minus)/(2*epsilon)
        differences.append(dict(epsilon=epsilon, derivative=value,
                                error=abs(value-reference)))
    fd_ratio = differences[0]['error']/differences[1]['error']
    # With F_charge=e F, division by 32*pi^2 and e^2=4*pi*alpha gives
    # the independent Drummond--Hathrell local one-loop coefficients.
    dh = [x/8 for x in (sp.Rational(-1, 18), sp.Rational(13, 45),
                        sp.Rational(-1, 45), sp.Rational(-4, 15))]
    dh_expected = [sp.Rational(-1, 144), sp.Rational(13, 360),
                   sp.Rational(-1, 360), sp.Rational(-1, 30)]
    tests = dict(
        exact_paired_volume_and_principal_trace=symbolic['passed'],
        Clifford_relations=bool(controls['clifford_error'] == 0.),
        exact_S4_heat_coefficient=bool(controls['S4_error'] < 2e-14),
        constant_magnetic_heat_control=bool(abs(controls['flat_magnetic_a6']) < 2e-14
                                           and abs(controls['flat_magnetic_a4']-4/3) < 2e-14),
        product_metric_jet_control=bool(controls['product_a6_error'] < 2e-13
                                       and controls['product_parallel_curvature_error'] < 2e-13),
        actual_Clifford_scalar_charge_controls=bool(actual['max_matrix_vs_scalar_scaled_error'] < 1e-10),
        actual_covariant_Bianchi_control=bool(actual['contracted_Bianchi_scaled_error'] < 1e-9),
        saved_lower_coefficients_recovered=bool(max(recovery.values()) < 2e-11),
        angular_radial_quadrature_refinement=bool(quadrature_error < 1e-8),
        third_jet_step_refinement=bool(step_error < 1e-8),
        independent_finite_metric_variation=bool(abs(fd_ratio-4) < .05
                                                 and differences[0]['error']/max(1., abs(reference)) < 1e-6),
        local_gauge_primary_coefficient_control=bool(dh == dh_expected))
    result = dict(name=NAME, passed=bool(all(tests.values())), tests=tests,
        source_sha256=sha(__file__), input_sha256=inputs, historical_source_evidence=history,
        parameters=dict(geometry=parameters, fixed_charge_q=1, heat_cutoff_s0_over_D0_squared=.01,
                        scheme_mu_H_times_D0=.7940236163832829),
        conventions={
            'operator': 'Rank-four massless Euclidean charged Dirac square on the actual prescribed paired smooth collars.',
            'heat': 'Tr exp(-sP) per Killing time ~ (4*pi*s)^(-2) sum_j A_(2j) s^j; reported coefficients exclude (4*pi)^(-2).',
            'variation': 'log(g_tt)/2 += epsilon*f; log(g_uu)/2 -= epsilon*f; magnetic connection fixed; exact pointwise volume preservation.',
            'measure': 'Both collars and complete azimuth: 4*pi*P^2*R^2 du dcos(theta), per unit Killing time.',
            'units': 'D0=1; deltaA2/D0, deltaA4*D0, deltaA6*D0^3 are the reported dimensionless numbers.',
            'curvature': 'Positive-sphere R_abcd=delta_ac delta_bd-delta_ad delta_bc; Vassilevich explicit Riemann signs translated.',
            'derivatives': 'Analytic second coordinate jets; fourth-order centered third jets and all four covariant derivative directions.',
            'IBP': 'Only integrated compact metric variations are used; all derivatives of f vanish at the collar boundaries.'},
        symbolic=symbolic, tensor_controls=controls, actual_tensor_controls=actual,
        numerical_refinement=dict(records=records, relative_last_quadrature_difference=quadrature_error,
            relative_step_difference=step_error, lower_coefficient_recovery=recovery,
            finite_variation=dict(order=32, complex_step_reference=float(reference),
                                  centered_differences=differences, error_ratio=fd_ratio)),
        reference_coefficient_estimate=final,
        fixed_geometry_charge_polynomial=dict(equation='deltaA6(q)=gravity+q^2*gauge',
            gravity=final['a6_gravity'], gauge=final['a6_gauge'],
            interpretation='Local heat coefficient only. Charge changes the Dirac operator; this does not scale its full vacuum response.'),
        local_primary_crosscheck=dict(coefficients_in_alpha_over_pi=[str(x) for x in dh],
            basis=['R F2', 'Ric_ab F_ac F_bc', 'R_abcd F_ab F_cd', '|div F|2']),
        short_time_status=dict(actual_U_bound=None, actual_remainder_constant_Cstar=None,
            formal_leading_signed_subtracted_heat_integral=final['a6']*.01/(16*np.pi**2),
            conditional_formula='If |deltaK(s)-H(s)-deltaA6*s/(16*pi^2)| <= Cstar*s^2/(16*pi^2), then U <= |deltaA6|*s0/(16*pi^2)+Cstar*s0^2/(32*pi^2).',
            missing='A proved quantitative remainder constant on 0<s<=s0, including localization/noncompact effects. The numerical A6 estimate is not an enclosure.'),
        limitations=[
            'This calculation estimates an integrated local heat coefficient; it supplies no sign or bound for the full renormalized source, U, or the infrared tail.',
            'Step and quadrature agreement are numerical checks, not interval enclosures or an asymptotic remainder theorem.',
            'The fixed q=1 geometry is not the nominal one-nanometer atomic candidate and is not a semiclassical solution.',
            'The local q^2 identity holds at fixed geometric shape. Large-q use of s0=.01 is not justified: q*F*s0 and higher heat terms require separate control.',
            'No Maxwell field equation is imposed on the pulled-back collar field when retaining div(F) terms.'],
        primary_sources=[
            dict(url='https://arxiv.org/pdf/hep-th/0306138', role='Vassilevich Eqs3.27-3.28 and4.29: Dirac bundle and integrated A6.'),
            dict(url='https://arxiv.org/pdf/gr-qc/0203034', role='Shore Eq4.10: independent Drummond-Hathrell gauge-coefficient normalization.'),
            dict(url='https://arxiv.org/pdf/1607.05152', role='Ludewig Theorem3.1: heat remainder existence, not an evaluated constant for this geometry.')],
        elapsed_seconds=perf_counter()-start)
    if not result['passed']:
        raise ArithmeticError(json.dumps(tests))
    return result


if __name__ == '__main__':
    output = run()
    destination = project_path(NAME+'_results.json')
    destination.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=output['passed'], checks=len(output['tests']),
        coefficient=output['reference_coefficient_estimate'],
        relative_refinement=output['numerical_refinement']['relative_last_quadrature_difference'],
        source_sha256=output['source_sha256'], results_sha256=sha(destination),
        elapsed_seconds=output['elapsed_seconds']), allow_nan=False))
