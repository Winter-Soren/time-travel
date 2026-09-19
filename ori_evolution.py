"""Local curved Ori Einstein-scalar evolution, not a chronology extension.

Second-order generalized harmonic equations (Lindblom et al., gr-qc/0512093,
Eqs. 4,7,11). Units 8*pi*G=c=1. z is periodic; chi=Phi(t,x,y)*exp(i*z)
represents two real canonical scalars. Metric is independent of z.
Boundary Taylor data are NOT constraint-preserving: only an early interior
with resolution and boundary-sensitivity checks is interpreted physically.
"""
import json
import hashlib
from pathlib import Path
import time

import numpy as np

from disk_basis import basis
from ori_disk_constraints import fields, longitudinal, seed_data, metric_geometry, divergence


H = np.array([0., 0., 0., -1.25])  # Prescribed covector, e=1/16.


def geometry(g, dg, ddg=None):
    """Batch tensors; dg[p,derivative,a,b], ddg[p,derivative,derivative,a,b]."""
    gi = np.linalg.inv(g)
    lower = .5*(dg.transpose(0, 2, 1, 3)+dg.transpose(0, 2, 3, 1)-dg)
    conn = np.einsum('pad,pdbc->pabc', gi, lower)
    dgi = -np.einsum('pam,plmn,pnb->plab', gi, dg, gi)
    dconn = np.einsum('plad,pdbc->plabc', dgi, lower)
    gcov = np.einsum('pbc,pabc->pa', gi, lower)
    dgcov = np.einsum('plbc,pabc->pla', dgi, lower)
    if ddg is not None:
        dlower = .5*(ddg.transpose(0, 1, 3, 2, 4)+ddg.transpose(0, 1, 3, 4, 2)-ddg)
        dconn += np.einsum('pad,pldbc->plabc', gi, dlower)
        dgcov += np.einsum('pbc,plabc->pla', gi, dlower)
    ricci = (np.einsum('pccab->pab', dconn)-np.einsum('pbcac->pab', dconn)+
             np.einsum('pd,pdab->pab', np.einsum('pccd->pd', conn), conn)-
             np.einsum('pcbd,pdac->pab', conn, conn))
    constraint = gcov+H
    covdc = dgcov-np.einsum('pbla,pb->pla', conn, constraint)
    reduced = ricci-.5*(covdc+covdc.transpose(0, 2, 1))
    return gi, conn, gcov, constraint, reduced, ricci


def d1(f, h, axis):
    if isinstance(h, tuple):
        du, inverse_jacobian, _ = h
        shape = [1]*f.ndim
        shape[axis] = len(inverse_jacobian)
        return inverse_jacobian.reshape(shape)*d1(f, du, axis)
    return (np.roll(f, 2, axis)-8*np.roll(f, 1, axis)+
            8*np.roll(f, -1, axis)-np.roll(f, -2, axis))/(12*h)


def d2(f, h, axis):
    if isinstance(h, tuple):
        du, inverse_jacobian, curvature = h
        shape = [1]*f.ndim
        shape[axis] = len(inverse_jacobian)
        return (inverse_jacobian.reshape(shape)**2*d2(f, du, axis)-
                curvature.reshape(shape)*d1(f, du, axis))
    return (-np.roll(f, 2, axis)+16*np.roll(f, 1, axis)-30*f+
            16*np.roll(f, -1, axis)-np.roll(f, -2, axis))/(12*h*h)


def coordinate_grid(n, half_width=.7, stretch=0.):
    if not (half_width > 0 and 0 <= stretch <= 4):
        raise ValueError('Need positive half-width and stretch in [0,4]')
    if stretch == 0:
        grid = np.linspace(-half_width, half_width, n)
        return grid, float(grid[1]-grid[0])
    u = np.linspace(-1., 1., n)
    x = half_width*np.sinh(stretch*u)/np.sinh(stretch)
    J = half_width*stretch*np.cosh(stretch*u)/np.sinh(stretch)
    return x, (float(u[1]-u[0]), 1/J, stretch**2*x/J**3)


def pack(g, v, phi, u):
    # The metric has ten independent components; keep redundant copies equal.
    g = .5*(g+g.swapaxes(-1, -2))
    v = .5*(v+v.swapaxes(-1, -2))
    return np.concatenate((g.reshape(*g.shape[:-2], 16), v.reshape(*v.shape[:-2], 16),
        phi.real[..., None], phi.imag[..., None], u.real[..., None], u.imag[..., None]), axis=-1)


def unpack(state):
    shape = state.shape[:-1]
    return (state[..., :16].reshape(*shape, 4, 4), state[..., 16:32].reshape(*shape, 4, 4),
            state[..., 32]+1j*state[..., 33], state[..., 34]+1j*state[..., 35])


def jets(state, h):
    """All needed derivatives except second time derivatives. Boundary values ignored."""
    shape = state.shape[:2]
    g, v, phi, u = unpack(state)
    dx, dy = d1(state, h, 0), d1(state, h, 1)
    dxx, dyy, dxy = d2(state, h, 0), d2(state, h, 1), d1(dx, h, 1)
    dg = np.zeros((*shape, 4, 4, 4))
    ddg = np.zeros((*shape, 4, 4, 4, 4))
    dg[..., 0, :, :] = v
    for i, first in ((1, dx), (2, dy)):
        dg[..., i, :, :] = first[..., :16].reshape(*shape, 4, 4)
        ddg[..., 0, i, :, :] = ddg[..., i, 0, :, :] = first[..., 16:32].reshape(*shape, 4, 4)
    for i, j, second in ((1, 1, dxx), (2, 2, dyy), (1, 2, dxy), (2, 1, dxy)):
        ddg[..., i, j, :, :] = second[..., :16].reshape(*shape, 4, 4)
    dp = np.stack((u, dx[..., 32]+1j*dx[..., 33], dy[..., 32]+1j*dy[..., 33], 1j*phi), axis=-1)
    ddp = np.zeros((*shape, 4, 4), dtype=complex)
    for i, first in ((1, dx), (2, dy)):
        ddp[..., 0, i] = ddp[..., i, 0] = first[..., 34]+1j*first[..., 35]
        ddp[..., i, 3] = ddp[..., 3, i] = 1j*dp[..., i]
    ddp[..., 0, 3] = ddp[..., 3, 0] = 1j*u
    ddp[..., 3, 3] = -phi
    for i, j, second in ((1, 1, dxx), (2, 2, dyy), (1, 2, dxy), (2, 1, dxy)):
        ddp[..., i, j] = second[..., 32]+1j*second[..., 33]
    return (g.reshape(-1, 4, 4), dg.reshape(-1, 4, 4, 4), ddg.reshape(-1, 4, 4, 4, 4),
            dp.reshape(-1, 4), ddp.reshape(-1, 4, 4))


def low_terms(g, dg):
    """Direct Lindblom Eq.4; full connection-derivative check stays in geometry."""
    gi = np.linalg.inv(g)
    lower = .5*(dg.transpose(0, 2, 1, 3)+dg.transpose(0, 2, 3, 1)-dg)
    gcov = np.einsum('pbc,pabc->pa', gi, lower)
    D = dg.transpose(0, 3, 1, 2)
    low = (np.einsum('pace,pbce->pab', D, gi[:, None]@D@gi[:, None])-
           np.einsum('pace,pbce->pab', lower, gi[:, None]@lower@gi[:, None])+
           np.einsum('pd,pdab->pab', gi@H, lower))
    return gi, gcov, low


def accelerations(g, dg, ddg, dp, ddp, monitor=None):
    spatial = g[:, 1:, 1:]
    # Sylvester's criterion checks the full symmetric spatial metric. Positive
    # inverse diagonals alone can miss a negative eigenvalue.
    minor2 = spatial[:, 0, 0]*spatial[:, 1, 1]-spatial[:, 0, 1]**2
    if (not np.all(np.isfinite(g)) or np.any(spatial[:, 0, 0] <= 0) or
            np.any(minor2 <= 0) or np.any(np.linalg.det(spatial) <= 0)):
        raise ValueError('Spatial metric is not positive definite, or metric nonfinite')
    gi, gcov, low = low_terms(g, dg)
    if np.max(gi[:, 0, 0]) >= 0 or not np.all(np.isfinite(gi)):
        raise ValueError('Spacelike time foliation lost, or metric nonfinite')
    if monitor is not None:
        alpha2 = -1/gi[:, 0, 0]
        beta = -gi[:, 0, 1:]/gi[:, 0, 0, None]
        spatial_inverse_diagonal = np.diagonal(gi[:, 1:, 1:], axis1=1, axis2=2)-gi[:, 0, 1:]**2/gi[:, 0, 0, None]
        if np.min(spatial_inverse_diagonal) <= 0:
            raise ValueError('Spatial inverse metric has nonpositive diagonal')
        speeds = abs(beta)+np.sqrt(alpha2[:, None]*spatial_inverse_diagonal)
        rate = max(np.max(speeds[:, :2]*monitor['inverse_jacobian']),
                   np.max(speeds[:, 2])*monitor['du'])  # k_z=1 retained exactly.
        monitor['rate'] = float(rate)
        monitor['step_rate'] = max(monitor['step_rate'], float(rate))
        monitor['step_speed'] = max(monitor['step_speed'], float(np.max(speeds[:, :2])))
    source = np.einsum('pa,pb->pab', dp.conj(), dp).real
    acc = (2*(low-source)-np.einsum('pcd,pcdab->pab', gi, ddg))/gi[:, 0, 0, None, None]
    gamma = np.einsum('pab,pb->pa', gi, gcov)
    scalar_acc = (np.einsum('pa,pa->p', gamma, dp)-np.einsum('pab,pab->p', gi, ddp))/gi[:, 0, 0]
    return acc, scalar_acc


def rhs(state, h, monitor=None):
    acc, scalar_acc = accelerations(*jets(state, h), monitor=monitor)
    _, v, _, u = unpack(state)
    return pack(v, acc.reshape(*state.shape[:2], 4, 4), u, scalar_acc.reshape(state.shape[:2]))


def background(x, y, t, epsilon=0.):
    """Exact vacuum or constant-intensity null-pair solution; same fixed H."""
    g = np.zeros((*x.shape, 4, 4))
    g[..., 0, 3] = g[..., 3, 0] = -1
    g[..., 1, 1] = g[..., 2, 2] = 1
    g[..., 1, 3] = g[..., 3, 1] = .105*x
    g[..., 2, 3] = g[..., 3, 2] = .145*y
    g[..., 3, 3] = (1/16-epsilon**2/2)*(x*x+y*y)-t
    v = np.zeros_like(g)
    v[..., 3, 3] = -1
    return pack(g, v, np.full(x.shape, epsilon, dtype=complex), np.zeros(x.shape, dtype=complex))


def initial_data(x, y, strong=False, data_path=None):
    """Batch reconstruction limits memory for a large high-degree disk solution."""
    saved = np.load(data_path or Path(__file__).with_name('ori_constraint_coefficients.npz'))
    prefix = 'strong_' if strong else ''
    xf, yf = np.asarray(x).ravel(), np.asarray(y).ravel()
    chunks = [_initial_data_chunk(xf[start:start+256], yf[start:start+256], saved, prefix)
              for start in range(0, len(xf), 256)]
    return np.concatenate(chunks).reshape(*x.shape, 36)


def _initial_data_chunk(x, y, saved, prefix):
    """Reconstruct CTT tensors and solve C_a=0 for initial g_0a velocities."""
    from ori_disk_constraints import scaled_basis
    degree, eps = int(saved[prefix+'degree']), float(saved[prefix+'epsilon'])
    assert float(saved['a']) == .02 and float(saved['e']) == 1/16 and float(saved['initial_t']) == -1
    radius, source_radius = float(saved['radius']), float(saved.get('source_radius', 1.))
    max_angular = int(saved.get(prefix+'max_angular', saved.get('max_angular', -1)))
    max_angular = None if max_angular < 0 else max_angular
    assert radius > 0 and 0 < source_radius <= radius and float(saved['wave_number']) == 1
    xf, yf = x.ravel(), y.ravel()
    if np.max(xf*xf+yf*yf) > radius**2*(1+1e-12):
        raise ValueError('Evolution grid extends outside constraint-data disk')
    data = seed_data(xf, yf, eps, source_radius=source_radius,
                     taper_kind=str(saved.get('taper_kind', 'polynomial')))
    values, first, second = fields(scaled_basis(degree, xf, yf, radius, max_angular=max_angular), saved[prefix+'coefficients'])
    psi, dl = np.exp(values[:, 0]), first[:, :, 0]
    L, _ = longitudinal(data, values[:, 1:], first[:, :, 1:], second[:, :, :, 1:])
    gamma = psi[:, None, None]**4*data['g']
    dgamma = psi[:, None, None, None]**4*(data['dg']+4*dl[:, :, None, None]*data['g'][:, None])
    Kcov = psi[:, None, None]**-2*(data['g']@(data['A']+L)@data['g'])+data['K'][:, None, None]*gamma/3
    D = data['det']
    dD = np.stack((2*(1/16-.105**2)*xf, 2*(1/16-.145**2)*yf, np.zeros_like(xf)), axis=-1)
    alpha = D**-.5
    dalpha = -.5*D[:, None]**-1.5*dD
    beta = np.stack((.105*xf, .145*yf, -np.ones_like(xf)), axis=-1)/D[:, None]
    dbeta = np.diag([.105, .145, 0.])[None]/D[:, None, None]-beta[:, None]*dD[:, :, None]/D[:, None, None]
    g = np.zeros((xf.size, 4, 4))
    dg = np.zeros((xf.size, 4, 4, 4))
    v = np.zeros_like(g)
    g[:, 1:, 1:] = gamma
    g[:, 0, 1:] = g[:, 1:, 0] = np.einsum('pij,pj->pi', gamma, beta)
    g[:, 0, 0] = -alpha**2+np.einsum('pi,pij,pj->p', beta, gamma, beta)
    dg[:, 1:, 1:, 1:] = dgamma
    dg[:, 1:, 0, 1:] = dg[:, 1:, 1:, 0] = (np.einsum('plij,pj->pli', dgamma, beta)+
                                                           np.einsum('pij,plj->pli', gamma, dbeta))
    dg[:, 1:, 0, 0] = -2*alpha[:, None]*dalpha+np.einsum('pi,plij,pj->pl', beta, dgamma, beta)+2*np.einsum('pli,pij,pj->pl', dbeta, gamma, beta)
    lie = np.einsum('pk,pkij->pij', beta, dgamma)+np.einsum('pkj,pik->pij', gamma, dbeta)+np.einsum('pik,pjk->pij', gamma, dbeta)
    v[:, 1:, 1:] = -2*alpha[:, None, None]*Kcov+lie
    dg[:, 0] = v
    gi, _, _, C, _, _ = geometry(g, dg)
    q = -C*np.array([2., 1., 1., 1.])[None]/gi[:, 0, 0, None]
    v[:, 0, :] = v[:, :, 0] = q
    dg[:, 0] = v
    assert np.max(abs(geometry(g, dg)[3])) < 2e-14
    u = alpha*psi**-6*data['Pi']+np.einsum('pi,pi->p', beta, data['C'])
    return pack(g, v, data['field'], u).reshape(*x.shape, 36)


def riemann(g, conn, ddg):
    """Covariant R_abcd; Ricci_bd=g^ac R_abcd, independently checked below."""
    return .5*(ddg.transpose(0, 3, 1, 2, 4)+ddg.transpose(0, 1, 3, 4, 2)-
               ddg.transpose(0, 3, 1, 4, 2)-ddg.transpose(0, 1, 3, 2, 4)) + (
        np.einsum('pmn,pmbc,pnad->pabcd', g, conn, conn)-
        np.einsum('pmn,pmbd,pnac->pabcd', g, conn, conn))


def axis_geodesic_velocity(g, dgdt, momentum=-.5):
    """Unit future velocity and t derivative from conserved p_z on the symmetry axis.

    The caller must check transverse geodesic acceleration. The rational branch
    avoids cancellation at g_zz=0 for momentum having the same sign as g_tz.
    """
    A, B, C = g[0, 0], g[0, 3], g[3, 3]
    dA, dB, dC = dgdt[0, 0], dgdt[0, 3], dgdt[3, 3]
    D, dD = B*B-A*C, 2*B*dB-dA*C-A*dC
    p = momentum
    if not np.all(np.isfinite([A, B, C, D, p])) or D <= 0 or p*p+C <= 0:
        raise ValueError('No regular future t-parametrized axis geodesic for this momentum')
    ut = np.sqrt((p*p+C)/D)
    dut = .5*ut*(dC/(p*p+C)-dD/D)
    if p*B > 0:
        numerator, denominator = -(A*p*p+B*B), D*(p+B*ut)
        uz = numerator/denominator
        dn = -(dA*p*p+2*B*dB)
        dd = dD*(p+B*ut)+D*(dB*ut+B*dut)
        duz = (dn-uz*dd)/denominator
    else:
        if C == 0:
            raise ValueError('Axis velocity has no finite limit on this horizon branch')
        uz = (p-B*ut)/C
        duz = -(dB*ut+B*dut+dC*uz)/C
    return np.array([ut, 0., 0., uz]), np.array([dut, 0., 0., duz])


def observer_tidal_norm(curvature, inverse_metric, observer):
    """Electric Riemann Frobenius norm, with spatiality and sign checks."""
    E = np.einsum('acbd,c,d->ab', curvature, observer, observer)
    orthogonality = float(np.max(abs(E@observer)))
    assert max(orthogonality, np.max(abs(E-E.T))) < 1e-10*(1+np.max(abs(E)))
    # The metric contraction avoids cancellation from a boosted spatial projector.
    E2 = float(np.einsum('ab,ab->', E, inverse_metric@E@inverse_metric))
    if E2 < -1e-12*(1+np.max(abs(E))**2):
        raise ValueError('Observer tidal squared norm is materially negative')
    return float(np.sqrt(max(0., E2))), orthogonality


def diagnostics(state, h, mask):
    """Independent ADM constraints from evolved spatial tensors and FD derivatives."""
    g4, dg4, ddg4, dp, ddp = jets(state, h)
    gi4, conn4, _, C, _, _ = geometry(g4, dg4)
    alpha = (-gi4[:, 0, 0])**-.5
    beta = -gi4[:, 0, 1:]/gi4[:, 0, 0, None]
    gamma, dgamma, ddgamma = g4[:, 1:, 1:], dg4[:, 1:, 1:, 1:], ddg4[:, 1:, 1:, 1:, 1:]
    gi, dgi, conn, _, R = metric_geometry(gamma, dgamma, ddgamma)
    # beta_i=g_0i; covariant derivative uses the 3D spatial connection.
    covdb = dg4[:, 1:, 0, 1:]-np.einsum('pkij,pk->pij', conn, g4[:, 0, 1:])
    Kcov = (-dg4[:, 0, 1:, 1:]+covdb+covdb.transpose(0, 2, 1))/(2*alpha[:, None, None])
    Kup, K = gi@Kcov@gi, np.einsum('pij,pij->p', gi, Kcov)
    shape = state.shape[:2]
    dKup, dK = np.zeros((len(g4), 3, 3, 3)), np.zeros((len(g4), 3))
    for axis in (0, 1):
        dKup[:, axis] = d1(Kup.reshape(*shape, 3, 3), h, axis).reshape(-1, 3, 3)
        dK[:, axis] = d1(K.reshape(shape), h, axis).ravel()
    Pi = (dp[:, 0]-np.einsum('pi,pi->p', beta, dp[:, 1:]))/alpha
    rho = .5*(abs(Pi)**2+np.einsum('pij,pi,pj->p', gi, dp[:, 1:].conj(), dp[:, 1:]).real)
    ham = R+K*K-np.einsum('pij,pij->p', Kcov, Kup)-2*rho
    mom = divergence(conn, Kup, dKup)-np.einsum('pij,pj->pi', gi, dK)+np.einsum('pij,pj->pi', gi, (Pi.conj()[:, None]*dp[:, 1:]).real)
    mnorm = np.sqrt(np.maximum(0, np.einsum('pij,pi,pj->p', gamma, mom, mom)))
    nup = np.column_stack((1/alpha, -beta/alpha[:, None]))
    positive = gi4+2*np.einsum('pa,pb->pab', nup, nup)
    cnorm = np.sqrt(np.maximum(0, np.einsum('pab,pa,pb->p', positive, C, C)))
    acc, _ = accelerations(g4, dg4, ddg4, dp, ddp)
    ddg4[:, 0, 0] = acc
    curvature = riemann(g4, conn4, ddg4)
    tidal = np.einsum('piajb,pa,pb->pij', curvature[:, 1:, :, 1:, :], nup, nup)
    tidal_norm = np.sqrt(np.maximum(0, np.einsum('pik,pjl,pij,pkl->p', gi, gi, tidal, tidal)))
    ricci = np.einsum('pac,pabcd->pbd', gi4, curvature)
    ricci_scalar = np.einsum('pab,pab->p', gi4, ricci)
    ricci_squared = np.einsum('pab,pab->p', ricci, gi4@ricci@gi4)
    source = np.einsum('pa,pb->pab', dp.conj(), dp).real
    source_trace = np.einsum('pab,pab->p', gi4, source)
    source_squared = np.einsum('pab,pab->p', source, gi4@source@gi4)
    sel = mask.ravel()
    center = (shape[0]//2)*shape[1]+shape[1]//2
    regular_norm = g4[center, 0, 0]+2*g4[center, 0, 3]+g4[center, 3, 3]
    regular_energy = (float(abs(dp[center, 0]+dp[center, 3])**2/(-regular_norm)+.5*source_trace[center])
                      if regular_norm < 0 else None)
    regular_tidal, regular_orthogonality = None, None
    if regular_norm < 0:
        observer = np.array([1., 0., 0., 1.])/np.sqrt(-regular_norm)
        regular_tidal, regular_orthogonality = observer_tidal_norm(curvature[center], gi4[center], observer)
    # ponytail: one conserved momentum suffices for the first free-fall comparison.
    momentum = -.5
    freefall, dfreefall = axis_geodesic_velocity(g4[center], dg4[center, 0], momentum)
    geodesic_acceleration = freefall[0]*dfreefall+np.einsum('abc,b,c->a', conn4[center], freefall, freefall)
    geodesic_residual = float(np.max(abs(geodesic_acceleration)))
    freefall_identity = max(abs(freefall@g4[center]@freefall+1), abs(g4[center, 3]@freefall-momentum))
    assert freefall_identity < 1e-10, freefall_identity
    if geodesic_residual > 1e-9*(1+np.max(abs(conn4[center]))*np.max(abs(freefall))**2):
        raise ValueError('Central axis does not satisfy the assumed geodesic equations')
    freefall_tidal, _ = observer_tidal_norm(curvature[center], gi4[center], freefall)
    freefall_energy = float(abs(dp[center]@freefall)**2+.5*source_trace[center])
    speed = np.maximum(abs(beta[:, 0])+alpha*np.sqrt(gi[:, 0, 0]), abs(beta[:, 1])+alpha*np.sqrt(gi[:, 1, 1]))
    return dict(max_H=float(np.max(abs(ham[sel]))), max_M=float(np.max(mnorm[sel])),
        max_gauge=float(np.max(cnorm[sel])), max_energy=float(np.max(rho[sel])),
        center_energy=float(rho[center]), center_gzz=float(g4[center, 3, 3]),
        center_tidal_norm=float(tidal_norm[center]), max_tidal_norm=float(np.max(tidal_norm[sel])),
        center_R=float(ricci_scalar[center]),
        center_Ricci_squared=float(ricci_squared[center]), center_source_R=float(source_trace[center]),
        center_source_Ricci_squared=float(source_squared[center]),
        max_einstein_component_residual=float(np.max(abs(ricci[sel]-source[sel]))),
        center_regular_energy=regular_energy, center_regular_direction_norm=float(regular_norm),
        center_regular_tidal_norm=regular_tidal, center_regular_tidal_orthogonality=regular_orthogonality,
        center_freefall_pz=momentum, center_freefall_ut=float(freefall[0]), center_freefall_uz=float(freefall[3]),
        center_freefall_energy=freefall_energy, center_freefall_tidal_norm=freefall_tidal,
        center_freefall_geodesic_residual=geodesic_residual, center_freefall_identity_residual=float(freefall_identity),
        center_gtt_inverse=float(gi4[center, 0, 0]), max_coordinate_speed=float(speed.max()),
        min_spatial_eigenvalue=float(np.linalg.eigvalsh(gamma).min()))


def evolve(n, end=-.8, boundary_order=2, strong=False, exact=None, cfl=.15,
           half_width=.7, stretch=0., data_path=None, adaptive=False, checkpoint_path=None):
    if n < 9 or n % 2 != 1 or not (-1 < end < 0) or boundary_order not in (1, 2) or not (0 < cfl <= .3):
        raise ValueError('Use odd n>=9, -1<end<0, Taylor order 1 or 2, and 0<cfl<=.3')
    start = time.perf_counter()
    # ponytail: checkpoints cover adaptive long runs; extend only if fixed-step resumes are needed.
    if checkpoint_path is not None and not adaptive:
        raise ValueError('Checkpointing requires adaptive=True')
    grid, h = coordinate_grid(n, half_width, stretch)
    du = h[0] if isinstance(h, tuple) else h
    x, y = np.meshgrid(grid, grid, indexing='ij')
    interior = x*x+y*y <= .2**2
    wall = np.ones(x.shape, dtype=bool)
    wall[2:-2, 2:-2] = False
    initial = initial_data(x, y, strong, data_path) if exact is None else background(x, y, -1, exact)
    # Boundary acceleration needs valid initial stencils; use the nearest
    # interior row. This is tested for interior sensitivity, not a CPBC.
    rows = np.clip(np.arange(n), 2, n-3)
    g0, v0, phi0, u0 = unpack(initial)
    _, a0, _, b0 = unpack(rhs(initial, h))
    a0, b0 = a0[rows[:, None], rows[None, :]], b0[rows[:, None], rows[None, :]]
    if boundary_order == 1:
        a0, b0 = np.zeros_like(a0), np.zeros_like(b0)
    def apply_wall(state, t):
        if exact is not None:
            state[wall] = background(x[wall], y[wall], t, exact)
        else:
            delta = t+1
            continued = pack(g0+delta*v0+.5*delta**2*a0, v0+delta*a0,
                             phi0+delta*u0+.5*delta**2*b0, u0+delta*b0)
            state[wall] = continued[wall]
        return state
    state = initial.copy()
    samples = [dict(t=-1., **diagnostics(state, h, interior))]
    fixed_steps = int(np.ceil((end+1)/(cfl*du)))
    fixed_dt = (end+1)/fixed_steps
    sample_times = (list(np.linspace(-1, end, 11)[1:]) if adaptive else
                    [-1+(fixed_steps//2)*fixed_dt, end])
    invj = h[1] if isinstance(h, tuple) else np.ones(n)
    monitor = dict(du=du, inverse_jacobian=np.stack(np.meshgrid(invj, invj, indexing='ij'), axis=-1).reshape(-1, 2),
                   rate=0., step_rate=0., step_speed=0.)
    t, steps, travel, max_courant, min_dt, max_dt = -1., 0, 0., 0., np.inf, 0.
    rejected_steps = 0
    previous_seconds = 0.
    checkpoint = Path(checkpoint_path) if checkpoint_path is not None else None
    if checkpoint is not None:
        case = dict(n=n, end=end, boundary_order=boundary_order, strong=strong, exact=exact,
            cfl=cfl, half_width=half_width, stretch=stretch,
            initial_sha256=hashlib.sha256(initial.tobytes()).hexdigest(),
            evolution_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        if checkpoint.exists():
            with np.load(checkpoint) as saved:
                meta = json.loads(str(saved['metadata']))
                if meta['case'] != case:
                    raise ValueError('Checkpoint case, initial state, or evolution code differs')
                state = saved['state'].copy()
            if state.shape != initial.shape or not np.all(np.isfinite(state)) or not (-1 < meta['t'] <= end):
                raise ValueError('Invalid checkpoint state or time')
            t, steps, travel, max_courant, min_dt, max_dt, rejected_steps = (
                meta[key] for key in ('t', 'steps', 'travel', 'max_courant', 'min_dt', 'max_dt', 'rejected_steps'))
            samples, previous_seconds = meta['samples'], meta['seconds']
            if steps < 1 or len(samples) < 2 or samples[-1]['t'] != t:
                raise ValueError('Invalid checkpoint progress metadata')
            sample_times = [value for value in sample_times if value > t+1e-13]
            print(json.dumps(dict(resumed_from=str(checkpoint), t=t, steps=steps)), flush=True)
    while t < end-1e-14:
        monitor['step_rate'], monitor['step_speed'] = 0., 0.
        k1 = rhs(apply_wall(state, t), h, monitor)
        initial_rate, initial_speed = monitor['rate'], monitor['step_speed']
        dt = min(.95*cfl*du/initial_rate, sample_times[0]-t, end-t) if adaptive else fixed_dt
        for attempt in range(16):
            monitor['step_rate'], monitor['step_speed'] = initial_rate, initial_speed
            k2 = rhs(apply_wall(state+.5*dt*k1, t+.5*dt), h, monitor)
            k3 = rhs(apply_wall(state+.5*dt*k2, t+.5*dt), h, monitor)
            k4 = rhs(apply_wall(state+dt*k3, t+dt), h, monitor)
            courant = dt*monitor['step_rate']/du
            if not adaptive or courant <= cfl*(1+1e-12):
                break
            rejected_steps += 1
            dt = min(.5*dt, .95*cfl*du/monitor['step_rate'])
        else:
            raise ValueError('Adaptive RK4 stage CFL failed after 16 attempts')
        state = apply_wall(state+dt*(k1+2*k2+2*k3+k4)/6, t+dt)
        t, steps = t+dt, steps+1
        travel += dt*monitor['step_speed']
        max_courant = max(max_courant, dt*monitor['step_rate']/du)
        min_dt, max_dt = min(min_dt, dt), max(max_dt, dt)
        if t >= sample_times[0]-1e-13:
            sample = dict(t=t, sampled_cartesian_travel_envelope=travel,
                          **diagnostics(state, h, interior))
            samples.append(sample)
            sample_times.pop(0)
            if checkpoint is not None:
                meta = dict(case=case, t=t, steps=steps, travel=travel, max_courant=max_courant,
                    min_dt=min_dt, max_dt=max_dt, rejected_steps=rejected_steps, samples=samples,
                    seconds=previous_seconds+time.perf_counter()-start)
                temporary = checkpoint.with_name(checkpoint.name+'.tmp')
                with temporary.open('wb') as stream:
                    np.savez_compressed(stream, state=state, metadata=json.dumps(meta))
                temporary.replace(checkpoint)
            if adaptive:
                print(json.dumps(dict(n=n, strong=strong, **sample)), flush=True)
    result = dict(n=n, dt=fixed_dt if not adaptive else None, min_dt=min_dt, max_dt=max_dt,
        steps=steps, adaptive=adaptive, cfl=cfl, stage_cfl_enforced=adaptive,
        rejected_steps=rejected_steps, spatial_stage_guard='Sylvester principal minors',
        max_stage_courant=max_courant,
        sampled_cartesian_travel_envelope=travel, nearest_wall_to_core=float(grid[-2]-.2),
        boundary_order=boundary_order, strong=strong,
        exact_epsilon=exact, interior_radius=.2, half_width=half_width, stretch=stretch,
        initial_data_file=str(data_path) if data_path else 'ori_constraint_coefficients.npz', samples=samples,
        seconds=previous_seconds+time.perf_counter()-start,
        checkpoint_file=str(checkpoint) if checkpoint is not None else None)
    if exact is not None:
        result['exact_error'] = float(np.max(abs(state[4:-4, 4:-4]-background(x, y, end, exact)[4:-4, 4:-4])))
    return result, state


def checkpoint_check():
    """Interrupt after a saved step; resumption must reproduce an uninterrupted run."""
    import contextlib
    import io
    import tempfile
    global rhs
    original_rhs = rhs
    with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
        path = Path(directory)/'resume.npz'
        options = dict(n=9, end=-.95, exact=.3, adaptive=True)
        reference, expected = evolve(**options)
        def interrupted(*args, **kwargs):
            if path.exists():
                raise RuntimeError('Intentional checkpoint-test interruption')
            return original_rhs(*args, **kwargs)
        try:
            rhs = interrupted
            try:
                evolve(**options, checkpoint_path=path)
            except RuntimeError as error:
                assert str(error) == 'Intentional checkpoint-test interruption', error
            else:
                raise AssertionError('Checkpoint interruption was not exercised')
        finally:
            rhs = original_rhs
        resumed, actual = evolve(**options, checkpoint_path=path)
        assert np.array_equal(actual, expected)
        for key in ('samples', 'steps', 'min_dt', 'max_dt', 'rejected_steps',
                    'max_stage_courant', 'sampled_cartesian_travel_envelope'):
            assert resumed[key] == reference[key], key
        try:
            evolve(**dict(options, cfl=.1), checkpoint_path=path)
        except ValueError as error:
            assert 'differs' in str(error), error
        else:
            raise AssertionError('Incompatible checkpoint was accepted')
    return True


def self_check():
    rng = np.random.default_rng(8427)
    g = np.diag([-1., 1., 1., 1.])[None]+.05*rng.normal(size=(3, 4, 4))
    g = .5*(g+g.transpose(0, 2, 1))
    dg = .1*rng.normal(size=(3, 4, 4, 4))
    dg = .5*(dg+dg.transpose(0, 1, 3, 2))
    ddg = .1*rng.normal(size=(3, 4, 4, 4, 4))
    ddg = .5*(ddg+ddg.transpose(0, 2, 1, 3, 4))
    ddg = .5*(ddg+ddg.transpose(0, 1, 2, 4, 3))
    gi, _, _, _, low, _ = geometry(g, dg)
    direct_error = float(np.max(abs(low_terms(g, dg)[2]-low)))
    assert direct_error < 2e-14, direct_error
    full = geometry(g, dg, ddg)[4]
    err = float(np.max(abs(full-low+.5*np.einsum('pcd,pcdab->pab', gi, ddg))))
    assert err < 2e-14, err
    connection, ricci = geometry(g, dg, ddg)[1], geometry(g, dg, ddg)[5]
    curvature = riemann(g, connection, ddg)
    curvature_error = float(np.max(abs(np.einsum('pac,pabcd->pbd', gi, curvature)-ricci)))
    assert curvature_error < 2e-14, curvature_error
    checks = dict(principal_identity_error=err, independent_ricci_error=curvature_error,
                  direct_low_error=direct_error)
    x, y = np.meshgrid(np.linspace(-.7, .7, 25), np.linspace(-.7, .7, 25), indexing='ij')
    for eps in (0., .3):
        state = background(x, y, -1, eps)
        a, b = accelerations(*jets(state, 1.4/24))
        actual = max(float(np.max(abs(a.reshape(25, 25, 4, 4)[4:-4, 4:-4]))),
                     float(np.max(abs(b.reshape(25, 25)[4:-4, 4:-4]))))
        assert actual < 1e-11, actual
        checks[f'exact_acceleration_{eps}'] = actual
    # Roundoff in redundant tensor entries must not seed antisymmetric modes.
    raw_g, raw_v, phi, u = unpack(state.copy())
    raw_g[..., 0, 1] += 1e-12
    raw_v[..., 0, 2] -= 2e-12
    phi, u = phi+.17j, u+.2+.3j
    rebuilt = pack(raw_g, raw_v, phi, u)
    actual_g, actual_v, actual_phi, actual_u = unpack(rebuilt)
    assert all(np.array_equal(value, value.swapaxes(-1, -2)) for value in (actual_g, actual_v))
    assert np.array_equal(actual_phi, phi) and np.array_equal(actual_u, u)
    assert all(np.array_equal(value, value.swapaxes(-1, -2)) for value in unpack(rhs(rebuilt, 1.4/24))[:2])
    checks['metric_and_rhs_symmetry_preserved'] = True
    # Independent exact pp wave: eta+F(dt+dz)^2, chi=eps*exp(i(t+z)).
    # Exercises nonzero mixed scalar time/z derivatives and a second gauge.
    x, y, eps, a = np.array([0., .2]), np.array([0., -.1]), .3, .07
    ell = np.array([1., 0., 0., 1.])
    ll = np.outer(ell, ell)
    F = a*(x*x-y*y)-eps**2*(x*x+y*y)/2
    metric = np.diag([-1., 1., 1., 1.])[None]+F[:, None, None]*ll
    first, second = np.zeros((2, 4, 4, 4)), np.zeros((2, 4, 4, 4, 4))
    first[:, 1] = ((2*a-eps**2)*x)[:, None, None]*ll
    first[:, 2] = ((-2*a-eps**2)*y)[:, None, None]*ll
    second[:, 1, 1] = (2*a-eps**2)*ll
    second[:, 2, 2] = (-2*a-eps**2)*ll
    field = eps*np.exp(.4j)
    dp = np.broadcast_to(1j*field*ell, (2, 4))
    ddp = np.broadcast_to(-field*ll, (2, 4, 4)).copy()
    ddp[:, 0, 0] = 0
    original_H = H.copy()
    try:
        H[:] = 0
        _, connection, _, C, _, ricci = geometry(metric, first, second)
        acc, scalar_acc = accelerations(metric, first, second, dp, ddp)
    finally:
        H[:] = original_H
    curvature = riemann(metric, connection, second)
    expected_tidal = np.diag([-a+eps**2/2, a+eps**2/2, 0.])
    benchmark_error = float(max(np.max(abs(ricci-eps**2*ll)), np.max(abs(C)), np.max(abs(acc)),
        np.max(abs(scalar_acc+field)), np.max(abs(curvature[0, 1:, 0, 1:, 0]-expected_tidal))))
    assert benchmark_error < 1e-13, benchmark_error
    checks['independent_pp_wave_error'] = benchmark_error
    mapping_errors = []
    for n in (33, 65, 129):
        grid, h = coordinate_grid(n, 2.1, 2.5)
        x, y = np.meshgrid(grid, grid, indexing='ij')
        cx, cy = .1+.7j, 1.1j
        f = np.exp(cx*x+cy*y)
        mask = (abs(x) <= 1.5) & (abs(y) <= 1.5)
        errors = [d1(f, h, 0)-cx*f, d2(f, h, 0)-cx*cx*f,
                  d1(d1(f, h, 0), h, 1)-cx*cy*f]
        mapping_errors.append(float(max(np.max(abs(error[mask])) for error in errors)))
    assert all(mapping_errors[i]/mapping_errors[i+1] > 10 for i in (0, 1)), mapping_errors
    checks['mapped_derivative_errors'] = mapping_errors
    # Normal observers become boosted at this regular null-pair horizon.
    # Their divergent energy alone must not be called a curvature singularity.
    x, y = np.meshgrid(np.linspace(-.04, .04, 9), np.linspace(-.04, .04, 9), indexing='ij')
    mask = (x == 0) & (y == 0)
    observer_errors = []
    freefall_errors = []
    for t in (-1., -.1):
        sample = diagnostics(background(x, y, t, .3), .01, mask)
        observer_errors.extend([abs(sample['center_energy']-.09/(-t)),
                                abs(sample['center_regular_energy']-.09/(2+t)),
                                abs(sample['center_regular_tidal_norm']-np.sqrt((.02**2+.09**2)/2)/(2+t))])
        uz = 1/(np.sqrt(.25-t)+.5)
        freefall_errors.extend([abs(sample['center_freefall_energy']-.09*uz**2),
            abs(sample['center_freefall_tidal_norm']-np.sqrt((.02**2+.09**2)/2)*uz**2),
            sample['center_freefall_geodesic_residual'], sample['center_freefall_identity_residual']])
    for t in (-1e-12, -1e-30, 0.):
        g = np.array([[0., 0., 0., -1.], [0., 1., 0., 0.], [0., 0., 1., 0.], [-1., 0., 0., -t]])
        dg = np.zeros((4, 4, 4))
        dg[0, 3, 3] = -1
        observer, derivative = axis_geodesic_velocity(g, dg[0])
        connection = geometry(g[None], dg[None])[1][0]
        acceleration = observer[0]*derivative+np.einsum('abc,b,c->a', connection, observer, observer)
        freefall_errors.extend([float(np.max(abs(acceleration))), abs(observer@g@observer+1),
            abs(g[3]@observer+.5), abs(observer[3]-1/(np.sqrt(.25-t)+.5))])
    assert max(freefall_errors) < 1e-12, freefall_errors
    checks['freefall_geodesic_and_horizon_identity_error'] = float(max(freefall_errors))
    assert max(observer_errors) < 1e-12, observer_errors
    checks['observer_energy_identity_error'] = float(max(observer_errors))
    # This indefinite matrix has positive inverse diagonals; the full guard
    # must reject it before using characteristic speeds.
    bad = np.diag([-1., 1., 1., 1.])[None]
    bad[0, 1:, 1:] = [[1., 2., 2.], [2., 1., -2.], [2., -2., 1.]]
    try:
        accelerations(bad, np.zeros((1, 4, 4, 4)), np.zeros((1, 4, 4, 4, 4)),
                      np.zeros((1, 4), dtype=complex), np.zeros((1, 4, 4), dtype=complex))
    except ValueError as error:
        assert 'positive definite' in str(error), error
    else:
        raise AssertionError('Indefinite spatial metric passed the stage guard')
    checks['indefinite_spatial_metric_rejected'] = True
    checks['checkpoint_resume_exact'] = checkpoint_check()
    return checks


def main():
    checks, runs, states = self_check(), {}, {}
    output = Path(__file__).with_name('ori_evolution_results.json')
    def run(name, **kwargs):
        record, state = evolve(**kwargs)
        runs[name], states[name] = record, state
        print(json.dumps(dict(name=name, **record)), flush=True)
        output.write_text(json.dumps(dict(scope='Early interior GH evolution; no chronology extension',
            checks=checks, runs=runs), indent=2)+'\n')
    for eps in (0., .3):
        run(f'exact_{eps}', n=25, exact=eps)
    for strong in (False, True):
        label = 'strong' if strong else 'weak'
        for n in (25, 49, 97):
            run(f'{label}_{n}', n=n, strong=strong)
        run(f'{label}_linear', n=97, strong=strong, boundary_order=1)
    run('weak_half_dt', n=49, cfl=.075)
    comparisons = {}
    for label in ('weak', 'strong'):
        distances = []
        for coarse, fine in ((25, 49), (49, 97)):
            grid = np.linspace(-.7, .7, coarse)
            x, y = np.meshgrid(grid, grid, indexing='ij')
            mask = x*x+y*y <= .2**2
            distances.append(float(np.max(abs(states[f'{label}_{coarse}'][mask]-
                                             states[f'{label}_{fine}'][::2, ::2][mask]))))
        grid = np.linspace(-.7, .7, 97)
        x, y = np.meshgrid(grid, grid, indexing='ij')
        mask = x*x+y*y <= .2**2
        boundary = float(np.max(abs(states[f'{label}_97'][mask]-states[f'{label}_linear'][mask])))
        final = runs[f'{label}_97']['samples'][-1]
        comparisons[label] = dict(nested_state_differences=distances,
            refinement_ratio=distances[0]/distances[1], boundary_state_difference=boundary,
            center_energy_ratio=final['center_energy']/runs[f'{label}_97']['samples'][0]['center_energy'])
    grid = np.linspace(-.7, .7, 49)
    x, y = np.meshgrid(grid, grid, indexing='ij')
    mask = x*x+y*y <= .2**2
    comparisons['time_step_difference'] = float(np.max(abs(states['weak_49'][mask]-states['weak_half_dt'][mask])))
    checks.update(exact_evolutions=all(runs[f'exact_{eps}']['exact_error'] < 2e-11 for eps in (0., .3)),
        spatial_convergence=all(10 < comparisons[label]['refinement_ratio'] < 22 for label in ('weak', 'strong')),
        constraint_convergence=all(runs[f'{label}_{coarse}']['samples'][-1][key] >
            8*runs[f'{label}_{fine}']['samples'][-1][key] for label in ('weak', 'strong')
            for coarse, fine in ((25, 49), (49, 97)) for key in ('max_H', 'max_M', 'max_gauge')),
        small_boundary_sensitivity=all(comparisons[label]['boundary_state_difference'] <
                                      comparisons[label]['nested_state_differences'][1]/20 for label in ('weak', 'strong')),
        time_error_subdominant=comparisons['time_step_difference'] < comparisons['weak']['nested_state_differences'][1]/10,
        weak_constraints=max(runs['weak_97']['samples'][-1][key] for key in ('max_H', 'max_M', 'max_gauge')) < 1e-7,
        strong_constraints=max(runs['strong_97']['samples'][-1][key] for key in ('max_H', 'max_M', 'max_gauge')) < 2e-4,
        positive_spatial_metric=all(sample['min_spatial_eigenvalue'] > 0 for run in runs.values() for sample in run['samples']))
    output.write_text(json.dumps(dict(scope='Early interior GH evolution; no chronology extension',
        checks=checks, comparisons=comparisons, runs=runs), indent=2)+'\n')
    np.savez_compressed(Path(__file__).with_name('ori_evolution_states.npz'), weak=states['weak_97'],
        strong=states['strong_97'], t=-.8, half_width=.7, interior_radius=.2)
    print(json.dumps(dict(checks=checks, comparisons=comparisons)), flush=True)
    assert all(value for value in checks.values() if isinstance(value, bool)), checks


if __name__ == '__main__':
    main()
