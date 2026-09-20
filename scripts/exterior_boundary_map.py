"""Connected flat two-mouth Weyl problem with a conservative P1 DG flux.

The finite-domain characteristic boundary is NOT a transparent exterior or a
physical escape channel.  Numerical jump loss is reported separately.  The
weighted unknown u=f/sqrt(rho) has measure rho*d_rho*d_z.  This module owns no
geometry, throat model, or vacuum-source update.
"""
from __future__ import annotations
from project_paths import project_path

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.special import jv

from incoming_modes import collar_mode


SX = np.array([[0., 1.], [1., 0.]], complex)
SY = np.array([[0., -1j], [1j, 0.]], complex)
SZ = np.diag([1., -1.]).astype(complex)
ID = np.eye(2, dtype=complex)


def pair_potential(rho, z, q=1., d=1.):
    """A_phi in the existing north-patch convention, away from the two foci."""
    return q/2*((z+d/2)/np.hypot(rho, z+d/2)
                -(z-d/2)/np.hypot(rho, z-d/2))


def _quadrature(order=4):
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x+1)/2, w/2
    # Duffy rule: weights sum to 1/2 on the reference triangle.
    bary, weights = [], []
    for s, ws in zip(x, w):
        for t, wt in zip(x, w):
            bary.append([1-s, s*(1-t), s*t])
            weights.append(ws*wt*s)
    return np.array(bary), np.array(weights), x, w


class ExteriorDG:
    """Strong-form weighted P1 upwind discretization; one LU, multiple inputs.

    An optional feedback triple (G, gamma, C) supplies a separately justified
    incoming lead map: rhs -> rhs + G diag(gamma) C u.  C acts on FULL DG
    coefficients, G has full coefficient rows.  No transparency is inferred.
    """

    def __init__(self, mesh, q=1, mu=.5, energy=.4, d=1.):
        if not all(math.isfinite(float(v)) for v in (q, mu, energy, d)) or d <= 0:
            raise ValueError('Finite coefficients and positive separation required')
        if abs(mu*2-round(mu*2)) > 1e-12 or round(mu*2) % 2 == 0:
            raise ValueError('The physical azimuthal label mu must be half-integral')
        if isinstance(mesh, dict):
            vertices, triangles, boundary = (mesh[k] for k in
                                             ('vertices', 'triangles', 'boundary_edges'))
        else:
            vertices, triangles, boundary = mesh
        self.vertices = np.asarray(vertices, float)
        self.triangles = np.asarray(triangles, int)
        self.q, self.mu, self.energy, self.d = q, mu, energy, d
        self.ndof = 6*len(self.triangles)
        self.boundary, self.interior, self.cells = [], [], []
        bary, qw, ex, ew = _quadrature()
        rows, cols, values = [], [], []

        def block(tl, tr, matrix):
            ii = 6*tl + np.arange(6)
            jj = 6*tr + np.arange(6)
            rows.extend(np.repeat(ii, 6))
            cols.extend(np.tile(jj, 6))
            values.extend(matrix.reshape(-1))

        edge_owners = {}
        forbidden = []
        for it, tri in enumerate(self.triangles):
            xy = self.vertices[tri]
            jac = np.linalg.det(np.stack((xy[1]-xy[0], xy[2]-xy[0])))
            if jac <= 0:
                raise ValueError('Triangles must be nondegenerate and counterclockwise')
            grad = np.linalg.inv(np.column_stack((np.ones(3), xy)))[1:].T
            points, weights = bary@xy, qw*jac
            rho, z = points.T
            weighted = weights*rho
            mass = np.einsum('p,pa,pb->ab', weighted, bary, bary)
            moments = np.einsum('p,pa->a', weighted, bary)
            plain_mass = np.einsum('p,pa,pb->ab', weights, bary, bary)
            azimuth = mu-pair_potential(rho, z, q, d)
            angular = np.einsum('p,pa,pb->ab', weights*azimuth, bary, bary)
            local = np.zeros((3, 2, 3, 2), complex)
            for a in range(3):
                for b in range(3):
                    local[a, :, b, :] = (moments[a]*(grad[b, 0]*SX+grad[b, 1]*SZ)
                                          +plain_mass[a, b]*SX/2
                                          +1j*angular[a, b]*SY
                                          -1j*energy*mass[a, b]*ID)
            block(it, it, local.reshape(6, 6))
            self.cells.append(dict(triangle=it, points=points, weights=weighted,
                                   basis=bary, mass=mass))
            for a, vertex in enumerate(tri):
                if abs(self.vertices[vertex, 0]) < 1e-13:
                    zz = self.vertices[vertex, 1]
                    axis_a = q if -d/2 < zz < d/2 else 0.
                    c = mu-axis_a
                    # Select the regular Frobenius branch.  P1 accurately
                    # resolves the integer 0/1 powers in the q=1 benchmark.
                    if abs(c) < .5-1e-12:
                        raise ValueError('Unclassified singular-axis branch')
                    forbidden.append(6*it+2*a+(1 if c > 0 else 0))
                    if abs(c) > .5+1e-12:
                        forbidden.append(6*it+2*a+(0 if c > 0 else 1))
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edge_owners.setdefault(tuple(sorted((int(tri[a]), int(tri[b])))),
                                       []).append((it, a, b))

        boundary_lookup = {tuple(sorted(map(int, item['vertices']))): item
                           for item in boundary}
        for key, owners in edge_owners.items():
            if len(owners) not in (1, 2):
                raise ValueError('Nonmanifold mesh edge')
            it, a, b = owners[0]
            p0, p1 = self.vertices[self.triangles[it, [a, b]]]
            tangent = p1-p0
            length = np.linalg.norm(tangent)
            normal = np.array([tangent[1], -tangent[0]])/length
            points = (1-ex[:, None])*p0+ex[:, None]*p1
            weights = ew*length*points[:, 0]
            basis = np.zeros((len(ex), 3))
            basis[:, a], basis[:, b] = 1-ex, ex
            an = normal[0]*SX+normal[1]*SZ
            pm, pp = (ID-an)/2, (ID+an)/2
            fm = np.einsum('p,pa,pb->ab', weights, basis, basis)
            block(it, it, np.kron(fm, pm))
            entry = dict(triangle=it, points=points, weights=weights, basis=basis,
                         normal=normal, pm=pm, pp=pp, vertices=key)
            if len(owners) == 2:
                jt, _, _ = owners[1]
                other_basis = np.zeros_like(basis)
                for a_local, global_vertex in enumerate(self.triangles[it]):
                    where = np.where(self.triangles[jt] == global_vertex)[0]
                    if len(where):
                        other_basis[:, where[0]] = basis[:, a_local]
                mixed = np.einsum('p,pa,pb->ab', weights, basis, other_basis)
                other_mass = np.einsum('p,pa,pb->ab', weights, other_basis, other_basis)
                block(it, jt, -np.kron(mixed, pm))
                block(jt, jt, np.kron(other_mass, pp))
                block(jt, it, -np.kron(mixed.T, pp))
                entry.update(other=jt, other_basis=other_basis)
                self.interior.append(entry)
            else:
                if key not in boundary_lookup:
                    raise ValueError('Unclassified boundary edge')
                supplied = boundary_lookup[key]
                if np.linalg.norm(normal-np.asarray(supplied['normal'])) > 1e-10:
                    raise ValueError('Mesh boundary outward normal is inconsistent')
                entry['kind'] = supplied['kind']
                self.boundary.append(entry)
        self.matrix_full = coo_matrix((values, (rows, cols)),
                                     shape=(self.ndof, self.ndof)).tocsc()
        self.active = np.setdiff1d(np.arange(self.ndof), np.unique(forbidden))
        self.matrix = self.matrix_full[self.active][:, self.active].tocsc()
        self.lu = splu(self.matrix)

    def input(self, profile):
        """Lift a callable profile(kind, points, normal) -> (npoints,2)."""
        rhs = np.zeros(self.ndof, complex)
        data, norm = [], 0.
        for edge in self.boundary:
            raw = np.asarray(profile(edge['kind'], edge['points'], edge['normal']), complex)
            if raw.shape != (len(edge['weights']), 2) or not np.all(np.isfinite(raw)):
                raise ValueError('Boundary profile must be a finite (npoints,2) array')
            incoming = raw@edge['pm'].T
            local = np.einsum('p,pa,ps->as', edge['weights'], edge['basis'], incoming)
            start = 6*edge['triangle']
            rhs[start:start+6] += local.reshape(6)
            norm += np.einsum('p,ps,ps->', edge['weights'], incoming.conj(), incoming).real
            data.append(incoming)
        return rhs, data, float(norm)

    def solve(self, rhs, feedback=None):
        """Solve fixed incoming data; optional finite-rank exterior feedback."""
        rhs = np.asarray(rhs, complex)
        reduced = rhs[self.active]
        result = self.lu.solve(reduced)
        if feedback is not None:
            G, gamma, C = map(np.asarray, feedback)
            if G.shape[0] != self.ndof or C.shape[1] != self.ndof:
                raise ValueError('Feedback must use full DG coefficient dimensions')
            lifted = self.lu.solve(G[self.active]*gamma[None, :])
            cr = C[:, self.active]
            result += lifted@np.linalg.solve(np.eye(len(gamma))-cr@lifted, cr@result)
        full = np.zeros((self.ndof,)+rhs.shape[1:], complex)
        full[self.active] = result
        return full

    def trace(self, solution, edge):
        return edge['basis']@solution[6*edge['triangle']:6*edge['triangle']+6].reshape(3, 2)

    def ledger(self, solution, incoming_data, rhs):
        """Discrete integration by parts, not an empirical fitted probability."""
        fluxes = {kind: 0. for kind in ('upper', 'lower', 'outer', 'axis')}
        input_norm, mismatch, actual_current, jump = 0., 0., 0., 0.
        for edge, incoming in zip(self.boundary, incoming_data):
            trace = self.trace(solution, edge)
            plus, minus = trace@edge['pp'].T, trace@edge['pm'].T
            norm = lambda v: float(np.einsum('p,ps,ps->', edge['weights'], v.conj(), v).real)
            fluxes[edge['kind']] += norm(plus)
            input_norm += norm(incoming)
            mismatch += norm(minus-incoming)
            actual_current += norm(plus)-norm(minus)
        for edge in self.interior:
            left = self.trace(solution, edge)
            j = edge['other']
            right = edge['other_basis']@solution[6*j:6*j+6].reshape(3, 2)
            diff = left-right
            jump += float(np.einsum('p,ps,ps->', edge['weights'], diff.conj(), diff).real)
        residual = self.matrix@solution[self.active]-rhs[self.active]
        algebraic = 2*np.vdot(solution[self.active], residual).real
        balance = sum(fluxes.values())+jump+mismatch-input_norm
        return dict(input_norm=input_norm, outgoing=fluxes, interior_jump_loss=jump,
                    incoming_mismatch_loss=mismatch, numerical_loss=jump+mismatch,
                    balance_residual=balance, algebraic_balance_term=float(algebraic),
                    identity_error=float(abs(balance-algebraic)),
                    physical_trace_net_current=actual_current,
                    linear_residual_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs[self.active]), 1e-300)))

    def projection(self, kind, profile, sign=1):
        """Unit-norm test trace and coefficient row; not a physical lead assertion."""
        row = np.zeros(self.ndof, complex)
        norm = 0.
        for edge in self.boundary:
            if edge['kind'] != kind:
                continue
            p = edge['pp'] if sign == 1 else edge['pm']
            value = np.asarray(profile(kind, edge['points'], edge['normal']))@p.T
            norm += np.einsum('p,ps,ps->', edge['weights'], value.conj(), value).real
            start = 6*edge['triangle']
            row[start:start+6] += np.einsum('p,pa,ps->as', edge['weights'], edge['basis'], value.conj()).reshape(6)
        if norm <= 0:
            raise ValueError('Projection profile has zero boundary norm')
        return row/math.sqrt(norm), float(norm)

    def l2_error(self, solution, exact):
        err, norm = 0., 0.
        for cell in self.cells:
            it = cell['triangle']
            numerical = cell['basis']@solution[6*it:6*it+6].reshape(3, 2)
            reference = exact(cell['points'])
            err += np.einsum('p,ps,ps->', cell['weights'], (numerical-reference).conj(), numerical-reference).real
            norm += np.einsum('p,ps,ps->', cell['weights'], reference.conj(), reference).real
        return float(math.sqrt(err/norm))

    def modal_port(self, kind, basis_function):
        """Facet-project and separately orthonormalize outgoing/incoming bases.

        basis_function(points) returns (plus,minus), each (2,npoints,nmodes).
        Cholesky orthogonalization keeps the FIRST column fixed except for a
        positive normalization, which preserves an explicitly first zero mode.
        Applying diagonal lead ratios afterwards is a convergent polygonal
        discretization, not exact spherical matching on each flat facet.
        """
        edges = [e for e in self.boundary if e['kind'] == kind]
        weights = np.concatenate([e['weights'] for e in edges])
        raw_plus, raw_minus = [], []
        for edge in edges:
            plus, minus = basis_function(edge['points'])
            plus, minus = np.moveaxis(plus, 0, 1), np.moveaxis(minus, 0, 1)
            raw_plus.append(np.einsum('ab,pbc->pac', edge['pp'], plus))
            raw_minus.append(np.einsum('ab,pbc->pac', edge['pm'], minus))

        def orthogonalize(raw):
            raw = np.concatenate(raw)
            gram = np.einsum('p,pai,paj->ij', weights, raw.conj(), raw)
            gram = (gram+gram.conj().T)/2
            factor = np.linalg.cholesky(gram)
            transform = np.linalg.solve(factor.conj().T, np.eye(len(factor)))
            basis = raw@transform
            after = np.einsum('p,pai,paj->ij', weights, basis.conj(), basis)
            return (basis, float(np.linalg.cond(gram)), float(np.max(abs(after-np.eye(len(after))))),
                    float(np.max(abs(gram-np.eye(len(gram))))))

        plus, plus_condition, plus_error, plus_deviation = orthogonalize(raw_plus)
        minus, minus_condition, minus_error, minus_deviation = orthogonalize(raw_minus)
        C = np.zeros((plus.shape[2], self.ndof), complex)
        G = np.zeros((self.ndof, minus.shape[2]), complex)
        start = 0
        for edge in edges:
            end = start+len(edge['weights'])
            dof = 6*edge['triangle']+np.arange(6)
            C[:, dof] += np.einsum('p,pa,psn->nas', edge['weights'], edge['basis'], plus[start:end].conj()).reshape(len(C), 6)
            G[dof] += np.einsum('p,pa,psn->asn', edge['weights'], edge['basis'], minus[start:end]).reshape(6, -1)
            start = end
        return dict(kind=kind, edges=edges, weights=weights, plus=plus, minus=minus,
                    G=G, C=C, gram_condition=max(plus_condition, minus_condition),
                    gram_error=max(plus_error, minus_error),
                    unnormalized_gram_deviation=max(plus_deviation, minus_deviation))

    def modal_incoming(self, ports, coefficients):
        """Expand all incoming modal coefficients into the ledger quadrature."""
        lookup = {}
        for port, coef in zip(ports, coefficients):
            values = port['minus']@coef
            start = 0
            for edge in port['edges']:
                end = start+len(edge['weights'])
                lookup[id(edge)] = values[start:end]
                start = end
        return [lookup.get(id(e), np.zeros((len(e['weights']), 2), complex)) for e in self.boundary]

    def cell_diagnostics(self, solution):
        values = solution.reshape(-1, 3, 2).mean(axis=1)
        centroids = self.vertices[self.triangles].mean(axis=1)
        density = centroids[:, 0]*np.sum(abs(values)**2, axis=1)
        currents = np.column_stack([centroids[:, 0]*np.einsum('pa,ab,pb->p', values.conj(), matrix, values).real
                                    for matrix in (SX, SZ)])
        return dict(vertices=self.vertices.tolist(), triangles=self.triangles.tolist(),
                    centroids=centroids.tolist(), rho_u_squared=density.tolist(),
                    weighted_meridional_current=currents.tolist(),
                    interpretation='Cell-centroid values of the actual computed field; not a trajectory or probability of escape.')


def mouth_profile(kind, energy=.4, d=1., outgoing=False):
    """Existing q=1 monopole reference spinor, projected on polygon normals.

    The lower outgoing reference is the physical negative-charge zero mode in
    the common gauge.  The lower incoming reference is auxiliary, not a zero
    mode of the negative mouth.  A finite collar model remains necessary.
    """
    center = d/2 if kind == 'upper' else -d/2
    sign = -1 if outgoing else 1

    def profile(edge_kind, points, normal):
        if edge_kind != kind:
            return np.zeros((len(points), 2), complex)
        rho, z = points.T
        reduced = collar_mode(1, 0, rho, z-center, energy, charge=sign).T
        return reduced/np.sqrt(rho[:, None])
    return profile


def _complex_json(a):
    a = np.asarray(a)
    return dict(real=a.real.tolist(), imag=a.imag.tolist())


def solve_case(h=.12, R=2., save_trace=False):
    from two_mouth_mesh import mesh
    started = time.perf_counter()
    model = ExteriorDG(mesh(b=.2, d=1., R=R, h=h))
    rhs_list, input_list = [], []
    for kind in ('upper', 'lower'):
        rhs, data, norm = model.input(mouth_profile(kind))
        rhs_list.append(rhs/math.sqrt(norm))
        input_list.append([v/math.sqrt(norm) for v in data])
    solutions = model.solve(np.column_stack(rhs_list))
    ledgers = [model.ledger(solutions[:, j], input_list[j], rhs_list[j]) for j in range(2)]
    projectors = np.array([model.projection(kind, mouth_profile(kind, outgoing=True))[0]
                           for kind in ('upper', 'lower')])
    coherent = projectors@solutions
    result = dict(h=h, R=R, triangles=len(model.triangles), active_dof=len(model.active),
                  ledgers=ledgers, mouth_coherent_projection=_complex_json(coherent),
                  seconds=time.perf_counter()-started)
    if save_trace:
        result['boundary_trace_map'] = []
        for edge in model.boundary:
            if edge['kind'] == 'axis':
                continue
            traces = np.stack([model.trace(solutions[:, j], edge)@edge['pp'].T for j in range(2)], axis=-1)
            result['boundary_trace_map'].append(dict(kind=edge['kind'], points=edge['points'].tolist(),
                        weights=edge['weights'].tolist(), normal=edge['normal'].tolist(),
                        outgoing=_complex_json(traces)))
    return result


def manufactured_case(h=.15):
    """Exact regular free cylindrical Bessel wave on the SAME two-hole domain."""
    from two_mouth_mesh import mesh
    energy, transverse = 1.3, .8
    kz = math.sqrt(energy**2-transverse**2)

    def exact(points):
        rho, z = points.T
        return np.exp(1j*kz*z[:, None])*np.column_stack((jv(0, transverse*rho),
                         1j*transverse/(energy+kz)*jv(1, transverse*rho)))

    model = ExteriorDG(mesh(b=.2, d=1., R=1.3, h=h), q=0, energy=energy)
    rhs, data, norm = model.input(lambda kind, points, normal: exact(points))
    solution = model.solve(rhs)
    return dict(h=h, triangles=len(model.triangles), error=model.l2_error(solution, exact),
                ledger=model.ledger(solution, data, rhs))


def free_outgoing_case(h=.12, nmax=10, R=2.):
    """Outgoing point Green spinor: independent full-free exterior map control."""
    from two_mouth_mesh import mesh
    from spherical_weyl_leads import outer_characteristics, outgoing_ratios
    energy, d = .4, 1.

    def exact(points):
        rho, z = points.T
        z = z-d/2
        r = np.hypot(rho, z)
        green = np.exp(1j*energy*r)/(4*np.pi*r)
        derivative = green*(1j*energy-1/r)
        return np.sqrt(2*np.pi)*np.column_stack((energy*green-1j*z/r*derivative,
                                                 -1j*rho/r*derivative))

    model = ExteriorDG(mesh(b=.2, d=d, R=R, h=h), q=0, energy=energy)
    outer = model.modal_port('outer', lambda points: outer_characteristics(.5,
                    np.arctan2(points[:, 0], points[:, 1]), R, nmax))
    gamma = outgoing_ratios(energy, R, .5, nmax)
    rhs, data, _ = model.input(lambda kind, points, normal:
                            np.zeros((len(points), 2), complex) if kind in ('outer','axis') else exact(points))
    solution = model.solve(rhs, (outer['G'], gamma, outer['C']))
    outgoing = outer['C']@solution
    incoming = gamma*outgoing
    feedback_data = model.modal_incoming([outer], [incoming])
    data = [a+b for a, b in zip(data, feedback_data)]
    total_rhs = rhs+outer['G']@incoming
    ledger = model.ledger(solution, data, total_rhs)
    net_outer = float(np.sum((1-abs(gamma)**2)*abs(outgoing)**2))
    raw_outer = ledger['outgoing']['outer']
    omitted = float(raw_outer-np.vdot(outgoing, outgoing).real)
    return dict(h=h, nmax=nmax, R=R, relative_l2_error=model.l2_error(solution, exact),
                outer_retained_net_flux=net_outer, outer_omitted_outgoing_norm=omitted,
                exact_full_outward_flux=energy**2/(2*np.pi), ledger=ledger,
                gram_error=outer['gram_error'], gram_condition=outer['gram_condition'])


def connected_lead_case(h=.08, R=2., nmax=8, energy=.4, save_field=False):
    """One physical upper zero-mode input, specified product/free lead closure.

    Nonzero mouth modes use a fixed-radius cylindrical monopole product lead;
    the exterior beyond R is replaced by free space.  Both are explicit model
    boundaries, not derived curved-mouth/exterior matching.  The actual pair
    potential is retained everywhere in the connected finite domain.
    """
    from two_mouth_mesh import mesh
    from spherical_weyl_leads import (outer_characteristics, outgoing_ratios,
                    mouth_characteristics, mouth_outgoing_ratios, zero_mode)
    started = time.perf_counter()
    b, d = .2, 1.
    model = ExteriorDG(mesh(b=b, d=d, R=R, h=h), energy=energy)

    def mouth_basis(port, points):
        theta = np.arctan2(points[:, 0], points[:, 1]-(d/2 if port == 'upper' else -d/2))
        plus, minus = mouth_characteristics(theta, b, port, nmax)
        zero = zero_mode(theta, b, port)[:, :, None]
        if port == 'upper':
            minus = np.concatenate((zero, minus), axis=2)
        else:
            plus = np.concatenate((zero, plus), axis=2)
        return plus, minus

    upper = model.modal_port('upper', lambda points: mouth_basis('upper', points))
    lower = model.modal_port('lower', lambda points: mouth_basis('lower', points))
    outer = model.modal_port('outer', lambda points: outer_characteristics(.5,
                          np.arctan2(points[:, 0], points[:, 1]), R, nmax))
    ports = [upper, lower, outer]
    beta = mouth_outgoing_ratios(energy, b, nmax)
    gamma = outgoing_ratios(energy, R, .5, nmax)
    G = np.column_stack((upper['G'][:, 1:], lower['G'], outer['G']))
    C = np.vstack((upper['C'], lower['C'][1:], outer['C']))
    ratios = np.concatenate((beta, beta, gamma))
    rhs = upper['G'][:, 0]
    solution = model.solve(rhs, (G, ratios, C))
    coefficients = [port['C']@solution for port in ports]
    reflected = [np.r_[1., beta*coefficients[0]], beta*coefficients[1][1:], gamma*coefficients[2]]
    data = model.modal_incoming(ports, reflected)
    total_rhs = sum((port['G']@inc for port, inc in zip(ports, reflected)), np.zeros(model.ndof, complex))
    ledger = model.ledger(solution, data, total_rhs)
    omitted = {port['kind']: float(ledger['outgoing'][port['kind']]-np.vdot(coef, coef).real)
               for port, coef in zip(ports, coefficients)}
    lower_zero = coefficients[1][0]
    lower_zero_flux = float(abs(lower_zero)**2)
    outer_flux = float(np.sum((1-abs(gamma)**2)*abs(coefficients[2])**2))
    mouth_nonzero_flux = float(np.sum((1-abs(beta)**2)*(abs(coefficients[0])**2+abs(coefficients[1][1:])**2)))
    reduced_balance = (lower_zero_flux+outer_flux+mouth_nonzero_flux+sum(omitted.values())
                       +ledger['numerical_loss']-1)
    result = dict(h=h, R=R, nmax=nmax, energy=energy, triangles=len(model.triangles),
                  active_dof=len(model.active), lower_zero_amplitude=_complex_json(lower_zero),
                  lower_zero_flux=lower_zero_flux, retained_free_outer_flux=outer_flux,
                  mouth_nonzero_net_flux=mouth_nonzero_flux,
                  omitted_outgoing_norm=omitted, ledger=ledger,
                  reduced_physical_plus_numerical_balance_residual=float(reduced_balance),
                  outgoing_coefficients={p['kind']:_complex_json(c) for p,c in zip(ports, coefficients)},
                  incoming_feedback_ratios=dict(mouth=_complex_json(beta), outer=_complex_json(gamma)),
                  modal_gram_error=max(p['gram_error'] for p in ports),
                  modal_gram_condition=max(p['gram_condition'] for p in ports),
                  unnormalized_gram_deviation={p['kind']:p['unnormalized_gram_deviation'] for p in ports},
                  seconds=time.perf_counter()-started)
    if save_field:
        result['field'] = model.cell_diagnostics(solution)
        result['boundary_trace_map'] = []
        for edge in model.boundary:
            if edge['kind'] == 'axis':
                continue
            trace = model.trace(solution, edge)
            result['boundary_trace_map'].append(dict(kind=edge['kind'], points=edge['points'].tolist(),
                weights=edge['weights'].tolist(), normal=edge['normal'].tolist(),
                outgoing=_complex_json(trace@edge['pp'].T), incoming=_complex_json(trace@edge['pm'].T)))
    return result


def checks():
    started = time.perf_counter()
    manufactured = [manufactured_case(h) for h in (.12, .08, .05)]
    free_green = [free_outgoing_case(h, 12) for h in (.12, .08, .05, .035)]
    mesh_refinement = [connected_lead_case(h, 2., 8) for h in (.12, .08, .05, .035)]
    angular_refinement = [mesh_refinement[2], connected_lead_case(.05, 2., 12)]
    radius_refinement = [angular_refinement[1], connected_lead_case(.05, 4., 12),
                         connected_lead_case(.05, 8., 12)]
    representative = connected_lead_case(.035, 4., 12, save_field=True)
    energy_samples = [connected_lead_case(.05, 4., 12, energy=e) for e in (.1, 1.)]
    energy_samples.insert(1, radius_refinement[1])
    connected = mesh_refinement+angular_refinement+radius_refinement+[representative]+energy_samples
    ledgers = [row['ledger'] for row in manufactured+free_green+connected]
    complex_value = lambda row: complex(row['lower_zero_amplitude']['real'], row['lower_zero_amplitude']['imag'])
    angular_difference = abs(complex_value(angular_refinement[1])-complex_value(angular_refinement[0]))
    radius_differences = [abs(complex_value(b)-complex_value(a))
                          for a,b in zip(radius_refinement[:-1], radius_refinement[1:])]
    random_identity = random_current_identity()
    tests = dict(
        exact_discrete_current_identity=max(r['identity_error'] for r in ledgers)<1e-10,
        sparse_linear_solve=max(r['linear_residual_relative'] for r in ledgers)<1e-10,
        manufactured_error_decreases=all(a['error']>b['error'] for a,b in zip(manufactured[:-1], manufactured[1:])),
        independent_free_green_error_decreases=all(a['relative_l2_error']>b['relative_l2_error'] for a,b in zip(free_green[:-1], free_green[1:])),
        independent_free_outgoing_flux_error_decreases=all(abs(a['outer_retained_net_flux']/a['exact_full_outward_flux']-1)>
                    abs(b['outer_retained_net_flux']/b['exact_full_outward_flux']-1) for a,b in zip(free_green[:-1],free_green[1:])),
        positive_numerical_losses=all(r['numerical_loss']>=0 for r in ledgers),
        connected_ledger_closes=max(abs(r['reduced_physical_plus_numerical_balance_residual']) for r in connected)<1e-10,
        numerical_loss_decreases_with_h=all(a['ledger']['numerical_loss']>b['ledger']['numerical_loss'] for a,b in zip(mesh_refinement[:-1],mesh_refinement[1:])),
        omitted_outgoing_norm_nonnegative=all(min(r['omitted_outgoing_norm'].values())>-1e-12 for r in connected),
        polygonal_modal_bases_orthonormal=max(r['modal_gram_error'] for r in connected)<1e-12,
        angular_refinement_small=angular_difference<2e-5,
        evanescent_mouth_modes_have_zero_net_flux=max(abs(r['mouth_nonzero_net_flux']) for r in connected)<1e-12,
        random_field_discrete_identity=random_identity['maximum_relative_error']<1e-12,
        connected_source_reaches_other_mouth=representative['lower_zero_flux']>0,
    )
    assert all(tests.values()), tests
    return dict(passed=True, tests=tests, parameters=dict(q=1, mu=.5, energy=.4, d=1., b=.2),
                representative=representative, manufactured=manufactured, free_green_controls=free_green,
                mesh_refinement=mesh_refinement, angular_refinement=angular_refinement,
                radius_refinement=radius_refinement, energy_samples=energy_samples,
                random_current_identity=random_identity,
                convergence=dict(angular_complex_amplitude_difference=angular_difference,
                    radius_complex_amplitude_differences=radius_differences,
                    representative_mesh_complex_amplitude_difference=abs(complex_value(representative)-complex_value(radius_refinement[1])),
                    outer_radius_convergence_certified=False,
                    interpretation='Finite differences are sensitivity diagnostics, not rigorous continuum error bounds. R dependence remains material.'),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                input_source_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest()
                       for name in ('two_mouth_mesh.py','spherical_weyl_leads.py','incoming_modes.py')},
                seconds=time.perf_counter()-started,
                equations=dict(unknown='u=f/sqrt(rho), inner product integral rho*u_dagger*u drho dz',
                   strong='sigma_x*u_rho+sigma_z*u_z+(sigma_x/2+i*sigma_y*(mu-Aphi))*u/rho-i*E*u=0',
                   face='strong correction P_minus(u_left-u_right); boundary P_minus(u-g)',
                   ledger='input = all_outgoing + interior_jump_loss + incoming_mismatch_loss',
                   effective_matrix='A_eff=A-G diag(gamma) C, solved by Woodbury over sparse LU',
                   closed_ledger='1=abs(T_lower_zero)^2+sum_outer[(1-abs(gamma)^2)*abs(c)^2]+omitted_outgoing+DG_loss',
                   axis='For q=1, mu=1/2: lower spin component vanishes on outside axis, upper on between-mouth axis.'),
                sources=[dict(url='https://arxiv.org/abs/1807.04726', use='Magnetic channel context; no full MMP geometry is solved here.'),
                         dict(url='https://arxiv.org/abs/1704.05481', use='Free spherical Weyl harmonics and first-order radial equations; independent lead checks are in the helper.')],
                scope=['Connected, flat, stationary single-Weyl-sector finite-domain diagnostic.',
                   'Actual pair potential is solved between two finite mouth cuts; q=1,d=1,b=.2 is distinct from the saved q=100000 physical-scale candidate.',
                   'One upper monopole zero-mode input and lower zero-mode output; higher mouth modes use an explicitly chosen semi-infinite fixed-radius product lead.',
                   'Outer modal map is outgoing free space after truncating the actual magnetic dipole field at R; radius convergence is not established.',
                   'Circle boundaries are polygonal: project to facet characteristics, jointly orthonormalize each side keeping zero first, then retain diagonal lead ratios. This is an extra geometric approximation that is refined.',
                   'The physical Green control uses analytic incoming data on both cuts and exact free outside; it is not the single-mouth production input.',
                   'Numerical loss is separate and must converge away before a physical scattering interpretation.',
                   'Omitted modal outgoing norm is an unresolved boundary remainder, not certified physical escape.',
                   'Raw complex traces and modal amplitudes retained; probabilities alone cannot determine a vacuum source.',
                   'Correct current balance and smooth controls do not prove absence of all DG spectral pollution or all channel errors.',
                   'No Casimir update, global exterior convergence, curved-mouth match, or time-travel conclusion.'])


def random_current_identity():
    """Integration-by-parts identity for arbitrary fields, independent of LU."""
    from two_mouth_mesh import mesh
    model = ExteriorDG(mesh(b=.2, d=1., R=1.3, h=.16))
    rng = np.random.default_rng(9081)
    errors = []
    zero_data = [np.zeros((len(e['weights']),2), complex) for e in model.boundary]
    for _ in range(3):
        solution = np.zeros(model.ndof, complex)
        solution[model.active] = rng.normal(size=len(model.active))+1j*rng.normal(size=len(model.active))
        ledger = model.ledger(solution, zero_data, np.zeros(model.ndof, complex))
        scale = sum(ledger['outgoing'].values())+ledger['numerical_loss']
        errors.append(ledger['identity_error']/scale)
    return dict(maximum_relative_error=max(errors), samples=len(errors), random_seed=9081)


if __name__ == '__main__':
    result = checks()
    project_path('exterior_boundary_map_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed', 'tests', 'seconds', 'source_sha256')}, indent=2))
