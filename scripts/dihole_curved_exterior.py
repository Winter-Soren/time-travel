"""Finite-domain Weyl response in the actual balanced-dihole coefficients.

For the canonical flat-measure spinor f, put u=sqrt(N/(P*rho))*f.
The stationary equation is

  [SX*(d_rho+1/(2*rho)) + SZ*d_z
   + i*SY*N*P*(mu-Aphysical)/rho - i*E*P/N] u = 0.

Thus the frozen ExteriorDG principal part, axis domain and boundary current
weights rho*ds remain valid.  Only angular volume terms and the imaginary-
frequency mass change.  The physical volume norm is int rho*(P/N)*|u|^2.
The half-disk and mouth cuts are polygonal approximations.  Zero incoming
outer data are a FINITE boundary condition, not transparent Melvin matching.
No throat state, return probability, or quantum stress is supplied here.

Geometry: Emparan, https://arxiv.org/abs/hep-th/9906160, equations 8--16,
using the existing balanced_dihole implementation and linked Q=1/5,D=1
parameters.  Their binary64 evaluation is numerical, not an interval bound.
"""
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.sparse import block_diag
from scipy.sparse.linalg import splu

from balanced_dihole import symbolic_fields
from dihole_stretched_collar import exact_target_parameters
from exterior_boundary_map import (ExteriorDG, ID, SY, _quadrature,
                                   mouth_profile, pair_potential)
from two_mouth_mesh import mesh


@lru_cache(maxsize=1)
def _field_function():
    return symbolic_fields()[0]


def coefficient_values(points, q=1):
    """Actual N,P,Aphysical,W at cylindrical points; horizon foci excluded."""
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.all(np.isfinite(points)):
        raise ValueError('points must be a finite array of shape (n,2).')
    if not math.isfinite(float(q)) or q <= 0 or int(q) != q:
        raise ValueError('q must be a positive integer.')
    p = exact_target_parameters()
    rho, z = points.T
    if np.any(rho < 0):
        raise ValueError('Require rho>=0.')
    upper, lower = np.hypot(rho, z-p['k']), np.hypot(rho, z+p['k'])
    if np.any((upper == 0) | (lower == 0)):
        raise ValueError('The horizon foci are not in this coefficient domain.')
    total = upper+lower
    r = p['M']+total/2
    # Algebraically (lower-upper)/(2*k), avoiding subtraction of distances.
    cosine = 2*z/total
    raw = _field_function()(r, cosine, p['M'], p['a'], p['B'])
    lapse, meridian, ahat = [np.broadcast_to(np.asarray(v, float), rho.shape)
                             for v in raw[:3]]
    potential = q*p['physical_A_factor']*ahat
    weight = meridian/lapse
    product = lapse*meridian
    if (not all(np.all(np.isfinite(v)) for v in
                (lapse, meridian, potential, weight, product))
            or np.any(lapse <= 0) or np.any(meridian <= 0)):
        raise ValueError('Nonpositive or nonfinite exterior metric coefficients.')
    return dict(N=lapse, P=meridian, Aphysical=potential, W=weight,
                NP=product, prolate_r=r, prolate_cosine=cosine)


def graded_mesh(Rtarget=8., baseR=1.6, h=.028, b=.05, d=1., r0=.65):
    """Radially stretch a base mesh outside r0, preserving the two cut meshes.

    F(r)=r+(Rtarget-baseR)*((r-r0)/(baseR-r0))**2 for r>r0.
    The triangle connectivity is unchanged; an inverted triangle is an error,
    never repaired by exchanging vertices.  This is a larger finite domain
    with coarse distant elements, not a transparent boundary construction.
    """
    if (not np.all(np.isfinite([Rtarget, baseR, h, b, d, r0]))
            or not (Rtarget >= baseR > r0 > d/2+b)):
        raise ValueError('Require finite Rtarget>=baseR>r0>d/2+b.')
    original = mesh(b=b, d=d, R=baseR, h=h)
    old = original['vertices']
    vertices = old.copy()
    radius = np.linalg.norm(old, axis=1)
    moved = radius > r0
    delta = Rtarget-baseR
    new_radius = radius[moved]+delta*((radius[moved]-r0)/(baseR-r0))**2
    vertices[moved] *= (new_radius/radius[moved])[:, None]
    triangles = original['triangles'].copy()
    corners = vertices[triangles]
    first, second = corners[:, 1]-corners[:, 0], corners[:, 2]-corners[:, 0]
    areas = (first[:, 0]*second[:, 1]-first[:, 1]*second[:, 0])/2
    if not np.all(np.isfinite(areas)) or np.any(areas <= 0):
        bad = np.flatnonzero(~np.isfinite(areas) | (areas <= 0))
        raise ValueError(f'Radial grading inverted/degenerate triangles {bad[:8].tolist()}; '
                         f'minimum signed area {np.min(areas):.17g}.')
    boundary = []
    mouth_vertices, outer_vertices = [], []
    normal_closure = np.zeros(2)
    divergence_area = 0.
    minimum_normal_margin = math.inf
    for source in original['boundary_edges']:
        edge = dict(source)
        ids = triangles[edge['triangle'], edge['local_edge']]
        ends = vertices[ids]
        tangent = ends[1]-ends[0]
        length = float(np.linalg.norm(tangent))
        normal = np.array([tangent[1], -tangent[0]])/length
        edge.update(vertices=ids.tolist(), normal=normal.tolist())
        boundary.append(edge)
        normal_closure += length*normal
        divergence_area += length*np.dot(ends.mean(axis=0), normal)/2
        kind = edge['kind']
        if kind in ('upper', 'lower'):
            mouth_vertices.extend(ids.tolist())
            center = np.array([0., (1 if kind == 'upper' else -1)*d/2])
            margin = -float(np.dot(ends.mean(axis=0)-center, normal))
        elif kind == 'outer':
            outer_vertices.extend(ids.tolist())
            margin = float(np.dot(ends.mean(axis=0), normal))
        elif kind == 'axis':
            margin = -float(normal[0])
        else:
            raise ValueError(f'Unexpected boundary kind {kind!r}.')
        minimum_normal_margin = min(minimum_normal_margin, margin)
    mouth_vertices = np.unique(mouth_vertices)
    if not np.array_equal(vertices[mouth_vertices], old[mouth_vertices]):
        raise AssertionError('The requested radial grading moved a mouth-cut vertex.')
    if minimum_normal_margin <= 0:
        raise ValueError('A recomputed boundary normal points into the exterior domain.')
    edge_lengths = np.stack([np.linalg.norm(corners[:, j]-corners[:, i], axis=1)
                             for i, j in ((0, 1), (1, 2), (2, 0))], axis=1)
    far = np.any(radius[triangles] > r0, axis=1)
    exact_area = math.pi*Rtarget**2/2-math.pi*b*b
    exact_weighted_volume = 2*Rtarget**3/3-4*b**3/3
    area = float(np.sum(areas))
    weighted_volume = float(np.sum(areas*corners[:, :, 0].mean(axis=1)))
    outer_error = float(np.max(abs(np.linalg.norm(vertices[np.unique(outer_vertices)], axis=1)-Rtarget)))
    geometry = dict(elements=len(triangles), vertices=len(vertices),
                    minimum_area=float(np.min(areas)), maximum_area=float(np.max(areas)),
                    area=area, exact_area=exact_area,
                    area_relative_error=abs(area/exact_area-1),
                    weighted_volume=weighted_volume, exact_weighted_volume=exact_weighted_volume,
                    weighted_volume_relative_error=abs(weighted_volume/exact_weighted_volume-1),
                    maximum_edge_length=float(np.max(edge_lengths)),
                    maximum_far_triangle_diameter=float(np.max(edge_lengths[far])),
                    outer_vertex_radius_error=outer_error,
                    boundary_normal_closure=float(np.linalg.norm(normal_closure)),
                    boundary_divergence_area_error=abs(divergence_area-area),
                    minimum_boundary_normal_margin=minimum_normal_margin)
    grading = dict(base_R=float(baseR), target_R=float(Rtarget), fixed_radius=float(r0),
                   base_h=float(h), radial_map='r+(Rtarget-baseR)*((r-r0)/(baseR-r0))^2 for r>r0',
                   maximum_radial_stretch=1+2*delta/(baseR-r0),
                   unchanged_inner_vertices=bool(np.array_equal(vertices[~moved], old[~moved])),
                   unchanged_mouth_vertices=True, unchanged_connectivity=True,
                   moved_vertices=int(np.count_nonzero(moved)),
                   base_geometry=original['geometry'],
                   scope='Larger distant elements limit angular/axial resolution; finite outer data remain approximate.')
    return dict(vertices=vertices, triangles=triangles, boundary_edges=boundary,
                parameters=dict(b=b, d=d, R=Rtarget, h=h), geometry=geometry, grading=grading)


def graded_mesh_checks():
    """Geometry-only R=8,16 controls; no coefficient evaluation or DG solve."""
    rows = []
    for target in (8., 16.):
        data = graded_mesh(Rtarget=target)
        v, t, g = data['vertices'], data['triangles'], data['geometry']
        edges = {tuple(sorted((int(tri[i]), int(tri[j])))) for tri in t
                 for i, j in ((0, 1), (1, 2), (2, 0))}
        tests = dict(positive_orientation=g['minimum_area'] > 0,
                     fixed_cuts=data['grading']['unchanged_mouth_vertices'],
                     fixed_near_vertices=data['grading']['unchanged_inner_vertices'],
                     connected_topology=len(np.unique(t))-len(edges)+len(t) == 1,
                     outward_normals=g['minimum_boundary_normal_margin'] > 0,
                     boundary_normal_closure=g['boundary_normal_closure'] < 2e-12*target,
                     area_divergence=g['boundary_divergence_area_error'] < 2e-12*target**2,
                     target_outer_vertices=g['outer_vertex_radius_error'] < 2e-12*target,
                     polygon_geometry=g['area_relative_error'] < .01
                                      and g['weighted_volume_relative_error'] < .01)
        rows.append(dict(R=target, tests={key: bool(value) for key, value in tests.items()},
                         passed=bool(all(tests.values())), geometry=g, grading=data['grading']))
    return dict(passed=all(row['passed'] for row in rows), cases=rows)


def curved_model(kappa=.1, h=.04, R=1.1, b=.05, q=1, mu=.5,
                 quadrature_order=4, triangulation=None):
    """Return (ExteriorDG instance, full rho*W mass CSC, diagnostics).

    The original order-four angular block is subtracted at its ORIGINAL
    quadrature before adding the actual curved block at quadrature_order.
    A higher requested order therefore does not leave a flat-potential
    quadrature defect.  Cell points/basis/weights/mass subsequently describe
    physical-norm quadrature; original principal data are retained separately.
    Every boundary and interior-face weight/projector is unchanged.  Optional
    triangulation is used directly, with its b,R,h metadata authoritative;
    it must have separation D=1.  Fields are evaluated on the new geometry,
    never interpolated from a base mesh.
    """
    if not math.isfinite(float(kappa)) or kappa < 0:
        raise ValueError('kappa must be finite and nonnegative.')
    if isinstance(quadrature_order, bool) or int(quadrature_order) != quadrature_order or quadrature_order < 2:
        raise ValueError('quadrature_order must be an integer >=2.')
    started = time.perf_counter()
    p = exact_target_parameters()
    if triangulation is None:
        triangulation = mesh(b=b, d=p['D'], R=R, h=h)
    else:
        mp = triangulation['parameters']
        if mp['d'] != p['D']:
            raise ValueError('Supplied mesh separation must equal the fixed balanced D=1.')
        b, R, h = (mp[key] for key in ('b', 'R', 'h'))
    model = ExteriorDG(triangulation, q=q, mu=mu, energy=0., d=p['D'])
    bary, qw, _, _ = _quadrature(int(quadrature_order))
    old_bary, old_qw, _, _ = _quadrature(4)
    xy = model.vertices[model.triangles]
    jac = np.linalg.det(np.stack((xy[:, 1]-xy[:, 0], xy[:, 2]-xy[:, 0]), axis=1))
    points = np.einsum('pa,tad->tpd', bary, xy)
    values = coefficient_values(points.reshape(-1, 2), q)
    fields = {key: value.reshape(points.shape[:2]) for key, value in values.items()}
    angular_blocks, mass_blocks = [], []
    for it, cell in enumerate(model.cells):
        old_points = old_bary@xy[it]
        old_angular = np.einsum('p,pa,pb->ab',
            old_qw*jac[it]*(mu-pair_potential(*old_points.T, q=q, d=p['D'])),
            old_bary, old_bary)
        plain_weights = qw*jac[it]
        actual_angular = np.einsum('p,pa,pb->ab',
            plain_weights*fields['NP'][it]*(mu-fields['Aphysical'][it]), bary, bary)
        weighted = plain_weights*points[it, :, 0]*fields['W'][it]
        local_mass = np.einsum('p,pa,pb->ab', weighted, bary, bary)
        angular_blocks.append(np.kron(actual_angular-old_angular, 1j*SY))
        mass_blocks.append(np.kron(local_mass, ID))
        cell['principal_quadrature'] = {key: cell[key] for key in
                                        ('points', 'weights', 'basis', 'mass')}
        cell.update(points=points[it], basis=bary, weights=weighted, mass=local_mass,
                    lapse=fields['N'][it], optical_weight=fields['W'][it])
    angular_correction = block_diag(angular_blocks, format='csc')
    mass = block_diag(mass_blocks, format='csc')
    model.matrix_full = (model.matrix_full+angular_correction+kappa*mass).tocsc()
    model.matrix = model.matrix_full[model.active][:, model.active].tocsc()
    model.lu = splu(model.matrix)
    model.energy, model.kappa = 1j*kappa, float(kappa)
    model.curved_weighted_mass = mass
    model.curved_parameters = p.copy()
    anti = angular_correction+angular_correction.conj().T
    mass_defect = mass-mass.conj().T
    diagnostics = dict(
        parameters={**p, 'kappa': float(kappa), 'h': float(h), 'R': float(R),
                    'b': float(b), 'q': int(q), 'mu': float(mu),
                    'quadrature_order': int(quadrature_order)},
        mesh=triangulation['geometry'], mesh_grading=triangulation.get('grading'),
        active_dof=int(len(model.active)),
        coefficient_ranges={key: [float(np.min(fields[key])), float(np.max(fields[key]))]
                            for key in ('N', 'P', 'NP', 'Aphysical', 'W')},
        angular_antihermitian_defect=float(np.max(abs(anti.data), initial=0.)),
        mass_hermitian_defect=float(np.max(abs(mass_defect.data), initial=0.)),
        minimum_cell_mass_eigenvalue=float(min(np.linalg.eigvalsh(c['mass'])[0]
                                               for c in model.cells)),
        norm='Integral rho*(P/N)*|u|^2 d_rho d_z; fixed-mu azimuthal normalization.',
        current='Boundary and interior-face weights remain rho*ds.',
        boundary='Finite polygonal cuts; outer incoming field is specified, not transparent.',
        coefficient_scope='Binary64 evaluation of linked exact Q=1/5,D=1 balanced constants.',
        build_seconds=time.perf_counter()-started)
    return model, mass, diagnostics


def attenuation_ledger(model, mass, kappa, solution, incoming_data, rhs):
    """DG current identity including 2*kappa times the physical bulk norm."""
    if float(kappa) != model.kappa:
        raise ValueError('Ledger kappa must match the assembled curved model.')
    solution = np.asarray(solution, complex)
    if solution.shape != (model.ndof,) or len(incoming_data) != len(model.boundary):
        raise ValueError('Supply one full solution and data on every boundary edge.')
    base = model.ledger(solution, incoming_data, rhs)
    bulk = float(np.vdot(solution, mass@solution).real)
    quadrature = 0.
    for cell in model.cells:
        start = 6*cell['triangle']
        value = cell['basis']@solution[start:start+6].reshape(3, 2)
        quadrature += float(np.einsum('p,ps,ps->', cell['weights'],
                                      value.conj(), value).real)
    attenuation = 2*kappa*bulk
    undamped = base.pop('balance_residual')
    base.pop('identity_error')
    base.update(bulk_weighted_norm=bulk, bulk_norm_quadrature=quadrature,
                bulk_norm_quadrature_difference=abs(bulk-quadrature),
                bulk_attenuation=attenuation, undamped_balance_residual=undamped,
                balance_with_attenuation_residual=undamped+attenuation,
                attenuation_identity_error=abs(undamped+attenuation-base['algebraic_balance_term']))
    return base


def self_check():
    """One small mouth-driven solve; no lead or continuum convergence claim."""
    started = time.perf_counter()
    model, mass, diagnostics = curved_model()
    # This helper's legacy collar profile requires positive real phase energy;
    # it is only prescribed boundary data, independent of the resolvent kappa.
    rhs, incoming, norm = model.input(mouth_profile('upper', energy=.4))
    if not norm > 0:
        raise AssertionError('The declared upper-mouth test datum has zero norm.')
    rhs /= math.sqrt(norm)
    incoming = [value/math.sqrt(norm) for value in incoming]
    solution = model.solve(rhs)
    ledger = attenuation_ledger(model, mass, model.kappa, solution, incoming, rhs)
    # Off-shell trial checks the full integration-by-parts identity, not just
    # cancellation of a solved linear system's residual.
    rng = np.random.default_rng(131)
    trial = np.zeros(model.ndof, complex)
    trial[model.active] = rng.normal(size=len(model.active))+1j*rng.normal(size=len(model.active))
    trial /= np.linalg.norm(trial)
    off_shell = attenuation_ledger(model, mass, model.kappa, trial, incoming, rhs)
    axis = coefficient_values(np.array([[0., -.9], [0., 0.], [0., .9]]))
    tests = dict(
        positive_lapse_and_weight=diagnostics['coefficient_ranges']['N'][0] > 0
                                  and diagnostics['coefficient_ranges']['W'][0] > 0,
        positive_physical_mass=diagnostics['minimum_cell_mass_eigenvalue'] > 0,
        angular_antihermitian=diagnostics['angular_antihermitian_defect'] < 2e-14,
        hermitian_physical_mass=diagnostics['mass_hermitian_defect'] < 2e-14,
        unchanged_regular_axis=np.max(abs(axis['NP']-1)) < 2e-13
                               and np.max(abs(axis['Aphysical']-[0., 1., 0.])) < 2e-13,
        linear_residual=ledger['linear_residual_relative'] < 2e-11,
        physical_norm_quadrature=ledger['bulk_norm_quadrature_difference'] < 2e-11,
        solved_attenuation_identity=ledger['attenuation_identity_error'] < 2e-11,
        off_shell_attenuation_identity=off_shell['attenuation_identity_error'] < 2e-11)
    result = dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                  diagnostics=diagnostics, ledger=ledger,
                  off_shell_identity_error=off_shell['attenuation_identity_error'],
                  input_scope='Normalized isolated-monopole reference profile (phase energy .4), not a matched throat state.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  seconds=time.perf_counter()-started)
    return result


if __name__ == '__main__':
    report = self_check()
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report['passed']:
        raise SystemExit('Curved exterior helper check failed.')
