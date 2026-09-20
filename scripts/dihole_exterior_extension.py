"""Compact exterior-only trial extensions, reflected across the equatorial plane.

New supports are rho in (.005,.5), z in (.60,1.20) and (-1.20,-.60).
Their squared polynomial envelopes and first derivatives vanish at every
support boundary. They are disjoint from the two mouth cuts, so their bridge
extension is zero. The returned data contain the frozen 48 global columns
first and the new component-first columns afterward.

Only these two supports are integrated. The old-old block of any resulting
Gram matrix is already included in the old full form and MUST NOT be added
again; the caller keeps just the new-new and old-new blocks. No physical
boundary condition, spectral convergence, or source conclusion is supplied.
"""
import json
import operator

import numpy as np
from numpy.polynomial.legendre import Legendre

from complete_angular_modes import _sector
from dihole_curved_exterior import coefficient_values
from dihole_sewn_form import (RHO_BOUNDS, SX4, SY4, SZ4, matrix_action,
                              scalar_trials, spin_values)
from dihole_stretched_collar import gauss


POSITIVE_Z_BOUNDS = (.60, 1.20)
OLD_DEGREES = (3, 4)
OLD_COLUMNS = 4*int(np.prod(OLD_DEGREES))


def _integer(value, name):
    try:
        if isinstance(value, bool):
            raise TypeError
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError(name+' must be a positive integer.') from exc
    if value < 1:
        raise ValueError(name+' must be a positive integer.')
    return value


def _degrees(degrees):
    if len(degrees) != 2:
        raise ValueError('degrees must contain rho and z polynomial counts.')
    return tuple(_integer(value, 'degree count') for value in degrees)


def extension_metadata(degrees=(2, 2)):
    """New columns only: component-first, positive support then reflected partner.

    new_index starts at zero; column starts at 48 in extension_data. These
    are not the final augmented-system indices when compact bridge columns
    already follow the global columns in the caller's old matrix.
    """
    degrees = _degrees(degrees)
    records = []
    for component in range(4):
        for support in ('positive_z', 'negative_z'):
            for rho_degree in range(degrees[0]):
                for z_degree in range(degrees[1]):
                    index = len(records)
                    records.append(dict(new_index=index, column=OLD_COLUMNS+index,
                        component=component, chirality='R' if component < 2 else 'L',
                        support=support, rho_degree=rho_degree, z_degree=z_degree))
    return records


def _one_coordinate(coordinate, bounds, count):
    lo, hi = bounds
    x = 2*(coordinate-lo)/(hi-lo)-1
    inside = (x > -1) & (x < 1)
    bump = np.where(inside, (1-x*x)**2, 0.)
    bump_derivative = np.where(inside, -4*x*(1-x*x), 0.)
    values, derivatives = [], []
    for degree in range(count):
        poly = Legendre.basis(degree)
        values.append(bump*poly(x))
        derivatives.append((bump_derivative*poly(x)+bump*poly.deriv()(x))*2/(hi-lo))
    return np.array(values).T, np.array(derivatives).T


def extension_scalar_trials(points, degrees=(2, 2)):
    """Return scalar value, d_rho, d_z; negative support is F(rho,-z)."""
    degrees = _degrees(degrees)
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.all(np.isfinite(points)):
        raise ValueError('points must be a finite array of shape (n,2).')
    radial, dr = _one_coordinate(points[:, 0], RHO_BOUNDS, degrees[0])
    families = []
    for sign in (1., -1.):
        axial, dz = _one_coordinate(sign*points[:, 1], POSITIVE_Z_BOUNDS, degrees[1])
        families.append(tuple(np.einsum('pi,pj->pij', left, right).reshape(len(points), -1)
                              for left, right in ((radial, axial), (dr, axial), (radial, sign*dz))))
    return tuple(np.concatenate((plus, minus), axis=1) for plus, minus in zip(*families))


def extension_nodes(order=24):
    """Exact rectangular support panels; positive weights are drho*dz only."""
    order = _integer(order, 'order')
    if order < 2:
        raise ValueError('order must be at least two.')
    nodes, weights = gauss(order)
    positive, positive_weights = [], []
    radial_panels = ((RHO_BOUNDS[0], .05), (.05, .1), (.1, .2), (.2, RHO_BOUNDS[1]))
    axial_panels = ((POSITIVE_Z_BOUNDS[0], .85), (.85, POSITIVE_Z_BOUNDS[1]))
    for rlo, rhi in radial_panels:
        rho = (rhi+rlo)/2+(rhi-rlo)*nodes/2
        rw = (rhi-rlo)*weights/2
        for zlo, zhi in axial_panels:
            z = (zhi+zlo)/2+(zhi-zlo)*nodes/2
            zw = (zhi-zlo)*weights/2
            positive.append(np.column_stack((np.repeat(rho, order), np.tile(z, order))))
            positive_weights.append((rw[:, None]*zw[None, :]).ravel())
    positive, positive_weights = np.vstack(positive), np.concatenate(positive_weights)
    reflected = positive*np.array([1., -1.])
    return np.vstack((positive, reflected)), np.r_[positive_weights, positive_weights]


def extension_data(mu=.5, order=24, degrees=(2, 2)):
    """Nine-key 137 form data: old globals then new exterior-only columns.

    Default 80 columns = 48 old globals + 32 new trials. f=fp=0 everywhere.
    weight contains rho*drho*dz; assemble supplies the physical P or 1/P.
    """
    mu = _sector(mu)+.5
    degrees = _degrees(degrees)
    points, weights = extension_nodes(order)
    rho = points[:, 0, None, None]
    fields = coefficient_values(points)
    old = scalar_trials(points, OLD_DEGREES)
    new = extension_scalar_trials(points, degrees)
    value, dr, dz = [np.concatenate((spin_values(a), spin_values(b)), axis=2)
                     for a, b in zip(old, new)]
    radial = -1j*(matrix_action(SX4, dr+value/(2*rho))+matrix_action(SZ4, dz))
    angular = matrix_action(SY4, value)*(fields['NP']*(mu-fields['Aphysical']))[:, None, None]/rho
    return dict(value=value, radial=radial, angular=angular, P=fields['P'], W=fields['W'],
                f=np.zeros(len(points)), fp=np.zeros(len(points)), weight=weights*points[:, 0])


def self_check():
    """Small support, reflection, derivative and quadrature controls; no eigensolve."""
    degrees = (2, 2)
    samples = np.array([[.03, .64], [.17, .78], [.36, 1.08], [.48, 1.17]])
    at = extension_scalar_trials(samples, degrees)
    reflected = extension_scalar_trials(samples*np.array([1., -1.]), degrees)
    reflection_errors = [float(np.max(abs(a[:, :4]-sign*b[:, 4:])))
                         for a, b, sign in zip(at, reflected, (1., 1., -1.))]
    step = 1e-4
    derivative_errors = []
    for coordinate in (0, 1):
        direction = np.eye(2)[coordinate]*step
        shifts = [extension_scalar_trials(samples+t*direction, degrees)[0] for t in (-2, -1, 1, 2)]
        derivative = (shifts[0]-8*shifts[1]+8*shifts[2]-shifts[3])/(12*step)
        derivative_errors.append(float(np.max(abs(derivative-at[coordinate+1]))))
    boundary = []
    for sign in (1., -1.):
        for rho in RHO_BOUNDS:
            boundary.extend((rho, sign*z) for z in np.linspace(*POSITIVE_Z_BOUNDS, 7))
        for z in POSITIVE_Z_BOUNDS:
            boundary.extend((rho, sign*z) for rho in np.linspace(*RHO_BOUNDS, 7))
    boundary_error = max(float(np.max(abs(value))) for value in extension_scalar_trials(boundary))
    theta = np.linspace(0., np.pi, 33)
    cuts = np.vstack([np.column_stack((.05*np.sin(theta), center+.05*np.cos(theta)))
                      for center in (-.5, .5)])
    cut_error = max(float(np.max(abs(value))) for value in extension_scalar_trials(cuts))

    rows, flat_grams = [], []
    for order in (8, 12):
        points, weights = extension_nodes(order)
        old = scalar_trials(points, OLD_DEGREES)[0]
        new = extension_scalar_trials(points)[0]
        weighted = weights*points[:, 0]
        flat_grams.append((new.T@(weighted[:, None]*new), old.T@(weighted[:, None]*new)))
        rows.append(dict(order=order, nodes=len(points), area=float(np.sum(weights)),
            rho_moment=float(weights@points[:, 0]), z_moment=float(weights@points[:, 1]),
            minimum_mouth_distance=float(min(np.min(np.hypot(points[:, 0], points[:, 1]-center))
                                             for center in (-.5, .5))),
            minimum_weight=float(np.min(weights))))
    quadrature_error = max(float(np.max(abs(a-b))) for a, b in zip(*flat_grams))
    data = extension_data(.5, order=12)
    points, _ = extension_nodes(12)
    old_values = spin_values(scalar_trials(points, OLD_DEGREES)[0])
    old_error = float(np.max(abs(data['value'][:, :, :OLD_COLUMNS]-old_values)))
    actual_cross = data['value'][:, 0, :12].conj().T@(
        (data['weight']*data['P'])[:, None]*data['value'][:, 0, OLD_COLUMNS:OLD_COLUMNS+8])
    labels = extension_metadata()
    expected_area = 2*(RHO_BOUNDS[1]-RHO_BOUNDS[0])*(POSITIVE_Z_BOUNDS[1]-POSITIVE_Z_BOUNDS[0])
    expected_rho = (RHO_BOUNDS[1]**2-RHO_BOUNDS[0]**2)*(POSITIVE_Z_BOUNDS[1]-POSITIVE_Z_BOUNDS[0])
    tests = dict(
        support_boundary_value_and_action_zero=boundary_error == 0.,
        zero_mouth_traces_and_actions=cut_error == 0.,
        reflected_values_and_derivatives=max(reflection_errors) == 0.,
        independent_derivatives=max(derivative_errors) < 2e-8,
        exact_rectangular_moments=max(max(abs(row['area']-expected_area), abs(row['rho_moment']-expected_rho),
                                        abs(row['z_moment'])) for row in rows) < 3e-14,
        polynomial_overlap_quadrature=quadrature_error < 3e-14,
        cuts_clear_and_weights_positive=min(row['minimum_mouth_distance'] for row in rows) > .05
                                        and min(row['minimum_weight'] for row in rows) > 0.,
        old_global_values_unchanged=old_error == 0.,
        actual_old_new_overlap=float(np.linalg.norm(actual_cross)) > 1e-8,
        component_first_counts=len(labels) == 32 and all(sum(r['component'] == c for r in labels) == 8 for c in range(4)),
        no_bridge_variation=np.all(data['f'] == 0.) and np.all(data['fp'] == 0.))
    return dict(passed=bool(all(tests.values())), tests={key: bool(value) for key, value in tests.items()},
                boundary_error=boundary_error, mouth_cut_error=cut_error,
                reflection_errors=reflection_errors, derivative_errors=derivative_errors,
                quadrature=rows, flat_gram_change=quadrature_error,
                actual_old_new_norm=float(np.linalg.norm(actual_cross)),
                old_value_error=old_error, default_total_columns=data['value'].shape[-1],
                supports=dict(rho=list(RHO_BOUNDS), positive_z=list(POSITIVE_Z_BOUNDS), negative_z=[-1.2, -.6]),
                scope='Exterior-only conforming trial extension. Keep new/cross Gram blocks only; no spectral/source convergence claim.')


if __name__ == '__main__':
    report = self_check()
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report['passed']:
        raise SystemExit('Exterior extension controls failed.')
