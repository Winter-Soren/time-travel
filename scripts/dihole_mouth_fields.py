"""Regular nonzero-trace fields on the frozen q=1, mu=1/2 sewn geometry.

The first four columns carry the monopole zero profile independently at each
mouth and in each Weyl block. Eight further columns carry the first angular
pair. Exterior fields are chi(r) U_side(theta) p(theta)/r, supported on
.05 <= r <= .1, with the existing right-cut phase. Their bridge extensions
are independent linear endpoint lifts. They are global H1 sections in the
declared spin/gauge bundle; normal derivatives need not match at the cuts.
This is a finite trial family, not new geometry, a boundary condition on the
physical field, or a source calculation. No azimuthal degeneracy is inserted.

Conventions: docs/GLOBAL_SEWN_FORM.md and scripts/dihole_azimuthal_bridge.py.
The actual magnetic/metric coefficients remain those of Emparan's balanced
dihole, https://arxiv.org/abs/hep-th/9906160, equations 8--16. The monopole
functions are a regular coordinate basis, not its exact angular eigenmodes.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from dihole_curved_exterior import coefficient_values
from dihole_enriched_bridge import _geometry
from dihole_radial_stream import bridge_chunks as old_bridge_chunks, RADIAL_PANELS
from dihole_sewn_form import (SX4, SY4, SZ4, RHO_BOUNDS, scalar_trials,
                              spin_values, matrix_action)
from dihole_stretched_collar import gauss, smooth_step, radius_profile, field_values
from project_paths import project_path


RC, ROUT, CENTER = .05, .1, .5
SY = np.array([[0., -1j], [1j, 0.]])
SIGMA = np.array([1., -1., -1., 1.])
FIELD_KEYS = ('value', 'radial', 'angular')


def column_metadata(include_pair=True):
    """New-column rows; indices0..3 are always the four zero-profile fields."""
    if not isinstance(include_pair, (bool, np.bool_)):
        raise ValueError('include_pair must be Boolean.')
    rows = []
    for profile in (('zero', 'upper1', 'lower1') if include_pair else ('zero',)):
        for side in ('left', 'right'):
            for chirality in ('R', 'L'):
                rows.append(dict(index=len(rows), family='mouth', profile=profile,
                                 side=side, chirality=chirality, mu=.5,
                                 l=0 if profile == 'zero' else 1))
    return rows


def profiles(theta):
    """Value, theta derivative and reference angular action: (points,2,3)."""
    theta = np.asarray(theta, float)
    if theta.ndim != 1 or not np.isfinite(theta).all() or np.any((theta < 0) | (theta > np.pi)):
        raise ValueError('theta must be a finite vector in [0,pi].')
    c, s = np.cos(theta), np.sin(theta)
    value = np.zeros((len(theta), 2, 3), complex)
    derivative, action = np.zeros_like(value), np.zeros_like(value)
    value[:, 0, 0] = 1/np.sqrt(2.)
    value[:, 0, 1] = np.sqrt(1.5)*c
    value[:, 1, 2] = .5j*np.sqrt(3.)*s
    derivative[:, 0, 1] = -np.sqrt(1.5)*s
    derivative[:, 1, 2] = .5j*np.sqrt(3.)*c
    action[:, 1, 1] = 1j*np.sqrt(1.5)*s
    action[:, 0, 2] = np.sqrt(3.)*c
    return value, derivative, action


def rotation(theta, side):
    """Exterior spinor = U_side times bridge-frame spinor, at the same angle."""
    if side not in ('left', 'right'):
        raise ValueError('side must be left or right.')
    sign = 1 if side == 'right' else -1
    c, s = np.cos(np.asarray(theta)/2), np.sin(np.asarray(theta)/2)
    return np.array([[c, -sign*s], [sign*s, c]], complex).transpose(2, 0, 1)


def taper(radius):
    """C-infinity compact taper and first radial derivative, flat at both ends."""
    step, derivative = smooth_step((np.asarray(radius, float)-RC)/(ROUT-RC))
    return 1-step, -derivative/(ROUT-RC)


def exterior_values(points, phase=0., include_pair=True):
    """New values and exact rho/z derivatives in the exterior flux frame.

    Only the actual exterior of both cuts is used by exterior_data. This
    analytic evaluation also accepts boundary points and the physical axis.
    """
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all() or np.any(points[:, 0] < 0):
        raise ValueError('Require finite (rho,z) points with rho>=0.')
    if not np.isfinite(phase):
        raise ValueError('phase must be finite.')
    metadata = column_metadata(include_pair)
    result = [np.zeros((len(points), 4, len(metadata)), complex) for _ in range(3)]
    for side in ('left', 'right'):
        sign = 1 if side == 'right' else -1
        rho, dz = points[:, 0], points[:, 1]-sign*CENTER
        radius = np.hypot(rho, dz)
        active = (radius <= ROUT) & (radius >= RC*(1-2e-13))
        at = np.flatnonzero(active)
        if not len(at):
            continue
        r = radius[at]
        theta = np.arctan2(rho[at], sign*dz[at])
        p, pt, _ = profiles(theta)
        U = rotation(theta, side)
        g = np.einsum('pab,pbj->paj', U, p)
        gt = np.einsum('pab,pbj->paj', U, pt-.5j*sign*np.einsum('ab,pbj->paj', SY, p))
        chi, chip = taper(r)
        h, hp = chi/r, chip/r-chi/r**2
        st, ct = np.sin(theta), np.cos(theta)
        value = h[:, None, None]*g
        dr = (st*hp)[:, None, None]*g+(ct*h/r)[:, None, None]*gt
        dzv = sign*((ct*hp)[:, None, None]*g-(st*h/r)[:, None, None]*gt)
        phase_factor = np.exp(1j*phase) if side == 'right' else 1.
        for row in metadata:
            if row['side'] != side:
                continue
            mode = ('zero', 'upper1', 'lower1').index(row['profile'])
            components = slice(0, 2) if row['chirality'] == 'R' else slice(2, 4)
            for target, raw in zip(result, (value, dr, dzv)):
                target[at, components, row['index']] = phase_factor*raw[:, :, mode]
    return tuple(result)


def bridge_values(s, theta, include_pair=True):
    """New bridge fields at one radial slice; physical angle spans both poles."""
    if not np.isfinite(s) or not -.1 <= s <= .1:
        raise ValueError('Bridge coordinate must lie in [-.1,.1].')
    p, _, reference = profiles(theta)
    if np.any((np.asarray(theta) <= 0) | (np.asarray(theta) >= np.pi)):
        raise ValueError('Differential data use interior polar nodes; endpoint values follow profiles().')
    function, parameters = _geometry()
    radius = float(radius_profile(min(.05, max(0., abs(s)-.05)))[0])
    fields = field_values(function, parameters, radius, np.cos(theta))
    correction = (fields['lapse']*fields['meridian_factor']*(.5-fields['unit_A_phi'])
                  -.5*np.cos(theta))/np.sin(theta)
    angular = (reference+correction[:, None, None]*np.einsum('ab,pbj->paj', SY, p))/radius
    metadata = column_metadata(include_pair)
    result = [np.zeros((len(theta), 4, len(metadata)), complex) for _ in range(3)]
    t = (s+.1)/.2
    for row in metadata:
        mode = ('zero', 'upper1', 'lower1').index(row['profile'])
        lift, slope = (t, 5.) if row['side'] == 'right' else (1-t, -5.)
        components = slice(0, 2) if row['chirality'] == 'R' else slice(2, 4)
        chiral = 1 if row['chirality'] == 'R' else -1
        j = row['index']
        result[0][:, components, j] = lift*p[:, :, mode]
        result[1][:, components, j] = -1j*chiral*np.array([1., -1.])*slope*p[:, :, mode]
        result[2][:, components, j] = chiral*lift*angular[:, :, mode]
    return dict(zip(FIELD_KEYS, result))


def bridge_chunks(phase=0., order=64, include_pair=True, slices_per_chunk=4):
    """Old312 plus new4/12 on exactly the frozen angular5/radial12 grid."""
    column_metadata(include_pair)
    x = gauss(order)[0]
    cmax = np.sqrt(1-(RHO_BOUNDS[0]/RC)**2)
    cosine = np.concatenate([(hi+lo)/2+(hi-lo)*x/2 for lo, hi in
                             ((-1., -cmax), (-cmax, cmax), (cmax, 1.))])
    theta = np.arccos(cosine)
    ss = np.concatenate([(hi+lo)/2+(hi-lo)*x/2 for lo, hi in RADIAL_PANELS])
    consumed = 0
    for old in old_bridge_chunks(.5, phase, 5, 12, order, (3, 4), slices_per_chunk):
        count = len(old['P'])//len(theta)
        new = [bridge_values(float(s), theta, include_pair) for s in ss[consumed:consumed+count]]
        result = dict(old)
        for key in FIELD_KEYS:
            result[key] = np.concatenate((old[key], np.concatenate([v[key] for v in new])), axis=2)
        consumed += count
        yield result
    if consumed != len(ss):
        raise ArithmeticError('Frozen radial streaming order changed.')


def exterior_data(phase=0., order=24, include_pair=True):
    """Return (data,coordinates): oldglobal48 + new4/12 on both round collars.

    Positive weights already equal rho*d_rho*d_z = r^2*dr*dcos(theta).
    At EACH radius, angular panels split exactly at rho=.005, the old
    global-field support edge. Existing exterior32 columns vanish here.
    """
    column_metadata(include_pair)
    x, w = gauss(order)
    rr, wr = (ROUT+RC)/2+(ROUT-RC)*x/2, (ROUT-RC)*w/2
    rows = []
    for side in ('left', 'right'):
        sign = 1 if side == 'right' else -1
        for r, weight in zip(rr, wr):
            cmax = np.sqrt(1-(RHO_BOUNDS[0]/r)**2)
            for panel, (lo, hi) in enumerate(((-1., -cmax), (-cmax, cmax), (cmax, 1.))):
                c, wc = (hi+lo)/2+(hi-lo)*x/2, (hi-lo)*w/2
                theta = np.arccos(c)
                rows.append((r*np.sin(theta), sign*(CENTER+r*c), weight*wc*r*r,
                             np.full(order, r), theta, np.full(order, sign), np.full(order, panel)))
    rho, z, weights, radius, theta, signs, panels = [np.concatenate(v) for v in zip(*rows)]
    points = np.column_stack((rho, z))
    fields = coefficient_values(points)
    old = [spin_values(v) for v in scalar_trials(points, (3, 4))]
    new = exterior_values(points, phase, include_pair)
    value, dr, dz = [np.concatenate((a, b), axis=2) for a, b in zip(old, new)]
    radial = -1j*(matrix_action(SX4, dr+value/(2*rho[:, None, None]))+matrix_action(SZ4, dz))
    angular = matrix_action(SY4, value)*(fields['NP']*(.5-fields['Aphysical'])/rho)[:, None, None]
    data = dict(value=value, radial=radial, angular=angular, P=fields['P'], W=fields['W'],
                f=np.zeros(len(points)), fp=np.zeros(len(points)), weight=weights)
    coords = dict(points=points, radius=radius, theta=theta, side=signs, polar_panel=panels)
    return data, coords


def self_check():
    """Focused field/domain checks; no exterior, determinant or heat solve."""
    from dihole_azimuthal_bridge import AzimuthalBridgeMatrices
    start = perf_counter()
    function, parameters = _geometry()
    model = AzimuthalBridgeMatrices(function, parameters, 1, 24, mu=.5)
    theta = np.r_[0., np.linspace(.13, np.pi-.13, 17), np.pi]
    p, pt, action = profiles(theta)
    f, g = model.profiles(theta)
    profile_error = max(np.max(abs(p[:, 0, :2]-f)), np.max(abs(p[:, 1, 2]-1j*g[:, 0])))
    inner = slice(1, -1)
    sx = np.array([[0., 1.], [1., 0.]])
    direct = (-1j*np.einsum('ab,pbj->paj', sx, pt[inner]+.5/np.tan(theta[inner])[:, None, None]*p[inner])
              +(.5/np.tan(theta[inner]))[:, None, None]*np.einsum('ab,pbj->paj', SY, p[inner]))
    angular_error = float(np.max(abs(direct-action[inner])))
    seam_error, seam_derivative_error, normal_error, derivative_error = 0., 0., 0., 0.
    phase = .7
    at = theta[inner]
    for side, sign in (('left', -1), ('right', 1)):
        points = np.column_stack((RC*np.sin(at), sign*(CENTER+RC*np.cos(at))))
        value, dr, dz = exterior_values(points, phase)
        U = rotation(at, side)
        U4 = np.zeros((len(at), 4, 4), complex)
        U4[:, :2, :2] = U4[:, 2:, 2:] = U
        phase_inverse = np.exp(-1j*phase) if sign == 1 else 1.
        recovered = phase_inverse*RC*np.einsum('pba,pbj->paj', U4.conj(), value)
        expected = bridge_values(sign*.1, at)['value']
        seam_error = max(seam_error, float(np.max(abs(recovered-expected))))
        tangent = RC*(np.cos(at)[:, None, None]*dr-sign*np.sin(at)[:, None, None]*dz)
        tangent += .5j*sign*matrix_action(np.kron(np.eye(2), SY), value)
        recovered_theta = phase_inverse*RC*np.einsum('pba,pbj->paj', U4.conj(), tangent)
        expected_theta = np.zeros_like(expected)
        for row in column_metadata():
            if row['side'] == side:
                comp = slice(0, 2) if row['chirality'] == 'R' else slice(2, 4)
                expected_theta[:, comp, row['index']] = pt[inner, :, ('zero', 'upper1', 'lower1').index(row['profile'])]
        seam_derivative_error = max(seam_derivative_error, float(np.max(abs(recovered_theta-expected_theta))))
        normal = -np.sin(at)[:, None, None]*SX4-sign*np.cos(at)[:, None, None]*SZ4
        transformed = np.einsum('pba,pbc,pcd->pad', U4.conj(), normal, U4)
        normal_error = max(normal_error, float(np.max(abs(transformed+sign*np.diag(SIGMA)))))
        sample = np.column_stack((.073*np.sin(at[::4]), sign*(CENTER+.073*np.cos(at[::4]))))
        _, exact_rho, exact_z = exterior_values(sample, phase)
        h = 2e-6
        for axis, exact in enumerate((exact_rho, exact_z)):
            delta = np.zeros(2); delta[axis] = h
            fm2, fm1, fp1, fp2 = [exterior_values(sample+j*delta, phase)[0] for j in (-2, -1, 1, 2)]
            fd = (fm2-8*fm1+8*fp1-fp2)/(12*h)
            derivative_error = max(derivative_error, float(np.max(abs(fd-exact))/(1+np.max(abs(exact)))))
    # Actual axis gauge/metric data fix the indicial matrices, not appearance.
    axes = np.array([[0., sign*(CENTER+r*c)] for sign in (-1, 1) for r in (.05, .073, .1) for c in (1., -1.)])
    axis_fields = coefficient_values(axes)
    axis_gauge_error = float(np.max(abs(axis_fields['Aphysical']-np.tile([0., 1.], 6))))
    axis_metric_error = float(np.max(abs(axis_fields['NP']-1)))
    polar_matrix_error, forbidden_error = 0., 0.
    for sign in (-1, 1):
        polar = np.array([[0., sign*(CENTER+.073)], [0., sign*(CENTER-.073)]])
        v = exterior_values(polar)[0]
        for j, charge in enumerate((.5, -.5)):
            indicial = SX4/2+1j*charge*SY4
            polar_matrix_error = max(polar_matrix_error, float(np.max(abs(indicial@v[j]))))
        small = np.array([1e-4, 5e-5, np.pi-1e-4, np.pi-5e-5])
        pp = np.column_stack((.073*np.sin(small), sign*(CENTER+.073*np.cos(small))))
        vv = exterior_values(pp)[0]
        for broad, narrow, components in ((0, 1, [1, 3]), (2, 3, [0, 2])):
            forbidden_error = max(forbidden_error, float(np.max(abs(vv[narrow, components]-vv[broad, components]/2))))
    outside = np.array([[ROUT*np.sin(t), sign*(CENTER+ROUT*np.cos(t))] for sign in (-1, 1) for t in at])
    outer = exterior_values(outside)
    outer_error = max(float(np.max(abs(v))) for v in outer)
    ext, coords = exterior_data(phase, order=16)
    volume_error = abs(float(ext['weight'].sum())-4*(ROUT**3-RC**3)/3)
    polar_old_zero = bool(np.all(ext['value'][coords['polar_panel'] != 1, :, :48] == 0))
    old_exact, scalar_exact, nested_error, chiral_error = True, True, 0., 0.
    for old, new in zip(old_bridge_chunks(.5, phase, 5, 12, 16, (3, 4)), bridge_chunks(phase, 16)):
        old_exact &= all(np.array_equal(old[k], new[k][:, :, :312]) for k in FIELD_KEYS)
        scalar_exact &= all(np.array_equal(old[k], new[k]) for k in old if k not in FIELD_KEYS)
        for k in FIELD_KEYS:
            a = new[k][:, :, 312:]
            sign = 1 if k == 'value' else -1
            chiral_error = max(chiral_error, float(np.max(abs(a[:, :2, ::2]-sign*a[:, 2:, 1::2]))))
    for large, small in zip(exterior_values(coords['points'], phase), exterior_values(coords['points'], phase, False)):
        nested_error = max(nested_error, float(np.max(abs(large[:, :, :4]-small))))
    tests = dict(reference_profiles=profile_error < 2e-14, reference_angular_action=angular_error < 2e-14,
                 seam_value=seam_error < 3e-13, seam_theta_derivative=seam_derivative_error < 5e-13,
                 normal_current=normal_error < 2e-14, cartesian_derivatives=derivative_error < 2e-8,
                 actual_axis_gauge=axis_gauge_error < 2e-13, actual_axis_metric=axis_metric_error < 2e-13,
                 indicial_domain=polar_matrix_error < 3e-13 and forbidden_error < 2e-10,
                 flat_outer_taper=outer_error < 1e-20, physical_collar_measure=volume_error < 1e-16,
                 old_polar_support=polar_old_zero, old_bridge_preserved=old_exact and scalar_exact,
                 nested_zero_fields=nested_error == 0, exact_chiral_pairing=chiral_error == 0)
    return dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                errors=dict(profile=float(profile_error), angular=angular_error, seam=seam_error,
                    seam_theta=seam_derivative_error, normal=normal_error, derivatives=derivative_error,
                    axis_gauge=axis_gauge_error, axis_metric=axis_metric_error, indicial=polar_matrix_error,
                    forbidden_component_halving=forbidden_error, outer=outer_error, volume=volume_error,
                    nested=nested_error, chiral=chiral_error), columns=column_metadata(), seconds=perf_counter()-start,
                scope=['Fixed q=1,mu=1/2 on the existing declared spin/gauge sewing; phase applies to both Weyls.',
                       'Four zero-profile traces precede eight first-pair traces. No extra parity or gauge-charge flip.',
                       'Regular monopole profiles and actual axis coefficients establish the local indicial domain.',
                       'Linear bridge lifts give H1 conformity; no matching normal derivatives or operator-domain claim.',
                       'Tests are floating controls, not continuum/source or quadrature-error bounds.'])


if __name__ == '__main__':
    result = self_check()
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    project_path('dihole_mouth_fields_checks.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))
    if not result['passed']:
        raise SystemExit('Mouth-field controls failed.')
