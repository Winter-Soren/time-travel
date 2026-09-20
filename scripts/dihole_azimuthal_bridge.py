"""Actual q=1 bridge at arbitrary half-odd-integer mu and physical Dirac mass.

Write a=mu-1/2.  In the reference upper-mouth gauge A=(1-cos(theta))/2,
K=-i*sigma_x*d_theta+sigma_y*(a+cos(theta)/2)/sin(theta).  Our profiles omit
sqrt(sin(theta)): the upper component is sqrt(sin)*f_l and the lower is
i*sqrt(sin)*g_l, with g_l=-(f_l'-a*csc(theta)*f_l)/sqrt(l*(l+1)).
The normalized associated Legendre f_l has its Condon--Shortley sign removed,
as in the frozen complete_angular_modes helper.  This fixes K_pair=+lambda*sx.

The spherical functions are a coordinate basis, not an approximation of the
actual angular coefficients.  We compress the actual N,P,A of the prescribed
bridge, with transverse correction [N*P*(mu-A)-(a+cos(theta)/2)]/(r*sin(theta)).
Both physical Weyls are present in coefficients(): plus=[R_up,L_down],
minus=[R_down,L_up].  Their off-diagonal physical mass is mass*P.  No endpoint
sewing sign is inserted here; use profiles at the specified endpoint angle.
No exterior solve, renormalized source, or self-consistent geometry is supplied.
"""
from functools import lru_cache
import hashlib
import json
import math
import operator
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import block_diag
from scipy.special import gammaln, lpmv

from balanced_dihole import symbolic_fields
from complete_angular_modes import _sector, mouth_characteristics, zero_mode
from dihole_bridge_scattering import collar_mask, complex_record, embed
from dihole_direct_source import direct_response
from dihole_massive_bridge import MassiveBridgeMatrices
from dihole_stretched_collar import exact_target_parameters, field_values, gauss, radius_profile
from project_paths import project_path
from spherical_weyl_leads import _mouth_arguments


class AzimuthalBridgeMatrices:
    """Finite regular angular compression; size counts paired modes, not lmax.

    profiles(theta) returns real (f,g) arrays shaped (len(theta),p/m), with
    dcos(theta) normalization.  plus_labels/minus_labels are their absolute
    l labels.  Only a=0 has an additional upper l=0 mode.  coefficients is
    full-Dirac embed-compatible; weyl_coefficients has dimensions p,m.
    """

    def __init__(self, function, parameters, size, order, mu=.5, mass=0.):
        try:
            if isinstance(size, bool) or isinstance(order, bool):
                raise TypeError
            size, order = operator.index(size), operator.index(order)
        except TypeError as exc:
            raise ValueError('size and order must be integers.') from exc
        self.a = _sector(mu)
        if size < 1:
            raise ValueError('size must include at least one paired angular mode.')
        first = max(1, abs(self.a))
        if order <= first+size:
            raise ValueError('Angular quadrature order must exceed the largest l plus one.')
        if isinstance(mass, bool) or not math.isfinite(float(mass)) or mass < 0:
            raise ValueError('Physical mass must be finite and nonnegative.')
        if parameters.get('q') != 1.:
            raise ValueError('This regular angular basis requires q=1.')
        self.function, self.parameters = function, dict(parameters)
        self.size, self.order, self.mu, self.mass = size, order, self.a+.5, float(mass)
        self.p, self.m = size+int(self.a == 0), size
        self.n = self.p+self.m
        self.plus = self.minus = self.n
        self.minus_labels = np.arange(first, first+size, dtype=int)
        self.plus_labels = (np.r_[0, self.minus_labels] if self.a == 0
                            else self.minus_labels.copy())
        self.paired_labels = self.minus_labels.copy()
        self.upper_labels, self.lower_labels = self.plus_labels, self.minus_labels
        self.cosine, self.weights = gauss(order)
        self.theta = np.arccos(self.cosine)
        self.upper, self.lower = self.profiles(self.theta)
        self.gram_error = float(max(
            np.max(abs(self.upper.T@(self.weights[:, None]*self.upper)-np.eye(self.p))),
            np.max(abs(self.lower.T@(self.weights[:, None]*self.lower)-np.eye(self.m)))))

    def profiles(self, theta):
        """Pole-safe scalar f,g; includes upper l=0 only in the a=0 sector."""
        theta, _, ep, em = _mouth_arguments(theta, 1., 'upper')
        plus, minus = mouth_characteristics(self.mu, theta, 1., 'upper',
                                            int(self.minus_labels[-1]))
        f = np.sum(ep[:, :, None]*minus, axis=0).real
        g = np.sum(em[:, :, None]*plus, axis=0).imag
        if self.a == 0:
            zeroth = np.sum(ep*zero_mode(self.mu, theta, 1., 'upper'), axis=0)
            f = np.column_stack((zeroth, f))
        return f, g

    @lru_cache(maxsize=4096)
    def weyl_at_radius(self, radius):
        """Return positive clock blocks Wp,Wm and the p by m transverse X."""
        radius = float(radius)
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError('Radius must be finite and positive.')
        # Match the frozen mu=.5 angular routine's coordinate evaluation.
        cosine = np.cos(self.theta)
        data = field_values(self.function, self.parameters, radius, cosine)
        lapse, meridian, potential = (np.asarray(data[key], float) for key in
                                     ('lapse', 'meridian_factor', 'unit_A_phi'))
        if (not all(np.all(np.isfinite(v)) for v in (lapse, meridian, potential))
                or np.any(lapse <= 0) or np.any(meridian <= 0)):
            raise ValueError('Actual angular metric and potential must be finite; N,P positive.')
        optical = meridian/lapse
        delta = (lapse*meridian*(self.mu-potential)-(self.a+.5*cosine))/np.sin(self.theta)
        cross = self.upper.T@((self.weights*delta/radius)[:, None]*self.lower)
        paired_rows = np.arange(self.m)+int(self.a == 0)
        cross[paired_rows, np.arange(self.m)] += np.sqrt(
            self.minus_labels*(self.minus_labels+1))/radius
        return (self.upper.T@((self.weights*optical)[:, None]*self.upper),
                self.lower.T@((self.weights*optical)[:, None]*self.lower), cross)

    @lru_cache(maxsize=4096)
    def mass_blocks_at_radius(self, radius):
        """Actual meridian-factor compressions multiplying fixed proper mass."""
        radius = float(radius)
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError('Radius must be finite and positive.')
        data = field_values(self.function, self.parameters, radius, self.cosine)
        meridian = np.asarray(data['meridian_factor'], float)
        if not np.all(np.isfinite(meridian)) or np.any(meridian <= 0):
            raise ValueError('The actual meridian factor must be finite and positive.')
        weighted = (self.weights*meridian)[:, None]
        return self.upper.T@(weighted*self.upper), self.lower.T@(weighted*self.lower)

    @lru_cache(maxsize=4096)
    def at_radius(self, radius):
        wp, wm, cross = self.weyl_at_radius(float(radius))
        pp, pm = self.mass_blocks_at_radius(float(radius))
        full = np.block([[cross, self.mass*pp], [self.mass*pm, -cross.conj().T]])
        return block_diag(wp, wm), block_diag(wm, wp), full

    @staticmethod
    def _radius_factor(s, epsilon):
        s, epsilon = float(s), float(epsilon)
        if not math.isfinite(s) or not math.isfinite(epsilon):
            raise ValueError('Coordinate and variation must be finite.')
        u = min(.05, max(0., abs(s)-.05))
        return float(radius_profile(u)[0]), math.exp(-epsilon*collar_mask(s))

    def coefficients(self, s, epsilon=0.):
        radius, factor = self._radius_factor(s, epsilon)
        wp, wm, cross = self.at_radius(radius)
        return factor*factor*wp, factor*factor*wm, factor*cross

    def weyl_coefficients(self, s, epsilon=0.):
        radius, factor = self._radius_factor(s, epsilon)
        wp, wm, cross = self.weyl_at_radius(radius)
        return factor*factor*wp, factor*factor*wm, factor*cross


def massless_permutations(p, m):
    """Index blockdiag(S_R,S_L) into the full helper's incoming/outgoing order."""
    total = p+m
    incoming = np.r_[np.arange(p), np.arange(total, total+m),
                     np.arange(p, total), np.arange(total+m, 2*total)]
    outgoing = np.r_[np.arange(m), np.arange(total, total+p),
                     np.arange(m, total), np.arange(total+p, 2*total)]
    return outgoing, incoming


def _angular_controls(model):
    """Independent associated-Legendre derivatives and regular endpoint powers."""
    at = np.linspace(.13, np.pi-.13, 37)
    cc, sn, aa = np.cos(at), np.sin(at), abs(model.a)
    ff, gg = [], []
    for ell in model.minus_labels:
        norm = (-1.)**aa*np.exp(.5*(np.log((2*ell+1)/2)
                              +gammaln(ell-aa+1)-gammaln(ell+aa+1)))
        poly = lpmv(aa, ell, cc)
        previous = lpmv(aa, ell-1, cc) if ell > aa else np.zeros_like(cc)
        derivative = (ell*cc*poly-(ell+aa)*previous)/(cc*cc-1)
        f, ft = norm*poly, -sn*norm*derivative
        ff.append(f)
        gg.append(-(ft-model.a*f/sn)/math.sqrt(ell*(ell+1)))
    independent_f, independent_g = np.array(ff).T, np.array(gg).T
    f, g = model.profiles(at)
    paired_f = f[:, int(model.a == 0):]
    associated_error = max(np.max(abs(paired_f-independent_f)), np.max(abs(g-independent_g)))
    step = 1e-4
    shifted = [model.profiles(at+t*step) for t in (-2, -1, 1, 2)]
    derivatives = tuple((a-8*b+8*c-d)/(12*step)
                        for a, b, c, d in zip(*shifted))
    ft = derivatives[0][:, int(model.a == 0):]
    gt = derivatives[1]
    lam = np.sqrt(model.minus_labels*(model.minus_labels+1))[None, :]
    action_error = max(np.max(abs(ft-model.a*paired_f/sn[:, None]+lam*g)),
                       np.max(abs(gt+(cc+model.a)[:, None]*g/sn[:, None]-lam*paired_f)))
    pole_angles = np.array([0., np.pi])
    fp, gp = model.profiles(pole_angles)
    expected_f, expected_g = np.zeros_like(fp), np.zeros_like(gp)
    cn = np.sqrt((2*model.minus_labels+1)/2)
    if model.a == 0:
        expected_f[:, 0] = 1/math.sqrt(2)
        expected_f[0, 1:], expected_f[1, 1:] = cn, (-1.)**model.minus_labels*cn
    elif model.a == 1:
        expected_g[1] = (-1.)**(model.minus_labels-1)*cn
    elif model.a == -1:
        expected_g[0] = -cn
    pole_error = max(np.max(abs(fp-expected_f)), np.max(abs(gp-expected_g)))
    # The two pole charts differ: g has powers |a|-1 and |a|+1, not |a| at both.
    powers_f = (aa, aa)
    powers_g = ((1, 1) if model.a == 0 else
                (aa+1, aa-1) if model.a > 0 else (aa-1, aa+1))
    power_error = 0.
    for end, direction, pf, pg in ((0., 1., powers_f[0], powers_g[0]),
                                   (np.pi, -1., powers_f[1], powers_g[1])):
        big = model.profiles(np.array([end+direction*1e-3]))
        small = model.profiles(np.array([end+direction*.5e-3]))
        for x, y, power in zip(big, small, (pf, pg)):
            power_error = max(power_error, float(np.max(abs(y/(x*2.**(-power))-1))))
    return dict(mu=model.mu, a=model.a, plus_labels=model.plus_labels.tolist(),
                minus_labels=model.minus_labels.tolist(), p=model.p, m=model.m,
                gram_error=model.gram_error, associated_legendre_error=float(associated_error),
                reference_partner_error=float(action_error), pole_limit_error=float(pole_error),
                pole_power_relative_error=power_error)


def _physical_block_control(model):
    p, m, n = model.p, model.m, model.n
    order = np.r_[np.arange(p), np.arange(n+p, 2*n), np.arange(p, n), np.arange(n, n+p)]
    current = np.diag(np.r_[np.ones(p), -np.ones(m)])
    rows = []
    for radius in (.025, .037, .05):
        wp, wm, x = model.weyl_at_radius(radius)
        pp, pm = model.mass_blocks_at_radius(radius)
        transverse = np.block([[np.zeros((p, p)), x], [x.conj().T, np.zeros((m, m))]])
        metric = block_diag(pp, pm)
        physical = np.block([[transverse, model.mass*metric], [model.mass*metric, -transverse]])
        wplus, wminus, full = model.at_radius(radius)
        reordered = np.block([[np.zeros((n, n)), full], [full.conj().T, np.zeros((n, n))]])
        old_clock = block_diag(wp, wm, wp, wm)
        rows.append(dict(radius=radius,
            physical_mass_block_error=float(np.max(abs(physical[np.ix_(order, order)]-reordered))),
            clock_block_error=float(np.max(abs(old_clock[np.ix_(order, order)]-block_diag(wplus, wminus)))),
            current_block_error=float(np.max(abs(block_diag(current, -current)[np.ix_(order, order)]
                                                -np.diag(np.r_[np.ones(n), -np.ones(n)])))),
            hermiticity_error=float(np.max(abs(reordered-reordered.conj().T))),
            clock_minimum=float(min(np.linalg.eigvalsh(wplus)[0], np.linalg.eigvalsh(wminus)[0])),
            mass_factor_minimum=float(min(np.linalg.eigvalsh(pp)[0], np.linalg.eigvalsh(pm)[0]))))
    varied = model.coefficients(-.08, .025)
    base = model.coefficients(-.08)
    f = math.exp(-.025*collar_mask(-.08))
    return dict(mu=model.mu, rows=rows, variation_scaling_error=float(max(
        np.max(abs(a-scale*b)) for a, b, scale in zip(varied, base, (f*f, f*f, f)))))


def _constant_mass_control(model):
    """Independent proper-mass dispersion, including the unpaired zero if present."""
    p, m, n = model.p, model.m, model.n
    lapse, meridian, mass, momentum, radius = .7, 2.3, 1.3, .4, .8
    weight = meridian/lapse
    cross = np.zeros((p, m))
    cross[np.arange(m)+int(model.a == 0), np.arange(m)] = np.sqrt(
        model.minus_labels*(model.minus_labels+1))/radius
    h0 = (np.diag(np.r_[np.full(p, momentum), np.full(m, -momentum)])
          +np.block([[np.zeros((p, p)), cross], [cross.T, np.zeros((m, m))]]))/weight
    full = np.block([[cross, mass*meridian*np.eye(p)], [mass*meridian*np.eye(m), -cross.T]])
    hm = (np.diag(np.r_[np.full(n, momentum), np.full(n, -momentum)])
          +np.block([[np.zeros((n, n)), full], [full.T, np.zeros((n, n))]]))/weight
    positive = np.sqrt(np.linalg.eigvalsh(h0)**2+(mass*lapse)**2)
    expected = np.sort(np.r_[-positive, positive])
    return dict(mu=model.mu, lapse=lapse, meridian=meridian, mass=mass,
                maximum_error=float(np.max(abs(np.linalg.eigvalsh(hm)-expected))))


def checks():
    started = perf_counter()
    names = ('balanced_dihole.py', 'dihole_stretched_collar.py', 'dihole_cylinder_transfer.py',
             'dihole_bridge_scattering.py', 'dihole_direct_source.py', 'dihole_massive_bridge.py',
             'complete_angular_modes.py', 'azimuthal_modes.py', 'spherical_weyl_leads.py', 'project_paths.py')
    inputs = {name: hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names}
    function, _ = symbolic_fields()
    parameters = exact_target_parameters()
    angular = [_angular_controls(AzimuthalBridgeMatrices(function, parameters, 3, 64, mu=a+.5))
               for a in range(-3, 4)]
    legacy = MassiveBridgeMatrices(function, parameters, 3, 48, mass=1.3)
    generalized = AzimuthalBridgeMatrices(function, parameters, 3, 48, mass=1.3)
    legacy_errors = dict(profiles=float(max(np.max(abs(generalized.upper-legacy.upper)),
                                           np.max(abs(generalized.lower-legacy.lower)))),
        coefficients=float(max(np.max(abs(a-b)) for s in (-.1, -.08, 0., .08, .1)
                               for eps in (0., .021)
                               for a, b in zip(generalized.coefficients(s, eps), legacy.coefficients(s, eps)))))
    physical, scattering, dispersion = [], [], []
    for mu in (-.5, .5, 1.5):
        model = AzimuthalBridgeMatrices(function, parameters, 2, 32, mu=mu, mass=1.3)
        zero = AzimuthalBridgeMatrices(function, parameters, 2, 32, mu=mu, mass=0.)
        physical.append(_physical_block_control(model))
        dispersion.append(_constant_mass_control(model))
        def left(s):
            wp, wm, x = zero.weyl_coefficients(s)
            return wm, wp, -x.conj().T
        right_s, _ = embed(zero.weyl_coefficients, zero.p, zero.m, .1j)
        left_s, _ = embed(left, zero.m, zero.p, .1j)
        joined, _ = embed(zero.coefficients, zero.n, zero.n, .1j)
        expected = block_diag(right_s, left_s)[np.ix_(*massless_permutations(zero.p, zero.m))]
        _, real_info = embed(model.coefficients, model.n, model.n, .2)
        imaginary_s, imaginary_info = embed(model.coefficients, model.n, model.n, .1j)
        scattering.append(dict(mu=mu, p=model.p, m=model.m,
            massless_permutation_error=float(np.max(abs(joined-expected))),
            real_energy=real_info, imaginary_energy=imaginary_info,
            imaginary_scattering=complex_record(imaginary_s)))
    direct_model = AzimuthalBridgeMatrices(function, parameters, 2, 32, mu=-.5, mass=1.3)
    direct, _ = direct_response(direct_model.coefficients, direct_model.n, direct_model.n, .1)
    eta = 1e-4
    finite = [direct_response(lambda s, e=e: direct_model.coefficients(s, e),
                             direct_model.n, direct_model.n, .1)[0]['direct_logdet']
              for e in (-eta, eta)]
    fd_error = abs(direct['direct_dlog']-(finite[1]-finite[0])/(2*eta))
    validation = []
    for kwargs in (dict(size=0), dict(order=3), dict(mu=0.), dict(mu=float('nan')),
                   dict(mass=-1.), dict(size=True), dict(parameters=dict(parameters, q=2.))):
        args = dict(function=function, parameters=parameters, size=2, order=32)
        args.update(kwargs)
        try:
            AzimuthalBridgeMatrices(**args)
        except ValueError:
            validation.append(True)
        else:
            validation.append(False)
    physical_rows = [row for control in physical for row in control['rows']]
    tests = dict(
        correct_channel_counts=all(r['p'] == 3+int(r['a'] == 0) and r['m'] == 3 for r in angular),
        angular_normalization=max(r['gram_error'] for r in angular) < 3e-12,
        independent_associated_legendre=max(r['associated_legendre_error'] for r in angular) < 5e-11,
        reference_differential_modes=max(r['reference_partner_error'] for r in angular) < 3e-8,
        regular_pole_limits=max(r['pole_limit_error'] for r in angular) < 3e-12,
        regular_pole_powers=max(r['pole_power_relative_error'] for r in angular) < 2e-5,
        legacy_mu_half_profiles=legacy_errors['profiles'] == 0.,
        actual_legacy_coefficients=legacy_errors['coefficients'] < 2e-11,
        physical_full_dirac_permutation=max(max(r['physical_mass_block_error'], r['clock_block_error'],
                                                r['current_block_error']) for r in physical_rows) < 2e-12,
        hermitian_transverse=max(r['hermiticity_error'] for r in physical_rows) < 2e-12,
        positive_clock_and_mass=min(min(r['clock_minimum'], r['mass_factor_minimum']) for r in physical_rows) > 0.,
        paired_variation=max(r['variation_scaling_error'] for r in physical) < 2e-12,
        massless_two_weyl_scattering=max(r['massless_permutation_error'] for r in scattering) < 2e-8,
        real_energy_unitarity=max(r['real_energy']['unitarity_error'] for r in scattering) < 2e-8,
        imaginary_passivity=min(r['imaginary_energy']['minimum_passivity_defect'] for r in scattering) > -2e-9,
        constant_proper_mass_dispersion=max(r['maximum_error'] for r in dispersion) < 3e-12,
        direct_metric_derivative=fd_error < 3e-6,
        invalid_input_rejection=all(validation),
        historical_inputs_preserved=all(hashlib.sha256(project_path(name).read_bytes()).hexdigest() == digest
                                        for name, digest in inputs.items()))
    return dict(complete=True, passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
        angular_controls=angular, legacy_mu_half_errors=legacy_errors,
        physical_block_controls=physical, scattering_controls=scattering,
        proper_mass_controls=dispersion,
        direct_control=dict(mu=-.5, kappa=.1, result=direct, metric_derivative_fd_error=fd_error),
        interface=dict(class_name='AzimuthalBridgeMatrices',
            constructor='(function, parameters, size, order, mu=.5, mass=0.)',
            profiles='profiles(theta) -> f[point,p],g[point,m]; normalization measure dcos(theta)',
            labels='plus_labels includes l=0 only for mu=.5; minus_labels has the paired absolute l labels',
            coefficients='full Dirac: plus=[R_up,L_down], minus=[R_down,L_up], each size n=p+m',
            weyl_coefficients='right Weyl: plus dimension p, minus dimension m'),
        scope=['Actual prescribed q1 bridge at finite angular size for each half-odd-integer mu; tested a=-3 through 3.',
               'size counts paired labels; a common size across sectors is not a complete fixed-l multiplet cutoff.',
               'Both physical Weyls are explicit. No extra Weyl factor belongs in a full determinant prefactor.',
               'Angular endpoint sewing is external: use the stated profiles and frame map, not an assumed parity.',
               'ODE and quadrature checks are numerical diagnostics, not continuum or roundoff enclosures.',
               'No exterior return, regulator sum, renormalized stress, self-consistent bridge, or chronology claim.'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), input_sha256=inputs,
        runtime_seconds=perf_counter()-started)


if __name__ == '__main__':
    result = checks()
    project_path('dihole_azimuthal_bridge_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=result['tests'],
        legacy_mu_half_errors=result['legacy_mu_half_errors'], runtime_seconds=result['runtime_seconds']), indent=2))
    if not result['passed']:
        raise SystemExit('Azimuthal bridge checks failed; diagnostic evidence retained.')
