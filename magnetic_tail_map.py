"""Outgoing spherical map with the actual flat point-pair magnetic tail.

The pair field is retained on R<=r<=L. Free outgoing data at L still omit
the field beyond L; far-radius and angular-channel refinements are separate.
This module does not change the declared mouth leads or solve gravity.
"""
import hashlib
import json
import operator
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from exterior_boundary_map import pair_potential, SX, SY, SZ
from spherical_weyl_leads import outer_characteristics, outgoing_ratios


class _AngularCoupling:
    """Precompute the normalized angular products, not a sampled radial field."""

    def __init__(self, q=1., d=1., mu=.5, nmax=12, angular_order=64):
        self.q, self.d, self.mu = map(float, (q, d, mu))
        if not np.all(np.isfinite([self.q, self.d, self.mu])) or self.d <= 0:
            raise ValueError('Require finite q,mu and positive d')
        try:
            self.nmax = operator.index(nmax)
            self.order = operator.index(angular_order)
        except TypeError as exc:
            raise ValueError('Channel cutoff and quadrature order must be integers') from exc
        if self.order < max(8, self.nmax+1):
            raise ValueError('Use angular_order>=max(8,nmax+1)')
        self.cosine, self.weights = np.polynomial.legendre.leggauss(self.order)
        self.theta = np.arccos(self.cosine)
        self.sine = np.sqrt((1-self.cosine)*(1+self.cosine))
        self.plus, self.minus = outer_characteristics(
            self.mu, self.theta, 1., self.nmax)
        first = int(round(abs(self.mu)+.5))
        self.channels = np.arange(first, self.nmax+1)
        self.free = np.diag(self.channels)
        self.kernel = np.einsum('apm,ab,bpn,p->pmn',
            self.plus.conj(), SY, self.minus, self.weights/self.sine)
        self.gram_error = max(float(np.max(abs(
            np.einsum('apm,apn,p->mn', basis.conj(), basis, self.weights)
            -np.eye(len(self.channels))))) for basis in (self.plus, self.minus))

    def matrix(self, radius):
        radius = float(radius)
        if not np.isfinite(radius) or radius <= self.d/2:
            raise ValueError('The full spherical expansion requires r>d/2')
        if self.q == 0:
            return self.free/radius
        potential = pair_potential(radius*self.sine, radius*self.cosine,
                                   q=self.q, d=self.d)
        J = np.einsum('p,pmn->mn', potential, self.kernel)
        return (self.free+1j*J)/radius


def _flux_diagnostics(Z):
    flux = np.eye(len(Z))-Z.conj().T@Z
    hermitian = (flux+flux.conj().T)/2
    eigenvalues = np.linalg.eigvalsh(hermitian)
    return dict(flux_metric_eigenvalues=eigenvalues.tolist(),
        minimum_flux_metric_eigenvalue=float(eigenvalues[0]),
        maximum_singular_value=float(np.linalg.svd(Z, compute_uv=False)[0]),
        flux_metric_hermiticity_error=float(np.max(abs(flux-flux.conj().T))),
        eigenvalues_clipped=False)


def tail_map(energy, R, L, q=1., d=1., mu=.5, nmax=12,
             angular_order=64, rtol=1e-10, atol=1e-12):
    """Return (complex Z, diagnostics), with incoming coefficients b=Z*a.

    Coordinates are centered halfway between the foci. E>0, L>R>d/2.
    The basis is exactly the free spherical characteristic convention of
    spherical_weyl_leads, with canonical radial coefficients u=Phi*a/r.
    Error tolerances are ODE controls, not certified global error bounds.
    """
    energy, R, L, rtol, atol = map(float, (energy, R, L, rtol, atol))
    angular = _AngularCoupling(q, d, mu, nmax, angular_order)
    if (not np.all(np.isfinite([energy, R, L, rtol, atol])) or energy <= 0
            or not L > R > angular.d/2 or not 0 < rtol < 1 or not 0 < atol < 1):
        raise ValueError('Require E>0, L>R>d/2 and finite positive ODE tolerances<1')
    count = len(angular.channels)
    initial = np.diag(outgoing_ratios(energy, L, mu, nmax))

    def rhs(radius, packed):
        Z = packed.reshape(count, count)
        C = angular.matrix(radius)
        return (C.conj().T-2j*energy*Z-Z@C@Z).reshape(-1)

    solution = solve_ivp(rhs, (L, R), initial.reshape(-1), method='DOP853',
                         rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError('Magnetic tail integration failed: '+solution.message)
    Z = solution.y[:, -1].reshape(count, count)
    if not np.all(np.isfinite(Z)):
        raise ArithmeticError('Nonfinite magnetic tail map')
    diagnostics = dict(energy=energy, R=R, L=L, q=float(q), d=float(d),
        mu=float(mu), nmax=angular.nmax, channels=angular.channels.tolist(),
        angular_order=angular.order, rtol=rtol, atol=atol,
        nfev=int(solution.nfev), accepted_steps=len(solution.t)-1,
        angular_gram_error=angular.gram_error, **_flux_diagnostics(Z))
    return Z, diagnostics


def _complex_record(value):
    value = np.asarray(value)
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def _operator_controls():
    angular = _AngularCoupling(nmax=3, angular_order=80)
    c11_errors = [abs(angular.matrix(r)[0, 0]-(1/r-1/(3*r*r)))
                  for r in (.75, 2., 10.)]
    C = angular.matrix(2.)
    identity = np.eye(3)
    generator = np.block([[.4j*identity, C], [C.conj().T, -.4j*identity]])
    metric = np.diag([1.]*3+[-1.]*3)

    # Independently apply the differential angular Weyl operator to Phi_minus.
    # This verifies the phase/normalization of J, not just its Hermitian assembly.
    theta, cosine, sine = angular.theta, angular.cosine, angular.sine
    step = 1e-5
    minus = lambda t: outer_characteristics(.5, t, 1., 3)[1]
    derivative = (-minus(theta+2*step)+8*minus(theta+step)
                  -8*minus(theta-step)+minus(theta-2*step))/(12*step)
    sigma_r = sine[:, None, None]*SX+cosine[:, None, None]*SZ
    sigma_theta = cosine[:, None, None]*SX-sine[:, None, None]*SZ
    potential = pair_potential(2*sine, 2*cosine)
    applied = (1j*np.einsum('pab,bpn->apn', sigma_r, angular.minus)
        -1j*np.einsum('pab,bpn->apn', sigma_theta, derivative)
        -.5j*np.einsum('ab,bpn,p->apn', SX, angular.minus, 1/sine)
        +np.einsum('ab,bpn,p->apn', SY, angular.minus, (.5-potential)/sine))
    projected = np.einsum('apm,apn,p->mn', angular.plus.conj(), applied,
                          angular.weights)
    # C=-i*K_plus_minus, with the dimensionless angular operator divided by r.
    angular_error = float(np.max(abs(projected-2j*C)))
    E, radius = .4, 2.
    z = -1/(1-2j*E*radius)
    free_n1_error = abs(-2j*E/(1-2j*E*radius)**2
                        -(1/radius-2j*E*z-z*z/radius))
    return dict(max_C11_exact_error=float(max(c11_errors)),
        C_at_R2=_complex_record(C), C12_real=float(C[0, 1].real),
        C21_real=float(C[1, 0].real),
        generator_flux_metric_error=float(np.max(abs(generator.conj().T@metric+metric@generator))),
        independent_differential_angular_projection_error=angular_error,
        free_n1_riccati_identity_error=float(free_n1_error))


def _fundamental_control():
    """Integrate the linear radial system independently of the Riccati RHS."""
    E, R, L, count = .4, 2., 5., 3
    angular = _AngularCoupling(nmax=count, angular_order=64)
    initial_Z = np.diag(outgoing_ratios(E, L, .5, count))
    initial = np.vstack((np.eye(count), initial_Z))

    def rhs(radius, packed):
        a, b = packed.reshape(2*count, count).reshape(2, count, count)
        C = angular.matrix(radius)
        return np.vstack((1j*E*a+C@b, C.conj().T@a-1j*E*b)).reshape(-1)

    solved = solve_ivp(rhs, (L, R), initial.reshape(-1), method='DOP853',
                       rtol=8e-13, atol=8e-15)
    if not solved.success:
        raise RuntimeError(solved.message)
    a, b = solved.y[:, -1].reshape(2, count, count)
    direct_Z = np.linalg.solve(a.T, b.T).T
    riccati_Z, _ = tail_map(E, R, L, nmax=count, rtol=8e-13, atol=8e-15)
    flux_initial = np.eye(count)-initial_Z.conj().T@initial_Z
    flux_final = a.conj().T@a-b.conj().T@b
    return dict(energy=E, R=R, L=L, nmax=count,
        maximum_map_difference=float(np.max(abs(direct_Z-riccati_Z))),
        conserved_flux_matrix_error=float(np.max(abs(flux_final-flux_initial))),
        final_plus_block_condition=float(np.linalg.cond(a)), nfev=int(solved.nfev))


def checks():
    """Small isolated controls and explicit, limited tail/truncation refinements."""
    E, R = .4, 2.
    fine = dict(rtol=3e-12, atol=3e-14)
    diagnostics, maps = {}, {}

    def evaluate(name, **kwargs):
        z, data = tail_map(E, kwargs.pop('R', R), **kwargs)
        diagnostics[name], maps[name] = data, z
        return z

    free = evaluate('free', L=32., q=0., **fine)
    free_error = float(np.max(abs(free-np.diag(outgoing_ratios(E, R, .5, 12)))))
    far_rows = []
    previous = None
    for L in (32., 64., 128., 256., 512.):
        current = evaluate('L'+str(int(L)), L=L, **fine)
        row = dict(L=L)
        if previous is not None:
            difference = current-previous
            row.update(maximum_entry_change=float(np.max(abs(difference))),
                       spectral_norm_change=float(np.linalg.norm(difference, 2)))
        far_rows.append(row)
        previous = current
    coarse = evaluate('ode_default', L=64.)
    accurate = maps['L64']
    tolerance_change = float(np.linalg.norm(coarse-accurate, 2))
    flux_change_bound = float((np.linalg.norm(coarse, 2)+np.linalg.norm(accurate, 2))*tolerance_change)
    quadrature = evaluate('angular32', L=64., angular_order=32, **fine)
    quadrature_change = float(np.max(abs(quadrature-accurate)))
    by_channels = {12: maps['L256']}
    for nmax in (8, 16, 20):
        by_channels[nmax] = evaluate('n'+str(nmax), L=256., nmax=nmax, **fine)
    channel_rows = []
    for lower, higher in ((8, 12), (12, 16), (16, 20)):
        difference = by_channels[lower]-by_channels[higher][:lower, :lower]
        channel_rows.append(dict(lower_nmax=lower, higher_nmax=higher,
            common_low_block_size=4,
            low_block_maximum_change=float(np.max(abs(difference[:4, :4]))),
            full_common_block_maximum_change=float(np.max(abs(difference)))))
    evaluate('R4_L512', R=4., L=512., **fine)
    operators = _operator_controls()
    fundamental = _fundamental_control()
    minimum_flux = min(row['minimum_flux_metric_eigenvalue'] for row in diagnostics.values())
    # An explicit numerical acceptance tolerance, not a certified continuum
    # error bound. Eigenvalues near zero are stored as computed, never clipped.
    contractivity_tolerance = 2e-10
    far_changes = [row['maximum_entry_change'] for row in far_rows[1:]]
    tests = dict(
        free_outgoing_map=free_error < 2e-10,
        exact_lowest_C_entry=operators['max_C11_exact_error'] < 2e-13,
        nonsymmetric_C_orientation=operators['C12_real'] < 0 < operators['C21_real'],
        generator_flux_metric=operators['generator_flux_metric_error'] < 1e-13,
        independent_angular_differentiation=operators['independent_differential_angular_projection_error'] < 2e-9,
        free_n1_riccati_identity=operators['free_n1_riccati_identity_error'] < 1e-14,
        independent_linear_radial_map=fundamental['maximum_map_difference'] < 2e-10,
        integrated_flux_metric=fundamental['conserved_flux_matrix_error'] < 2e-9,
        angular_quadrature=quadrature_change < 2e-10,
        ode_tolerance_refinement=tolerance_change < 2e-9,
        far_radius_changes_decrease=all(a > b for a, b in zip(far_changes[:-1], far_changes[1:])),
        far_radius_final_observed_change=far_changes[-1] < 1e-5,
        retained_low_channel_block=channel_rows[1]['low_block_maximum_change'] < 2e-10 and channel_rows[2]['low_block_maximum_change'] < 2e-10,
        numerical_contractivity_without_clipping=minimum_flux >= -contractivity_tolerance,
        angular_basis_normalization=max(row['angular_gram_error'] for row in diagnostics.values()) < 2e-12)
    tests = {key: bool(value) for key, value in tests.items()}
    assert all(tests.values()), tests
    root = Path(__file__).parent
    inputs = ('exterior_boundary_map.py', 'spherical_weyl_leads.py')
    return dict(passed=True, tests=tests, operator_controls=operators,
        independent_linear_system=fundamental, free_map_maximum_error=free_error,
        refinement=dict(energy=E, base_R=R,
            angular_order_32_to_64_maximum_change=quadrature_change,
            ode_default_to_refined_spectral_change=tolerance_change,
            computed_flux_metric_difference_upper_bound=flux_change_bound,
            far_radius=far_rows, channel_cutoff=channel_rows),
        numerical_contractivity=dict(minimum_observed_flux_eigenvalue=minimum_flux,
            allowed_numerical_tolerance=contractivity_tolerance, clipping=False,
            interpretation='Exact finite-channel outgoing maps are contractive. The tolerance only accepts bounded observed numerical negativity; ODE tolerances and refinement differences are not certified true-error bounds.'),
        cases=diagnostics,
        reference_maps=[dict(parameters=diagnostics[key], Z=_complex_record(maps[key]))
                        for key in ('L512', 'R4_L512')],
        input_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in inputs},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        equations=dict(
            basis='u=[Phi_plus*a+Phi_minus*b]/r; Phi=sqrt(2*pi)*eta; integral_-1^1 Phi_s^dagger Phi_s dc=I. Use outer_characteristics at radius=1.',
            magnetic='Aphi=(q/2)*[(r*c+d/2)/sqrt(r^2+d^2/4+r*d*c)-(r*c-d/2)/sqrt(r^2+d^2/4-r*d*c)].',
            angular='J=integral Phi_plus^dagger*sigma_y*Phi_minus*Aphi/sqrt(1-c^2) dc; C=[diag(n)+i*J]/r. No extra 2*pi or radial normalization belongs in J.',
            radial='a_prime=i*E*a+C*b; b_prime=C^dagger*a-i*E*b. C is generally nonsymmetric.',
            riccati='Z=b*a^(-1); Z_prime=C^dagger-2*i*E*Z-Z*C*Z; Z(L)=diag(free outgoing_ratios(E,L)).',
            current='j=a^dagger*a-b^dagger*b; radial generator satisfies G^dagger*diag(I,-I)+diag(I,-I)*G=0.',
            exact_C11='For mu=1/2 and r>d/2, C11=1/r-q*d/(3*r^2), since integral Aphi dc=2*q*d/(3*r).',
            flux_metric_comparison='||Q(Z1)-Q(Z2)||_2 <= (||Z1||_2+||Z2||_2)*||Z1-Z2||_2, Q=I-Z^dagger Z. This bounds the difference of computed metrics, not the continuum error.'),
        scope=[
            'Flat, stationary point-pair exterior in one fixed half-integral azimuthal sector, positive energy, and R>d/2. Units are the caller\'s common length and inverse-length energy units.',
            'The actual magnetic potential is evaluated throughout R<=r<=L; only the initialization at L replaces the remaining exterior by free space.',
            'The reported L sequence measures sensitivity to that far initialization. Its decreasing differences do not constitute an analytic infinity-error bound.',
            'Channel refinements establish the reported low block only. The highest common channels retain larger truncation sensitivity, reported separately; no uniform whole-operator convergence is asserted.',
            'Z is in the continuous spherical coefficient basis, before polygonal projection or Gram coordinate conversion. No DG map, singular-value clipping, geometry correction, or mouth-lead choice is performed here.',
            'This calculation does not upgrade the declared product mouth leads to a physical curved throat, and does not determine a vacuum source, stable medium, or backward-time passage.'])


if __name__ == '__main__':
    result = checks()
    Path(__file__).with_name('magnetic_tail_map_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in (
        'passed', 'tests', 'operator_controls', 'independent_linear_system',
        'refinement', 'numerical_contractivity', 'source_sha256')}, indent=2))
