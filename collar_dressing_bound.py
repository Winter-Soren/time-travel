"""Analytic local dressed-state residual bounds in the flat point-pair collar.

Reads frozen incoming/angular-mode records; never reruns those benchmarks.
The independent q=1 test constructs the corrector from a closed angular CDF.
No continuum, propagation, curved-interface, or vacuum-source bound is implied.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp


def _number(value):
    value = mp.mpf(str(value))
    if not mp.isfinite(value):
        raise ValueError('Parameters must be finite real numbers')
    return value


def _bound_quantities(q, k, b, d, energy):
    if (not isinstance(q, int) or isinstance(q, bool) or q <= 0
            or not isinstance(k, int) or isinstance(k, bool) or not 0 <= k < q):
        raise ValueError('Require integer q>0 and integer 0<=k<q')
    b, d, energy = map(_number, (b, d, energy))
    if not 0 < b < d:
        raise ValueError('Require 0<b<d')
    eps = b/d
    logM = -q*mp.log1p(-eps*eps)
    variance = mp.mpf((k+1)*(q-k))/((q+1)**2*(q+2))
    alpha = q*b/(d*d)
    gap0 = mp.sqrt(q+1)/b
    perturbation = q*b/(4*(d-b)**2)
    if not perturbation < gap0/2:
        raise ValueError('The sufficient isolated-gap condition delta<gap0/2 failed')
    gap = gap0-perturbation
    root_variance_M = mp.sqrt(variance)*mp.exp(logM/2)
    sigma = alpha*root_variance_M/(1-eps)**3
    rho = q/(d*d)*(1+eps+2*eps*eps)*root_variance_M/(1-eps)**5
    G = alpha/(1-eps*eps)
    chi_second = rho+G*sigma
    eta = sigma/gap
    if not eta < 1:
        raise ValueError('The sufficient positive-flux condition eta_bound<1 failed')
    C = q/(2*d*d*(1-eps)**3)
    eta_prime = (chi_second+sigma/b+C*eta)/gap
    log_flux_factor = -mp.log1p(-eta*eta)/2
    flux_factor = mp.exp(log_flux_factor)
    residual = flux_factor*(eta_prime+2*abs(energy)*eta
        +eta*eta_prime*mp.sqrt(1+eta*eta)/(1-eta*eta))
    return dict(epsilon=eps, log_M=logM, beta_variance=variance,
        gap0=gap0, angular_perturbation_bound=perturbation, gap_lower_bound=gap,
        gap_perturbation_ratio=perturbation/gap0,
        sigma_bound=sigma, rho_bound=rho, G_bound=G,
        chi_second_bound=chi_second, radial_operator_C_bound=C,
        eta_bound=eta, eta_prime_bound=eta_prime,
        flux_deficit_upper=eta*eta,
        flux_factor_minus_one_upper=mp.expm1(log_flux_factor),
        energy=energy, flux_normalized_residual_bound=residual)


def analytic_bound(q, k, b, d, energy, precision=70):
    """Evaluate the sufficient bound; decimal offsets retain tiny flux changes."""
    if not isinstance(precision, int) or precision < 50:
        raise ValueError('Use at least 50 decimal digits')
    with mp.workdps(precision):
        values = _bound_quantities(q, k, b, d, energy)
        record = {key: float(value) for key, value in values.items()}
        record['decimal_values'] = {key: mp.nstr(values[key], 45) for key in (
            'eta_bound', 'eta_prime_bound', 'flux_deficit_upper',
            'flux_factor_minus_one_upper', 'flux_normalized_residual_bound')}
        return record


def _direct_q1_check():
    """Closed-CDF corrector, angular inverse, and Cartesian Dirac checks.

    This constructs t and t_b without using any bound coefficient. Quadrature
    uses rational combinations of CDF jets, avoiding endpoint fractional powers.
    """
    from dipole_collar_modes import correction

    with mp.workdps(65):
        b, energy = mp.mpf('.2'), mp.mpf('.4')

        def integral(tau, radius):
            D = mp.sqrt((1+radius)**2-4*radius*tau)
            return (tau+1+radius-D-(1-radius)*mp.log(2/(D+1-radius)))/(1+radius)

        def cdf(tau, radius):
            return integral(tau, radius)/integral(1, radius)

        def density(tau, radius):
            D = mp.sqrt((1+radius)**2-4*radius*tau)
            return (D+1+radius)/((1+radius)*(D+1-radius)*integral(1, radius))

        def corrector_jets(tau, radius):
            p = density(tau, radius)
            pb = mp.diff(lambda r: density(tau, r), radius)
            J = mp.diff(lambda r: cdf(tau, r), radius)/2
            Jb = mp.diff(lambda r: cdf(tau, r), radius, 2)/2
            Q = -radius*J
            Qb_effective = -J-radius*Jb-Q*pb/(2*p)
            return p, pb, Q, Qb_effective

        def h(theta, radius):
            tau = mp.sin(theta/2)**2
            return mp.sqrt(density(tau, radius)*mp.sqrt(tau*(1-tau)))

        def t(theta, radius):
            tau = mp.sin(theta/2)**2
            J = mp.diff(lambda r: cdf(tau, r), radius)/2
            return -radius*J/h(theta, radius)

        quadratures = []
        for order in (24, 48):
            nodes, weights = mp.gauss_quadrature(order, 'legendre')
            terms = []
            for node, weight in zip(nodes, weights):
                tau = (node+1)/2
                p, pb, Q, T = corrector_jets(tau, b)
                denominator = p*tau*(1-tau)
                terms.append([weight*x/2 for x in (
                    Q*Q/denominator, Q*T/denominator, T*T/denominator,
                    pb*pb/(4*p), p)])
            norm_t2, inner, norm_tp2, norm_hp2, norm_h2 = [
                mp.fsum(row[j] for row in terms) for j in range(5)]
            factor = 1/mp.sqrt(1-norm_t2)
            factor_prime = factor*inner/(1-norm_t2)
            residual2 = (factor_prime**2+factor**2*norm_tp2
                +2*factor*factor_prime*inner+factor_prime**2*norm_t2
                +4*energy**2*factor**2*norm_t2)
            quadratures.append(dict(order=order, norm_t=mp.sqrt(norm_t2),
                norm_t_prime=mp.sqrt(norm_tp2), inner_t_t_prime=inner,
                norm_h_prime=mp.sqrt(norm_hp2), norm_h_squared=norm_h2,
                flux_factor=factor, flux_factor_prime=factor_prime,
                flux_normalized_residual=mp.sqrt(residual2),
                normalized_flux=factor**2*(norm_h2-norm_t2)))
        coarse, fine = quadratures
        refinement = max(abs(coarse[key]/fine[key]-1) for key in (
            'norm_t', 'norm_t_prime', 'inner_t_t_prime', 'flux_normalized_residual'))
        factor, factor_prime = fine['flux_factor'], fine['flux_factor_prime']
        sx = mp.matrix([[0, 1], [1, 0]])
        sy = mp.matrix([[0, -1j], [1j, 0]])
        sz = mp.diag([1, -1])
        angular_errors, cdf_errors, cartesian_errors = [], [], []
        for tau in map(mp.mpf, ('.17', '.53', '.83')):
            theta = 2*mp.asin(mp.sqrt(tau))
            cdf_errors.append(abs(mp.diff(lambda x: integral(x, b), tau)
                -mp.exp(2*correction(1, b, 1, tau)[0])))
            rho, z = b*mp.sin(theta), b*mp.cos(theta)
            Aphi = ((z+1)/mp.sqrt(rho*rho+(z+1)**2)-z/b)/2
            V = (mp.mpf('.5')-Aphi)/mp.sin(theta)
            hb = mp.diff(lambda r: h(theta, r), b)
            inverse = (-1j/b)*(mp.diff(lambda x: t(x, b), theta)+V*t(theta, b))
            angular_errors.append(abs(inverse-1j*hb))

            def cartesian_spinor(rho_value, z_value):
                radius = mp.sqrt(rho_value**2+z_value**2)
                angle = mp.atan2(rho_value, z_value)
                rotation = mp.matrix([[mp.cos(angle/2), -mp.sin(angle/2)],
                                      [mp.sin(angle/2), mp.cos(angle/2)]])
                # Only the independently integrated first normalizer jet is
                # needed by the first-order Dirac operator at radius b.
                normalizer = factor+factor_prime*(radius-b)
                return (mp.exp(1j*energy*radius)/mp.sqrt(radius)*normalizer
                    *rotation*mp.matrix([h(angle, radius), t(angle, radius)]))

            f = cartesian_spinor(rho, z)
            dr = mp.matrix([mp.diff(lambda x: cartesian_spinor(x, z)[j], rho) for j in (0, 1)])
            dz = mp.matrix([mp.diff(lambda x: cartesian_spinor(rho, x)[j], z) for j in (0, 1)])
            actual = -1j*sx*dr-1j*sz*dz+sy*(mp.mpf('.5')-Aphi)/rho*f-energy*f
            rotation = mp.matrix([[mp.cos(theta/2), -mp.sin(theta/2)],
                                  [mp.sin(theta/2), mp.cos(theta/2)]])
            expected = mp.exp(1j*energy*b)/mp.sqrt(b)*rotation*mp.matrix([
                -1j*factor_prime*h(theta, b),
                1j*(factor*mp.diff(lambda r: t(theta, r), b)+factor_prime*t(theta, b))
                -2*energy*factor*t(theta, b)])
            cartesian_errors.append(mp.norm(actual-expected)/mp.norm(f))
        bound = _bound_quantities(1, 0, b, 1, energy)
        return dict(q=1, k=0, b=float(b), d=1., energy=float(energy),
            quadratures=[{key: int(value) if key == 'order' else float(value)
                          for key, value in row.items()} for row in quadratures],
            quadrature_relative_difference=float(refinement),
            cdf_density_identity_error=float(max(cdf_errors)),
            cdf_endpoint_error=float(max(abs(cdf(0, b)), abs(cdf(1, b)-1))),
            angular_inverse_error=float(max(angular_errors)),
            full_cartesian_dirac_identity_error=float(max(cartesian_errors)),
            actual_over_bound={key: float(fine[actual]/bound[key]) for key, actual in (
                ('sigma_bound', 'norm_h_prime'), ('eta_bound', 'norm_t'),
                ('eta_prime_bound', 'norm_t_prime'),
                ('flux_normalized_residual_bound', 'flux_normalized_residual'))},
            bound={key: float(value) for key, value in bound.items()})


def checks():
    root = Path(__file__).parent
    frozen = {}
    hashes = {}
    for name in ('incoming_modes', 'dipole_collar_modes'):
        path = root/(name+'_results.json')
        data = json.loads(path.read_text(encoding='utf-8'))
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[name+'.py'] = hashlib.sha256((root/(name+'.py')).read_bytes()).hexdigest()
        if not data['passed'] or data['source_sha256'] != hashes[name+'.py']:
            raise ValueError('Frozen source/result mismatch: '+name)
        frozen[name] = data
    incoming, angular = frozen['incoming_modes'], frozen['dipole_collar_modes']
    saved = incoming['saved_collar']
    if not saved['passed'] or angular['saved_collar_parameters'] != saved['parameters']:
        raise ValueError('Saved collar parameters differ or their controls failed')
    if angular['frozen_input_sha256'] != hashes['incoming_modes_results.json']:
        raise ValueError('Angular result used a different incoming-mode record')
    parameters = saved['parameters']
    q, b, d = parameters['q'], parameters['isotropic_b'], parameters['separation_d']
    energies = dict(zero=0., **saved['energy_scales'])
    rows, precision_errors = [], []
    for old in angular['actual_collar']:
        k = old['k']
        row = dict(q=q, k=k, b=b, d=d, **analytic_bound(q, k, b, d, 0))
        row['saved_undressed_radial_residual'] = old['radial_residual_norm']
        row['energy_comparisons'] = []
        for label, energy in energies.items():
            bound = analytic_bound(q, k, b, d, energy)
            residual = bound['flux_normalized_residual_bound']
            row['energy_comparisons'].append(dict(label=label, energy=energy,
                residual_bound=residual,
                bound_over_saved_undressed_residual=residual/old['radial_residual_norm'],
                bound_over_energy=None if energy == 0 else residual/abs(energy)))
        with mp.workdps(60):
            coarse = _bound_quantities(q, k, b, d, 0)
        with mp.workdps(80):
            fine = _bound_quantities(q, k, b, d, 0)
            precision_errors.append(float(max(abs(coarse[key]/fine[key]-1) for key in (
                'eta_bound', 'eta_prime_bound', 'flux_factor_minus_one_upper',
                'flux_normalized_residual_bound'))))
        rows.append(row)
    direct = _direct_q1_check()
    rejected = []
    for args in ((0, 0, .2, 1, 0), (1, 1, .2, 1, 0), (1, 0, 1, 1, 0),
                 (1, 0, .9, 1, 0), (1, 0, .55, 1, 0), (1, 0, .2, 1, 'nan')):
        try:
            analytic_bound(*args)
        except ValueError:
            rejected.append(True)
        else:
            rejected.append(False)
    tests = dict(
        frozen_inputs_match=True,
        sufficient_domain_and_gap_flux_guards=all(rejected),
        direct_closed_cdf_density_and_endpoints=max(direct['cdf_density_identity_error'], direct['cdf_endpoint_error']) < 1e-50,
        direct_angular_inverse=direct['angular_inverse_error'] < 1e-50,
        direct_cartesian_dirac_residual=direct['full_cartesian_dirac_identity_error'] < 1e-45,
        independent_corrector_quadrature=direct['quadrature_relative_difference'] < 1e-27,
        direct_normalized_flux=abs(direct['quadratures'][-1]['normalized_flux']-1) < 1e-14,
        direct_corrector_and_residual_below_bounds=all(0 < value < 1 for value in direct['actual_over_bound'].values()),
        saved_radial_residual_below_sigma_bound=all(row['saved_undressed_radial_residual'] <= row['sigma_bound'] for row in rows),
        high_precision_evaluation=max(precision_errors) < 1e-48,
        tiny_flux_factor_offset_retained=all(row['flux_factor_minus_one_upper'] > 0 for row in rows))
    assert all(tests.values()), tests
    return dict(passed=True, tests=tests, actual_collar=rows, independent_q1=direct,
        parameters=parameters, sampled_energy_scales=energies,
        high_precision_relative_differences=precision_errors,
        input_sha256=hashes,
        upstream_closure_input_sha256=saved['input_sha256'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        units='re=1 and hbar=c=1. b,d in re; E, sigma, eta_prime and residual in inverse re; eta and flux offsets dimensionless. Saved exterior energy numbers are comparisons in the flat model, without a new curved tetrad matching.',
        equations=dict(
            angular='H=-i*sigma_z*d_b+K; K=(-i/b)*[[0,d_theta+V],[d_theta-V,0]]; h_theta=V*h; integral h^2 dtheta=1.',
            radial_source='h_b=(g-mean(g))*h; sigma=||h_b||; J(theta)=integral_0^theta h*h_b ds; t=-b*J/h.',
            inverse='K*(0,t)=i*(h_b,0), with t in the negative-chirality complement.',
            gap='gap0=sqrt(q+1)/b; delta=q*b/[4*(d-b)^2]<gap0/2; Delta=gap0-delta.',
            moment_bounds='eps=b/d; M=(1-eps^2)^(-q); v=(k+1)*(q-k)/[(q+1)^2*(q+2)]; alpha=q*b/d^2; sigmaB=alpha*sqrt(v*M)/(1-eps)^3; rhoB=(q/d^2)*(1+eps+2eps^2)*sqrt(v*M)/(1-eps)^5.',
            derivative_bound='G=alpha/(1-eps^2); chiSecondB=rhoB+G*sigmaB; etaB=sigmaB/Delta; C=q/[2*d^2*(1-eps)^3]; etaPrimeB=(chiSecondB+sigmaB/b+C*etaB)/Delta.',
            raw_residual='For exp(i*E*b)*(h,t), (H-E)state=(0,i*t_b-2*E*t); raw outward flux=1-||t||^2.',
            normalized_residual='F=(1-||t||^2)^(-1/2); residual=(-i*Fprime*h, i*(F*t_b+Fprime*t)-2*E*F*t).',
            bound='[etaPrimeB+2*abs(E)*etaB+etaB*etaPrimeB*sqrt(1+etaB^2)/(1-etaB^2)]/sqrt(1-etaB^2), requiring etaB<1.',
            q1_cdf='D=sqrt((1+b)^2-4*b*tau); I=[tau+1+b-D-(1-b)*log(2/(D+1-b))]/(1+b); CDF=I/I(1,b); P=d_tau I/I(1,b); J=.5*d_b CDF.'),
        proof_outline=[
            'The exact zero-mode count and the bounded angular perturbation give the isolated complementary gap under delta<gap0/2; solve the angular inverse on negative chirality.',
            'The normalized Beta-density tilt and radial derivative bounds give ||h_b||<=sigmaB and ||h_bb||<=chiSecondB.',
            'Differentiate K*t=i*h_b using K_b=-K/b plus a multiplication operator bounded by C, yielding ||t_b||<=etaPrimeB.',
            'Cauchy-Schwarz gives abs(Fprime)<=etaB*etaPrimeB/(1-etaB^2)^(3/2). Apply the triangle inequality to the normalized residual and ||(h,t)||<=sqrt(1+etaB^2).'],
        scope=[
            'Analytic sufficient local bound for a constructed flux-normalized dressed state on an upper-centered sphere in the flat point-pair field, with 0<b<d and the stated gap/flux guards.',
            'Reported decimal bounds use high-precision arithmetic, not outward-rounded interval certification. Saved input decimals retain their original precision.',
            'The q=1 corrector is constructed independently from its closed normalized CDF, including its radial derivative; Cartesian Dirac differentiation checks signs and frame conversion.',
            'The four physical-q rows are analytic bound evaluations, not numerical inversion of the high-q angular operator.',
            'The sampled energies and weighted RMS are comparison scales, not a continuum or ultraviolet cutoff and not a unique assignment of loop energy to angular k.',
            'A local residual bound does not bound accumulated transport error without a propagation estimate and boundary data. No physical return/escape amplitude, vacuum stress, global state, curved-mouth matching, stable medium, or backward date passage is established.'])


if __name__ == '__main__':
    result = checks()
    Path(__file__).with_name('collar_dressing_bound_results.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in ('passed', 'tests', 'actual_collar', 'independent_q1')}, indent=2))
