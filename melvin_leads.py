"""Regular-axis transverse Melvin spectrum and oriented imaginary-energy leads.

The saved classical dihole approaches this environment only at infinity.
These leads do not supply a finite-cut return map or a supporting quantum state.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eig_banded
from scipy.integrate import quad
from scipy.optimize import brentq

from balanced_dihole import parameters


def coefficients(mu, beta, flux):
    if not np.all(np.isfinite([mu, beta, flux])) or beta <= 0 or flux <= 0:
        raise ValueError('Finite half-integer mu and positive beta, flux required.')
    if abs(mu-(round(mu-.5)+.5)) > 1e-12:
        raise ValueError('Smooth rotating-frame spinors have half-integer mu.')
    return beta*(2*mu-flux), beta**2*(mu-flux)


def radial_matrix(mu, spin, size, omega, beta, flux):
    """Exact polynomial Ritz matrix in the regular radial-oscillator domain.

    Three padding rows retain excursions of x^3 outside the retained basis.
    The critical inverse-square term at |mu-spin/2|=0 is treated analytically.
    """
    a, c = coefficients(mu, beta, flux)
    if spin not in (-1, 1) or not isinstance(size, int) or size < 4 or not math.isfinite(omega) or omega <= 0:
        raise ValueError('Require spin +/-1, integer size>=4 and finite omega>0.')
    nu = abs(mu-spin/2)
    j = np.arange(size+3, dtype=float)
    off = -np.sqrt((j[:-1]+1)*(j[:-1]+nu+1))
    x = np.diag(2*j+nu+1)+np.diag(off, 1)+np.diag(off, -1)
    x2 = x@x
    h = (np.diag(2*omega*(2*j+nu+1)+a*(2*mu+spin))
         +(a*a+c*(2*mu+3*spin)-omega*omega)*x/omega
         +2*a*c*x2/omega**2+c*c*(x2@x)/omega**3)
    return h[:size, :size]


def spectrum(mu, spin, size, omega, beta, flux, count=4):
    h = radial_matrix(mu, spin, size, omega, beta, flux)
    if not isinstance(count, int) or not 1 <= count <= size:
        raise ValueError('Invalid retained eigenvalue count.')
    band = np.zeros((4, size))
    for k in range(4):
        band[k, :size-k] = np.diag(h, -k)
    values, vectors = eig_banded(band, lower=True, select='i', select_range=(0, count-1))
    # Full differential-operator residual: polynomial H reaches only 3 more modes.
    larger = radial_matrix(mu, spin, size+3, omega, beta, flux)
    extended = np.pad(vectors, ((0, 3), (0, 0)))
    residual = np.linalg.norm(larger@extended-extended*values, axis=0)
    return values, vectors, residual


def zero_profile(mu, beta, flux):
    """Normalize exact upper zero mode; log scaling handles distant annuli."""
    a, c = coefficients(mu, beta, flux)
    if not .5 <= mu < flux:
        raise ValueError('Positive-flux zero modes require .5<=mu<flux.')
    peak2 = 2*mu/(math.sqrt(a*a-4*c*mu)-a)
    peak = math.sqrt(peak2)
    log_density = lambda r: 2*mu*math.log(r)+a*r*r+c*r**4/2
    offset = log_density(peak)
    width = 1/math.sqrt(-2*(-mu/peak2+a+3*c*peak2))
    # Integrate in widths about the peak so a narrow distant ring is not missed.
    lower = -peak/width
    density = lambda x: math.exp(log_density(peak+width*x)-offset) if x > lower else 0.
    norm, error = quad(density, lower, 16., epsabs=2e-12, epsrel=2e-12, points=[0.])
    def cdf(x):
        return quad(density, lower, x, epsabs=2e-12, epsrel=2e-12, points=[0.] if x>0 else None)[0]/norm
    quantiles = [peak+width*brentq(lambda x: cdf(x)-p, lower+1e-12, 16.) for p in (.05,.5,.95)]
    omitted = quad(density, 16., np.inf, epsabs=1e-13, epsrel=1e-11)[0]/norm
    return dict(mu=mu, peak=peak, width=width, radial_quantiles_05_50_95=quantiles,
                log_norm=offset+math.log(width*norm), quadrature_relative_error=error/norm,
                numerical_omitted_tail=omitted)


def outgoing_pair(kappa, mass, end=1):
    """Unit Euclidean-norm fixed mass-pair vector at axial +/- infinity.

    H=p*sigma3+mass*sigma1. p=end*i*sqrt(kappa^2+mass^2).
    This is resolvent data, not a real-frequency transmission probability.
    """
    if not np.all(np.isfinite([kappa, mass])) or kappa <= 0 or mass <= 0 or end not in (-1, 1):
        raise ValueError('Require finite kappa,mass>0 and end +/-1.')
    scale = max(kappa, mass)
    kk, mm = kappa/scale, mass/scale
    ratio = -1j*mm/(kk+math.hypot(kk, mm))
    vector = np.array([1., ratio] if end == 1 else [ratio, 1.], complex)
    return vector/np.linalg.norm(vector)


def gap_comparison():
    """Interval check of a global analytic oscillator lower comparison.

    Exact calibration Q=1/5,D=1,q=1. Negative-mu upper partner and positive
    mu>=21.5 lower partner potentials increase outward in mu (RESEARCH110).
    Thus 23 sectors cover every half-integer, including all positive partners
    of the 21 exact zeros. Fixed scales need no numerical optimization here.
    """
    import mpmath as mp
    iv, old = mp.iv, mp.iv.dps
    iv.dps = 70
    try:
        q, k = iv.mpf(1)/5, iv.mpf(1)/2
        mass = (iv.sqrt(k*k+4*q*k)-k)/2
        a = iv.sqrt(k*k-mass*mass)
        field = mass/((k+a)*(mass+k))
        beta, flux = field*field/4, 1/(q*field)
        p = parameters()
        assert abs(p['B']-float(field.mid)) < 1e-14 and p['Q'] == .2 and p['D'] == 1.
        assert flux.a > iv.mpf('20.5') and flux.b < iv.mpf('21.5')
        scales = ('0.3518','0.3462','0.3284','0.3106','0.2928','0.2750',
                  '0.2573','0.2396','0.2220','0.2045','0.1871','0.1698',
                  '0.1525','0.1354','0.1185','0.1017','0.0851','0.0687',
                  '0.0524','0.0364','0.0205','0.0049','0.3471')
        rows = []
        for j, scale in enumerate(scales):
            mu, spin = iv.mpf(2*j-1)/2, 1 if j == 0 else -1
            nu, omega = abs((2*j-1)/2-spin/2), iv.mpf(scale)
            aa, cc = beta*(2*mu-flux), beta*beta*(mu-flux)
            c0, c2 = aa*(2*mu+spin), aa*aa+cc*(2*mu+3*spin)
            c4, c6 = 2*aa*cc, cc*cc
            disc = 4*c4*c4-12*c6*(c2-omega*omega)
            candidates = [iv.mpf(0)]
            if disc.a >= 0:
                for sign in (-1, 1):
                    x = (-2*c4+sign*iv.sqrt(disc))/(6*c6)
                    # Including negative critical points only weakens the bound.
                    candidates.append((c2-omega*omega)*x+c4*x*x+c6*x*x*x)
            else:
                assert disc.b < 0, 'Critical-point discrimination not enclosed.'
            bound = 2*omega*(nu+1)+c0+min(value.a for value in candidates)
            assert bound.a > iv.mpf('0.60'), (j, bound)
            rows.append(dict(mu=(2*j-1)/2, spin=spin, oscillator_scale=scale,
                             squared_gap_lower_interval=str(bound)))
        return dict(passed=True, digits=70, uniform_squared_lower='0.60',
                    exact_calibration='Q=1/5,D=1,q=1', flux_interval=str(flux), rows=rows)
    finally:
        iv.dps = old


def run():
    p = parameters()
    beta, flux = p['B']**2/4, p['q']/(p['Q']*p['B'])
    zeros = [zero_profile(float(mu), beta, flux) for mu in np.arange(.5, flux, 1.)]
    rows = []
    for mu in np.arange(-.5, math.ceil(flux)+1, 1.):
        mu = float(mu)
        a, c = coefficients(mu, beta, flux)
        has_zero = .5 <= mu < flux
        omega = (max(.07, (mu+1)/zeros[int(mu-.5)]['peak']**2) if has_zero
                 else max(abs(a), math.sqrt(abs(c)), .07))
        controls = []
        for size, frequency in ((160, omega), (224, omega), (224, omega*1.1)):
            partners = []
            for spin in (1, -1):
                ev, _, residual = spectrum(mu, spin, size, frequency, beta, flux, 4)
                partners.append(dict(spin=spin, squared_eigenvalues=ev.tolist(),
                                     full_radial_residual_norms=residual.tolist()))
            controls.append(dict(size=size, omega=frequency, partners=partners))
        index = int(has_zero)
        upper, lower = [np.array(v['squared_eigenvalues']) for v in controls[1]['partners']]
        rows.append(dict(mu=mu, has_exact_zero=has_zero, omega=omega, controls=controls,
            first_three_squared_masses=lower[:3].tolist(),
            raw_upper_zero=float(upper[0]) if has_zero else None,
            paired_squared_mass_error=float(np.max(abs(upper[index:index+3]-lower[:3]))),
            size_change=max(float(np.max(abs(np.array(c0['squared_eigenvalues'])-c1['squared_eigenvalues'])))
                            for c0, c1 in zip(controls[0]['partners'], controls[1]['partners'])),
            scale_change=max(float(np.max(abs(np.array(c0['squared_eigenvalues'])-c1['squared_eigenvalues'])))
                             for c0, c1 in zip(controls[1]['partners'], controls[2]['partners']))))
    # At flux=mu, the cubic coefficient vanishes and the oscillator is exact.
    oscillator_errors = []
    for mu in (.5, 5.5):
        for spin in (-1, 1):
            omega = beta*mu
            ev, _, residual = spectrum(mu, spin, 16, omega, beta, mu)
            exact = 2*omega*(2*np.arange(4)+abs(mu-spin/2)+1)+omega*(2*mu+spin)
            oscillator_errors.append(float(np.max(abs(ev-exact))))
    lead_rows = []
    for kappa in (.001, .1, 1., 10.):
        for mass in (math.sqrt(rows[1]['first_three_squared_masses'][0]),
                     math.sqrt(rows[-2]['first_three_squared_masses'][0])):
            for end in (-1, 1):
                vector = outgoing_pair(kappa, mass, end)
                momentum = end*1j*math.hypot(kappa, mass)
                matrix = np.array([[momentum, mass], [mass, -momentum]])
                residual = np.linalg.norm(matrix@vector-1j*kappa*vector)
                lead_rows.append(dict(kappa=kappa, mass=mass, end=end,
                    vector_real=vector.real.tolist(), vector_imag=vector.imag.tolist(),
                    eigen_residual=float(residual),
                    normal_current=float(end*(abs(vector[0])**2-abs(vector[1])**2))))
    refined = [partner for row in rows for control in row['controls'][1:]
               for partner in control['partners']]
    comparison = gap_comparison()
    tests = dict(all_sector_interval_gap=comparison['passed'],
        exact_oscillator_limit=max(oscillator_errors)<1e-12,
        all_21_zero_modes=len(zeros)==21 and max(abs(row['raw_upper_zero']) for row in rows if row['has_exact_zero'])<1e-9,
        independent_partner_spectra=max(row['paired_squared_mass_error'] for row in rows)<1e-8,
        basis_size_control=max(row['size_change'] for row in rows)<1e-7,
        basis_scale_control=max(row['scale_change'] for row in rows)<1e-7,
        full_domain_residual=max(max(partner['full_radial_residual_norms']) for partner in refined)<2e-8,
        positive_partner_spectra=min(row['first_three_squared_masses'][0] for row in rows)>0,
        exact_zero_quadrature=max(row['quadrature_relative_error'] for row in zeros)<1e-9,
        zero_profile_tail=max(row['numerical_omitted_tail'] for row in zeros)<1e-20,
        oriented_resolvent_eigenvectors=max(row['eigen_residual'] for row in lead_rows)<1e-12,
        oriented_resolvent_current=min(row['normal_current'] for row in lead_rows)>0)
    result = dict(passed=all(tests.values()), tests=tests, parameters=p, beta=beta, flux=flux,
        sectors=rows, zero_modes=zeros, oriented_leads=lead_rows, gap_comparison=comparison,
        oscillator_limit_errors=oscillator_errors,
        smallest_computed_positive_mass=math.sqrt(min(row['first_three_squared_masses'][0] for row in rows)),
        scope=['Massless charged Dirac field on the asymptotic pure Melvin environment of the saved classical dihole.',
            'The raw squared zero Ritz values are retained; no negative numerical values are clipped into physical zeros.',
            'Ritz values are variational upper bounds in exact arithmetic. Full L2 residuals include leakage beyond the retained oscillator basis; floating arithmetic is not enclosed.',
            'A small residual alone does not identify an eigenvalue index. Exact zero modes, partner spectra, basis size and scale changes are separately checked.',
            'A separate analytic oscillator comparison, with 70-digit interval evaluation and pointwise sector monotonicity, gives squared transverse gap>0.60 for the exact Q=1/5,D=q=1 Melvin calibration.',
            'Zero channels have no transverse gap and travel axially. The other modes have transverse masses; no finite-cut throat connection, quantum stress, stability or time travel is shown.',
            'The imaginary-energy vectors specify outgoing resolvent data, not probabilities or real-frequency waves.',
            'The actual dihole has a gravitational long-range tail and requires propagation before these data reach a finite interface.'],
        sources=['https://arxiv.org/pdf/1508.07307 (exact Dirac operator22-25, not later small-radius spectrum)',
                 'https://arxiv.org/pdf/hep-th/9906160 (classical geometry and asymptotic environment)'],
        input_sha256=hashlib.sha256(Path('balanced_dihole.py').read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if not result['passed']:
        raise AssertionError(tests)
    Path(__file__).with_name('melvin_leads_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests,
        smallest_computed_positive_mass=result['smallest_computed_positive_mass'],
        maximum_size_change=max(row['size_change'] for row in rows),
        maximum_scale_change=max(row['scale_change'] for row in rows),
        maximum_full_residual=max(max(partner['full_radial_residual_norms']) for partner in refined),
        outermost_zero=zeros[-1]), indent=2))
    return result


if __name__ == '__main__':
    run()
