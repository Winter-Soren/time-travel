"""Full-angular proper-Dirac geometry on an auxiliary sewn collar buffer.

No angular-sector restriction or full-field smoothing value is asserted here.
The exact collar profile is constant for s <= 0, is the saved smooth opening
on [0,.05], and continues as R=s for s >= .05.  The auxiliary Dirichlet
region is (-.020,.070) x S2.  All lengths are in baseline D units.
"""
import math

import numpy as np

from dihole_channel_bounds import calibration, field_function, lower, upper


def _magnitude(value):
    return max(abs(lower(value)), abs(upper(value)))


def _geometry(radius, x, parameters, raw, sqrt):
    """Regularized log jets; neither a polar division nor a 1/R cancellation."""
    k, mass, ag, _, _ = parameters
    f, weight, flux = raw[:3]
    fR, _, fluxR = raw[3:6]
    fx, _, fluxx = raw[6:9]
    fRR, _, _ = raw[9:12]
    _, _, fluxRx = raw[12:15]
    fxx, _, _ = raw[15:18]
    distance = sqrt(radius*radius+4*k*k+4*k*radius*x)
    dR = (radius+2*k*x)/distance
    dx = 2*k*radius/distance
    dRR = 4*k*k*(1-x*x)/distance**3
    dxx = -(2*k*radius)**2/distance**3
    r = mass+(radius+distance)/2
    c = (distance-radius)/(2*k)
    rR, cR = (1+dR)/2, (dR-1)/(2*k)
    rx, cx = dx/2, dx/(2*k)
    sigma = r*r-ag*ag*c*c
    sigmaR = 2*r*rR-2*ag*ag*c*cR
    sigmax = 2*r*rx-2*ag*ag*c*cx
    sigmaRR = 2*rR*rR+r*dRR-2*ag*ag*cR*cR-ag*ag*c*dRR/k
    sigmaxx = 2*rx*rx+r*dxx-2*ag*ag*cx*cx-ag*ag*c*dxx/k
    ldR, ldx = dR/distance, dx/distance
    ldRR = dRR/distance-ldR*ldR
    ldxx = dxx/distance-ldx*ldx
    lsR, lsx = sigmaR/sigma, sigmax/sigma
    lsRR = sigmaRR/sigma-lsR*lsR
    lsxx = sigmaxx/sigma-lsx*lsx
    # n0=log(N/R), p0=log(P*R); L0=n0-p0=2log(d/sigma).
    n0R = (fR/f)/2+ldR-lsR
    n0RR = (fRR/f-(fR/f)**2)/2+ldRR-lsRR
    nx = (fx/f)/2+ldx-lsx
    nxx = (fxx/f-(fx/f)**2)/2+ldxx-lsxx
    L0R, L0RR = 2*(ldR-lsR), 2*(ldRR-lsRR)
    Lx, Lxx = 2*(ldx-lsx), 2*(ldxx-lsxx)
    return dict(N=sqrt(f/weight), P=sqrt(f*weight), NP=f, W=weight,
                flux=flux, flux_R=fluxR, flux_x=fluxx, flux_Rx=fluxRx,
                n0_R=n0R, n0_RR=n0RR, n_x=nx, n_xx=nxx,
                L0_R=L0R, L0_RR=L0RR, L_x=Lx, L_xx=Lxx)


def scalar_curvature(radius, x, rate_squared, second, data):
    """R_h of P^2(ds^2+R^2 dtheta^2)+R^2 sin^2(theta)/N^2 dphi^2."""
    nR, nx = data['n0_R'], data['n_x']
    LR, LRR, Lx, Lxx = (data[name] for name in ('L0_R', 'L0_RR', 'L_x', 'L_xx'))
    angular = -x*(2*nx+Lx)+(1-x*x)*(Lxx-nx*nx)
    radial = LRR+LR/radius-nR*nR
    return 2*((1+angular)/radius**2+rate_squared*radial+second*LR)/data['P']**2


def formula_controls():
    """Independent unregularized warped-product identity and exact simple metrics."""
    from mpmath.ctx_iv import MPIntervalContext
    iv = MPIntervalContext()
    iv.dps = 60
    parameters = [float(v.mid) for v in calibration(iv)]
    function = field_function(np)
    errors = []
    for radius, x, rate, second in ((.025, -.8, 0., 0.), (.037, .2, .4, 25.),
                                    (.05, .9, 1., 0.), (.069, -.95, 1., 0.)):
        raw = function(radius, x, *parameters)
        data = _geometry(radius, x, parameters, raw, np.sqrt)
        nR = 1/radius+data['n0_R']
        pR = -1/radius+data['n0_R']-data['L0_R']
        nRR = -1/radius**2+data['n0_RR']
        pRR = 1/radius**2+data['n0_RR']-data['L0_RR']
        nx, nxx = data['n_x'], data['n_xx']
        px, pxx = nx-data['L_x'], nxx-data['L_xx']
        sine2 = 1-x*x
        ntheta2 = -x*nx+sine2*nxx
        ptheta2 = -x*px+sine2*pxx
        ns, ps = rate*nR, rate*pR
        nss, pss = rate*rate*nRR+second*nR, rate*rate*pRR+second*pR
        base_scalar = -2*(second/radius+pss+rate*ps/radius+ptheta2/radius**2)/data['P']**2
        warped_laplacian = (second/radius+rate*rate/radius**2-nss
                            -3*rate*ns/radius+ns*ns
                            +(-1-ntheta2+2*x*nx+sine2*nx*nx)/radius**2)/data['P']**2
        direct = base_scalar-2*warped_laplacian
        regular = scalar_curvature(radius, x, rate*rate, second, data)
        errors.append(float(abs(direct-regular)/max(1., abs(direct), abs(regular))))
    radius, rate, second = .3, .4, .7
    flat = dict(P=1., n0_R=-1/radius, n_x=0., L0_R=-2/radius,
                L0_RR=2/radius**2, L_x=0., L_xx=0.)
    flat_error = abs(scalar_curvature(radius, .2, rate*rate, second, flat)
                     -(2*(1-rate*rate)/radius**2-4*second/radius))
    poles = []
    for radius in (.025, .047, .07):
        for x in (-1., 1.):
            raw = np.asarray(function(radius, x, *parameters), float)
            poles.append(float(max(abs(raw[0]-1), abs(raw[3]), abs(raw[5]),
                                   abs(raw[2]-(1-x)/2))))
    controls = dict(warped_product_relative_errors=errors,
                    flat_spherical_scalar_error=float(flat_error),
                    balanced_polar_value_derivative_errors=poles)
    controls['passed'] = bool(max(errors)<1e-10 and flat_error<1e-10 and max(poles)<1e-9)
    return controls


def buffer_bounds(iv, radial=24, angular=32):
    """Directed compact bounds; C(q)=C0+Cq*abs(q) is uniform in frequency."""
    if radial < 1 or angular < 1:
        raise ValueError('Positive interval cover dimensions required.')
    parameters = calibration(iv)
    function = field_function(iv)
    Rmin, Rmax = iv.mpf(1)/40, iv.mpf(7)/100
    r0 = Rmin
    rate, speed_squared, second = iv.mpf([0, 1]), iv.mpf([0, 1]), iv.mpf([0, 160])
    bounds = {name:[math.inf, -math.inf] for name in
              ('N', 'P', 'W', 'NP', 'R_h', 'volume_ratio', 'log_lapse_gradient_squared',
               'radial_volume_connection', 'radial_spin_connection', 'radial_laplacian_drift',
               'n_R', 'p_R')}
    rows = []
    arx = 0.
    for i in range(radial):
        radius = iv.mpf([Rmin+(Rmax-Rmin)*i/radial, Rmin+(Rmax-Rmin)*(i+1)/radial])
        for j in range(angular):
            x = iv.mpf([-1+iv.mpf(2)*j/angular, -1+iv.mpf(2)*(j+1)/angular])
            data = _geometry(radius, x, parameters, function(radius, x, *parameters), iv.sqrt)
            if data['N'].a <= 0 or data['P'].a <= 0:
                raise ArithmeticError('Cover failed to separate positive metric factors; refine it.')
            sine2 = 1-x*x
            nR, pR = 1/radius+data['n0_R'], -1/radius+data['n0_R']-data['L0_R']
            gradient = (speed_squared*nR**2+sine2*data['n_x']**2/radius**2)/data['P']**2
            values = dict(N=data['N'], P=data['P'], W=data['W'], NP=data['NP'],
                          R_h=scalar_curvature(radius, x, speed_squared, second, data),
                          volume_ratio=data['P']**2*radius**2/(data['N']*r0**2),
                          log_lapse_gradient_squared=gradient,
                          radial_volume_connection=rate*(data['n0_R']/2-data['L0_R']-1/(2*radius)),
                          radial_spin_connection=iv.sqrt(sine2)*(data['n_x']-data['L_x'])/(2*radius),
                          radial_laplacian_drift=rate*(1/radius-data['n0_R']), n_R=nR, p_R=pR)
            for name, value in values.items():
                bounds[name][0] = min(bounds[name][0], lower(value))
                bounds[name][1] = max(bounds[name][1], upper(value))
            arx = max(arx, _magnitude(data['flux_Rx']))
            rows.append((radius, data))
    magnetic = 0.
    for radius, data in rows:
        # a_R vanishes at both polar axes. Integrating a_Rx to the nearer
        # pole proves |a_R|/sin(theta)<=sup|a_Rx|, including endpoint limits.
        value = iv.sqrt(data['flux_x']**2/radius**4+iv.mpf(arx)**2/radius**2)/data['W']
        magnetic = max(magnetic, upper(value))
    nmin, nmax = map(iv.mpf, bounds['N'])
    pmin, pmax = map(iv.mpf, bounds['P'])
    absbound = lambda name:iv.mpf(max(abs(v) for v in bounds[name]))
    C0 = max(iv.mpf(0), -iv.mpf(bounds['R_h'][0]))/4+iv.mpf(bounds['log_lapse_gradient_squared'][1])/4
    connection = absbound('radial_volume_connection')+absbound('radial_spin_connection')
    cutoff_first, cutoff_second = iv.mpf(1600), iv.mpf(4800000)
    cutoff_gradient = cutoff_first/pmin
    cutoff_laplacian = (cutoff_second+cutoff_first*absbound('radial_laplacian_drift'))/pmin**2
    density_min, density_max = map(iv.mpf, bounds['volume_ratio'])
    inverse_metric_min = min(1/pmax**2, r0**2/(pmax**2*Rmax**2), r0**2*nmin**2/Rmax**2)
    comparison = density_min*inverse_metric_min/density_max
    controls = formula_controls()
    result = dict(domain_s_over_D=[-.020,.070], domain_length_over_D=.09,
                  central_length_over_D=.1,
                  source_s_over_D=[0.,.05], jet_buffer_s_over_D=[-.005,.055],
                  R_lower=lower(Rmin), R_upper=upper(Rmax),
                  R_lower_D=lower(Rmin), R_upper_D=upper(Rmax), reference_sphere_radius_over_D=.025,
                  radius_rate_abs_upper=1., radius_second_upper_D=160.,
                  radial_tiles=radial, angular_tiles=angular, boxes=radial*angular,
                  intervals=bounds, N_lower=lower(nmin), N_upper=upper(nmax),
                  P_lower=lower(pmin), P_upper=upper(pmax),
                  scalar_curvature_lower_D2=bounds['R_h'][0],
                  scalar_curvature_negative_upper_D2=max(0., -bounds['R_h'][0]),
                  magnetic_unit_flux_upper_D2=magnetic,
                  magnetic_norm_per_abs_q_upper_D2=magnetic,
                  log_lapse_gradient_squared_upper_D2=bounds['log_lapse_gradient_squared'][1],
                  polar_flux_Rx_upper_D=arx,
                  lichnerowicz_C0_upper_D2=upper(C0), lichnerowicz_Cq_upper_D2=magnetic,
                  canonical_radial_connection_upper_D=upper(connection),
                  radial_connection_upper_D=upper(connection),
                  radial_laplacian_drift_upper_D=upper(absbound('radial_laplacian_drift')),
                  cutoff_first_upper_D=1600., cutoff_second_upper_D2=4800000.,
                  cutoff_gradient_upper_D=upper(cutoff_gradient),
                  cutoff_laplacian_upper_D2=upper(cutoff_laplacian),
                  scalar_volume_ratio_lower=lower(density_min), scalar_volume_ratio_upper=upper(density_max),
                  scalar_laplacian_comparison_lower=lower(comparison),
                  controls=controls, geometry_controls=controls,
                  scope=['Full polar sphere and all angular states; no mu-band derivative estimate.',
                         'P_kappa >= covariant spatial Laplacian - C0 - abs(q)*Cq as a compact Dirichlet form.',
                         'Local Dirichlet comparison is not an upper bound for the global heat kernel.',
                         'No Hilbert-Schmidt smoothing norm, physical source, or continuum heat value is evaluated.'])
    if not result['controls']['passed']:
        raise ArithmeticError('Independent buffer formula control failed.')
    return result
