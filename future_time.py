"""Ideal finite-acceleration roundtrip, not a propulsion or survival model."""
import json
from pathlib import Path

import mpmath as mp
from scipy.integrate import quad


def proton_ring_checks():
    """Lorentz orbit plus a separate, orbit-averaged radiation-loss check."""
    import math
    import numpy as np
    from scipy.constants import c, e, epsilon_0, m_p, hbar
    from scipy.integrate import solve_ivp
    gamma, field, day = 10., 1., 86400.
    momentum = math.sqrt(gamma**2-1)  # momentum in units m*c
    beta = momentum/gamma
    radius = m_p*c*momentum/(e*field)
    omega = e*field/(gamma*m_p)
    period = 2*math.pi/omega
    acceleration = (beta*c)**2/radius
    power = e**2*gamma**4*acceleration**2/(6*math.pi*epsilon_0*c**3)
    kappa = e**4*field**2/(6*math.pi*epsilon_0*m_p**3*c**3)
    assert math.isclose(e*beta*c*field, gamma*m_p*acceleration, rel_tol=1e-14)
    assert math.isclose(power, m_p*c*c*kappa*(gamma**2-1), rel_tol=1e-14)

    # ponytail: ten turns check the ideal orbit; this is not a one-day lattice/stability simulation.
    def orbit_rhs(s, state):
        px, py = state[2:4]
        instantaneous_gamma = math.sqrt(1+px*px+py*py)
        factor = gamma/instantaneous_gamma
        return (factor*px/momentum, factor*py/momentum,
                factor*py, -factor*px, 1/instantaneous_gamma)
    times = np.linspace(0, 20*math.pi, 501)
    orbit = solve_ivp(orbit_rhs, (0, times[-1]), (1, 0, 0, -momentum, 0),
                     t_eval=times, method='DOP853', rtol=2e-12, atol=2e-13, max_step=.1)
    assert orbit.success
    position_error = float(np.max(np.abs(orbit.y[:2]-np.array([np.cos(times), -np.sin(times)]))))
    gamma_error = float(np.max(np.abs(np.sqrt(1+np.sum(orbit.y[2:4]**2, axis=0))-gamma)))
    clock_error = float(np.max(np.abs(orbit.y[4]-times/gamma)))
    assert max(position_error, gamma_error, clock_error)<1e-9

    with mp.workdps(60):
        k = mp.mpf(str(kappa)); u0 = mp.acoth(gamma)
        final_gamma = mp.coth(u0+k*day)
        proper_time = (mp.log(mp.cosh(u0+k*day))-mp.log(mp.cosh(u0)))/k
    loss = solve_ivp(lambda s, y: (-kappa*day*(y[0]**2-1), day/y[0]),
                     (0, 1), (gamma, 0), method='DOP853', rtol=2e-12, atol=2e-13)
    assert loss.success
    loss_errors = (abs(loss.y[0, -1]-float(final_gamma)),
                   abs(loss.y[1, -1]-float(proper_time)))
    assert loss_errors[0]<1e-11 and loss_errors[1]<1e-7
    assert float(final_gamma)<gamma and float(proper_time)>day/gamma
    final_radius = m_p*c*math.sqrt(float(final_gamma)**2-1)/(e*field)
    return dict(
        model='One proton, perpendicular ideal uniform magnetic field B=1T, nominal gamma10. Circular storage interval, not a complete injection/extraction trajectory',
        constants=dict(source='CODATA2022 via installed scipy.constants; verified against https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
                       c=c, elementary_charge=e, proton_mass_kg=m_p, vacuum_permittivity=epsilon_0),
        equations='R=gamma*m*v/(eB), omega=eB/(gamma*m), a_lab=v^2/R, a_proper=gamma^2*a_lab, d_tau=d_t/gamma; B does no work',
        nominal=dict(gamma=gamma, field_tesla=field, speed_over_c=beta,
            bending_radius_metres=radius, orbit_diameter_metres=2*radius,
            circumference_metres=2*math.pi*radius, revolution_frequency_hz=1/period,
            lab_acceleration_m_per_s2=acceleration, proper_acceleration_m_per_s2=gamma**2*acceleration,
            lab_storage_seconds=day, proton_proper_seconds=day/gamma,
            payload_kinetic_energy_GeV=(gamma-1)*m_p*c*c/(e*1e9),
            payload_kinetic_energy_joules=(gamma-1)*m_p*c*c),
        radiation=dict(formula='Free-space classical transverse Lienard power P=e^2*gamma^4*a_lab^2/(6*pi*epsilon0*c^3)',
            power_per_proton_watts=power, energy_per_turn_eV=power*period/e,
            energy_per_day_eV_if_gamma_held=power*day/e,
            effective_replacement_voltage_per_turn=power*period/e,
            fraction_of_kinetic_energy_per_day_if_gamma_held=power*day/((gamma-1)*m_p*c*c),
            quantum_recoil_parameter=gamma*beta*e*hbar*field/(m_p*m_p*c*c),
            compensation_scope='Constant gamma assumes ideal energy replacement, e.g. RF. Effective voltage is only energy/e per turn, not an RF hardware specification. Free-space power omits chamber shielding/edge fields and additional radiation from control fields'),
        no_energy_replacement=dict(
            model='Orbit-averaged leading classical damping in uniform B; slow spiral, not full radiation-reaction/lattice dynamics',
            equations='gamma_dot=-kappa*(gamma^2-1); gamma=coth(acoth(gamma0)+kappa*t); tau=log(cosh(acoth(gamma0)+kappa*t)/cosh(acoth(gamma0)))/kappa',
            kappa_per_second=kappa, gamma_after_one_day=float(final_gamma),
            proper_seconds_after_one_day=float(proper_time),
            extra_proper_seconds_over_constant_gamma=float(proper_time)-day/gamma,
            orbit_radius_change_metres=final_radius-radius,
            independent_ODE_gamma_error=loss_errors[0], independent_ODE_clock_error_seconds=loss_errors[1]),
        orbit_checks=dict(turns=10, maximum_position_error_in_radius_units=position_error,
            maximum_gamma_error=gamma_error, maximum_clock_error_in_inverse_omega_units=clock_error),
        scope='Bending radius is not total facility size. Uniform B does not focus vertical perturbations. No one-day confinement, magnet lattice, vacuum lifetime, injection, extraction, deceleration or measured proton clock is demonstrated. Both clocks accumulate time continuously; the proton is not frozen or sent to the past. A proton is not a grain or organic sample; extension needs charge-to-mass and internal-force tests',
        sources=['https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/xray_production/InvSynFracInt.html (magnetic rigidity and radiation scope)',
                 'https://arxiv.org/pdf/2107.09131 (Prat, equations8-9, classical radiation power and per-turn loss)',
                 'https://www.nature.com/articles/268301a0 (Bailey et al.,1977: circular-orbit time-dilation experiment using muons, not this one-day proton scenario)'])


def proton_focusing_checks(nominal):
    """Exact local Maxwell field, conservative aperture bound, and disturbed orbits."""
    import math
    import numpy as np
    import sympy as sp
    from scipy.integrate import solve_ivp
    r, z, n = sp.symbols('r z n', real=True)
    psi = z*(1+n/2)-n*(z*r*r-sp.Rational(2, 3)*z**3)/2
    br, bz = sp.diff(psi, r), sp.diff(psi, z)
    assert sp.simplify(sp.diff(r*br, r)/r+sp.diff(bz, z)) == 0
    assert sp.simplify(sp.diff(br, z)-sp.diff(bz, r)) == 0
    assert sp.simplify((-r*sp.diff(bz, r)/bz).subs({r:1, z:0})-n) == 0
    ap = r*(1+n/2)/2-n*r**3/8+n*r*z*z/2
    assert sp.simplify(sp.diff(r*ap, r)/r-bz) == 0
    assert sp.simplify(-sp.diff(ap, z)-br) == 0
    ell = -sp.Rational(1, 2)+n/8
    potential = (ell/r-ap)**2
    assert all(sp.simplify(sp.diff(potential, v).subs({r:1, z:0})) == 0 for v in (r, z))
    assert sp.simplify(sp.hessian(potential, (r, z)).subs({r:1, z:0})-sp.diag(2*(1-n), 2*n)) == sp.zeros(2)
    radius, gamma0 = nominal['bending_radius_metres'], nominal['gamma']
    momentum = math.sqrt(gamma0**2-1)
    index, kick = .4, 1e-5
    with mp.workdps(60):
        R, nn = mp.mpf(str(radius)), mp.mpf('.4')
        aperture = mp.mpf('.01')/R
        ell0 = -mp.mpf('.5')+nn/8
        def azimuthal_potential(rr, zz):
            return rr*(1+nn/2)/2-nn*rr**3/8+nn*rr*zz**2/2
        def V(rr, zz):
            return (ell0/rr-azimuthal_potential(rr, zz))**2
        # M=ell/r-Aphi is negative throughout the aperture, so M=0 supplies no missing extrema.
        assert ell0<0 and (1+nn/2)/2-nn*(1+aperture)**2/8>0
        candidates = [(rr, zz) for rr in (1-aperture, 1+aperture)
                      for zz in (-aperture, mp.mpf(0), aperture)]
        aa, bb = 3*nn/8, -((1+nn/2)/2+nn*aperture**2/2)
        for vv in ((-bb-mp.sqrt(bb*bb+4*aa*ell0))/(2*aa),
                   (-bb+mp.sqrt(bb*bb+4*aa*ell0))/(2*aa)):
            rr = mp.sqrt(vv)
            if 1-aperture <= rr <= 1+aperture:
                candidates.extend(((rr, -aperture), (rr, aperture)))
        barrier = min(V(rr, zz) for rr, zz in candidates)
        initial_level = V(1+mp.mpf('.001')/R, mp.mpf('.0005')/R)+2*mp.mpf(str(kick))**2
        assert 1<initial_level<barrier
        barrier_record = dict(aperture_half_width_metres=.01,
            boundary_minimum_above_reference=mp.nstr(barrier-1, 30),
            initial_level_above_reference=mp.nstr(initial_level-1, 30),
            positive_margin=mp.nstr(barrier-initial_level, 30),
            boundary_to_initial_excess_ratio=mp.nstr((barrier-1)/(initial_level-1), 20))

    cases = {}
    samples = np.linspace(0, 200*math.pi, 8001)
    for field_index in (0., index):
        def aphi(rr, zz):
            return rr*(1+field_index/2)/2-field_index*rr**3/8+field_index*rr*zz*zz/2
        rr, zz, canonical_l = 1+.001/radius, .0005/radius, -.5+field_index/8
        initial = np.array([rr, 0, zz, momentum*kick,
                            momentum*(canonical_l/rr-aphi(rr, zz)), momentum*kick, 0])
        initial_gamma = math.sqrt(1+np.dot(initial[3:6], initial[3:6]))
        def rhs(s, state):
            x, y, z, px, py, pz, clock = state
            factor = gamma0/math.sqrt(1+px*px+py*py+pz*pz)
            bx, by = -field_index*z*x, -field_index*z*y
            bz = 1+field_index/2-field_index*(x*x+y*y)/2+field_index*z*z
            return (factor*px/momentum, factor*py/momentum, factor*pz/momentum,
                    factor*(py*bz-pz*by), factor*(pz*bx-px*bz), factor*(px*by-py*bx),
                    factor/gamma0)
        result = solve_ivp(rhs, (0, samples[-1]), initial, t_eval=samples,
                           method='DOP853', rtol=2e-12, atol=2e-13, max_step=.1)
        assert result.success
        radii = np.hypot(result.y[0], result.y[1])
        vertical = result.y[2]
        energy = np.sum(result.y[3:6]**2, axis=0)/momentum**2
        angular = ((result.y[0]*result.y[4]-result.y[1]*result.y[3])/momentum
                   +radii*aphi(radii, vertical))
        energy_error = float(np.max(np.abs(energy-energy[0])))
        angular_error = float(np.max(np.abs(angular-canonical_l)))
        clock_error = float(np.max(np.abs(result.y[6]-samples/initial_gamma)))
        assert max(energy_error, angular_error, clock_error)<1e-9
        radial_mm, vertical_mm = (radii-1)*radius*1000, vertical*radius*1000
        if field_index:
            assert max(np.max(np.abs(radial_mm)), np.max(np.abs(vertical_mm)))<10
            assert energy_error<float(barrier-initial_level)*1e-3
        else:
            expected = zz+gamma0/initial_gamma*kick*samples
            assert np.max(np.abs(vertical-expected))<1e-11
            assert np.max(np.abs(vertical_mm))>10
        short = samples <= 20*math.pi
        cases[str(field_index)] = dict(field_index=field_index, reference_turns=100,
            gamma_initial=initial_gamma, maximum_radial_displacement_mm=float(np.max(np.abs(radial_mm))),
            maximum_vertical_displacement_mm=float(np.max(np.abs(vertical_mm))),
            energy_invariant_error=energy_error, canonical_angular_momentum_error=angular_error,
            clock_error_in_inverse_omega_units=clock_error,
            plot=dict(reference_turn=(samples[short]/(2*math.pi)).tolist(),
                      radial_mm=radial_mm[short].tolist(), vertical_mm=vertical_mm[short].tolist()))
    return dict(
        field='Dimensionless cylindrical r,z scaled by R: B/B0=(-n*z*r radial,0 azimuthal,1+n/2-n*r^2/2+n*z^2 axial); exact divB=curlB=0 locally, no electric field',
        scalar_potential='Psi/(B0*R)=z*(1+n/2)-n*(z*r^2-2*z^3/3)/2',
        vector_potential='Aphi/(B0*R)=r*(1+n/2)/2-n*r^3/8+n*r*z^2/2',
        nominal_field_index=index, radial_gradient_tesla_per_metre=-index*nominal['field_tesla']/radius,
        linear_tunes=dict(radial=math.sqrt(1-index), vertical=math.sqrt(index)),
        confinement='At fixed canonical ell=Lz/(e*B0*R^2)=-1/2+n/8, V=(ell/r-Aphi)^2 has Hessian diag(2*(1-n),2*n) at(1,0). Energy gives p_r^2/p0^2+p_z^2/p0^2+V=constant. A complete edge-extremum check yields a boundary level higher than the initial level, preventing escape from the1cm rectangle for this conservative model',
        aperture_bound=barrier_record,
        initial_condition='Radial offset+1mm, vertical offset+.5mm, p_r/p0=p_z/p0=1e-5. Azimuthal momentum is chosen to preserve the reference canonical angular momentum; do not renormalize to gamma10 afterward. Actual initial gamma is reported per case',
        cases=cases,
        scope='Ideal static local field and classical point proton. The energy barrier is an all-time conservative confinement argument for the specified invariants, not a proof for RF/radiation, field errors, collisions or finite matter. The polynomial field is a local target, not a global finite-energy magnet design. One hundred reference turns are about65.6microseconds, not a one-day storage test. Injection/extraction and a clock observable remain unresolved',
        source='Kim and Semertzidis, PhysRevAccelBeams25,024001(2022), https://arxiv.org/abs/2110.10003 (weak-focusing dynamics); the exact local field and aperture-level calculation here are our chosen benchmark')


def proton_slow_loss_checks(ring):
    """Adiabatic circular-orbit drift and leading linear LL transverse envelopes."""
    import math
    import numpy as np
    import sympy as sp
    from scipy.integrate import solve_ivp
    nominal = ring['nominal']
    radius, gamma0 = nominal['bending_radius_metres'], nominal['gamma']
    kappa = ring['no_energy_replacement']['kappa_per_second']
    omega0 = 2*math.pi*nominal['revolution_frequency_hz']
    day, index = 86400., .4
    rr, nn = sp.symbols('r n', positive=True)
    ff = 1+nn/2-nn*rr**2/2
    hh = rr*ff
    ell = -(1+nn/2)*rr**2/2+3*nn*rr**4/8
    assert sp.simplify(sp.diff(ell, rr)+rr*sp.diff(hh, rr)) == 0
    def quantities(r):
        f, fp = 1+index/2-index*r*r/2, -index*r
        h, hp = r*f, f+r*fp
        gamma = math.sqrt(1+(gamma0*gamma0-1)*h*h)
        gp = (gamma0*gamma0-1)*h*hp/gamma
        n = index*r*r/f
        np_ = index*(2*r*f-r*r*fp)/(f*f)
        drift = -kappa*gamma*f*f*h/hp
        omega = omega0*gamma0*f/gamma
        radial, vertical = omega*math.sqrt(1-n), omega*math.sqrt(n)
        bx = kappa*f*f*((1-n)/gamma-gamma*(1-2*n)/(1-n))
        bz = kappa*f*f*n/gamma
        logwrp, logwzp = fp/f-gp/gamma-np_/(2*(1-n)), fp/f-gp/gamma+np_/(2*n)
        return gamma, drift, n, radial, vertical, bx, bz, logwrp, logwzp
    def rhs(t, y):
        gamma, drift, n, wr, wz, bx, bz, dx, dz = quantities(1+y[0])
        return (day*drift, day/gamma, -.5*day*(bx+dx*drift), -.5*day*(bz+dz*drift))
    samples = np.linspace(0, 1, 101)
    result = solve_ivp(rhs, (0, 1), np.zeros(4), method='DOP853', rtol=2e-12,
                       atol=2e-14, t_eval=samples, max_step=.05)
    assert result.success
    rfinal = 1+result.y[0, -1]
    final = quantities(rfinal)
    with mp.workdps(60):
        kk, gg, ni = mp.mpf(str(kappa)), mp.mpf(str(gamma0)), mp.mpf('.4')
        def geom(r):
            f = 1+ni/2-ni*r*r/2
            h, hp = r*f, 1+ni/2-3*ni*r*r/2
            gamma = mp.sqrt(1+(gg*gg-1)*h*h)
            return f, h, hp, gamma
        def inverse_drift(r):
            f, h, hp, gamma = geom(r)
            return hp/(kk*gamma*f*f*h)
        rcheck = mp.findroot(lambda r: mp.quad(inverse_drift, [r, 1])-day, (.9999, 1))
        taucheck = mp.quad(lambda r: inverse_drift(r)/geom(r)[3], [rcheck, 1])
        edge_days = mp.quad(inverse_drift, [1-mp.mpf('.01')/mp.mpf(str(radius)), 1])/day
        radius_error = abs(mp.mpf(str(result.y[0, -1]))-(rcheck-1))*mp.mpf(str(radius))
        clock_error = abs(result.y[1, -1]-float(taucheck))
        assert radius_error<mp.mpf('1e-12') and clock_error<1e-8
    ratios = np.exp(result.y[2:, -1])
    initial_radial = math.hypot(.001, radius*1e-5/math.sqrt(1-index))
    initial_vertical = math.hypot(.0005, radius*1e-5/math.sqrt(index))
    radial_envelope = np.abs(result.y[0])*radius+initial_radial*np.exp(result.y[2])
    vertical_envelope = initial_vertical*np.exp(result.y[3])
    # In uniform B, the moving-orbit radial envelope is constant to this order:
    # b_x=-kappa*(gamma-1/gamma), d(log omega)/dt=+kappa*(gamma-1/gamma).
    gx, kap = sp.symbols('gamma kappa', positive=True)
    assert sp.simplify(kap*(1/gx-gx)+kap*(gx*gx-1)/gx) == 0
    ffinal = 1+index/2-index*rfinal*rfinal/2
    hfinal = rfinal*ffinal
    gp = (gamma0*gamma0-1)*hfinal*(1+index/2-3*index*rfinal*rfinal/2)/final[0]
    assert abs(gp*final[1]+kappa*ffinal**2*(final[0]**2-1))<1e-22
    assert 0<final[2]<1 and final[1]<0 and final[0]<gamma0
    return dict(
        model='Slow radiation-driven motion along nearby circular equilibria, plus leading linear Landau-Lifshitz/WKB transverse envelopes. No RF or full day-long particle tracking',
        equations='r=R/R0,f=1+n/2-n*r^2/2,h=r*f,gamma=sqrt(1+(gamma0^2-1)*h^2); gamma_dot=-kappa0*f^2*(gamma^2-1); r_dot=-kappa0*gamma*f^2*h/h_prime; tau_dot=1/gamma',
        transverse='N=n*r^2/f,Omega=e*B0*f/(m*gamma); omega_r=Omega*sqrt(1-N),omega_z=Omega*sqrt(N). b_r=kappa0*f^2*((1-N)/gamma-gamma*(1-2*N)/(1-N)); b_z=kappa0*f^2*N/gamma. A_i/A_i0=sqrt(omega_i0/omega_i)*exp(-integral b_i dt/2). Derived to first order in radiation/displacement with the LL field-gradient term retained',
        parameters=dict(gamma_initial=gamma0, field_index=index, lab_seconds=day, kappa_per_second=kappa),
        one_day=dict(radius_change_metres=float(result.y[0, -1]*radius), gamma=final[0],
                     proper_seconds=float(result.y[1, -1]), local_field_index=final[2],
                     radial_amplitude_ratio=float(ratios[0]), vertical_amplitude_ratio=float(ratios[1]),
                     estimated_max_radial_envelope_metres=float(np.max(radial_envelope)),
                     estimated_max_vertical_envelope_metres=float(np.max(vertical_envelope)),
                     circular_center_days_to_inward_one_cm_edge=float(edge_days),
                     drift_to_vertical_frequency_ratio=abs(final[1]/rfinal)/final[4]),
        verification=dict(independent_radius_quadrature_error_metres=float(radius_error),
                          independent_clock_quadrature_error_seconds=clock_error,
                          energy_chain_rule=True, canonical_angular_momentum_torque=True,
                          uniform_field_WKB_cancellation=True),
        plot=dict(lab_hours=(24*samples).tolist(), radius_shift_mm=(1000*radius*result.y[0]).tolist(),
                  radial_envelope_mm=(1000*radial_envelope).tolist(), vertical_envelope_mm=(1000*vertical_envelope).tolist()),
        scope='Table starts on the nominal circular gamma10 orbit; it is not the exact clock of the previous kicked orbit. Envelopes apply leading linear offsets/kicks around the slowly moving circle. Tiny timescale ratios support averaging but are not a certified accumulated-error bound. Radiation changes canonical angular momentum, so the previous exact fixed-invariant aperture proof no longer applies. No field errors, collisions, RF, injection/extraction, global magnet design or organic survival are established',
        sources=['https://arxiv.org/pdf/2107.09131 (Prat equations8-9: transverse radiation power)',
                 'https://arxiv.org/pdf/physics/9911027 (Spohn equation8: Landau-Lifshitz including field gradients; equation10: uniform magnetic field). The slow-orbit and envelope application here is our calculation'])


def proton_rf_checks(ring):
    """Localized effective RF voltage: mean balance, linear map, nonlinear averaged bucket."""
    import math
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.constants import c, e, m_p
    nominal = ring['nominal']
    with mp.workdps(65):
        gamma, index, harmonic = mp.mpf(10), mp.mpf('.4'), 1
        total_energy = gamma*mp.mpf(str(m_p))*mp.mpf(str(c))**2/mp.mpf(str(e))
        beta2 = 1-1/gamma**2
        loss = mp.mpf(str(ring['radiation']['energy_per_turn_eV']))
        voltage = mp.mpf('.001')  # eV gain for one proton; effective peak gap voltage1mV.
        alpha_c, eta = 1/(1-index), 1/(1-index)-1/gamma**2
        phase = mp.pi-mp.asin(loss/voltage)
        a = 2*mp.pi*harmonic*eta/(beta2*total_energy)
        b = voltage*mp.cos(phase)
        damping = loss/(beta2*total_energy)*(2+beta2-index/(1-index))
        assert loss<voltage and damping>0 and -4+2*damping<a*b<0
        assert a*voltage*mp.cos(mp.pi-phase)>0  # Rising branch fails the stability sign.
        matrix = mp.matrix([[1, a], [b, 1-damping+a*b]])
        assert abs(mp.det(matrix)-(1-damping))<mp.mpf('1e-60')
        trace, modulus = matrix[0, 0]+matrix[1, 1], mp.sqrt(1-damping)
        angle = mp.acos(trace/(2*modulus))
        revolution = mp.mpf(str(nominal['revolution_frequency_hz']))
        synchrotron_hz = angle*revolution/(2*mp.pi)
        lower_saddle = voltage*(-2*mp.cos(phase)+(mp.pi-2*phase)*mp.sin(phase))
        half_height = mp.sqrt(2*lower_saddle/a)
        dispersion_radius = mp.mpf(str(nominal['bending_radius_metres']))*alpha_c*half_height/(beta2*total_energy)
        # Differentiate the exact circular family independently of the quoted slip/loss formulas.
        def family(r):
            f = 1+index/2-index*r*r/2
            g = mp.sqrt(1+(gamma*gamma-1)*(r*f)**2)
            return g/gamma/f, loss*g/gamma*f*(g*g-1)/(gamma*gamma-1), g/gamma*total_energy
        energy_prime = mp.diff(lambda r: family(r)[2], 1)
        time_slope = mp.diff(lambda r: family(r)[0], 1)/energy_prime
        loss_slope = mp.diff(lambda r: family(r)[1], 1)/energy_prime
        assert abs(2*mp.pi*harmonic*time_slope-a)<mp.mpf('1e-60')
        assert abs(loss_slope-damping)<mp.mpf('1e-60')
        tune_scale, energy_scale = mp.sqrt(a*voltage), mp.sqrt(voltage/a)
        fixed_point = dict(effective_peak_voltage_volts=float(voltage), harmonic=harmonic,
            frequency_hz=float(revolution*harmonic), synchronous_phase_degrees=float(phase*180/mp.pi),
            mean_gain_and_loss_eV_per_turn=float(loss), time_based_slip_factor=float(eta),
            phase_slope_rad_per_eV=float(a), rf_slope_eV_per_rad=float(b),
            loss_energy_derivative_per_turn=mp.nstr(damping, 25),
            synchrotron_frequency_hz=float(synchrotron_hz),
            longitudinal_amplitude_damping_years=float(-1/(revolution*mp.log(modulus))/(365.25*86400)),
            one_day_linear_damping_multiplier=mp.nstr(mp.exp(revolution*86400*mp.log(modulus)), 25),
            averaged_bucket_half_height_eV=float(half_height),
            bucket_linear_dispersion_metres=float(dispersion_radius))
    phase, rate, escale, damp = map(float, (phase, tune_scale, energy_scale, damping))
    initial = np.array([.2, 100/escale])  # phase relative to synchronous arrival; energy deviation stored separately.
    def rhs(s, y):
        return (y[1], 2*math.sin(y[0]/2)*math.cos(phase+y[0]/2))
    # One real-parameter synchrotron cycle: a discrete thin-gap map versus its conservative average.
    turns = round(2*math.pi/rate)
    phase_error, energy_error = initial
    discrete_samples = [(0, phase_error, energy_error)]
    stride = max(1, turns//400)
    for turn in range(turns):
        phase_error += rate*energy_error
        energy_error += rate*2*math.sin(phase_error/2)*math.cos(phase+phase_error/2)
        if (turn+1)%stride == 0 or turn+1 == turns:
            discrete_samples.append((turn+1, phase_error, energy_error))
    end = turns*rate
    solution = solve_ivp(rhs, (0, end), initial, dense_output=True, method='DOP853',
                         rtol=2e-12, atol=2e-14, max_step=.025)
    assert solution.success
    discrete = np.array(discrete_samples)
    map_error = float(np.max(np.abs(discrete[:, 1:].T-solution.sol(discrete[:, 0]*rate))))
    assert map_error<2e-6
    samples = np.linspace(0, end, 401); values = solution.sol(samples)
    hamiltonian = (.5*values[1]**2+np.cos(phase+values[0])-math.cos(phase)
                   +values[0]*math.sin(phase))
    invariant_error = float(np.max(np.abs(hamiltonian-hamiltonian[0])))
    assert invariant_error<1e-10 and 0<hamiltonian[0]<float(lower_saddle/voltage)
    return dict(
        model='Localized thin RF gap plus existing circular-family radiation mean; linear phase-slip/loss expansion and nonlinear sinusoidal averaged bucket. Effective voltage includes transit-time factor; no cavity Maxwell solution or full transverse/RF transport',
        sign_convention='Positive phase means later laboratory arrival. Voltage is positive along clockwise proton motion, gain=e*Vhat*sin(phi). eta_T=dlnT/dlnp=alpha_c-1/gamma^2>0. This is opposite the frequency-based slip sign used in Tecker. Stable synchronous phase is pi-asin(U0/eVhat)',
        map='phi_next=phi+2*pi*h*(T(E)/T0-1); deltaE_next=deltaE-U(E)+e*Vhat*sin(phi_next). Store deltaE separately; U0 is about one double-precision spacing of the9.38GeV total energy',
        linear_stability='M=[[1,a],[b,1-D+a*b]], a=2*pi*h*eta_T/(beta^2*E0),b=eVhat*cos(phi_s),D=U_prime(E0). detM=1-D; Schur conditionsD>0 and -4+2D<a*b<0. Verified with65-digit arithmetic',
        radiation='U(r)=U0*(gamma(r)/gamma0)*f(r)*(gamma(r)^2-1)/(gamma0^2-1); exact circular-family derivative checked independently',
        fixed_point=fixed_point,
        averaged_hamiltonian='H=a*deltaE^2/2+eVhat*(cos(phi)-cos(phi_s)+(phi-phi_s)*sin(phi_s)). With energy-dependent loss,H_prime=-a*deltaE*(U(E0+deltaE)-U0)<=0 whereU increases. H is not an exact discrete-map invariant',
        numerical_check=dict(discrete_turns=turns, compared_samples=len(discrete_samples),
            conservative_map_versus_averaged_sampled_state_error=map_error,
            error_coordinates=f'Maximum absolute difference over phase error in radians and energy deviation divided by {escale:.15g} eV',
            averaged_dimensionless_H_error=invariant_error,
            scope='Discrete and ODE comparison uses linear phase slip and the nonlinear RF sine, not full nonlinear T(E); it removes energy-dependent damping to isolate averaging error. Synchronous mean loss is still balanced. Tiny radiative damping is checked separately in the high-precision linear matrix, not claimed resolved by this floating-point turn loop'),
        plot=dict(seconds=(samples/rate/nominal['revolution_frequency_hz']).tolist(),
                  phase_error_degrees=(values[0]*180/math.pi).tolist(), energy_deviation_eV=(escale*values[1]).tolist()),
        scope='Mean energy replacement and longitudinal focusing in an ideal effective-gap model. This does not prove combined transverse confinement, control-noise tolerance, stochastic-emission behavior, field geometry, injection/extraction, one-day apparatus operation or finite-matter survival. The energy-offset bucket is not an injection system',
        sources=['https://arxiv.org/pdf/1601.04901 (Tecker: synchrotron phase stability and averaged Hamiltonian)',
                 'https://arxiv.org/pdf/1601.05230 (Jensen: RF cavity fields and effective accelerating voltage)'])


def proton_kicker_checks(ring):
    """Finite static vacuum fringe scattering and outward aperture departure."""
    import cmath
    import math
    import numpy as np
    import sympy as sp
    from scipy.integrate import solve_ivp
    nominal = ring['nominal']
    radius, gamma0 = nominal['bending_radius_metres'], nominal['gamma']
    field0, omega0 = nominal['field_tesla'], 2*math.pi*nominal['revolution_frequency_hz']
    momentum, index = math.sqrt(gamma0*gamma0-1), .4
    aperture, target, length, fringe = .01, .012, 1., .05
    tune, wave_number = math.sqrt(1-index), math.sqrt(1-index)/radius

    # One analytic tanh end satisfies both vacuum magnetostatic identities.
    u, v = sp.symbols('u v', real=True)
    denominator = sp.cosh(2*u)+sp.cos(2*v)
    bs, bz = sp.sin(2*v)/denominator, sp.sinh(2*u)/denominator
    assert sp.simplify(sp.diff(bs, u)+sp.diff(bz, v)) == 0
    assert sp.simplify(sp.diff(bs, v)-sp.diff(bz, u)) == 0
    def profile(s, z=0.):
        w = complex(s, z)
        return .5*(cmath.tanh((w+length/2)/fringe)-cmath.tanh((w-length/2)/fringe))
    for ss, zz in ((-.52, .005), (.48, -.01), (0., .01)):
        def end(x):
            den = math.cosh(2*x/fringe)+math.cos(2*zz/fringe)
            return complex(math.sinh(2*x/fringe)/den, math.sin(2*zz/fringe)/den)
        assert abs(profile(ss, zz)-.5*(end(ss+length/2)-end(ss-length/2)))<1e-13
    assert .01<math.pi*fringe/2
    extent = length/2+8*fringe
    edge = math.asin(extent/radius)
    start, end = -edge, math.pi/(2*tune)+.05
    # On the nominal circular orbit the straight magnet coordinate is s=R*sin(theta).
    circular_factor = radius/length*quad(
        lambda theta: profile(radius*math.sin(theta)).real*math.cos(tune*theta),
        -edge, edge, epsabs=1e-13, epsrel=1e-12, points=[-.5/radius, .5/radius])[0]
    straight_factor = (math.sin(wave_number*length/2)/(wave_number*length/2)
                       *(math.pi*fringe*wave_number/2)/math.sinh(math.pi*fringe*wave_number/2))
    integrated_profile = quad(lambda s: profile(s).real, -extent, extent,
                              epsabs=1e-13, epsrel=1e-12, points=[-.5, .5])[0]
    tail_fraction = fringe/length*math.log1p(math.exp(-16))
    assert abs(integrated_profile/length-(1-tail_fraction))<1e-13
    peak = -field0*tune*target/(length*circular_factor)  # Negative vertical field kicks clockwise protons outward.
    samples = np.unique(np.concatenate((np.linspace(start, edge, 161), np.linspace(edge, end, 481))))
    post = samples>=edge

    def trajectory(strength, refinement=1.):
        def rhs(theta, state):
            x, y, z, px, py, pz, clock, linear_x, linear_dx = state
            factor = gamma0/math.sqrt(1+px*px+py*py+pz*pz)
            extra = strength*peak/field0*profile(-radius*y, radius*z)
            bx, by = -index*z*x, -index*z*y-extra.imag
            bz = 1+index/2-index*(x*x+y*y)/2+index*z*z+extra.real
            forcing = -radius*strength*peak/field0*profile(radius*math.sin(theta)).real
            return (factor*px/momentum, factor*py/momentum, factor*pz/momentum,
                    factor*(py*bz-pz*by), factor*(pz*bx-px*bz), factor*(px*by-py*bx),
                    factor/gamma0, linear_dx, -tune*tune*linear_x+forcing)
        def departure(theta, state):
            return radius*(math.hypot(state[0], state[1])-1)-aperture
        departure.direction = 1
        initial = (math.cos(start), -math.sin(start), 0.,
                   -momentum*math.sin(start), -momentum*math.cos(start), 0., 0., 0., 0.)
        first = solve_ivp(rhs, (start, edge), initial, dense_output=True, method='DOP853',
                          rtol=2e-12*refinement**2, atol=2e-13*refinement**2,
                          max_step=fringe/(4*radius)*refinement)
        assert first.success
        second = solve_ivp(rhs, (edge, end), first.y[:, -1], dense_output=True, events=departure,
                           method='DOP853', rtol=2e-12*refinement**2, atol=2e-13*refinement**2,
                           max_step=.02*refinement)
        assert second.success
        values = np.empty((9, len(samples)))
        before = samples<=edge
        values[:, before], values[:, ~before] = first.sol(samples[before]), second.sol(samples[~before])
        radial = radius*(np.hypot(values[0], values[1])-1)
        gamma_error = float(np.max(np.abs(np.sqrt(1+np.sum(values[3:6]**2, axis=0))-gamma0)))
        clock_error = float(np.max(np.abs(values[6]-(samples-start)/gamma0))/omega0)
        thin = np.where(samples>=0, -strength*peak*length/(field0*tune)*np.sin(tune*samples), 0.)
        finite_error = float(np.max(np.abs(radial[post]-values[7, post])))
        thin_error = float(np.max(np.abs(radial[post]-thin[post])))
        assert gamma_error<1e-9 and clock_error<1e-15 and np.max(np.abs(values[2]))<1e-14
        crossing = None
        if len(second.t_events[0]):
            tt, st = second.t_events[0][0], second.y_events[0][0]
            pr = st[0]*st[3]+st[1]*st[4]
            ps = st[1]*st[3]-st[0]*st[4]
            crossing = dict(lab_seconds=float((tt-start)/omega0),
                            proton_proper_seconds=float(st[6]/omega0),
                            outward_angle_radians=math.atan2(pr, ps),
                            radial_displacement_metres=float(radius*(math.hypot(st[0], st[1])-1)),
                            position_metres=(radius*st[:3]).tolist(), momentum_over_mc=st[3:6].tolist())
            assert crossing['outward_angle_radians']>0
        record = dict(strength_fraction=strength, sampled_maximum_radial_metres=float(np.max(radial)),
                      maximum_gamma_error=gamma_error, clock_error_seconds=clock_error,
                      post_gap_finite_linear_position_error_metres=finite_error,
                      post_gap_thin_position_error_metres=thin_error,
                      relative_finite_linear_error=finite_error/(strength*target) if strength else None,
                      departure=crossing)
        return record, radial, values[7], thin

    main, radial, linear, thin = trajectory(1.)
    weak, weak_radial, weak_linear, weak_thin = trajectory(.01)
    control, control_radial, _, _ = trajectory(0.)
    fine, fine_radial, _, _ = trajectory(1., .5)
    assert main['departure'] is not None and weak['departure'] is None
    assert control['departure'] is None and np.max(np.abs(control_radial))<1e-9
    assert weak['relative_finite_linear_error']<main['relative_finite_linear_error']/5
    assert main['relative_finite_linear_error']<.01
    refinement_error = float(np.max(np.abs(radial-fine_radial)))
    crossing_error = abs(main['departure']['lab_seconds']-fine['departure']['lab_seconds'])
    assert refinement_error<1e-9 and crossing_error<1e-13
    linear_quadrature = target*np.sin(tune*samples)
    quadrature_error = float(np.max(np.abs(linear[post]-linear_quadrature[post])))
    assert quadrature_error<1e-9
    field_samples = np.linspace(-extent, extent, 241)
    return dict(
        model='Finite static vacuum-fringe kicker added to the existing n=.4 magnetic ring; a partial-turn Lorentz scattering check from a nominal gamma10 orbit',
        field='Local s=-Y,Z=z in metres, w=s+i*z; b(w)=Bk/2*(tanh((w+ell/2)/a)-tanh((w-ell/2)/a)); global deltaB=(0,-Im(b),Re(b)). Cauchy-Riemann identities give divB=curlB=0 for abs(z)<pi*a/2. Electric field is zero during this static flat-top calculation',
        parameters=dict(nominal_gamma=gamma0, field_index=index, aperture_departure_metres=aperture,
                        linear_target_amplitude_metres=target, magnet_length_parameter_metres=length,
                        fringe_scale_metres=fringe, vertical_analytic_half_width_metres=math.pi*fringe/2,
                        vertical_field_coefficient_tesla=peak, integrated_vertical_field_tesla_metres=peak*length,
                        initial_nominal_angle_radians=start, final_nominal_angle_radians=end,
                        elapsed_lab_seconds=(end-start)/omega0),
        thin_prediction=dict(minimum_angle_to_reach_one_cm_radians=tune*aperture/radius,
                             minimum_field_integral_tesla_metres=field0*tune*aperture,
                             quarter_betatron_distance_metres=math.pi*radius/(2*tune),
                             qualification='The minimum only grazes1cm at its turning point. The chosen12mm amplitude gives an outward crossing; physical septum/chamber clearance is unspecified'),
        finite_fringe=dict(straight_reference_response_factor=straight_factor,
                           circular_reference_response_factor=circular_factor,
                           entrance_and_exit_integral_tail_fraction=tail_fraction,
                           linear_ODE_vs_circular_quadrature_error_metres=quadrature_error),
        cases=dict(nominal_strength=main, one_percent_strength=weak, no_kick=control, refined_nominal=fine),
        verification=dict(vacuum_Maxwell_identities=True, actual_complex_field_matches_real_formula=True,
                          sampled_main_refinement_difference_metres=refinement_error,
                          departure_time_refinement_difference_seconds=crossing_error),
        plot=dict(lab_microseconds=((samples-start)/omega0*1e6).tolist(),
                  radial_mm=(1000*radial).tolist(), finite_linear_radial_mm=(1000*linear).tolist(),
                  thin_radial_mm=(1000*thin).tolist(), no_kick_radial_mm=(1000*control_radial).tolist(),
                  weak_radial_mm=(1000*weak_radial).tolist(), weak_finite_linear_radial_mm=(1000*weak_linear).tolist(),
                  field_s_metres=field_samples.tolist(),
                  midplane_vertical_field_millitesla=[1000*peak*profile(s).real for s in field_samples],
                  longitudinal_field_at_z5mm_millitesla=[1000*peak*profile(s, .005).imag for s in field_samples]),
        scope='This is aperture departure, not completed extraction. The original polynomial ring field has no matched exterior or septum; the kicker is a local analytic field, not a finite coil/return-current design. A static flat-top scattering calculation does not model switch-on/off, induced electric fields, injection history, or a complete storage pulse. RF, radiation, stochastic recoil and existing kicked-orbit offsets are absent. No combined long-term confinement, field-free exit, deceleration, capture or reunion is established. Analytic fringe tails are not manually cut off; the finite entrance omits a quantified small prior tail',
        sources=['https://arxiv.org/pdf/1103.1583 (Barnes et al., sections2-3: kicker timing and septum roles)',
                 'https://arxiv.org/pdf/1103.1062 (Barnes et al., section4, equation2: magnetic field-integral deflection)',
                 'https://arxiv.org/pdf/physics/0302069 (Papaphilippou et al., sectionII, equations1-2: vacuum scalar potential and required fringe components). The tanh field and ring application here are our benchmark'])


def proton_outgoing_candidate_checks(ring, kicker):
    """One harmonic transition field seeded by, but not globally joined to, the kicker state."""
    import cmath
    import math
    import numpy as np
    import sympy as sp
    from scipy.integrate import solve_ivp
    from scipy.special import erf
    departure = kicker['cases']['nominal_strength']['departure']
    nominal, constants = ring['nominal'], ring['constants']
    radius, field0 = nominal['bending_radius_metres'], nominal['field_tesla']
    c, charge, mass = constants['c'], constants['elementary_charge'], constants['proton_mass_kg']
    origin = np.array(departure['position_metres'])
    initial_p = np.array(departure['momentum_over_mc'])
    entry_radius = math.hypot(origin[0], origin[1])
    radial = np.array([origin[0]/entry_radius, origin[1]/entry_radius, 0.])
    tangent = np.array([origin[1]/entry_radius, -origin[0]/entry_radius, 0.])
    basis = np.column_stack((tangent, radial, np.array([0., 0., 1.])))
    assert np.max(np.abs(basis.T@basis-np.eye(3)))<1e-14 and abs(np.linalg.det(basis)-1)<1e-14
    local_p = basis.T@initial_p
    index, scale, checkpoint, half_width, half_height = .4, 1., 5., .5, .01
    entry_field = field0*(1+index/2-index*entry_radius**2/(2*radius**2))
    gradient = -index*field0*entry_radius/radius**2
    qmc = charge/(mass*c)

    # The potential is harmonic because it is linear in x and analytic in s+i*z.
    ss, xx, zz, aa, bd, gg = sp.symbols('s x z a Bd g', real=True)
    wp, wm = (ss+sp.I*zz)/aa, (ss-sp.I*zz)/aa
    ap = aa*sp.sqrt(sp.pi)/2*sp.erf(wp)
    am = aa*sp.sqrt(sp.pi)/2*sp.erf(wm)
    potential = (bd+gg*xx)*(ap-am)/(2*sp.I)
    cp, cm = sp.exp(-wp**2), sp.exp(-wm**2)
    symbolic_field = sp.Matrix([(bd+gg*xx)*(cp-cm)/(2*sp.I),
                                gg*(ap-am)/(2*sp.I), (bd+gg*xx)*(cp+cm)/2])
    coordinates = (ss, xx, zz)  # Right-handed order: clockwise tangent, outward, vertical.
    assert all(sp.simplify(sp.diff(potential, axis)-symbolic_field[i]) == 0
               for i, axis in enumerate(coordinates))
    jacobian = symbolic_field.jacobian(coordinates)
    assert sp.simplify(sp.trace(jacobian)) == 0
    assert sp.simplify(jacobian-jacobian.T) == sp.zeros(3)
    entry = {ss:0, xx:0, zz:0}
    expected_jet = sp.Matrix([[0, 0, 0], [0, 0, gg], [0, gg, 0]])
    assert sp.simplify(jacobian.subs(entry)-expected_jet) == sp.zeros(3)
    assert symbolic_field.subs(entry) == sp.Matrix([0, 0, bd])
    rr, rd, nn, bb = sp.symbols('R rd n B0', positive=True)
    original_field = sp.Matrix([-nn*bb*zz*ss/rr**2, -nn*bb*zz*(rd+xx)/rr**2,
                                 bb*(1+nn/2-nn*((rd+xx)**2+ss**2)/(2*rr**2)+nn*zz**2/rr**2)])
    assert sp.simplify(original_field.jacobian(coordinates).subs(entry)
                       -expected_jet.subs(gg, -nn*bb*rd/rr**2)) == sp.zeros(3)
    assert sp.simplify(original_field.subs(entry)-sp.Matrix(
        [0, 0, bb*(1+nn/2-nn*rd**2/(2*rr**2))])) == sp.zeros(3, 1)

    def field(s, x, z, original=False):
        if original:
            return np.array([-index*field0*z*s/radius**2,
                             -index*field0*z*(entry_radius+x)/radius**2,
                             field0*(1+index/2-index*((entry_radius+x)**2+s*s)/(2*radius**2)
                                     +index*z*z/radius**2)])
        w = complex(s, z)/scale
        analytic_c = cmath.exp(-w*w)
        analytic_a = scale*math.sqrt(math.pi)/2*erf(w)
        return np.array([(entry_field+gradient*x)*analytic_c.imag,
                         gradient*analytic_a.imag, (entry_field+gradient*x)*analytic_c.real])
    assert np.max(np.abs(field(0, 0, 0)-field(0, 0, 0, True)))<1e-15
    with mp.workdps(45):
        field_formula_error = 0.
        for s, x, z in ((0., .01, .005), (.7, -.03, -.01), (4., -.1, .008)):
            w = mp.mpc(s, z)/scale
            cc, av = mp.exp(-w*w), scale*mp.sqrt(mp.pi)/2*mp.erf(w)
            expected = np.array([float((entry_field+gradient*x)*cc.imag),
                                 float(gradient*av.imag), float((entry_field+gradient*x)*cc.real)])
            field_formula_error = max(field_formula_error, float(np.max(np.abs(field(s, x, z)-expected))))
        assert field_formula_error<1e-13
    samples = np.linspace(0, checkpoint, 401)

    def trajectory(original=False, refinement=1., vertical_offset=0.):
        # Reparameterize Lorentz motion by increasing s: dt/ds=gamma/(c*p_s), d_tau/ds=1/(c*p_s).
        initial = np.r_[0., vertical_offset, local_p, 0., 0.]
        gamma_initial = math.sqrt(1+np.dot(local_p, local_p))
        def rhs(s, state):
            x, z, ps, px, pz, ct, cclock = state
            p = np.array([ps, px, pz])
            gamma = math.sqrt(1+np.dot(p, p))
            force = qmc*np.cross(p, field(s, x, z, original))/ps
            return (px/ps, pz/ps, *force, gamma/ps, 1/ps)
        def side(s, state):
            return half_width-abs(state[0])
        def height(s, state):
            return half_height-abs(state[1])
        def forward(s, state):
            return state[2]-1.
        for event in (side, height, forward):
            event.terminal, event.direction = True, -1
        solution = solve_ivp(rhs, (0, checkpoint), initial, dense_output=True, method='DOP853',
                             events=(side, height, forward), max_step=.025*refinement,
                             rtol=2e-12*refinement**2, atol=2e-13*refinement**2)
        assert solution.success and solution.t[-1] == checkpoint
        assert all(len(events) == 0 for events in solution.t_events)
        values = solution.sol(samples)
        gamma_values = np.sqrt(1+np.sum(values[2:5]**2, axis=0))
        gamma_error = float(np.max(np.abs(gamma_values-gamma_initial)))
        clock_error = float(np.max(np.abs(values[6]-values[5]/gamma_initial))/c)
        assert gamma_error<1e-10 and clock_error<1e-18
        locations = origin[:, None]+basis@np.vstack((samples, values[:2]))
        radial_excess = np.hypot(locations[0], locations[1])-radius
        field_norm = np.array([np.linalg.norm(field(s, x, z, original))
                               for s, x, z in zip(samples, values[0], values[1])])
        record = dict(lab_seconds=float(values[5, -1]/c), proton_proper_seconds=float(values[6, -1]/c),
                      exit_position_metres=locations[:, -1].tolist(),
                      exit_momentum_over_mc=(basis@values[2:5, -1]).tolist(),
                      exit_local_sxz_metres=[checkpoint, float(values[0, -1]), float(values[1, -1])],
                      exit_radius_above_nominal_metres=float(radial_excess[-1]),
                      exit_field_norm_tesla=float(field_norm[-1]),
                      sampled_max_abs_x_metres=float(np.max(np.abs(values[0]))),
                      sampled_max_abs_z_metres=float(np.max(np.abs(values[1]))),
                      minimum_forward_momentum_over_mc=float(np.min(values[2])),
                      initial_gamma=gamma_initial, maximum_gamma_error=gamma_error,
                      proper_clock_identity_error_seconds=clock_error)
        return record, values, radial_excess, field_norm

    main, values, radial_excess, norms = trajectory()
    control, old_values, old_radial, old_norms = trajectory(original=True)
    fine, fine_values, _, _ = trajectory(refinement=.5)
    off_axis, _, _, _ = trajectory(vertical_offset=.001)
    assert main['exit_field_norm_tesla']<1e-9 and control['exit_field_norm_tesla']>.9*field0
    assert main['exit_radius_above_nominal_metres']>.2 and control['exit_radius_above_nominal_metres']<.02
    refinement_error = float(np.max(np.abs(values[:2]-fine_values[:2])))
    time_error = abs(main['lab_seconds']-fine['lab_seconds'])
    assert refinement_error<1e-9 and time_error<1e-16
    return dict(
        model='Candidate local outgoing transition field, using one harmonic potential throughout its finite test region and the saved kicker departure as initial data. This is not a global joining of the old and new field maps',
        potential='Right-handed local(s,x,z): w=(s+i*z)/a, A=a*sqrt(pi)/2*erf(w), C=exp(-w^2), Psi=(Bd+g*x)*Im(A); B=( (Bd+g*x)*Im(C), g*Im(A), (Bd+g*x)*Re(C) ). E=0. No field is switched or manually cut off in this calculation',
        parameters=dict(fringe_scale_metres=scale, checkpoint_s_metres=checkpoint,
                        entry_vertical_field_tesla=entry_field, entry_radial_gradient_tesla_per_metre=gradient,
                        valid_s_interval_metres=[-.05, 5.05], valid_x_interval_metres=[-half_width, half_width],
                        valid_z_interval_metres=[-half_height, half_height],
                        origin_global_metres=origin.tolist(), local_to_global_basis=basis.tolist()),
        matching=dict(field_and_first_spatial_derivative_at_entry=True,
                      original_d2Bz_ds2_tesla_per_metre2=-index*field0/radius**2,
                      candidate_d2Bz_ds2_tesla_per_metre2=-2*entry_field/scale**2,
                      qualification='Matching B and its first derivative at one event is not equality on an open region. Harmonic uniqueness forbids an exactly unchanged polynomial field on an open current-free region from acquiring a different continuation in the same connected vacuum. A physical insertion therefore requires a redesigned common field/source geometry and rechecking the preceding orbit'),
        cases=dict(candidate=main, unchanged_ring_control=control, refined_candidate=fine,
                   candidate_with_one_mm_vertical_offset=off_axis),
        verification=dict(harmonic_potential_and_vacuum_Maxwell_identities=True,
                          actual_field_vs_high_precision_error_tesla=field_formula_error,
                          no_test_region_boundary_or_forward_momentum_event=True,
                          sampled_position_refinement_error_metres=refinement_error,
                          checkpoint_lab_time_refinement_error_seconds=time_error),
        plot=dict(s_metres=samples.tolist(), lab_nanoseconds=(values[5]/c*1e9).tolist(),
                  candidate_local_x_mm=(1000*values[0]).tolist(), control_local_x_mm=(1000*old_values[0]).tolist(),
                  candidate_radius_above_nominal_mm=(1000*radial_excess).tolist(),
                  control_radius_above_nominal_mm=(1000*old_radial).tolist(),
                  candidate_field_tesla=norms.tolist(), control_field_tesla=old_norms.tolist()),
        scope='A weak-field outgoing checkpoint in a replacement local field, not completed extraction, a finite coil design, or an exactly field-free region. The saved state is initial data; the prior kicker/storage solution has not been integrated through a common global field. No current sheet crosses the beam in this test because the potential is harmonic throughout the declared box. Sources, return flux, septum, vacuum chamber and global finite-energy field remain unspecified. RF, radiation, stochastic recoil, switching, deceleration, capture and reunion are absent',
        sources=['https://arxiv.org/pdf/physics/0302069 (Papaphilippou, Wei and Talman, sectionII equations1-2: current-free scalar potential and Laplace equation)',
                 'https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.18.064001 (Muratori, Jones and Wolski: analytical Maxwell-consistent fringe models for initial beam-dynamics studies). The erf potential, local matching and trajectory here are our candidate calculation'])


def _finite_coil_loop_sum(rr, zz, coil_radius, half_separation, current,
                          sqrt, Kfun, Efun, pi, permeability):
    """Common circular-loop field and vector potential, with injectable scalar arithmetic."""
    aa = coil_radius
    ap = br = bz = 0
    if rr == 0:
        for center in (-half_separation, half_separation):
            bz += permeability*current*aa*aa/(2*(aa*aa+(zz-center)**2)**1.5)
        return ap, br, bz
    for center in (-half_separation, half_separation):
        dz = zz-center
        dd, qq = (aa+rr)**2+dz*dz, (aa-rr)**2+dz*dz
        if qq == 0:
            raise ValueError('The filament source is singular; no trajectory may cross it')
        mm = 4*aa*rr/dd
        kk, ee, root = Kfun(mm), Efun(mm), sqrt(dd)
        ap += permeability*current*root/(4*pi*rr)*((2-mm)*kk-2*ee)
        br += permeability*current*dz/(2*pi*rr*root)*(-kk+(aa*aa+rr*rr+dz*dz)/qq*ee)
        bz += permeability*current/(2*pi*root)*(kk+(aa*aa-rr*rr-dz*dz)/qq*ee)
    return ap, br, bz


def proton_finite_coil_checks(ring):
    """One finite-position current source, actual storage, and a local escape barrier."""
    import math
    import numpy as np
    from scipy.constants import c, e, m_p, mu_0
    from scipy.integrate import solve_ivp
    from scipy.special import ellipk, ellipe

    nominal = ring['nominal']
    radius, gamma0 = nominal['bending_radius_metres'], nominal['gamma']
    field0, omega0 = nominal['field_tesla'], 2*math.pi*nominal['revolution_frequency_hz']
    momentum = math.sqrt(gamma0*gamma0-1)
    coil_radius, half_separation, aperture, kick = 50., 40., .01, 1e-5
    assert math.isclose(m_p*c*momentum, e*field0*radius, rel_tol=1e-14)

    # Closed circular filament loops at z=+/-h, with equal same-direction currents.
    # scipy's elliptic argument is m=k^2, not k. No polynomial field or cutoff is added.
    def raw_loop_sum(rr, zz, current, sqrt, Kfun, Efun, pi, permeability):
        return _finite_coil_loop_sum(rr, zz, coil_radius, half_separation, current,
                                     sqrt, Kfun, Efun, pi, permeability)

    with mp.workdps(60):
        R, B0, permeability = map(lambda x: mp.mpf(str(x)), (radius, field0, mu_0))
        def high_field(rr, zz, current):
            return raw_loop_sum(rr, zz, current, mp.sqrt, mp.ellipk, mp.ellipe, mp.pi, permeability)
        current = B0/high_field(R, mp.mpf(0), mp.mpf(1))[2]
        index = -R*mp.diff(lambda rr: high_field(rr, mp.mpf(0), current)[2], R)/B0
        alpha0 = high_field(R, mp.mpf(0), current)[0]/(B0*R)
        ell = alpha0-1
        def high_alpha(rr, zz):
            return high_field(R*rr, R*zz, current)[0]/(B0*R)
        def potential(rr, zz):
            return (ell/rr-high_alpha(rr, zz))**2
        one, zero = mp.mpf(1), mp.mpf(0)
        gradient = [mp.diff(lambda rr: potential(rr, zero), one),
                    mp.diff(lambda zz: potential(one, zz), zero)]
        hessian = [mp.diff(lambda rr: potential(rr, zero), one, 2),
                   mp.diff(lambda zz: potential(one, zz), zero, 2),
                   mp.diff(lambda rr: mp.diff(lambda zz: potential(rr, zz), zero), one)]
        assert 0<index<1 and max(abs(v) for v in gradient)<mp.mpf('1e-45')
        assert max(abs(hessian[0]-2*(1-index)), abs(hessian[1]-2*index),
                   abs(hessian[2]))<mp.mpf('1e-45')

        # Bound the complete aperture boundary, not just a midplane radial slice.
        # |D_v^j(1/d)| <= j!/d^(j+1) for a unit meridional direction v follows
        # from the Legendre expansion. Integral |cos(phi)|=4 per loop.
        physical_width = mp.mpf(str(aperture))
        width = physical_width/R
        source_radius, source_height = map(lambda x: mp.mpf(str(x)), (coil_radius, half_separation))
        assert source_radius-R>physical_width and source_height>physical_width
        dmin = mp.sqrt((source_radius-R-physical_width)**2+(source_height-physical_width)**2)
        coefficient = 2*permeability*current*coil_radius/mp.pi
        derivative_bounds = [mp.factorial(j)*abs(ell)/(1-width)**(j+1)
                             +coefficient*mp.factorial(j)*R**(j-1)/(B0*dmin**(j+1))
                             for j in range(4)]
        M0, M1, M2, M3 = derivative_bounds
        third_bound = 6*M1*M2+2*M0*M3
        lower_coefficient = min(index, 1-index)-third_bound*mp.sqrt(2)*width/6
        boundary_excess = lower_coefficient*width**2
        initial_r, initial_z = one+mp.mpf('.001')/R, mp.mpf('.0005')/R
        initial_phi = ell/initial_r-high_alpha(initial_r, initial_z)
        initial_excess = initial_phi**2+2*mp.mpf(str(kick))**2-1
        assert lower_coefficient>0 and 0<initial_excess<boundary_excess
        bound = dict(
            aperture_half_width_metres=aperture,
            geometry='The meridional box |r-R|<=1 cm and |z|<=1 cm is revolved about the z axis into a toroidal tube; it is not a 1 cm Cartesian box',
            numerical_evaluation='Coefficients are evaluated with 60-digit arithmetic. These decimal evaluations are not an outward-rounded interval certificate; the continuum bound follows from the stated analytic inequalities and positive margin',
            minimum_distance_to_any_filament_metres=mp.nstr(dmin, 30),
            directional_derivative_bounds_of_M=[mp.nstr(v, 30) for v in derivative_bounds],
            third_directional_derivative_bound_of_V=mp.nstr(third_bound, 30),
            boundary_potential_excess_lower_bound=mp.nstr(boundary_excess, 30),
            initial_energy_level_excess=mp.nstr(initial_excess, 30),
            positive_margin=mp.nstr(boundary_excess-initial_excess, 30),
            boundary_to_initial_excess_ratio=mp.nstr(boundary_excess/initial_excess, 20),
            proof='For u=(rho-1,zeta), V=1+1/2*u^T*H*u+remainder, |remainder|<=C3*|u|^3/6. Every box boundary has b<=|u|<=sqrt(2)*b. Hence V-1>=[min(n,1-n)-C3*sqrt(2)*b/6]*b^2. The source-distance derivative bound holds throughout the entire box. Conserved energy and canonical angular momentum prevent this initially interior state from crossing its boundary')
        current_float, index_float, ell_float = float(current), float(index), float(ell)
        initial_phi_float = float(initial_phi)
        field_residual = mp.mpf(0)
        for rr, zz in ((R, zero), (R+mp.mpf('.008'), mp.mpf('.005')),
                       (mp.mpf(15), mp.mpf(12)), (mp.mpf(60), mp.mpf(15))):
            ap, br, bz = high_field(rr, zz, current)
            da_r = mp.diff(lambda r: high_field(r, zz, current)[0], rr)
            da_z = mp.diff(lambda z: high_field(rr, z, current)[0], zz)
            div = (mp.diff(lambda r: high_field(r, zz, current)[1], rr)+br/rr
                   +mp.diff(lambda z: high_field(rr, z, current)[2], zz))
            curl = (mp.diff(lambda z: high_field(rr, z, current)[1], zz)
                    -mp.diff(lambda r: high_field(r, zz, current)[2], rr))
            field_residual = max(field_residual, abs(da_r+ap/rr-bz)/B0,
                                 abs(da_z+br)/B0, abs(div)*R/B0, abs(curl)*R/B0)
        assert field_residual<mp.mpf('1e-45')

    def field(rr, zz):
        return raw_loop_sum(rr, zz, current_float, math.sqrt, ellipk, ellipe, math.pi, mu_0)

    # Independent source integral checks every used component, including A_phi.
    source_error = 0.
    for rr, zz in ((radius, 0.), (radius+.008, .005), (15., 12.), (60., 15.), (200., -10.)):
        numeric = np.zeros(3)
        for center in (-half_separation, half_separation):
            dz = zz-center
            def integrand(phi, component):
                cs = math.cos(phi)
                distance = math.sqrt(coil_radius**2+rr**2-2*coil_radius*rr*cs+dz*dz)
                return (coil_radius*cs/distance if component == 0 else
                        coil_radius*cs*dz/distance**3 if component == 1 else
                        (coil_radius**2-coil_radius*rr*cs)/distance**3)
            for component in range(3):
                numeric[component] += mu_0*current_float/(4*math.pi)*quad(
                    lambda phi: integrand(phi, component), 0, 2*math.pi,
                    epsabs=1e-13, epsrel=2e-12)[0]
        errors = np.abs(numeric-np.array(field(rr, zz)))/np.array([field0*radius, field0, field0])
        source_error = max(source_error, float(np.max(errors)))
    assert source_error<2e-12
    assert abs(field(radius, 0.)[2]/field0-1)<1e-14 and field(radius, 0.)[1] == 0
    axis_ap, axis_br, axis_bz = field(0., 0.)
    assert axis_ap == axis_br == 0 and axis_bz>0

    def rhs(theta, state):
        x, y, z, px, py, pz, clock = state
        rr = math.hypot(x, y)
        _, br, bz = field(radius*rr, radius*z)
        bx, by, bz = br*x/(rr*field0), br*y/(rr*field0), bz/field0
        factor = gamma0/math.sqrt(1+px*px+py*py+pz*pz)
        return (factor*px/momentum, factor*py/momentum, factor*pz/momentum,
                factor*(py*bz-pz*by), factor*(pz*bx-px*bz), factor*(px*by-py*bx),
                factor/gamma0)

    cases, plot = {}, {}
    for name, turns in (('nominal_circle', 10), ('perturbed_fixed_canonical_momentum', 20)):
        initial = np.array([1., 0., 0., 0., -momentum, 0., 0.])
        if name != 'nominal_circle':
            initial = np.array([1+.001/radius, 0., .0005/radius,
                                momentum*kick, momentum*initial_phi_float, momentum*kick, 0.])
        initial_gamma = math.sqrt(1+np.dot(initial[3:6], initial[3:6]))
        times = np.linspace(0, 2*math.pi*turns, 801)
        def radial_maximum(theta, state):
            return state[0]*state[3]+state[1]*state[4]
        def vertical_maximum(theta, state):
            return state[5]
        radial_maximum.direction = vertical_maximum.direction = -1
        solution = solve_ivp(rhs, (0, times[-1]), initial, t_eval=times,
                             events=(radial_maximum, vertical_maximum) if name != 'nominal_circle' else None,
                             method='DOP853', rtol=2e-12, atol=2e-13, max_step=.08)
        assert solution.success
        states = solution.y
        radii = np.hypot(states[0], states[1])
        energy = np.sum(states[3:6]**2, axis=0)/momentum**2
        alpha = np.array([field(radius*rr, radius*zz)[0]/(field0*radius)
                          for rr, zz in zip(radii, states[2])])
        canonical = (states[0]*states[4]-states[1]*states[3])/momentum+radii*alpha
        energy_error = float(np.max(np.abs(energy-energy[0])))
        canonical_error = float(np.max(np.abs(canonical-ell_float)))
        clock_error = float(np.max(np.abs(states[6]-times/initial_gamma))/omega0)
        radial_mm, vertical_mm = (radii-1)*radius*1000, states[2]*radius*1000
        assert energy_error<1e-11 and canonical_error<1e-11 and clock_error<1e-16
        record = dict(reference_turns=turns, elapsed_lab_seconds=times[-1]/omega0,
                      initial_gamma=initial_gamma, elapsed_proton_proper_seconds=float(states[6, -1]/omega0),
                      maximum_radial_displacement_mm=float(np.max(np.abs(radial_mm))),
                      maximum_vertical_displacement_mm=float(np.max(np.abs(vertical_mm))),
                      energy_invariant_error=energy_error, canonical_angular_momentum_error=canonical_error,
                      proper_clock_error_seconds=clock_error)
        if name == 'nominal_circle':
            position_error = float(np.max(np.abs(states[:3]-np.array([np.cos(times), -np.sin(times), np.zeros_like(times)]))))
            assert position_error<1e-9
            record['analytic_circle_position_error_in_radius_units'] = position_error
        else:
            measured = [2*math.pi/float(np.median(np.diff(ev))) for ev in solution.t_events]
            expected = [math.sqrt(1-index_float), math.sqrt(index_float)]
            assert min(map(len, solution.t_events))>=10
            assert max(abs(a-b) for a, b in zip(measured, expected))<1e-5
            assert max(np.max(np.abs(radial_mm)), np.max(np.abs(vertical_mm)))<10
            assert energy_error<float(bound['positive_margin'])*1e-3
            record['measured_tunes_from_peak_spacing'] = dict(radial=measured[0], vertical=measured[1])
            plot.update(reference_turn=(times/(2*math.pi)).tolist(),
                        radial_mm=radial_mm.tolist(), vertical_mm=vertical_mm.tolist())
        cases[name] = record

    radial_samples = np.linspace(radius-.02, 200., 301)
    radial_fields = np.array([field(rr, 0.) for rr in radial_samples])
    radial_potential = (ell_float*radius/radial_samples-radial_fields[:, 0]/(field0*radius))**2
    plot.update(midplane_radius_metres=radial_samples.tolist(),
                midplane_Bz_tesla=radial_fields[:, 2].tolist(),
                midplane_effective_potential=radial_potential.tolist())
    return dict(
        model='Two explicitly located circular filament current loops provide the same static vacuum field for storage, the aperture barrier, and the exterior profile. E=0; no old field map, switched source, or aperture cutoff is spliced in',
        parameters=dict(coil_radius_metres=coil_radius, coil_z_metres=[-half_separation, half_separation],
                        ampere_turns_per_loop=current_float, equal_same_direction_currents=True,
                        vacuum_permeability_SI=mu_0, reference_orbit_radius_metres=radius,
                        reference_gamma=gamma0, reference_Bz_tesla=field0,
                        field_index=index_float, alpha_at_reference=float(alpha0),
                        reference_canonical_angular_momentum=ell_float,
                        linear_tunes=dict(radial=math.sqrt(1-index_float), vertical=math.sqrt(index_float))),
        field_formula='For each loop, D=(a+r)^2+(z-zc)^2, Q=(a-r)^2+(z-zc)^2, m=4*a*r/D. Aphi=mu0*NI*sqrt(D)/(4*pi*r)*[(2-m)K(m)-2E(m)]; Br=mu0*NI*(z-zc)/(2*pi*r*sqrt(D))*[-K+(a^2+r^2+(z-zc)^2)*E/Q]; Bz=mu0*NI/(2*pi*sqrt(D))*[K+(a^2-r^2-(z-zc)^2)*E/Q]. Sum both loops; axis uses its regular limit',
        invariants='rho=r/R, zeta=z/R, alpha=Aphi/(B0*R), p0=e*B0*R. ell=Pphi/(p0*R)=rho*(p_phi/p0+alpha) is constant. (p_r/p0)^2+(p_z/p0)^2+V=(gamma^2-1)/(gamma0^2-1), V=(ell/rho-alpha)^2. At ell=alpha(1,0)-1 the clockwise circle has V=1, gradV=0, Hessian diag(2*(1-n),2*n)',
        aperture_bound=bound,
        initial_condition='The perturbed case starts at radial+1 mm and vertical+0.5 mm, with p_r/p0=p_z/p0=1e-5 and p_phi chosen to retain nominal ell. Its actual gamma is reported; momentum is not renormalized afterward',
        cases=cases,
        verification=dict(independent_Biot_Savart_and_vector_potential_maximum_normalized_error=source_error,
                          high_precision_vacuum_Maxwell_and_potential_maximum_normalized_residual=mp.nstr(field_residual, 8),
                          effective_potential_stationarity_and_Hessian_checked=True,
                          source_bound_covers_complete_aperture_boundary=True,
                          nominal_orbit_and_perturbed_Lorentz_runs_pass=True),
        plot=plot,
        escape_obstruction='This specified conservative state cannot leave the toroidal tube defined by the 1 cm meridional box in this unchanging coil field, despite the field decaying at infinity. Escape requires invalidating at least one fixed-field, fixed-energy, or fixed-canonical-angular-momentum hypothesis. Field controls can change the effective potential V itself and can induce electric fields, requiring a separate dynamical model. A radial potential plot alone is not a global three-dimensional escape threshold. No outgoing trajectory is established here',
        scope='Mathematical finite-position filament source and classical test-particle result, not an engineered magnet or completed extraction. Filament fields are singular at the loops and their magnetic self energy diverges; finite coil radius does not imply finite field energy. Conductor thickness, current density, mechanical forces, cooling and source preparation are unspecified. The barrier assumes this unchanging field and fixed conservative invariants; it does not prove confinement under RF, radiation, collisions, field errors or extraction controls. Old polynomial/fringe calculations remain separate benchmarks and are not joined to this field. No injection, outgoing control, capture, deceleration, reunion or measured proton clock is included',
        sources=['https://ntrs.nasa.gov/api/citations/20010038494/downloads/20010038494.pdf (Simpson, Lane, Immer and Youngquist, Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop, equations1-5: source integral, elliptic vector potential and curl; indexed primary text accessed, direct PDF fetch returned403)',
                 'https://cds.cern.ch/record/370628/files/1.pdf (Schmueser, Basic Course on Accelerator Optics, equations1.6-1.8: weak focusing; indexed primary text accessed). Coil selection, canonical barrier, and numerical trajectory checks are our calculation'])


def proton_finite_coil_escape_checks(ring, coils):
    """Conditional momentum rotation followed by escape or return in the same coil field."""
    import math
    import numpy as np
    from scipy.constants import c, m_p
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    from scipy.special import ellipk, ellipe

    params = coils['parameters']
    radius, field0 = params['reference_orbit_radius_metres'], params['reference_Bz_tesla']
    gamma0 = params['reference_gamma']
    coil_radius, height = params['coil_radius_metres'], params['coil_z_metres'][1]
    current, permeability = params['ampere_turns_per_loop'], params['vacuum_permeability_SI']
    assert params['coil_z_metres'][0] == -height
    assert radius == ring['nominal']['bending_radius_metres'] and gamma0 == ring['nominal']['gamma']
    momentum, beta = math.sqrt(gamma0*gamma0-1), math.sqrt(1-1/gamma0**2)
    omega0 = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    tail_start, checkpoint = 200., 1000.

    def field(rr, zz):
        return _finite_coil_loop_sum(rr, zz, coil_radius, height, current,
                                     math.sqrt, ellipk, ellipe, math.pi, permeability)

    with mp.workdps(60):
        R, B0, aa, hh, mu = map(lambda v: mp.mpf(str(v)),
                               (radius, field0, coil_radius, height, permeability))
        def high_field(rr, current_value):
            return _finite_coil_loop_sum(rr, mp.mpf(0), aa, hh, current_value,
                                         mp.sqrt, mp.ellipk, mp.ellipe, mp.pi, mu)
        high_current = B0/high_field(R, mp.mpf(1))[2]
        assert abs(high_current/mp.mpf(str(current))-1)<mp.mpf('1e-14')
        def F(rr):
            return rr/R*(high_field(rr, high_current)[0]/(B0*R)-1)
        alpha0 = high_field(R, high_current)[0]/(B0*R)
        candidate = mp.findroot(lambda rr: rr*high_field(rr, high_current)[2]/(B0*R)-1,
                                (mp.mpf(40), mp.mpf(43)))
        lower = F(candidate)  # Any evaluated point is a lower bound on the supremum.
        coefficient = 2*mu*high_current*aa/mp.pi
        target_width = mp.mpf('1e-10')
        # Every cell has an analytic upper bound. This is not a sampled-grid maximum.
        cells, accepted, evaluations, maximum_depth = [(R, mp.mpf(str(tail_start)), 0)], [], 0, 0
        while cells:
            left, right, depth = cells.pop()
            center, half_width = (left+right)/2, (right-left)/2
            ap, _, bz = high_field(center, high_current)
            value = center/R*(ap/(B0*R)-1)
            derivative = (center*bz/(B0*R)-1)/R
            distance = mp.sqrt(hh*hh+max(left-aa, aa-right, mp.mpf(0))**2)
            second_bound = (2*coefficient/distance**2+2*right*coefficient/distance**3)/(B0*R**2)
            upper = value+abs(derivative)*half_width+second_bound*half_width**2/2
            evaluations += 1
            if upper <= lower+target_width:
                accepted.append(upper)
                maximum_depth = max(maximum_depth, depth)
            else:
                assert depth<30
                cells.extend(((left, center, depth+1), (center, right, depth+1)))
        upper = max(accepted)
        T = mp.mpf(str(tail_start))
        assert T>aa
        tail_upper = T/R*(coefficient/(B0*R*(T-aa))-1)
        assert tail_upper<lower<=upper<lower+target_width
        angle_lower, angle_upper = mp.acos(alpha0-lower), mp.acos(alpha0-upper)
        assert 0<angle_lower<angle_upper<mp.pi/2 and alpha0<1
        requirement = dict(
            candidate_stationary_radius_metres=mp.nstr(candidate, 30),
            ell_supremum_lower_bound=mp.nstr(lower, 30), ell_supremum_upper_bound=mp.nstr(upper, 30),
            critical_angle_lower_degrees=mp.nstr(angle_lower*180/mp.pi, 30),
            critical_angle_upper_degrees=mp.nstr(angle_upper*180/mp.pi, 30),
            critical_angle_lower_radians=mp.nstr(angle_lower, 30),
            critical_angle_upper_radians=mp.nstr(angle_upper, 30),
            rotation_impulse_magnitude_over_p0_bounds=[mp.nstr(2*mp.sin(v/2), 30) for v in (angle_lower, angle_upper)],
            finite_interval_metres=[radius, tail_start], accepted_interval_count=len(accepted),
            interval_evaluations=evaluations, maximum_bisection_depth=maximum_depth,
            entire_outer_tail_F_upper_bound=mp.nstr(tail_upper, 30),
            upper_constraint_G_lower_bound=1.,
            finite_interval_proof='F(r)=(r/R)*(Aphi/(B0*R)-1), Fprime=(r*Bz/(B0*R)-1)/R. On a cell[l,u], dmin=sqrt(h^2+dist(a,[l,u])^2), CA=2*mu0*NI*a/pi, and |Fsecond|<=(2*CA/dmin^2+2*u*CA/dmin^3)/(B0*R^2). Thus the entire cell lies below F(center)+|Fprime(center)|*halfwidth+bound*halfwidth^2/2. Adaptive bisection covers every point of[R,200m]; no uniqueness assumption about stationary points is needed',
            tail_proof='For r>=T>a, |Aphi|<=CA/(r-a). Therefore F(r)<=(r/R)*(CA/(B0*R*(r-a))-1). The right side has derivative -1/R-CA*a/(B0*R^2*(r-a)^2)<0, so its value at T bounds the entire infinite tail',
            upper_constraint_proof='For positive loop current and r>0, pair source elements with opposite cos(phi). The positive-cos element is nearer, making their combined contribution to Aphi positive. Thus G=(r/R)*(Aphi/(B0*R)+1)>=r/R>=1. For clockwise rotations0<theta<=pi/2, ell=alpha0-cos(theta)<=alpha0<1, so the upper constraint is never active',
            precision_scope='Analytic derivative and tail bounds cover continuous intervals. Their coefficients are evaluated at60 decimal digits, without outward-rounded interval arithmetic; this is not a formal interval-arithmetic certificate. The two tested rotations lie far outside the narrow threshold enclosure')

    alpha_reference = field(radius, 0.)[0]/(field0*radius)
    def F_float(rr):
        return rr/radius*(field(rr, 0.)[0]/(field0*radius)-1)
    def rhs(time_scaled, state):
        x, y, z, px, py, pz, clock = state
        rr = math.hypot(x, y)
        _, br, bz = field(radius*rr, radius*z)
        bx, by, bz = br*x/(rr*field0), br*y/(rr*field0), bz/field0
        factor = gamma0/math.sqrt(1+px*px+py*py+pz*pz)
        return (factor*px/momentum, factor*py/momentum, factor*pz/momentum,
                factor*(py*bz-pz*by), factor*(pz*bx-px*bz), factor*(px*by-py*bx),
                factor/gamma0)
    def far_event(time_scaled, state):
        return math.hypot(state[0], state[1])-checkpoint/radius
    def return_event(time_scaled, state):
        return math.hypot(state[0], state[1])-1
    def turning_event(time_scaled, state):
        return state[0]*state[3]+state[1]*state[4]
    far_event.terminal = return_event.terminal = True
    far_event.direction, return_event.direction, turning_event.direction = 1, -1, -1

    cases, plots = {}, {}
    for name, degrees in (('below_threshold_5_degrees', 5.), ('above_threshold_10_degrees', 10.)):
        angle = math.radians(degrees)
        ell = alpha_reference-math.cos(angle)
        initial = np.array([1., 0., 0., momentum*math.sin(angle), -momentum*math.cos(angle), 0., 0.])
        solution = solve_ivp(rhs, (0, 100*math.pi), initial, dense_output=True,
                             events=(far_event, return_event, turning_event), method='DOP853',
                             rtol=2e-12, atol=2e-13, max_step=.05)
        assert solution.success and solution.status == 1
        is_above = degrees>float(angle_upper*180/mp.pi)
        assert bool(len(solution.t_events[0])) == is_above
        assert bool(len(solution.t_events[1])) != is_above
        samples = np.linspace(0, solution.t[-1], 501)
        values = solution.sol(samples)
        radii = np.hypot(values[0], values[1])
        field_values = np.array([field(radius*rr, radius*zz) for rr, zz in zip(radii, values[2])])
        alpha = field_values[:, 0]/(field0*radius)
        radial_momentum = (values[0]*values[3]+values[1]*values[4])/(radii*momentum)
        azimuthal_momentum = (values[0]*values[4]-values[1]*values[3])/(radii*momentum)
        canonical = radii*(azimuthal_momentum+alpha)
        energy = np.sum(values[3:6]**2, axis=0)/momentum**2
        energy_error = float(np.max(np.abs(energy-1)))
        canonical_error = float(np.max(np.abs(canonical-ell)))
        clock_error = float(np.max(np.abs(values[6]-samples/gamma0))/omega0)
        reduced_error = float(np.max(np.abs(radial_momentum**2+(ell/radii-alpha)**2-1)))
        assert max(energy_error, canonical_error, reduced_error)<1e-9 and clock_error<1e-16
        assert np.max(np.abs(values[2])) == 0 and np.max(np.abs(values[5])) == 0
        final_radius = float(radii[-1]*radius)
        lab_time, proper_time = float(samples[-1]/omega0), float(values[6, -1]/omega0)
        record = dict(
            post_control_rotation_degrees=degrees, post_control_gamma=gamma0,
            canonical_angular_momentum=ell,
            post_control_momentum_over_mc=initial[3:6].tolist(),
            required_momentum_change_over_p0=[math.sin(angle), 1-math.cos(angle), 0.],
            required_momentum_change_magnitude_kg_m_per_s=2*m_p*c*momentum*math.sin(angle/2),
            outcome='outgoing checkpoint reached' if is_above else 'outer turning point then inward return to nominal radius',
            endpoint_radius_metres=final_radius, elapsed_lab_seconds=lab_time,
            elapsed_proton_proper_seconds=proper_time,
            endpoint_position_metres=(values[:3, -1]*radius).tolist(),
            endpoint_momentum_over_mc=values[3:6, -1].tolist(),
            endpoint_field_tesla=float(np.hypot(field_values[-1, 1], field_values[-1, 2])),
            endpoint_radial_momentum_over_p0=float(radial_momentum[-1]),
            maximum_energy_invariant_error=energy_error,
            maximum_canonical_angular_momentum_error=canonical_error,
            maximum_reduced_radial_equation_error=reduced_error, proper_clock_error_seconds=clock_error)
        if is_above:
            assert ell>float(upper) and len(solution.t_events[2]) == 0
            assert np.min(radial_momentum)>0 and abs(final_radius-checkpoint)<1e-8
            independent_time = quad(lambda rr: 1/(beta*c*math.sqrt(
                1-(ell*radius/rr-field(rr, 0.)[0]/(field0*radius))**2)), radius, checkpoint,
                epsabs=1e-17, epsrel=2e-11, points=[float(candidate), 100., tail_start])[0]
            quadrature_error = abs(lab_time-independent_time)
            assert quadrature_error<1e-14
            record['independent_radial_flight_time_seconds'] = independent_time
            record['Lorentz_vs_radial_quadrature_time_error_seconds'] = quadrature_error
        else:
            assert ell<float(lower) and len(solution.t_events[2]) == 1
            predicted_turn = brentq(lambda rr: F_float(rr)-ell, radius, float(candidate), xtol=1e-11)
            turning_state = solution.y_events[2][0]
            actual_turn = radius*math.hypot(turning_state[0], turning_state[1])
            turn_error = abs(actual_turn-predicted_turn)
            assert turn_error<1e-7 and radial_momentum[-1]<0
            record.update(outer_turning_radius_metres=actual_turn,
                          independent_effective_potential_turning_radius_metres=predicted_turn,
                          turning_radius_error_metres=turn_error)
        cases[name] = record
        plots[name] = dict(lab_microseconds=(samples/omega0*1e6).tolist(),
                           x_metres=(values[0]*radius).tolist(), y_metres=(values[1]*radius).tolist(),
                           radius_metres=(radii*radius).tolist(),
                           radial_momentum_over_p0=radial_momentum.tolist(),
                           field_magnitude_tesla=np.hypot(field_values[:, 1], field_values[:, 2]).tolist())

    return dict(
        model='Conditional post-control momentum rotation at the nominal circular orbit, followed by exact Lorentz motion in the unchanged two-loop field. The initial rotation is prescribed, not generated by a modeled electromagnetic control',
        initial_state='At(r=R,z=0), p_r=p0*sin(theta), p_phi=-p0*cos(theta), p_z=0 with unchanged gamma10. ell=alpha0-cos(theta). Midplane reflection symmetry preserves z=p_z=0 exactly in this ideal model',
        escape_condition='For direct monotonically outward midplane motion, |ell/rho-alpha(rho,0)|<1 at every rho>=1. Equivalently sup[rho*(alpha-1)]<ell<inf[rho*(alpha+1)]. Equality is marginal and is not an escaping trajectory. The threshold here is for this one-parameter family, not the minimum over arbitrary controls or three-dimensional paths',
        source_parameters=params, requirements=requirement,
        outgoing_checkpoint_metres=checkpoint,
        guaranteed_minimum_midplane_distance_to_filaments_metres=height,
        cases=cases, plot=plots,
        scope='The above-threshold initial data escape the original storage barrier in the same global filament field; no field map is spliced or manually switched. A finite outward checkpoint and the analytic no-further-barrier condition are established. The checkpoint field remains nonzero. No finite control, physical kick, source preparation, septum, injection, deceleration, capture, reunion or clock readout has been demonstrated. The zero-vertical-motion invariant is exact only for the specified symmetric initial data; robustness to out-of-plane errors is untested. RF, radiation reaction and collisions are absent. Ideal filament currents retain their source singularities and divergent self energy; finite conductor engineering remains unspecified',
        sources=coils['sources'],
        derivation_attribution='The conserved-momentum escape condition, analytic finite-interval and tail bounds, threshold enclosure, and trajectories are our calculations using the previously checked loop source')


def _retarded_ring_family_control(center, permeability, light_speed, duration, peak_current,
                                  source_pairs, azimuth_order):
    """Retarded fields of weighted coaxial ring pairs; each row is radius, positive z, weight.

    Source-pair weights are normalized to one for one total current per side.
    Every azimuth integral is split at its exact retarded waveform fronts.
    """
    import math
    import numpy as np
    nodes, weights = np.polynomial.legendre.leggauss(azimuth_order)
    pairs = np.asarray(source_pairs, dtype=float)
    assert pairs.ndim == 2 and pairs.shape[1] == 3
    assert np.min(pairs[:, :2])>0 and np.min(pairs[:, 2])>=0
    assert abs(float(np.sum(pairs[:, 2]))-1)<1e-13
    radii = np.repeat(pairs[:, 0, None], 2, axis=1)
    heights = pairs[:, 1, None]*np.array([-1., 1.])[None, :]
    source_weights = pairs[:, 2]
    prefactor = permeability/(4*math.pi)

    def waveform(times, amplitude):
        phase = math.pi*np.clip(times/duration, 0., 1.)
        active = (times>0)&(times<duration)
        sine, cosine = np.sin(phase), np.cos(phase)
        return amplitude*sine**4*active, amplitude*4*math.pi/duration*sine**3*cosine*active

    def pair_sum(values):
        # Sum the reflected pair first: exact z=0 radial-field cancellation is preserved.
        return float(np.dot(source_weights, values[:, 0]+values[:, 1]))

    def control(time, position, amplitude=peak_current):
        if amplitude == 0 or time<=0:
            return np.zeros(3), np.zeros(3), np.zeros(3)
        dx, dy, zz = np.asarray(position)-center
        rr = math.hypot(dx, dy)
        dz = zz-heights
        if rr == 0:
            distances = np.hypot(radii, dz)
            current, derivative = waveform(time-distances/light_speed, amplitude)
            bz = permeability/2*pair_sum(radii**2*(
                current/distances**3+derivative/(light_speed*distances**2)))
            return np.zeros(3), np.zeros(3), np.array([0., 0., bz])
        dmin, dmax = np.hypot(rr-radii, dz), np.hypot(rr+radii, dz)
        low, high = np.maximum(dmin, light_speed*(time-duration)), np.minimum(dmax, light_speed*time)
        active = high>low
        if not np.any(active):
            return np.zeros(3), np.zeros(3), np.zeros(3)
        base = rr*rr+radii*radii+dz*dz
        denominator = 2*rr*radii
        angle_low = np.where(low == dmin, 0., np.arccos(np.clip((base-low*low)/denominator, -1., 1.)))
        angle_high = np.where(high == dmax, math.pi,
                              np.arccos(np.clip((base-high*high)/denominator, -1., 1.)))
        width = np.where(active, angle_high-angle_low, 0.)
        angles = (angle_low+angle_high)[..., None]/2+width[..., None]/2*nodes
        integration_weights = width[..., None]*weights  # Both mirror halves of each ring.
        cosine = np.cos(angles)
        distances = np.sqrt(base[..., None]-denominator[..., None]*cosine)
        current, derivative = waveform(time-distances/light_speed, amplitude)
        kernel = current/distances**3+derivative/(light_speed*distances**2)
        def azimuth_sum(values):
            return np.sum(integration_weights*values, axis=-1)
        ap = prefactor*pair_sum(azimuth_sum(radii[..., None]*cosine*current/distances))
        ep = -prefactor*pair_sum(azimuth_sum(radii[..., None]*cosine*derivative/distances))
        br = prefactor*pair_sum(azimuth_sum(radii[..., None]*cosine*dz[..., None]*kernel))
        bz = prefactor*pair_sum(azimuth_sum(
            (radii[..., None]**2-radii[..., None]*rr*cosine)*kernel))
        radial = np.array([dx/rr, dy/rr, 0.])
        azimuthal = np.array([-dy/rr, dx/rr, 0.])
        return ap*azimuthal, ep*azimuthal, br*radial+np.array([0., 0., bz])
    return control


def proton_retarded_control_checks(ring, coils, pulse_parameters=None, screening=False,
                                  surface_poloidal_order=16, launch_conditions=None,
                                  launch_resolution='screen'):
    """Continuous causal control and escape; optional parameters preserve the original defaults.

    Screening uses the same equations at lower numerical cost and finite time horizons.
    It does not certify nonescape when a checkpoint is not reached.
    """
    import math
    import numpy as np
    from scipy.constants import c, e, m_p
    from scipy.integrate import solve_ivp, quad_vec
    from scipy.special import ellipk, ellipe

    source = coils['parameters']
    radius, field0, gamma0 = (source['reference_orbit_radius_metres'], source['reference_Bz_tesla'],
                              source['reference_gamma'])
    mu, storage_current = source['vacuum_permeability_SI'], source['ampere_turns_per_loop']
    storage_a, storage_h = source['coil_radius_metres'], source['coil_z_metres'][1]
    omega0 = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    momentum = math.sqrt(gamma0**2-1)
    control_a, control_h, duration, peak_current = 5., 5., 80e-9, -1e7
    start_azimuth, checkpoint = .5, 5000.
    control_b = 0.
    if pulse_parameters is not None:
        allowed = {'source_pulse_duration_seconds', 'peak_impressed_current_per_loop_amperes',
                   'starting_azimuth_radians', 'control_loop_radius_metres',
                   'control_loop_half_separation_metres', 'control_torus_minor_radius_metres'}
        assert set(pulse_parameters)<=allowed
        duration = float(pulse_parameters.get('source_pulse_duration_seconds', duration))
        peak_current = float(pulse_parameters.get('peak_impressed_current_per_loop_amperes', peak_current))
        start_azimuth = float(pulse_parameters.get('starting_azimuth_radians', start_azimuth))
        control_a = float(pulse_parameters.get('control_loop_radius_metres', control_a))
        control_h = float(pulse_parameters.get('control_loop_half_separation_metres', control_h))
        assert duration>0 and math.isfinite(peak_current) and math.isfinite(start_azimuth)
        control_b = float(pulse_parameters.get('control_torus_minor_radius_metres', control_b))
        assert 0<control_a<math.inf and 0<control_h<math.inf
        assert 0<=control_b<min(control_a, control_h)
    finite_torus = control_b>0
    assert isinstance(surface_poloidal_order, int) and surface_poloidal_order>=4 and surface_poloidal_order%2 == 0
    source_distance_lower_bound = control_h-control_b
    localized_geometry = control_a != 5. or control_h != 5. or finite_torus
    control_resolution_scale = min(1., control_a/5., control_h/5.)
    center = np.array([radius, 0., 0.])
    prefactor = mu/(4*math.pi)

    def waveform(t, amplitude):
        values = np.asarray(t)
        phase = math.pi*np.clip(values/duration, 0., 1.)
        sine, cosine = np.sin(phase), np.cos(phase)
        active = (values>0)&(values<duration)
        return amplitude*sine**4*active, amplitude*4*math.pi/duration*sine**3*cosine*active

    def source_distances(position):
        dx, dy, zz = position-center
        local_r = math.hypot(dx, dy)
        minimum = min(math.hypot(local_r-control_a, zz-zc) for zc in (-control_h, control_h))
        maximum = max(math.hypot(local_r+control_a, zz-zc) for zc in (-control_h, control_h))
        if finite_torus:
            minimum = min(abs(math.hypot(local_r-control_a, zz-zc)-control_b)
                          for zc in (-control_h, control_h))
            maximum += control_b
        return minimum, maximum

    def make_control(order, poloidal_order=None):
        if finite_torus:
            count = surface_poloidal_order if poloidal_order is None else poloidal_order
            theta = 2*math.pi*np.arange(count)/count
            pairs = np.column_stack((control_a+control_b*np.cos(theta),
                                     control_h+control_b*np.sin(theta), np.full(count, 1./count)))
            return _retarded_ring_family_control(center, mu, c, duration, peak_current, pairs, order)
        nodes, weights = np.polynomial.legendre.leggauss(order)
        def control(t, position, amplitude=peak_current):
            if amplitude == 0 or t<=0:
                return np.zeros(3), np.zeros(3), np.zeros(3)
            dx, dy, zz = position-center
            rr = math.hypot(dx, dy)
            ap = ep = br = bz = 0.
            for zc in (-control_h, control_h):
                dz = zz-zc
                if rr == 0:
                    distance = math.hypot(control_a, dz)
                    current, derivative = waveform(t-distance/c, amplitude)
                    bz += mu*control_a**2/2*(current/distance**3+derivative/(c*distance**2))
                    continue
                dmin, dmax = math.hypot(rr-control_a, dz), math.hypot(rr+control_a, dz)
                low, high = max(dmin, c*(t-duration)), min(dmax, c*t)
                if high<=low:
                    continue
                # These are exact retarded waveform-support limits on the source, not a field cutoff.
                # Keep geometric endpoints exact: acos(1-roundoff) otherwise adds spurious arc noise.
                angle_low = (0. if low == dmin else math.acos(np.clip(
                    (rr*rr+control_a**2+dz*dz-low*low)/(2*control_a*rr), -1., 1.)))
                angle_high = (math.pi if high == dmax else math.acos(np.clip(
                    (rr*rr+control_a**2+dz*dz-high*high)/(2*control_a*rr), -1., 1.)))
                angles = (angle_low+angle_high)/2+(angle_high-angle_low)/2*nodes
                integration_weights = (angle_high-angle_low)*weights  # Both mirror halves of the loop.
                cosine = np.cos(angles)
                distances = np.sqrt(rr*rr+control_a**2+dz*dz-2*control_a*rr*cosine)
                currents, derivatives = waveform(t-distances/c, amplitude)
                kernel = currents/distances**3+derivatives/(c*distances**2)
                ap += prefactor*np.dot(integration_weights, control_a*cosine*currents/distances)
                ep -= prefactor*np.dot(integration_weights, control_a*cosine*derivatives/distances)
                br += prefactor*np.dot(integration_weights, control_a*cosine*dz*kernel)
                bz += prefactor*np.dot(integration_weights, (control_a**2-control_a*rr*cosine)*kernel)
            radial = np.array([dx/rr, dy/rr, 0.]) if rr else np.zeros(3)
            azimuthal = np.array([-dy/rr, dx/rr, 0.]) if rr else np.zeros(3)
            return ap*azimuthal, ep*azimuthal, br*radial+np.array([0., 0., bz])
        return control

    control = make_control(32)
    accurate_control = make_control(64, 2*surface_poloidal_order if finite_torus else None)
    field_probes = ((70e-9, center+np.array([3., 1., .5])),
                    (110e-9, center+np.array([-8., -15., 1.])))
    if pulse_parameters is not None:
        field_probes = ((.55*duration, center+np.array([3., 1., .5])),
                        (.9*duration, center+np.array([-8., -15., 1.])))
    if localized_geometry:
        field_probes = tuple((time, center+(position-center)*control_resolution_scale)
                             for time, position in field_probes)
    def background(position):
        rr = math.hypot(position[0], position[1])
        ap, br, bz = _finite_coil_loop_sum(rr, position[2], storage_a, storage_h, storage_current,
                                           math.sqrt, ellipk, ellipe, math.pi, mu)
        return ap, np.array([br*position[0]/rr, br*position[1]/rr, bz])
    initial = np.array([math.cos(start_azimuth), math.sin(start_azimuth), 0.,
                        momentum*math.sin(start_azimuth), -momentum*math.cos(start_azimuth), 0., 0., 0.])
    launch_mode = launch_conditions is not None
    launch_values = dict(radial_offset_metres=0., vertical_offset_metres=0.,
                        radial_direction_radians=0., vertical_direction_radians=0.)
    if launch_mode:
        assert set(launch_conditions)<=set(launch_values)
        launch_values.update({key:float(value) for key,value in launch_conditions.items()})
        assert all(math.isfinite(value) for value in launch_values.values())
        dr, dz, ar, az = (launch_values[key] for key in launch_values)
        assert radius+dr>0 and abs(ar)<math.pi/2 and abs(az)<math.pi/2
        if any(launch_values.values()):
            er = np.array([math.cos(start_azimuth),math.sin(start_azimuth),0.])
            ep = np.array([-math.sin(start_azimuth),math.cos(start_azimuth),0.])
            pp = momentum*(math.sin(ar)*math.cos(az)*er-
                math.cos(ar)*math.cos(az)*ep+np.array([0.,0.,math.sin(az)]))
            initial[:3] = ((radius+dr)*er+np.array([0.,0.,dz]))/radius
            initial[3:6] = pp

    def make_rhs(field_function, amplitude=peak_current, cleared=False):
        def rhs(s, state):
            position, pp = state[:3]*radius, state[3:6]
            gamma = math.sqrt(1+np.dot(pp, pp))
            _, magnetic = background(position)
            if cleared:
                electric = np.zeros(3)
            else:
                _, electric, extra = field_function(s/omega0, position, amplitude)
                magnetic += extra
            force = gamma0*(electric/(c*field0)+np.cross(pp, magnetic)/(gamma*field0))
            work = gamma0*np.dot(pp, electric)/(gamma*c*field0)
            return np.concatenate((gamma0/gamma*pp/momentum, force, [1/gamma, work]))
        return rhs
    def clear_event(s, state):
        return s/omega0-source_distances(state[:3]*radius)[1]/c-duration
    clear_event.direction, clear_event.terminal = 1, True
    def far_event(s, state):
        return radius*math.hypot(state[0], state[1])-checkpoint
    far_event.direction, far_event.terminal = 1, True


    if launch_mode:
        # Reuse the same retarded fields and Lorentz/work RHS, without a midplane restriction.
        profiles = {
            'screen': (24, 8, 2e-10, 2e-11, .06*control_resolution_scale, .2),
            'main': (32, 16, 2e-12, 2e-13, .015*control_resolution_scale, .15),
            'source': (64, 32, 2e-12, 2e-13, .015*control_resolution_scale, .15),
            'fine': (64, 32, 5e-13, 5e-14, .0075*control_resolution_scale, .075)}
        assert finite_torus and launch_resolution in profiles
        azimuth_order, poloidal_order, rtol, atol, pulse_step, coast_step = profiles[launch_resolution]
        fields = make_control(azimuth_order, poloidal_order)
        ceiling, checkpoint_ceiling = max(5e-6, 50*duration), 30e-6

        def signed_surface_distance(position):
            local_r = math.hypot(position[0]-center[0], position[1]-center[1])
            return min(math.hypot(local_r-control_a, position[2]-zc)-control_b
                       for zc in (-control_h, control_h))

        def storage_distance(position):
            rr = math.hypot(position[0], position[1])
            return min(math.hypot(rr-storage_a, position[2]-zc) for zc in (-storage_h, storage_h))

        def collision_event(s, state):
            return signed_surface_distance(state[:3]*radius)
        collision_event.direction, collision_event.terminal = -1, True

        def preclear_checkpoint(s, state):
            return far_event(s, state)
        preclear_checkpoint.direction, preclear_checkpoint.terminal = 1, False
        initial_gamma = math.sqrt(1+float(np.dot(initial[3:6], initial[3:6])))
        assert abs(initial_gamma-gamma0)<5e-14
        assert signed_surface_distance(initial[:3]*radius)>0
        pulse = solve_ivp(make_rhs(fields), (0, omega0*ceiling), initial,
                          events=(clear_event, collision_event, preclear_checkpoint),
                          dense_output=True, method='DOP853', rtol=rtol, atol=atol, max_step=pulse_step)
        cleared = bool(len(pulse.t_events[0]))
        hit_source = bool(len(pulse.t_events[1]))
        preclear_crossing = bool(len(pulse.t_events[2]))
        coast = None
        if not pulse.success:
            outcome = 'pulse solver failure'
        elif hit_source:
            outcome = 'control surface intersection before clearance'
        elif not cleared:
            outcome = 'causal clearance not reached within stated horizon'
        elif preclear_crossing:
            outcome = 'cylindrical checkpoint crossed before source clearance'
        elif pulse.t[-1]/omega0>=checkpoint_ceiling:
            outcome = 'source clearance occurred after checkpoint-time horizon'
        elif math.hypot(*pulse.y[:2,-1])*radius>=checkpoint:
            outcome = 'already beyond cylindrical checkpoint at source clearance'
        else:
            coast = solve_ivp(make_rhs(fields, cleared=True),
                (pulse.t[-1], omega0*checkpoint_ceiling), pulse.y[:,-1],
                events=(far_event, collision_event), dense_output=True, method='DOP853',
                rtol=rtol, atol=atol, max_step=coast_step)
            outcome = ('coast solver failure' if not coast.success else
                       'control surface intersection after clearance' if len(coast.t_events[1]) else
                       'outgoing cylindrical checkpoint reached' if len(coast.t_events[0]) else
                       'no outgoing cylindrical checkpoint within stated horizon')
        node_times = pulse.t if coast is None else np.concatenate((pulse.t, coast.t[1:]))
        node_values = pulse.y if coast is None else np.concatenate((pulse.y, coast.y[:,1:]), axis=1)
        final_time, final_state = float(node_times[-1]), node_values[:,-1]

        def trace_sample(times):
            values = np.empty((8, len(times)))
            before = times<=pulse.t[-1]
            if np.any(before):
                values[:,before] = pulse.sol(times[before])
            if np.any(~before):
                assert coast is not None
                values[:,~before] = coast.sol(times[~before])
            return values

        sample_times = np.linspace(0, final_time, 1001)
        sample_values = trace_sample(sample_times)
        all_times = np.unique(np.concatenate((node_times, sample_times)))
        # The normalized RHS has this constant coordinate-speed scale, equal to c in the model.
        speed_scale = radius*omega0*gamma0/momentum
        assert abs(speed_scale/c-1)<1e-13
        speed_bound = max(c, speed_scale)
        # Refine only loose geometric intervals, using existing dense output rather than changing physics.
        for clearance_round in range(13):
            all_values = trace_sample(all_times)
            positions = all_values[:3].T*radius
            surface_distances = np.array([signed_surface_distance(pos) for pos in positions])
            storage_distances = np.array([storage_distance(pos) for pos in positions])
            half_travel = speed_bound*np.diff(all_times)/(2*omega0)
            surface_intervals = np.minimum(surface_distances[:-1],surface_distances[1:])-half_travel
            storage_intervals = np.minimum(storage_distances[:-1],storage_distances[1:])-half_travel
            unresolved = (surface_intervals<=0)|(storage_intervals<=0)
            if not np.any(unresolved) or clearance_round == 12 or hit_source:
                break
            midpoints = (all_times[:-1][unresolved]+all_times[1:][unresolved])/2
            refined_times = np.unique(np.concatenate((all_times,midpoints)))
            if len(refined_times) == len(all_times):
                break
            all_times = refined_times
        momenta = all_values[3:6].T
        gammas = np.sqrt(1+np.sum(momenta**2,axis=1))
        work_error = float(np.max(np.abs(gammas-gamma0-all_values[7])))
        surface_lower = float(np.min(surface_intervals))
        storage_lower = float(np.min(storage_intervals))
        safe_continuation = surface_lower>0 and storage_lower>0
        clock = all_values[6]/omega0
        assert np.all(np.isfinite(all_values)) and np.min(gammas)>=1
        assert np.all(np.diff(clock)>=0)
        verification_ok = work_error<1e-9 and safe_continuation
        if outcome == 'outgoing cylindrical checkpoint reached' and not verification_ok:
            outcome = 'checkpoint reached but work or continuous-clearance check failed'

        def state_record(s, state):
            pos, pp = state[:3]*radius, state[3:6]
            rr, spherical = math.hypot(*pos[:2]), float(np.linalg.norm(pos))
            norm = float(np.linalg.norm(pp))
            ap, magnetic = background(pos)
            ell = ((state[0]*state[4]-state[1]*state[3])/momentum+
                   rr/radius*ap/(field0*radius))
            direction = pp/norm
            return dict(lab_seconds=s/omega0, proton_proper_seconds=state[6]/omega0,
                position_metres=pos.tolist(), momentum_over_mc=pp.tolist(),
                momentum_unit_direction=direction.tolist(),
                cylindrical_radius_metres=rr, spherical_radius_metres=spherical,
                gamma=math.sqrt(1+norm*norm), original_source_canonical_ell=ell,
                cylindrical_radial_direction_cosine=float(np.dot(pos[:2],direction[:2])/rr),
                spherical_radial_direction_cosine=float(np.dot(pos,direction)/spherical),
                storage_field_tesla=float(np.linalg.norm(magnetic)),
                electric_work_over_mc2=float(state[7]))

        initial_record = state_record(0., initial)
        final_record = state_record(final_time, final_state)
        clearance_record = state_record(float(pulse.t[-1]), pulse.y[:,-1]) if cleared else None
        if clearance_record is not None:
            clearance_record['minimum_continuous_source_retarded_seconds'] = (
                clearance_record['lab_seconds']-source_distances(pulse.y[:3,-1]*radius)[1]/c)
        gamma_drift = ell_drift = None
        if coast is not None:
            after_times = np.linspace(coast.t[0],coast.t[-1],301)
            after_values = coast.sol(after_times)
            after_records = [state_record(float(s), state) for s,state in zip(after_times,after_values.T)]
            gamma_drift = max(abs(row['gamma']-clearance_record['gamma']) for row in after_records)
            ell_drift = max(abs(row['original_source_canonical_ell']-
                                clearance_record['original_source_canonical_ell']) for row in after_records)
            for s,state in zip(after_times[1::50],after_values.T[1::50]):
                assert clear_event(s,state)>0
                assert max(np.linalg.norm(value) for value in fields(s/omega0,state[:3]*radius)) == 0
            assert max(gamma_drift,ell_drift)<1e-8

        far_bound = None
        if outcome == 'outgoing cylindrical checkpoint reached':
            position = final_state[:3]*radius
            pp = final_state[3:6]
            spherical = float(np.linalg.norm(position))
            pp_norm = float(np.linalg.norm(pp))
            outward_cosine = float(np.dot(position,pp)/(spherical*pp_norm))
            source_sphere = math.hypot(storage_a,storage_h)
            field_coefficient = mu*abs(storage_current)*storage_a
            momentum_SI = m_p*c*pp_norm
            assert spherical>source_sphere and outward_cosine>0
            subtraction = 2*e*field_coefficient/(momentum_SI*(spherical-source_sphere))
            margin = outward_cosine**2-subtraction
            with mp.workdps(60):
                pos_hp = list(map(lambda value:mp.mpf(str(value)),position))
                pp_hp = list(map(lambda value:mp.mpf(str(value)),pp))
                r_hp = mp.sqrt(sum(value*value for value in pos_hp))
                pnorm_hp = mp.sqrt(sum(value*value for value in pp_hp))
                w_hp = sum(x*p for x,p in zip(pos_hp,pp_hp))/(r_hp*pnorm_hp)
                mu_hp,I_hp,a_hp,h_hp,q_hp,m_hp,c_hp = map(lambda value:mp.mpf(str(value)),
                    (mu,abs(storage_current),storage_a,storage_h,e,m_p,c))
                subtraction_hp = 2*q_hp*mu_hp*I_hp*a_hp/(
                    m_hp*c_hp*pnorm_hp*(r_hp-mp.sqrt(a_hp*a_hp+h_hp*h_hp)))
                margin_hp = w_hp*w_hp-subtraction_hp
                assert abs(float(margin_hp)-margin)<2e-14
            far_bound = dict(spherical_checkpoint_radius_metres=spherical,
                storage_source_enclosing_radius_metres=source_sphere,
                storage_field_bound_coefficient_tesla_square_metres=field_coefficient,
                momentum_kg_metres_per_second=momentum_SI,
                outward_spherical_direction_cosine=outward_cosine,
                integrated_magnetic_turning_subtraction=subtraction,
                positive_radial_direction_squared_margin=margin,
                margin_60_digit_evaluation=mp.nstr(margin_hp,30),
                sufficient_outward_continuation_condition_satisfied=bool(margin>0),
                equations='For spherical radius r>s, |B|<=C/(r-s)^2, C=mu0*abs(Istorage)*astorage '
                    'for the two finite loops, s=sqrt(a_storage^2+h_storage^2). At fixed momentum p '
                    'and w=rhat dot vhat>0, d(w^2)/dr>=-2*abs(q)*|B|/p. Thus '
                    'w(r)^2>=w0^2-2*abs(q)*C/[p*(r0-s)]>0 prevents every future radial turn',
                scope='A sufficient three-dimensional bound after complete control-source clearance, '
                    'using the prescribed static storage field. It requires spherical radius and '
                    'outward spherical direction even though the checkpoint event is cylindrical. '
                    'It proves outward continuation, not aiming, capture, material survival or a '
                    'complete transport apparatus; evaluated trajectory data are not interval certified')
        plot = dict(lab_microseconds=(sample_times/omega0*1e6).tolist(),
            x_metres=(sample_values[0]*radius).tolist(),
            y_metres=(sample_values[1]*radius).tolist(),
            z_metres=(sample_values[2]*radius).tolist(),
            cylindrical_radius_metres=(np.hypot(sample_values[0],sample_values[1])*radius).tolist(),
            momentum_over_mc=sample_values[3:6].T.tolist(),
            proton_proper_microseconds=(sample_values[6]/omega0*1e6).tolist())
        return dict(
            parameters=dict(control_loop_radius_metres=control_a,
                control_loop_half_separation_metres=control_h,control_torus_minor_radius_metres=control_b,
                source_pulse_duration_seconds=duration,peak_impressed_current_per_loop_amperes=peak_current,
                starting_azimuth_radians=start_azimuth),
            launch_conditions=launch_values, initial=initial_record, pulse_clear=clearance_record,
            endpoint=final_record, outcome=outcome,
            pulse_solver_success=bool(pulse.success),
            coast_solver_success=None if coast is None else bool(coast.success),
            source_surface_intersection=bool(hit_source or (coast is not None and len(coast.t_events[1]))),
            preclear_cylindrical_checkpoint_crossing=preclear_crossing,
            numerical_profile=dict(name=launch_resolution,azimuth_order=azimuth_order,
                poloidal_order=poloidal_order,rtol=rtol,atol=atol,
                pulse_maximum_dimensionless_step=pulse_step,coast_maximum_dimensionless_step=coast_step,
                pulse_rhs_evaluations=int(pulse.nfev),coast_rhs_evaluations=0 if coast is None else int(coast.nfev)),
            horizons=dict(source_clearance_lab_seconds=ceiling,
                outgoing_checkpoint_lab_seconds=checkpoint_ceiling,
                outgoing_checkpoint_cylindrical_radius_metres=checkpoint),
            verification=dict(initial_gamma_error=abs(initial_gamma-gamma0),
                maximum_gamma_minus_initial_minus_work_residual=work_error,
                post_pulse_gamma_drift=gamma_drift,post_pulse_canonical_ell_drift=ell_drift,
                sampled_maximum_absolute_z_metres=float(np.max(np.abs(positions[:,2]))),
                continuous_control_surface_clearance_lower_bound_metres=surface_lower,
                continuous_storage_filament_clearance_lower_bound_metres=storage_lower,
                sampled_minimum_signed_control_surface_separation_metres=float(np.min(surface_distances)),
                sampled_minimum_storage_filament_distance_metres=float(np.min(storage_distances)),
                coordinate_speed_upper_bound_metres_per_second=speed_bound,
                clearance_derivation='Signed exterior torus separation and storage-filament distance '
                    'are 1-Lipschitz. On each time interval, the nearest endpoint is at most Delta t/2 '
                    'away, hence d(t)>=min(d_i,d_j)-vmax*Delta t/2. A positive minimum over intervals '
                    'excludes source crossing without assuming z=0',
                clearance_scope='The geometric/time bound applies to any relativistic path consistent '
                    'with these nodal positions. Numerical nodal errors are not interval certified; '
                    'selected source/timestep refinements are reported separately'),
            far_three_dimensional_escape_bound=far_bound,
            plot=plot,
            scope='The 5 km event is a cylindrical radius, with actual z retained. No midplane assertion '
                'or one-dimensional barrier proof is used. Arrival point, direction and time still '
                'need a receiver/capture model; reaching an outward checkpoint is not captured delivery')

    def integrate(order, step_factor, poloidal_order=None):
        fields = make_control(order, poloidal_order)
        ceiling = 5e-6 if pulse_parameters is None else max(5e-6, 50*duration)
        pulse = solve_ivp(make_rhs(fields), (0, omega0*ceiling), initial, events=clear_event,
                          dense_output=True, method='DOP853', rtol=2e-12*step_factor**2,
                          atol=2e-13*step_factor**2,
                          max_step=.015*step_factor*control_resolution_scale)
        assert pulse.success and pulse.status == 1
        assert pulse.t[-1]<omega0*30e-6 and radius*math.hypot(*pulse.y[:2,-1])<checkpoint, (
            'Full audit requires clearance before the 5 km checkpoint and 30 us horizon')
        # This reduction is exact after the entire retarded source has passed, not a field splice.
        coast = solve_ivp(make_rhs(fields, cleared=True), (pulse.t[-1], omega0*30e-6), pulse.y[:, -1],
                          events=far_event, dense_output=True, method='DOP853',
                          rtol=2e-12*step_factor**2, atol=2e-13*step_factor**2,
                          max_step=.15*step_factor)
        assert coast.success and coast.status == 1
        def sample(times):
            values = np.empty((8, len(times)))
            before = times<=pulse.t[-1]
            if np.any(before):
                values[:, before] = pulse.sol(times[before])
            if np.any(~before):
                values[:, ~before] = coast.sol(times[~before])
            return values
        return pulse, coast, sample

    if screening:
        # Same retarded fields and Lorentz/work equations; bounded, coarser screening only.
        fields = make_control(24, max(4, surface_poloidal_order//2) if finite_torus else None)
        ceiling = max(5e-6, 50*duration)
        def preclear_checkpoint_event(s, state):
            return far_event(s, state)
        preclear_checkpoint_event.direction = 1
        preclear_checkpoint_event.terminal = False
        pulse = solve_ivp(make_rhs(fields), (0, omega0*ceiling), initial,
                          events=(clear_event, preclear_checkpoint_event), dense_output=True, method='DOP853',
                          rtol=2e-10, atol=2e-11, max_step=.06*control_resolution_scale)
        state = pulse.y[:, -1]
        position = state[:3]*radius
        rr = math.hypot(position[0], position[1])
        gamma = math.sqrt(1+np.dot(state[3:6], state[3:6]))
        ap, _ = background(position)
        ell = ((state[0]*state[4]-state[1]*state[3])/momentum+rr/radius*ap/(field0*radius))
        radial = (state[0]*state[3]+state[1]*state[4])/(rr/radius*momentum)
        result = dict(parameters=dict(source_pulse_duration_seconds=duration,
                    peak_impressed_current_per_loop_amperes=peak_current,
                    starting_azimuth_radians=start_azimuth),
                    pulse_solver_success=bool(pulse.success), causal_clearance_reached=bool(pulse.status==1),
                    stop_lab_seconds=float(pulse.t[-1]/omega0), radius_metres=rr, gamma=gamma,
                    original_source_canonical_ell=ell, radial_momentum_over_original_p0=radial,
                    work_residual=float(abs(gamma-gamma0-state[7])),
                    checkpoint_crossed_before_clearance=bool(len(pulse.t_events[1])),
                    first_preclear_checkpoint_lab_seconds=(float(pulse.t_events[1][0]/omega0)
                                                           if len(pulse.t_events[1]) else None),
                    screening_source_order=24, screening_rtol=2e-10)
        if finite_torus:
            result['parameters']['control_torus_minor_radius_metres'] = control_b
            result['screening_poloidal_source_order'] = max(4, surface_poloidal_order//2)
        if localized_geometry:
            result['parameters'].update(control_loop_radius_metres=control_a,
                control_loop_half_separation_metres=control_h)
            result['screening_horizons'] = dict(causal_clearance_lab_seconds=ceiling,
                outgoing_checkpoint_lab_seconds=30e-6, outgoing_checkpoint_metres=checkpoint,
                maximum_dimensionless_pulse_step=.06*control_resolution_scale)
        if not pulse.success or pulse.status!=1:
            result['outcome']='solver failure' if not pulse.success else 'causal clearance not reached within screening time'
            return result
        if len(pulse.t_events[1]):
            result['outcome']='outgoing checkpoint crossed before causal clearance'
            return result
        if rr>=checkpoint:
            result['outcome']='already beyond outgoing checkpoint at causal clearance'
            return result
        if pulse.t[-1]>=omega0*30e-6:
            result['outcome']='causal clearance reached after checkpoint-time horizon'
            return result
        coast = solve_ivp(make_rhs(fields, cleared=True), (pulse.t[-1], omega0*30e-6), state,
                          events=far_event, method='DOP853', rtol=2e-10, atol=2e-11, max_step=.2)
        result.update(coast_solver_success=bool(coast.success),
                      outgoing_checkpoint_reached=bool(coast.status==1),
                      endpoint_radius_metres=float(radius*math.hypot(coast.y[0,-1],coast.y[1,-1])),
                      endpoint_lab_seconds=float(coast.t[-1]/omega0))
        result['outcome']=('coast solver failure' if not coast.success else
                           'outgoing checkpoint reached' if coast.status==1 else
                           'no outgoing checkpoint within screening time')
        return result

    if finite_torus:
        poloidal_nodes, poloidal_weights = np.polynomial.legendre.leggauss(2*surface_poloidal_order)
        theta = math.pi*(poloidal_nodes+1)
        direct_rings = [(control_a+control_b*math.cos(angle),
                         sign*(control_h+control_b*math.sin(angle)), weight/2)
                        for angle, weight in zip(theta, poloidal_weights) for sign in (-1, 1)]
    else:
        direct_rings = [(control_a, zc, 1.) for zc in (-control_h, control_h)]
    source_integral_error = 0.
    for time, position in field_probes:
        def direct_integrand(phi):
            answer = np.zeros(9)
            for ring_radius, zc, current_weight in direct_rings:
                tangent = np.array([-ring_radius*math.sin(phi), ring_radius*math.cos(phi), 0.])
                point = center+np.array([ring_radius*math.cos(phi), ring_radius*math.sin(phi), zc])
                displacement = position-point
                distance = float(np.linalg.norm(displacement))
                current, derivative = waveform(time-distance/c, peak_current)
                answer += prefactor*current_weight*np.concatenate((tangent*current/(distance*control_a),
                    -tangent*derivative/(c*distance),
                    np.cross(tangent, displacement)*(current/distance**3+derivative/(c*distance**2))))
            return answer
        independent, estimate = quad_vec(direct_integrand, 0, 2*math.pi, epsabs=1e-11, epsrel=1e-10)
        ap, electric, magnetic = control(time, position)
        error = float(np.max(np.abs(independent-np.concatenate((ap/control_a, electric/c, magnetic)))))
        source_integral_error = max(source_integral_error, error)
        assert error<2e-9 and estimate<2e-9

    def curl(jacobian):
        return np.array([jacobian[2, 1]-jacobian[1, 2], jacobian[0, 2]-jacobian[2, 0],
                         jacobian[1, 0]-jacobian[0, 1]])
    maxwell_residuals, potential_residuals = [], []
    for step in (.002, .001):
        maxwell_error = potential_error = 0.
        for time, position in field_probes:
            coordinates = np.concatenate(([c*time], position))
            def field_values(q):
                ap, electric, magnetic = accurate_control(q[0]/c, q[1:])
                return np.concatenate((ap, electric/c, magnetic))
            current_values = field_values(coordinates)
            derivatives = np.empty((9, 4))
            for axis in range(4):
                unit = np.eye(4)[axis]*step
                derivatives[:, axis] = (-field_values(coordinates+2*unit)+8*field_values(coordinates+unit)
                    -8*field_values(coordinates-unit)+field_values(coordinates-2*unit))/(12*step)
            potential_error = max(potential_error,
                float(np.max(np.abs(curl(derivatives[:3, 1:])-current_values[6:]))),
                float(np.max(np.abs(derivatives[:3, 0]+current_values[3:6]))))
            maxwell_error = max(maxwell_error, control_a*max(
                abs(float(np.trace(derivatives[3:6, 1:]))), abs(float(np.trace(derivatives[6:, 1:]))),
                float(np.max(np.abs(curl(derivatives[3:6, 1:])+derivatives[6:, 0]))),
                float(np.max(np.abs(curl(derivatives[6:, 1:])-derivatives[3:6, 0])))))
        assert maxwell_error<1e-7 and potential_error<1e-8
        maxwell_residuals.append(maxwell_error)
        potential_residuals.append(potential_error)

    pulse, coast, sample = integrate(32, 1.)
    high_pulse, high_coast, high_sample = integrate(64, 1.)
    fine_pulse, fine_coast, fine_sample = integrate(64, .5, 2*surface_poloidal_order if finite_torus else None)
    if finite_torus:
        poloidal_pulse, poloidal_coast, poloidal_sample = integrate(64, 1., 2*surface_poloidal_order)
        common_times = np.linspace(0, min(coast.t[-1], high_coast.t[-1], fine_coast.t[-1],
                                         poloidal_coast.t[-1]), 1001)
        step_reference_values = poloidal_sample(common_times)
        step_reference_coast = poloidal_coast
    else:
        common_times = np.linspace(0, min(coast.t[-1], high_coast.t[-1], fine_coast.t[-1]), 1001)
        step_reference_values = high_sample(common_times)
        step_reference_coast = high_coast
    main_values, high_values, fine_values = sample(common_times), high_sample(common_times), fine_sample(common_times)
    source_position_error = float(np.max(np.abs(main_values[:3]-high_values[:3]))*radius)
    time_position_error = float(np.max(np.abs(step_reference_values[:3]-fine_values[:3]))*radius)
    source_momentum_error = float(np.max(np.abs(main_values[3:6]-high_values[3:6])))
    time_momentum_error = float(np.max(np.abs(step_reference_values[3:6]-fine_values[3:6])))
    assert source_position_error<1e-6 and time_position_error<1e-4
    assert source_momentum_error<1e-8 and time_momentum_error<1e-6

    clear_state, final_state = pulse.y[:, -1], coast.y[:, -1]
    gamma_clear = math.sqrt(1+np.dot(clear_state[3:6], clear_state[3:6]))
    clear_radius = radius*math.hypot(clear_state[0], clear_state[1])
    clear_ap, _ = background(clear_state[:3]*radius)
    ell_clear = ((clear_state[0]*clear_state[4]-clear_state[1]*clear_state[3])/momentum
                 +clear_radius/radius*clear_ap/(field0*radius))
    pulse_plot_end = 1.1*pulse.t[-1] if pulse_parameters is None else min(1.1*pulse.t[-1], coast.t[-1])
    pulse_times = np.linspace(0, pulse_plot_end, 601)
    route_times = np.linspace(0, coast.t[-1], 501)
    pulse_values, route_values = sample(pulse_times), sample(route_times)
    all_values = np.concatenate((pulse_values, route_values), axis=1)
    gammas = np.sqrt(1+np.sum(all_values[3:6]**2, axis=0))
    work_error = float(np.max(np.abs(gammas-gamma0-all_values[7])))
    assert work_error<1e-10
    assert np.max(np.abs(all_values[2])) == 0 and np.max(np.abs(all_values[5])) == 0
    coast_times = np.linspace(pulse.t[-1], coast.t[-1], 501)
    coast_values = coast.sol(coast_times)
    coast_gamma = np.sqrt(1+np.sum(coast_values[3:6]**2, axis=0))
    coast_radii = np.hypot(coast_values[0], coast_values[1])
    coast_alpha = np.array([background(value[:3]*radius)[0]/(field0*radius) for value in coast_values.T])
    coast_ell = ((coast_values[0]*coast_values[4]-coast_values[1]*coast_values[3])/momentum
                 +coast_radii*coast_alpha)
    coast_radial = (coast_values[0]*coast_values[3]+coast_values[1]*coast_values[4])/(coast_radii*momentum)
    gamma_error = float(np.max(np.abs(coast_gamma-gamma_clear)))
    ell_error = float(np.max(np.abs(coast_ell-ell_clear)))
    assert gamma_error<1e-9 and ell_error<1e-9 and np.min(coast_radial)>0
    assert max(np.linalg.norm(v) for v in control(0., initial[:3]*radius)) == 0
    for s, value in zip(coast_times[1::50], coast_values.T[1::50]):
        assert clear_event(s, value)>0
        assert max(np.linalg.norm(v) for v in control(s/omega0, value[:3]*radius)) == 0
    minimum_source_distance = min(source_distances(value[:3]*radius)[0] for value in all_values.T)
    assert minimum_source_distance>=source_distance_lower_bound
    geometry_verification = {}
    if localized_geometry:
        from scipy.optimize import brentq, minimize_scalar
        distance_times = np.unique(np.concatenate((pulse.t, coast.t)))
        local_offsets = np.array([math.hypot(*(value[:2]*radius-center[:2]))-control_a
                                  for value in sample(distance_times).T])
        crossings = np.flatnonzero(local_offsets[:-1]*local_offsets[1:]<0)
        if len(crossings):
            index = int(crossings[0])
            def offset(s):
                value = sample(np.array([s]))[:, 0]
                return math.hypot(*(value[:2]*radius-center[:2]))-control_a
            near_time = brentq(offset, distance_times[index], distance_times[index+1], xtol=1e-14)
            near_distance = source_distances(sample(np.array([near_time]))[:3, 0]*radius)[0]
            assert abs(near_distance-source_distance_lower_bound)<1e-12
            distance_proof = ('Exact midplane symmetry gives distance>=h everywhere. The sampled numerical '
                'trajectory brackets local cylindrical radius=a, so continuity attains that lower bound; '
                'the dense-output root verifies it. This uses the computed trajectory, not interval-certified ODE errors')
            bracket = dict(lab_seconds=(distance_times[index:index+2]/omega0).tolist(),
                           local_radius_minus_a_metres=local_offsets[index:index+2].tolist())
        else:
            # A numerical minimum along dense output if the absolute geometric lower bound is not attained.
            candidates = [(abs(value), float(s)) for value, s in zip(local_offsets, distance_times)]
            for index in range(1, len(distance_times)-1):
                if abs(local_offsets[index])<=min(abs(local_offsets[index-1]), abs(local_offsets[index+1])):
                    def squared_offset(s):
                        value = sample(np.array([s]))[:, 0]
                        return (math.hypot(*(value[:2]*radius-center[:2]))-control_a)**2
                    minimum = minimize_scalar(squared_offset,
                        bounds=(distance_times[index-1], distance_times[index+1]), method='bounded')
                    assert minimum.success
                    candidates.append((math.sqrt(minimum.fun), float(minimum.x)))
            _, near_time = min(candidates)
            near_distance = source_distances(sample(np.array([near_time]))[:3, 0]*radius)[0]
            distance_proof = ('Numerical minimization along the computed dense-output trajectory; not a global '
                'numerical minimum certificate. Exact midplane symmetry independently gives distance>=h everywhere')
            bracket = None
        geometry_verification['continuous_trajectory_source_distance'] = dict(
            minimum_distance_metres=near_distance, occurrence_lab_seconds=near_time/omega0,
            analytic_global_lower_bound_metres=source_distance_lower_bound, lower_bound_attainment_bracket=bracket,
            qualification=distance_proof)
        geometry_verification['maximum_dimensionless_pulse_steps'] = [
            .015*control_resolution_scale, .015*control_resolution_scale, .0075*control_resolution_scale]

    zero = solve_ivp(make_rhs(control, amplitude=0.), (0, coast.t[-1]), initial, t_eval=route_times,
                     method='DOP853', rtol=2e-12, atol=2e-13, max_step=.1)
    assert zero.success
    expected_circle = np.array([np.cos(start_azimuth-route_times), np.sin(start_azimuth-route_times),
                                np.zeros_like(route_times)])
    zero_position_error = float(np.max(np.abs(zero.y[:3]-expected_circle))*radius)
    zero_clock_error = float(np.max(np.abs(zero.y[6]-route_times/gamma0))/omega0)
    assert zero_position_error<1e-7 and zero_clock_error<1e-16 and np.max(np.abs(zero.y[7])) == 0

    # Apply the actual post-pulse energy and ell, not the earlier fixed-gamma rotation threshold.
    with mp.workdps(60):
        R, B0, aa, hh, current, permeability = map(lambda v: mp.mpf(str(v)),
            (radius, field0, storage_a, storage_h, storage_current, mu))
        q = mp.sqrt((mp.mpf(str(gamma_clear))**2-1)/(mp.mpf(str(gamma0))**2-1))
        ell = mp.mpf(str(ell_clear))
        coefficient = 2*permeability*current*aa/mp.pi
        left_edge = mp.mpf(str(clear_radius))
        tail_start = max(mp.mpf(200), left_edge)
        assert R<left_edge and q>0
        pending, accepted = ([(left_edge, tail_start)] if left_edge<tail_start else []), []
        while pending:
            left, right = pending.pop()
            midpoint, half = (left+right)/2, (right-left)/2
            ap, _, bz = _finite_coil_loop_sum(midpoint, mp.mpf(0), aa, hh, current,
                                               mp.sqrt, mp.ellipk, mp.ellipe, mp.pi, permeability)
            value, derivative = midpoint/R*(ap/(B0*R)-q), (midpoint*bz/(B0*R)-q)/R
            distance = mp.sqrt(hh*hh+max(left-aa, aa-right, mp.mpf(0))**2)
            second = (2*coefficient/distance**2+2*right*coefficient/distance**3)/(B0*R**2)
            upper = value+abs(derivative)*half+second*half**2/2
            if upper<ell-mp.mpf('1e-5'):
                accepted.append(upper)
            else:
                assert right-left>mp.mpf('1e-10')
                pending.extend(((left, midpoint), (midpoint, right)))
        outer_upper = tail_start/R*(coefficient/(B0*R*(tail_start-aa))-q)
        lower_constraint_margin = ell-max(accepted+[outer_upper])
        upper_constraint_margin = q*left_edge/R-ell  # Aphi>0 for positive storage-loop current.
        assert min(lower_constraint_margin, upper_constraint_margin)>0
        escape_bound = dict(post_pulse_momentum_magnitude_over_original_p0=mp.nstr(q, 25),
            finite_interval_metres=[clear_radius, float(tail_start)], accepted_interval_count=len(accepted),
            lower_constraint_margin=mp.nstr(lower_constraint_margin, 25),
            upper_constraint_margin_lower_bound=mp.nstr(upper_constraint_margin, 25),
            entire_outer_tail_Fq_upper_bound=mp.nstr(outer_upper, 25),
            proof='With actual post-pulse q=|p|/p0 and ell, require sup[rho*(alpha-q)]<ell<inf[rho*(alpha+q)]. The same source-distance second-derivative bound encloses every finite interval from the clearance radius to200m. The decreasing analytic tail bound uses-q instead of-1. Positive Aphi gives rho*(alpha+q)>=q*rho. Thus the already outward midplane trajectory has no later radial turning barrier in the original storage field',
            precision_scope='Bounds are analytic in the continuous spatial variable; numerical coefficients and the rounded simulated clearance state are evaluated with60-digit arithmetic, without outward-rounded interval certification. The invariant/refinement errors are much smaller than the reported positive margins')

    pulse_fields = [control(s/omega0, value[:3]*radius) for s, value in zip(pulse_times, pulse_values.T)]
    distances = np.array([source_distances(value[:3]*radius) for value in pulse_values.T])
    pulse_gamma = np.sqrt(1+np.sum(pulse_values[3:6]**2, axis=0))
    route_gamma = np.sqrt(1+np.sum(route_values[3:6]**2, axis=0))
    current_samples = waveform(pulse_times/omega0, peak_current)[0]
    electric_norm = np.array([np.linalg.norm(values[1]) for values in pulse_fields])
    magnetic_norm = np.array([np.linalg.norm(values[2]) for values in pulse_fields])
    energy_gain = (gamma_clear-gamma0)*m_p*c*c
    result = dict(
        model='Nominal circular proton motion acted on continuously by the unchanged storage coils plus a finite-duration retarded electromagnetic field from two displaced closed current loops; no prescribed momentum rotation',
        parameters=dict(storage_source=source, control_loop_radius_metres=control_a,
            control_loop_centers_metres=[[radius, 0., -control_h], [radius, 0., control_h]],
            peak_impressed_current_per_loop_amperes=peak_current, source_pulse_duration_seconds=duration,
            starting_azimuth_radians=start_azimuth, outgoing_checkpoint_metres=checkpoint,
            source_quadrature_orders=[32, 64], initial_gamma=gamma0),
        source_and_fields='Each closed loop has uniform impressed I(t)=Ipeak*sin(pi*t/T)^4 for0<t<T, zero otherwise; this compact waveform is C3, not C-infinity. The distributional current is divergence-free with rho_charge=0, so Phi=0. A=mu0/(4pi)*integral I(t-d/c)*dl/d; E=-mu0/(4pi)*integral Iprime(t-d/c)*dl/d; B=mu0/(4pi)*integral [I/d^3+Iprime/(c*d^2)]*dl cross displacement. Gauss integration is split at the exact retarded support boundaries on each source loop',
        Maxwell_basis='The closed-loop current obeys charge conservation. The retarded vector potential therefore obeys the Lorenz condition and the sourced wave equation. Away from the filaments these imply all vacuum Maxwell equations; E and B above are its time derivative and curl. The displaced control axis makes the total field nonaxisymmetric about the storage center',
        pulse_clear=dict(lab_seconds=pulse.t[-1]/omega0, proton_proper_seconds=clear_state[6]/omega0,
            position_metres=(clear_state[:3]*radius).tolist(), momentum_over_mc=clear_state[3:6].tolist(),
            radius_metres=clear_radius, gamma=gamma_clear, original_source_canonical_ell=ell_clear,
            minimum_continuous_source_retarded_seconds=pulse.t[-1]/omega0-source_distances(clear_state[:3]*radius)[1]/c,
            criterion='For each horizontal loop, dmax=sqrt((local_cylindrical_radius+a)^2+(z-zc)^2). Clearance uses t-max(dmax)/c=T over both entire loops. Along any later subluminal trajectory, d(t-d/c)/dt>=1-|v|/c>0 for each fixed source point, so no later signal reentry occurs. At the proton and along its future timelike path the total field then reduces to the unchanged storage field; the emitted radiation shell can still exist elsewhere'),
        endpoint=dict(lab_seconds=coast.t[-1]/omega0, proton_proper_seconds=final_state[6]/omega0,
            position_metres=(final_state[:3]*radius).tolist(), momentum_over_mc=final_state[3:6].tolist(),
            radius_metres=radius*math.hypot(final_state[0], final_state[1]),
            gamma=float(coast_gamma[-1]), radial_momentum_over_original_p0=float(coast_radial[-1]),
            original_storage_field_tesla=float(np.linalg.norm(background(final_state[:3]*radius)[1]))),
        electric_work=dict(equation='d(gamma)/dt=e*v dot E/(m*c^2); the independent ledger state integrates this actual retarded electric field. Magnetic force does no work',
            gain_in_gamma=gamma_clear-gamma0, gain_joules=energy_gain, gain_MeV=energy_gain/(e*1e6),
            maximum_gamma_minus_initial_minus_work_residual=work_error,
            sampled_peak_electric_field_V_per_m=float(np.max(electric_norm)),
            sampled_peak_control_magnetic_field_tesla=float(np.max(magnetic_norm))),
        post_pulse_escape_bound=escape_bound,
        verification=dict(independent_Cartesian_source_integral_maximum_error_tesla=source_integral_error,
            four_point_Maxwell_residual_times_control_radius_tesla=maxwell_residuals,
            four_point_potential_field_residual_tesla=potential_residuals,
            source_quadrature_sampled_position_difference_metres=source_position_error,
            source_quadrature_sampled_momentum_over_mc_difference=source_momentum_error,
            half_step_sampled_position_difference_metres=time_position_error,
            half_step_sampled_momentum_over_mc_difference=time_momentum_error,
            source_quadrature_endpoint_time_difference_seconds=abs(coast.t[-1]-high_coast.t[-1])/omega0,
            half_step_endpoint_time_difference_seconds=abs(step_reference_coast.t[-1]-fine_coast.t[-1])/omega0,
            post_pulse_maximum_gamma_drift=gamma_error, post_pulse_maximum_canonical_ell_drift=ell_error,
            sampled_minimum_distance_to_control_filaments_metres=minimum_source_distance,
            analytic_minimum_midplane_distance_to_control_filaments_metres=source_distance_lower_bound,
            midplane_and_zero_vertical_momentum_preserved=True,
            zero_current_circle_position_error_metres=zero_position_error,
            zero_current_circle_clock_error_seconds=zero_clock_error,
            exactly_zero_control_fields_before_arrival_and_after_clearance=True, **geometry_verification),
        plot=dict(pulse=dict(lab_nanoseconds=(pulse_times/omega0*1e9).tolist(), gamma=pulse_gamma.tolist(),
            radius_metres=(np.hypot(pulse_values[0], pulse_values[1])*radius).tolist(),
            source_current_MA=(current_samples/1e6).tolist(), electric_work_over_mc2=pulse_values[7].tolist(),
            control_E_over_c_tesla=(electric_norm/c).tolist(), control_B_tesla=magnetic_norm.tolist(),
            min_retarded_nanoseconds=((pulse_times/omega0-distances[:, 1]/c)*1e9).tolist(),
            max_retarded_nanoseconds=((pulse_times/omega0-distances[:, 0]/c)*1e9).tolist()),
            route=dict(lab_microseconds=(route_times/omega0*1e6).tolist(),
            x_metres=(route_values[0]*radius).tolist(), y_metres=(route_values[1]*radius).tolist(),
            radius_metres=(np.hypot(route_values[0], route_values[1])*radius).tolist(), gamma=route_gamma.tolist(),
            proton_proper_microseconds=(route_values[6]/omega0*1e6).tolist(),
            zero_current_radius_metres=(np.hypot(zero.y[0], zero.y[1])*radius).tolist())),
        scope=('A continuous Lorentz trajectory in a specified retarded-field source model, not an engineered pulsed magnet. The10 MA current waveform is impressed; distributed drives, voltage, conductor thickness, source preparation, mechanical response and emitted-energy budget remain unspecified. Filament fields are singular and their self energy diverges. No instantaneous magnetic switching, omitted induced electric field, or field-map splice is used. The proton starts on a nominal circle; prior-day storage, physical injection, RF and radiation reaction are absent. Exact midplane symmetry supplies the tested route; out-of-plane robustness, vacuum hardware, capture, deceleration, reunion, clock readout, material samples and humans are untested. Reaching5km does not mean an exactly field-free checkpoint').replace('The10 MA current waveform', 'The10 MA current waveform' if pulse_parameters is None else f'The {abs(peak_current)/1e6:g} MA current waveform'),
        sources=['https://www.damtp.cam.ac.uk/user/tong/em/el5.pdf (David Tong, Electromagnetic Radiation, section6.1 equation6.7 and section6.1.3: retarded potentials and Lorenz gauge from current conservation; primary PDF read)',
                 coils['sources'][0]],
        derivation_attribution='The pulsed loop geometry, support-split quadrature, continuous trajectory, electric-work ledger and post-pulse escape bound are our calculations')
    if finite_torus:
        poloidal_position_error = float(np.max(np.abs(high_values[:3]-step_reference_values[:3]))*radius)
        poloidal_momentum_error = float(np.max(np.abs(high_values[3:6]-step_reference_values[3:6])))
        assert poloidal_position_error<1e-6 and poloidal_momentum_error<1e-8
        result['model'] = ('Continuous Lorentz motion in the original storage-coil field plus the exact retarded '
            'field of two displaced torus-surface impressed azimuthal currents')
        result['parameters'].update(control_torus_minor_radius_metres=control_b,
            surface_poloidal_quadrature_orders=[surface_poloidal_order, 2*surface_poloidal_order])
        result['source_and_fields'] = ('On each torus r_s=a+b*cos(theta), z_s=+/-h+b*sin(theta), '
            'dS=b*r_s*dtheta*dphi and K_phi(t)=I(t)/(2*pi*b). Thus K*dS=(I*dtheta/(2*pi))*dl, '
            'where dl=r_s*e_phi*dphi. Each weighted circular ring uses the retarded I and Iprime at its '
            'own source distance. The surface current is independent of azimuth, so its distributional '
            'divergence vanishes and rho_charge=Phi=0 is consistent. A=mu0/(4*pi) integral K(t-d/c)dS/d; '
            'E=-partial_t A and B=curl A. The compact sin^4 pulse is C3. Uniform periodic poloidal '
            'quadrature is combined with exact source-front splitting of each azimuth integral. '
            'Reflected top/bottom pairs preserve the exact midplane symmetry; reflecting the bottom '
            'poloidal parameter leaves its uniform surface-current integral unchanged')
        result['Maxwell_basis'] = ('This fixed closed surface supports an azimuthal divergence-free '
            'impressed current. Its retarded potential satisfies the Lorenz condition and sourced wave '
            'equation, hence the vacuum Maxwell equations at the proton. No instantaneous field or '
            'missing induced electric term is introduced. Material redistribution is not a field equation '
            'solved here; the prescribed surface current is the source definition')
        result['pulse_clear']['criterion'] = ('For each torus centered at zc, '
            'dmax=sqrt((local_cylindrical_radius+a)^2+(z-zc)^2)+b exactly over its entire surface. '
            'Clearance solves t-max(dmax)/c=T over both tori. Every later subluminal trajectory '
            'has d(t-d/c)/dt>=1-|v|/c>0 for every fixed source point, so all control fields remain '
            'zero at the proton after this event, while radiation may still exist elsewhere')
        verify = result['verification']
        verify['poloidal_quadrature_sampled_position_difference_metres'] = poloidal_position_error
        verify['poloidal_quadrature_sampled_momentum_over_mc_difference'] = poloidal_momentum_error
        verify['poloidal_quadrature_endpoint_time_difference_seconds'] = abs(
            high_coast.t[-1]-poloidal_coast.t[-1])/omega0
        verify['azimuth_poloidal_step_refinements'] = [
            [32, surface_poloidal_order, 1.], [64, surface_poloidal_order, 1.],
            [64, 2*surface_poloidal_order, 1.], [64, 2*surface_poloidal_order, .5]]
        verify['maximum_dimensionless_pulse_steps'] = [
            row[2]*.015*control_resolution_scale for row in verify['azimuth_poloidal_step_refinements']]
        verify['sampled_minimum_distance_to_control_surface_metres'] = verify.pop(
            'sampled_minimum_distance_to_control_filaments_metres')
        verify['analytic_minimum_midplane_distance_to_control_surface_metres'] = verify.pop(
            'analytic_minimum_midplane_distance_to_control_filaments_metres')
        verify['continuous_trajectory_source_distance']['qualification'] = (
            'The exact distance to each torus surface is abs(sqrt((r_local-a)^2+(z-zc)^2)-b). '
            'Exact midplane symmetry and h>b imply distance>=h-b globally. A bracketed numerical '
            'crossing of r_local=a attains that bound; this uses the computed continuous trajectory, '
            'without interval-certified ODE errors')
        result['scope'] = ('A specified divergence-free impressed surface current and a test-proton '
            'trajectory, not a material-response calculation or engineered conductor. The extended '
            'source is included in the Lorentz fields and its continuous wavefront-clearance event. '
            'Its current is imposed uniformly in poloidal angle; skin currents, resistivity, drive '
            'distribution, mechanical motion, heating, feeds and radiation reaction are not solved. '
            'Earlier thin-wire stored-energy and stress estimates are not an exact source-energy or '
            'mechanical solution for this surface distribution. The original storage loops remain '
            'ideal filaments. The exact midplane route is not an out-of-plane robustness or human '
            'transport result; prior storage, injection, RF, capture, deceleration and reunion are absent')
        result['derivation_attribution'] = ('The torus surface-current reduction, exact continuous source '
            'distance extrema, refined retarded-source trajectory, work ledger and escape bound are '
            'our calculations using the same primary retarded-potential framework')
    return result


def proton_control_source_feasibility_checks(control_result):
    """Finite far-zone radiation of the saved impressed currents, retaining all retardation."""
    import math
    import numpy as np
    from scipy.constants import c

    params = control_result['parameters']
    assert params.get('control_torus_minor_radius_metres', 0.) == 0, (
        'This radiation helper is for filament pairs; a torus surface requires its full poloidal source integral')
    mu = params['storage_source']['vacuum_permeability_SI']
    a, centers = params['control_loop_radius_metres'], np.array(params['control_loop_centers_metres'])
    current, duration = params['peak_impressed_current_per_loop_amperes'], params['source_pulse_duration_seconds']
    assert np.array_equal(centers[0, :2], centers[1, :2])
    midpoint = np.mean(centers, axis=0)
    height = abs(centers[1, 2]-centers[0, 2])/2
    assert a>0 and duration>0
    cached_rules = {}
    def rule(order):
        if order not in cached_rules:
            cached_rules[order] = np.polynomial.legendre.leggauss(order)
        return cached_rules[order]

    def line_sums(x, cosine_theta, pulse_duration, source_order):
        """Each loop's dimensionless integral of cos(phi)*gprime, with g=sin^4(pi*x)."""
        x = np.atleast_1d(x)
        transverse_delay = a/(c*pulse_duration)*math.sqrt(max(0., 1-cosine_theta**2))
        axial_delay = height/(c*pulse_duration)*cosine_theta
        if transverse_delay == 0:
            return np.zeros((2, len(x)))  # The line integral is exactly zero on the axis.
        nodes, weights = rule(source_order)
        contributions = []
        for shift in (-axial_delay, axial_delay):
            cosine_low = np.clip(-(x+shift)/transverse_delay, -1., 1.)
            cosine_high = np.clip((1-x-shift)/transverse_delay, -1., 1.)
            angle_low, angle_high = np.arccos(cosine_high), np.arccos(cosine_low)
            width = angle_high-angle_low
            cosine_phi = np.cos((angle_low+angle_high)[:, None]/2+width[:, None]/2*nodes)
            phase = math.pi*np.clip(x[:, None]+shift+transverse_delay*cosine_phi, 0., 1.)
            derivative = 4*math.pi*np.sin(phase)**3*np.cos(phase)
            # Doubled mirror-half integral: the Gauss half-width and mirror factor cancel.
            contributions.append(width*np.sum(weights*cosine_phi*derivative, axis=1))
        return np.array(contributions)

    def time_integrals(cosine_theta, pulse_duration, time_order, source_order):
        transverse = a/(c*pulse_duration)*math.sqrt(max(0., 1-cosine_theta**2))
        axial = height/(c*pulse_duration)*cosine_theta
        if transverse == 0:
            return 0., 0.
        # Every onset/end of a retarded source arc is an integration boundary.
        boundaries = np.unique([-sign*axial+value for sign in (-1, 1)
                                for value in (-transverse, transverse, 1-transverse, 1+transverse)])
        left, right = boundaries[:-1], boundaries[1:]
        nodes, weights = rule(time_order)
        times = ((left+right)[:, None]/2+(right-left)[:, None]/2*nodes).ravel()
        time_weights = ((right-left)[:, None]/2*weights).ravel()
        pair = line_sums(times, cosine_theta, pulse_duration, source_order)
        coherent = float(np.dot(time_weights, (pair[0]+pair[1])**2))
        self_terms = np.array([np.dot(time_weights, value**2) for value in pair])
        assert coherent>=0 and np.min(self_terms)>=0
        assert abs(self_terms[0]-self_terms[1])<1e-10*(1+max(self_terms))
        return coherent, float(np.sum(self_terms))

    def radiation(pulse_duration, angular_order, time_order, source_order):
        nodes, weights = rule(angular_order)
        totals = np.zeros(2)
        for cosine_theta, weight in zip((nodes+1)/2, weights/2):
            totals += weight*np.array(time_integrals(cosine_theta, pulse_duration, time_order, source_order))
        # Azimuthal integral 2pi and equal hemispheres are included; S=(a*I0/T)*line_sum.
        factor = mu*a*a*current*current/(4*math.pi*c*pulse_duration)
        return factor*totals

    orders = ((32, 24, 24), (64, 48, 32), (128, 96, 48))
    refinements = [radiation(duration, *order) for order in orders]
    radiated, incoherent = map(float, refinements[-1])
    refinement_difference = abs(float(refinements[-1][0]-refinements[-2][0]))
    assert radiated>0 and max(abs(float(row[0])-radiated)/radiated for row in refinements)<1e-10
    second_derivative_integral = quad(lambda x: (2*math.pi**2*(math.cos(2*math.pi*x)
                                                                -math.cos(4*math.pi*x)))**2, 0, 1)[0]
    assert abs(second_derivative_integral/(4*math.pi**4)-1)<1e-13
    def dipole_energy(pulse_duration):
        moment_amplitude = 2*math.pi*a*a*current
        return mu*moment_amplitude**2/(6*math.pi*c**3*pulse_duration**3)*4*math.pi**4

    scan = []
    for multiplier in (1., 2., 5., 10., 30., 100.):
        pulse_duration = duration*multiplier
        energy = radiated if multiplier == 1 else float(radiation(pulse_duration, 64, 48, 32)[0])
        dipole = dipole_energy(pulse_duration)
        scan.append(dict(multiplier=multiplier, duration_seconds=pulse_duration,
                         energy_joules=energy, dipole_energy_joules=dipole,
                         exact_to_dipole_ratio=energy/dipole))
    assert 0<scan[0]['exact_to_dipole_ratio']<1
    if duration == 80e-9 and a == 5. and height == 5.:
        assert scan[0]['exact_to_dipole_ratio']<.5
    assert .98<scan[3]['exact_to_dipole_ratio']<1
    assert 0<1-scan[-1]['exact_to_dipole_ratio']<2e-4
    assert all(scan[i]['exact_to_dipole_ratio']<scan[i+1]['exact_to_dipole_ratio'] for i in range(len(scan)-1))

    # The power profile uses retarded-time cuts centered at the pair midpoint.
    extent = math.hypot(a, height)/(c*duration)
    times = np.linspace(-extent, 1+extent, 1001)
    powers = np.zeros_like(times)
    nodes, weights = rule(128)
    for cosine_theta, weight in zip((nodes+1)/2, weights/2):
        values = line_sums(times, cosine_theta, duration, 48)
        powers += weight*(values[0]+values[1])**2
    powers *= mu*a*a*current*current/(4*math.pi*c*duration**2)
    integrated_power = float(np.trapezoid(powers, times*duration))
    assert np.min(powers)>=0 and abs(integrated_power/radiated-1)<1e-9
    peak_index = int(np.argmax(powers))
    angles = np.linspace(0, math.pi/2, 129)
    angular_energy = np.array([time_integrals(math.cos(angle), duration, 48, 32)[0] for angle in angles])
    angular_energy *= mu*a*a*current*current/(16*math.pi**2*c*duration)
    assert angular_energy[0] == 0 and np.min(angular_energy)>=0

    maximum_current_slew = 3*math.sqrt(3)*math.pi*abs(current)/(4*duration)
    light_circuit = 2*math.pi*a/c
    return dict(
        model='Far-zone energy emitted by the existing two prescribed control-loop currents in free space, with full retardation around each loop and coherent interference between loops. No small-source or magnetic-dipole approximation is used for the main result',
        parameters=dict(loop_radius_metres=a, half_axial_separation_metres=height,
                        loop_centers_metres=centers.tolist(), peak_current_per_loop_amperes=current,
                        source_pulse_duration_seconds=duration, vacuum_permeability_SI=mu, light_speed_SI=c),
        derivation=dict(
            field='u=t-r/c, S(n,u)=sum_loops integral Iprime(u+n dot x_source/c) dl. Then r*E_rad=-mu0*S/(4pi) and B_rad=n cross E_rad/c. There is no extra 1/c in the electric-field prefactor',
            transversality='n dot S=c*sum_loops integral d[I(u+n dot x_source/c)]=0 around closed loops; the leading radiation field is already transverse',
            energy='dU/(du*dOmega)=mu0*|S|^2/(16*pi^2*c). With n=(sqrt(1-nu^2),0,nu), only Sy remains: a*integral cos(phi)*[Iprime(u+a*sqrt(1-nu^2)*cos(phi)/c+h*nu/c)+Iprime(u+a*sqrt(1-nu^2)*cos(phi)/c-h*nu/c)]dphi. Both loop fields are added before squaring',
            translation='A common source translation adds n dot displacement/c to the retarded argument. Integration over all u removes this shift, leaving total radiated energy unchanged. An instantaneous angle-integrated power profile does depend on the chosen asymptotic retarded-time cuts',
            finiteness='The line integrals contain bounded Iprime on finite loops, and each angular waveform has compact retarded-time support. Thus this far-zone energy is finite despite the filament model having a divergent near-field self energy'),
        radiation=dict(energy_joules=radiated, energy_gigajoules=radiated/1e9,
                       incoherent_two_loop_sum_joules=incoherent,
                       interference_cross_term_joules=radiated-incoherent,
                       single_loop_energy_joules=incoherent/2,
                       dipole_energy_joules=dipole_energy(duration),
                       exact_to_dipole_ratio=radiated/dipole_energy(duration),
                       sampled_peak_far_zone_power_watts=float(powers[peak_index]),
                       sampled_peak_retarded_time_seconds=float(times[peak_index]*duration),
                       power_profile_reference_midpoint_metres=midpoint.tolist(),
                       retarded_time_support_seconds=[float(times[0]*duration), float(times[-1]*duration)],
                       qualification='This is radiated energy for these impressed currents with no added shielding or other compensating sources. Far-zone peak power is not drive-terminal peak power. The cross term is an interference contribution to a positive radiation energy, not a negative-energy material'),
        source_timing=dict(loop_circumference_light_time_seconds=light_circuit,
                           source_duration_over_circumference_light_time=duration/light_circuit,
                           maximum_abs_dI_dt_amperes_per_second=maximum_current_slew,
                           slew_extrema_source_times_seconds=[duration/3, 2*duration/3],
                           radius_over_cT=a/(c*duration), half_separation_over_cT=height/(c*duration),
                           ka_at_2pi_over_T=2*math.pi*a/(c*duration),
                           ka_at_4pi_over_T=4*math.pi*a/(c*duration),
                           interpretation=('An 80 ns globally prescribed current changes faster than a signal can traverse the 104.8 ns loop circumference. A single lumped current/inductance circuit is therefore not established. This timing fact does not prohibit previously coordinated distributed drives; those drives and finite conductor fields require a separate model' if duration == 80e-9 and a == 5. else
                           f'The {duration*1e9:g} ns pulse duration is {duration/light_circuit:.6g} circumference light-travel times. This ratio alone does not establish a uniform current or a lumped circuit. The prescribed currents still require a causal distributed drive and finite-conductor model')),
        dipole_limit=dict(equations='m(t)=2*pi*a^2*I(t); Pdip=mu0*|msecond|^2/(6*pi*c^3). For g(x)=sin(pi*x)^4 on [0,1], integral(gsecond)^2 dx=4*pi^4, giving Udip=mu0*(2*pi*a^2*I0)^2*4*pi^4/(6*pi*c^3*T^3)',
                          normalization_integral=second_derivative_integral,
                          qualification='The duration scan changes only the prescribed waveform length at fixed geometry/current. It tests the long-pulse radiation limit and does not demonstrate extraction for those alternative pulses'),
        pulse_duration_scan=scan,
        verification=dict(quadrature_orders_angular_time_source=[list(order) for order in orders],
                          refined_energy_joules=[float(row[0]) for row in refinements],
                          last_refinement_absolute_difference_joules=refinement_difference,
                          sampled_power_integral_joules=integrated_power,
                          sampled_power_integral_difference_joules=abs(integrated_power-radiated),
                          integrand_nonnegative=True, equal_single_loop_time_integrals_checked=True,
                          long_pulse_dipole_limit_checked=True,
                          precision_scope='Refinement agreement checks numerical convergence; it is not an outward-rounded interval certificate. Angular and temporal integrals cover their entire finite domains, split at the exact source-waveform fronts'),
        plot=dict(power=dict(retarded_time_nanoseconds=(times*duration*1e9).tolist(), power_watts=powers.tolist()),
                  angular=dict(polar_angle_degrees=(angles*180/math.pi).tolist(),
                               energy_per_steradian_joules=angular_energy.tolist(),
                               symmetry='Axisymmetric about the control-loop axis, with equal northern and southern hemispheres')),
        scope=f'The pulse radiates about {radiated/1e9:.3g} GJ in this free-space source model; any proposed realization retaining these fields must account for that outgoing energy. This is one energy requirement, not a complete source-energy or efficiency budget. The saved one-proton electric work is a local exchange and does not measure the global power supply. No finite wire radius, current density, distributed driver, terminal voltage, self-inductance, mechanical containment, heat, recovery or sustainable pulse circuit is established. Filament near-field self energy remains divergent. Additional shields, conductors or source currents would change the field problem and require the trajectory to be checked again',
        sources=['https://www.damtp.cam.ac.uk/user/tong/em/el5.pdf (David Tong, equation 6.7: retarded potential; section 6.1.3: current conservation and Lorenz gauge; section 6.2.4 equations 6.24-6.25: magnetic-dipole radiation. Primary PDF read)',
                 'The full retarded loop integration, interference decomposition, duration scan and source audit are our calculations'])



def proton_lower_radiation_control_checks(ring, coils, original_control, original_radiation):
    """Bounded search using actual retarded Lorentz trajectories, followed by one full audit."""
    import math
    from scipy.constants import c

    omega = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    # This records the completed adaptive screen, not a Cartesian grid or a global optimization.
    requests = [
        (160e-9, -5e6, 0.), (320e-9, -5e6, 0.),
        (160e-9, -7.5e6, 0.), (160e-9, -1e7, 0.),
        (320e-9, -7.5e6, 0.), (320e-9, -1e7, 0.),
        (640e-9, -2.5e6, 0.), (640e-9, -5e6, 0.),
        (640e-9, -7.5e6, 0.), (640e-9, -1e7, 0.),
        (640e-9, -5e6, -.6), (640e-9, -5e6, -.3),
        (640e-9, -5e6, .3), (640e-9, -5e6, .6),
        (1280e-9, -5e6, 0.), (1280e-9, -7.5e6, 0.),
        (1280e-9, -6e6, 0.), (1280e-9, -6.5e6, 0.),
        (1280e-9, -7e6, 0.), (1280e-9, -5e6, -.3),
        (1280e-9, -5e6, .3)]
    reference_energies = {}
    attempts = []
    for index, (duration, amplitude, phase_offset) in enumerate(requests, 1):
        parameters = dict(source_pulse_duration_seconds=duration,
                          peak_impressed_current_per_loop_amperes=amplitude,
                          starting_azimuth_radians=omega*(duration/2+math.sqrt(50)/c)+phase_offset)
        record = proton_retarded_control_checks(ring, coils, parameters, screening=True)
        if duration not in reference_energies:
            radiation_parameters = dict(original_control['parameters'],
                                        source_pulse_duration_seconds=duration)
            reference = proton_control_source_feasibility_checks({'parameters': radiation_parameters})
            reference_energies[duration] = reference['radiation']['energy_joules']
        energy = reference_energies[duration]*(amplitude/original_control['parameters'][
            'peak_impressed_current_per_loop_amperes'])**2
        record.update(attempt_id=index, phase_offset_radians=phase_offset,
                      radiation_energy_joules=energy)
        attempts.append(record)
    eligible = [record for record in attempts if record['outcome']=='outgoing checkpoint reached']
    assert eligible, 'No screened candidate reached the checkpoint within the stated acceptance horizons'
    chosen = min(eligible, key=lambda record: record['radiation_energy_joules'])
    selected = proton_retarded_control_checks(ring, coils, chosen['parameters'])
    source_audit = proton_control_source_feasibility_checks(selected)
    selected_energy = source_audit['radiation']['energy_joules']
    assert abs(selected_energy/chosen['radiation_energy_joules']-1)<1e-12
    assert selected['pulse_clear']['radius_metres']<5000
    assert selected['pulse_clear']['lab_seconds']<selected['endpoint']['lab_seconds']<30e-6
    assert selected['endpoint']['radial_momentum_over_original_p0']>0
    assert float(selected['post_pulse_escape_bound']['lower_constraint_margin'])>0
    assert float(selected['post_pulse_escape_bound']['upper_constraint_margin_lower_bound'])>0
    original_energy = original_radiation['radiation']['energy_joules']
    assert selected_energy<original_energy/1000
    assert max(record['work_residual'] for record in attempts)<1e-9
    outcomes = {}
    for record in attempts:
        outcomes[record['outcome']] = outcomes.get(record['outcome'], 0)+1
        record['full_audit_performed'] = record['attempt_id']==chosen['attempt_id']
    return dict(
        model='A bounded pulse-duration, current-amplitude and initial-orbit-phase search in the same retarded control-loop and storage-coil field, followed by full checks of the selected continuous trajectory',
        search=dict(
            attempt_count=len(attempts), selected_attempt_id=chosen['attempt_id'],
            selection='Lowest radiated energy among the screened candidates that reached 5 km within the stated horizons; only the selected candidate receives the full Maxwell, refinement, work and continuous escape-bound audit',
            phase_rule='Initial azimuth=omega0*(T/2+sqrt(a^2+h^2)/c)+listed phase offset, with a=h=5 m. The source current starts at t=0 for every trial; initial phase changes the proton timing without switching or splicing a field',
            fixed_source='The same two control loops of radius 5 m at (R,0,+/-5 m), the same original storage coils, gamma=10 initial circular momentum, and the same compact sin^4 waveform',
            screening='Source quadrature order 24, DOP853 rtol=2e-10 and atol=2e-11. Clearance is sought until max(5 us,50*T), with dimensionless maximum step .06. A subsequent coast is followed until 30 us, with maximum step .2',
            acceptance='The proton must causally clear the entire control source before reaching 5 km and before the 30 us checkpoint horizon, then cross 5 km outward. Selected full verification also proves no later midplane radial barrier using its actual post-pulse energy and canonical momentum',
            outcome_qualification='A time-limited failure to reach the checkpoint is not a proof of permanent trapping or impossibility. Unclear pulses, failed ODE solves, and checkpoint-time failures are distinct outcomes. No failed invariant-bound candidate was discarded: the one selected full audit passed',
            radiation_computation='At each tested duration, reuse the exact all-angle retarded source-energy integral at the original current amplitude, then scale by Ipeak^2. This scaling follows linear Maxwell fields; total radiation is independent of the initial proton phase. The selected energy is independently recomputed at its own amplitude',
            outcomes=outcomes),
        attempted_cases=attempts,
        comparison=dict(original_radiated_energy_joules=original_energy,
                        selected_radiated_energy_joules=selected_energy,
                        radiation_reduction_factor=original_energy/selected_energy,
                        radiation_reduction_fraction=1-selected_energy/original_energy,
                        current_amplitude_ratio=abs(selected['parameters'][
                            'peak_impressed_current_per_loop_amperes']/original_control['parameters'][
                            'peak_impressed_current_per_loop_amperes']),
                        pulse_duration_ratio=selected['parameters']['source_pulse_duration_seconds']/
                            original_control['parameters']['source_pulse_duration_seconds']),
        selected_control=selected,
        selected_source_feasibility=source_audit,
        scope='This verifies a lower-radiation outgoing trajectory for the specified impressed-current model and exact midplane initial state. It is not a global optimum, hardware design or proof that a 6 MA pulse can be generated. Finite conductors, causal drive distribution, stresses, heating and full source-energy accounting remain unmodeled; the radiation figure is one calculated outgoing-energy term. Different feeds or shields would require new total-field and trajectory calculations. The proton is not a material grain or human, and there is no injection, long-duration storage, RF/radiation-reaction combination, capture, deceleration, reunion or clock-readout demonstration',
        sources=source_audit['sources'])

def proton_finite_conductor_checks(lower_radiation_result, wire_radii_metres=None):
    """Quasistatic thin-conductor estimates for the selected control pair alone."""
    selected = lower_radiation_result['selected_control']
    radiation = lower_radiation_result['selected_source_feasibility']
    params = selected['parameters']
    centers = params['control_loop_centers_metres']
    assert centers[0][:2] == centers[1][:2]
    with mp.workdps(60):
        a, current, duration, mu, c = map(lambda value: mp.mpf(str(value)), (
            params['control_loop_radius_metres'],
            params['peak_impressed_current_per_loop_amperes'],
            params['source_pulse_duration_seconds'],
            params['storage_source']['vacuum_permeability_SI'],
            radiation['parameters']['light_speed_SI']))
        separation = abs(mp.mpf(str(centers[1][2]))-mp.mpf(str(centers[0][2])))
        assert min(a, duration, mu, separation)>0 and current != 0

        def mutual(a1, a2, distance):
            parameter = 4*a1*a2/((a1+a2)**2+distance**2)
            return mu*mp.sqrt(a1*a2/parameter)*(
                (2-parameter)*mp.ellipk(parameter)-2*mp.ellipe(parameter))

        def self_external(radius, wire_radius):
            return mu*radius*(mp.log(8*radius/wire_radius)-2)

        mutual_value = mutual(a, a, separation)
        neumann_value = mu*a*a/2*mp.quad(
            lambda phi: mp.cos(phi)/mp.sqrt(2*a*a+separation**2-2*a*a*mp.cos(phi)),
            [0, mp.pi, 2*mp.pi])
        mutual_relative_error = abs(neumann_value/mutual_value-1)
        assert mutual_relative_error<mp.mpf('1e-50')
        assert mutual(a*.9, a*1.1, separation) == mutual(a*1.1, a*.9, separation)
        mutual_radial_derivative = mp.diff(lambda a1: mutual(a1, a, separation), a)
        mutual_axial_derivative = mp.diff(lambda distance: mutual(a, a, distance), separation)
        # Independent Biot-Savart field per ampere at the other loop's centerline.
        def distance_cube(phi):
            return (2*a*a+separation**2-2*a*a*mp.cos(phi))**mp.mpf('1.5')
        partner_bz_per_amp = mu/(4*mp.pi)*mp.quad(
            lambda phi: (a*a-a*a*mp.cos(phi))/distance_cube(phi), [0, mp.pi, 2*mp.pi])
        partner_br_per_amp = mu/(4*mp.pi)*mp.quad(
            lambda phi: a*separation*mp.cos(phi)/distance_cube(phi), [0, mp.pi, 2*mp.pi])
        radial_error = abs(mutual_radial_derivative/(2*mp.pi*a*partner_bz_per_amp)-1)
        axial_error = abs(mutual_axial_derivative/(-2*mp.pi*a*partner_br_per_amp)-1)
        assert max(radial_error, axial_error)<mp.mpf('1e-50')
        assert mutual_value>0 and mutual_radial_derivative>0 and mutual_axial_derivative<0
        slew = 3*mp.sqrt(3)*mp.pi*abs(current)/(4*duration)
        saved_slew = mp.mpf(str(radiation['source_timing']['maximum_abs_dI_dt_amperes_per_second']))
        assert abs(slew/saved_slew-1)<mp.mpf('1e-14')
        g = lambda x: mp.sin(mp.pi*x)**4
        gp = lambda x: 4*mp.pi*mp.sin(mp.pi*x)**3*mp.cos(mp.pi*x)
        assert abs(abs(gp(mp.mpf(1)/3))-3*mp.sqrt(3)*mp.pi/4)<mp.mpf('1e-55')
        rise = mp.quad(lambda x: 2*g(x)*gp(x), [0, mp.mpf('.5')])
        fall = mp.quad(lambda x: 2*g(x)*gp(x), [mp.mpf('.5'), 1])
        assert max(abs(rise-1), abs(fall+1))<mp.mpf('1e-50')
        mutual_tension = current**2*mutual_radial_derivative/(2*mp.pi)
        rows, hoop_errors = [], []
        wire_radii = ('.001', '.01', '.05', '.1') if wire_radii_metres is None else map(str, wire_radii_metres)
        for wire_radius in map(mp.mpf, wire_radii):
            aspect = wire_radius/a
            assert 0<aspect<=mp.mpf('.02' if wire_radii_metres is None else '.1')
            assert 2*wire_radius<separation
            inductance = self_external(a, wire_radius)
            equal_current_inductance = inductance+mutual_value
            assert inductance>mutual_value>0
            energy = equal_current_inductance*current**2
            emf = equal_current_inductance*slew
            surface_field = mu*abs(current)/(2*mp.pi*wire_radius)
            pressure = surface_field**2/(2*mu)
            self_tension = mu*current**2/(4*mp.pi)*(mp.log(8*a/wire_radius)-1)
            total_tension = self_tension+mutual_tension
            area = mp.pi*wire_radius**2
            # Vary both major radii at fixed wire radii and separation: twice the single-loop force.
            energy_derivative = mp.diff(
                lambda radius: (self_external(radius, wire_radius)+
                                mutual(radius, radius, separation))*current**2, a)
            error = abs(energy_derivative/(4*mp.pi*total_tension)-1)
            assert error<mp.mpf('1e-50')
            hoop_errors.append(error)
            next_order = mu*a*((1+mp.mpf('.75')*aspect**2)*mp.log(8*a/wire_radius)
                               -2-mp.mpf('1.5')*aspect**2)
            correction = next_order/inductance-1
            assert abs(correction-mp.mpf('.75')*aspect**2)<mp.mpf('1e-55')
            rows.append(dict(
                wire_radius_metres=float(wire_radius), wire_radius_over_loop_radius=float(aspect),
                external_self_inductance_henries=float(inductance),
                equal_current_quasistatic_inductance_per_loop_henries=float(equal_current_inductance),
                control_pair_peak_quasistatic_magnetic_energy_joules=float(energy),
                quasistatic_peak_inductive_emf_per_loop_volts=float(emf),
                leading_local_self_surface_field_tesla=float(surface_field),
                leading_local_self_magnetic_pressure_pascals=float(pressure),
                self_hoop_tension_newtons=float(self_tension),
                mutual_hoop_tension_newtons=float(mutual_tension),
                pair_only_hoop_tension_newtons=float(total_tension),
                assumed_solid_load_area_square_metres=float(area),
                pair_only_equivalent_hoop_stress_pascals=float(total_tension/area),
                next_order_surface_self_inductance_henries=float(next_order),
                next_order_self_inductance_relative_correction=float(correction),
                mutual_finite_wire_expansion_parameter=float(wire_radius**2/(a*separation))))
        for key in ('external_self_inductance_henries',
                    'control_pair_peak_quasistatic_magnetic_energy_joules',
                    'quasistatic_peak_inductive_emf_per_loop_volts',
                    'leading_local_self_magnetic_pressure_pascals',
                    'pair_only_equivalent_hoop_stress_pascals'):
            assert all(rows[j][key]>rows[j+1][key]>0 for j in range(len(rows)-1))
        return dict(
            model='Thin-wire surface-current and quasistatic estimates for the selected two control loops, with exact filament mutual inductance. These are control-pair contributions, not total apparatus loads or exact transient source energy',
            parameters=dict(loop_major_radius_metres=float(a), loop_centers_metres=centers,
                axial_center_separation_metres=float(separation),
                signed_peak_current_per_loop_amperes=[float(current), float(current)],
                source_pulse_duration_seconds=float(duration),
                maximum_abs_dI_dt_amperes_per_second=float(slew),
                voltage_extrema_source_times_seconds=[float(duration/3), float(2*duration/3)],
                peak_current_and_energy_source_time_seconds=float(duration/2),
                vacuum_permeability_SI=float(mu),
                loop_circumference_light_time_seconds=float(2*mp.pi*a/c),
                pulse_duration_over_circumference_light_time=float(duration*c/(2*mp.pi*a)),
                ka_at_2pi_over_T=float(2*mp.pi*a/(c*duration)),
                ka_at_4pi_over_T=float(4*mp.pi*a/(c*duration))),
            equations=dict(
                self_inductance='Lext=mu0*a*[ln(8*a/b)-2], with loop major radius a and conductor minor radius b; internal inductance is omitted in the assumed surface-current limit',
                mutual_inductance='m=4*a1*a2/[(a1+a2)^2+d^2]; M=mu0*sqrt(a1*a2/m)*[(2-m)*K(m)-2*E(m)], where K and E take the elliptic parameter m, not the modulus',
                energy_and_emf='Wpair=.5*L*I1^2+.5*L*I2^2+M*I1*I2=(L+M)*I^2 for equal signed currents. The quasistatic per-loop flux derivative is V=(L+M)*Iprime; this is a loop EMF estimate, not a modeled feed-terminal voltage',
                surface_pressure='Bself_surface approximately mu0*abs(I)/(2*pi*b), Pself approximately Bself_surface^2/(2*mu0): leading local compressive pressure of a straight surface-current wire segment',
                hoop_load='At fixed currents, b, other loop radius and d, the coenergy derivative gives Qradial=.5*I^2*dL/da+I^2*dM/da1. Hoop tension=Qradial/(2*pi); self term=mu0*I^2*[ln(8*a/b)-1]/(4*pi). Divide self+mutual tension by pi*b^2 only for a full solid cross-section carrying this unopposed pair load',
                finite_radius_indicator='Dengler surface-current next term: Lnext=mu0*a*{[1+3*b^2/(4*a^2)]ln(8*a/b)-2-3*b^2/(2*a^2)}, so (Lnext-Lext)/Lext=3*b^2/(4*a^2). This is an asymptotic correction indicator, not a certified remainder bound'),
            mutual=dict(filament_mutual_inductance_henries=float(mutual_value),
                elliptic_parameter_m=float(4*a*a/(4*a*a+separation**2)),
                derivative_one_major_radius_henries_per_metre=float(mutual_radial_derivative),
                derivative_separation_henries_per_metre=float(mutual_axial_derivative),
                pair_only_axial_attraction_per_loop_newtons=float(-current**2*mutual_axial_derivative),
                geometry_scope='Exact for circular filament centerlines. Finite conductor redistribution is uncomputed; b^2/(a*d) is a small geometric expansion parameter for mutual corrections, not a numerical error bound'),
            rows=rows,
            verification=dict(
                neumann_vs_elliptic_mutual_relative_error=mp.nstr(mutual_relative_error, 8),
                mutual_radial_derivative_vs_Biot_Savart_force_relative_error=mp.nstr(radial_error, 8),
                mutual_axial_derivative_vs_Biot_Savart_force_relative_error=mp.nstr(axial_error, 8),
                maximum_pair_coenergy_derivative_vs_hoop_tension_relative_error=mp.nstr(max(hoop_errors), 8),
                normalized_quasistatic_rising_field_work=float(rise),
                normalized_quasistatic_falling_field_work=float(fall),
                positive_inductance_matrix_all_rows=True,
                precision_scope='60-digit arithmetic checks the stated formulas and independent integrals. It does not certify the finite-conductor, current-distribution, retardation or mechanical approximations'),
            resistive_heating=dict(joules=None, reason='No conductivity, temperature dependence, skin depth, current distribution or material is specified; resistive heating is not estimated'),
            scope='The surface-current limit is an assumption, not a derived skin-depth calculation. Wire size supplies an outer radius and an assumed full solid load area. Local pressure and equivalent hoop stress describe different force directions; neither is a yield or survival calculation. Supports, hollow sections, material response, stress waves and original storage-field forces can change actual stresses. Positive pair magnetic energy excludes original storage-field cross energy and is not a universal drive-energy floor. Ideal quasistatic rising/falling work is stored and returned within that approximation; it is not dissipated energy and cannot be added to exact radiation as an established apparatus budget. The pulse spans about12.2 circumference light-travel times and contains higher frequencies; exact transient stored energy and actual terminal voltage are not established. Feeds, capacitance, displacement currents, redistribution, breakdown, source impedance and recovery remain unmodeled. Real conductors or feeds require a new field and trajectory check'.replace('about12.2 circumference', 'about12.2 circumference' if wire_radii_metres is None else f'about {float(duration*c/(2*mp.pi*a)):.3g} circumference'),
            sources=[
                'https://www.aemjournal.org/index.php/AEM/article/view/331/pdf ; Dengler, Advanced Electromagnetics 5(1),1-8(2016), DOI10.7716/aem.v5i1.331, Eqs2-3, sections3.2/4.1, pages1/2/5: current energy, SI Neumann integral, thin circular-loop inductance and surface-current correction. Published primary PDF read',
                'https://dlmf.nist.gov/19.34 ; NIST DLMF19.34.1-19.34.6: exact coaxial-circle mutual integral and elliptic forms, in cgs; the SI Neumann prefactor is verified independently here',
                'https://w3.pppl.gov/~hji/paper/Myers16-pop.pdf ; Myers et al., Physics of Plasmas23,112102(2016), DOI10.1063/1.4966691, page4 Eqs6-8: large-aspect-ratio hoop force and inductance, with zero internal inductance for surface current. The magnetic force formula is used without importing plasma material properties',
                'https://web.mit.edu/6.013_book/www/chapter11/11.7.html ; Haus and Melcher, section11.7: magnetoquasistatic fixed-current coenergy and coupled electrical terminal pairs',
                'https://ocw.mit.edu/courses/6-013-electromagnetics-and-applications-spring-2009/e163397146247bc4558e1acfa8ff74d2_MIT6_013S09_chap06.pdf ; MIT6.013 Chapter6, Eqs6.3.30/6.4.15: magnetic surface pressure. Primary course PDF read'])


def proton_localized_control_checks(ring, coils, baseline_control, baseline_conductor,
                                    attempted_cases=None):
    """Eighteen bounded geometry/pulse screens and one full continuous-trajectory audit."""
    import copy
    import math
    from scipy.constants import c
    from scipy.special import ellipk, ellipe

    omega = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    reference_parameters = baseline_control['selected_control']['parameters']
    baseline_row = next(row for row in baseline_conductor['rows']
                        if row['wire_radius_over_loop_radius'] == .02)
    reference_energy = baseline_row['control_pair_peak_quasistatic_magnetic_energy_joules']
    # The first four retain the established starting phase. Later cases use the new source delay.
    requests = []
    for a, current in ((2.5, -6e6), (2.5, -4e6), (1.25, -6e6), (1.25, -4e6)):
        requests.append(dict(control_loop_radius_metres=a, control_loop_half_separation_metres=a,
            peak_impressed_current_per_loop_amperes=current,
            source_pulse_duration_seconds=reference_parameters['source_pulse_duration_seconds'],
            starting_azimuth_radians=reference_parameters['starting_azimuth_radians']))
    geometric_requests = [
        (2.5, -5e6, 1280e-9, 0.), (1.25, -5e6, 1280e-9, 0.),
        (1.25, -5.5e6, 1280e-9, 0.), (1.25, -4.5e6, 1280e-9, 0.),
        (2.5, -6e6, 640e-9, 0.), (1.25, -6e6, 640e-9, 0.),
        (1.25, -5e6, 640e-9, 0.), (1.25, -4e6, 640e-9, 0.),
        (1.25, -4e6, 2560e-9, 0.), (1.25, -5e6, 2560e-9, 0.),
        (1.25, -4e6, 1920e-9, 0.), (1.25, -5e6, 1920e-9, 0.),
        (1.25, -5.25e6, 1280e-9, 0.), (1.25, -5e6, 1280e-9, .3)]
    for a, current, duration, phase_offset in geometric_requests:
        requests.append(dict(control_loop_radius_metres=a, control_loop_half_separation_metres=a,
            peak_impressed_current_per_loop_amperes=current, source_pulse_duration_seconds=duration,
            starting_azimuth_radians=omega*(duration/2+math.sqrt(2)*a/c)+phase_offset))
    assert len(requests) == 18
    attempts = ([proton_retarded_control_checks(ring, coils, parameters, screening=True)
                 for parameters in requests] if attempted_cases is None else copy.deepcopy(attempted_cases))
    assert len(attempts) == len(requests)
    permeability = coils['parameters']['vacuum_permeability_SI']
    for index, (record, parameters) in enumerate(zip(attempts, requests), 1):
        assert record['parameters'] == parameters
        a, h, current = (parameters['control_loop_radius_metres'],
                          parameters['control_loop_half_separation_metres'],
                          parameters['peak_impressed_current_per_loop_amperes'])
        b, m = .02*a, a*a/(a*a+h*h)
        mutual = permeability*a/math.sqrt(m)*((2-m)*ellipk(m)-2*ellipe(m))
        energy = (permeability*a*(math.log(8*a/b)-2)+mutual)*current**2
        # Homologous sources at fixed b/a provide an independent scaling check.
        scaled = reference_energy*a/reference_parameters['control_loop_radius_metres']*(
            current/reference_parameters['peak_impressed_current_per_loop_amperes'])**2
        assert abs(energy/scaled-1)<2e-14
        record.update(attempt_id=index, ranking_wire_radius_metres=b,
            estimated_control_pair_peak_quasistatic_magnetic_energy_joules=energy)
    eligible = [record for record in attempts if record['outcome'] == 'outgoing checkpoint reached']
    if not eligible:
        return dict(model='Bounded localized retarded-control search', search=dict(attempt_count=len(attempts),
            outcome='No candidate reached the outgoing checkpoint within the stated acceptance horizons'),
            attempted_cases=attempts, selected_control=None,
            scope='A bounded null result does not prove permanent trapping or impossibility')
    chosen = min(eligible, key=lambda record:
                 record['estimated_control_pair_peak_quasistatic_magnetic_energy_joules'])
    selected = proton_retarded_control_checks(ring, coils, chosen['parameters'])
    radiation = proton_control_source_feasibility_checks(selected)
    pair = dict(selected_control=selected, selected_source_feasibility=radiation)
    a = selected['parameters']['control_loop_radius_metres']
    conductors = proton_finite_conductor_checks(pair, [aspect*a for aspect in (.0002, .002, .01, .02)])
    fixed_conductor = proton_finite_conductor_checks(pair, [.1])
    selected_row = next(row for row in conductors['rows'] if row['wire_radius_over_loop_radius'] == .02)
    fixed_row = fixed_conductor['rows'][0]
    energy = selected_row['control_pair_peak_quasistatic_magnetic_energy_joules']
    assert abs(energy/chosen['estimated_control_pair_peak_quasistatic_magnetic_energy_joules']-1)<2e-14
    assert energy<reference_energy
    assert selected['pulse_clear']['lab_seconds']<selected['endpoint']['lab_seconds']<30e-6
    assert selected['pulse_clear']['radius_metres']<5000
    assert selected['endpoint']['radial_momentum_over_original_p0']>0
    assert min(float(selected['post_pulse_escape_bound'][key]) for key in
               ('lower_constraint_margin', 'upper_constraint_margin_lower_bound'))>0
    source_distance = selected['verification']['continuous_trajectory_source_distance']
    clearance_rows = []
    for label, row in (('fixed aspect b/a=.02', selected_row), ('fixed wire radius b=.1 m', fixed_row)):
        b = row['wire_radius_metres']
        lower = source_distance['analytic_global_lower_bound_metres']-b
        actual = source_distance['minimum_distance_metres']-b
        separation = conductors['parameters']['axial_center_separation_metres']
        assert min(lower, actual, separation-2*b)>0
        clearance_rows.append(dict(case=label, conductor_minor_radius_metres=b,
            minimum_trajectory_distance_to_outer_conductor_metres=actual,
            analytic_midplane_lower_bound_to_outer_conductor_metres=lower,
            partner_conductor_surface_gap_metres=separation-2*b,
            minimum_occurrence_lab_seconds=source_distance['occurrence_lab_seconds']))
    outcomes = {}
    for record in attempts:
        outcomes[record['outcome']] = outcomes.get(record['outcome'], 0)+1
        record['full_audit_performed'] = record['attempt_id'] == chosen['attempt_id']
    compared_keys = (
        'control_pair_peak_quasistatic_magnetic_energy_joules',
        'quasistatic_peak_inductive_emf_per_loop_volts',
        'leading_local_self_surface_field_tesla',
        'leading_local_self_magnetic_pressure_pascals',
        'pair_only_equivalent_hoop_stress_pascals')
    return dict(
        model='Smaller displaced control loops, the original storage coils, and continuous exact-retarded-source Lorentz motion; stored-energy comparison at equal conductor aspect ratio',
        search=dict(attempt_count=len(attempts), new_case_limit=18, selected_attempt_id=chosen['attempt_id'],
            outcomes=outcomes,
            selection='Lowest estimated control-pair quasistatic stored magnetic energy at b/a=.02 among these screened successes. Only the chosen candidate receives the full trajectory, source and conductor audits; this is not a global optimum',
            phase_rule='Cases1-4 retain the established 5 m-pair starting phase. Other cases use phi0=omega0*[T/2+sqrt(a^2+h^2)/c]+offset, with a=h and offset zero except case18 (+.3 rad). Every impressed source starts at t=0',
            geometry='Control-loop centers (R,0,+/-h), radius a, with a=h=2.5 or1.25 m. The control field is nonaxisymmetric about the storage axis; both fields are computed in the same coordinates without source-field splicing',
            screening='Source order24, DOP853 rtol2e-10/atol2e-11, pulse max step .06*min(1,a/5,h/5) in nominal orbital angle; coast max step .2. The source-resolution time step decreases with loop size',
            acceptance='Clear all continuous control-source retarded times by max(5 us,50*T), before first5km passage and before30us; then reach5km outward by30us. A failed bounded checkpoint is not a claim of permanent trapping',
            ranking='Wpair=(Lext+M)Ipeak^2, with b=.02*a, Lext=mu0*a[ln(8a/b)-2] and exact filament mutual inductance; compare only this pair contribution, not the full apparatus budget',
            outcome_qualification='All18 pulse and coast solves completed. Failed checkpoint cases are retained. No invariant-bound failure was discarded: the one selected full audit passed'),
        attempted_cases=attempts,
        comparison=dict(baseline_control_pair=baseline_row,
            selected_at_fixed_aspect=selected_row, selected_at_fixed_wire_radius=fixed_row,
            baseline_over_selected_stored_energy=reference_energy/energy,
            selected_over_baseline_at_fixed_aspect={key:selected_row[key]/baseline_row[key] for key in compared_keys},
            selected_over_baseline_at_fixed_wire_radius={key:fixed_row[key]/baseline_row[key] for key in compared_keys},
            baseline_radiated_energy_joules=baseline_control['selected_source_feasibility']['radiation']['energy_joules'],
            selected_radiated_energy_joules=radiation['radiation']['energy_joules'],
            interpretation='At fixed b/a, shrinking the wire as well as the loop raises surface field and equivalent hoop stress. Keeping b=.1m tests that tradeoff separately: it changes the aspect to.08 for a1.25m and is not used to rank the18 source cases. Its next-order self-inductance indicator is not a certified remainder bound'),
        selected_control=selected, selected_source_feasibility=radiation,
        selected_conductor=conductors, fixed_radius_conductor=fixed_conductor,
        conductor_clearance=dict(rows=clearance_rows,
            qualification='Distances are to a hypothetical circular conductor envelope around each ideal filament. The computed centerline-source trajectory is unchanged; finite-conductor current redistribution and feed fields have not been modeled. The saved pulse-clear event belongs to the filament sources, not the larger hypothetical conductor-envelope wavefront. Exact midplane symmetry gives a source-distance lower bound globally. The reported attained distance is separately checked along the numerical continuous trajectory',
            fixed_radius_approximation='The explicit fixed-b row permits b/a<=.1, while old default conductor checks still require b/a<=.02. The .08 case uses the published moderate-aspect asymptotic correction as an indicator; no error certificate, material survival, or exact transient circuit inference follows'),
        scope='This is one lower stored-energy source model and outgoing proton trajectory, not a buildable magnet or time machine. Uniform loop currents remain impressed, the trajectory uses ideal filaments, and radiation and quasistatic pair energy are different quantities. Storage-field cross energy/forces, feeds, capacitance, actual terminal voltage, current distribution, materials, support dynamics, heating and energy recovery are not a completed apparatus model. No conductor alteration has been substituted into the Lorentz fields. The exact midplane path lacks an out-of-plane robustness test; no prior-day storage, RF/radiation-reaction combination, injection, capture, deceleration, reunion or human transport is established',
        sources=conductors['sources']+radiation['sources'])


def proton_finite_torus_control_checks(ring, coils, localized_result, screening_cases=None):
    """Extend the fixed selected pulse to two finite torus-surface currents, without retuning."""
    import copy
    import math
    import numpy as np

    baseline = localized_result['selected_control']
    baseline_parameters = baseline['parameters']
    a = baseline_parameters['control_loop_radius_metres']
    centers = baseline_parameters['control_loop_centers_metres']
    h = abs(centers[1][2]-centers[0][2])/2
    fixed = {key:baseline_parameters[key] for key in (
        'control_loop_radius_metres', 'peak_impressed_current_per_loop_amperes',
        'source_pulse_duration_seconds', 'starting_azimuth_radians')}
    fixed['control_loop_half_separation_metres'] = h
    requests = [dict(fixed, control_torus_minor_radius_metres=b) for b in (.025, .1)]
    screens = ([proton_retarded_control_checks(ring, coils, parameters, screening=True)
                for parameters in requests] if screening_cases is None else copy.deepcopy(screening_cases))
    assert len(screens) == 2
    cases = []
    for parameters, screen in zip(requests, screens):
        assert screen['parameters'] == parameters
        b = parameters['control_torus_minor_radius_metres']
        current = parameters['peak_impressed_current_per_loop_amperes']
        # Independent current normalization and dipole moment of the defined surface distribution.
        moment_errors = []
        for order in (8, 16, 32):
            theta = 2*math.pi*np.arange(order)/order
            radii = a+b*np.cos(theta)
            assert np.min(radii)>0 and b<h
            moment = math.pi*current*float(np.mean(radii*radii))
            expected_moment = math.pi*current*(a*a+b*b/2)
            error = abs(moment/expected_moment-1)
            assert error<1e-14
            moment_errors.append(error)
        case = dict(minor_radius_metres=b, minor_over_major_radius=b/a,
            screening=screen, outcome=screen['outcome'],
            source_current=dict(
                signed_peak_surface_current_amperes_per_metre=current/(2*math.pi*b),
                signed_peak_current_per_meridional_cut_amperes=current,
                surface_area_per_torus_square_metres=4*math.pi**2*a*b,
                signed_peak_pair_dipole_moment_amperes_square_metres=
                    2*math.pi*current*(a*a+b*b/2),
                quadrature_current_moment_relative_errors=moment_errors),
            control=None)
        if screen['outcome'] == 'outgoing checkpoint reached':
            control = proton_retarded_control_checks(ring, coils, parameters, surface_poloidal_order=16)
            assert control['parameters']['control_torus_minor_radius_metres'] == b
            assert control['pulse_clear']['lab_seconds']<control['endpoint']['lab_seconds']<30e-6
            assert control['pulse_clear']['radius_metres']<5000
            assert control['endpoint']['radial_momentum_over_original_p0']>0
            assert min(float(control['post_pulse_escape_bound'][key]) for key in (
                'lower_constraint_margin', 'upper_constraint_margin_lower_bound'))>0
            assert abs(control['verification']['continuous_trajectory_source_distance'][
                'minimum_distance_metres']-(h-b))<1e-12
            case['control'] = control
            case['outcome'] = 'fully audited outgoing checkpoint and no later midplane radial barrier'
            case['comparison_to_filament'] = dict(
                pulse_clear_lab_seconds_difference=
                    control['pulse_clear']['lab_seconds']-baseline['pulse_clear']['lab_seconds'],
                endpoint_lab_seconds_difference=
                    control['endpoint']['lab_seconds']-baseline['endpoint']['lab_seconds'],
                endpoint_proper_seconds_difference=
                    control['endpoint']['proton_proper_seconds']-baseline['endpoint']['proton_proper_seconds'],
                endpoint_position_difference_metres=float(np.linalg.norm(
                    np.array(control['endpoint']['position_metres'])-
                    np.array(baseline['endpoint']['position_metres']))),
                post_pulse_gamma_difference=
                    control['pulse_clear']['gamma']-baseline['pulse_clear']['gamma'],
                post_pulse_canonical_ell_difference=
                    control['pulse_clear']['original_source_canonical_ell']-
                    baseline['pulse_clear']['original_source_canonical_ell'],
                interpretation='The finite source changes the force history and endpoint direction. '
                    'Differences from the filament trajectory are model changes, not numerical error '
                    'estimates. The positive escape margin is not a field or timing tolerance')
        cases.append(case)
    outcomes = {}
    for case in cases:
        outcomes[case['outcome']] = outcomes.get(case['outcome'], 0)+1
    return dict(
        model='Two fixed finite torus-surface impressed control currents acting on the same initial proton, with the selected filament pulse parameters held unchanged',
        fixed_parameters=fixed, torus_center_positions_metres=centers,
        tested_minor_radii_metres=[.025, .1], retuning_performed=False,
        source_definition=dict(
            surface='x=R+(a+b*cos(theta))*cos(phi), y=(a+b*cos(theta))*sin(phi), z=+/-h+b*sin(theta)',
            area_and_current='dS=b*(a+b*cos(theta))*dtheta*dphi; K=I(t)/(2*pi*b)*e_phi. Hence K*dS=[I(t)*dtheta/(2*pi)]*dl, dl=(a+b*cos(theta))*e_phi*dphi; no additional radius factor multiplies the ring current',
            continuity='Only the azimuthal component is nonzero and it is independent of phi; the surface and distributional three-dimensional current divergence vanish. rho_charge=Phi=0 is consistent with the specified impressed closed source',
            waveform='Each torus carries Ipeak*sin(pi*t/T)^4 on0<t<T, zero otherwise. Source positions are static. The current is prescribed throughout the surface; causal distributed generation is not solved',
            exact_distances='For either torus centered at zc: dmin=abs(sqrt((r_local-a)^2+(z-zc)^2)-b), dmax=sqrt((r_local+a)^2+(z-zc)^2)+b. Min/max over both tori give the exact continuous source extrema, independent of quadrature nodes',
            reference_radius='a is the torus major radius and b its minor radius; these are control-source dimensions about the displaced axis, not the storage-orbit radius'),
        numerical_protocol=dict(
            screening='Exactly two unchanged-pulse cases; azimuth24/poloidal8, DOP853 rtol2e-10/atol2e-11. Retain bounded failures without retuning',
            full_refinements='Separate azimuth32->64 at poloidal16, then poloidal16->32 at azimuth64, then halve the time step at azimuth64/poloidal32. Independent source check uses Cartesian azimuth integration and Gauss rather than periodic poloidal quadrature',
            checks='Vacuum Maxwell and potential derivatives, actual electric-work ledger, exact continuous source clearance, zero-current circular motion, post-clearance energy/canonical momentum, and the spatially continuous radial escape bound for each actual endpoint',
            baseline_scope='The original22 result blocks remain unchanged. The source field, wavefront event and trajectory are recomputed here; a geometric envelope is no longer substituted for the field source',
            finite_horizons='Continuous source clearance by max(5us,50T), before first5km crossing and30us; then5km outward by30us. A bounded miss would not prove permanent trapping'),
        cases=cases, outcomes=outcomes,
        scope='This is a specified finite surface-current source, not a conductor-response or drive design. Uniform poloidal current weighting need not equal an induced skin-current distribution. No resistivity, current redistribution, support motion, mechanical stress, heating, feeds, source impedance or radiation reaction is solved. Earlier thin-wire energy/hoop estimates cannot be promoted to an exact cost or survival result for this new current distribution. The original storage coils remain ideal filaments; midplane robustness, prior storage, RF, injection, capture, deceleration, reunion and human transport remain outside this calculation',
        sources=baseline['sources'],
        derivation_attribution='The surface-current reduction, exact torus distance extrema, fixed-pulse two-case checks and continuous retarded Lorentz trajectories are our calculations; the retarded-potential framework is the cited primary source')



def proton_torus_launch_perturbation_checks(ring, coils, torus,
                                          screening_cases=None, audited_cases=None):
    """Bounded launch sensitivities in the fixed finite source, with two refined trajectories.

    Cached records may be supplied for focused, deterministic aggregation. Omitting them
    reruns only this launch study, not earlier source searches or storage calculations.
    """
    import copy
    import math
    import numpy as np

    fixed = dict(torus['fixed_parameters'], control_torus_minor_radius_metres=.1)
    nominal = next(row['control'] for row in torus['cases']
                   if row['minor_radius_metres'] == .1)
    omega = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    combined = dict(radial_offset_metres=.001, vertical_offset_metres=.001,
                    radial_direction_radians=1e-5, vertical_direction_radians=1e-5)
    requests = [('baseline', 0., {}), ('combined_3d', 0., combined)]
    for dt in (-.1e-9, .1e-9, -1e-9, 1e-9, -25e-9, 25e-9, 50e-9, 37.5e-9):
        requests.append((f'timing_{dt}', dt, {}))
    for key, amplitudes in (
            ('radial_offset_metres', (.0005, .001)),
            ('vertical_offset_metres', (.0005, .001)),
            ('radial_direction_radians', (5e-6, 1e-5)),
            ('vertical_direction_radians', (5e-6, 1e-5))):
        for amplitude in amplitudes:
            for sign in (-1, 1):
                requests.append((f'{key}_{sign*amplitude}', 0., {key:sign*amplitude}))
    assert len(requests) == 26 and len({row[0] for row in requests}) == 26
    defaults = dict(radial_offset_metres=0., vertical_offset_metres=0.,
                    radial_direction_radians=0., vertical_direction_radians=0.)
    success = 'outgoing cylindrical checkpoint reached'

    reference_endpoint = copy.deepcopy(nominal['endpoint'])
    ref_position = np.array(reference_endpoint['position_metres'])
    ref_momentum = np.array(reference_endpoint['momentum_over_mc'])
    ref_direction = ref_momentum/np.linalg.norm(ref_momentum)
    reference_endpoint['momentum_unit_direction'] = ref_direction.tolist()

    def arrival_comparison(record):
        if record['outcome'] != success:
            return None
        endpoint = record['endpoint']
        delta = np.array(endpoint['position_metres'])-ref_position
        direction = np.array(endpoint['momentum_unit_direction'])
        angle = math.atan2(np.linalg.norm(np.cross(direction, ref_direction)),
                           float(np.dot(direction, ref_direction)))
        return dict(endpoint_position_difference_metres=float(np.linalg.norm(delta)),
            position_delta_metres=delta.tolist(), direction_angle_radians=angle,
            lab_seconds_difference=endpoint['lab_seconds']-reference_endpoint['lab_seconds'],
            proper_seconds_difference=endpoint['proton_proper_seconds']-
                reference_endpoint['proton_proper_seconds'],
            interpretation='Differences at each trajectory crossing of the same cylindrical radius. '
                'They combine prescribed launch changes (zero for the baseline) with finite numerical '
                'error, assessed separately; no receiver acceptance or error bound is specified')

    def normalize(record, case_id, dt, launch, profile, keep_plot=False):
        row = copy.deepcopy(record)
        expected = dict(fixed, starting_azimuth_radians=fixed['starting_azimuth_radians']+omega*dt)
        assert row['parameters'] == expected
        assert row['launch_conditions'] == dict(defaults, **launch)
        assert row['numerical_profile']['name'] == profile
        assert row['horizons']['outgoing_checkpoint_cylindrical_radius_metres'] == 5000.
        assert row['horizons']['outgoing_checkpoint_lab_seconds'] == 30e-6
        assert row['verification']['initial_gamma_error']<5e-14
        row['case_id'] = case_id
        row['phase_equivalent_timing_error_seconds'] = dt
        row['arrival_comparison'] = arrival_comparison(row)
        # Shared mathematical qualifications are kept once at study level.
        row.pop('scope', None)
        row['verification'].pop('clearance_derivation', None)
        row['verification'].pop('clearance_scope', None)
        bound = row['far_three_dimensional_escape_bound']
        if bound is not None:
            if 'certified_outward_continuation' in bound:
                bound['sufficient_outward_continuation_condition_satisfied'] = (
                    bound.pop('certified_outward_continuation'))
            bound.pop('equations', None)
            bound.pop('scope', None)
        if keep_plot:
            if len(row['plot']['lab_microseconds']) == 1001:
                row['plot'] = {key:values[::5] for key,values in row['plot'].items()}
            assert all(len(values) == 201 for values in row['plot'].values())
        else:
            row.pop('plot', None)
        return row

    if screening_cases is None:
        screening_cases = []
        for case_id, dt, launch in requests:
            parameters = dict(fixed, starting_azimuth_radians=fixed['starting_azimuth_radians']+omega*dt)
            row = proton_retarded_control_checks(ring, coils, parameters,
                                                 launch_conditions=launch, launch_resolution='screen')
            row['case_id'] = case_id
            screening_cases.append(row)
    cached = {row['case_id']:row for row in screening_cases}
    assert len(screening_cases) == len(cached) == len(requests)
    assert set(cached) == {row[0] for row in requests}
    screens = [normalize(cached[case_id], case_id, dt, launch, 'screen')
               for case_id,dt,launch in requests]
    lookup = {row['case_id']:row for row in screens}
    outcomes = {}
    for row in screens:
        outcomes[row['outcome']] = outcomes.get(row['outcome'], 0)+1

    selected = ('combined_3d', 'timing_3.75e-08')
    request_map = {case_id:(dt, launch) for case_id,dt,launch in requests}
    if audited_cases is None:
        audited_cases = []
        for case_id in selected:
            dt, launch = request_map[case_id]
            parameters = dict(fixed, starting_azimuth_radians=fixed['starting_azimuth_radians']+omega*dt)
            for profile in ('main', 'source', 'fine'):
                row = proton_retarded_control_checks(ring, coils, parameters,
                    launch_conditions=launch, launch_resolution=profile)
                row['case_id'] = case_id
                audited_cases.append(row)
    # Accept flattened cached records or the finalized grouped schema for focused reproduction.
    if audited_cases and 'main' in audited_cases[0]:
        audited_cases = [group[key] for group in audited_cases
                         for key in ('main', 'source_refined', 'fine')]
    audit_lookup = {(row['case_id'], row['numerical_profile']['name']):row for row in audited_cases}
    assert len(audited_cases) == len(audit_lookup) == 6
    assert set(audit_lookup) == {(case_id,profile) for case_id in selected
                                for profile in ('main','source','fine')}

    def numerical_difference(first, second):
        result = {}
        for event in ('pulse_clear', 'endpoint'):
            left, right = first[event], second[event]
            result[event] = dict(
                position_metres=float(np.linalg.norm(np.array(left['position_metres'])-
                                                     np.array(right['position_metres']))),
                momentum_over_mc=float(np.linalg.norm(np.array(left['momentum_over_mc'])-
                                                       np.array(right['momentum_over_mc']))),
                lab_seconds=abs(left['lab_seconds']-right['lab_seconds']),
                proper_seconds=abs(left['proton_proper_seconds']-right['proton_proper_seconds']),
                gamma=abs(left['gamma']-right['gamma']),
                canonical_ell=abs(left['original_source_canonical_ell']-
                                  right['original_source_canonical_ell']))
        return result

    audits = []
    for case_id in selected:
        dt,launch = request_map[case_id]
        records = [normalize(audit_lookup[(case_id,profile)],case_id,dt,launch,profile,
                             keep_plot=profile == 'fine') for profile in ('main','source','fine')]
        for row in records:
            assert row['outcome'] == success
            assert not row['source_surface_intersection']
            assert row['verification']['maximum_gamma_minus_initial_minus_work_residual']<1e-9
            assert row['verification']['continuous_control_surface_clearance_lower_bound_metres']>0
            assert row['verification']['continuous_storage_filament_clearance_lower_bound_metres']>0
            assert row['far_three_dimensional_escape_bound'][
                'sufficient_outward_continuation_condition_satisfied']
        source_change = numerical_difference(records[0], records[1])
        step_change = numerical_difference(records[1], records[2])
        for differences in (source_change, step_change):
            for event in differences.values():
                assert event['position_metres']<1e-5
                assert event['lab_seconds']<1e-14
                assert event['gamma']<1e-9
        audits.append(dict(case_id=case_id, phase_equivalent_timing_error_seconds=dt,
            main=records[0], source_refined=records[1], fine=records[2],
            convergence=dict(source_quadrature_change=source_change,
                timestep_and_tolerance_change=step_change,
                screen_to_fine=numerical_difference(lookup[case_id],records[2]),
                scope='Numerical refinements at identical physical launch/source parameters. '
                    'The 10 micrometre, 1e-14 second and 1e-9 gamma checks are numerical validation '
                    'criteria, not physical alignment tolerances or rigorous error enclosures')))

    def central_difference(plus, minus, step):
        assert plus['outcome'] == minus['outcome'] == success
        pp,mm = plus['endpoint'],minus['endpoint']
        return dict(position_metres_per_input_unit=((np.array(pp['position_metres'])-
                np.array(mm['position_metres']))/(2*step)).tolist(),
            direction_per_input_unit=((np.array(pp['momentum_unit_direction'])-
                np.array(mm['momentum_unit_direction']))/(2*step)).tolist(),
            lab_seconds_per_input_unit=(pp['lab_seconds']-mm['lab_seconds'])/(2*step),
            proper_seconds_per_input_unit=(pp['proton_proper_seconds']-
                mm['proton_proper_seconds'])/(2*step),
            gamma_per_input_unit=(pp['gamma']-mm['gamma'])/(2*step),
            canonical_ell_per_input_unit=(pp['original_source_canonical_ell']-
                mm['original_source_canonical_ell'])/(2*step))

    sensitivities = []
    for key,unit,steps in (
            ('phase_equivalent_timing_error_seconds','second',(.1e-9,1e-9)),
            ('radial_offset_metres','metre',(.0005,.001)),
            ('vertical_offset_metres','metre',(.0005,.001)),
            ('radial_direction_radians','radian',(5e-6,1e-5)),
            ('vertical_direction_radians','radian',(5e-6,1e-5))):
        levels = []
        for step in steps:
            prefix = 'timing' if key.startswith('phase') else key
            value = central_difference(lookup[f'{prefix}_{step}'],lookup[f'{prefix}_{-step}'],step)
            levels.append(dict(symmetric_step_in_input_units=step,central_difference=value))
        small = np.array(levels[0]['central_difference']['position_metres_per_input_unit'])
        large = np.array(levels[1]['central_difference']['position_metres_per_input_unit'])
        sensitivities.append(dict(input=key,input_unit=unit,levels=levels,
            position_derivative_relative_step_change=float(np.linalg.norm(small-large)/
                                                           max(np.linalg.norm(small),1e-30)),
            scope='Two-amplitude central differences of the coarse trajectory solve. '
                'These estimate local arrival sensitivity; they do not prove a tolerance region '
                'or establish every derivative component above numerical noise'))
    reflection_errors = []
    for key,steps in (('vertical_offset_metres',(.0005,.001)),
                      ('vertical_direction_radians',(5e-6,1e-5))):
        for step in steps:
            positive = np.array(lookup[f'{key}_{step}']['endpoint']['position_metres'])
            negative = np.array(lookup[f'{key}_{-step}']['endpoint']['position_metres'])
            reflection_errors.append(float(np.linalg.norm(positive-negative*np.array([1,1,-1]))))
    assert max(reflection_errors)<1e-6
    passing = [row for row in screens if row['outcome'] == success]
    for row in passing:
        assert row['far_three_dimensional_escape_bound'][
            'sufficient_outward_continuation_condition_satisfied']

    return dict(
        model='Bounded timing and three-dimensional launch variations in the unchanged finite torus control and original storage fields',
        fixed_source_parameters=fixed, source_retuning_performed=False,
        reference=dict(source='finite_torus_control.cases minor_radius_metres=.1, previous fully audited trajectory',
            crossing='First outward crossing of cylindrical radius sqrt(x^2+y^2)=5000 metres after complete source clearance',
            endpoint=reference_endpoint,
            scope='A reference arrival event on a cylinder, not a target aperture, field-free detector or captured endpoint'),
        launch_definition=dict(
            position='x0=(R+dr)*e_r(phi0)+dz*e_z',
            momentum='p_r=p0*sin(ar)*cos(az); p_phi=-p0*cos(ar)*cos(az); p_z=p0*sin(az), p0=m*c*sqrt(10^2-1)',
            timing='phi0=phi_nominal+omega0*delta_t. Positive delta_t means later unperturbed orbital arrival at a fixed azimuth while the source pulse stays fixed at t=0. This is not an explicitly simulated source-trigger delay; resetting after a positive trigger delay gives the opposite phase shift',
            nominal_angular_frequency_radians_per_second=omega,
            gamma='All launch perturbations preserve initial gamma=10; pulse electric work may subsequently change it'),
        numerical_protocol=dict(screening_case_count=26,refined_case_count=2,
            refined_runs_per_case=3,persisted_trajectory_evaluations=32,
            count_scope='26 screening cases and six refined records; repeated verification and legacy regression do not add physical cases',
            screening='Azimuth24/poloidal8, DOP853 rtol2e-10/atol2e-11. All 26 defined cases are retained',
            refinements='For the fixed combined 3D launch and +37.5ns case: main azimuth32/poloidal16; source azimuth64/poloidal32 at the same time settings; fine retains64/32, halves maximum steps and tightens tolerances',
            source_checks='The retarded torus field kernel and source parameters are unchanged from finite_torus_control. Its Maxwell/potential checks remain applicable. New trajectory checks use actual electric work, all-source causal clearance, signed continuous separation, and full three-dimensional post-clearance invariants',
            horizons='Source clearance by64us and before the cylindrical checkpoint; outgoing cylindrical5km checkpoint by30us. A time-horizon miss is not proof of permanent trapping',
            mathematical_scope='Bounds use numerical nodal/end states; no outward-rounded interval integration or universal perturbation cell is claimed'),
        screening_cases=screens,outcomes=outcomes,audited_cases=audits,sensitivities=sensitivities,
        checks=dict(maximum_vertical_reflection_position_error_metres=max(reflection_errors),
            maximum_screen_initial_gamma_error=max(row['verification']['initial_gamma_error'] for row in screens),
            minimum_successful_three_dimensional_escape_margin=min(row['far_three_dimensional_escape_bound'][
                'positive_radial_direction_squared_margin'] for row in passing)),
        sampled_timing_transition=dict(passing_phase_equivalent_seconds=37.5e-9,
            missed_phase_equivalent_seconds=50e-9, missed_outcome=lookup['timing_5e-08']['outcome'],
            scope='A passing and a missed sample bound a change in this finite-horizon test. '
                'No exact threshold, monotonicity or failure of all larger offsets is proved'),
        continuous_clearance=dict(
            distance='Signed exterior control-surface separation is min_zc(sqrt((r_local-a)^2+(z-zc)^2)-b). Storage distance is the minimum exact filament distance',
            bound='Both distances are 1-Lipschitz. On each sampled interval, d(t)>=min(d_i,d_j)-vmax*Delta t/2. Loose intervals are subdivided using existing dense output; a positive minimum excludes crossing conditional on the numerical nodal positions',
            causal_clearance='t-max_source_distance/c=T, with exact torus maximum sqrt((r_local+a)^2+(z-zc)^2)+b and actual z. Its derivative is at least1-|v|/c>0, so the control signal cannot re-enter this future subluminal path'),
        three_dimensional_escape_bound=dict(
            equations='Storage sources lie inside s=sqrt(a_storage^2+h_storage^2); |B|<=C/(r-s)^2 with C=mu0*abs(Istorage)*astorage for the pair. After clearance p is constant. For spherical w=rhat dot vhat>0, d(w^2)/dr>=-2*abs(q)*|B|/p, so w0^2-2*abs(q)*C/[p*(r0-s)]>0 prevents any later spherical radial turn',
            scope='Uses spherical radius and outward direction evaluated at a cylindrical checkpoint. No midplane barrier or exact z=0 assumption is used. The sufficient condition is evaluated on numerical endpoint data and establishes continued escape in the prescribed static source, not delivery'),
        scope='An unchanged prescribed divergence-free surface current produces these test-particle trajectories. Numerical source clearance and escape do not show receiver alignment, capture, deceleration, reunion, material response or practical source generation. Initial displacements/directions are prescribed launch data, not a modeled injector. Finite differences are local sampled sensitivities; surviving samples do not prove generic stability',
        sources=torus['sources'],
        derivation_attribution='Launch normalization, signed-distance interval bounds, full-3D far-field escape inequality and sampled transport sensitivities are our calculations using the previously cited retarded Maxwell and Lorentz framework')


def _electrostatic_torus_pair_kernel(inner, outer, minor, charge, permittivity,
                                     azimuth_order, poloidal_order, poloidal_weighting=None):
    """Static +/-Q torus surfaces in local coordinates, sharing the z axis and z=0 plane.

    Equal theta/phi charge measure dQ=+/-Q*dtheta*dphi/(4*pi^2) implies
    sigma(theta)=+/-Q/[4*pi^2*b*(a+b*cos(theta))], not uniform charge per area.
    These are prescribed fixed charges, not solved equipotential conductors.
    Optional 'cosine_squared' multiplies the normalized angular measure by
    2*cos(theta)^2; the default arithmetic is unchanged.
    """
    import math
    import numpy as np
    assert 0<minor<inner<outer and outer-inner>2*minor
    assert azimuth_order>=8 and poloidal_order>=4
    theta,phi = np.meshgrid(2*math.pi*np.arange(poloidal_order)/poloidal_order,
                           2*math.pi*np.arange(azimuth_order)/azimuth_order,indexing='ij')
    er = np.column_stack((np.cos(phi).ravel(),np.sin(phi).ravel(),np.zeros(phi.size)))
    axial = np.column_stack((np.zeros(phi.size),np.zeros(phi.size),minor*np.sin(theta).ravel()))
    first = (inner+minor*np.cos(theta).ravel())[:,None]*er+axial
    second = first+(outer-inner)*er
    delta = second-first
    coefficient = charge/(4*math.pi*permittivity)
    angular_weights=None
    if poloidal_weighting is not None:
        if poloidal_weighting!='cosine_squared':raise ValueError('Unknown receiver angular measure')
        angular_weights=2*np.cos(theta).ravel()**2
        assert abs(float(np.mean(angular_weights))-1)<1e-14

    def field(position):
        point = np.asarray(position)
        r1,r2 = point-first,point-second
        d1sq,d2sq = np.sum(r1*r1,axis=1),np.sum(r2*r2,axis=1)
        d1,d2 = np.sqrt(d1sq),np.sqrt(d2sq)
        # Rationalize paired inverse distances rather than subtract nearly equal monopoles.
        difference_sq = np.sum(delta*(second+first-2*point),axis=1)
        inverse_difference = difference_sq/(d1*d2*(d1+d2))
        inverse_cube_difference = difference_sq*(d1sq+d1*d2+d2sq)/(
            (d1+d2)*d1**3*d2**3)
        if angular_weights is not None:
            potential=coefficient*np.mean(angular_weights*inverse_difference)
            electric=coefficient*np.mean(angular_weights[:,None]*(delta/d1[:,None]**3+
                                        r2*inverse_cube_difference[:,None]),axis=0)
            return potential,electric
        potential = coefficient*np.mean(inverse_difference)
        electric = coefficient*np.mean(delta/d1[:,None]**3+
                                        r2*inverse_cube_difference[:,None],axis=0)
        return potential,electric

    return field


def _electrostatic_receiver_definition(ring, coils, torus):
    """Fix one receiver from the previous arrival; its design charge is never retuned by a solve."""
    import math
    import numpy as np
    from scipy.constants import c,e,m_p,epsilon_0
    nominal = next(row['control'] for row in torus['cases'] if row['minor_radius_metres'] == .1)
    point = np.array(nominal['endpoint']['position_metres'])
    direction = np.array(nominal['endpoint']['momentum_over_mc'])
    direction /= np.linalg.norm(direction)
    up = np.array([0.,0.,1.])
    first = up-np.dot(up,direction)*direction
    first /= np.linalg.norm(first)
    frame = np.column_stack((first,np.cross(direction,first),direction))
    assert np.max(np.abs(frame.T@frame-np.eye(3)))<1e-14 and np.linalg.det(frame)>0
    pulse = dict(torus['fixed_parameters'],control_torus_minor_radius_metres=.1)
    radius = coils['parameters']['reference_orbit_radius_metres']
    phase = pulse['starting_azimuth_radians']
    initial = radius*np.array([math.cos(phase),math.sin(phase),0.])
    local_initial = frame.T@(initial-point)
    unit = _electrostatic_torus_pair_kernel(10.,20.,.1,1.,epsilon_0,64,32)
    at_center = float(unit(np.zeros(3))[0])
    at_launch = float(unit(local_initial)[0])
    design_gamma = 1.01
    voltage = (nominal['endpoint']['gamma']-design_gamma)*m_p*c*c/e
    charge = voltage/(at_center-at_launch)
    return dict(center_metres=point.tolist(),frame_columns=frame.tolist(),
        axis_unit_direction=direction.tolist(),inner_major_radius_metres=10.,
        outer_major_radius_metres=20.,minor_radius_metres=.1,
        inner_total_charge_coulombs=charge,outer_total_charge_coulombs=-charge,
        vacuum_permittivity_SI=epsilon_0,design_gamma=design_gamma,
        design_potential_rise_volts=voltage,
        central_potential_volts=charge*at_center,initial_potential_volts=charge*at_launch,
        reference_incoming_gamma=nominal['endpoint']['gamma'],
        pulse_parameters=pulse)


def proton_electrostatic_receiver_trace(ring, coils, torus, profile='main', receiver_design=None):
    """One continuous launch-to-receiver Lorentz trajectory in all prescribed sources."""
    import math
    import numpy as np
    from scipy.constants import c,e,m_p
    from scipy.integrate import solve_ivp
    from scipy.optimize import minimize_scalar
    from scipy.special import ellipk,ellipe

    design = (_electrostatic_receiver_definition(ring,coils,torus)
              if receiver_design is None else dict(receiver_design))
    center,frame = np.array(design['center_metres']),np.array(design['frame_columns'])
    pulse_parameters = design['pulse_parameters']
    source = coils['parameters']
    radius,field0,gamma0 = (source['reference_orbit_radius_metres'],
                           source['reference_Bz_tesla'],source['reference_gamma'])
    mu,current = source['vacuum_permeability_SI'],source['ampere_turns_per_loop']
    storage_a,storage_h = source['coil_radius_metres'],source['coil_z_metres'][1]
    ca,ch,cb = (pulse_parameters['control_loop_radius_metres'],
                pulse_parameters['control_loop_half_separation_metres'],
                pulse_parameters['control_torus_minor_radius_metres'])
    duration = pulse_parameters['source_pulse_duration_seconds']
    amplitude = pulse_parameters['peak_impressed_current_per_loop_amperes']
    phase = pulse_parameters['starting_azimuth_radians']
    omega = 2*math.pi*ring['nominal']['revolution_frequency_hz']
    p0 = math.sqrt(gamma0*gamma0-1)
    control_center = np.array([radius,0.,0.])
    profiles = {'main':(32,16,2e-12,2e-13,.00375,.15),
                'source':(64,32,2e-12,2e-13,.00375,.15),
                'fine':(64,32,5e-13,5e-14,.001875,.075)}
    assert profile in profiles
    azimuth,poloidal,rtol,atol,pulse_step,coast_step = profiles[profile]
    theta = 2*math.pi*np.arange(poloidal)/poloidal
    pairs = np.column_stack((ca+cb*np.cos(theta),ch+cb*np.sin(theta),
                             np.full(poloidal,1./poloidal)))
    control = _retarded_ring_family_control(control_center,mu,c,duration,amplitude,pairs,azimuth)
    receiver = _electrostatic_torus_pair_kernel(10.,20.,.1,
        design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI'],azimuth,poloidal,
        design.get('poloidal_weighting'))
    initial = np.array([math.cos(phase),math.sin(phase),0.,
                       p0*math.sin(phase),-p0*math.cos(phase),0.,0.,0.,0.])
    initial_potential = float(receiver(frame.T@(initial[:3]*radius-center))[0])

    def local(position):
        return frame.T@(position-center)

    def background(position):
        rr = math.hypot(*position[:2])
        _,br,bz = _finite_coil_loop_sum(rr,position[2],storage_a,storage_h,current,
                                       math.sqrt,ellipk,ellipe,math.pi,mu)
        return np.array([br*position[0]/rr,br*position[1]/rr,bz])

    def source_gaps(position):
        rr = math.hypot(position[0]-radius,position[1])
        control_gap = min(math.hypot(rr-ca,position[2]-height)-cb for height in (-ch,ch))
        coordinate = local(position)
        transverse = math.hypot(*coordinate[:2])
        receiver_gaps = [math.hypot(transverse-a,coordinate[2])-.1 for a in (10.,20.)]
        rr_storage = math.hypot(*position[:2])
        storage_gap = min(math.hypot(rr_storage-storage_a,position[2]-height)
                          for height in (-storage_h,storage_h))
        return np.array([control_gap,*receiver_gaps,storage_gap])

    def rhs(cleared):
        def evaluate(s,state):
            position,pp = state[:3]*radius,state[3:6]
            gamma = math.sqrt(1+np.dot(pp,pp))
            _,electric_local = receiver(local(position))
            electric_receiver = frame@electric_local
            if cleared:
                electric_control,extra = np.zeros(3),np.zeros(3)
            else:
                _,electric_control,extra = control(s/omega,position)
            electric = electric_receiver+electric_control
            magnetic = background(position)+extra
            force = gamma0*(electric/(c*field0)+np.cross(pp,magnetic)/(gamma*field0))
            factor = gamma0/(gamma*c*field0)
            return np.concatenate((gamma0/gamma*pp/p0,force,
                [1/gamma,factor*np.dot(pp,electric_control),factor*np.dot(pp,electric_receiver)]))
        return evaluate

    def clear_event(s,state):
        x=state[:3]*radius
        rr=math.hypot(x[0]-radius,x[1])
        dmax=max(math.hypot(rr+ca,x[2]-height) for height in (-ch,ch))+cb
        return s/omega-dmax/c-duration
    clear_event.direction,clear_event.terminal=1,True

    def collision(s,state):
        return float(np.min(source_gaps(state[:3]*radius)[:3]))
    collision.direction,collision.terminal=-1,True

    def plane(s,state):
        return float(local(state[:3]*radius)[2])
    plane.direction,plane.terminal=1,False

    def downstream(s,state):
        return plane(s,state)-200.
    downstream.direction,downstream.terminal=1,True

    def reflected(s,state):
        return plane(s,state)+200.
    reflected.direction,reflected.terminal=-1,True

    assert np.min(source_gaps(initial[:3]*radius))>0
    pulse = solve_ivp(rhs(False),(0,omega*64e-6),initial,events=(clear_event,collision),
        dense_output=True,method='DOP853',rtol=rtol,atol=atol,max_step=pulse_step)
    assert pulse.success and len(pulse.t_events[0]),'Control clearance not achieved; retain failing run'
    assert not len(pulse.t_events[1])
    assert pulse.t[-1]/omega<30e-6
    coast = solve_ivp(rhs(True),(pulse.t[-1],omega*30e-6),pulse.y[:,-1],
        events=(plane,downstream,reflected,collision),dense_output=True,
        method='DOP853',rtol=rtol,atol=atol,max_step=coast_step)
    outcome = ('coast solver failure' if not coast.success else
        'source surface intersection' if len(coast.t_events[3]) else
        'transmitted through receiver and reached downstream plane' if len(coast.t_events[1]) else
        'reflected to upstream plane' if len(coast.t_events[2]) else
        'no downstream or reflected checkpoint by thirty microseconds')
    finish = float(coast.t[-1])
    nodes=np.unique(np.concatenate((pulse.t,coast.t,np.linspace(0,finish,1201))))

    def values(times):
        result=np.empty((9,len(times)))
        before=times<=pulse.t[-1]
        if np.any(before):result[:,before]=pulse.sol(times[before])
        if np.any(~before):result[:,~before]=coast.sol(times[~before])
        return result

    speed_bound=max(c,radius*omega*gamma0/p0)
    assert abs(speed_bound/c-1)<1e-13
    for _ in range(13):
        sampled=values(nodes)
        gaps=np.array([source_gaps(x) for x in sampled[:3].T*radius])
        bounds=np.minimum(gaps[:-1],gaps[1:])-speed_bound*np.diff(nodes)[:,None]/(2*omega)
        unresolved=np.any(bounds<=0,axis=1)
        if not np.any(unresolved) or outcome=='source surface intersection':break
        nodes=np.unique(np.concatenate((nodes,(nodes[:-1][unresolved]+nodes[1:][unresolved])/2)))
    sampled=values(nodes)
    gammas=np.sqrt(1+np.sum(sampled[3:6]**2,axis=0))
    potentials=np.array([float(receiver(local(x))[0]) for x in sampled[:3].T*radius])
    energy_unit=m_p*c*c
    work_residual=gammas-gamma0-sampled[7]-sampled[8]
    potential_residual=sampled[8]+e*(potentials-initial_potential)/energy_unit
    static_residual=gammas-gamma0+e*(potentials-initial_potential)/energy_unit-sampled[7]
    assert np.all(np.isfinite(sampled)) and np.min(gammas)>=1
    assert np.all(np.diff(sampled[6])>=0)
    if coast.success and not len(coast.t_events[3]):
        assert np.min(bounds)>0
        assert max(np.max(np.abs(work_residual)),np.max(np.abs(potential_residual)),
                   np.max(np.abs(static_residual)))<1e-8

    def state_record(s,state):
        position,pp=state[:3]*radius,state[3:6]
        coordinate=local(position)
        potential,electric=receiver(coordinate)
        gamma=math.sqrt(1+float(np.dot(pp,pp)))
        return dict(lab_seconds=float(s/omega),proton_proper_seconds=float(state[6]/omega),
            position_metres=position.tolist(),momentum_over_mc=pp.tolist(),
            receiver_coordinates_metres=coordinate.tolist(),
            receiver_transverse_offset_metres=math.hypot(*coordinate[:2]),
            receiver_axial_momentum_over_mc=float(np.dot(frame[:,2],pp)),
            gamma=gamma,beta=math.sqrt(1-1/(gamma*gamma)),
            kinetic_energy_eV=(gamma-1)*energy_unit/e,
            receiver_potential_volts=float(potential),
            receiver_electric_field_volts_per_metre=float(np.linalg.norm(electric)),
            control_work_joules=float(state[7]*energy_unit),
            receiver_work_joules=float(state[8]*energy_unit),
            gamma_plus_potential_over_mc2=float(gamma+e*potential/energy_unit))

    minimum_index=int(np.argmin(gammas))
    lo,hi=nodes[max(0,minimum_index-1)],nodes[min(len(nodes)-1,minimum_index+1)]
    minimum=minimize_scalar(lambda s:math.sqrt(1+float(np.dot(values(np.array([s]))[3:6,0],
                                                               values(np.array([s]))[3:6,0]))),
        bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
    candidates=[(float(nodes[minimum_index]),sampled[:,minimum_index])]
    if minimum.success:
        candidates.append((float(minimum.x),values(np.array([minimum.x]))[:,0]))
    min_s,min_state=min(candidates,key=lambda row:np.dot(row[1][3:6],row[1][3:6]))
    crossing=(state_record(float(coast.t_events[0][0]),coast.y_events[0][0])
              if len(coast.t_events[0]) else None)
    minimum_record=state_record(min_s,min_state)
    if crossing is not None and crossing['gamma']<minimum_record['gamma']:
        minimum_record=crossing
    near=np.array([local(x)[2] for x in sampled[:3].T*radius])
    eligible=np.where(np.abs(near)<250)[0]
    local_start=nodes[eligible[0]] if len(eligible) else pulse.t[-1]
    event_times=[float(pulse.t[-1]),min_s,*coast.t_events[0],finish]
    plot_times=np.unique(np.concatenate((np.linspace(0,finish,201),np.linspace(local_start,finish,401),event_times)))
    plot_values=values(plot_times)
    plot_local=np.array([local(x) for x in plot_values[:3].T*radius])
    plot_gamma=np.sqrt(1+np.sum(plot_values[3:6]**2,axis=0))
    return dict(profile=profile,design=design,outcome=outcome,
        pulse_clear=state_record(float(pulse.t[-1]),pulse.y[:,-1]),
        receiver_plane_crossing=crossing,
        minimum_gamma_state=minimum_record,
        endpoint=state_record(finish,coast.y[:,-1]),
        initial=state_record(0.,initial),
        numerical_profile=dict(azimuth_order=azimuth,poloidal_order=poloidal,rtol=rtol,atol=atol,
            pulse_maximum_dimensionless_step=pulse_step,coast_maximum_dimensionless_step=coast_step,
            pulse_rhs_evaluations=pulse.nfev,coast_rhs_evaluations=coast.nfev),
        verification=dict(maximum_total_work_residual_gamma=float(np.max(np.abs(work_residual))),
            maximum_receiver_work_plus_potential_residual_gamma=float(np.max(np.abs(potential_residual))),
            maximum_static_energy_minus_control_work_residual_gamma=float(np.max(np.abs(static_residual))),
            continuous_source_gap_lower_bounds_metres=np.min(bounds,axis=0).tolist(),
            source_gap_order=['control torus pair','inner receiver torus','outer receiver torus','storage coil pair'],
            gap_scope='1-Lipschitz signed surface/filament distances with vmax*Delta t/2 endpoint-interval bounds; conditional on numerical nodal positions, not interval-certified integration',
            pulse_clear_event_residual_seconds=float(clear_event(pulse.t[-1],pulse.y[:,-1])),
            all_receiver_fields_present_from_initial_time=True),
        plot=dict(lab_microseconds=(plot_times/omega*1e6).tolist(),
            receiver_axial_metres=plot_local[:,2].tolist(),
            receiver_transverse_metres=np.hypot(plot_local[:,0],plot_local[:,1]).tolist(),
            gamma=plot_gamma.tolist(),
            kinetic_energy_MeV=((plot_gamma-1)*energy_unit/e/1e6).tolist(),
            proton_proper_microseconds=(plot_values[6]/omega*1e6).tolist(),
            x_metres=(plot_values[0]*radius).tolist(),y_metres=(plot_values[1]*radius).tolist(),
            z_metres=(plot_values[2]*radius).tolist()),
        scope='Full initial launch uses the receiver field as well as the unchanged retarded control and static storage fields. Receiver charge was designed from an earlier arrival estimate and not retuned. A slowed plane crossing or reflection is not permanent capture; no energy sink or source switching is modeled')


def proton_electrostatic_receiver_checks(ring, coils, torus, launch_study, trajectories=None):
    """Audit one fixed-charge finite receiver and three resolutions of its full incoming route."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c,e,m_p

    design=_electrostatic_receiver_definition(ring,coils,torus)
    charge,eps=design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI']
    coefficient=charge/(4*math.pi*eps)
    field_scale=abs(coefficient)/10**2
    coarse=_electrostatic_torus_pair_kernel(10,20,.1,charge,eps,32,16)
    fine_field=_electrostatic_torus_pair_kernel(10,20,.1,charge,eps,64,32)
    center=np.array(design['center_metres']);frame=np.array(design['frame_columns'])
    radius=coils['parameters']['reference_orbit_radius_metres']
    phase=design['pulse_parameters']['starting_azimuth_radians']
    initial_position=radius*np.array([math.cos(phase),math.sin(phase),0.])
    initial_local=frame.T@(initial_position-center)
    probes=[np.zeros(3),np.array([0.,0.,10.]),np.array([2.,1.,-30.]),
            np.array([.1,.2,-2.]),np.array([60.,40.,50.]),initial_local]
    probe_rows=[]
    for point in probes:
        potential,electric=fine_field(point)
        gradient=np.empty(3);jacobian=np.empty((3,3))
        for j in range(3):
            complex_point=point.astype(complex)
            complex_point[j]+=1j*1e-25
            pot,vec=fine_field(complex_point)
            gradient[j]=pot.imag/1e-25;jacobian[:,j]=vec.imag/1e-25
        low_pot,low_field=coarse(point)
        row=dict(receiver_coordinates_metres=point.tolist(),
            potential_volts=float(potential),electric_magnitude_volts_per_metre=float(np.linalg.norm(electric)),
            potential_gradient_error_volts_per_metre=float(np.linalg.norm(electric+gradient)),
            gradient_error_over_characteristic_source_field=float(np.linalg.norm(electric+gradient)/field_scale),
            divergence_relative_to_derivative_scale=float(abs(np.trace(jacobian))/max(1,np.linalg.norm(jacobian))),
            curl_relative_to_derivative_scale=float(np.max(np.abs(jacobian-jacobian.T))/max(1,np.linalg.norm(jacobian))),
            source_refinement_potential_difference_volts=float(abs(low_pot-potential)),
            source_refinement_field_difference_volts_per_metre=float(np.linalg.norm(low_field-electric)))
        assert max(row['gradient_error_over_characteristic_source_field'],
                   row['divergence_relative_to_derivative_scale'],
                   row['curl_relative_to_derivative_scale'])<1e-12
        assert row['source_refinement_field_difference_volts_per_metre']<1e-4
        probe_rows.append(row)
    zero=_electrostatic_torus_pair_kernel(10,20,.1,0.,eps,32,16)
    assert zero(probes[1])[0] == 0 and np.max(np.abs(zero(probes[1])[1])) == 0
    with mp.workdps(50):
        qh,eh,bh=map(lambda value:mp.mpf(str(value)),(charge,eps,.1))
        kh=qh/(4*mp.pi*eh)
        def axis_potential(z):
            zh=mp.mpf(str(z))
            def integrand(theta):
                ri,ro=10+bh*mp.cos(theta),20+bh*mp.cos(theta)
                dz=zh-bh*mp.sin(theta)
                di,do=mp.sqrt(ri*ri+dz*dz),mp.sqrt(ro*ro+dz*dz)
                return (ro*ro-ri*ri)/(di*do*(di+do))
            return kh*mp.quad(integrand,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi)
        axis_checks=[]
        for z in (0.,10.,1000.,1e6):
            expected=axis_potential(z)
            actual=float(fine_field(np.array([0.,0.,z]))[0])
            relative=abs(actual/float(expected)-1)
            assert relative<2e-14
            axis_checks.append(dict(axial_metres=z,potential_volts=actual,
                independent_50_digit_axis_integral=mp.nstr(expected,30),
                relative_difference=relative))
    at_origin=probes[0].astype(complex)
    electric_jacobian=np.empty((3,3))
    for j in range(3):
        point=at_origin.copy();point[j]+=1j*1e-25
        electric_jacobian[:,j]=fine_field(point)[1].imag/1e-25
    hessian=-electric_jacobian
    assert hessian[2,2]<0 and min(hessian[0,0],hessian[1,1])>0
    assert abs(np.trace(hessian))/np.linalg.norm(hessian)<1e-13
    surfaces=[]
    for a,sign in ((10.,1),(20.,-1)):
        density_peak=charge/(4*math.pi**2*.1*(a-.1))
        surfaces.append(dict(major_radius_metres=a,minor_radius_metres=.1,
            total_charge_coulombs=sign*charge,
            maximum_absolute_prescribed_surface_charge_coulombs_per_square_metre=density_peak,
            maximum_absolute_normal_electric_field_jump_volts_per_metre=density_peak/eps,
            jump_scope='The jump across the prescribed charged sheet is sigma/epsilon0. This is not the electric field on either side or a conductor voltage specification'))
    normalization=[]
    for order in (8,16,32):
        theta=2*math.pi*np.arange(order)/order
        for a,sign in ((10.,1),(20.,-1)):
            density=sign*charge/(4*math.pi**2*.1*(a+.1*np.cos(theta)))
            integrated=4*math.pi**2*float(np.mean(density*.1*(a+.1*np.cos(theta))))
            normalization.append(abs(integrated/(sign*charge)-1))
    assert max(normalization)<1e-14

    if trajectories is None:
        trajectories=[proton_electrostatic_receiver_trace(ring,coils,torus,profile)
                      for profile in ('main','source','fine')]
    if isinstance(trajectories,dict):
        trajectories=[trajectories[key] for key in ('main','source_refined','fine')]
    traces={row['profile']:copy.deepcopy(row) for row in trajectories}
    assert set(traces)=={'main','source','fine'} and len(trajectories)==3
    for profile,row in traces.items():
        assert row['design']==design
        assert row['verification']['all_receiver_fields_present_from_initial_time']
        crossing=row['receiver_plane_crossing']
        if crossing is not None and crossing['gamma']<row['minimum_gamma_state']['gamma']:
            row['minimum_gamma_state']=copy.deepcopy(crossing)
        row['minimum_gamma_scope']='Lowest gamma among evaluated sample/local-search/crossing candidates, not a certified global minimum; no optimizer-convergence flag was retained in the original cached records'
        if profile=='fine':
            plot=row['plot']
            by_time={t:{key:values[i] for key,values in plot.items()}
                     for i,t in enumerate(plot['lab_microseconds'])}
            for event in ('pulse_clear','receiver_plane_crossing','minimum_gamma_state','endpoint'):
                state=row[event]
                if state is None:continue
                time=state['lab_seconds']*1e6
                by_time[time]=dict(lab_microseconds=time,
                    receiver_axial_metres=state['receiver_coordinates_metres'][2],
                    receiver_transverse_metres=state['receiver_transverse_offset_metres'],
                    gamma=state['gamma'],kinetic_energy_MeV=state['kinetic_energy_eV']/1e6,
                    proton_proper_microseconds=state['proton_proper_seconds']*1e6,
                    x_metres=state['position_metres'][0],y_metres=state['position_metres'][1],
                    z_metres=state['position_metres'][2])
            row['plot']={key:[by_time[t][key] for t in sorted(by_time)] for key in plot}
        if profile!='fine':row.pop('plot',None)
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
        for key in ('maximum_total_work_residual_gamma',
                    'maximum_receiver_work_plus_potential_residual_gamma',
                    'maximum_static_energy_minus_control_work_residual_gamma'):
            assert row['verification'][key]<1e-8
    def differences(first,second):
        result={}
        for event in ('pulse_clear','receiver_plane_crossing','endpoint'):
            left,right=first[event],second[event]
            if left is None or right is None:
                result[event]=dict(same_event_presence=left is right)
                continue
            result[event]=dict(
                position_metres=float(np.linalg.norm(np.array(left['position_metres'])-
                                                     np.array(right['position_metres']))),
                lab_seconds=abs(left['lab_seconds']-right['lab_seconds']),
                proper_seconds=abs(left['proton_proper_seconds']-right['proton_proper_seconds']),
                gamma=abs(left['gamma']-right['gamma']))
        return result
    source_change=differences(traces['main'],traces['source'])
    step_change=differences(traces['source'],traces['fine'])
    for comparison in (source_change,step_change):
        for event in comparison.values():
            if 'position_metres' in event:
                assert event['position_metres']<1e-4 and event['gamma']<1e-8
                assert event['lab_seconds']<1e-12
    assert len({row['outcome'] for row in traces.values()})==1
    result=traces['fine'];crossing=result['receiver_plane_crossing']
    if crossing is not None:
        total=result['initial']['gamma']+e*result['initial']['receiver_potential_volts']/(m_p*c*c)
        total+=crossing['control_work_joules']/(m_p*c*c)
        gamma_predicted=total-e*crossing['receiver_potential_volts']/(m_p*c*c)
        assert abs(gamma_predicted-crossing['gamma'])<1e-8
    else:
        gamma_predicted=None
    matching=[]
    axis=frame[:,2];projector=np.eye(3)-np.outer(axis,axis)
    for row in launch_study['sensitivities']:
        derivative=np.array(row['levels'][0]['central_difference']['position_metres_per_input_unit'])
        projected=float(np.linalg.norm(projector@derivative))
        scale=9.9/projected if projected else None
        tested=max(level['symmetric_step_in_input_units'] for level in row['levels'])
        matching.append(dict(input=row['input'],input_unit=row['input_unit'],
            previous_receiver_plane_projected_position_sensitivity_metres_per_unit=projected,
            largest_previously_tested_amplitude=tested,
            linear_bore_equivalent_amplitude_within_previous_sample_range=(
                scale if scale is not None and scale<=tested else None),
            scope='Project the prior arrival sensitivity without the receiver field onto the fixed receiver plane. '
                'This is an alignment estimate for the 9.9m open bore, not the acceptance of the new electric lens; '
                'values outside the previously sampled launch range are not extrapolated'))
    kinetic_voltage=(design['reference_incoming_gamma']-1)*m_p*c*c/e
    return dict(
        model='One static net-neutral finite charged-surface receiver, present continuously throughout the original launch',
        design=design,receiver_frame='Columns of frame_columns are local transverse1,transverse2,axial unit vectors in lab coordinates; xlocal=frame.T@(xlab-center)',
        source_definition=dict(
            torus_surfaces='xlocal=((a+b*cos(theta))*cos(phi),(a+b*cos(theta))*sin(phi),b*sin(theta)), a=10m or20m,b=.1m',
            charge_measure='dQ=+/-Q*dtheta*dphi/(4*pi^2), sigma=+/-Q/[4*pi^2*b*(a+b*cos(theta))]. Both static surfaces have J=0; the signed total charge and dipole moment of the pair vanish',
            interpretation='Prescribed fixed-charge surfaces, not solved conducting electrodes. Surface potentials need not be constant; charge preparation, supports, dielectric response and field emission are not modeled',
            fields='Phi=(1/(4*pi*epsilon0))*integral dQ/|x-xs|; E=-gradPhi,Breceiver=0. The neutral pair uses rationalized inverse-distance differences to avoid cancellation. Original retarded current sources and static storage B are included from t=0',
            far_axis='Phi(z)=Q*(a_outer^2-a_inner^2)/(8*pi*epsilon0*|z|^3)+O(|z|^-5). Net-charge and dipole cancellation reduce but do not remove the upstream field',
            finite_surface_scope='The smooth charged sheets have finite extent and bounded density; their continuum Coulomb self energy is finite. This calculation does not evaluate that energy or a total apparatus budget; original storage coils remain ideal filaments'),
        source_audit=dict(surfaces=surfaces,
            maximum_charge_normalization_relative_error=max(normalization),
            characteristic_source_field_volts_per_metre=field_scale,
            point_checks=probe_rows,independent_axis_potential_checks=axis_checks,
            central_potential_hessian_volts_per_square_metre=hessian.tolist(),
            saddle_scope='The electrostatic center has positive transverse and negative axial potential curvature with zero trace. This diagnoses an electrostatic saddle; it is not a general theorem forbidding all combined electromagnetic traps',
            roundoff_scope='At the exact zero-field center, cancellation leaves a small absolute field residual. Its gradient check is normalized to the stated source field scale, with the absolute residual retained. High precision and numerical quadrature checks are not outward-rounded certificates'),
        energy_accounting=dict(
            reference_stopping_potential_rise_volts=kinetic_voltage,
            chosen_design_potential_rise_volts=design['design_potential_rise_volts'],
            design_gamma_is_achieved_constraint=False,
            exact_work_equations='Delta(gamma mc^2)=Wcontrol+Wreceiver; Wreceiver=-e[Phi(x)-Phi(x0)]. Thus gamma+e*Phi/(mc^2)-Wcontrol/(mc^2) is constant even during the pulse. After control clearance, gamma mc^2+ePhi alone is conserved',
            receiver_plane_gamma_from_actual_energy=gamma_predicted,
            interpretation='Receiver work first transfers kinetic energy into interaction potential energy and then returns it as the particle descends the static potential. No irreversible energy sink, detector absorption or charge switching is included'),
        trajectories=dict(main=traces['main'],source_refined=traces['source'],fine=traces['fine']),
        convergence=dict(source_quadrature_change=source_change,timestep_and_tolerance_change=step_change,
            criteria='Same physical source and launch in all three runs. Position1e-4m, time1e-12s and gamma1e-8 are numerical refinement checks, not target alignment tolerances or rigorous error bounds'),
        receiver_matching_context=dict(clear_central_bore_radius_metres=9.9,
            reference_arrival=next(row['control']['endpoint'] for row in torus['cases']
                                   if row['minor_radius_metres']==.1),
            prior_sensitivity_estimates=matching,
            new_acceptance_region_established=False,
            qualification='The prior cylindrical-checkpoint launch study supplied a nominal position/direction and local matching estimates. The receiver changes the global field and optics; only the nominal launch is recomputed here. A slowed crossing does not determine capture or acceptance for perturbed launches'),
        outcome=result['outcome'],capture_demonstrated=False,
        numerical_protocol='Exactly one fixed physical receiver, with main32/16 angular source orders, source64/32, and fine64/32 plus halved steps/tighter tolerances. Center-plane crossing and +/-200m exit/return events are measured in the fixed receiver frame, with a30us lab deadline. No charge, alignment or launch retuning',
        scope='The explicit finite static receiver slows the simulated proton without removing its energy permanently. It then transmits or reflects according to the computed trajectory. No old magnetic-only escape bound is applied in this new electric field. Capture needs a separately modeled time-dependent or dissipative stage, and a realizable charge/support/drive model remains missing. Radiation reaction, material response, detector capture, deceleration to rest, storage-duration preparation, reunion and human transport are outside this calculation',
        sources=[
            dict(url='https://www.damtp.cam.ac.uk/user/tong/em/el1.pdf',
                 support='Sections2.1.4 and2.2.3, Eq2.21: surface charges and Coulomb Green function; sections2.2.5 and2.3: electrostatic saddle/equilibrium and potential energy'),
            dict(url='https://arxiv.org/pdf/1909.07493',
                 support='GBAR electrostatic decelerator, section2 and Figure1: a pulsed drift tube prevents downstream reacceleration. This is a distinct low-energy apparatus, not evidence that the present multi-GV charged surfaces can be built')],
        derivation_attribution='Finite torus charge normalization, neutral-pair cancellation, fixed receiver design, full continuous Maxwell/Lorentz trajectory and matching estimates are our calculations. The cited papers supply established field theory and a distinct switching mechanism')

def _receiver_discharge_waveform(retarded_time, onset, duration, charge):
    """C4 prescribed endpoint charge and its first two time derivatives."""
    import numpy as np
    value=np.asarray(retarded_time,dtype=float)
    u=np.clip((value-onset)/duration,0.,1.)
    v=np.minimum(u,1-u)
    polynomial=v**5*(126+v*(-420+v*(540+v*(-315+70*v))))
    fraction=np.where(u<=.5,1-polynomial,polynomial)
    first=-630*charge/duration*u**4*(1-u)**4
    second=-2520*charge/duration**2*u**3*(1-u)**3*(1-2*u)
    return charge*fraction,first,second


def _switched_torus_receiver_kernel(design, onset, duration,
                                    azimuth_order=16,poloidal_order=8,radial_order=16,
                                    source_batch_size=None):
    """Exact retarded charge/current fields, quadrature-evaluated inside the open bore.

    Equal endpoint measures +/-q(t)dtheta dphi/(4pi^2) are joined by
    J=-qdot(t) Delta integral_0^1 delta(x-x_inner-s Delta)ds, averaged
    over theta,phi. The implied surface density is q/[4pi^2 b(a+b cos theta)].
    The radial current has an integrable edge enhancement, not a solved
    conductor response. Returns Phi,A,E,B in the fixed receiver frame.
    Optional design['poloidal_weighting']='cosine_squared' uses the bounded
    two-branch radial density induced by the measure 2*cos(theta)^2.
    Optional source batching changes only summation grouping, not source nodes.
    """
    import math
    import numpy as np
    from scipy.constants import c
    from numpy.polynomial.legendre import leggauss
    inner,outer,minor=(design['inner_major_radius_metres'],
        design['outer_major_radius_metres'],design['minor_radius_metres'])
    charge,eps=design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI']
    assert 0<minor<inner<outer and outer-inner>2*minor and duration>0
    k=1/(4*math.pi*eps);mu_over_4pi=k/c**2
    theta=2*math.pi*np.arange(poloidal_order)/poloidal_order
    phi=2*math.pi*np.arange(azimuth_order)/azimuth_order
    angular_weights=None
    weighting=design.get('poloidal_weighting')
    if weighting is not None:
        if weighting!='cosine_squared':raise ValueError('Unknown receiver angular measure')
        angular_weights=np.repeat(2*np.cos(theta)**2,azimuth_order)
        assert abs(float(np.mean(angular_weights))-1)<1e-14
    radial=np.repeat(inner+minor*np.cos(theta),azimuth_order)
    height=np.repeat(minor*np.sin(theta),azimuth_order)
    direction=np.column_stack((np.tile(np.cos(phi),poloidal_order),
        np.tile(np.sin(phi),poloidal_order),np.zeros(azimuth_order*poloidal_order)))
    first=radial[:,None]*direction;first[:,2]=height
    delta=(outer-inner)*direction;second=first+delta
    nodes,weights=leggauss(radial_order)
    static=_electrostatic_torus_pair_kernel(inner,outer,minor,charge,eps,
                                            azimuth_order,poloidal_order,weighting)
    def distances(position):
        rr=math.hypot(*position[:2]);zz=float(position[2])
        gap=math.hypot(rr-min(outer,max(inner,rr)),zz)-minor
        far=math.hypot(rr+outer,zz)+minor
        return gap,far
    def field(time,position):
        x=np.asarray(position,dtype=float)
        rr=math.hypot(*x[:2])
        if rr>=inner-minor:
            raise ValueError('Switched receiver radial quadrature is restricted to the open central bore')
        gap,far=distances(x)
        if time-gap/c<=onset:
            pot,vec=static(x)
            return float(pot),np.zeros(3),np.asarray(vec,dtype=float),np.zeros(3)
        if time-far/c>=onset+duration:
            return 0.,np.zeros(3),np.zeros(3),np.zeros(3)
        R1=x-first;R2=x-second
        d1=np.linalg.norm(R1,axis=1);d2=np.linalg.norm(R2,axis=1)
        q1,dq1,_=_receiver_discharge_waveform(time-d1/c,onset,duration,charge)
        q2,dq2,_=_receiver_discharge_waveform(time-d2/c,onset,duration,charge)
        if angular_weights is not None:
            q1=q1*angular_weights;dq1=dq1*angular_weights
            q2=q2*angular_weights;dq2=dq2*angular_weights
        potential=k*np.mean(q1/d1-q2/d2)
        electric=k*np.mean((q1/d1**3+dq1/(c*d1**2))[:,None]*R1-
                           (q2/d2**3+dq2/(c*d2**2))[:,None]*R2,axis=0)
        # For bore points, distance increases monotonically along each outward
        # radial segment. Split at both causal waveform fronts before Gauss integration.
        low=np.maximum(d1,c*(time-onset-duration))
        high=np.minimum(d2,c*(time-onset))
        active=high>low
        if not np.any(active):
            return float(potential),np.zeros(3),electric,np.zeros(3)
        along=direction@x
        impact=np.sum(x[:2]**2)-along**2+(x[2]-height)**2
        lo_source=along+np.sqrt(np.maximum(0.,low**2-impact))
        hi_source=along+np.sqrt(np.maximum(0.,high**2-impact))
        slo=np.where(low<=d1,0.,(lo_source-radial)/(outer-inner))
        shi=np.where(high>=d2,1.,(hi_source-radial)/(outer-inner))
        slo=np.clip(slo,0.,1.);shi=np.clip(shi,0.,1.)
        def current_fields(selection):
            midpoint=(slo[selection]+shi[selection])/2
            width=(shi[selection]-slo[selection])/2
            s=midpoint[:,None]+width[:,None]*nodes
            points=first[selection,None,:]+s[:,:,None]*delta[selection,None,:]
            R=x-points;distance=np.linalg.norm(R,axis=2)
            _,dq,ddq=_receiver_discharge_waveform(time-distance/c,onset,duration,charge)
            w=width[:,None]*weights[None,:]/len(first)
            if angular_weights is not None:w=w*angular_weights[selection,None]
            A=-mu_over_4pi*np.sum((w*dq/distance)[:,:,None]*delta[selection,None,:],axis=(0,1))
            E=mu_over_4pi*np.sum((w*ddq/distance)[:,:,None]*delta[selection,None,:],axis=(0,1))
            B=-mu_over_4pi*np.sum(
                (w*(dq/distance**3+ddq/(c*distance**2)))[:,:,None]*
                np.cross(delta[selection,None,:],R),axis=(0,1))
            return A,E,B
        if source_batch_size is None:
            A,current_E,magnetic=current_fields(active)
        else:
            if not isinstance(source_batch_size,int) or source_batch_size<=0:
                raise ValueError('Source batch size must be a positive integer or None')
            indices=np.flatnonzero(active)
            A=np.zeros(3);current_E=np.zeros(3);magnetic=np.zeros(3)
            for first_index in range(0,len(indices),source_batch_size):
                av,ev,bv=current_fields(indices[first_index:first_index+source_batch_size])
                A+=av;current_E+=ev;magnetic+=bv
        electric+=current_E
        return float(potential),A,electric,magnetic
    return field,distances


def _switched_receiver_definition(receiver,duration_seconds=100e-9):
    """Fix the one discharge from the saved actual fine center crossing, without retuning."""
    import copy
    import math
    initial=copy.deepcopy(receiver['trajectories']['fine']['receiver_plane_crossing'])
    assert initial is not None
    from scipy.constants import c
    design=copy.deepcopy(receiver['design'])
    local=initial['receiver_coordinates_metres']
    r=math.hypot(*local[:2])
    inner,outer,b=(design['inner_major_radius_metres'],
        design['outer_major_radius_metres'],design['minor_radius_metres'])
    gap=math.hypot(r-min(outer,max(inner,r)),local[2])-b
    onset=initial['lab_seconds']-gap/c
    if not math.isfinite(duration_seconds) or duration_seconds<=0:
        raise ValueError('Receiver discharge duration must be finite and positive')
    return dict(receiver=design,initial=initial,onset_seconds=onset,
        duration_seconds=duration_seconds,
        source_initial_nearest_distance_metres=gap,
        downstream_axial_metres=200.,upstream_axial_metres=-200.,
        maximum_elapsed_after_center_seconds=100e-6)


def proton_switched_receiver_trace(ring,coils,receiver,profile='main',duration_seconds=100e-9):
    """Continue the same incoming proton after its causally unchanged center crossing."""
    import math
    import numpy as np
    from scipy.constants import c,e,m_p
    from scipy.integrate import solve_ivp
    from scipy.special import ellipk,ellipe
    definition=_switched_receiver_definition(receiver,duration_seconds)
    design=definition['receiver'];start=definition['initial']
    center=np.array(design['center_metres']);frame=np.array(design['frame_columns'])
    onset,T=definition['onset_seconds'],definition['duration_seconds']
    t0=start['lab_seconds'];energy=m_p*c*c
    profiles={'main':(16,8,16,2e-11,2e-12,.025,1.),
              'source':(32,16,32,2e-11,2e-12,.025,1.),
              'fine':(32,16,32,3e-13,3e-14,.0125,.5),
              'fine_poloidal':(32,64,32,3e-13,3e-14,.0125,.5)}
    assert profile in profiles
    az,pol,rad,rtol,atol,step,coast_step=profiles[profile]
    source_batch_size=128 if profile=='fine_poloidal' else None
    field,distances=_switched_torus_receiver_kernel(design,onset,T,az,pol,rad,source_batch_size)
    source=coils['parameters']
    mu,current=source['vacuum_permeability_SI'],source['ampere_turns_per_loop']
    a,h=source['coil_radius_metres'],source['coil_z_metres'][1]
    pulse=design['pulse_parameters'];ring_radius=source['reference_orbit_radius_metres']
    control_a,control_h,control_b=(pulse['control_loop_radius_metres'],
        pulse['control_loop_half_separation_metres'],pulse['control_torus_minor_radius_metres'])
    pp=frame.T@np.array(start['momentum_over_mc'])
    initial=np.concatenate((np.array(start['receiver_coordinates_metres']),pp,[0.,0.]))
    gamma_initial=math.sqrt(1+float(np.dot(pp,pp)))
    factor=e*T/(m_p*c)
    def background(local):
        x=center+frame@local
        rr=math.hypot(*x[:2])
        _,br,bz=_finite_coil_loop_sum(rr,x[2],a,h,current,
                                     math.sqrt,ellipk,ellipe,math.pi,mu)
        return frame.T@np.array([br*x[0]/rr,br*x[1]/rr,bz])
    def gaps(local):
        x=center+frame@local;rr=math.hypot(x[0]-ring_radius,x[1])
        control=min(math.hypot(rr-control_a,x[2]-zz)-control_b for zz in(-control_h,control_h))
        rr=math.hypot(*x[:2])
        storage=min(math.hypot(rr-a,x[2]-zz) for zz in(-h,h))
        return np.array([distances(local)[0],control,storage])
    def rhs(cleared):
        def evaluate(s,state):
            position,pp=state[:3],state[3:6]
            gamma=math.sqrt(1+float(np.dot(pp,pp)))
            if cleared:E=np.zeros(3);B=np.zeros(3)
            else:_,_,E,B=field(t0+T*s,position)
            B=B+background(position)
            return np.concatenate((c*T*pp/gamma,
                factor*(E+c*np.cross(pp,B)/gamma),
                [1/gamma,factor*np.dot(pp,E)/gamma]))
        return evaluate
    def clear_event(s,state):
        return (t0+T*s-distances(state[:3])[1]/c-onset-T)/T
    clear_event.direction,clear_event.terminal=1,True
    def collision(s,state):return float(np.min(gaps(state[:3])))
    collision.direction,collision.terminal=-1,True
    def bore_exit(s,state):
        return design['inner_major_radius_metres']-design['minor_radius_metres']-math.hypot(*state[:2])
    bore_exit.direction,bore_exit.terminal=-1,True
    def downstream(s,state):return state[2]-definition['downstream_axial_metres']
    downstream.direction,downstream.terminal=1,True
    def upstream(s,state):return state[2]-definition['upstream_axial_metres']
    upstream.direction,upstream.terminal=-1,True
    def axial_turn(s,state):return state[5]
    axial_turn.direction,axial_turn.terminal=(0 if profile=='fine_poloidal' else -1),False
    limit=definition['maximum_elapsed_after_center_seconds']/T
    pulse=solve_ivp(rhs(False),(0,limit),initial,method='DOP853',dense_output=True,
        events=(clear_event,collision,bore_exit,downstream,upstream,axial_turn),
        rtol=rtol,atol=atol,max_step=step)
    cleared=bool(len(pulse.t_events[0]))
    coast=None
    if pulse.success and cleared:
        coast=solve_ivp(rhs(True),(pulse.t[-1],limit),pulse.y[:,-1],method='DOP853',
            dense_output=True,events=(downstream,upstream,collision,bore_exit),
            rtol=rtol,atol=atol,max_step=coast_step)
    last=coast if coast is not None else pulse
    if not last.success:outcome='integration failure'
    elif (coast is not None and len(coast.t_events[2])) or len(pulse.t_events[1]):
        outcome='source-support intersection'
    elif (coast is not None and len(coast.t_events[3])) or len(pulse.t_events[2]):
        outcome='left the explicitly modeled central-bore region'
    elif (coast is not None and len(coast.t_events[0])) or len(pulse.t_events[3]):
        outcome='reached downstream plane'
    elif (coast is not None and len(coast.t_events[1])) or len(pulse.t_events[4]):
        outcome='returned to upstream plane'
    else:outcome='no checkpoint within fixed hundred-microsecond continuation'
    finish=float(last.t[-1])
    def values(times):
        array=np.atleast_1d(times)
        out=np.empty((8,len(array)))
        before=array<=pulse.t[-1]
        if np.any(before):out[:,before]=pulse.sol(array[before])
        if np.any(~before):out[:,~before]=coast.sol(array[~before])
        return out
    nodes=np.unique(np.concatenate((pulse.t,last.t,np.linspace(0,finish,801))))
    for _ in range(12):
        sampled=values(nodes);gap_values=np.array([gaps(x) for x in sampled[:3].T])
        bounds=np.minimum(gap_values[:-1],gap_values[1:])-c*T*np.diff(nodes)[:,None]/2
        bad=np.any(bounds<=0,axis=1)
        if not np.any(bad) or outcome=='source-support intersection':break
        nodes=np.unique(np.concatenate((nodes,(nodes[:-1][bad]+nodes[1:][bad])/2)))
    sampled=values(nodes)
    gammas=np.sqrt(1+np.sum(sampled[3:6]**2,axis=0))
    work_error=gammas-gamma_initial-sampled[7]
    assert np.all(np.isfinite(sampled)) and np.all(np.diff(sampled[6])>=0)
    if outcome in ('reached downstream plane','returned to upstream plane'):
        assert np.min(bounds)>0 and np.max(np.abs(work_error))<1e-8
    def state_record(s,state):
        local,pp=state[:3],state[3:6]
        gamma=math.sqrt(1+float(np.dot(pp,pp)))
        time=t0+T*s
        potential,A,E,B=field(time,local)
        return dict(lab_seconds=float(time),
            proton_proper_seconds=float(start['proton_proper_seconds']+T*state[6]),
            position_metres=(center+frame@local).tolist(),
            momentum_over_mc=(frame@pp).tolist(),receiver_coordinates_metres=local.tolist(),
            receiver_transverse_offset_metres=math.hypot(*local[:2]),
            receiver_axial_momentum_over_mc=float(pp[2]),gamma=gamma,
            beta=math.sqrt(1-1/gamma**2),kinetic_energy_eV=(gamma-1)*energy/e,
            receiver_potential_volts=float(potential),
            receiver_electric_field_volts_per_metre=float(np.linalg.norm(E)),
            receiver_magnetic_field_tesla=float(np.linalg.norm(B)),
            added_receiver_work_joules=float(state[7]*energy),
            total_receiver_work_joules=float(start['receiver_work_joules']+state[7]*energy),
            control_work_joules=start['control_work_joules'])
    clear_state=state_record(float(pulse.t[-1]),pulse.y[:,-1]) if cleared else None
    endpoint=state_record(finish,last.y[:,-1])
    min_index=int(np.argmin(gammas))
    minimum=state_record(float(nodes[min_index]),sampled[:,min_index])
    turns=[state_record(float(s),state) for s,state in zip(pulse.t_events[5],pulse.y_events[5])]
    if profile=='fine_poloidal':
        for row,s,state in zip(turns,pulse.t_events[5],pulse.y_events[5]):
            row['axial_momentum_direction_change']=('increasing' if rhs(False)(s,state)[5]>0 else 'decreasing')
    if turns:
        minimum=min([minimum,*turns],key=lambda row:row['gamma'])
    events=[0.,finish,float(nodes[min_index]),*pulse.t_events[5]]
    if cleared:events.append(float(pulse.t[-1]))
    plot_times=np.unique(np.concatenate((np.linspace(0,finish,401),
        np.linspace(0,min(finish,pulse.t[-1]),201),events)))
    plot_states=values(plot_times)
    plot_gamma=np.sqrt(1+np.sum(plot_states[3:6]**2,axis=0))
    postclear=nodes>=pulse.t[-1] if cleared else np.zeros(len(nodes),dtype=bool)
    gamma_coast_error=(float(np.max(np.abs(gammas[postclear]-clear_state['gamma'])))
                       if cleared else None)
    return dict(profile=profile,definition=definition,outcome=outcome,
        initial=state_record(0.,initial),receiver_pulse_clear=clear_state,
        minimum_gamma_sample_or_turn=minimum,axial_turns=turns,endpoint=endpoint,
        numerical_profile=dict(azimuth_order=az,poloidal_order=pol,radial_order=rad,
            rtol=rtol,atol=atol,maximum_switch_step_seconds=step*T,
            maximum_coast_step_seconds=coast_step*T,
            switch_rhs_evaluations=pulse.nfev,coast_rhs_evaluations=(coast.nfev if coast else 0),
            **({'source_batch_size':source_batch_size} if source_batch_size is not None else {})),
        verification=dict(maximum_added_work_residual_gamma=float(np.max(np.abs(work_error))),
            postclear_gamma_variation=gamma_coast_error,
            continuous_source_gap_lower_bounds_metres=np.min(bounds,axis=0).tolist(),
            source_gap_order=['full receiver radial-current capsule','original control tori','original storage filaments'],
            gap_scope='Signed 1-Lipschitz source distances and c*Delta t/2 interval reduction. Conditional on numerical endpoint positions, not interval-certified integration',
            clear_event_residual_seconds=(float(clear_event(pulse.t[-1],pulse.y[:,-1])*T) if cleared else None),
            maximum_transverse_offset_metres=float(np.max(np.hypot(sampled[0],sampled[1]))),
            solver_success=bool(pulse.success and (coast is None or coast.success))),
        plot=dict(lab_microseconds=((t0+T*plot_times)*1e6).tolist(),
            receiver_axial_metres=plot_states[2].tolist(),
            receiver_transverse_metres=np.hypot(plot_states[0],plot_states[1]).tolist(),
            gamma=plot_gamma.tolist(),kinetic_energy_MeV=((plot_gamma-1)*energy/e/1e6).tolist(),
            proton_proper_microseconds=((start['proton_proper_seconds']+T*plot_states[6])*1e6).tolist(),
            receiver_work_joules=(start['receiver_work_joules']+energy*plot_states[7]).tolist(),
            x_metres=(center[0]+(frame@plot_states[:3])[0]).tolist(),
            y_metres=(center[1]+(frame@plot_states[:3])[1]).tolist(),
            z_metres=(center[2]+(frame@plot_states[:3])[2]).tolist()))


def proton_switched_receiver_checks(ring,coils,receiver,trajectories=None):
    """One fixed causal discharge: source identities, three numerical resolutions, and outcome."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c,e,m_p
    from numpy.polynomial.legendre import leggauss
    definition=_switched_receiver_definition(receiver)
    design=definition['receiver'];initial=definition['initial']
    onset,T=definition['onset_seconds'],definition['duration_seconds']
    Q,eps=design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI']
    inner,outer,b=(design['inner_major_radius_metres'],
        design['outer_major_radius_metres'],design['minor_radius_metres'])
    k=1/(4*math.pi*eps);energy=m_p*c*c
    field,distances=_switched_torus_receiver_kernel(design,onset,T,32,16,32)
    coarse,_=_switched_torus_receiver_kernel(design,onset,T,16,8,16)
    high,_=_switched_torus_receiver_kernel(design,onset,T,64,32,48)
    q=lambda u:_receiver_discharge_waveform(onset+u*T,onset,T,Q)
    transferred=quad(lambda u:-float(q(u)[1])*T,0,1,epsabs=1e-11,epsrel=1e-12)[0]
    assert abs(transferred/Q-1)<1e-12
    symmetry=max(abs(float(q(u)[0]+q(1-u)[0])-Q) for u in np.linspace(0,1,33))
    assert symmetry/Q<1e-12
    peak_current=630*Q/T/256
    assert abs(-float(q(.5)[1])/peak_current-1)<1e-12
    assert all(float(q(u)[j])==0. for u in (-1.,2.) for j in (1,2))
    assert float(q(-1.)[0])==Q and float(q(2.)[0])==0

    # Independent weak continuity check against a polynomial test function.
    theta=2*math.pi*np.arange(16)/16;phi=2*math.pi*np.arange(32)/32
    er=np.column_stack((np.tile(np.cos(phi),16),np.tile(np.sin(phi),16),np.zeros(512)))
    first=np.repeat(inner+b*np.cos(theta),32)[:,None]*er
    first[:,2]=np.repeat(b*np.sin(theta),32)
    delta=(outer-inner)*er;second=first+delta
    sn,sw=leggauss(8);s=(sn+1)/2;sw=sw/2
    points=first[:,None,:]+s[None,:,None]*delta[:,None,:]
    def psi(p):
        x,y,z=np.moveaxis(p,-1,0)
        return x*x+3*y*y+2*z*z+.01*x**4+.02*y*y*z*z
    def grad_psi(p):
        x,y,z=np.moveaxis(p,-1,0)
        return np.stack((2*x+.04*x**3,6*y+.04*y*z*z,4*z+.04*y*y*z),axis=-1)
    endpoint_difference=float(np.mean(psi(second)-psi(first)))
    line_integral=float(np.mean(np.sum(sw[None,:]*np.sum(
        grad_psi(points)*delta[:,None,:],axis=2),axis=1)))
    continuity_error=abs(line_integral/endpoint_difference-1)
    assert continuity_error<1e-13

    field_scale=k*Q/inner**2;derivative_scale=field_scale/inner
    def derivative(fn,step):
        return (fn(-2*step)-8*fn(-step)+8*fn(step)-fn(2*step))/(12*step)
    def curl(jac):
        return np.array([jac[2,1]-jac[1,2],jac[0,2]-jac[2,0],jac[1,0]-jac[0,1]])
    probes=[]
    for point,delay in (([.1,.2,3.],45e-9),([.2,-.1,5.],90e-9),
                        ([.1,.2,11.],140e-9)):
        x=np.array(point);time=onset+delay
        potential,A,E,B=field(time,x)
        low=coarse(time,x);upper=high(time,x)
        derivative_rows=[]
        for step in (.002,.001):
            dt=step/c
            gradients=[[],[],[],[]]
            for j in range(3):
                basis=np.eye(3)[j]
                for index in range(4):
                    gradients[index].append(derivative(
                        lambda h:field(time,x+h*basis)[index],step))
            grad_phi=np.array(gradients[0])
            ja,je,jb=(np.array(gradients[j]).T for j in (1,2,3))
            dphi=derivative(lambda h:field(time+h,x)[0],dt)
            da=derivative(lambda h:field(time+h,x)[1],dt)
            de=derivative(lambda h:field(time+h,x)[2],dt)
            db=derivative(lambda h:field(time+h,x)[3],dt)
            row=dict(spatial_step_metres=step,time_step_seconds=dt,
                E_from_potential_relative=float(np.linalg.norm(E+grad_phi+da)/field_scale),
                B_from_potential_relative=float(c*np.linalg.norm(B-curl(ja))/field_scale),
                divergence_E_relative=float(abs(np.trace(je))/derivative_scale),
                divergence_B_relative=float(c*abs(np.trace(jb))/derivative_scale),
                Faraday_relative=float(np.linalg.norm(curl(je)+db)/derivative_scale),
                Ampere_relative=float(c*np.linalg.norm(curl(jb)-de/c**2)/derivative_scale),
                Lorenz_gauge_relative=float(c*abs(np.trace(ja)+dphi/c**2)/field_scale))
            assert max(value for key,value in row.items() if key.endswith('_relative'))<1e-8
            derivative_rows.append(row)
        changes=[]
        for other in (low,upper):
            change=dict(potential_volts=float(abs(other[0]-potential)),
                A_tesla_metres=float(np.linalg.norm(other[1]-A)),
                electric_volts_per_metre=float(np.linalg.norm(other[2]-E)),
                magnetic_tesla=float(np.linalg.norm(other[3]-B)))
            assert max(change['electric_volts_per_metre'],
                       c*change['magnetic_tesla'])/field_scale<1e-9
            changes.append(change)
        probes.append(dict(receiver_coordinates_metres=point,seconds_after_source_onset=delay,
            potential_volts=float(potential),electric_vector_volts_per_metre=E.tolist(),
            magnetic_vector_tesla=B.tolist(),potential_and_vacuum_Maxwell_checks=derivative_rows,
            source_changes=dict(main_to_source=changes[0],source_to_higher_order=changes[1])))
    axis_checks=[]
    for z,delay in ((3.,45e-9),(5.,90e-9),(11.,140e-9)):
        time=onset+delay
        def integrand(theta):
            zz=z-b*math.sin(theta)
            result=0.
            for a,sign in ((inner,1),(outer,-1)):
                distance=math.hypot(a+b*math.cos(theta),zz)
                qq,dq,_=_receiver_discharge_waveform(time-distance/c,onset,T,Q)
                result+=sign*(qq/distance**3+dq/(c*distance**2))*zz
            return k*result/(2*math.pi)
        exact=quad(integrand,0,2*math.pi,epsabs=1e-4,epsrel=1e-12)[0]
        values=field(time,np.array([0.,0.,z]))
        error=abs(float(values[2][2])-exact)/field_scale
        assert error<1e-11
        assert max(np.linalg.norm(values[1]),np.linalg.norm(values[3]))<1e-10
        axis_checks.append(dict(axial_metres=z,seconds_after_source_onset=delay,
            independent_axis_E_volts_per_metre=exact,
            relative_difference_over_source_scale=error))
    gap,far=distances(np.array(initial['receiver_coordinates_metres']))
    assert abs(gap-definition['source_initial_nearest_distance_metres'])<1e-14
    onset_residual=initial['lab_seconds']-gap/c-onset
    assert abs(onset_residual)<1e-20
    p=np.array(initial['position_metres'])
    radius=coils['parameters']['reference_orbit_radius_metres']
    pulse=design['pulse_parameters']
    local_r=math.hypot(p[0]-radius,p[1])
    control_far=max(math.hypot(local_r+pulse['control_loop_radius_metres'],p[2]-z)
        for z in (-pulse['control_loop_half_separation_metres'],
                   pulse['control_loop_half_separation_metres']))+pulse['control_torus_minor_radius_metres']
    old_control_margin=initial['lab_seconds']-control_far/c-pulse['source_pulse_duration_seconds']
    assert old_control_margin>0

    if trajectories is None:
        trajectories=[proton_switched_receiver_trace(ring,coils,receiver,profile)
                      for profile in ('main','source','fine')]
    if isinstance(trajectories,dict):
        trajectories=[trajectories[key] for key in ('main','source_refined','fine')]
    traces={row['profile']:copy.deepcopy(row) for row in trajectories}
    assert len(trajectories)==3 and set(traces)=={'main','source','fine'}
    for name,row in traces.items():
        assert row['definition']==definition
        assert row['verification']['solver_success']
        assert row['verification']['maximum_added_work_residual_gamma']<1e-8
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
        if row['receiver_pulse_clear'] is not None:
            assert row['verification']['postclear_gamma_variation']<1e-10
            assert abs(row['verification']['clear_event_residual_seconds'])<1e-15
        if name!='fine':row.pop('plot',None)
    def changes(left,right):
        output={}
        for event in ('receiver_pulse_clear','endpoint'):
            a,b=left[event],right[event]
            if a is None or b is None:
                output[event]=dict(same_event_presence=a is b);continue
            output[event]=dict(position_metres=float(np.linalg.norm(
                np.array(a['position_metres'])-b['position_metres'])),
                lab_seconds=abs(a['lab_seconds']-b['lab_seconds']),
                proper_seconds=abs(a['proton_proper_seconds']-b['proton_proper_seconds']),
                gamma=abs(a['gamma']-b['gamma']))
        return output
    source_change=changes(traces['main'],traces['source'])
    step_change=changes(traces['source'],traces['fine'])
    for comparison in (source_change,step_change):
        for row in comparison.values():
            if 'position_metres' in row:
                assert row['position_metres']<1e-4 and row['lab_seconds']<1e-12
                assert row['gamma']<1e-8
    assert len({row['outcome'] for row in traces.values()})==1
    final=traces['fine'];end=final['endpoint'];clear=final['receiver_pulse_clear']
    if clear:
        at_clear=np.array(clear['receiver_coordinates_metres'])
        # Monotonic source retarded times keep all receiver fields zero at this
        # particle and every subsequent timelike continuation, not everywhere in space.
        after=field(clear['lab_seconds']+1e-12,at_clear)
        assert after[0]==0 and all(np.max(np.abs(v))==0 for v in after[1:])
    initial_launch=receiver['trajectories']['fine']['initial']
    full_work=(end['gamma']-initial_launch['gamma']-
        (end['total_receiver_work_joules']+end['control_work_joules'])/energy)
    assert abs(full_work)<1e-8
    multipole=Q*(outer**2-inner**2)/2
    enclosing=20.2
    exterior_floor=3*k*multipole**2/(10*enclosing**5)
    source_times=np.linspace(onset-.1*T,onset+1.1*T,241)
    charges,currents,_=_receiver_discharge_waveform(source_times,onset,T,Q)
    charge_transfer_peak_time=onset+T/2
    static_end=receiver['trajectories']['fine']['endpoint']
    # A deliberately loose Lienard envelope on this nonradiating trajectory.
    # Whole-interval source-distance bounds, rather than sampled field maxima,
    # control both the retarded fields and an independent gamma upper bound.
    radiation=None
    if clear is not None:
        receiver_gap,_,storage_gap=final['verification']['continuous_source_gap_lower_bounds_metres']
        peak_second=2520*Q/T**2*(3/14)**3/math.sqrt(7)
        second_u=(1-1/math.sqrt(7))/2
        assert abs(abs(float(q(second_u)[2]))/peak_second-1)<1e-11
        delta=outer-inner
        Ebound=2*k*(Q/receiver_gap**2+peak_current/(c*receiver_gap))
        Ebound+=k/c**2*peak_second*delta/receiver_gap
        Bbound=k/c**2*delta*(peak_current/receiver_gap**2+peak_second/(c*receiver_gap))
        storage=coils['parameters']
        Bstorage=storage['vacuum_permeability_SI']*abs(storage['ampere_turns_per_loop'])*storage['coil_radius_metres']/storage_gap**2
        switch_interval=clear['lab_seconds']-initial['lab_seconds']
        coast_interval=end['lab_seconds']-clear['lab_seconds']
        gamma_bound=initial['gamma']+e*Ebound*switch_interval/(m_p*c)
        prefactor=e**4/(6*math.pi*eps*m_p**2*c**3)
        radiation_bound=prefactor*gamma_bound**2*(Ebound+c*(Bbound+Bstorage))**2*switch_interval
        radiation_bound+=prefactor*clear['gamma']**2*(c*Bstorage)**2*coast_interval
        assert radiation_bound>0 and radiation_bound/(end['kinetic_energy_eV']*e)<1e-7
        radiation=dict(
            peak_absolute_q_second_derivative_C_per_s2=peak_second,
            receiver_E_envelope_volts_per_metre=Ebound,
            receiver_B_envelope_tesla=Bbound,storage_B_envelope_tesla=Bstorage,
            gamma_envelope_during_switch=gamma_bound,
            integrated_proton_radiation_upper_bound_joules=radiation_bound,
            upper_bound_over_retained_kinetic_energy=radiation_bound/(end['kinetic_energy_eV']*e),
            field_bounds='d is the continuous receiver-source gap and delta=a_outer-a_inner. |E|<=2ke(Q/d^2+Imax/(cd))+ke*qddmax*delta/(c^2*d); |Breceiver|<=ke*delta/c^2*(Imax/d^2+qddmax/(cd)); |Bstorage|<=mu0*Istorage*a_storage/dstorage^2 for the two loops',
            radiation_equation='P=e^4*gamma^2/(6*pi*epsilon0*mp^2*c^3)*[(E+v cross B)^2-(v dot E)^2/c^2] <=prefactor*gamma^2*(Ebound+c*Bbound_total)^2. dgamma/dt<=e*Ebound/(mp*c) supplies the independent gamma envelope',
            scope='Conditional on the numerical continuous source-gap bounds, this bounds Lienard radiation evaluated on the prescribed nonradiating path from the reused receiver-center state to the recorded upstream endpoint, excluding the original launch segment. It is not a self-force trajectory-error bound or control-source radiation. The loose bound exceeds the momentary axial-turn kinetic energy, so that tiny energy is not a physically precise stopping or capture prediction')
    result=dict(
        model='One fixed causal discharge of the finite charged receiver, with explicit radial transfer current and retarded Maxwell fields',
        definition=definition,
        source_definition=dict(
            endpoint_charges='rho=q(t)<delta(x-x1)-delta(x-x2)>, with <...>=integral dtheta dphi/(4pi^2) and x_i=((a_i+b*cos(theta))*cos(phi),(a_i+b*cos(theta))*sin(phi),b*sin(theta))',
            endpoint_surface_density='sigma_i(theta,t)=+/-q(t)/[4*pi^2*b*(a_i+b*cos(theta))]; these are prescribed charge sheets, not equipotential conducting electrodes',
            transfer_current='J=-qdot(t)<Delta*integral_0^1 delta(x-x1-s*Delta) ds>, Delta=(a_outer-a_inner)*e_r. Its distributional divergence cancels rhodot exactly',
            current_support='The annular capsule obtained by rotating dist((r,z),[a_inner,a_outer] x {0})<=b about the receiver axis, including both charged boundary tori',
            waveform='q=Q[1-126u^5+420u^6-540u^7+315u^8-70u^9] for u=(t-ton)/T in(0,1), Q before and0 after. Stable symmetric evaluation is used nearu=1; q is C4',
            transfer_distance_metres=outer-inner,
            peak_outward_transfer_current_amperes=peak_current,
            peak_transfer_current_lab_seconds=charge_transfer_peak_time,
            current_density_limitation='In the interior annulus Jr=-qdot/[2*pi^2*r*sqrt(b^2-z^2)]. Its edge divergence is integrable, but J^2 has a logarithmically divergent integral for a finite positive resistivity. This ideal source needs a regularized current distribution before any resistive material model or engineering claim',
            timing_limitation='The waveform is prescribed at all distributed source points by synchronized preprogrammed control. No instantaneous command propagation is assumed or demonstrated; feed, carrier motion, material response and preparation remain unmodeled',
            source_transfer_light_time_seconds=(outer-inner)/c,
            field_equations='Phi=ke integral[rho]/d; A=ke/c^2 integral[J]/d; E=ke integral([rho]R/d^3+[rhodot]R/(c*d^2)-[Jdot]/(c^2*d)); B=ke/c^2 integral([J]/d^3+[Jdot]/(c*d^2)) cross R. Brackets denote actual retarded time',
            magnetic_prefactor='Receiver mu=1/(epsilon0*c^2) enforces exact Maxwell normalization; the unchanged independently prescribed static storage B retains its earlier saved mu0',
            kernel_domain='The numerical radial-front splitting is valid in the open central cylinder r<a_inner-b. Every reported trajectory remains in that bore; no exterior kernel or arbitrary field cutoff is used'),
        source_audit=dict(transferred_charge_coulombs=transferred,
            total_charge_conserved=True,waveform_symmetry_relative_error=symmetry/Q,
            weak_continuity_polynomial_relative_error=continuity_error,
            vacuum_field_scale_volts_per_metre=field_scale,
            potential_and_Maxwell_probes=probes,independent_axis_checks=axis_checks,
            source_refinement_orders=dict(main=[16,8,16],source=[32,16,32],point_only_higher=[64,32,48]),
            numerical_scope='Finite-difference field identities, independent axial quadrature, and angular/radial refinement are numerical checks of the exact retarded construction, not outward-rounded interval certificates'),
        causal_reuse=dict(
            exact_saved_state_key='electrostatic_receiver.trajectories.fine.receiver_plane_crossing',
            onset_condition_residual_seconds=onset_residual,
            original_control_all_source_clearance_margin_seconds=old_control_margin,
            proof='For each fixed source point, d(t-|x(t)-xs|/c)/dt>=1-|v|/c>0. At the reused actual center state all retarded receiver source times<=ton by the full-support minimum distance. The complete earlier gamma10 launch is therefore unchanged. All position, momentum, absolute proper/lab times and work are carried forward; no fresh low-energy tracer or unequal-field splice is introduced',
            inherited_accuracy='The incoming state has the prior full-route numerical error and matching limitations. The present refinement tests the continuation conditional on that fixed state and timing',
            after_receiver_clearance='The outgoing electromagnetic disturbance remains elsewhere, but receiver fields vanish at the proton and along its future timelike continuation. Original static storage B alone remains there'),
        trajectories=dict(main=traces['main'],source_refined=traces['source'],fine=final),
        convergence=dict(source_quadrature_change=source_change,timestep_and_tolerance_change=step_change,
            criteria='Same fixed source, onset and incoming state. Position1e-4m, eventtime1e-12s and gamma1e-8 are numerical comparisons, not receiver acceptance or physical tolerances'),
        energy_accounting=dict(
            equation='gamma(t)-gamma(center)=W_receiver_after_center/(mp*c^2); include the inherited receiver and magnetic-control electric work for the complete launch. Gamma+ePhi/(mp*c^2) is not conserved during switching',
            complete_launch_work_residual_gamma=full_work,
            launch_kinetic_energy_eV=initial_launch['kinetic_energy_eV'],
            final_kinetic_energy_eV=end['kinetic_energy_eV'],
            launch_to_final_kinetic_energy_removed_eV=initial_launch['kinetic_energy_eV']-end['kinetic_energy_eV'],
            static_receiver_downstream_kinetic_energy_eV=static_end['kinetic_energy_eV'],
            static_receiver_downstream_gamma=static_end['gamma'],
            comparison_scope='The static run transmits to+200m; the switched run outcome is measured independently and may return to-200m. These endpoints are different events, not a same-target efficiency comparison. Single-particle electric work does not determine total source dissipation or radiation'),
        source_energy_lower_bound=dict(
            quadrupole_C_m2=multipole,enclosing_radius_metres=enclosing,
            initial_field_energy_outside_sphere_lower_bound_joules=exterior_floor,
            equation='Q2=Q(a_outer^2-a_inner^2)/2; Uoutside(R)>=3*ke*Q2^2/(10*R^5) forR>a_outer+b',
            scope='Positive orthogonal multipole contributions give an exact lower bound for the initial prescribed receiver field outside this sphere. This is not a total-energy estimate, emitted-radiation calculation, driver-energy demand or universal receiver bound'),
        proton_radiation_diagnostic=radiation,
        outcome=final['outcome'],
        downstream_transport_achieved=final['outcome']=='reached downstream plane',
        retained_kinetic_energy_reduction_demonstrated=bool(clear is not None and
            end['gamma']<initial_launch['gamma'] and final['verification']['postclear_gamma_variation']<1e-10),
        capture_demonstrated=False,
        plot=dict(source=dict(lab_microseconds=(source_times*1e6).tolist(),
            remaining_inner_charge_coulombs=charges.tolist(),
            outward_transfer_current_amperes=(-currents).tolist()),
            static_reference=copy.deepcopy(receiver['trajectories']['fine']['plot'])),
        scope='This is one untuned prescribed-source experiment. Retained energy removal, transmission and permanent capture are distinct outcomes. The full induction/radiation fields can reverse the proton and restore part of its kinetic energy. No dissipative trap, radiation reaction, finite-current-density conductor, support or feed solution, timing acceptance, capture, final rest, reunion or human transport is established',
        sources=[
            dict(url='https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node145.html',
                 support='Fitzpatrick, Retarded Potentials, Eqs1746-1747: causal scalar and vector potentials for conserved charge/current distributions'),
            dict(url='https://farside.ph.utexas.edu/teaching/jk1/lectures/node32.html',
                 support='Multipole expansion used in the independently audited positive exterior-energy lower bound'),
            dict(url='https://farside.ph.utexas.edu/teaching/em/lectures/node56.html',
                 support='Electrostatic field energy; source-field energy is distinct from one test-particle work'),
            dict(url='https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node156.html',
                 support='Fitzpatrick Eq1944: relativistic Larmor/Lienard radiation. Substituting the Lorentz force and applying the source-gap bounds gives the conditional diagnostic here'),
            dict(url='https://arxiv.org/pdf/1909.07493',
                 support='GBAR section2/Figure1 supplies a distinct experimental example of timed electrostatic deceleration; it does not validate the present enormous source parameters')],
        derivation_attribution='The finite torus transfer current, waveform, causal reuse, retarded continuation and source lower bound are our calculations or independent team derivations from established Maxwell theory, not apparatus results from the cited experiments')
    return result


def _regularized_receiver_definition(receiver,coils):
    """Change only the normalized poloidal charge measure; keep the previous Q and geometry."""
    import copy
    import math
    import numpy as np
    design=copy.deepcopy(receiver['design'])
    design['poloidal_weighting']='cosine_squared'
    design['original_uniform_design_gamma']=design.pop('design_gamma')
    field=_electrostatic_torus_pair_kernel(10.,20.,.1,
        design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI'],64,32,
        design['poloidal_weighting'])
    radius=coils['parameters']['reference_orbit_radius_metres']
    phase=design['pulse_parameters']['starting_azimuth_radians']
    launch=radius*np.array([math.cos(phase),math.sin(phase),0.])
    frame=np.array(design['frame_columns']);center=np.array(design['center_metres'])
    design['central_potential_volts']=float(field(np.zeros(3))[0])
    design['initial_potential_volts']=float(field(frame.T@(launch-center))[0])
    design['design_potential_rise_volts']=design['central_potential_volts']-design['initial_potential_volts']
    design['charge_choice']='The previous uniform-measure Q is held fixed. Reweighting changes the static field and achieved gamma; no charge, launch or source retuning is performed'
    return design


def _regularized_receiver_source_checks(definition,probe_delays_seconds=None,
                                       spatial_steps_metres=(.002,.001),
                                       quadrature_orders=((16,8,16),(32,16,32),(64,32,48)),
                                       source_batch_size=None):
    """Weighted source normalization, branch-resolved finite current norm, and Maxwell checks."""
    import math
    import numpy as np
    from scipy.constants import c
    from numpy.polynomial.legendre import leggauss
    design=definition['receiver'];assert design['poloidal_weighting']=='cosine_squared'
    ai,ao,b=(design['inner_major_radius_metres'],design['outer_major_radius_metres'],
             design['minor_radius_metres'])
    Q,eps=design['inner_total_charge_coulombs'],design['vacuum_permittivity_SI']
    onset,T=definition['onset_seconds'],definition['duration_seconds']
    delays=((45e-9,90e-9,140e-9) if probe_delays_seconds is None else tuple(probe_delays_seconds))
    assert len(delays)==3
    assert len(spatial_steps_metres)==2 and all(h>0 for h in spatial_steps_metres)
    assert len(quadrature_orders)==3 and all(len(orders)==3 for orders in quadrature_orders)
    k=1/(4*math.pi*eps);scale=k*Q/ai**2
    def logarithms(D):
        return (math.log1p(2*D/(ai-D))+4*math.log((ao-D)/(ai+D))+
                math.log1p(2*D/(ao-D)))
    C=quad(lambda u:(1-u*u)*logarithms(b*math.sqrt(max(0.,1-u*u))),
           -1,1,epsabs=1e-13,epsrel=1e-13)[0]/(2*math.pi**3*b)
    # A separate two-dimensional quadrature squares the branch sum before integration.
    zn,zw=leggauss(128);rn,rw=leggauss(32)
    direct=0.
    for z,wz in zip(b*zn,b*zw):
        D=math.sqrt(b*b-z*z)
        for lo,hi,count in ((ai-D,ai+D,1),(ai+D,ao-D,2),(ao-D,ao+D,1)):
            r=(lo+hi)/2+(hi-lo)*rn/2
            density=D*count/(2*math.pi**2*r*b*b)
            direct+=wz*(hi-lo)/2*float(np.sum(rw*2*math.pi*r*density**2))
    assert C>0 and abs(direct/C-1)<1e-9
    time_norm=4410/2431*Q*Q/T
    measured_time_norm=quad(lambda u:float(_receiver_discharge_waveform(
        onset+u*T,onset,T,Q)[1])**2*T,0,1,epsabs=.01,epsrel=1e-12)[0]
    assert abs(measured_time_norm/time_norm-1)<1e-11
    peak_current=630*Q/(256*T)
    max_density=max(peak_current/(math.pi**2*b*(ai+b)),
                    peak_current/(2*math.pi**2*b*(ai-b)))
    flux_checks=[]
    for z in (-.09,-.05,0.,.05,.09):
        D=math.sqrt(b*b-z*z)
        at_radius=(ai+ao)/2
        density=D/(math.pi**2*at_radius*b*b)
        flux_checks.append(dict(z_metres=z,central_annulus_J_over_current_per_m2=density))
    flux=quad(lambda z:2*math.pi*((ai+ao)/2)*math.sqrt(max(0.,b*b-z*z))/
              (math.pi**2*((ai+ao)/2)*b*b),-b,b,epsabs=1e-12)[0]
    assert abs(flux-1)<1e-12
    theta=2*math.pi*np.arange(32)/32
    weights=2*np.cos(theta)**2
    assert abs(float(np.mean(weights))-1)<1e-14 and np.min(weights)>=0
    phi=2*math.pi*np.arange(64)/64
    directions=np.column_stack((np.tile(np.cos(phi),32),np.tile(np.sin(phi),32),np.zeros(2048)))
    first=np.repeat(ai+b*np.cos(theta),64)[:,None]*directions
    first[:,2]=np.repeat(b*np.sin(theta),64)
    delta=(ao-ai)*directions;second=first+delta
    angular=np.repeat(weights,64)
    sn,sw=leggauss(8);s=(sn+1)/2;sw=sw/2
    points=first[:,None,:]+s[None,:,None]*delta[:,None,:]
    def psi(p):
        x,y,z=np.moveaxis(p,-1,0)
        return x*x+3*y*y+2*z*z+.01*x**4+.02*y*y*z*z
    def grad(p):
        x,y,z=np.moveaxis(p,-1,0)
        return np.stack((2*x+.04*x**3,6*y+.04*y*z*z,4*z+.04*y*y*z),axis=-1)
    endpoints=float(np.mean(angular*(psi(second)-psi(first))))
    line=float(np.mean(angular*np.sum(sw[None,:]*np.sum(grad(points)*delta[:,None,:],axis=2),axis=1)))
    assert abs(line/endpoints-1)<1e-13
    # The paired l=2 multipole is unaffected by the common minor-radius weighting.
    q2=float(np.mean(angular*(Q*(3*first[:,2]**2-np.sum(first*first,axis=1))/2-
                                       Q*(3*second[:,2]**2-np.sum(second*second,axis=1))/2)))
    assert abs(q2/(Q*(ao*ao-ai*ai)/2)-1)<1e-14
    static=_electrostatic_torus_pair_kernel(ai,ao,b,Q,eps,64,32,'cosine_squared')
    static_rows=[]
    for z in (0.,5.,100.,1000.):
        def integrand(t):
            zz=z-b*math.sin(t)
            di=math.hypot(ai+b*math.cos(t),zz);do=math.hypot(ao+b*math.cos(t),zz)
            return (k*Q*2*math.cos(t)**2*((ao-ai)*(ao+ai+2*b*math.cos(t)))/
                    (di*do*(di+do))/(2*math.pi))
        expected=quad(integrand,0,2*math.pi,epsabs=1e-4,epsrel=1e-12)[0]
        value=float(static(np.array([0.,0.,z]))[0])
        assert abs(value/expected-1)<1e-12
        static_rows.append(dict(axial_metres=z,potential_volts=value,
            independent_axis_quadrature_volts=expected,relative_difference=abs(value/expected-1)))
    low,_=_switched_torus_receiver_kernel(design,onset,T,*quadrature_orders[0],source_batch_size)
    field,distances=_switched_torus_receiver_kernel(design,onset,T,*quadrature_orders[1],source_batch_size)
    high,_=_switched_torus_receiver_kernel(design,onset,T,*quadrature_orders[2],source_batch_size)
    def derivative(fn,h):
        return (fn(-2*h)-8*fn(-h)+8*fn(h)-fn(2*h))/(12*h)
    def curl(j):
        return np.array([j[2,1]-j[1,2],j[0,2]-j[2,0],j[1,0]-j[0,1]])
    probes=[]
    for point,delay in zip(([.1,.2,3.],[.2,-.1,5.],[.1,.2,11.]),delays):
        x=np.array(point);t=onset+delay
        pot,A,E,B=field(t,x);checks=[]
        for h in spatial_steps_metres:
            dt=h/c;derivatives=[[],[],[],[]]
            for j in range(3):
                for index in range(4):
                    derivatives[index].append(derivative(lambda z:field(t,x+z*np.eye(3)[j])[index],h))
            gp=np.array(derivatives[0]);ja,je,jb=(np.array(derivatives[i]).T for i in (1,2,3))
            dp=derivative(lambda z:field(t+z,x)[0],dt)
            da=derivative(lambda z:field(t+z,x)[1],dt)
            de=derivative(lambda z:field(t+z,x)[2],dt)
            db=derivative(lambda z:field(t+z,x)[3],dt)
            row=dict(spatial_step_metres=h,
                potential_E_relative=float(np.linalg.norm(E+gp+da)/scale),
                potential_B_relative=float(c*np.linalg.norm(B-curl(ja))/scale),
                divergence_E_relative=float(abs(np.trace(je))/(scale/ai)),
                divergence_B_relative=float(c*abs(np.trace(jb))/(scale/ai)),
                Faraday_relative=float(np.linalg.norm(curl(je)+db)/(scale/ai)),
                Ampere_relative=float(c*np.linalg.norm(curl(jb)-de/c**2)/(scale/ai)),
                Lorenz_gauge_relative=float(c*abs(np.trace(ja)+dp/c**2)/scale))
            assert max(v for name,v in row.items() if name.endswith('_relative'))<1e-8
            checks.append(row)
        changes=[]
        for fn in (low,high):
            pv,av,ev,bv=fn(t,x)
            change=dict(potential_volts=abs(pv-pot),electric_volts_per_metre=float(np.linalg.norm(ev-E)),
                        magnetic_tesla=float(np.linalg.norm(bv-B)))
            assert max(change['electric_volts_per_metre'],c*change['magnetic_tesla'])/scale<1e-9
            changes.append(change)
        probes.append(dict(receiver_coordinates_metres=point,seconds_after_onset=delay,
            electric_vector_volts_per_metre=E.tolist(),magnetic_vector_tesla=B.tolist(),
            checks=checks,main_to_source_change=changes[0],source_to_higher_change=changes[1]))
    for z,delay in zip((3.,5.,11.),delays):
        t=onset+delay
        def integrand(theta):
            zz=z-b*math.sin(theta);value=0.
            for a,sign in ((ai,1),(ao,-1)):
                d=math.hypot(a+b*math.cos(theta),zz)
                q,dq,_=_receiver_discharge_waveform(t-d/c,onset,T,Q)
                value+=sign*(q/d**3+dq/(c*d*d))*zz
            return k*2*math.cos(theta)**2*value/(2*math.pi)
        expected=quad(integrand,0,2*math.pi,epsabs=1e-4,epsrel=1e-12)[0]
        values=field(t,np.array([0.,0.,z]))
        assert abs(values[2][2]-expected)/scale<1e-11
        assert max(np.linalg.norm(values[1]),np.linalg.norm(values[3]))<1e-10
    return dict(
        weighting='w(theta)=2*cos(theta)^2, mean1 under dtheta/(2*pi); each endpoint retains total+/-q(t)',
        surface_charge_density='sigma_i=+/-q(t)*cos(theta)^2/[2*pi^2*b*(a_i+b*cos(theta))]',
        radial_density='D=sqrt(b^2-z^2); Jr=-qdot*D/(2*pi^2*r*b^2)*[1_[ai-D,ao-D](r)+1_[ai+D,ao+D](r)] for|z|<b, and0 otherwise. Edge bands have one branch, overlap has two',
        support='The same full annular-capsule closure as before. The density tends to0 at|z|=b and is bounded; jumps at the charged surfaces supply the required normal current',
        normalized_charge_error=abs(float(np.mean(weights))-1),
        weak_continuity_relative_error=abs(line/endpoints-1),
        radial_flux_over_negative_qdot=flux,
        current_squared_geometry_per_metre=C,
        independent_volume_quadrature_geometry_per_metre=direct,
        squared_current_time_integral_A2_seconds=time_norm,
        space_time_current_squared_A2_seconds_per_metre=C*time_norm,
        peak_transfer_current_amperes=peak_current,
        supremum_current_density_A_per_m2=max_density,
        peak_density_equation='The overlap supremum atz=0,r->ai+b from above is Imax/[pi^2*b*(ai+b)]; the inner single-branch supremum is Imax/[2*pi^2*b*(ai-b)]. The larger is retained',
        endpoint_surface_peaks=[dict(major_radius_metres=a,
            maximum_absolute_surface_charge_C_per_m2=Q/(2*math.pi**2*b*(a-b)),
            maximum_normal_E_jump_volts_per_metre=Q/(2*math.pi**2*b*(a-b)*eps))
            for a in (ai,ao)],
        surface_peak_scope='Maxima occur at the inward side theta=pi. The normal-field jump is sigma/epsilon0, not a one-sided field, conductor voltage or material-tolerance prediction',
        current_norm_equation='C=1/(2*pi^3*b^4) integral_-b^b D^2{ln[(ai+D)/(ai-D)]+4ln[(ao-D)/(ai+D)]+ln[(ao+D)/(ao-D)]}dz; integralJ^2dV=C*qdot^2; integralqdot^2dt=(4410/2431)Q^2/T',
        resistive_interpretation='Multiplying the finite space-time integral by an assumed constant resistivity in ohm metres gives an Ohmic diagnostic in joules. No resistivity, material response, drive circuit, Joule heat prediction or total energy cost is supplied; the impressed current is not a self-consistent conductor solution',
        density_samples=flux_checks,static_axis_checks=static_rows,
        dynamic_Maxwell_checks=probes,field_scale_volts_per_metre=scale,
        unchanged_quadrupole_C_m2=q2,
        initial_field_outside_20_2m_lower_bound_joules=3*k*q2*q2/(10*20.2**5),
        exterior_energy_scope='Exact positive multipole lower bound for this prescribed initial field; common minor-radius weighting cancels from the paired quadrupole. Not total field energy, radiation or a universal receiver bound',
        scope='The previous current-density divergence is removed in this explicit source family. Finite source-current norm and source conservation do not establish material feasibility, stable control, downstream transport or capture')


def proton_regularized_receiver_checks(ring,coils,torus,previous_receiver,
                                      incoming_trajectories=None,discharge_trajectories=None):
    """One bounded-current source, full new incoming launch, and one untuned 100ns discharge."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c,e,m_p
    design=_regularized_receiver_definition(previous_receiver,coils)
    if incoming_trajectories is None:
        incoming_trajectories=[proton_electrostatic_receiver_trace(
            ring,coils,torus,name,receiver_design=design) for name in ('main','source','fine')]
    def profiles(records):
        if isinstance(records,dict):
            records=[records[name] for name in ('main','source_refined','fine')]
        assert len(records)==3
        result={row['profile']:copy.deepcopy(row) for row in records}
        assert set(result)=={'main','source','fine'}
        return result
    incoming=profiles(incoming_trajectories)
    for name,row in incoming.items():
        assert row['design']==design
        assert row['receiver_plane_crossing'] is not None
        assert row['verification']['all_receiver_fields_present_from_initial_time']
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
        for key in ('maximum_total_work_residual_gamma',
                    'maximum_receiver_work_plus_potential_residual_gamma',
                    'maximum_static_energy_minus_control_work_residual_gamma'):
            assert row['verification'][key]<1e-8
        row['minimum_gamma_state']=min((row['minimum_gamma_state'],row['receiver_plane_crossing']),
                                      key=lambda state:state['gamma'])
        row['minimum_gamma_scope']='Lowest sampled/locally-refined/crossing candidate, not a certified global minimum'
        if name=='fine':
            plot=row['plot']
            by_time={t:{key:values[i] for key,values in plot.items()}
                     for i,t in enumerate(plot['lab_microseconds'])}
            for event in ('pulse_clear','receiver_plane_crossing','minimum_gamma_state','endpoint'):
                state=row[event]
                if state is None:continue
                time=state['lab_seconds']*1e6
                by_time[time]=dict(lab_microseconds=time,
                    receiver_axial_metres=state['receiver_coordinates_metres'][2],
                    receiver_transverse_metres=state['receiver_transverse_offset_metres'],
                    gamma=state['gamma'],kinetic_energy_MeV=state['kinetic_energy_eV']/1e6,
                    proton_proper_microseconds=state['proton_proper_seconds']*1e6,
                    x_metres=state['position_metres'][0],y_metres=state['position_metres'][1],
                    z_metres=state['position_metres'][2])
            row['plot']={key:[by_time[t][key] for t in sorted(by_time)] for key in plot}
        else:row.pop('plot',None)
    new_receiver=dict(design=design,trajectories=dict(
        main=incoming['main'],source_refined=incoming['source'],fine=incoming['fine']))
    definition=_switched_receiver_definition(new_receiver)
    source_audit=_regularized_receiver_source_checks(definition)
    if discharge_trajectories is None:
        discharge_trajectories=[proton_switched_receiver_trace(ring,coils,new_receiver,name)
                                for name in ('main','source','fine')]
    outgoing=profiles(discharge_trajectories)
    for name,row in outgoing.items():
        assert row['definition']==definition
        assert row['verification']['solver_success']
        assert row['verification']['maximum_added_work_residual_gamma']<1e-8
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
        if row['receiver_pulse_clear'] is not None:
            assert row['verification']['postclear_gamma_variation']<1e-10
            assert abs(row['verification']['clear_event_residual_seconds'])<1e-15
        if name!='fine':row.pop('plot',None)
    def differences(left,right,events):
        values={}
        for event in events:
            a,b=left[event],right[event]
            if a is None or b is None:
                assert a is b;continue
            values[event]=dict(position_metres=float(np.linalg.norm(
                np.array(a['position_metres'])-b['position_metres'])),
                lab_seconds=abs(a['lab_seconds']-b['lab_seconds']),
                proper_seconds=abs(a['proton_proper_seconds']-b['proton_proper_seconds']),
                gamma=abs(a['gamma']-b['gamma']))
            assert values[event]['position_metres']<1e-4
            assert values[event]['lab_seconds']<1e-12 and values[event]['gamma']<1e-8
        return values
    convergence={}
    for label,records,events in (
        ('incoming',incoming,('pulse_clear','receiver_plane_crossing','endpoint')),
        ('discharge',outgoing,('receiver_pulse_clear','endpoint'))):
        assert len({row['outcome'] for row in records.values()})==1
        convergence[label]=dict(
            source_quadrature_change=differences(records['main'],records['source'],events),
            timestep_and_tolerance_change=differences(records['source'],records['fine'],events))
    seed=incoming['fine']['receiver_plane_crossing'];final=outgoing['fine'];end=final['endpoint']
    energy=m_p*c*c
    complete_work=(end['gamma']-incoming['fine']['initial']['gamma']-
        (end['total_receiver_work_joules']+end['control_work_joules'])/energy)
    assert abs(complete_work)<1e-8
    field,distances=_switched_torus_receiver_kernel(design,definition['onset_seconds'],
        definition['duration_seconds'],32,16,32)
    local=np.array(seed['receiver_coordinates_metres'])
    onset_error=seed['lab_seconds']-distances(local)[0]/c-definition['onset_seconds']
    assert abs(onset_error)<1e-20
    # Recheck the old extraction-control wave, whose clearance must survive reflection.
    pulse=design['pulse_parameters'];position=np.array(seed['position_metres'])
    rr=math.hypot(position[0]-coils['parameters']['reference_orbit_radius_metres'],position[1])
    far=max(math.hypot(rr+pulse['control_loop_radius_metres'],position[2]-z)
        for z in (-pulse['control_loop_half_separation_metres'],pulse['control_loop_half_separation_metres']))
    far+=pulse['control_torus_minor_radius_metres']
    control_margin=seed['lab_seconds']-far/c-pulse['source_pulse_duration_seconds']
    assert control_margin>0
    clear=final['receiver_pulse_clear']
    if clear:
        values=field(clear['lab_seconds']+1e-12,np.array(clear['receiver_coordinates_metres']))
        assert values[0]==0 and all(np.max(np.abs(v))==0 for v in values[1:])
    # Keep the full incoming prefix and actual discharged continuation in a single
    # plot, excluding the separate static receiver's unphysical-for-this-run exit.
    common=list(incoming['fine']['plot'])
    merged={}
    for plot,before in ((incoming['fine']['plot'],True),(final['plot'],False)):
        for i,t in enumerate(plot['lab_microseconds']):
            if before and t>seed['lab_seconds']*1e6:continue
            merged[t]={key:plot[key][i] for key in common}
    continuous={key:[merged[t][key] for t in sorted(merged)] for key in common}
    assert continuous['lab_microseconds'][0]==0
    assert all(np.diff(continuous['lab_microseconds'])>0)
    assert all(np.diff(continuous['proton_proper_microseconds'])>=0)
    xi=np.linspace(-1,1,201)
    shape=2*np.sqrt(np.maximum(0.,1-xi*xi))/math.pi
    source_t=np.linspace(definition['onset_seconds']-.1e-7,definition['onset_seconds']+1.1e-7,241)
    charges,currents,_=_receiver_discharge_waveform(source_t,definition['onset_seconds'],1e-7,
                                                  design['inner_total_charge_coulombs'])
    return dict(
        model='Regularized finite receiver source: bounded radial transfer current, recomputed original gamma10 launch, one fixed100ns causal discharge',
        design=design,source_audit=source_audit,
        incoming=dict(trajectories=dict(main=incoming['main'],source_refined=incoming['source'],fine=incoming['fine']),
            scope='The changed weighted static charge is present from the original launch at t=0. This recomputes the entire incoming control/storage/receiver trajectory; the old uniform-charge center state is not reused. The static continuation after center is retained only as a separate comparison'),
        discharge=dict(definition=definition,trajectories=dict(main=outgoing['main'],
            source_refined=outgoing['source'],fine=final),
            scope='All discharge resolutions use the new fine incoming state and one fixed onset. They test the continuation conditional on that state; incoming refinement is reported separately. This is not an interval-certified whole-route error bound'),
        convergence=convergence,
        numerical_protocol='Incoming32/16,64/32 angular orders with a third run halving steps/tightening tolerances. Discharge16/8/16,32/16/32 azimuth/poloidal/radial orders, then halved steps/tighter tolerances. Position1e-4m,eventtime1e-12s,gamma1e-8 are numerical comparison criteria. One fixed source; no timing, duration, charge or launch scan',
        causal_connection=dict(
            onset_residual_seconds=onset_error,
            original_control_all_source_clearance_margin_seconds=control_margin,
            proof='Every fixed-source retarded time t-distance/c increases along a subluminal worldline. Onset is fixed by the full-support minimum distance at the new actual center state, so the earlier weighted-static trajectory is unchanged. After all-source discharge clearance, receiver fields cannot re-enter a future timelike path; the radiated wave remains elsewhere',
            preparation='Endpoint charges plus weighted radial impressed current obey distributional continuity. Synchronized preprogrammed source timing is assumed; no instantaneous command, feed design, current-carrier dynamics, constitutive material or causal preparation mechanism is demonstrated'),
        energy_accounting=dict(
            complete_launch_work_residual_gamma=complete_work,
            launch_kinetic_energy_eV=incoming['fine']['initial']['kinetic_energy_eV'],
            actual_new_center_kinetic_energy_eV=seed['kinetic_energy_eV'],
            final_kinetic_energy_eV=end['kinetic_energy_eV'],
            launch_to_final_kinetic_energy_removed_eV=incoming['fine']['initial']['kinetic_energy_eV']-end['kinetic_energy_eV'],
            old_uniform_center_gamma=previous_receiver['trajectories']['fine']['receiver_plane_crossing']['gamma'],
            new_weighted_center_gamma=seed['gamma'],
            static_weighted_endpoint_gamma=incoming['fine']['endpoint']['gamma'],
            equation='Delta gamma from original launch equals inherited magnetic-control electric work plus all weighted receiver work divided by mp*c^2. During the incoming static-receiver segment, gamma+e*Phi/(mp*c^2)-W_control/(mp*c^2) is conserved; after control clearance W_control is constant. During receiver discharge use the actual electric-work integral',
            scope='Single-proton work, finite integralJ^2, initial receiver field energy and emitted source radiation are separate quantities. The static comparison exits downstream; the discharged case may return upstream, so their endpoint comparison is not same-target efficiency'),
        outcome=final['outcome'],
        current_density_divergence_removed=True,
        downstream_transport_achieved=final['outcome']=='reached downstream plane',
        capture_demonstrated=False,
        plot=dict(continuous=continuous,weighted_static_comparison=copy.deepcopy(incoming['fine']['plot']),
            current_flux_profiles=dict(z_over_b=xi.tolist(),regularized_normalized_flux_density=shape.tolist(),
                old_uniform_formula='1/(pi*sqrt(1-xi^2)) for|xi|<1; integrable endpoint divergence, no finite endpoint sample',
                scope='Unit-integral radial-flux profile through a cylinder in the central overlap; radial edge bands carry a single branch'),
            source=dict(lab_microseconds=(source_t*1e6).tolist(),remaining_inner_charge_coulombs=charges.tolist(),
                        outward_transfer_current_amperes=(-currents).tolist())),
        scope='The divergent source-current norm is repaired within a precisely prescribed Maxwell source family. This does not show Ohmic self-consistency, practical source generation, capture, controlled rest, launch-error acceptance or reunion. The old unweighted-source radiation bound is not applied to this changed trajectory; no radiation reaction, source radiation budget or physically precise sub-eV turning-point prediction is established',
        sources=[
            dict(url='https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node145.html',
                 support='Retarded scalar/vector potentials, Eqs1746-1747. Weighted conserved sources and their trajectory here are our calculations'),
            dict(url='https://farside.ph.utexas.edu/teaching/jk1/lectures/node32.html',
                 support='Multipole expansion: paired quadrupole remains unchanged under the common normalized poloidal weight'),
            dict(url='https://farside.ph.utexas.edu/teaching/em/lectures/node56.html',
                 support='Electrostatic field energy used for the positive exterior multipole lower bound')],
        derivation_attribution='The2cos^2 source choice, branch density, finite current norm, new incoming match and discharge are team derivations checked against independent integration; no cited paper claims this apparatus')


def _receiver_duration_source_checks(definition):
    """Resolve the selected short pulse without requiring coarse screening nodes to pass."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c
    relative=copy.deepcopy(definition)
    relative['onset_seconds']=0.
    T=relative['duration_seconds'];design=relative['receiver']
    delays=(math.hypot(10,3)/c+.4*T,
            (math.hypot(10,5)+math.hypot(20,5))/(2*c)+.5*T,
            math.hypot(20,11)/c+.6*T)
    orders=((32,64,32),(32,128,32),(64,128,48))
    audit=_regularized_receiver_source_checks(relative,delays,(.0005,.00025),orders,128)
    for probe in audit['dynamic_Maxwell_checks']:
        probe['lower_to_reference_change']=probe.pop('main_to_source_change')
        probe['reference_to_higher_change']=probe.pop('source_to_higher_change')
    scale=audit['field_scale_volts_per_metre']
    center_delay=definition['initial']['lab_seconds']-definition['onset_seconds']
    points=[([.1,.2,3.],delays[0]),([.2,-.1,5.],delays[1]),([.1,.2,11.],delays[2]),
            ([.0002,0.,.2],center_delay+.3e-9),([.0002,0.,.2],center_delay+33.5e-9)]
    functions={pol:_switched_torus_receiver_kernel(design,0.,T,32,pol,32,128)[0]
               for pol in (8,16,64,128)}
    shifted=_switched_torus_receiver_kernel(design,definition['onset_seconds'],T,32,128,32,128)[0]
    unbatched=_switched_torus_receiver_kernel(design,0.,T,32,16,32)[0]
    probes=[]
    for position,time in points:
        x=np.array(position);reference=functions[128](time,x)
        rows=[]
        for pol,field in functions.items():
            value=field(time,x)
            electric=float(np.linalg.norm(value[2]-reference[2]))
            magnetic=float(np.linalg.norm(value[3]-reference[3]))
            rows.append(dict(poloidal_order=pol,electric_difference_volts_per_metre=electric,
                magnetic_difference_tesla=magnetic,
                maximum_E_or_cB_difference_over_source_scale=max(electric,c*magnetic)/scale))
        assert rows[2]['maximum_E_or_cB_difference_over_source_scale']<1e-9
        value=shifted(time+definition['onset_seconds'],x)
        time_shift_error=max(float(np.linalg.norm(value[2]-reference[2])),
                             c*float(np.linalg.norm(value[3]-reference[3])))/scale
        a,b=unbatched(time,x),functions[16](time,x)
        batching_error=max(float(np.linalg.norm(a[2]-b[2])),
                           c*float(np.linalg.norm(a[3]-b[3])))/scale
        assert time_shift_error<1e-9 and batching_error<1e-11
        probes.append(dict(receiver_coordinates_metres=position,seconds_after_onset=time,
            reference_electric_norm_volts_per_metre=float(np.linalg.norm(reference[2])),
            quadrature_rows=rows,time_origin_shift_relative_error=time_shift_error,
            summation_batching_relative_error=batching_error))
    return dict(identities_and_current_norm=audit,
        accepted_Maxwell_quadrature_orders=[list(x) for x in orders],
        spatial_derivative_steps_metres=[.0005,.00025],source_batch_size=128,
        poloidal_front_probes=probes,
        qualification='The8-duration grid is screening. Its8-poloidal and the16-poloidal source samples do not resolve every ns wavefront to the strict field-check tolerance. The selected trajectory uses64 poloidal nodes; separate128-node probes and independent front-split integration test it. No coarse-check assertion is weakened',
        conditioning='Finite-difference Maxwell checks translate the common time origin to onset0, avoiding subtraction of18microsecond offsets in picosecond derivatives. The actual absolute-time fields are separately compared at the same events. Batching only regroups the conserved-source integral',
        resource_failures='One unbatched64-poloidal trajectory and one parallel high-order field audit ended on NumPy temporary-array allocation failures before producing new trajectory results. They are computational failures, not reflection, transmission or source-model outcomes')


def proton_receiver_duration_scan_checks(ring,coils,regularized,
                                        screening_records=None,selected_records=None):
    """Eight fixed durations, one selected source refinement, and honest transport limits."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c,e,m_p
    receiver=dict(design=copy.deepcopy(regularized['design']),
                  trajectories=regularized['incoming']['trajectories'])
    durations=(100e-9,50e-9,20e-9,10e-9,5e-9,2e-9,1e-9,3.5e-9)
    selected_T=3.5e-9
    if screening_records is None:
        screening_records=[copy.deepcopy(regularized['discharge']['trajectories']['main'])]
        screening_records.extend(proton_switched_receiver_trace(ring,coils,receiver,'main',T)
                                 for T in durations[1:])
    screens=copy.deepcopy(screening_records)
    assert len(screens)==len(durations)
    source=regularized['source_audit'];Q=receiver['design']['inner_total_charge_coulombs']
    C=source['current_squared_geometry_per_metre']
    def cost(T):
        factor=100e-9/T;peak=630*Q/(256*T)
        return dict(peak_transfer_current_amperes=peak,
            peak_space_integral_current_squared_A2_per_metre=C*peak**2,
            space_time_current_squared_A2_seconds_per_metre=C*(4410/2431)*Q**2/T,
            supremum_current_density_A_per_m2=source['supremum_current_density_A_per_m2']*factor,
            current_and_integrated_J_squared_ratio_to_100ns=factor,
            peak_J_squared_ratio_to_100ns=factor**2)
    for T,row in zip(durations,screens):
        actual_T=row['definition']['duration_seconds']
        assert math.isclose(actual_T,T,rel_tol=1e-15,abs_tol=0.)
        assert row['profile']=='main' and row['definition']==_switched_receiver_definition(receiver,actual_T)
        row['duration_seconds']=actual_T;row['source_cost']=cost(actual_T);row.pop('plot',None)
        assert row['verification']['solver_success']
        assert row['verification']['maximum_added_work_residual_gamma']<1e-8
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
    downstream=[row for row in screens if row['outcome']=='reached downstream plane']
    assert downstream
    selected_screen=min(downstream,key=lambda row:row['endpoint']['kinetic_energy_eV'])
    assert selected_screen['duration_seconds']==selected_T
    if selected_records is None:
        selected_records=[copy.deepcopy(selected_screen)]
        selected_records.extend(proton_switched_receiver_trace(ring,coils,receiver,name,selected_T)
                                for name in ('source','fine','fine_poloidal'))
    if isinstance(selected_records,dict):selected_records=list(selected_records.values())
    traces={row['profile']:copy.deepcopy(row) for row in selected_records}
    assert set(traces)=={'main','source','fine','fine_poloidal'}
    definition=_switched_receiver_definition(receiver,selected_T)
    for name,row in traces.items():
        row.pop('duration_seconds',None);row.pop('source_cost',None)
        assert row['definition']==definition and row['verification']['solver_success']
        assert row['verification']['maximum_added_work_residual_gamma']<1e-8
        assert min(row['verification']['continuous_source_gap_lower_bounds_metres'])>0
        assert row['receiver_pulse_clear'] is not None
        assert row['verification']['postclear_gamma_variation']<1e-10
        assert abs(row['verification']['clear_event_residual_seconds'])<1e-15
        if name!='fine_poloidal':row.pop('plot',None)
    final=traces['fine_poloidal'];end=final['endpoint'];seed=definition['initial']
    assert final['outcome']=='reached downstream plane'
    source_audit=_receiver_duration_source_checks(definition)
    def differences(left,right):
        result={}
        for event in ('receiver_pulse_clear','endpoint'):
            a,b=left[event],right[event]
            result[event]=dict(position_metres=float(np.linalg.norm(
                np.array(a['position_metres'])-b['position_metres'])),
                lab_seconds=abs(a['lab_seconds']-b['lab_seconds']),
                proper_seconds=abs(a['proton_proper_seconds']-b['proton_proper_seconds']),
                gamma=abs(a['gamma']-b['gamma']),
                kinetic_energy_eV=abs(a['kinetic_energy_eV']-b['kinetic_energy_eV']))
        return result
    convergence=dict(screen_to_source=differences(traces['main'],traces['source']),
        source_to_smaller_steps=differences(traces['source'],traces['fine']),
        poloidal16_to64=differences(traces['fine'],final))
    independent=dict(
        receiver_pulse_clear=dict(lab_seconds=1.853671450748213e-5,
            proton_proper_seconds=1.9162307449522943e-6,
            position_metres=[-4985.5380797785965,-380.83478050809674,-1.4057273506351177e-14],
            gamma=1.0000804909386072,kinetic_energy_eV=75522.4011470908),
        endpoint=dict(lab_seconds=7.113412519412374e-5,
            proton_proper_seconds=5.450940815738009e-5,
            position_metres=[-5184.437099330031,-402.2979815575685,-5.227520101970059e-12],
            gamma=1.0000804909386072,kinetic_energy_eV=75522.40114709076))
    independent_change=differences(final,independent)
    for event in independent_change.values():
        assert event['position_metres']<1e-4 and event['lab_seconds']<1e-12
        assert event['gamma']<1e-8 and event['kinetic_energy_eV']<1.
    independent_audit=dict(reference_events=independent,differences=independent_change,
        source='Read-only independent reviewer integration, completed separately from this production helper; these retained reference events are not regenerated by cached reproduction',
        method='Cartesian position and mechanical momentum in SI, Radau rtol3e-12;24 azimuths offset0.237 grid spacings,32 Gauss-theta nodes per interval split at both endpoint leading/trailing wavefronts,32 Gauss nodes per active radial interval, four azimuths per memory batch. Independent256-node Cartesian storage-loop integral. The same source, seed, waveform and onset are used',
        time_coordinate='s=(t-t_center)/T; source phase=(t_center-onset)+T*s-distance/c; active maxdt=.05T and postclear maxdt100ns',
        comparison_scope='Numerical cross-check at the clear and downstream events, not an interval-certified whole-route error bound or physical receiver acceptance tolerance')
    # Step comparison at fixed source order. Coarse source changes are retained,
    # not passed through the former100ns event-time criterion.
    for event in convergence['source_to_smaller_steps'].values():
        assert event['position_metres']<1e-4 and event['lab_seconds']<1e-12
        assert event['gamma']<1e-8
    design=definition['receiver'];r=math.hypot(*seed['receiver_coordinates_metres'][:2])
    z=seed['receiver_coordinates_metres'][2];minor=design['minor_radius_metres']
    windows=[]
    for radius in (design['inner_major_radius_metres'],design['outer_major_radius_metres']):
        near=math.hypot(r-radius,z)-minor;far=math.hypot(r+radius,z)+minor
        first=definition['onset_seconds']+near/c;last=definition['onset_seconds']+far/c
        windows.append(dict(major_radius_metres=radius,nearest_distance_metres=near,
            farthest_distance_metres=far,first_onset_after_center_seconds=first-seed['lab_seconds'],
            last_onset_after_center_seconds=last-seed['lab_seconds'],
            last_end_after_center_seconds=last+selected_T-seed['lab_seconds']))
    spread=windows[-1]['last_onset_after_center_seconds']-windows[0]['first_onset_after_center_seconds']
    assert abs(windows[0]['first_onset_after_center_seconds'])<1e-20 and spread>30e-9
    norm=m_p*c*c;incoming=regularized['incoming']['trajectories']['fine']
    work_residual=end['gamma']-incoming['initial']['gamma']-(
        end['total_receiver_work_joules']+end['control_work_joules'])/norm
    assert abs(work_residual)<1e-8
    # Reuse precisely the incoming prefix up to the seed; the100ns and static
    # outgoing comparisons are separate trajectories, not part of this path.
    keys=list(incoming['plot']);merged={}
    for plot,before in ((incoming['plot'],True),(final['plot'],False)):
        for i,t in enumerate(plot['lab_microseconds']):
            if before and t>seed['lab_seconds']*1e6:continue
            merged[t]={key:plot[key][i] for key in keys}
    continuous={key:[merged[t][key] for t in sorted(merged)] for key in keys}
    assert continuous['lab_microseconds'][0]==0 and all(np.diff(continuous['lab_microseconds'])>0)
    assert all(np.diff(continuous['proton_proper_microseconds'])>=0)
    source_times=definition['onset_seconds']+np.linspace(-.1,1.1,121)*selected_T
    charge,current,_=_receiver_discharge_waveform(source_times,definition['onset_seconds'],selected_T,Q)
    final_records={name:row for name,row in traces.items() if name!='fine_poloidal'}
    return dict(
        model='Regularized retarded receiver: bounded discharge-duration scan with the original gamma10 launch, unchanged charge/geometry and causal onset',
        screening=dict(records=screens,durations_seconds=list(durations),
            initial_durations_seconds=list(durations[:4]),extension_durations_seconds=list(durations[4:7]),
            final_bracket_midpoint_seconds=selected_T,
            reused_100ns_baseline=True,number_of_new_duration_screens=7,
            selection='Lowest final kinetic energy among the screened downstream outcomes. Only that candidate is refined; this bounded grid is not an optimum or a fully converged ranking of all durations',
            turning_event_scope='Legacy screening profiles record positive-to-negative axial momentum crossings only. The selected fine_poloidal profile records both directions; absence of a screened forward reversal is not evidence of monotonic motion',
            outcome_counts={outcome:sum(row['outcome']==outcome for row in screens)
                            for outcome in sorted({row['outcome'] for row in screens})},
            deadline_after_center_seconds=definition['maximum_elapsed_after_center_seconds']),
        selected=dict(duration_seconds=selected_T,definition=definition,
            trajectory=final,refinement_records=final_records,source_cost=cost(selected_T),
            convergence=convergence,independent_SI_audit=independent_audit,
            numerical_scope='Source refinement and tighter-step changes are recorded separately. The coarse screen is not an accuracy reference.64-poloidal production is additionally checked by128-poloidal field probes and independent front-split SI integration; numerical differences are not receiver acceptance tolerances or interval error bounds'),
        source_audit=source_audit,
        causal_timing=dict(fixed_center_source_windows=windows,
            full_support_onset_arrival_spread_seconds=spread,
            first_and_last_refer_to='Support-closure arrival bounds at the fixed actual center event; current weights vanish at isolated edge points but have the same support closure',
            relation_to_turns='Each screen retains its actual axial turning events. Shortening a shared waveform does not remove the roughly34ns inner-to-outer front spread. Charge, induction, current and magnetic terms are all included; the outcome cannot be inferred from a static outer-ring attraction alone',
            onset_scope='This changes duration, not intervention before arrival. Onset equals actual seed time minus nearest full-source distance/c. Thus every earlier source-retarded time precedes onset, and the complete changed-static incoming path can be reused exactly conditional on its numerical seed',
            preparation='Source clocks are preprogrammed and synchronized. Continuity of the prescribed charge/current is exact, but no distributed trigger, material/current-carrier dynamics or drive circuit is constructed'),
        energy_accounting=dict(full_launch_work_residual_gamma=work_residual,
            initial_launch_kinetic_energy_eV=incoming['initial']['kinetic_energy_eV'],
            initial_continuation_kinetic_energy_eV=seed['kinetic_energy_eV'],
            final_kinetic_energy_eV=end['kinetic_energy_eV'],
            launch_to_final_energy_removed_eV=incoming['initial']['kinetic_energy_eV']-end['kinetic_energy_eV'],
            initial_source_exterior_field_lower_bound_joules=source['initial_field_outside_20_2m_lower_bound_joules'],
            source_cost_scaling='At fixed charge and geometry, peak current and integrated space-time J^2 scale as1/T; peak spatial J^2 scales as1/T^2. Initial electrostatic field energy is unchanged. No resistivity, Joule-heat number, emitted source energy or total drive budget is inferred',
            equation='Full Delta gamma=(inherited control electric work+all receiver electric work)/(mp*c^2). Time-dependent scalar potential is not a conserved energy; after all-source clearance the remaining static magnetic field does no work'),
        receiver_matching=dict(final_transverse_offset_metres=end['receiver_transverse_offset_metres'],
            final_plane_axial_metres=definition['downstream_axial_metres'],
            open_bore_radius_metres=design['inner_major_radius_metres']-minor,
            conclusion='Downstream passage at the recorded offset, allowing the actual intervening reversals. The target plane is not a detector or capture device; no perturbed-launch acceptance or stable storage at the receiver is proved'),
        outcome=final['outcome'],capture_demonstrated=False,controlled_rest_demonstrated=False,
        plot=dict(continuous=continuous,source=dict(lab_microseconds=(source_times*1e6).tolist(),
            remaining_inner_charge_coulombs=charge.tolist(),outward_transfer_current_amperes=(-current).tolist())),
        scope='A test-particle trajectory in prescribed conserved finite sources. The receiver source is bounded in current density but remains unsupported by a material response or preparation model. Radiation reaction, precise turning-point energies, source radiation, beam acceptance, final capture, reunion and human transport are outside this calculation. The old magnetic-only far-escape bound is not applied',
        sources=copy.deepcopy(regularized['sources']),
        derivation_attribution='Duration screens, source-wavefront refinements and costs are our calculations; the primary references establish field equations, not this apparatus')


def _receiver_radiation_temporal_amplitude(x):
    """Real midpoint-phase Beta(5,5) transform, with its removable zero limit."""
    import math
    from scipy.special import spherical_jn
    y=abs(float(x))/2
    if y<.25:
        value=term=1.
        for m in range(1,9):
            term*=-y*y/(2*m*(2*m+9))
            value+=term
        return value
    return float(945*spherical_jn(4,y)/y**4)


def _receiver_radiation_tail_bound(design,duration,cutoff):
    """Analytic positive-frequency tail bound; independent of numerical quadrature."""
    import math
    from scipy.constants import c
    Q=design['inner_total_charge_coulombs'];L=design['outer_major_radius_metres']-design['inner_major_radius_metres']
    mu=1/(design['vacuum_permittivity_SI']*c*c)
    assert cutoff*L/(c*duration)>=2
    coefficients=(1,40,2160,100800,2822400)
    total=sum(value*cutoff**(-9-2*j)/(9+2*j)*
        (8*math.log(cutoff*L/(c*duration))+4+8/(9+2*j))
        for j,value in enumerate(coefficients))
    return mu*Q*Q*c/(8*math.pi**2*duration)*945**2*2**10*total


def _receiver_radiation_spatial_amplitude(k,beta,design,poloidal_order=32):
    """Exact radial/azimuth Bessel reduction, numerical normalized poloidal average."""
    import math
    import numpy as np
    from scipy.special import j0
    angle=np.atleast_1d(beta)
    alpha=2*math.pi*(np.arange(poloidal_order)+.5)/poloidal_order
    weights=2*np.cos(alpha)**2/poloidal_order
    b=design['minor_radius_metres']
    ri=design['inner_major_radius_metres']+b*np.cos(alpha)
    ro=design['outer_major_radius_metres']+b*np.cos(alpha)
    height=b*np.sin(alpha)
    result=np.empty(len(angle))
    for first in range(0,len(angle),128):
        beta_part=angle[first:first+128]
        nu=np.cos(beta_part);v=k*np.sin(beta_part)[:,None]
        nonzero=v!=0
        radial=np.divide(j0(v*ri)-j0(v*ro),v,out=np.zeros((len(v),poloidal_order)),where=nonzero)
        small=np.abs(v*ro)<.1
        if np.any(small):
            series=np.zeros_like(radial)
            for m in range(1,6):
                series+=(-1)**(m+1)*v**(2*m-1)*(ro**(2*m)-ri**(2*m))/(4**m*math.factorial(m)**2)
            radial=np.where(small,series,radial)
        result[first:first+len(beta_part)]=np.sum(weights*np.cos(k*nu[:,None]*height)*radial,axis=1)
    return result


def _receiver_radiation_spectrum(design,duration,frequency_order=8,
                                  angular_multiplier=1,poloidal_order=32,
                                  cutoff=80.,panel_width=.25):
    """Stream the full isolated-source spectral integral, retaining all interference."""
    import math
    import numpy as np
    from scipy.constants import c
    from scipy.special import roots_legendre
    assert design['poloidal_weighting']=='cosine_squared'
    Q=design['inner_total_charge_coulombs'];eps=design['vacuum_permittivity_SI']
    mu=1/(eps*c*c);prefactor=mu*Q*Q/(8*math.pi**2*c*duration**3)
    outer=design['outer_major_radius_metres'];b=design['minor_radius_metres']
    low_band_x=2*c*duration/(outer+b)
    edges=np.unique(np.r_[np.arange(0,cutoff,panel_width),low_band_x,cutoff])
    fn,fw=roots_legendre(frequency_order)
    observer_grids={};panel_energies=[];xs=[];densities=[]
    for left,right in zip(edges[:-1],edges[1:]):
        frequencies=(left+right)/2+(right-left)*fn/2
        density=[]
        for x in frequencies:
            k=x/(c*duration)
            needed=32+k*(outer+b)
            order=angular_multiplier*2**int(math.ceil(math.log2(needed)))
            if order not in observer_grids:
                z,w=roots_legendre(order)
                beta=(z+1)*math.pi/4
                observer_grids[order]=(beta,w*math.pi/4)
            beta,wb=observer_grids[order]
            H=_receiver_radiation_spatial_amplitude(k,beta,design,poloidal_order)
            G=float(2*np.sum(wb*np.sin(beta)*np.cos(beta)**2*H**2))
            A=_receiver_radiation_temporal_amplitude(x)
            value=prefactor*x*x*A*A*G
            assert math.isfinite(value) and value>=0
            density.append(value)
        panel_energies.append(float((right-left)/2*np.dot(fw,density)))
        xs.extend(frequencies.tolist());densities.extend(density)
    total=math.fsum(panel_energies)
    low_band=math.fsum(value for value,right in zip(panel_energies,edges[1:]) if right<=low_band_x)
    return dict(frequency_order=frequency_order,angular_multiplier=angular_multiplier,
        poloidal_order=poloidal_order,observer_orders=sorted(observer_grids),
        source_parameters=dict(inner_major_radius_metres=design['inner_major_radius_metres'],
            outer_major_radius_metres=outer,minor_radius_metres=b,
            endpoint_charge_magnitude_coulombs=Q,pulse_duration_seconds=duration,
            vacuum_permittivity_SI=eps,poloidal_weighting=design['poloidal_weighting']),
        frequency_cutoff_x=cutoff,maximum_frequency_panel_width_x=panel_width,
        frequency_nodes_evaluated=len(xs),observer_source_batch_size=128,
        truncated_radiated_energy_joules=total,positive_low_band_energy_joules=low_band,
        plot=dict(frequency_x=xs,frequency_hertz=(np.array(xs)/(2*math.pi*duration)).tolist(),
            spectral_energy_joules_per_x=densities,
            spectral_energy_joules_per_hertz=(np.array(densities)*2*math.pi*duration).tolist(),
            cumulative_frequency_x=edges.tolist(),
            cumulative_energy_joules=np.r_[0.,np.cumsum(panel_energies)].tolist()))


def proton_receiver_radiation_spectrum_checks(duration_scan,quadrature_results=None):
    """Full conserved receiver spectrum with separate numerical and analytic error checks."""
    import copy
    import math
    import numpy as np
    from scipy.constants import c
    from scipy.special import roots_legendre
    definition=duration_scan['selected']['definition'];design=copy.deepcopy(definition['receiver'])
    T=definition['duration_seconds'];Q=design['inner_total_charge_coulombs']
    ai,ao,b=(design['inner_major_radius_metres'],design['outer_major_radius_metres'],design['minor_radius_metres'])
    mu=1/(design['vacuum_permittivity_SI']*c*c);cutoff=80.
    # Characteristic-function normalization and its phase-independent real part.
    temporal=[]
    for x in (0.,1e-6,.1,.5,1.,5.,10.,30.,80.,200.):
        expected=quad(lambda u:630*u**4*(1-u)**4*math.cos(x*(u-.5)),0,1,
                      epsabs=1e-13,epsrel=1e-12)[0]
        value=_receiver_radiation_temporal_amplitude(x)
        assert abs(value-expected)<1e-11
        temporal.append(dict(frequency_x=x,amplitude=value,independent_time_quadrature=expected,
                             absolute_difference=abs(value-expected)))
    # Exact finite-polynomial identity behind the temporal tail envelope.
    p=np.polynomial.polynomial.polymul([1,0,-45,0,105],[1,0,-45,0,105])
    q=np.polynomial.polynomial.polymul([0,10,0,-105],[0,10,0,-105])
    p[:len(q)]+=q
    assert np.array_equal(p,[1,0,10,0,135,0,1575,0,11025])
    spatial=[]
    alpha=2*math.pi*(np.arange(64)+.5)/64;phi=2*math.pi*(np.arange(128)+.5)/128
    weights=2*np.cos(alpha)**2/64
    ri=ai+b*np.cos(alpha);ro=ao+b*np.cos(alpha);z=b*np.sin(alpha)
    for k,beta in ((1e-5,.3),(.03,.7),(.3,1.1),(1.,.5)):
        projection=k*math.sin(beta)*np.cos(phi)
        radial=(ao-ai)*np.exp(-1j*projection[:,None]*(ri+ro)/2)*np.sinc(projection[:,None]*(ao-ai)/(2*math.pi))
        phase=np.exp(-1j*k*math.cos(beta)*z)
        vector=-np.mean(np.cos(phi)*np.sum(weights*phase*radial,axis=1))
        H=float(_receiver_radiation_spatial_amplitude(k,[beta],design,64)[0])
        assert abs(vector-1j*H)<1e-11
        spatial.append(dict(wavenumber_per_metre=k,observer_polar_angle_radians=beta,
            bessel_amplitude_metres=H,independent_Cartesian_amplitude_real_metres=float(vector.real),
            independent_Cartesian_amplitude_imaginary_metres=float(vector.imag),
            absolute_difference_metres=float(abs(vector-1j*H))))
    z,w=roots_legendre(128);beta=(z+1)*math.pi/4;k=1e-6
    H=_receiver_radiation_spatial_amplitude(k,beta,design,64)
    G=float(2*np.sum(w*math.pi/4*np.sin(beta)*np.cos(beta)**2*H**2))
    low_k_ratio=G/((ao*ao-ai*ai)**2*k*k/60)
    assert abs(low_k_ratio-1)<1e-8
    g_second_norm=quad(lambda u:(2520*u*u*(1-u)**2*(3-14*u+14*u*u))**2,
                       0,1,epsabs=1e-9,epsrel=1e-12)[0]
    assert abs(g_second_norm/(272160/143)-1)<1e-12
    z,w=roots_legendre(256);x=100*(z+1)
    fourier_norm=float(100*np.sum(w*x**4*np.array([
        _receiver_radiation_temporal_amplitude(value)**2 for value in x])))
    parseval_tail=945**2*2**10*sum(value*200.**(-5-2*j)/(5+2*j)
        for j,value in enumerate((1,40,2160,100800,2822400)))
    assert -1e-8<=math.pi*g_second_norm-fourier_norm<=parseval_tail+1e-8
    long_pulse=mu*Q*Q*(ao*ao-ai*ai)**2/(480*math.pi*c**3)*g_second_norm
    profiles=((8,1,32),(12,1,32),(12,2,64))
    if quadrature_results is None:
        quadrature_results=[_receiver_radiation_spectrum(design,T,*profile,cutoff) for profile in profiles]
    rows=copy.deepcopy(quadrature_results)
    assert len(rows)==3
    source_parameters=dict(inner_major_radius_metres=ai,outer_major_radius_metres=ao,
        minor_radius_metres=b,endpoint_charge_magnitude_coulombs=Q,
        pulse_duration_seconds=T,vacuum_permittivity_SI=design['vacuum_permittivity_SI'],
        poloidal_weighting=design['poloidal_weighting'])
    for row,profile in zip(rows,profiles):
        assert row['source_parameters']==source_parameters
        assert (row['frequency_order'],row['angular_multiplier'],row['poloidal_order'])==profile
        assert row['frequency_cutoff_x']==cutoff
        assert row['truncated_radiated_energy_joules']>0
    reference=rows[-1]['truncated_radiated_energy_joules']
    frequency_change=abs(rows[1]['truncated_radiated_energy_joules']-rows[0]['truncated_radiated_energy_joules'])
    spatial_change=abs(reference-rows[1]['truncated_radiated_energy_joules'])
    assert max(frequency_change,spatial_change)/reference<1e-8
    omega=2*c/(ao+b);moment=Q*(ao*ao-ai*ai)/2
    lower=mu*moment**2*omega**5/(2400*math.pi**2*c**3)*math.cos(2*b/(ao+b))**2*math.cos(omega*T/2)**2
    assert rows[-1]['positive_low_band_energy_joules']>=lower>3.676e8
    tail=_receiver_radiation_tail_bound(design,T,cutoff)
    assert tail>0 and tail/reference<1e-6
    # Total initial energy is bounded analytically, without a singular self-energy quadrature.
    with mp.workdps(50):
        ia,oa,bb,qq,ee=map(lambda v:mp.mpf(str(v)),
            (ai,ao,b,Q,design['vacuum_permittivity_SI']))
        def self_bound(a):
            A=a-bb
            return mp.mpf(3)/4*(mp.asinh(A/bb)/A+mp.asinh(bb/A)/bb)
        b_inner,b_outer=self_bound(ia),self_bound(oa)
        cross_lower=1/(ia+oa+2*bb)
        initial_upper=qq**2/(8*mp.pi*ee)*(b_inner+b_outer-2*cross_lower)
        # A deliberately coarser outward decimal, well above the50-digit value.
        upper_rounded=float(mp.ceil(initial_upper/1000)*1000)
    initial_lower=duration_scan['energy_accounting']['initial_source_exterior_field_lower_bound_joules']
    lower_rounded=math.floor(initial_lower/1000)*1000.
    assert 0<lower_rounded<initial_lower<float(initial_upper)<upper_rounded<reference
    source_work_lower=reference-upper_rounded
    source_work_upper=reference+tail-lower_rounded
    plot=rows[-1].pop('plot')
    for row in rows[:-1]:row.pop('plot',None)
    return dict(
        model='Full radiation spectrum of the isolated prescribed regularized receiver discharge in vacuum; unchanged accepted3.5ns source',
        parameters=dict(inner_major_radius_metres=ai,outer_major_radius_metres=ao,
            minor_radius_metres=b,endpoint_charge_magnitude_coulombs=Q,pulse_duration_seconds=T,
            source_onset_seconds=definition['onset_seconds'],vacuum_permeability_SI=mu,
            vacuum_permittivity_SI=design['vacuum_permittivity_SI'],light_speed_SI=c),
        equations=dict(temporal='qdot_hat=-Q*exp[i*omega*(onset+T/2)]*A(x), x=omega*T; A=945*j4(x/2)/(x/2)^4, A(0)=1. j4 is spherical Bessel, not ordinary J4',
            spatial='H=<2cos(alpha)^2*cos(k*b*sin(alpha)*cos(beta))*[J0(k*ri*sin(beta))-J0(k*ro*sin(beta))]/[k*sin(beta)]>alpha; ri/o=ai/o+b*cos(alpha)',
            angular='G(k)=2*integral_0^(pi/2)sin(beta)*cos(beta)^2*H^2 d beta; all intra/inter-surface interference remains inside H^2',
            radiation='W=mu0*Q^2/(8*pi^2*c*T^3)*integral_0^infinity x^2*A(x)^2*G(x/(c*T)) dx. Positive frequencies for real fields, no extra sinusoidal time-average factor',
            source='rho=q(t)<delta_i-delta_o>; J=-qdot(t)<er*integral_ri^ro delta(x-r*er-b*sin(alpha)*ez)dr>, normalized weighted angle average. Continuity incorporates the endpoint charges; adding separate charge radiation would double count'),
        radiation=dict(evaluated_spectrum_joules=reference,evaluated_spectrum_gigajoules=reference/1e9,
            omitted_positive_tail_upper_bound_joules=tail,tail_bound_over_evaluated_spectrum=tail/reference,
            cutoff_dimensionless_frequency_x=cutoff,cutoff_angular_frequency_per_second=cutoff/T,
            cutoff_frequency_hertz=cutoff/(2*math.pi*T),
            independently_derived_low_band_lower_bound_joules=lower,
            evaluated_low_band_joules=rows[-1]['positive_low_band_energy_joules'],
            interpretation='The finite-band value estimates total emitted energy with a separately bounded nonnegative tail. Quadrature refinement is numerical evidence, not an outward-rounded certificate; do not present the finite-band estimate plus tail as a rigorous whole-spectrum interval without also certifying that quadrature'),
        tail_bound=dict(temporal_squared_envelope='A(x)^2<=945^2*2^10*x^-10*(1+40/x^2+2160/x^4+100800/x^6+2822400/x^8)',
            spatial_envelope='|H|<=min(L,2/[k*sin(beta)]), L=ao-ai. Splitting the angular integral at nu0=sqrt(1-4/(kL)^2) gives G<=[8*log(kL)+4]/k^2 for kL>=2',
            integrated_bound='mu0*Q^2*c*945^2*2^10/(8*pi^2*T) times sum_j c_j*X^(-9-2j)/(9+2j)*[8log(XL/(cT))+4+8/(9+2j)], c=(1,40,2160,100800,2822400)',
            source_condition_cutoff_times_L_over_cT=cutoff*(ao-ai)/(c*T),
            scope='Analytic inequality over every omitted frequency. It includes the full spatial transform rather than a long-wavelength multipole approximation; ordinary floating evaluation of bound constants is not interval arithmetic'),
        verification=dict(temporal_transform_checks=temporal,spatial_transform_checks=spatial,
            low_frequency_G_to_quadrupole_limit_ratio=low_k_ratio,
            temporal_Parseval_check=dict(exact_g_second_squared_integral=272160/143,
                numerical_time_integral=g_second_norm,
                evaluated_positive_frequency_x4_A2_integral_to_200=fourier_norm,
                exact_full_positive_frequency_integral=math.pi*g_second_norm,
                omitted_positive_tail_upper_bound=parseval_tail),
            strict_long_pulse_limit_W_times_T5_joule_seconds5=long_pulse,
            long_pulse_scope='Analytic T->infinity normalization benchmark. The electrically large3.5ns source is evaluated with its full transform, not this multipole approximation',
            quadrature_results=rows,frequency_refinement_change_joules=frequency_change,
            angular_and_poloidal_refinement_change_joules=spatial_change,
            maximum_refinement_change_over_spectrum=max(frequency_change,spatial_change)/reference,
            protocol='One frequency at a time,128 observation-angle rows per source batch. Frequency8->12 Gauss per panel isolates frequency error; doubled observer orders and32->64 poloidal nodes then test spatial quadrature. No large Cartesian field grid or proton trajectory rerun'),
        energy_accounting=dict(initial_exterior_field_energy_lower_bound_joules=duration_scan['energy_accounting']['initial_source_exterior_field_lower_bound_joules'],
            total_initial_field_energy_upper_bound_joules=float(initial_upper),
            conservatively_rounded_total_initial_upper_bound_joules=upper_rounded,
            conservatively_rounded_initial_lower_bound_joules=lower_rounded,
            initial_self_inverse_distance_upper_bounds_per_metre=[float(b_inner),float(b_outer)],
            initial_cross_inverse_distance_lower_bound_per_metre=float(cross_lower),
            initial_energy_bound_equations='A=a-b; B(a)=3/4*[asinh(A/b)/A+asinh(b/A)/b]; U0<=ke*Q^2/2*[B(ai)+B(ao)-2/(ai+ao+2b)]',
            initial_energy_bound_derivation='The normalized poloidal weight is1+cos2alpha, with autocorrelation1+.5cos2delta<=1.5. Exact same-torus distance squared is4b^2sin^2(delta_alpha/2)+4r*rprime*sin^2(delta_phi/2), withr,rprime>=a-b. Folding angles and using sin(u)>=2u/pi on[0,pi/2] gives B(a). Opposite-torus distances never exceed ai+ao+2b, giving the cross-term lower bound before its negative energy sign',
            initial_energy_bound_scope='Analytic inequalities for the full weighted charged-sheet source.50-digit evaluation and conservative decimal rounding are recorded; no exact totalU0, conductor equilibrium or material energy is computed. Independently audited source-team bound, reproduced here',
            conditional_source_to_EM_work_lower_estimate_joules=source_work_lower,
            conditional_source_to_EM_work_upper_estimate_joules=source_work_upper,
            conditional_source_work_sign='Positive for the evaluated spectrum and proven stored-energy bounds',
            total_initial_receiver_field_energy_computed=False,driver_work_computed=False,
            identity='For the isolated receiver, -integral J dot E_receiver dVdt=Wrad-U0, where U0 is total initial electrostatic field energy, not its exterior lower bound',
            scope='The work estimates combine numerical finite-band radiation, its analytic positive tail, and analytic lower/upper bounds on totalU0. Because the finite-band quadrature is not interval-certified, this is a conditional numerical source-work bracket, not a rigorous standalone work theorem. It concerns the imposed receiver current and its own Maxwell field. Driver material/carrier changes, losses, energy recovery and other-source cross terms are not solved. Feeds, shields or additional currents may alter radiation and require a new total-source calculation'),
        plot=plot,
        scope='Retarded Maxwell radiation from the same impressed conserved source used for the accepted trajectory. This does not establish material response, current generation, energy recovery, a total apparatus budget, capture or time travel',
        sources=[dict(url='https://arxiv.org/pdf/1801.03263',support='Wang et al., Eq2.4, printedpage3: exact transverse current Fourier transform of the Maxwell far field'),
            dict(url='https://dlmf.nist.gov/10.49',support='NIST DLMF10.49.2: exact elementary spherical-Bessel formulas, used to derive the temporal tail envelope'),
            dict(url='https://arxiv.org/pdf/physics/0506053',support='Jackson, Eq1: Poynting energy bookkeeping; radiation is distinct from source work and stored field energy')],
        derivation_attribution='The weighted receiver transform, Beta(5,5) waveform, integrated pulse spectrum and sharper tail bound are our derivations. SOURCE_FEASIBILITY.md supplied the independently reviewed exact spatial/Fourier expression and low-band benchmark')


def proton_staged_rf_deceleration_checks(ring):
    """Synchronous averaged ramp with signed induction work; no beam capture model."""
    import math
    import numpy as np
    import sympy as sp
    from scipy.integrate import solve_ivp

    # The old static weak-focusing proof is not used for this time-dependent field.
    r, z, t = sp.symbols('r z t', real=True)
    rs, bstart, cs = sp.symbols('R Bstart c', positive=True)
    bdot = sp.symbols('Bdot', real=True)
    ns = sp.Rational(2, 5)
    ap = r*(1+ns/2)/2-ns*r**3/8+ns*r*z*z/2
    br, bz = -ns*r*z, 1+ns/2-ns*r*r/2+ns*z*z
    eramp = -bdot*rs*ap
    brr, bzz = (bstart+bdot*t)*br, (bstart+bdot*t)*bz
    maxwell = [sp.diff(r*brr, r)/(rs*r)+sp.diff(bzz, z)/rs,
        -sp.diff(eramp, z)/rs+sp.diff(brr, t),
        sp.diff(r*eramp, r)/(rs*r)+sp.diff(bzz, t),
        (sp.diff(brr, z)-sp.diff(bzz, r))/rs-sp.diff(eramp, t)/cs**2,
        sp.diff(r*ap, r)/r-bz, -sp.diff(ap, z)-br]
    assert all(sp.simplify(value) == 0 for value in maxwell)
    zeta_symbolic = sp.simplify(ap.subs({r:1, z:0}))
    assert zeta_symbolic == sp.Rational(11, 20)
    gs = sp.symbols('gamma', positive=True)
    us = sp.sqrt(gs*gs-1)
    assert sp.simplify(sp.diff(gs/us, gs, 2)-3*gs/us**5) == 0
    assert sp.simplify(sp.diff(1/us, gs, 2)-(2*gs*gs+1)/us**5) == 0

    with mp.workdps(65):
        constants = ring['constants']
        c, charge, mass, eps = [mp.mpf(str(constants[key])) for key in
            ('c', 'elementary_charge', 'proton_mass_kg', 'vacuum_permittivity')]
        radius = mp.mpf(str(ring['nominal']['bending_radius_metres']))
        circumference = 2*mp.pi*radius
        rest_eV = mass*c*c/charge
        gi, kinetic_final, decrement = mp.mpf(10), mp.mpf('5.3e6'), mp.mpf(250)
        gf = 1+kinetic_final/rest_eV
        assert gi > gf > 1 and decrement > 0
        ui, uf = mp.sqrt(gi*gi-1), mp.sqrt(gf*gf-1)
        zeta = mp.mpf(11)/20
        field_initial = mass*c*ui/(charge*radius)
        field_slope = -decrement/(radius*circumference)
        scale = circumference*rest_eV/(c*decrement)
        lab = scale*(ui-uf)
        proper = scale*(mp.acosh(gi)-mp.acosh(gf))
        turns = rest_eV*(gi-gf)/decrement
        kinetic_initial = (gi-1)*rest_eV
        delta_kinetic = kinetic_final-kinetic_initial
        induction_turn = -zeta*decrement
        induction_total = zeta*delta_kinetic
        induction_Ephi = -field_slope*radius*zeta
        assert field_slope < 0 < induction_Ephi
        assert abs(-charge*induction_Ephi*circumference/charge-induction_turn) < mp.mpf('1e-55')

        # Free-space Lienard diagnostic of the smooth reference, including pdot.
        radiation_coefficient = charge/(3*eps*radius)
        momentum_dot = -decrement*charge/circumference
        parallel_power = charge**2*momentum_dot**2/(6*mp.pi*eps*mass**2*c**3)
        transverse_total = rest_eV/decrement*radiation_coefficient*(ui**5-uf**5)/5
        parallel_total = parallel_power*lab/charge
        radiation_total = transverse_total+parallel_total
        rf_total = (1-zeta)*delta_kinetic+radiation_total
        ledger = induction_total+rf_total-radiation_total-delta_kinetic
        assert abs(ledger) < mp.mpf('1e-50')

        def state(gamma):
            u = mp.sqrt(gamma*gamma-1)
            beta = u/gamma
            elapsed = scale*(ui-u)
            clock = scale*(mp.acosh(gi)-mp.acosh(gamma))
            frequency = beta*c/circumference
            transverse = radiation_coefficient*gamma*u**3
            radiation = transverse+parallel_power/(charge*frequency)
            rf_work = -(1-zeta)*decrement+radiation
            assert rf_work < 0 and abs(induction_turn+rf_work-radiation+decrement) < mp.mpf('1e-55')
            return dict(lab_seconds=float(elapsed), proper_seconds=float(clock),
                kinetic_energy_MeV=float((gamma-1)*rest_eV/1e6), gamma=float(gamma),
                speed_over_c=float(beta), momentum_GeV_over_c=float(rest_eV*u/1e9),
                field_tesla=float(mass*c*u/(charge*radius)), rf_frequency_hz=float(frequency),
                cumulative_turns=float(rest_eV*(gi-gamma)/decrement),
                induction_work_per_turn_eV=float(induction_turn), rf_work_per_turn_eV=float(rf_work),
                radiation_loss_per_turn_eV=float(radiation),
                rf_effective_decelerating_crest_voltage=float(-rf_work))

        # Independent energy-parameter quadratures and a lab-time energy/clock ODE.
        gamma_breaks = [gf, mp.mpf('1.1'), mp.mpf(2), gi]
        integral_lab = scale*mp.quad(lambda g:g/mp.sqrt(g*g-1), gamma_breaks)
        integral_proper = scale*mp.quad(lambda g:1/mp.sqrt(g*g-1), gamma_breaks)
        integral_radiation = rest_eV/decrement*mp.quad(
            lambda g:radiation_coefficient*g*(g*g-1)**mp.mpf('1.5'), gamma_breaks)
        integral_errors = [abs(integral_lab-lab), abs(integral_proper-proper),
                           abs(integral_radiation-transverse_total)]
        assert max(integral_errors) < mp.mpf('1e-48')

        def turn_clock_bracket(step, is_proper, direct_check=False):
            number = rest_eV*(gi-gf)/step
            full = int(mp.floor(number)); fraction = number-full
            h = step/rest_eV
            gn = gi-h*full
            clock_scale = circumference*rest_eV/(c*step)
            def orbit_period(g):
                return circumference/(c*mp.sqrt(g*g-1))*(1 if is_proper else g)
            def elapsed(g):
                return clock_scale*((mp.acosh(gi)-mp.acosh(g)) if is_proper
                                    else (ui-mp.sqrt(g*g-1)))
            def derivative(g):
                return h*circumference/c*(g if is_proper else 1)/(g*g-1)**mp.mpf('1.5')
            # Composite trapezoid error is integral theta(1-theta)F''/2.
            # F''>0 and theta(1-theta)/2<=1/8; no sum over the real turn count.
            left_lower = elapsed(gn)+(orbit_period(gi)-orbit_period(gn))/2+fraction*orbit_period(gn)
            width = (derivative(gn)-derivative(gi))/8
            right_lower = left_lower+orbit_period(gn)-orbit_period(gi)+fraction*(orbit_period(gf)-orbit_period(gn))
            continuum = elapsed(gf)
            assert width > 0 and left_lower+width < continuum < right_lower
            result = dict(continuum_seconds=float(continuum),
                start_energy_clock_enclosure_seconds=[float(left_lower), float(left_lower+width)],
                end_energy_clock_enclosure_seconds=[float(right_lower), float(right_lower+width)],
                start_energy_bias_enclosure_seconds=[float(left_lower-continuum), float(left_lower+width-continuum)],
                end_energy_bias_enclosure_seconds=[float(right_lower-continuum), float(right_lower+width-continuum)],
                enclosure_width_seconds=float(width),
                monotonicity_absolute_error_bound_seconds=float(orbit_period(gf)-orbit_period(gi)))
            if direct_check:
                direct_left = mp.fsum(orbit_period(gi-h*j) for j in range(full))+fraction*orbit_period(gn)
                direct_right = mp.fsum(orbit_period(gi-h*j) for j in range(1, full+1))+fraction*orbit_period(gf)
                assert left_lower <= direct_left <= left_lower+width
                assert right_lower <= direct_right <= right_lower+width
                result.update(direct_start_energy_seconds=float(direct_left), direct_end_energy_seconds=float(direct_right),
                    integer_turns=full, final_turn_fraction=float(fraction))
            return result

        step_bounds = dict(lab=turn_clock_bracket(decrement, False),
                           proper=turn_clock_bracket(decrement, True))
        for key, bound in (('lab', 2.756e-6), ('proper', 3.034e-6)):
            assert max(abs(value) for name in ('start_energy_bias_enclosure_seconds',
                'end_energy_bias_enclosure_seconds') for value in step_bounds[key][name]) < bound
            step_bounds[key]['conservative_absolute_error_bound_seconds'] = bound
        test_step = rest_eV*(gi-gf)/mp.mpf('1024.3')
        direct_bounds = dict(lab=turn_clock_bracket(test_step, False, True),
                             proper=turn_clock_bracket(test_step, True, True))
        stage_energies = [kinetic_initial, mp.mpf('1e9'), mp.mpf('1e8'), mp.mpf('2e7'), kinetic_final]
        stages = [state(1+kinetic/rest_eV) for kinetic in stage_energies]
        plot_gammas = [mp.sqrt(1+(ui-(ui-uf)*mp.mpf(j)/160)**2) for j in range(161)]
        plot_gammas += [1+kinetic/rest_eV for kinetic in stage_energies[1:-1]]
        plot_states = [state(g) for g in sorted(set(plot_gammas), reverse=True)]
        plot = {key:[row[key] for row in plot_states] for key in plot_states[0]}
        # Float conversion can put the analytically zero initial clock below zero by 1ulp.
        plot['lab_seconds'][0] = plot['proper_seconds'][0] = plot['cumulative_turns'][0] = 0.
        source_radiation_error = abs(stages[0]['radiation_loss_per_turn_eV']/ring['radiation']['energy_per_turn_eV']-1)
        assert source_radiation_error < 2e-14
        updiff, lab_float, proper_float = float(ui-uf), float(lab), float(proper)
        gammaf_float, turns_float = float(gf), float(turns)
        maximum_turn_fraction = float(decrement/(rest_eV*(gf*gf-1)/gf))
        result = dict(
            model='Ideal synchronous circular ramp from already-stored gamma10 proton to5.3MeV; effective RF gaps plus an explicit local induced electric field. This is a separate deceleration benchmark, not a continuation of the receiver trajectory',
            parameters=dict(constants=dict(constants), bending_radius_metres=float(radius), circumference_metres=float(circumference),
                weak_focusing_index=.4, initial_gamma=10., final_kinetic_energy_eV=float(kinetic_final),
                net_energy_removed_per_turn_eV=float(decrement), rf_harmonic=1, rest_energy_eV=float(rest_eV)),
            equations=dict(energy='dE/dN=-D; dN/dt=beta*c/C; dp/dt=-D_J/C; p=e*R*B0',
                clocks='t=C*(pc_i-pc_f)/(c*D); tau=C*(m*c^2)*(acosh(gamma_i)-acosh(gamma_f))/(c*D), using one consistent energy unit',
                rf='f_RF=beta*c/C; signed work W_RF/e=-(1-zeta)*D+Urad is negative, while the effective peak amplitude Vhat=(1-zeta)*D-Urad is positive, in volts for a proton when D,Urad are in eV. Phase phi_s=-pi/2 gives W_RF/e=Vhat*sin(phi_s)',
                signed_work='W_ind=-zeta*D; W_RF=-(1-zeta)*D+Urad; W_radiation=-Urad; their sum=-D',
                induction_fraction='zeta=Aphi(R,0)/(B0*R)=1/2+n/8=11/20 for this specified field; it is not a universal RF/induction partition'),
            reference=dict(lab_elapsed_seconds=lab_float, proton_proper_elapsed_seconds=proper_float,
                lab_minus_proper_seconds=float(lab-proper), lab_to_proper_ratio=float(lab/proper),
                total_turns=turns_float, complete_turns=int(mp.floor(turns)), last_turn_fraction=float(turns-mp.floor(turns)),
                initial_kinetic_energy_eV=float(kinetic_initial), final_gamma=gammaf_float,
                maximum_fractional_momentum_change_per_turn=maximum_turn_fraction),
            induction=dict(field_definition='Aphi=B0(t)*R*[r*(1+n/2)/2-n*r^3/8+n*r*z^2/2], with r,z in units R; B=curl A and E_ind=-partial_t A',
                B_shape='Br/B0=-n*r*z; Bz/B0=1+n/2-n*r^2/2+n*z^2; Bphi=0',
                field_initial_tesla=float(field_initial), field_final_tesla=stages[-1]['field_tesla'],
                field_slope_tesla_per_second=float(field_slope), Ephi_at_reference_volts_per_metre=float(induction_Ephi),
                flux_through_reference_disk_weber_initial=float(2*mp.pi*radius**2*field_initial*zeta),
                work_per_turn_eV=float(induction_turn),
                sign='B0 decreases, hence E_ind points in positive azimuth; the positive proton moves clockwise and qE_ind dot v is negative',
                scope='During the linear ramp this polynomial field satisfies the local vacuum Maxwell equations exactly. It is not a finite-source magnet design. Ramp start/end transitions, source energy, finite gaps and off-reference dynamics are not supplied'),
            energy_accounting=dict(kinetic_energy_change_eV=float(delta_kinetic), induction_work_eV=float(induction_total),
                effective_rf_work_eV=float(rf_total), radiated_energy_eV=float(radiation_total),
                signed_balance_residual_eV=float(abs(ledger)),
                scope='Signed work on one prescribed reference proton. Radiation compensation is included once in RF work. This is not wall-plug power, energy recovery, magnet/cavity stored energy or a complete apparatus budget'),
            radiation=dict(transverse_per_turn_formula_eV='e*gamma*(gamma^2-1)^(3/2)/(3*epsilon0*R)',
                transverse_integral_formula_eV='(m*c^2/e)/D * e/(15*epsilon0*R) * [(gamma_i^2-1)^(5/2)-(gamma_f^2-1)^(5/2)]',
                transverse_total_eV=float(transverse_total), longitudinal_total_eV=float(parallel_total),
                longitudinal_power_watts=float(parallel_power), total_eV=float(radiation_total),
                fraction_of_removed_kinetic_energy=float(radiation_total/(-delta_kinetic)),
                scope='Free-space classical Lienard radiation of the smooth circular reference, with longitudinal acceleration pdot/(m*gamma^3). Actual gap/fringe radiation, chamber effects and a self-force trajectory are not computed'),
            finite_turn_step=dict(**step_bounds, theorem='For F(n)=C/(beta*c) or C/(beta*gamma*c), gamma=gamma_i-D*n/(m*c^2), F increases and F_doubleprime>0. Over complete unit-turn intervals, 0<=trapezoid-integral<=Delta F_prime/8. The final fractional turn is handled explicitly. This encloses start-energy and end-energy clock sums without iterating33.8million turns',
                conventions='Start-energy means pre-loss energy is used during the revolution, followed by the loss; end-energy means the loss occurs first and the revolution uses post-loss energy. The final turn has fractional weight r and energy decrement r*D',
                scope='An ideal instantaneous-kick quadrature benchmark, not an error bound on the actual distributed-induction and localized-RF trajectory. The inequalities are analytic; displayed constants are65-digit evaluations rounded to floats, not outward-rounded interval certificates. Conservative absolute bounds are2.756microseconds lab and3.034microseconds proper. The last fractional turn is mathematical and does not demonstrate extraction at that phase'),
            stages=stages, plot=plot,
            phase_and_capture=dict(phase_perturbations_solved=False, transverse_ramp_confinement_solved=False, capture_solved=False,
                effective_phase_program='At the averaged reference, oscillator phase may be2*pi*N(t)-pi/2 with positive effective amplitude112.5 minus the numerical per-turn radiation loss in eV, expressed in volts for this proton; exact arrivals still require a finite-gap/chirp model',
                scope='Decelerating crest has zero linear voltage slope and supplies no claimed phase focusing. Existing fixed-energy RF and static transverse checks do not prove stability along this ramp. Bunch emittance, cooling, finite gap transit, field errors, magnet transitions and final trapping remain unsolved'),
            verification=dict(local_Maxwell_symbolic_residuals=[0]*len(maxwell), electric_divergence_zero_by_axisymmetry=True,
                turn_clock_second_derivative_identities_verified=True,
                independent_energy_integral_lab_error_seconds=float(integral_errors[0]),
                independent_energy_integral_proper_error_seconds=float(integral_errors[1]),
                independent_radiation_integral_error_eV=float(integral_errors[2]),
                initial_radiation_relative_difference_from_saved_ring=source_radiation_error,
                direct_1024_turn_fractional_quadrature_checks=direct_bounds),
            scope='Reference clocks and a signed energy/field/frequency ramp, not an implemented decelerator, continuous finite-source transport, capture demonstration, human system or backward time travel. Both clocks run throughout; the17s difference accrues during this already-stored reference interval',
            sources=[dict(url='https://arxiv.org/pdf/1601.04901',support='Tecker: Eqs1-9 Lorentz energy and induction; Eqs31-38 RF synchronism and fixed-radius magnetic ramp. The induction/RF partition here is derived from our specified field, not copied from a no-induction approximation'),
                dict(url='https://arxiv.org/pdf/1601.05230',support='Jensen: section3.6, Eqs21-22 define effective accelerating voltage and its finite-transit reduction; hardware voltage is not fixed by energy per turn alone'),
                dict(url='https://arxiv.org/pdf/2107.09131',support='Prat: classical free-space synchrotron radiation power and bending-magnet scope')],
            derivation_attribution='The clocks, explicit n=.4 induction partition, radiation integral and finite-turn quadrature enclosures are calculations for this ideal ramp. No experimental machine specification is inferred')

    def rhs(s, y):
        gamma = y[0]
        u = math.sqrt(gamma*gamma-1)
        return (-updiff*u/gamma, lab_float/gamma)
    integration = solve_ivp(rhs, (0., 1.), (10., 0.), method='DOP853',
                            rtol=2e-13, atol=(2e-14, 2e-14), max_step=.005)
    assert integration.success
    ode_errors = (abs(integration.y[0, -1]-gammaf_float),
                  abs(integration.y[1, -1]-proper_float))
    assert ode_errors[0] < 2e-12 and ode_errors[1] < 2e-11
    result['verification'].update(independent_lab_time_ODE_gamma_error=ode_errors[0],
        independent_lab_time_ODE_proper_error_seconds=ode_errors[1],
        independent_lab_time_ODE_evaluations=int(integration.nfev))
    assert all(np.isfinite(result['plot'][key]).all() for key in result['plot'])
    assert np.all(np.diff(result['plot']['lab_seconds']) > 0)
    assert np.all(np.diff(result['plot']['field_tesla']) < 0)
    assert np.all(np.diff(result['plot']['rf_frequency_hz']) < 0)
    return result


def _rf_ramp_spread_model(ring, ramp):
    """Canonical guiding-circle model and a separate linear betatron norm bound."""
    import math
    import numpy as np
    constants = ring['constants']
    c, charge, mass, eps = [constants[key] for key in
        ('c', 'elementary_charge', 'proton_mass_kg', 'vacuum_permittivity')]
    radius = ramp['parameters']['bending_radius_metres']
    circumference = ramp['parameters']['circumference_metres']
    rest = ramp['parameters']['rest_energy_eV']
    decrement = ramp['parameters']['net_energy_removed_per_turn_eV']
    ui = math.sqrt(99.)
    gf = ramp['reference']['final_gamma']; uf = math.sqrt(gf*gf-1)
    k = c*decrement/(circumference*rest)
    field0 = ramp['induction']['field_initial_tesla']
    duration = ramp['reference']['lab_elapsed_seconds']
    radiation_coefficient = charge/(3*eps*radius)
    parallel_power = ramp['radiation']['longitudinal_power_watts']
    canonical_scale = 1e-6
    aperture = .01

    def reference(t):
        u = ui-k*np.asarray(t)
        gamma = np.sqrt(1+u*u)
        field = field0*u/ui
        frequency = c*u/(circumference*gamma)
        radiation = radiation_coefficient*gamma*u**3+parallel_power/(charge*frequency)
        return u, gamma, field, frequency, radiation

    def circular(t, canonical_error_scaled):
        us, gs, field, fs, radiation_s = reference(t)
        delta_ell = canonical_scale*np.asarray(canonical_error_scaled)*field0/field
        # Invert ell(r)=-.6r^2+.15r^4 on the branch near r=1.
        q = delta_ell/.15
        discriminant = 1+q
        assert np.all(discriminant > 0), 'guiding-circle inversion left its local branch'
        r2minus1 = -q/(1+np.sqrt(discriminant))
        rminus1 = r2minus1/(1+np.sqrt(1+r2minus1))
        rho = 1+rminus1
        fshape = 1-.2*r2minus1
        u = us*rho*fshape
        gamma = np.sqrt(1+u*u)
        # Rationalized frequency difference avoids cancellation at the reference.
        slip = r2minus1*(.2*(1+fshape)+us*us*fshape*fshape)/(gamma*(gamma+fshape*gs))
        radiation = radiation_coefficient*gamma*u**3/rho+parallel_power/(charge*(fs*(1-slip)))
        return rho, rminus1, u, gamma, 2*math.pi*fs*slip, radiation-radiation_s

    def rhs(phase_s, t, state):
        count = len(state)//2
        psi, canonical = state[:count], state[count:]
        _, _, _, _, radiation_s = reference(t)
        peak_voltage = (112.5-radiation_s)/(-math.sin(phase_s))
        rho, dr, u, gamma, phase_rate, radiation_difference = circular(t, canonical)
        rf_difference = 2*peak_voltage*np.cos(phase_s+psi/2)*np.sin(psi/2)
        canonical_rate = -(rf_difference-radiation_difference)/(2*math.pi*field0*radius**2*canonical_scale)
        return np.concatenate((phase_rate, canonical_rate))

    def initial_state(phase, relative_momentum_error):
        # p/p_s=r*f(r). Solve only this monotone near-reference branch.
        from scipy.optimize import brentq
        rho = brentq(lambda r:r*(1.2-.2*r*r)-(1+relative_momentum_error), .99, 1.01,
                     xtol=5e-15, rtol=1e-14)
        delta_ell = (rho-1)*(rho+1)*(-.45+.15*rho*rho)
        return float(phase), delta_ell/canonical_scale

    def transverse(t, direction):
        tune = math.sqrt(.6 if direction == 'radial' else .4)
        us, gamma, _, _, _ = reference(t)
        omega0 = tune*c*ui/(radius*10.)
        x0, slope0 = .0001, 1e-6
        xdot0 = c*ui/10*slope0
        initial_norm = math.hypot(x0, (xdot0-k*x0/(2*ui))/omega0)
        def primitive(u):
            return -3*np.sqrt(1+u*u)/(2*u*u)+.5*np.arcsinh(1/u)
        exponent = k*radius/(8*tune*c)*(primitive(ui)-primitive(us))
        return initial_norm*np.sqrt(ui/us)*np.exp(exponent), exponent

    return dict(c=c, charge=charge, mass=mass, eps=eps, radius=radius,
        circumference=circumference, rest=rest, decrement=decrement, ui=ui, uf=uf, k=k,
        field0=field0, duration=duration, canonical_scale=canonical_scale, aperture=aperture,
        reference=reference, circular=circular, rhs=rhs, initial_state=initial_state, transverse=transverse)


def _rf_ramp_spread_trajectory(ring, ramp, crest=False, rtol=3e-10):
    """Resolve synchrotron motion over the real ramp; never speed up the ramp."""
    import math
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    model = _rf_ramp_spread_model(ring, ramp)
    phase_s = -math.pi/2 if crest else 7*math.pi/6
    seeds = [(.02, 1e-7)] if crest else [(p, d) for p in (-.02, .02) for d in (-1e-7, 1e-7)]
    initial = [model['initial_state'](*seed) for seed in seeds]
    state0 = np.array([pair[0] for pair in initial]+[pair[1] for pair in initial])
    count = len(seeds)

    def guarantee_inside(t, state):
        offset = model['radius']*model['circular'](t, state[count:])[1]
        envelope = model['transverse'](t, 'radial')[0]
        return model['aperture']-float(np.max(np.abs(offset)))-float(envelope)

    def guarantee_outside(t, state):
        offset = model['radius']*model['circular'](t, state[count:])[1]
        envelope = model['transverse'](t, 'radial')[0]
        return model['aperture']-float(np.max(np.abs(offset)))+float(envelope)

    guarantee_inside.direction = guarantee_outside.direction = -1
    guarantee_inside.terminal = not crest
    guarantee_outside.terminal = crest
    started = time.perf_counter()
    integration = solve_ivp(lambda t, y:model['rhs'](phase_s, t, y), (0., model['duration']), state0,
        method='DOP853', rtol=rtol, atol=np.array([rtol*.01]*count+[rtol*.1]*count),
        max_step=.0005, dense_output=True, events=(guarantee_inside, guarantee_outside))
    assert integration.success
    elapsed = time.perf_counter()-started
    times = integration.t
    values = integration.y
    rho, offset_fraction, momentum, gamma, _, _ = model['circular'](times, values[count:])
    offset = model['radius']*offset_fraction
    radial_envelope = model['transverse'](times, 'radial')[0]
    vertical_envelope = model['transverse'](times, 'vertical')[0]
    inside_times, outside_times = [events.tolist() for events in integration.t_events]
    outcome = ('guaranteed_outside_in_reduced_model' if outside_times else
               'inside_guarantee_ended' if inside_times else 'four_sample_envelopes_within_aperture')
    sample_times = np.unique(np.r_[np.linspace(0, times[-1], 1201), inside_times, outside_times])
    sample_values = integration.sol(sample_times)
    sampled = model['circular'](sample_times, sample_values[count:])
    us, gs, _, fs, _ = model['reference'](times)
    energy_error = model['rest']*(gamma-gs)
    per_case = []
    for j, seed in enumerate(seeds):
        margin = model['aperture']-np.abs(offset[j])-radial_envelope
        per_case.append(dict(initial_phase_error_radians=seed[0], initial_relative_momentum_error=seed[1],
            maximum_sampled_absolute_guiding_offset_metres=float(np.max(np.abs(offset[j]))),
            minimum_sampled_radial_envelope_margin_metres=float(np.min(margin)),
            maximum_sampled_absolute_phase_error_radians=float(np.max(np.abs(values[j]))),
            maximum_sampled_absolute_energy_error_eV=float(np.max(np.abs(energy_error[j]))),
            final_phase_error_radians=float(values[j, -1]), final_canonical_error_scaled=float(values[count+j, -1]),
            final_guiding_offset_metres=float(offset[j, -1]), final_energy_error_eV=float(energy_error[j, -1])))
    record = dict(crest_control=crest, phase_s_radians=phase_s, rtol=rtol, outcome=outcome,
        final_lab_seconds=float(times[-1]), solver_evaluations=int(integration.nfev),
        accepted_steps=len(times)-1, execution_seconds=elapsed,
        inside_guarantee_end_times_seconds=inside_times, guaranteed_outside_times_seconds=outside_times,
        minimum_vertical_envelope_margin_metres=float(np.min(model['aperture']-vertical_envelope)),
        cases=per_case,
        plot=dict(lab_seconds=sample_times.tolist(), phase_error_radians=sample_values[:count].tolist(),
            guiding_offset_metres=(model['radius']*sampled[1]).tolist(),
            radial_betatron_envelope_metres=np.asarray(model['transverse'](sample_times, 'radial')[0]).tolist(),
            vertical_betatron_envelope_metres=np.asarray(model['transverse'](sample_times, 'vertical')[0]).tolist()),
        scope='Finite sampled corners in the nonlinear turn-averaged guiding-circle model plus a separate nominal linear betatron envelope. Numerical extrema are sampled at adaptive solver nodes. Envelope overlap alone is not a proven wall hit; the crest outside bound implies departure within this reduced model. No whole-box, full-gap or complete-ramp Lorentz acceptance theorem')
    if not crest and rtol <= 3e-12 and times[-1] == model['duration']:
        checkpoints = []
        for nominal_turn in (0, ramp['reference']['complete_turns']-20):
            g = 10-model['decrement']*nominal_turn/model['rest']
            u = math.sqrt(g*g-1)
            time_at_turn = (model['ui']-u)/model['k']
            point = integration.sol(time_at_turn)
            checkpoints.append(dict(nominal_turn=nominal_turn, lab_seconds=time_at_turn,
                representative_case_index=3, phase_error_radians=float(point[3]),
                canonical_error_scaled=float(point[count+3])))
        record['lorentz_window_checkpoints'] = checkpoints
    return record, integration, model


def _rf_ramp_lorentz_window(ring, ramp, checkpoint, fine=False):
    """True-ramp guide Lorentz motion between explicitly ideal thin RF gaps."""
    import math
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    model = _rf_ramp_spread_model(ring, ramp)
    c, radius, rest = model['c'], model['radius'], model['rest']
    nominal_turn = checkpoint['nominal_turn']
    tstart = checkpoint['lab_seconds']
    psi0, canonical0 = checkpoint['phase_error_radians'], checkpoint['canonical_error_scaled']
    phase_s = 7*math.pi/6
    us0, gs0, field_start, _, _ = model['reference'](tstart)
    rho, _, u0, gamma0, _, _ = model['circular'](tstart, canonical0)
    phase_rate, canonical_rate = model['rhs'](phase_s, tstart, np.array((psi0, canonical0)))
    bdot = ramp['induction']['field_slope_tesla_per_second']
    delta_ell = model['canonical_scale']*canonical0*model['field0']/field_start
    ell_rate = model['canonical_scale']*canonical_rate*model['field0']/field_start-bdot*delta_ell/field_start
    guiding_radial_velocity = radius*ell_rate/(-1.2*rho+.6*rho**3)

    def betatron_seed(direction):
        tune = math.sqrt(.6 if direction == 'radial' else .4)
        omega_i = tune*c*model['ui']/(10*radius)
        omega = tune*c*us0/(gs0*radius)
        x0, xdot0 = .0001, c*model['ui']/10*1e-6
        a, b = x0, (xdot0-model['k']*x0/(2*model['ui']))/omega_i
        phase = 2*math.pi*tune*nominal_turn
        scale = math.sqrt(model['ui']/us0)
        value = scale*(a*math.cos(phase)+b*math.sin(phase))
        velocity = model['k']/(2*us0)*value+scale*omega*(-a*math.sin(phase)+b*math.cos(phase))
        return value, velocity

    xbeta, vx_beta = betatron_seed('radial')
    zbeta, vz_beta = betatron_seed('vertical')
    pr0, pz0 = gamma0*(guiding_radial_velocity+vx_beta)/c, gamma0*vz_beta/c
    pp0 = -math.sqrt(u0*u0-pr0*pr0-pz0*pz0)
    # State: cylindrical radius,z,p_r/(mc),p_phi/(mc),p_z/(mc),delta_t,delta_tau,W_ind/eV.
    initial = np.array((radius*rho+xbeta, zbeta, pr0, pp0, pz0, 0., 0., 0.))
    initial_gamma = math.sqrt(1+pr0*pr0+pp0*pp0+pz0*pz0)
    assert abs(initial_gamma-gamma0) < 3e-15
    force_factor = model['charge']/(model['mass']*c)

    def local_reference(delta_t):
        u = us0-model['k']*delta_t
        gamma = math.sqrt(1+u*u)
        field = field_start+bdot*delta_t
        # The start index is an integer; use a rationalized LOCAL turn advance.
        advance = model['rest']/model['decrement']*model['k']*delta_t*(us0+u)/(gs0+gamma)
        frequency = c*u/(model['circumference']*gamma)
        rad = model['charge']/(3*model['eps']*radius)*gamma*u**3
        rad += ramp['radiation']['longitudinal_power_watts']/(model['charge']*frequency)
        return field, advance, rad

    def rhs(theta, y):
        rr, zz, pr, pp, pz, dt, tau, work = y
        gamma = math.sqrt(1+pr*pr+pp*pp+pz*pz)
        vr, vp, vz = c*pr/gamma, c*pp/gamma, c*pz/gamma
        field = field_start+bdot*dt
        rh, zh = rr/radius, zz/radius
        br = -.4*field*rh*zh
        bz = field*(1.2-.2*rh*rh+.4*zh*zh)
        ap = .6*rh-.05*rh**3+.2*rh*zh*zh
        ep = -bdot*radius*ap
        assert pp < 0 and rr > 0
        time_factor = -rr/vp
        derivative = np.array((vr, vz,
            pp*vp/rr+force_factor*vp*bz,
            -pr*vp/rr+force_factor*(ep+vz*br-vr*bz),
            -force_factor*vp*br, 1., 1/gamma, ep*vp))
        return derivative*time_factor

    turns = 16
    theta0, theta_end = -psi0, -psi0+2*math.pi*turns
    gap_angles = [2*math.pi*j for j in range(math.floor(theta0/(2*math.pi))+1,
        math.floor(theta_end/(2*math.pi))+1)]
    boundaries = gap_angles+[theta_end]
    sample_angles = np.unique(np.r_[np.linspace(theta0, theta_end, 1025), gap_angles])
    all_angles, all_states, all_rf_work = [theta0], [initial.copy()], [0.]
    plot_angles, plot_states, plot_rf = [theta0], [initial.copy()], [0.]
    start, state, rf_work = theta0, initial, 0.
    kicks = []
    started = time.perf_counter()
    evaluations = 0
    for stop in boundaries:
        solve = solve_ivp(rhs, (start, stop), state, method='DOP853',
            rtol=3e-12 if fine else 2e-10,
            atol=(1e-11, 1e-11, 2e-14, 2e-14, 2e-14, 1e-18, 1e-18, 1e-9),
            max_step=.025 if fine else .05, dense_output=True)
        assert solve.success
        evaluations += solve.nfev
        all_angles.extend(solve.t[1:].tolist())
        all_states.extend(solve.y[:, 1:].T.copy())
        all_rf_work.extend([rf_work]*(len(solve.t)-1))
        here = sample_angles[(sample_angles > start)&(sample_angles <= stop)]
        plot_angles.extend(here.tolist())
        plot_states.extend(solve.sol(here).T.copy())
        plot_rf.extend([rf_work]*len(here))
        state = solve.y[:, -1].copy()
        if stop in gap_angles:
            _, advance, rad = local_reference(state[5])
            gain = 2*(112.5-rad)*math.sin(2*math.pi*advance+phase_s)
            gamma_before = math.sqrt(1+sum(state[2:5]**2))
            change = gain/rest
            pp_square = state[3]**2+2*gamma_before*change+change*change
            assert pp_square > 0
            state[3] = -math.sqrt(pp_square)
            rf_work += gain
            kicks.append(dict(lab_seconds=tstart+state[5], work_eV=gain,
                rf_phase_radians=2*math.pi*advance+phase_s))
            all_angles.append(stop); all_states.append(state.copy()); all_rf_work.append(rf_work)
            # At a gap, retain the post-kick momentum at the same sampled position/time.
            if plot_angles[-1] == stop:
                plot_states[-1], plot_rf[-1] = state.copy(), rf_work
        start = stop
    elapsed = time.perf_counter()-started
    all_states = np.array(all_states).T
    plot_states = np.array(plot_states).T
    gammas = np.sqrt(1+np.sum(all_states[2:5]**2, axis=0))
    work_residual = rest*(gammas-initial_gamma)-all_states[7]-np.asarray(all_rf_work)
    assert np.max(np.abs(work_residual)) < .01
    assert tstart+state[5] < model['duration']
    comparison = solve_ivp(lambda t, y:model['rhs'](phase_s, t, y),
        (tstart, tstart+state[5]), (psi0, canonical0), method='DOP853',
        rtol=2e-12, atol=1e-13, dense_output=True)
    assert comparison.success
    average = comparison.sol(tstart+plot_states[5])
    average_rho = model['circular'](tstart+plot_states[5], average[1])[0]
    field = field_start+bdot*plot_states[5]
    rh, zh = plot_states[0]/radius, plot_states[1]/radius
    ap = .6*rh-.05*rh**3+.2*rh*zh*zh
    physical_canonical_error = (rh*plot_states[3]/model['ui']+field/model['field0']*(rh*ap+.45))/model['canonical_scale']
    projected_rho = model['circular'](tstart+plot_states[5], physical_canonical_error)[0]
    phase_actual = np.array([2*math.pi*local_reference(t)[1] for t in plot_states[5]])-np.asarray(plot_angles)
    return dict(nominal_start_turn=nominal_turn, start_lab_seconds=tstart, turns=turns,
        representative_case_index=3, fine=fine, execution_seconds=elapsed, solver_evaluations=evaluations,
        initial_total_gamma=initial_gamma, elapsed_lab_seconds=float(state[5]),
        elapsed_proper_seconds=float(state[6]), induction_work_eV=float(state[7]), rf_work_eV=float(rf_work),
        final_total_gamma=float(gammas[-1]), maximum_work_residual_eV=float(np.max(np.abs(work_residual))),
        minimum_sampled_radial_aperture_margin_metres=float(np.min(model['aperture']-np.abs(all_states[0]-radius))),
        minimum_sampled_vertical_aperture_margin_metres=float(np.min(model['aperture']-np.abs(all_states[1]))),
        maximum_projected_guiding_difference_from_average_metres=float(np.max(np.abs(radius*(projected_rho-average_rho)))),
        maximum_phase_difference_from_average_radians=float(np.max(np.abs(phase_actual-average[0]))),
        gap_count=len(kicks), kicks=kicks,
        plot=dict(angle_radians=np.asarray(plot_angles).tolist(), local_lab_seconds=plot_states[5].tolist(),
            radial_offset_metres=(plot_states[0]-radius).tolist(), vertical_offset_metres=plot_states[1].tolist(),
            projected_guiding_offset_metres=(radius*(projected_rho-1)).tolist(),
            averaged_guiding_offset_metres=(radius*(average_rho-1)).tolist(), phase_error_radians=phase_actual.tolist()),
        scope='Sixteen real revolutions at the unchanged physical ramp rate. Guide-field segments solve relativistic Lorentz motion with the exact local induction E/B. RF is an ideal zero-length longitudinal energy kick preserving transverse mechanical momentum; no finite Maxwell cavity, RF fringe focusing or gap-radiation calculation is provided. No radiation reaction is included in these windows. The initialization uses the separate adiabatic transverse approximation, and these two local comparisons do not validate full-ramp Lorentz confinement')


def proton_rf_ramp_spread_checks(ring, ramp, records=None):
    """Four finite spreads in an averaged ramp, with explicitly local gap checks."""
    import copy
    import math
    import numpy as np
    import sympy as sp
    model = _rf_ramp_spread_model(ring, ramp)
    if records is None:
        main = _rf_ramp_spread_trajectory(ring, ramp)[0]
        fine = _rf_ramp_spread_trajectory(ring, ramp, rtol=3e-12)[0]
        crest = _rf_ramp_spread_trajectory(ring, ramp, crest=True)[0]
        crest_fine = _rf_ramp_spread_trajectory(ring, ramp, crest=True, rtol=3e-12)[0]
        windows = [_rf_ramp_lorentz_window(ring, ramp, point, fine=resolution)
            for point in fine['lorentz_window_checkpoints'] for resolution in (False, True)]
    else:
        main, fine, crest, crest_fine, windows = [copy.deepcopy(records[key])
            for key in ('main', 'fine', 'crest', 'crest_fine', 'windows')]
    assert len(main['cases']) == len(fine['cases']) == 4 and len(windows) == 4
    expected_seeds = [(p, d) for p in (-.02, .02) for d in (-1e-7, 1e-7)]
    for trace, tolerance, is_crest in ((main, 3e-10, False), (fine, 3e-12, False),
        (crest, 3e-10, True), (crest_fine, 3e-12, True)):
        assert trace['rtol'] == tolerance and trace['crest_control'] == is_crest
        expected = [(.02, 1e-7)] if is_crest else expected_seeds
        assert [(case['initial_phase_error_radians'], case['initial_relative_momentum_error'])
                for case in trace['cases']] == expected
        assert trace['phase_s_radians'] == (-math.pi/2 if is_crest else 7*math.pi/6)
    assert main['outcome'] == fine['outcome'] == 'four_sample_envelopes_within_aperture'
    assert crest['outcome'] == crest_fine['outcome'] == 'guaranteed_outside_in_reduced_model'
    assert main['final_lab_seconds'] == fine['final_lab_seconds'] == model['duration']

    # Check the energy identity underlying the canonical torque equation.
    r, u, k = sp.symbols('r u k', positive=True)
    fshape = sp.Rational(6, 5)-r*r/5
    ap = 3*r/5-r**3/20
    ell = r*ap-r*r*fshape
    gshape = r*fshape
    assert sp.simplify(sp.diff(ell, r)+r*sp.diff(gshape, r)) == 0
    assert sp.simplify(gshape-sp.diff(gshape, r)*ell/sp.diff(ell, r)-ap) == 0
    assert sp.simplify(-sp.diff(gshape, r)/sp.diff(ell, r)-1/r) == 0
    numerator = 3*k*k/(4*u*u)-k*k/(2*(1+u*u))
    assert sp.simplify(numerator-k*k*(u*u+3)/(4*u*u*(1+u*u))) == 0
    primitive = -3*sp.sqrt(1+u*u)/(2*u*u)+sp.asinh(1/u)/2
    assert sp.simplify(sp.diff(primitive, u)-(u*u+3)/(u**3*sp.sqrt(1+u*u))) == 0
    probe_errors = []
    for fraction in (0., .5, .999):
        time = fraction*model['duration']
        us, gs, field, fs, rad_s = model['reference'](time)
        for psi in (-.02, 0., .02):
            for canonical in (-1., 0., 1.):
                rho, _, momentum, gamma, _, rad_difference = model['circular'](time, canonical)
                derivative = model['rhs'](7*math.pi/6, time, np.array((psi, canonical)))
                delta_ell = model['canonical_scale']*canonical*model['field0']/field
                bdot = ramp['induction']['field_slope_tesla_per_second']
                delta_ell_dot = model['canonical_scale']*derivative[1]*model['field0']/field-bdot*delta_ell/field
                rdot = delta_ell_dot/(-1.2*rho+.6*rho**3)
                udot = -model['k']*rho*(1.2-.2*rho*rho)+us*(1.2-.6*rho*rho)*rdot
                energy_rate = model['rest']*momentum/gamma*udot
                velocity = model['c']*momentum/gamma
                actual_frequency = velocity/(model['circumference']*rho)
                induction_rate = bdot*model['radius']*(.6*rho-.05*rho**3)*velocity
                rf_work = 2*(112.5-rad_s)*math.sin(7*math.pi/6+psi)
                field_rate = induction_rate+actual_frequency*(rf_work-rad_s-rad_difference)
                probe_errors.append(abs(energy_rate-field_rate)/max(1., abs(energy_rate)))
    assert max(probe_errors) < 2e-13
    assert np.array_equal(model['rhs'](7*math.pi/6, 0., np.zeros(2)), np.zeros(2))

    # Full-interval bound for the separate nominal linear transverse equation.
    transverse = {}
    with mp.workdps(65):
        km, rm, cm, ui, uf = [mp.mpf(str(model[key])) for key in ('k', 'radius', 'c', 'ui', 'uf')]
        integral = mp.quad(lambda v:(v*v+3)/(v**3*mp.sqrt(1+v*v)), [uf, 1, ui])
        primitive_m = lambda v:-3*mp.sqrt(1+v*v)/(2*v*v)+mp.asinh(1/v)/2
        assert abs(integral-(primitive_m(ui)-primitive_m(uf))) < mp.mpf('1e-55')
        for direction, tune_squared in (('radial', mp.mpf('.6')), ('vertical', mp.mpf('.4'))):
            tune = mp.sqrt(tune_squared)
            exponent = km*rm*integral/(8*tune*cm)
            envelope, production_exponent = model['transverse'](model['duration'], direction)
            assert abs(float(exponent)-production_exponent) < 1e-17
            transverse[direction] = dict(tune=float(tune), final_amplitude_bound_metres=float(envelope),
                full_interval_norm_exponent=float(exponent), excess_norm_factor=float(mp.expm1(exponent)),
                maximum_reference_frequency_adiabatic_parameter=float(km*rm/(tune*cm*uf*uf*mp.sqrt(1+uf*uf))))

        rest, charge, mass, eps = [mp.mpf(str(model[key])) for key in ('rest', 'charge', 'mass', 'eps')]
        circumference, decrement = mp.mpf(str(model['circumference'])), mp.mpf(250)
        alpha = mp.mpf(5)/3
        gf = mp.mpf(str(ramp['reference']['final_gamma']))
        peak_gamma = mp.sqrt(3/alpha)
        def frozen(gamma):
            um = mp.sqrt(gamma*gamma-1)
            freq = cm*um/(circumference*gamma)
            rad = charge*gamma*um**3/(3*eps*rm)
            rad += mp.mpf(str(ramp['radiation']['longitudinal_power_watts']))/(charge*freq)
            voltage = 2*(mp.mpf('112.5')-rad)
            eta = alpha-1/(gamma*gamma)
            omega = mp.sqrt(eta/gamma*charge*voltage*mp.sqrt(3)/2/(2*mp.pi*mass*rm*rm))
            return freq, omega
        def frozen_map_difference(gamma):
            freq, omega = frozen(gamma)
            mu = omega/freq
            return 2*mp.asin(mu/2)-mu
        accumulated_map_difference = rest/decrement*mp.quad(frozen_map_difference,
            [gf, mp.mpf('1.1'), peak_gamma, 2, 5, 10])
        frozen_rows = []
        for gamma in (mp.mpf(10), peak_gamma, gf):
            freq, omega = frozen(gamma)
            frozen_rows.append(dict(gamma=float(gamma), kinetic_energy_MeV=float((gamma-1)*rest/1e6),
                kinetic_curvature_factor=float((alpha-1/gamma**2)/gamma),
                synchrotron_frequency_hz=float(omega/(2*mp.pi)), linear_phase_advance_per_turn=float(omega/freq)))

    assert main['plot']['lab_seconds'] == fine['plot']['lab_seconds']
    comparison = dict(
        maximum_common_display_sample_phase_difference_radians=float(np.max(np.abs(
            np.array(main['plot']['phase_error_radians'])-np.array(fine['plot']['phase_error_radians'])))),
        maximum_common_display_sample_guiding_difference_metres=float(np.max(np.abs(
            np.array(main['plot']['guiding_offset_metres'])-np.array(fine['plot']['guiding_offset_metres'])))),
        maximum_final_phase_difference_radians=max(abs(a['final_phase_error_radians']-b['final_phase_error_radians'])
            for a, b in zip(main['cases'], fine['cases'])),
        maximum_final_guiding_difference_metres=max(abs(a['final_guiding_offset_metres']-b['final_guiding_offset_metres'])
            for a, b in zip(main['cases'], fine['cases'])),
        maximum_adaptive_node_peak_estimate_change_metres=max(abs(a['maximum_sampled_absolute_guiding_offset_metres']-
            b['maximum_sampled_absolute_guiding_offset_metres']) for a, b in zip(main['cases'], fine['cases'])),
        crest_inside_time_difference_seconds=abs(crest['inside_guarantee_end_times_seconds'][0]-crest_fine['inside_guarantee_end_times_seconds'][0]),
        crest_outside_time_difference_seconds=abs(crest['guaranteed_outside_times_seconds'][0]-crest_fine['guaranteed_outside_times_seconds'][0]))
    assert comparison['maximum_common_display_sample_phase_difference_radians'] < 2e-8
    assert comparison['maximum_common_display_sample_guiding_difference_metres'] < 4e-9
    assert comparison['crest_outside_time_difference_seconds'] < 1e-10
    window_comparisons = []
    for coarse, resolved in zip(windows[::2], windows[1::2]):
        assert coarse['fine'] is False and resolved['fine'] is True
        assert coarse['nominal_start_turn'] == resolved['nominal_start_turn']
        assert coarse['gap_count'] == resolved['gap_count'] == 16
        assert coarse['plot']['angle_radians'] == resolved['plot']['angle_radians']
        position_difference = max(float(np.max(np.abs(np.array(coarse['plot'][name])-np.array(resolved['plot'][name]))))
            for name in ('radial_offset_metres', 'vertical_offset_metres'))
        assert position_difference < 2e-9
        window_comparisons.append(dict(nominal_start_turn=coarse['nominal_start_turn'],
            maximum_common_angle_position_difference_metres=position_difference,
            elapsed_lab_difference_seconds=abs(coarse['elapsed_lab_seconds']-resolved['elapsed_lab_seconds']),
            elapsed_proper_difference_seconds=abs(coarse['elapsed_proper_seconds']-resolved['elapsed_proper_seconds']),
            guiding_model_difference_metres=resolved['maximum_projected_guiding_difference_from_average_metres']))

    return dict(
        status='verified bounded reduced-model study; local Lorentz/thin-gap checks completed',
        parameters=dict(net_energy_removed_per_turn_eV=250., actual_ramp_lab_seconds=model['duration'],
            phase_s_radians=7*math.pi/6, phase_s_degrees=210., effective_peak_voltage_formula='Vhat(t)=2*(112.5-Urad_s[eV]) volts, with signed RF work Vhat*sin(phi_s+psi) eV',
            initial_phase_errors_radians=[-.02, .02], initial_relative_momentum_errors=[-1e-7, 1e-7],
            transverse_initial_position_metres=dict(radial=.0001, vertical=.0001),
            transverse_initial_slopes=dict(radial=1e-6, vertical=1e-6),
            radial_aperture_half_width_metres=.01, vertical_aperture_half_width_metres=.01,
            fixed_lorentz_representative_initial_phase=.02, fixed_lorentz_representative_initial_relative_momentum=1e-7),
        model=dict(canonical='L=e*B0*R^2*ell(r), ell=-.6*r^2+.15*r^4; deltaL=L-L_s and saved canonical error is deltaL/(e*B_i*R^2*1e-6)',
            equations='psi_dot=omega_s-omega(r,t); deltaL_dot=-e*[DeltaW_RF-DeltaUrad]/(2*pi), with W in eV. p/(mc)=u_s*r*(1.2-.2*r^2). The actual linear B0(t) is retained',
            induction_identity='ell_prime=-r*(r*f)_prime makes dE/dt=e*R*Bdot*a_phi*v+f_rev*(W_RF-Urad)*e exactly in the guiding-circle model. The induction contribution is not dropped or counted twice',
            phase_convention='psi=theta_s-theta for clockwise orbital angle theta; positive psi means later RF-gap arrival. phi_s=210degrees has negative voltage slope and restores small longitudinal errors',
            averaging_scope='The RF torque is averaged over turns and the guiding state is assigned to an instantaneous circular orbit. The real24.111s ramp is not accelerated. Radiation uses leading transverse free-space loss plus the negligible nominal longitudinal term; finite-gap effects and guiding/betatron coupling are not closed'),
        main=main, fine=fine, crest=crest, crest_fine=crest_fine, windows=windows,
        transverse_linear_bound=dict(**transverse,
            equation='xddot+(gammadot/gamma)*xdot+omega_beta^2*x=0, omega_beta=nu*c*u/(R*gamma), u=u_i-k*t',
            transformation='x=S*y, S=sqrt(u_i/u), theta_dot=omega_beta. y_theta_theta+(1+delta)*y=0; delta=[k^2*(u^2+3)/(4*u^2*gamma^2)]/omega_beta^2',
            bound='sqrt(y^2+y_theta^2) grows by at most exp(kappa), kappa=k*R/(8*nu*c)*integral (u^2+3)/(u^3*sqrt(1+u^2))du fromu(t) tou_i. Initial norm includes the Sdot*x0 correction',
            scope='A full-interval norm bound for the stated decoupled nominal linear oscillator, not a bound on omitted RF fringe forces, transverse-longitudinal coupling or full nonlinear Lorentz trajectories'),
        frozen_longitudinal_diagnostics=dict(rows=frozen_rows,
            curvature_maximum_gamma=float(peak_gamma),
            accumulated_frozen_map_vs_ODE_phase_difference_radians=float(accumulated_map_difference),
            definition='Kinetic curvature is proportional to eta/gamma, eta=5/3-1/gamma^2; its maximum is atgamma=sqrt(1.8). The map diagnostic integrates[2*asin(mu/2)-mu]dN, mu=omega_syn/f_rev',
            scope='Frozen linear focusing and discrete-kick phase diagnostics only. The accumulated phase difference is not a global error bound for the changing nonlinear ramp; endpoint-only curvature tests would miss the intermediate maximum'),
        verification=dict(symbolic_canonical_energy_identity=True, symbolic_transverse_residual_and_primitive=True,
            maximum_pointwise_relative_energy_rate_error=max(probe_errors),
            actual_reference_is_exact_equilibrium_of_error_ODE=True, tolerance_comparison=comparison,
            lorentz_window_comparisons=window_comparisons,
            maximum_local_Lorentz_work_residual_eV=max(row['maximum_work_residual_eV'] for row in windows)),
        display_scope='The full-ramp plot contains1201 common display samples and does not resolve every synchrotron oscillation. Peak estimates use the much denser adaptive solver nodes. They are sampled extrema, not certified suprema; differing nodes explain micrometre changes in peak estimates despite nanometre common-time convergence',
        outcome_scope='All four declared off-crest sample envelopes stay inside the1cm radial/vertical aperture in the reduced model. This is not a whole-distribution acceptance theorem. The crest sample loses its inside guarantee and subsequently has |guiding offset|-envelope>1cm, forcing departure within that same reduced model. Two16-turn Lorentz/thin-gap windows remain inside locally but show a late~0.575mm discrepancy from averaging. No full-ramp finite-cavity confinement, injection, extraction, cooling, trapping or physical capture is established',
        sources=[dict(url='https://arxiv.org/pdf/1601.04901',support='Tecker: longitudinal synchronism, dispersion and phase stability. The positive time-slip convention and explicit induction/canonical equations here are derived for our chosen field'),
            dict(url='https://arxiv.org/pdf/1601.05230',support='Jensen: effective RF voltage requires a field and transit-time model; an energy kick alone is not a finite cavity design')],
        derivation_attribution='Canonical nonlinear averaging, the actual-ramp transverse norm bound and the finite-window comparison are calculations for this specified ideal field and source approximation, independently audited at the equation level')


def _rf_event_native_code(ring, ramp):
    """Small streaming C implementation of the reduced, actual-turn RF map."""
    model = _rf_ramp_spread_model(ring, ramp)
    constants = dict(CLIGHT=model['c'], RADIUS=model['radius'], CIRC=model['circumference'],
        REST=model['rest'], D=model['decrement'], RADCOEF=model['charge']/(3*model['eps']*model['radius']),
        PARALLEL=ramp['radiation']['longitudinal_power_watts']/model['charge'])
    declarations = '\n'.join('static const real '+name+' = '+repr(value)+'L;'
                             for name, value in constants.items())
    return r'''
#define __USE_MINGW_ANSI_STDIO 1
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#ifdef EXTENDED
typedef long double real;
#define SQRT sqrtl
#define SIN sinl
#define ATAN atanl
#define ABS fabsl
#define LOG logl
#define MANTISSA LDBL_MANT_DIG
#else
typedef double real;
#define SQRT sqrt
#define SIN sin
#define ATAN atan
#define ABS fabs
#define LOG log
#define MANTISSA DBL_MANT_DIG
#endif
''' + declarations + r'''
static real pi, ui, initial_v, initial_delta, initial_w, initial_ns;
static real rf_correction=0, rf_compensation=0;
static real max_offset=0, max_mid_offset=0, max_phase=0, max_ledger=0;
static long long kicks=0, evaluations=0, segments=0;
static int method=4, subdivisions=1;
typedef struct {real gs, us, g, u, x, slip, rad, rad_s, freq;} Geometry;

static Geometry geometry(real v, const real *y) {
    Geometry a;
    const real q=y[1]/REST;
    a.gs=(10-(D/REST)*v)-(D/REST)*y[0];
    a.us=SQRT((a.gs-1)*(a.gs+1));
    a.g=a.gs+q;
    a.u=SQRT(a.us*a.us+q*(2*a.gs+q));
    const real dp=q*(2*a.gs+q)/(a.us*(a.u+a.us));
    a.x=dp/(real).6;
    for(int j=0;j<2;j++)
        a.x-=(a.x*((real).6-a.x*((real).6+(real).2*a.x))-dp)
            /((real).6-a.x*((real)1.2+(real).6*a.x));
    const real r=1+a.x, r2m=a.x*(2+a.x), shape=1-(real).2*r2m;
    if(!(a.gs>1 && a.g>1 && ABS(a.x)<(real).01)) {
        fprintf(stderr,"Local guiding-circle branch left: v=%.17g x=%.17g gamma=%.17g\n",
                (double)v,(double)a.x,(double)a.g); exit(3);
    }
    a.slip=r2m*((real).2*(1+shape)+a.us*a.us*shape*shape)/(a.g*(a.g+shape*a.gs));
    const real fs=CLIGHT*a.us/(CIRC*a.gs);
    a.freq=fs*(1-a.slip);
    a.rad=RADCOEF*a.g*a.u*a.u*a.u/r+PARALLEL/a.freq;
    a.rad_s=RADCOEF*a.gs*a.us*a.us*a.us+PARALLEL/fs;
    return a;
}

static void derivative(real v, const real *y, real *dy) {
    Geometry a=geometry(v,y);
    real r2m=a.x*(2+a.x);
    dy[0]=a.slip/(1-a.slip);
    dy[1]=D*(dy[0]+(real).45-(real).5*r2m+(real).05*r2m*r2m)-a.rad;
    dy[2]=-CIRC*(1+a.x)*y[1]/(CLIGHT*REST*a.u*a.gs);
    dy[3]=-D*((real).5*r2m-(real).05*r2m*r2m);
    dy[4]=a.rad;
    evaluations++;
}

static void compensated_add(real *value, real *carry, real increment) {
    real adjusted=increment-*carry, total=*value+adjusted;
    *carry=(total-*value)-adjusted; *value=total;
}

static void advance(real v, const real *input, const real *carry_in,
                    real length, real *output, real *carry_out) {
    real h=length/subdivisions, a[5], b[5], c[5], d[5], temp[5];
    for(int j=0;j<5;j++){output[j]=input[j];carry_out[j]=carry_in[j];}
    for(int n=0;n<subdivisions;n++) {
        real start=v+n*h;
        derivative(start,output,a);
        for(int j=0;j<5;j++)temp[j]=output[j]+h*a[j]/2;
        derivative(start+h/2,temp,b);
        if(method==2) {
            for(int j=0;j<5;j++)compensated_add(output+j,carry_out+j,h*b[j]);
        } else {
            for(int j=0;j<5;j++)temp[j]=output[j]+h*b[j]/2;
            derivative(start+h/2,temp,c);
            for(int j=0;j<5;j++)temp[j]=output[j]+h*c[j];
            derivative(start+h,temp,d);
            for(int j=0;j<5;j++)
                compensated_add(output+j,carry_out+j,h*(a[j]+2*b[j]+2*c[j]+d[j])/6);
        }
    }
}

static real ledger(real v, const real *y) {
    return (y[1]-initial_w)-D*(y[0]-initial_delta)
        -(real)112.5*((v-initial_v)-kicks)-y[3]-rf_correction+y[4];
}

static void track(real v, const real *y, real gain, int at_gap) {
    Geometry a=geometry(v,y);
    real offset=ABS(RADIUS*a.x), phase=ABS(2*pi*y[0]), error=ABS(ledger(v,y));
    if(offset>max_offset)max_offset=offset;
    if(phase>max_phase)max_phase=phase;
    if(error>max_ledger)max_ledger=error;
    if(at_gap) {
        real midpoint[5]; for(int j=0;j<5;j++)midpoint[j]=y[j]; midpoint[1]-=gain/2;
        real mid=ABS(RADIUS*geometry(v,midpoint).x);
        if(mid>max_mid_offset)max_mid_offset=mid;
    }
}

static void emit(int kind, real v, const real *y, real gain) {
    Geometry a=geometry(v,y);
    real ns=v+y[0], refclock=CIRC*REST/(CLIGHT*D)*LOG((10+ui)/(a.gs+a.us));
    real lab=CIRC*REST/(CLIGHT*D)*(ui-a.us);
    real midpoint[5];for(int j=0;j<5;j++)midpoint[j]=y[j];midpoint[1]-=gain/2;
    real mid_offset=RADIUS*geometry(v,midpoint).x;
    printf("{\"kind\":%d,\"orbital_turn\":%.17g,\"nominal_turn\":%.17g,"
           "\"lab_seconds\":%.17g,\"phase_error_radians\":%.17g,\"energy_error_eV\":%.17g,"
           "\"guiding_offset_metres\":%.17g,\"midpoint_energy_error_eV\":%.17g,"
           "\"midpoint_guiding_offset_metres\":%.17g,\"proper_clock_error_seconds\":%.17g,"
           "\"proper_seconds\":%.17g,\"rf_gap_work_eV\":%.17g,\"work_residual_eV\":%.17g}\n",
           kind,(double)v,(double)ns,(double)lab,(double)(2*pi*y[0]),(double)y[1],
           (double)(RADIUS*a.x),(double)midpoint[1],(double)mid_offset,(double)y[2],
           (double)(refclock+y[2]),(double)gain,(double)ledger(v,y));
}

int main(int argc, char **argv) {
    if(argc!=9){fprintf(stderr,"Expected start_N stop_N phase dp method subdivisions stride dense_tail\n");return 2;}
    initial_ns=(real)strtold(argv[1],0); real stop=(real)strtold(argv[2],0);
    real phase=(real)strtold(argv[3],0), dp=(real)strtold(argv[4],0);
    method=atoi(argv[5]);subdivisions=atoi(argv[6]);long long stride=atoll(argv[7]);
    long long dense_tail=atoll(argv[8]);
    if(!(stop>initial_ns && stride>0 && subdivisions>0 && (method==2||method==4)))return 2;
    pi=4*ATAN((real)1);ui=SQRT((real)99);
    initial_delta=phase/(2*pi);initial_v=initial_ns-initial_delta;
    real gs=10-D/REST*initial_ns, us=SQRT((gs-1)*(gs+1)), u=us*(1+dp);
    initial_w=REST*us*us*dp*(2+dp)/(SQRT(1+u*u)+gs);
    real y[5]={initial_delta,initial_w,0,0,0},carry[5]={0,0,0,0,0};
    real v=initial_v; long long next=(long long)floorl((long double)v)+1;
    long long first_next=next;
    emit(0,v,y,0); track(v,y,0,0); clock_t loop_started=clock();
    for(;;) {
        real length=(real)next-v, trial[5], trial_carry[5];
        advance(v,y,carry,length,trial,trial_carry);
        if((real)next+trial[0]>=stop) {
            // Solve the actual fixed source-time endpoint, without an extra kick.
            real fraction=length*(stop-(v+y[0]))/(((real)next+trial[0])-(v+y[0]));
            for(int n=0;n<5;n++) {
                advance(v,y,carry,fraction,trial,trial_carry);
                Geometry a=geometry(v+fraction,trial);
                real residual=(v-stop)+fraction+trial[0];
                fraction-=residual*(1-a.slip);
            }
            advance(v,y,carry,fraction,trial,trial_carry);v+=fraction;
            for(int j=0;j<5;j++){y[j]=trial[j];carry[j]=trial_carry[j];}
            segments++;track(v,y,0,0);emit(2,v,y,0);break;
        }
        v=(real)next;
        for(int j=0;j<5;j++){y[j]=trial[j];carry[j]=trial_carry[j];}
        int diagnostic=((next-first_next)%16==0);
        if(diagnostic)track(v,y,0,0);
        Geometry a=geometry(v,y);
        real gain=2*((real)112.5-a.rad_s)*SIN(7*pi/6+2*pi*y[0]);
        compensated_add(y+1,carry+1,gain);
        compensated_add(&rf_correction,&rf_compensation,gain+(real)112.5);
        kicks++;segments++;if(diagnostic)track(v,y,gain,1);
        if(next==first_next || (next-first_next)%stride==0 || stop-(real)next<dense_tail)
            emit(1,v,y,gain);
        next++;
    }
    printf("{\"summary\":true,\"native_loop_clock_seconds\":%.17g,\"mantissa_bits\":%d,\"gaps\":%lld,\"drift_segments\":%lld,"
           "\"rhs_evaluations\":%lld,\"actual_revolutions\":%.17g,\"source_endpoint_residual_turns\":%.17g,"
           "\"maximum_pre_post_guiding_offset_metres\":%.17g,\"maximum_midpoint_guiding_offset_metres\":%.17g,"
           "\"maximum_phase_error_radians\":%.17g,\"maximum_work_residual_eV\":%.17g,"
           "\"induction_work_eV\":%.17g,\"rf_work_eV\":%.17g,\"radiated_energy_eV\":%.17g}\n",
           (double)(clock()-loop_started)/CLOCKS_PER_SEC,MANTISSA,kicks,segments,evaluations,(double)(v-initial_v),(double)((v-stop)+y[0]),
           (double)max_offset,(double)max_mid_offset,(double)max_phase,(double)max_ledger,
           (double)(-(real)137.5*(v-initial_v)+y[3]),(double)(-(real)112.5*kicks+rf_correction),(double)y[4]);
    return 0;
}
'''


def _rf_event_native_trace(ring, ramp, *, stop_nominal_turns=None, start_nominal_turns=0.,
        method=4, subdivisions=1, extended=False, sample_stride=8192, dense_tail=128):
    """Compile a temporary native recurrence; no persistent binaries or fast ramp."""
    import hashlib
    import shutil
    import subprocess
    import tempfile
    import time
    compiler = shutil.which('gcc') or shutil.which(r'C:\MinGW\bin\gcc.exe')
    if compiler is None:
        raise RuntimeError('This focused event check requires an installed C99 GCC compiler')
    stop = ramp['reference']['total_turns'] if stop_nominal_turns is None else stop_nominal_turns
    code = _rf_event_native_code(ring, ramp)
    flags = ['-std=c99', '-O3', '-fno-fast-math', '-ffp-contract=off']
    if extended:
        flags.append('-DEXTENDED')
    with tempfile.TemporaryDirectory(prefix='rf_actual_turns_') as directory:
        source, executable = Path(directory)/'event.c', Path(directory)/'event.exe'
        source.write_text(code, encoding='utf-8')
        compiled = subprocess.run([compiler, *flags, str(source), '-o', str(executable), '-lm'],
                                  text=True, capture_output=True, check=True)
        started = time.perf_counter()
        result = subprocess.run([str(executable), repr(start_nominal_turns), repr(stop), '.02', '1e-7',
            str(method), str(subdivisions), str(sample_stride), str(dense_tail)],
            text=True, capture_output=True, check=True)
        elapsed = time.perf_counter()-started
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    summary = rows.pop()
    assert summary.pop('summary') is True and rows[-1]['kind'] == 2
    assert summary['mantissa_bits'] == (64 if extended else 53)
    assert abs(summary['source_endpoint_residual_turns']) < 1e-8
    return dict(parameters=dict(start_nominal_turns=start_nominal_turns, stop_nominal_turns=stop,
        initial_phase_error_radians=.02, initial_relative_momentum_error=1e-7,
        drift_order=method, drift_subdivisions=subdivisions, extended_precision=extended,
        sample_stride_turns=sample_stride, dense_last_turns=dense_tail, extrema_stride_turns=16),
        implementation=dict(compiler=compiler, flags=flags, native_code_sha256=hashlib.sha256(code.encode()).hexdigest(),
            compilation_warnings=compiled.stderr, execution_seconds=elapsed,
            arithmetic='Compensated state/work additions, phase stored separately from the integer turn; extended run uses C long double'),
        summary=summary, initial=rows[0], endpoint=rows[-1], records=rows)


def _rf_event_short_independent_checks(ring, ramp):
    """Independent short event solves using a fresh circle root and local Ns."""
    import math
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    m = _rf_ramp_spread_model(ring, ramp)
    phase0, dp = .02, 1e-7
    rows = []
    for start_ns in (0., float(ramp['reference']['complete_turns']-32)):
        stop_ns = start_ns+16.25
        native = _rf_event_native_trace(ring, ramp, start_nominal_turns=start_ns,
            stop_nominal_turns=stop_ns, sample_stride=1, dense_tail=0)
        delta0 = phase0/(2*math.pi)
        v0 = start_ns-delta0
        gs0 = 10-m['decrement']/m['rest']*start_ns
        us0 = math.sqrt((gs0-1)*(gs0+1))
        w0 = m['rest']*us0**2*dp*(2+dp)/(math.sqrt(1+us0**2*(1+dp)**2)+gs0)
        state = np.array((0., w0, 0.))
        def rhs(local_turn, y):
            gs = (10-m['decrement']/m['rest']*start_ns)-m['decrement']/m['rest']*y[0]
            us = math.sqrt((gs-1)*(gs+1))
            gamma = gs+y[1]/m['rest']
            u = math.sqrt((gamma-1)*(gamma+1))
            r = brentq(lambda radius:radius*(1.2-.2*radius**2)-u/us,
                       .99, 1.01, xtol=5e-15, rtol=1e-14)
            ratio = gamma/(gs*(1.2-.2*r*r))
            f_actual = m['c']*u/(m['circumference']*r*gamma)
            radiation = m['charge']*gamma*u**3/(3*m['eps']*m['radius']*r)
            radiation += ramp['radiation']['longitudinal_power_watts']/(m['charge']*f_actual)
            return (ratio, m['decrement']*(ratio-(.6*r*r-.05*r**4))-radiation,
                    -m['circumference']*r*y[1]/(m['c']*m['rest']*u*gs))
        def endpoint(local_turn, y):
            return y[0]-(stop_ns-start_ns)
        endpoint.terminal = True; endpoint.direction = 1
        local = 0.; count = 0; evaluations = 0
        next_gap = math.floor(v0)+1
        while True:
            target = next_gap-v0
            solve = solve_ivp(rhs, (local, target), state, method='DOP853', rtol=3e-13,
                atol=(2e-13, 2e-10, 1e-22), max_step=.1, events=endpoint)
            assert solve.success
            evaluations += solve.nfev
            local, state = solve.t[-1], solve.y[:, -1].copy()
            if len(solve.t_events[0]):
                break
            gs = (10-m['decrement']/m['rest']*start_ns)-m['decrement']/m['rest']*state[0]
            us = math.sqrt((gs-1)*(gs+1))
            fs = m['c']*us/(m['circumference']*gs)
            radiation_s = m['charge']*gs*us**3/(3*m['eps']*m['radius'])
            radiation_s += ramp['radiation']['longitudinal_power_watts']/(m['charge']*fs)
            phase = phase0+2*math.pi*(state[0]-local)
            state[1] += 2*(112.5-radiation_s)*math.sin(7*math.pi/6+phase)
            count += 1; next_gap += 1
        errors = dict(phase_radians=abs(phase0+2*math.pi*(state[0]-local)-native['endpoint']['phase_error_radians']),
            energy_eV=abs(state[1]-native['endpoint']['energy_error_eV']),
            proper_clock_seconds=abs(state[2]-native['endpoint']['proper_clock_error_seconds']))
        assert count == native['summary']['gaps'] and errors['phase_radians'] < 2e-9
        assert errors['energy_eV'] < 2e-5 and errors['proper_clock_seconds'] < 1e-15
        rows.append(dict(start_nominal_turn=start_ns, stop_nominal_turn=stop_ns,
            independent_method='DOP853 with local reference-turn state and a fresh Brent circle root at every RHS; exact event RF phase',
            independent_evaluations=evaluations, gap_count=count, errors=errors, native=native))
    return rows


def _rf_saved_half_energy_processing(ring, ramp, spread):
    """Reprocess saved Lorentz windows at exact midpoint-of-energy-kick sections."""
    import math
    import numpy as np
    from scipy.integrate import solve_ivp
    m = _rf_ramp_spread_model(ring, ramp)
    radius, ui, rest = m['radius'], m['ui'], m['rest']
    field0, scale = m['field0'], m['canonical_scale']
    rows, curves = [], {}
    for window in spread['windows']:
        p = window['plot']
        times = window['start_lab_seconds']+np.asarray(p['local_lab_seconds'])
        angles = np.asarray(p['angle_radians'])
        r = 1+np.asarray(p['radial_offset_metres'])/radius
        z = np.asarray(p['vertical_offset_metres'])/radius
        rg = 1+np.asarray(p['projected_guiding_offset_metres'])/radius
        field = np.asarray(m['reference'](times)[2])
        ap = .6*r-.05*r**3+.2*r*z*z
        physical_p = field/field0*(-.6*rg*rg+.15*rg**4)
        indices = np.where((np.abs(angles/(2*math.pi)-np.round(angles/(2*math.pi))) < 2e-13)
                           & (np.arange(len(angles)) > 0))[0]
        assert len(indices) == len(window['kicks']) == 16
        previous_post = physical_p[0]
        midpoint_canonical = []
        for index, kick in zip(indices, window['kicks']):
            pre = ui/r[index]*(previous_post-field[index]/field0*r[index]*ap[index])
            post = ui/r[index]*(physical_p[index]-field[index]/field0*r[index]*ap[index])
            square = (pre*pre+post*post)/2-(kick['work_eV']/(2*rest))**2
            assert pre < 0 and post < 0 and square > 0
            middle = -math.sqrt(square)
            midpoint_canonical.append((r[index]*middle/ui+
                field[index]/field0*(r[index]*ap[index]+.45))/scale)
            previous_post = physical_p[index]
        t = times[indices]
        phase = np.asarray(p['phase_error_radians'])[indices]
        canonical = np.asarray(midpoint_canonical)
        solve = solve_ivp(lambda t, y:m['rhs'](7*math.pi/6, t, y),
            (t[0], t[-1]), (phase[0], canonical[0]), method='DOP853',
            rtol=3e-13, atol=1e-14, t_eval=t)
        assert solve.success and len(solve.t) == 16
        residual = radius*(m['circular'](t, canonical)[0]-m['circular'](t, solve.y[1])[0])
        phase_residual = phase-solve.y[0]
        key = (window['nominal_start_turn'], bool(window['fine']))
        curves[key] = (residual, phase_residual)
        rows.append(dict(nominal_start_turn=key[0], fine=key[1], centered_sections=16, compared_intervals=15,
            first_lab_seconds=float(t[0]), last_lab_seconds=float(t[-1]),
            first_centered_phase_radians=float(phase[0]), first_centered_canonical_error_scaled=float(canonical[0]),
            maximum_guiding_residual_metres=float(np.max(np.abs(residual))),
            final_guiding_residual_metres=float(residual[-1]),
            maximum_phase_residual_radians=float(np.max(np.abs(phase_residual))),
            lab_seconds=t.tolist(), guiding_residual_metres=residual.tolist(), phase_residual_radians=phase_residual.tolist(),
            prior_raw_maximum_guiding_residual_metres=window['maximum_projected_guiding_difference_from_average_metres']))
    changes = []
    for turn in sorted({key[0] for key in curves}):
        coarse, fine = curves[(turn, False)], curves[(turn, True)]
        changes.append(dict(nominal_start_turn=turn,
            maximum_same_gap_guiding_residual_change_metres=float(np.max(np.abs(coarse[0]-fine[0]))),
            maximum_same_gap_phase_residual_change_radians=float(np.max(np.abs(coarse[1]-fine[1])))))
    return dict(windows=rows, coarse_fine_comparisons=changes,
        equations='E_mid=(E_pre+E_post)/2 at fixed transverse mechanical momentum; pphi_mid=-sqrt[(pphi_pre^2+pphi_post^2)/2-(DeltaE/(2mc^2))^2] in units mc. Recover each pre-gap Pphi from the prior post-gap value because these axisymmetric guide windows omit radiation reaction',
        scope='Reuses saved radiation-free Lorentz/thin-gap windows. The local average is restarted at the first midpoint-of-energy-kick section; sixteen sections give fifteen subsequent intervals. This removes the leading sampling/micromotion offset, not transverse motion, finite-gap effects or full-ramp Lorentz uncertainty',
        independent_derivation='The section reconstruction and independent numerical calculation were supplied by mass_audit; this helper reproduces them from the saved source data')


def _rf_event_matched_average(ring, ramp, event_trace, rtol=3e-12):
    """Compare both raw and consistently centered initializations at gap times."""
    import math
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    m = _rf_ramp_spread_model(ring, ramp)
    gaps = [row for row in event_trace['records'] if row['kind'] == 1]
    first, initial = gaps[0], event_trace['initial']
    def canonical(time, offset):
        us = m['reference'](time)[0]
        x = offset/m['radius']; r2m = x*(2+x)
        return float(us/m['ui']*r2m*(-.3+.15*r2m)/m['canonical_scale'])
    raw_initial = (initial['phase_error_radians'], canonical(0., initial['guiding_offset_metres']))
    tfirst = first['lab_seconds']
    prelude = solve_ivp(lambda t,y:m['rhs'](7*math.pi/6,t,y), (0.,tfirst), raw_initial,
        method='DOP853', rtol=3e-13, atol=1e-14)
    assert prelude.success
    raw_at_first = prelude.y[:, -1]
    matched = (first['phase_error_radians'], canonical(tfirst, first['midpoint_guiding_offset_metres']))
    state0 = (raw_at_first[0], matched[0], raw_at_first[1], matched[1])
    started = time.perf_counter()
    tolerance_scale = rtol/3e-12
    atols = [value*tolerance_scale for value in (3e-14,3e-14,3e-13,3e-13)]
    solve = solve_ivp(lambda t,y:m['rhs'](7*math.pi/6,t,y), (tfirst,m['duration']), state0,
        method='DOP853', rtol=rtol, atol=atols,
        max_step=.0005, dense_output=True)
    assert solve.success and solve.t[-1] == m['duration']
    seconds = time.perf_counter()-started
    times = np.asarray([row['lab_seconds'] for row in gaps])
    states = solve.sol(times)
    offsets = m['radius']*m['circular'](times, states[2:])[1]
    raw = np.asarray([row['guiding_offset_metres'] for row in gaps])
    mid = np.asarray([row['midpoint_guiding_offset_metres'] for row in gaps])
    phases = np.asarray([row['phase_error_radians'] for row in gaps])
    difference = mid-offsets[1]
    raw_difference = raw-offsets[0]
    phase_difference = phases-states[1]
    whole = dict(lab_seconds=times.tolist(), orbital_turns=[row['orbital_turn'] for row in gaps],
        phase_error_radians=phases.tolist(), raw_guiding_offset_metres=raw.tolist(),
        midpoint_guiding_offset_metres=mid.tolist(), averaged_matched_guiding_offset_metres=offsets[1].tolist(),
        averaged_raw_initial_guiding_offset_metres=offsets[0].tolist(),
        midpoint_minus_matched_average_metres=difference.tolist(),
        midpoint_minus_matched_average_phase_radians=phase_difference.tolist())
    last_turn = gaps[-1]['orbital_turn']
    late_indices = np.where(np.asarray(whole['orbital_turns']) >= last_turn-2047)[0]
    late = {key:np.asarray(values)[late_indices].tolist() for key, values in whole.items()}
    final = solve.y[:, -1]
    final_offsets = m['radius']*m['circular'](m['duration'], final[2:])[1]
    return dict(method='DOP853, actual lab-time canonical guiding-circle ODE; two initializations integrated together',
        rtol=rtol, atol=atols, max_step_seconds=.0005,
        execution_seconds=seconds, evaluations=int(solve.nfev), accepted_steps=len(solve.t)-1,
        first_midpoint_initialization=dict(lab_seconds=tfirst, phase_error_radians=matched[0],
            canonical_error_scaled=matched[1], energy_error_eV=first['midpoint_energy_error_eV']),
        maximum_sampled_midpoint_guiding_difference_metres=float(np.max(np.abs(difference))),
        maximum_sampled_midpoint_phase_difference_radians=float(np.max(np.abs(phase_difference))),
        maximum_sampled_raw_postkick_difference_metres=float(np.max(np.abs(raw_difference))),
        final_matched_average_phase_radians=float(final[1]),
        final_matched_average_guiding_offset_metres=float(final_offsets[1]),
        final_raw_initial_average_phase_radians=float(final[0]),
        final_raw_initial_average_guiding_offset_metres=float(final_offsets[0]),
        plot=whole, late_plot=late,
        scope='Averages and map share source time at each sampled midpoint-of-energy-kick section. The matched average starts at the first centered event; the separate raw-initialized average illustrates initialization plus micromotion effects. Sparse whole-ramp samples do not resolve every synchrotron oscillation. The last2048 consecutive gap sections are resolved. This comparison tests a reduced longitudinal model, not transverse acceptance or a finite RF cavity')


def proton_rf_event_ramp_checks(ring, ramp, spread, records=None):
    """One actual-ramp event trajectory with drift, roundoff and section checks."""
    import copy
    import hashlib
    import math
    import numpy as np
    import sympy as sp
    m = _rf_ramp_spread_model(ring, ramp)
    if records is None:
        pilots = [_rf_event_native_trace(ring, ramp, stop_nominal_turns=100000.,
            method=order, extended=extended, sample_stride=10000, dense_tail=4)
            for order, extended in ((2,False), (4,False), (2,True))]
        shorts = _rf_event_short_independent_checks(ring, ramp)
        main = _rf_event_native_trace(ring, ramp, method=2, dense_tail=2048)
        fine = _rf_event_native_trace(ring, ramp, method=4, dense_tail=2048)
        extended = _rf_event_native_trace(ring, ramp, method=4, extended=True, dense_tail=2048)
        average = _rf_event_matched_average(ring, ramp, extended)
        average_coarse = _rf_event_matched_average(ring, ramp, extended, rtol=3e-10)
    else:
        pilots, shorts, main, fine, extended, average, average_coarse = [copy.deepcopy(records[key])
            for key in ('pilots','short_checks','main','fine','extended','matched_average','matched_average_coarse')]
    processing = _rf_saved_half_energy_processing(ring, ramp, spread)
    r = sp.symbols('r', positive=True)
    fshape = sp.Rational(6,5)-r*r/5
    ap = 3*r/5-r**3/20
    ell, gshape = r*ap-r*r*fshape, r*fshape
    assert sp.simplify(sp.diff(ell,r)+r*sp.diff(gshape,r)) == 0
    assert sp.simplify(gshape-sp.diff(gshape,r)*ell/sp.diff(ell,r)-ap) == 0
    assert sp.simplify(r*ap-(sp.Rational(11,20)+(r*r-1)/2-(r*r-1)**2/20)) == 0
    code_hash = hashlib.sha256(_rf_event_native_code(ring, ramp).encode()).hexdigest()
    assert len(pilots) == 3 and len(shorts) == 2
    for trace in pilots+[row['native'] for row in shorts]+[main,fine,extended]:
        assert trace['implementation']['native_code_sha256'] == code_hash
        assert trace['parameters']['initial_phase_error_radians'] == .02
        assert trace['parameters']['initial_relative_momentum_error'] == 1e-7
        assert trace['parameters']['extrema_stride_turns'] == 16
        assert all(np.isfinite(row[key]) for row in trace['records'] for key in row)
        assert trace['initial']['kind'] == 0 and trace['endpoint']['kind'] == 2
    for trace, order, precision in ((main,2,False),(fine,4,False),(extended,4,True)):
        assert trace['parameters']['start_nominal_turns'] == 0.
        assert trace['parameters']['stop_nominal_turns'] == ramp['reference']['total_turns']
        assert trace['parameters']['drift_order'] == order
        assert trace['parameters']['drift_subdivisions'] == 1
        assert trace['parameters']['extended_precision'] == precision
        assert trace['summary']['gaps'] == math.floor(trace['endpoint']['orbital_turn'])-math.floor(trace['initial']['orbital_turn'])
        expected_turns = (ramp['reference']['total_turns']+
            (trace['initial']['phase_error_radians']-trace['endpoint']['phase_error_radians'])/(2*math.pi))
        assert abs(trace['summary']['actual_revolutions']-expected_turns) < 1e-8
        assert abs(trace['endpoint']['lab_seconds']-m['duration']) < 1e-12
        assert abs(trace['endpoint']['proper_seconds']-ramp['reference']['proton_proper_elapsed_seconds']-
                   trace['endpoint']['proper_clock_error_seconds']) < 1e-12
        assert trace['summary']['maximum_work_residual_eV'] < 1e-5
        assert np.all(np.diff([row['lab_seconds'] for row in trace['records']]) > 0)
        for row in trace['records']:
            if row['kind'] == 1:
                assert row['orbital_turn'] == round(row['orbital_turn'])
                assert abs(row['midpoint_energy_error_eV']-(row['energy_error_eV']-row['rf_gap_work_eV']/2)) < 1e-9
    comparisons = []
    for label, a, b in (('midpoint_drift_vs_RK4',main,fine),('double_vs_extended_RK4',fine,extended)):
        ga, gb = [[row for row in trace['records'] if row['kind'] == 1] for trace in (a,b)]
        assert [row['orbital_turn'] for row in ga] == [row['orbital_turn'] for row in gb]
        errors = {key:max(abs(x[key]-y[key]) for x,y in zip(ga,gb)) for key in
            ('phase_error_radians','energy_error_eV','midpoint_guiding_offset_metres','proper_clock_error_seconds','lab_seconds')}
        assert errors['midpoint_guiding_offset_metres'] < 3e-9
        assert errors['phase_error_radians'] < 1e-8
        assert errors['energy_error_eV'] < .001
        comparisons.append(dict(comparison=label, matching_gap_sections=len(ga), maximum_same_gap_difference=errors,
            endpoint_difference={key:abs(a['endpoint'][key]-b['endpoint'][key]) for key in
                ('phase_error_radians','energy_error_eV','guiding_offset_metres','proper_seconds')}))
    assert average['rtol'] == 3e-12 and average_coarse['rtol'] == 3e-10
    assert average['plot']['lab_seconds'] == average_coarse['plot']['lab_seconds']
    mean_errors = {key:float(np.max(np.abs(np.array(average['plot'][key])-np.array(average_coarse['plot'][key]))))
        for key in ('averaged_matched_guiding_offset_metres','midpoint_minus_matched_average_phase_radians')}
    assert mean_errors['averaged_matched_guiding_offset_metres'] < 3e-9
    assert mean_errors['midpoint_minus_matched_average_phase_radians'] < 3e-8
    assert len(average['late_plot']['orbital_turns']) == 2048
    assert np.all(np.diff(average['late_plot']['orbital_turns']) == 1)
    old_reference = spread['fine']['cases'][3]
    prior_errors = dict(phase_radians=abs(average['final_raw_initial_average_phase_radians']-
            old_reference['final_phase_error_radians']),
        guiding_offset_metres=abs(average['final_raw_initial_average_guiding_offset_metres']-
            old_reference['final_guiding_offset_metres']))
    assert prior_errors['phase_radians'] < 2e-10 and prior_errors['guiding_offset_metres'] < 3e-11
    end = extended['endpoint']
    last_gap = [row for row in extended['records'] if row['kind'] == 1][-1]
    return dict(status='verified bounded event-resolved longitudinal study; one declared sample over the real ramp',
        parameters=dict(initial_phase_error_radians=.02, initial_relative_momentum_error=1e-7,
            actual_source_ramp_seconds=m['duration'], net_reference_energy_decrement_eV=250.,
            source_phase_radians=7*math.pi/6, source_peak_voltage='2*(112.5-Urad_s[eV]) volts',
            reference_final_kinetic_energy_eV=5.3e6, declared_guiding_aperture_half_width_metres=.01,
            representative_count=1, rf_harmonic=1),
        equations=dict(section='Clockwise orbital turns v=theta/(2*pi); canonical gap variables are (t,-E). Instantaneous circular reduction of K_theta=e*r_phys*Aphi-r_phys*sqrt(E^2/c^2-m^2*c^2-pr^2-pz^2). No transverse action is evolved',
            circle='p/(mc)=u_s*r*(1.2-.2*r^2); F=r*a_phi=.6*r^2-.05*r^4; L=e*B0*R^2*(-.6*r^2+.15*r^4)',
            drift='delta=N_s(t)-v, w=E-E_s(N_s) in eV; Rf=f_s/f_actual=gamma/[gamma_s*(1.2-.2*r^2)]; ddelta/dv=Rf-1; dw/dv=D*(Rf-F)-Urad',
            induction='Between gaps dE/dv=-D*F-Urad: the induction field removes energy and is included exactly in this reduced drift. This follows from Bdot=-D/(R*C) with consistent volts/eV conventions for a proton',
            kick='At every integer v: w_plus=w_minus+Vhat(t)*sin(phi_s+2*pi*delta). Neither clock advances across the ideal gap. Source phase is evaluated at the actual event time',
            clocks='t=C*M*(u_i-u_s)/(c*D); d[tau-tau_s(N_s)]/dv=-C*r*w/(c*M*u*gamma_s), with M=mc^2 in eV',
            endpoint='Start v=-delta0 so N_s=0. Stop on N_s=v+delta=N_total by a partial final drift; physical revolutions=N_total+delta0-delta_final. The final particle energy includes w_final',
            centering='At a gap E_mid=(E_minus+E_plus)/2. Project this energy onto the circular family; arithmetic averaging of Pphi is not the exact section processor. Initialize the comparison average at the same first midpoint section'),
        pilots=pilots, short_checks=shorts, main=main, fine=fine, extended=extended,
        matched_average=average, matched_average_coarse=average_coarse, saved_window_processing=processing,
        accepted=dict(method='extended-precision RK4 drift with exact ideal RF event kicks',
            actual_revolutions=extended['summary']['actual_revolutions'], gap_encounters=extended['summary']['gaps'],
            final_partial_orbital_turn=end['orbital_turn']-last_gap['orbital_turn'],
            final_lab_seconds=end['lab_seconds'], final_proper_seconds=end['proper_seconds'],
            final_energy_error_eV=end['energy_error_eV'], final_kinetic_energy_eV=5.3e6+end['energy_error_eV'],
            maximum_sampled_midpoint_vs_matched_average_metres=average['maximum_sampled_midpoint_guiding_difference_metres'],
            maximum_sampled_midpoint_vs_matched_average_phase_radians=average['maximum_sampled_midpoint_phase_difference_radians'],
            maximum_sampled_raw_postkick_vs_raw_initial_average_metres=average['maximum_sampled_raw_postkick_difference_metres'],
            maximum_sampled_midpoint_guiding_offset_metres=max(extended['summary']['maximum_midpoint_guiding_offset_metres'],
                max(abs(x) for x in average['plot']['midpoint_guiding_offset_metres'])),
            native_every16_turn_pre_post_guiding_maximum_metres=extended['summary']['maximum_pre_post_guiding_offset_metres']),
        verification=dict(symbolic_canonical_induction_identity=True, finite_source_rate_unchanged=True,
            native_refinements=comparisons, matched_average_refinement=mean_errors,
            prior_averaged_representative_endpoint_comparison=prior_errors,
            extended_work_ledger_maximum_residual_eV=extended['summary']['maximum_work_residual_eV'],
            source_endpoint_error_seconds=abs(end['lab_seconds']-m['duration']),
            maximum_short_independent_energy_error_eV=max(row['errors']['energy_eV'] for row in shorts),
            source_code_sha256=code_hash),
        work_accounting=dict(induction_work_eV=extended['summary']['induction_work_eV'],
            rf_work_eV=extended['summary']['rf_work_eV'], radiated_energy_eV=extended['summary']['radiated_energy_eV'],
            equation='Delta E=W_ind+sum W_RF-Urad. Compensated corrections to the large nominal work sums are accumulated separately to test this identity without subtracting GeV-scale quantities',
            scope='Leading smooth-orbit radiation diagnostic retained as a dissipative term. No singular thin-gap radiation or actual cavity power/source-energy calculation is supplied'),
        interpretation='Discrete kicks can be resolved efficiently at the actual ramp rate in this longitudinal reduction. Matched midpoint sections remove most of the raw kick/average offset; a nonzero accumulated difference remains and is well above the measured drift, precision and mean-reference changes',
        limitations='One fixed initial sample, not a beam-acceptance theorem. Relative to a physical fixed-position gap kick, the circularly reduced canonical Pphi jump differs at second order in kick energy; the displacement relative to the new circle is first order, while its radial action/energy is quadratic. That transverse motion is omitted here. The separately derived prior transverse envelope is not promoted to a full event-map confinement proof. Extrema are sampled: native checks use every16th pre/post/mid gap, while saved late data contain every gap. The terminal point is generally between gaps and is not relabeled as a midpoint section. No finite Maxwell RF cavity, ramp-end transition, injection, cooling, trapping or capture is established',
        sources=copy.deepcopy(spread['sources']),
        derivation_attribution='The induction-aware event equations, stable clock variables, exact half-energy processor and native implementation are calculations for this stated ideal model; mass_audit independently checked the canonical equations and supplied the saved-window processor')


def run_checks():
    with mp.workdps(60):
        c, g0, day, mass = map(mp.mpf, ('299792458', '9.80665', '86400', '1e-15'))
        cases = {}
        for name in ('one_g_one_lab_day', 'tenfold_one_lab_day'):
            h = (mp.asinh(g0*day/(4*c)) if name.startswith('one_g')
                 else mp.findroot(lambda h: mp.sinh(h)/h-10, (4, 5)))
            a = 4*c*mp.sinh(h)/day
            tau = 4*c*h/a
            gamma = mp.cosh(h)
            distance = 4*c*c/a*mp.sinh(h/2)**2
            energy = 2*mp.sinh(h/2)**2*mass*c*c
            # ponytail: prescribed acceleration jumps; finite-jerk propulsion is a separate model.
            legs = ((0, 1), (h, -1), (0, -1), (-h, 1))
            times = [quad(lambda s: float(mp.cosh(start+sign*h*s)), 0, 1,
                          epsabs=1e-11, epsrel=1e-12)[0]*float(tau/4)
                     for start, sign in legs]
            distances = [quad(lambda s: float(mp.sinh(start+sign*h*s)), 0, 1,
                              epsabs=1e-11, epsrel=1e-12)[0]*float(c*tau/4)
                         for start, sign in legs]
            assert abs(sum(times)/float(day)-1)<1e-12
            assert abs(sum(distances)/float(distance))<1e-12
            assert abs(sum(distances[:2])/float(distance)-1)<1e-12
            assert abs(4*c/a*mp.sinh(h)-day)<mp.mpf('1e-48')
            assert day>tau>0 and 0<mp.tanh(h)<1 and energy>0
            assert abs(gamma**2*(1-mp.tanh(h)**2)-1)<mp.mpf('1e-50')
            if name.startswith('one_g'):
                assert abs(a-g0)<mp.mpf('1e-50')
            else:
                assert abs(day/tau-10)<mp.mpf('1e-50')
            cases[name] = {key: mp.nstr(value, 30) for key, value in dict(
                lab_elapsed_seconds=day, probe_elapsed_seconds=tau,
                age_difference_seconds=day-tau, lab_to_probe_time_ratio=day/tau,
                proper_acceleration_m_per_s2=a, proper_acceleration_in_g=a/g0,
                peak_rapidity=h, peak_gamma=gamma, peak_speed_over_c=mp.tanh(h),
                maximum_distance_metres=distance,
                peak_payload_kinetic_energy_joules=energy).items()}
            cases[name]['independent_integral_lab_time_error_seconds'] = abs(sum(times)-float(day))
        ring = proton_ring_checks()
        kicker = proton_kicker_checks(ring)
        finite_coils = proton_finite_coil_checks(ring)
        retarded_control = proton_retarded_control_checks(ring, finite_coils)
        source_radiation = proton_control_source_feasibility_checks(retarded_control)
        lower_radiation = proton_lower_radiation_control_checks(
            ring, finite_coils, retarded_control, source_radiation)
        conductors = proton_finite_conductor_checks(lower_radiation)
        localized = proton_localized_control_checks(ring, finite_coils, lower_radiation, conductors)
        torus = proton_finite_torus_control_checks(ring, finite_coils, localized)
        launch_study = proton_torus_launch_perturbation_checks(ring, finite_coils, torus)
        receiver = proton_electrostatic_receiver_checks(ring, finite_coils, torus, launch_study)
        regularized = proton_regularized_receiver_checks(ring, finite_coils, torus, receiver)
        duration_scan = proton_receiver_duration_scan_checks(ring, finite_coils, regularized)
        staged_rf = proton_staged_rf_deceleration_checks(ring)
        rf_spread = proton_rf_ramp_spread_checks(ring, staged_rf)
        return dict(
            model='Flat-spacetime roundtrip: start/end at rest at the same inertial station. Four equal proper-time legs with signed proper acceleration +a,-a,-a,+a',
            equations='h=a*tau/(4c); t=4c*sinh(h)/a; xmax=2c^2*(cosh(h)-1)/a; gamma_max=cosh(h); K_max=(gamma_max-1)*m*c^2',
            derivation='Integrate dt/dtau=cosh(theta), dx/dtau=c*sinh(theta), dtheta/dtau=+/-a/c. Rapidity follows0->h->0->-h->0, so the trajectory reunites at rest',
            proton_ring=ring,
            proton_focusing=proton_focusing_checks(ring['nominal']),
            proton_slow_loss=proton_slow_loss_checks(ring),
            proton_rf=proton_rf_checks(ring),
            proton_kicker=kicker,
            proton_outgoing_candidate=proton_outgoing_candidate_checks(ring, kicker),
            proton_finite_coil=finite_coils,
            proton_finite_coil_escape=proton_finite_coil_escape_checks(ring, finite_coils),
            proton_retarded_control=retarded_control,
            source_feasibility=source_radiation,
            lower_radiation_control=lower_radiation,
            finite_conductor=conductors,
            localized_control=localized,
            finite_torus_control=torus,
            torus_launch_perturbations=launch_study,
            electrostatic_receiver=receiver,
            switched_receiver=proton_switched_receiver_checks(ring, finite_coils, receiver),
            regularized_receiver=regularized,
            receiver_duration_scan=duration_scan,
            receiver_radiation_spectrum=proton_receiver_radiation_spectrum_checks(duration_scan),
            staged_rf_deceleration=staged_rf,
            rf_ramp_spread=rf_spread,
            rf_event_ramp=proton_rf_event_ramp_checks(ring, staged_rf, rf_spread),
            payload_mass_kg=str(mass),
            mass_scope='1 picogram is a mass normalization, not a specified solid grain. Clock ratio, required acceleration and route size are independent of test mass; payload kinetic energy scales linearly with mass',
            cases=cases,
            date_example='If the inertial station calls departure13 September2027, reunion is14 September2027 after86400 station seconds. In the tenfold case only2.4hours of probe proper time pass. This does not return to12 September2027',
            scope='Exact ideal special-relativistic trajectory, not a demonstrated flight. Acceleration is finite but switches discontinuously. No propulsion, fuel, field confinement, radiation, collisions, bonding or finite-size survival model; gravity of Earth/Sun and station motion are neglected. K_max is only peak payload kinetic energy, not total mission energy',
            source='Chou et al., Optical Clocks and Relativity, Science329,1630-1633(2010), doi10.1126/science.1192720; https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=905055. Experiments support time dilation; this four-leg calculation is our application, not their flight experiment')


if __name__ == '__main__':
    result = run_checks()
    Path(__file__).with_name('future_time_results.json').write_text(
        json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result['cases'], indent=2))
