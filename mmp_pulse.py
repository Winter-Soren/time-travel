"""Linear MMP/JT pulse, with a separate first-interval boundary/CFT response test."""
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

WEAK_SOURCE_SCOPE = 'Declared alpha sets this observer-specific weak-source screen. Order1 fails the assumed small-source hierarchy, not a universal JT-validity, horizon or collapse criterion. Small values do not prove global stability. Omitted4D terms, trace corrections, fluctuations, finite-mouth/exterior geometry and preparation remain untested'


def time_offset_checks():
    """Kinematic helical-cylinder family and a necessary static-throat source test."""
    L, delta, q = sp.symbols('L delta q', positive=True)
    ratio = delta/L
    metric = sp.diag(-1, 1)
    boost = sp.Matrix([[1, -ratio], [-ratio, 1]])/sp.sqrt(1-ratio**2)
    assert sp.simplify(boost.T*metric*boost-metric) == sp.zeros(2)
    shifted = sp.simplify(boost*sp.Matrix([delta, L]))
    assert shifted[0] == 0
    assert sp.simplify(shifted[1]**2-(L**2-delta**2)) == 0
    rho = -sp.pi*q/(6*(L**2-delta**2))
    stress = sp.simplify(boost.T*(rho*sp.eye(2))*boost)
    null_jacobian = sp.Matrix([[1, 1], [1, -1]])/2  # t=(u+v)/2, x=(u-v)/2.
    null_stress = sp.simplify(null_jacobian.T*stress*null_jacobian)
    expected = sp.diag(-sp.pi*q/(12*(L+delta)**2), -sp.pi*q/(12*(L-delta)**2))
    assert sp.simplify(null_stress-expected) == sp.zeros(2)
    invariant = sp.factor(sp.trace(stress*metric*stress*metric))
    assert sp.simplify(invariant-2*rho**2) == 0
    t, x = sp.symbols('t x', real=True)
    clock = t-ratio*x
    assert sp.simplify(clock.subs({t:t+delta, x:x+L}, simultaneous=True)-clock) == 0
    assert sp.simplify(-sp.diff(clock, t)**2+sp.diff(clock, x)**2+1-ratio**2) == 0
    # Any time-independent dilaton in this fixed conformal metric has equal null constraints.
    chi = sp.symbols('chi', real=True); phi = sp.Function('phi')(chi)
    du = lambda f: sp.diff(f, chi)/2
    dv = lambda f: -sp.diff(f, chi)/2
    assert sp.simplify(du(du(phi))-sp.tan(chi)*du(phi)
                       -dv(dv(phi))-sp.tan(chi)*dv(phi)) == 0
    tauu, tauv = expected[0, 0]+q/(48*sp.pi), expected[1, 1]+q/(48*sp.pi)
    difference = sp.factor(tauu-tauv)
    assert sp.simplify(difference-sp.pi*q*L*delta/(3*(L**2-delta**2)**2)) == 0
    tau, jp, jm = sp.symbols('tau j_plus j_minus', real=True)
    du_full = lambda f: (sp.diff(f, tau)+sp.diff(f, chi))/2
    dv_full = lambda f: (sp.diff(f, tau)-sp.diff(f, chi))/2
    forced = -(jp+jm)*(1+chi*sp.tan(chi))-(jp-jm)*tau*sp.tan(chi)
    forced_residuals = [sp.simplify(sp.trigsimp(f)) for f in (
        du_full(du_full(forced))-sp.tan(chi)*du_full(forced)+jp,
        dv_full(dv_full(forced))+sp.tan(chi)*dv_full(forced)+jm,
        du_full(dv_full(forced))+forced/(2*sp.cos(chi)**2))]
    assert forced_residuals == [0, 0, 0]
    # The invariant changes only if the invariant loop length changes, not under a mere boost.
    ell, rapidity = sp.symbols('ell rapidity', positive=True)
    assert sp.simplify((L**2-delta**2).subs({L:ell*sp.cosh(rapidity), delta:ell*sp.sinh(rapidity)})-ell**2) == 0
    samples = []
    length = math.pi+.1
    for a in (0., .1, .5, .9, .99, .999):
        long, short = length*(1+a), length*(1-a)
        plus, minus = 1/6-2*math.pi**2/(3*long**2), 1/6-2*math.pi**2/(3*short**2)
        samples.append(dict(offset_over_loop_length=a, shorter_return_gap=short,
            proper_circumference_over_L=math.sqrt(1-a*a),
            short_chiral_vacuum_magnitude_ratio=1/(1-a)**2,
            rest_frame_energy_magnitude_ratio=1/(1-a*a),
            invariant_stress_squared_ratio=1/(1-a*a)**2,
            lab_energy_magnitude_ratio=(1+a*a)/(1-a*a)**2,
            full_throat_null_plus_over_alpha=plus, full_throat_null_minus_over_alpha=minus,
            static_source_mismatch_over_alpha=plus-minus))
    return dict(
        model='Separate fixed-parameter flat1+1 quotient family (t,x)~(t+Delta,x+L), c=hbar=1, with q complex massless fermions in the antiperiodic rest-cylinder vacuum. Not a dynamical formation simulation or a constructed MMP extension',
        units='The JT comparison and sample table use L=L_physical/ell0,Delta=Delta_physical/ell0 and dimensionless throat null coordinates. Restoring a physical2D stress requires dividing the dimensionless stress by ell0^2. The independent flat-cylinder formulas also hold in physical natural units if all lengths and coordinate components are restored consistently',
        invariant_loop='ell_loop^2=L^2-Delta^2. For abs(Delta)<L, T=t-(Delta/L)*x is single-valued and has timelike gradient; a global time exists. At equality the compact orbit is null; above it the orbit is timelike',
        vacuum='In the rest frame rho=p=-pi*q/(6*(L^2-Delta^2)). With u=t+x,v=t-x: f_uu=-pi*q/[12*(L+Delta)^2], f_vv=-pi*q/[12*(L-Delta)^2]',
        stress_invariant='T_ab*T^ab=pi^2*q^2/[18*(L^2-Delta^2)^2] in the flat quotient. At fixed finite L it diverges as abs(Delta)->L from below; this is not merely a boosted observer energy',
        boost_control='Holding ell_loop fixed while L=ell_loop*cosh(theta),Delta=ell_loop*sinh(theta) preserves the invariant vacuum stress. Delta/L approaches1 only with unbounded L,Delta, and no finite null identification is reached',
        static_throat_check='Adding the identical q/(48*pi) throat Weyl term to both chiral stresses leaves T_uu-T_vv=pi*q*L*Delta/[3*(L^2-Delta^2)^2]. Any static phi(chi) in the old fixed AdS2 conformal chart has equal uu/vv JT constraints, so the shifted vacuum cannot source that unchanged static throat for Delta!=0',
        source_normalization='Total JT sources j_plus/minus=(2*G/r_e^2)*(f_plus/minus+q/(48*pi)); alpha=G*q/(4*pi*r_e^2). Table full_throat_null entries are j_plus/minus divided by alpha, not raw2D stress divided by alpha. Sample L=pi+0.1',
        forced_response='If unequal total null sources j_plus,j_minus are artificially held constant, phi=beta*(1+chi*tan(chi))-d*tau*tan(chi), beta=-(j_plus+j_minus),d=j_plus-j_minus, solves all three local JT equations. The center stays constant while d/dtau(bR-bL)=-2*d. Homogeneous modes bounded in time cannot remove this drift. This is a fixed-source diagnostic, not a closed evolving vacuum/mouth solution',
        forced_symbolic_residuals=[str(r) for r in forced_residuals],
        interpretation='A physical mouth offset changes the composed sewing/return maps. Relabeling mouth clocks while transforming all sewing maps consistently changes no invariant return time. A nonzero imposed offset with maps held fixed is a different physical ansatz, not automatically a coordinate relabeling',
        samples=samples,
        scope='The vacuum result is state-dependent and the fixed-L scan is a family of backgrounds, not an evolution through a chronology horizon. Real MMP field-line delays differ; the fastest exterior causal route can differ from a chosen fermion loop. Divergent short-period modes also leave the low-energy approximation. Do not identify this toy divergence with a proved4D chronology threshold or protection theorem. Other states, geometry, counterflow, exterior dynamics and backreaction are not solved. This rejects retaining the old static source/throat unchanged while adding this offset',
        sources=['https://arxiv.org/pdf/1807.04726 (MMP: antiperiodic cylinder energy, local anomaly and causal background)',
                 'https://doi.org/10.1103/PhysRevLett.61.1446 (Morris-Thorne-Yurtsever: physical wormhole time-offset proposal; not this quantum source solution)',
                 'https://arxiv.org/pdf/gr-qc/9603012 (Kay-Radzikowski-Wald: separate compactly generated horizon theorem for linear scalar fields; not automatically applicable to this fermion quotient family)'])


def stationary_balance_checks():
    """Sharp stationary NS mean-source bound, explicit states, and conditional throat matching."""
    import mpmath as mp
    from scipy.optimize import brentq
    long, short, q = sp.symbols('P_long P_short q', positive=True)
    probability = (1-short**2/long**2)/24
    f_long = -sp.pi*q/(12*long**2)
    f_short = 2*sp.pi*(q*probability-q/24)/short**2
    assert sp.simplify(f_short-f_long) == 0
    beta_over_alpha = -16*sp.pi*(f_long+q/(48*sp.pi))/q
    assert sp.simplify(beta_over_alpha-((2*sp.pi/long)**2-1)/3) == 0
    p = sp.symbols('p', nonnegative=True)
    # Basis: vacuum, particle, antiparticle, neutral pair in the lowest NS modes.
    density = sp.diag(1-p, 0, 0, p)
    hamiltonian, charge_op = sp.diag(0, sp.Rational(1, 2), sp.Rational(1, 2), 1), sp.diag(0, 1, -1, 0)
    assert sp.trace(density) == 1 and density*hamiltonian == hamiltonian*density
    assert sp.trace(density*hamiltonian) == p
    assert sp.trace(density*charge_op**2) == 0
    assert sp.simplify(sp.trace(density*hamiltonian**2)-p**2-p*(1-p)) == 0
    count = sp.symbols('count', integer=True, nonnegative=True)
    assert sp.simplify((1-p)*count**2+p*(count+1)**2-(count+p)**2-p*(1-p)) == 0
    length, channels = math.pi+.1, 100000
    fixed = []
    for delta in (0., .1, 1., 2., math.pi-.1, 3.1):
        pl, ps = length+delta, length-delta
        prob = (1-(ps/pl)**2)/24
        assert 0 <= prob < 1/24
        common_j = (1-(2*math.pi/pl)**2)/6
        beta = -2*common_j
        h = channels*prob
        # Two admissible preparations, same mean: independent pairs or adjacent total-pair counts.
        fraction = h-math.floor(h)
        variance_independent, variance_correlated = channels*prob*(1-prob), fraction*(1-fraction)
        source_factor = 16*math.pi**2/(channels*ps**2)
        fixed.append(dict(delta=delta, pair_probability=prob, short_excitation_L0=h,
            beta_max_over_alpha=beta, minimum_compensation_j_over_alpha=16*math.pi**2*prob/ps**2,
            independent_zero_mode_j_sd_over_alpha=source_factor*math.sqrt(variance_independent),
            correlated_zero_mode_j_sd_over_alpha=source_factor*math.sqrt(variance_correlated)))
    assert fixed[-2]['beta_max_over_alpha'] == 0 and fixed[-1]['beta_max_over_alpha'] < 0
    x, w = sp.symbols('x w', nonnegative=True); z = x+w
    polynomial = sp.expand(4*z**3-(z+3)*(z+x)**2)
    assert sp.simplify(polynomial.subs(w, 1)+x*(8*x+5)) == 0
    assert sp.simplify(polynomial.subs(w, 3)-27*(x+2)) == 0
    # Match ell=4*r_e/(pi*beta) with beta_max. Solve for w=z-x to preserve the return gap.
    def matched(offset, exterior=.1):
        xx = (offset+exterior)/math.pi
        def residual(ww):
            return (3*ww**2*(ww-1)+xx*ww*(7*ww-12)+4*xx**2*(ww-3))/(1+xx)**2
        ww = brentq(residual, 1., 3., xtol=1e-14, rtol=1e-14)
        zz = xx+ww
        with mp.workdps(60):
            xxm = (mp.mpf(str(offset))+mp.mpf(str(exterior)))/mp.pi
            exact = mp.findroot(lambda v: 3*v*v*(v-1)+xxm*v*(7*v-12)+4*xxm**2*(v-3),
                                (mp.mpf('1'), mp.mpf('3')))
            root_error = abs(mp.mpf(str(ww))-exact)
            assert root_error < mp.mpf('2e-13')
        assert 1 <= ww <= 3 and zz >= 1 and zz > xx
        assert abs(1/zz-((2*zz/(zz+xx))**2-1)/3) < 1e-13
        # All lengths here are divided by ell_v; this evaluates a lower bound, not a new spacetime.
        return dict(physical_offset_over_ell_v=offset, exterior_delay_over_ell_v=exterior,
            ell_over_ell_v=zz, beta_over_alpha=1/zz,
            shorter_return_over_ell_v=math.pi*ww+2*exterior,
            offset_over_physical_loop=offset/(math.pi*zz+exterior),
            independent_matching_root_error=float(root_error),
            throat_quantum_control_q_over_16z=channels/(16*zz))
    assert matched(0., 0.)['ell_over_ell_v'] == 1
    matching = [matched(value) for value in (0., .1, 1., 10., 100., 1000., 1e6)]
    varying = []
    weights, delays = np.array([.2, .3, .5]), np.array([.1, .15, .3])
    for row in matching:
        zz, offset = row['ell_over_ell_v'], row['physical_offset_over_ell_v']
        upper = float((np.sum(weights*(2*math.pi/(math.pi+(delays+offset)/zz))**2)-1)/3)
        fast_cap = ((2*math.pi/(math.pi+(.1+offset)/zz))**2-1)/3
        assert upper <= fast_cap+1e-14
        if upper>0:
            assert math.pi*zz > .1+offset
        varying.append(dict(physical_offset_over_ell_v=offset, weighted_beta_cap_over_alpha=upper,
                            fastest_path_beta_cap_over_alpha=fast_cap))
    return dict(
        model='Stationary finite-energy density matrices in the positive-energy free NS fermion sector on the spacelike helical quotient; leading local JT equations. Not driven state preparation or a global4D solution',
        spectrum='For P_long=L+abs(Delta),P_short=L-abs(Delta)>0, f_i=2*pi*(<L0_i>-q/24)/P_i^2. Each L0_i>=0. Stationarity removes nonzero stress Fourier modes even for non-Gaussian/mixed/cross-chiral-entangled states',
        sharp_bound='For equal chiral mean stresses f, f>=-pi*q/(12*P_long^2). beta=-2*j, j=(2G/r_e^2)*(f+q/(48pi)); beta<=alpha/3*((2pi/P_long)^2-1). Some balanced stationary state has beta>0 iff P_long<2pi, with mixed states allowed',
        state='Keep long chirality in vacuum. For each short complex fermion, use (1-p)|0><0|+p|pair><pair| with p=(1-(P_short/P_long)^2)/24 and pair=b^dagger_(1/2)*d^dagger_(1/2)|0>, L0=1. Neutral, even parity, stationary, positive density matrix; finite-mode change preserves the local Hadamard property at each fixed spacelike quotient. No physical preparation protocol supplied',
        fixed_length_parameters=dict(L=length, D=.1, q=channels, support_threshold_delta=math.pi-.1,
                                     null_identification_delta=length),
        fixed_length_samples=fixed,
        fluctuation_scope='Only variance of the whole-loop chiral excitation zero mode is computed. Independent pair mixtures give Var L0=q*p*(1-p). Alternatively mix symmetric neutral states with floor(q*p) and ceil(q*p) total pairs, giving Var L0=frac(q*p)*(1-frac(q*p)). They share the same mean. Neither is the locally smeared stress noise or a metric-stability test; preparation-dependent noise cannot serve as a universal obstruction',
        physical_matching='Conditional leading symmetric mouth matching: ell_v=4*r_e/(pi*alpha), z=ell/ell_v, x=(d+S)/(pi*ell_v), S=abs(physical clock offset) in light-travel length. Optimal balanced source gives (z+3)*(z+x)^2=4*z^3. Unique z=x+w with1<=w<3 (w=1 atx0); asx grows, w=3-27/(4x)+O(x^-2)',
        causal_bound='Already beta>0 requires pi*ell>d+S, so shortest common-delay return pi*ell+d-S>2*d. The matched optimal branch gives return/ell_v=pi*w+2*d/ell_v, strictly positive. Lengthening restores mean-source support but not a past-directed loop',
        matched_common_delay_samples=matching,
        variable_channel_bound='For normalized positive weights and exterior field-line delays d_i>=d_fast>=0, total flux balance and L0 positivity give beta/alpha<=1/3*(sum w_i*(2pi/(pi+(d_i+S)/ell))^2-1)<=1/3*((2pi/(pi+(d_fast+S)/ell))^2-1). Hence beta>0 requires pi*ell>d_fast+S and fastest exterior return>2*d_fast. Does not assume equal channel delays or per-channel flux balance',
        variable_channel_parameters=dict(weights=weights.tolist(), delays_over_ell_v=delays.tolist(), fastest_exterior_over_ell_v=.1),
        variable_channel_samples=varying,
        scope='Sharp stationary mean-source result within this specified free sector and leading throat geometry, not universal chronology protection. Dynamical states/geometry, other fields, interactions, actual exteriors, physical preparation, local fluctuations and transmission remain open. The matched scale assumes the old symmetric overlap matching can be used conditionally; it is not a solved offset4D spacetime. Large-z samples are algebraic only: q/(16z)=r_e^3/(G*ell) must be much greater than1 for gravitational fluctuation control. Low-energy fermion and other MMP assumptions need separate checks',
        sources=['https://arxiv.org/pdf/hep-th/9211071 (Cappelli et al. sections3.1-3.3: NS charged-fermion spectrum and Casimir shift; used as CFT input, not a wormhole result)',
                 'https://arxiv.org/pdf/math-ph/0412028 (Fewster-Hollands section5.1: positive L0)',
                 'https://arxiv.org/pdf/1807.04726 (MMP anomaly, leading matching5.39-5.40, quantum control5.43 and unequal field-line delays5.44). Bounds, states and offset matching here are our applications'])


def _plateau_chiral_state(x,period,plateau_width,ramp_width,high,low):
    """Smooth periodic diffeomorphic NS vacuum; return V' , V'' and f/q."""
    if abs(x)>period/2: x=(x+period/2)%period-period/2
    vacuum=-math.pi/(12*period**2)
    z=(abs(x)-plateau_width/2)/ramp_width
    if z<=0: v,vp,vpp=high,0.,0.
    elif z>=1: v,vp,vpp=low,0.,0.
    else:
        logit=1/z-1/(1-z); small=math.exp(-abs(logit))
        step=small/(1+small) if logit>=0 else 1/(1+small)
        variance=small/(1+small)**2
        b=1/z**2+1/(1-z)**2
        bp=-2/z**3+2/(1-z)**3
        v=low+(high-low)*(1-step)
        vp=-(high-low)*variance*b/ramp_width*(1 if x>0 else -1)
        vpp=-(high-low)*variance*((1-2*step)*b*b+bp)/ramp_width**2
    f=vacuum*v*v-(vpp/v-1.5*(vp/v)**2)/(24*math.pi)
    return v,vp,f


def nonstationary_support_checks(balance):
    """Screen a shaped NS state against unchanged support over a finite throat strip."""
    length=balance['fixed_length_parameters']['L']
    beta=balance['fixed_length_samples'][0]['beta_max_over_alpha']
    required=-(1+3*beta)/(48*math.pi)  # f/q, in the existing throat optical chart.
    maximum_width=2*math.pi/math.sqrt(1+3*beta)
    assert math.isclose(required,-math.pi/(12*length**2),rel_tol=2e-15)
    assert math.isclose(maximum_width,length,rel_tol=2e-15)
    offset=.9*length; period=length+offset; epsilon=.9; harmonic=2
    root=math.sqrt(1-epsilon**2); frequency=2*math.pi*harmonic/period
    vacuum=-math.pi/(12*period**2)

    def state(x):
        weight=1+epsilon*math.cos(frequency*x)
        return vacuum*(harmonic**2-(harmonic**2-1)*root**2/weight**2)

    def from_covariance(x):
        g=1+epsilon*math.cos(frequency*x)
        gp=-epsilon*frequency*math.sin(frequency*x)
        gpp=-epsilon*frequency**2*math.cos(frequency*x)
        schwarzian=-gpp/g+.5*(gp/g)**2
        return vacuum*(root/g)**2-schwarzian/(24*math.pi)

    sample=np.linspace(-period/2,period/2,401)
    covariance_error=max(abs(state(x)-from_covariance(x)) for x in sample)
    assert covariance_error < 2e-14
    winding=quad(lambda x: root/(1+epsilon*math.cos(frequency*x)),0,period,
                 epsabs=1e-12,epsrel=2e-12)[0]
    mean=math.pi/(12*period**2)*((harmonic**2-1)/root-harmonic**2)
    average=quad(state,0,period,epsabs=1e-12,epsrel=2e-12)[0]/period
    weighted=quad(lambda x: (1+epsilon*math.cos(frequency*x))*state(x),
                  0,period,epsabs=1e-12,epsrel=2e-12)[0]
    weighted_bound=-math.pi/(12*period)*(harmonic**2-(harmonic**2-1)*root)
    assert abs(winding/period-1)<2e-13 and abs(average-mean)<2e-13
    assert abs(weighted-weighted_bound)<2e-13 and mean>=vacuum
    ratio=(period/length)**2
    edge_cos=(root*math.sqrt((harmonic**2-1)/(harmonic**2-ratio))-1)/epsilon
    width=period/(math.pi*harmonic)*math.acos(edge_cos)
    assert abs(state(width/2)/required-1)<2e-14
    assert state(0)<required and width<math.pi
    # A smooth plateau tests a stronger candidate than the sinusoidal family.
    plateau_width=math.pi; ramp_width=.05; high=period/length
    low=(period-high*(plateau_width+ramp_width))/(period-plateau_width-ramp_width)
    assert low>0
    def plateau(x):
        return _plateau_chiral_state(x,period,plateau_width,ramp_width,high,low)
    edges=[0.,plateau_width/2,plateau_width/2+ramp_width,period/2]
    def whole_integral(function):
        return 2*sum(quad(function,left,right,epsabs=1e-10,epsrel=2e-12)[0]
                     for left,right in zip(edges,edges[1:]))
    plateau_energy=whole_integral(lambda x: plateau(x)[2])
    by_parts=whole_integral(lambda x: vacuum*plateau(x)[0]**2+
                           (plateau(x)[1]/plateau(x)[0])**2/(48*math.pi))
    assert abs(whole_integral(lambda x: plateau(x)[0])/period-1)<2e-13
    assert abs(plateau_energy/by_parts-1)<2e-12
    assert abs(plateau(0)[2]/required-1)<2e-14
    short=length-offset; pair_probability=(1-(short/length)**2)/24
    short_stress=2*math.pi*(pair_probability-1/24)/short**2
    assert abs(short_stress/required-1)<2e-14
    plateau_samples=np.unique(np.r_[np.linspace(-period/2,period/2,241),
        np.linspace(plateau_width/2,plateau_width/2+ramp_width,121),
        np.linspace(-plateau_width/2-ramp_width,-plateau_width/2,121)])
    # The normalized cos² sampler has integral1 and integral w'^2/w=4pi²/W².
    sampler_mass=quad(lambda x: 2/length*math.cos(math.pi*x/length)**2,
                      -length/2,length/2)[0]
    derivative_integral=quad(lambda x: 8*math.pi**2/length**3*math.sin(math.pi*x/length)**2,
                             -length/2,length/2)[0]
    assert abs(sampler_mass-1)<2e-14
    assert math.isclose(derivative_integral/(48*math.pi),-required,rel_tol=2e-14)
    rows=[]
    for rho in (4.,10.,100.,1000.):
        span=2*math.atan(rho)
        duration=maximum_width-span
        assert 0<duration<span and 2*span>maximum_width
        rows.append(dict(rho_cut=rho,throat_optical_span=span,
                         maximum_unchanged_support_duration=duration,
                         null_crossing_time=span,unchanged_support_for_full_crossing=False))
    # ponytail: fixed positive-energy free NS sector and optical chart only;
    # evolving geometry and its coupled state require different equations.
    return dict(
        model='Nonstationary states of the same fixed-circle free NS CFT; positive-energy vacuum representation, c=q. This is a source feasibility screen, not formation or a dynamic spacetime',
        parameters=dict(original_optical_period=length,beta_over_alpha=beta,
                        offset=offset,offset_over_original_period=.9,long_period=period,
                        short_period=length-offset,required_chiral_stress_over_q=required),
        weighted_circle_bound='For smooth periodic w>0: integral(w*f)>=-q/(48*pi)*integral(w_prime^2/w)-pi*q/(12*integral(1/w)). Derived from circle covariance and positive L0; not a new published theorem',
        local_bound='For a nonnegative sampler supported within a proper arc: integral(w*f)>=-q/(12*pi)*integral((sqrt(w))_prime^2). For normalized cos² of width W, lower bound=-q*pi/(12*W²); use smooth approximation at endpoints',
        state=dict(epsilon=epsilon,harmonic=harmonic,
                   diffeomorphism='V_prime=sqrt(1-epsilon²)/(1+epsilon*cos(2*pi*n*x/P)); V(x+P)=V(x)+P',
                   stress='f/q=-pi/(12P²)*[n²-(n²-1)*(1-epsilon²)/(1+epsilon*cos(2*pi*n*x/P))²]',
                   peak_negative_over_required=state(0)/required,
                   qualifying_chiral_window=width,mean_stress_over_q=mean,
                   vacuum_stress_over_q=vacuum,weighted_integral_over_q=weighted,
                   weighted_lower_bound_over_q=weighted_bound,
                   excess_L0_per_flavor=(harmonic**2-1)*(1/root-1)/24),
        plateau_candidate=dict(guaranteed_local_window=plateau_width,transition_width=ramp_width,
            high_diffeomorphism_rate=high,low_diffeomorphism_rate=low,
            construction='V_prime=high on |x|<=W/2, low beyond W/2+delta; C-infinity exp(-1/z) step joins them. Its integral over P is P. On the plateau the Schwarzian vanishes and f equals the old required stress',
            energy_over_q=plateau_energy,L0_per_flavor=plateau_energy*period/(2*math.pi)+1/24,
            excess_energy_over_old_chiral_vacuum_magnitude=(plateau_energy-vacuum*period)/(math.pi/(12*length)),
            short_chiral_neutral_pair_probability=pair_probability,
            short_chiral_stress_over_q=short_stress,energy_integration_relative_error=abs(plateau_energy/by_parts-1),
            local_window_pass=True,whole_throat_for_crossing_pass=False,
            scope='Admissible mean-stress state in the prescribed free sector; large transition stresses and positive mean energy require responding geometry. No preparation protocol or survival conclusion'),
        rectangle='For chi in[-a,a] and tau in[0,T], either chiral coordinate covers width W=T+2a. Unchanged old supporting stress everywhere requires W<=L0. If W>=P_long, the whole-period mean bound already excludes it because P_long>L0',
        maximum_chiral_support_width=maximum_width,cases=rows,
        verification=dict(covariance_max_absolute_error=covariance_error,
                          winding_relative_error=abs(winding/period-1),
                          mean_absolute_error=abs(average-mean),
                          weighted_bound_absolute_error=abs(weighted-weighted_bound)),
        plot=dict(chiral_coordinate=plateau_samples.tolist(),stress_over_required_magnitude=[plateau(x)[2]/abs(required) for x in plateau_samples]),
        conclusion='Reject changing the state while retaining the original supporting stress across the entire fixed throat for a full light crossing. Local negative pulses remain allowed; the sinusoidal family is an example, not the universe of states',
        scope='No conclusion about a particular ray, which can keep one null coordinate fixed, or passage with evolving geometry. No collapse theorem, physical preparation, interacting charged-mode correction, 4D completion, seconds/calendar conversion, or new-physics novelty claim. Do not extrapolate this spacelike quotient through offset=L0',
        sources=['https://arxiv.org/pdf/math-ph/0412028 (Fewster-Hollands covariance, positive-energy representation and Theorem4.1; circle and rectangle formulas are our corollaries)',
                 'https://arxiv.org/html/1807.04726v3 (MMP state-dependent chiral stress and separate geometric anomaly)'])


def plateau_gravity_response_checks(support, parameters):
    """Local leading-JT response; changing mouth/state sewing is deliberately unresolved."""
    from scipy.optimize import brentq
    p,c=support['parameters'],support['plateau_candidate']
    beta=p['beta_over_alpha']; b0=math.pi*beta/2; half=math.pi/2
    ramp=c['transition_width']; edge=half+ramp; alpha=parameters['alpha']
    def source(u):
        if abs(u)<=half: return 0.
        f=_plateau_chiral_state(u,p['long_period'],c['guaranteed_local_window'],ramp,
                               c['high_diffeomorphism_rate'],c['low_diffeomorphism_rate'])[2]
        return 8*math.pi*(f-p['required_chiral_stress_over_q'])
    sol=solve_ivp(lambda u,y:(y[1],y[2],-y[1]-source(u)),(half,edge),(0.,0.,0.),
                  method='DOP853',rtol=2e-12,atol=2e-13,max_step=ramp/80,dense_output=True)
    assert sol.success
    ae,pe,qe=sol.y[:,-1]; constant=source(edge)
    def jet(u):
        assert -half<=u<=half+.25
        if u<=half: return np.zeros(3)
        if u<=edge: return sol.sol(u)
        d=u-edge; sn,cs=math.sin(d),math.cos(d)
        return np.array([ae+(pe+constant)*sn+qe*2*math.sin(d/2)**2-constant*d,
                         (pe+constant)*cs+qe*sn-constant,
                         -(pe+constant)*sn+qe*cs])
    thresholds=[]; grid=np.linspace(half,half+.25,501)
    for fraction in (.01,.1,.5,1.):
        f=lambda u: jet(u)[0]+fraction*b0
        brackets=[(a,b) for a,b in zip(grid,grid[1:]) if f(a)>0>=f(b)]
        assert brackets
        root=brentq(f,*brackets[0],xtol=5e-15)
        thresholds.append(dict(fractional_right_coefficient_loss=fraction,
                               formal_global_time=root-half,A_jet=jet(root).tolist()))
    stop=thresholds[-1]['formal_global_time']; errors=[]
    for u in (half+.01,half+.025,edge,half+.1,half+stop):
        intervals=[half]+([edge] if u>edge else [])+[u]
        kernels=(lambda d:2*math.sin(d/2)**2,math.sin,math.cos)
        exact=np.array([-sum(quad(lambda x:source(x)*kernel(u-x),a,b,
                     epsabs=2e-11,epsrel=2e-12)[0] for a,b in zip(intervals,intervals[1:]))
                     for kernel in kernels])
        errors.append(float(np.max(abs(jet(u)-exact)/(1+abs(exact)))))
    assert max(errors)<2e-10
    cut=100.; a=math.atan(cut); times=np.linspace(0,stop,121)
    chi_grid=np.linspace(-a,a,81); boundary=[]; min_area=math.inf; residuals=[]
    for tau in times:
        row=[]
        for chi in chi_grid:
            u=float(tau+chi); A,Ap,App=jet(u); Appp=-Ap-source(u)
            tangent=math.tan(chi); sec2=1+tangent*tangent
            phi=beta*(1+chi*tangent)+A*tangent+Ap
            gu=beta*(tangent+chi*sec2)+A*sec2
            pu=App+Ap*tangent+gu/2; pv=-gu/2
            puu=Appp+App*tangent+Ap*sec2+beta*(sec2+chi*sec2*tangent)/2+A*sec2*tangent/2
            pvv=beta*(sec2+chi*sec2*tangent)/2+A*sec2*tangent/2
            puv=-(Ap*sec2+beta*(sec2+chi*sec2*tangent)+A*sec2*tangent)/2
            scale=1+abs(puu)+abs(pvv)+abs(phi*sec2)
            residuals.append(max(abs(puu-tangent*pu-beta/2+source(u)),
                                 abs(pvv+tangent*pv-beta/2),abs(puv+phi*sec2/2))/scale)
            row.append(phi)
        min_area=min(min_area,1+alpha*min(row))
        boundary.append(dict(global_time=float(tau),
            right_formal_coefficient_ratio=1+jet(float(tau+half))[0]/b0,
            left_formal_coefficient_ratio=1.,right_fixed_cut_phi_over_alpha=row[-1],
            left_fixed_cut_phi_over_alpha=row[0]))
    assert max(residuals)<2e-13 and min_area>0
    return dict(
        model='Local leading-JT area response to the previously constructed prescribed quantum state. All three local equations hold, but exterior quantum transport and moving-mouth feedback are NOT closed',
        parameters=dict(alpha=alpha,beta_over_alpha=beta,baseline_formal_mouth_coefficient_over_alpha=b0,
                        rho_cut=cut,chi_cut=a),
        equations="Phi/alpha=beta_hat*(1+chi*tan chi)+A(u)*tan chi+A_prime(u); u=tau+chi; A'''+A'=-8pi*(f_long/q-f_req/q). Short source remains at its old value",
        initial_data='A=A_prime=A_second=0 throughout the width-pi plateau. At tau=0 the finite throat strip has exactly its old leading area and time derivative. This supplies local throat data, not a complete exterior initial slice',
        first_finite_right_boundary_source_change=half-a,
        thresholds=thresholds,formal_positive_right_coefficient_gate=stop,
        minimum_sampled_total_area_over_extremal_area=min_area,
        verification=dict(convolution_max_scaled_error=max(errors),local_JT_max_scaled_residual=max(residuals),
                          sampled_times=len(times),sampled_chi=len(chi_grid)),
        plot=boundary,
        decision='The shaped state requires moving-mouth/state feedback as soon as it changes boundary data. Holding its original incoming profile fixed formally drives b_R to zero; that is a diagnostic gate, not a predicted physical collapse time',
        scope='No feedback-free validity time follows from the old delays without exterior initial data. Do not feed this state into the static-past classical-carrier solver unchanged. No new topology, closed spacetime, horizon theorem, safe transit, seconds/calendar conversion or novelty claim',
        sources=['https://arxiv.org/html/1807.04726v3 (MMP leading throat equations and normalization)',
                 'https://arxiv.org/pdf/1402.6334 (Almheiri-Polchinski local dilaton equations; this state-specific response is our calculation)'])


def _exterior_quantum_seed(s, amplitude, delay=.1):
    """Initial excess exterior chiral flux from a compact smooth circle diffeomorphism."""
    z=s/delay
    if not 0<z<1: return 0.
    w=z*(1-z); b=math.exp(4-1/w)
    if b == 0: return 0.  # Flat-end derivatives are below floating-point range here.
    b1=b*(1-2*z)/w**2
    b2=b*(1-6*w+6*w*w)/w**4
    b3=b*(1-2*z)*(1-10*w+24*w*w-12*w**3)/w**6
    v=1+amplitude*b1
    assert v>0
    vp=amplitude*b2/delay; vpp=amplitude*b3/delay**2
    vacuum=-math.pi/(12*(math.pi+delay)**2)
    return 8*math.pi*vacuum*(2*amplitude*b1+(amplitude*b1)**2)-(vpp/v-1.5*(vp/v)**2)/3


def quantum_seed_feedback_checks():
    """Evolve an initially exterior quantum packet with the responding reduced loop."""
    delay=.1; length=math.pi+delay; final=2*length
    samples=np.unique(np.r_[np.linspace(0,final,201),np.linspace(0,delay,41),
                             np.linspace(length,length+delay,41)])
    w0=(3-math.sqrt(3))/6
    maximum_bprime=math.exp(4-1/w0)*math.sqrt(1-4*w0)/w0**2
    def squared_terms(z,amplitude):
        if not 0<z<1: return 0.,0.
        w=z*(1-z); bump=math.exp(4-1/w)
        if bump == 0: return 0.,0.
        bp=bump*(1-2*z)/w**2; bpp=bump*(1-6*w+6*w*w)/w**4
        return bp*bp,(bpp/(1+amplitude*bp))**2
    cases={}; values={}; errors=[]
    for amplitude in (0.,.001,.01,.05,.1):
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        record,vector=nonlinear_pulse_solution(0.,loops=2,samples=samples,initial_chiral=seed)
        integral1=quad(lambda z:squared_terms(z,amplitude)[0],0,1,epsabs=1e-11)[0]
        integral2=quad(lambda z:squared_terms(z,amplitude)[1],0,1,epsabs=1e-9,epsrel=2e-12)[0]
        by_parts=amplitude**2/(6*delay)*integral2-(2*math.pi/length)**2*amplitude**2*delay/6*integral1
        error=abs(record['initial_quantum_packet_energy_over_alpha']-by_parts)
        assert error<2e-11 and 1-amplitude*maximum_bprime>0 and by_parts>=0
        assert all(probe['complete_return_delay'] is None or probe['complete_return_delay']>0
                   for probe in record['launched_probes'])
        record['diffeomorphism_amplitude']=amplitude
        record['minimum_initial_diffeomorphism_derivative']=1-amplitude*maximum_bprime
        cases[str(amplitude)]=record; values[amplitude]=vector; errors.append(error)
    assert cases['0.0']['maximum_sampled_clock_rate_change']<1e-12
    refinements=[]
    for amplitude in (.01,.1):
        fine,vector=nonlinear_pulse_solution(0.,.5,loops=2,samples=samples,
                           initial_chiral=lambda s:_exterior_quantum_seed(s,amplitude,delay))
        difference=abs(vector-values[amplitude])
        error=float(np.max(difference)); scaled=float(np.max(difference/(1+abs(vector))))
        clock_error=float(np.max(difference[:,:8])); ledger_error=float(np.max(difference[:,14:]))
        # Near a small mouth-clock rate A_uu can be large. Keep absolute clock
        # and energy tests, and compare every state component with its own scale.
        assert scaled<2e-8 and clock_error<2e-8 and ledger_error<2e-8
        check=dict(amplitude=amplitude,maximum_full_state_difference=error,
                   maximum_component_scaled_difference=scaled,maximum_clock_state_difference=clock_error,
                   maximum_integral_state_difference=ledger_error,
                   refined_ledger_error=fine['maximum_reduced_energy_ledger_error'])
        if amplitude == .1:
            _,finer=nonlinear_pulse_solution(0.,.25,loops=2,samples=samples,
                              initial_chiral=lambda s:_exterior_quantum_seed(s,amplitude,delay))
            next_error=float(np.max(abs(finer-vector)/(1+abs(finer))))
            assert next_error<scaled
            check['further_refined_component_scaled_difference']=next_error
        refinements.append(check)
    return dict(
        model='Responding leading-JT mouths and retarded free-CFT state in a fixed-delay exterior. Quantum packet is initial data in the exterior, with synchronized initial mouths; no imposed time offset or continuing input/absorber',
        parameters=dict(exterior_delay=delay,original_loop_time=length,final_exterior_time=final,
                        maximum_abs_bump_derivative=maximum_bprime),
        initial_state='V(s)=s+amplitude*D*exp(4-1/[z*(1-z)]), z=s/D in(0,1), identity outside, extended by degree-one periodicity. Only one chirality is changed. f=V_prime²*f0-q*Schwarzian(V)/(24pi), with f0=-pi*q/(12L²)',
        initial_geometry='Old leading throat area and synchronized clocks; pulse is entirely in the initially exterior segment. Negative-time characteristic lookup encodes this initial slice, not a physical stationary preparation history. Exterior gravitational constraints are not solved',
        transport='For0<s<D: J_R=1/6+[Q0+delta(s)+Sch_r/3]/r_prime²; leftseed0. Evaluate delta in exterior time s. Later incoming fields come from opposite outgoing fields at s-D, including clock Jacobians and Schwarzian; initial excitation therefore recirculates',
        energy='B=-b0*(KR+KL-1); X=integral_(s-D)^s Q_out_excess + integral_s^D delta for s<D. B+X=integral_0^D delta. Initial energy is supplied in the state, not generated for free. Multiply normalized ledger byq/(8pi) for dimensionless exterior-clock energy',
        verification=dict(initial_energy_max_absolute_error=max(errors),step_refinements=refinements,
                          classical_path_baseline='Before extension, epsilon=.01/one-loop record and full state were captured; default-path rerun reproduced both exactly'),
        cases=cases,
        decision='These admissible quantum packets generate clock differences but no tested probe returns before launch. Stronger packets need not improve passage: some probes have not exited by the finite simulation end. That is not proof of permanent trapping',
        scope='Self-consistent only within this stated reduced fixed-D mean-field model. No exterior metric response, interacting charged sector, finite-mouth stress matching, quantum-noise control, preparation mechanism, topology creation, physical matter transit or calendar mapping. Clock offsets/rates are coordinates, not alone a chronology test',
        sources=['https://arxiv.org/pdf/math-ph/0412028 (CFT diffeomorphism covariance and positive energy)',
                 'https://arxiv.org/html/1807.04726v3 (leading MMP throat and quantum-state scope)',
                 'https://arxiv.org/pdf/1606.01857 (JT boundary clocks; initial-state extension and feedback test are our reduced-model calculation)'])


def causal_foliation_checks(feedback):
    """Separate the exterior-time solver guard from actual closed-return gates."""
    delay=feedback['parameters']['exterior_delay']; k=delay/math.pi
    cases={}
    for key,record in feedback['cases'].items():
        plot=record['plot']; time=np.array(plot['exterior_time'])
        clocks=np.array([time+np.array(plot[name]) for name in ('right_clock_shift','left_clock_shift')])
        # Two saved times differ by one floating-point ulp; their reconstructed
        # clock values can round equal. Nondecreasing endpoints suffice here.
        assert np.all(np.diff(time)>0) and np.all(np.diff(clocks,axis=1)>=0)
        assert min(plot['right_clock_rate']+plot['left_clock_rate'])>0
        # For monotone clocks, endpoint rectangles enclose EVERY point of each
        # interval. No interpolant or additional pulse integration is needed.
        gap_lower=np.array([math.pi+clocks[side,:-1]-clocks[1-side,1:] for side in (0,1)])
        margin=float(gap_lower.min())
        assert margin>0
        # The D-shifted closed-return gate is larger than its one-way gate.
        # Bound it directly too, splitting wherever either sampled clock jumps
        # to its next enclosing interval. Negative arguments are initial data.
        breaks=np.unique(np.r_[time,(time+delay)[time+delay<time[-1]]])
        def lower(s,side):
            if s<0: return s
            return clocks[side,max(0,np.searchsorted(time,s,side='right')-1)]
        def upper(s,side):
            if s<=0: return s
            return clocks[side,min(len(time)-1,np.searchsorted(time,s,side='left'))]
        loop_lower=min(math.pi+lower(a,side)-upper(b-delay,1-side)
                       for a,b in zip(breaks[:-1],breaks[1:]) for side in (0,1))
        assert loop_lower>0
        cases[key]=dict(one_way_clock_margin_lower_bound=margin,
                        complete_return_clock_margin_lower_bound=float(loop_lower),
                        scope='Endpoint enclosure conditional on monotone computed clocks; not a rigorous ODE error bound or a statement after the saved time interval')
    # Exact affine clock controls, not sourced gravitational solutions:
    # r(s)=s+offset, l(s)=s. The closed loop is L-offset, not pi-offset.
    examples=[]
    for offset in (0.,math.pi+.05,math.pi+delay,math.pi+delay+.05):
        one_way=math.pi-offset; complete=one_way+delay
        throat_lambda=one_way+k*math.pi
        exterior_lambda=(1+k)*delay-k*offset
        gradient_norm=-(1+k)**2+(offset/math.pi)**2
        assert math.isclose(throat_lambda,complete,abs_tol=1e-14)
        assert math.isclose(exterior_lambda,k*complete,abs_tol=1e-14)
        assert math.isclose(gradient_norm,-((math.pi+delay)**2-offset**2)/math.pi**2,abs_tol=1e-14)
        examples.append(dict(imposed_clock_offset=offset,one_way_exterior_elapsed=one_way,
                             complete_return_elapsed=complete,throat_lambda_elapsed=throat_lambda,
                             exterior_lambda_elapsed=exterior_lambda,time_gradient_norm=gradient_norm))
    assert examples[1]['one_way_exterior_elapsed']<0<examples[1]['complete_return_elapsed']
    assert examples[1]['time_gradient_norm']<0 and examples[-1]['time_gradient_norm']>0
    # An actual retarded transport recurrence in lambda passes through the
    # region rejected by an exterior-time one-leg lookup. Compare characteristic
    # tracking around the affine loop against its independent exterior-time map.
    test=examples[1]; step=test['throat_lambda_elapsed']+test['exterior_lambda_elapsed']
    cycle=delay+math.pi-test['imposed_clock_offset']
    arrival=0.; errors=[]
    for count in range(1,21):
        arrival+=test['throat_lambda_elapsed']
        arrival+=test['exterior_lambda_elapsed']
        errors.append(abs(arrival/(1+k)-count*cycle))
    assert max(errors)<2e-13 and step>0
    # Independently check the equal-lambda slice construction for changing,
    # positive-rate clocks, not just the affine counterexample.
    from scipy.optimize import brentq
    maps=(lambda s:s+.2*math.sin(s),lambda s:s-.1*math.sin(.7*s))
    rates=(lambda s:1+.2*math.cos(s),lambda s:1-.07*math.cos(.7*s))
    skew=[]; slice_rates=[]; relation_errors=[]
    for lam in np.linspace(0,12,61):
        exterior=[brentq(lambda s:s+k*f(s)-lam,-1,13,xtol=5e-14) for f in maps]
        throat=[f(s) for f,s in zip(maps,exterior)]
        ds=exterior[1]-exterior[0]; dt=throat[1]-throat[0]
        relation_errors.append(abs(ds+k*dt)); skew.append(abs(dt)/math.pi)
        slice_rates.extend(v/(1+k*v) for v in (rates[i](exterior[i]) for i in (0,1)))
        assert abs(dt)<math.pi and abs(ds)<delay
    assert max(relation_errors)<1e-13 and min(slice_rates)>0
    return dict(
        model='Causal gates and computational foliation for the existing two-mouth, pi-throat, fixed-D exterior graph; no new matter source or geometry formation',
        gates='For launch at s, R->L complete return C_R=l_inverse(r(s)+pi)+D-s. Its sign is the sign of M_R=pi+r(s)-l(s-D); swap R/L for M_L. One-way sign is pi+r(s)-l(s), a different gate',
        solver_limitation='nonlinear_pulse_solution requires outgoing emission in COMPLETED exterior-time history; its final one-way guard also requires transit>D/2. It cannot explore a negative one-way exterior delay even if the full return is positive. Successful runs do not test chronology protection',
        foliation='lambda_i(s)=s+(D/pi)*y_i(s), derivative1+(D/pi)*y_i_prime>0. On each throat ray Delta_lambda=Delta_s+D. On an exterior ray i@s -> j@(s+D), Delta_lambda=(D/pi)*[pi+y_j(s+D)-y_i(s)]. Thus both edge types are retarded if both elementary full-return gates are positive throughout the domain',
        affine_global_time='For l=s,r=s+S and throat x in[0,pi] fromLtoR: lambda=(1+k)*tau-S*x/pi. Exterior x in[0,D] fromRtoL: lambda=(1+k)*s+k*S*(1-x/D). These match at both mouths; conformal gradient norm=-(1+k)^2+(S/pi)^2<0 iff |S|<pi+D. Equality gives a closed null loop; larger |S| gives a timelike loop in this imposed affine quotient',
        dynamic_slices='At equal lambda the inverse boundary maps give ds=s_L-s_R=-(D/pi)*(l(s_L)-r(s_R)). Positive full-return gates imply abs(dtau)<pi and abs(ds)<D. Linear interpolation between equal-lambda endpoints gives spacelike slices in both conformal strips. Endpoint derivatives are positive, so slices are nested. This constructs a piecewise smooth causal foliation of the covered glued-strip domain; strict margins allow local smoothing',
        dynamic_scope='Only the ideal fixed-D glued strips and covered inverse-map domains are addressed; no4D exterior geometry or source existence theorem. The new coupled lambda feedback solver is not yet implemented; affine controls have imposed clocks rather than solved gravitational sources',
        parameters=dict(exterior_delay=delay,k=k),computed_interval_cases=cases,
        affine_examples=examples,transport_check=dict(circuits=20,maximum_exterior_arrival_error=max(errors)),
        changing_clock_slice_check=dict(samples=61,maximum_normalized_throat_skew=max(skew),
                                       minimum_throat_slice_endpoint_rate=min(slice_rates),
                                       maximum_equal_lambda_relation_error=max(relation_errors)),
        decision='Saved quantum packets remain below even the one-way gate over their computed interval. Retire exterior-clock failure as a physical obstruction; implement coupled feedback in lambda before using this code to approach a chronology threshold',
        sources=['https://arxiv.org/abs/1606.01857 (JT boundary time reparametrizations; gate and lambda construction here are our derivation)',
                 'https://arxiv.org/abs/1807.04726 (parent long wormhole has no ambient causality violation; fixed-D graph is not a4D construction)'])


def lambda_feedback_checks(feedback):
    """Verify actual coupled evolution in the causal lambda foliation."""
    cases={}; overlap=[]; states={}
    for amplitude in (0.,.01,.05):
        key=str(amplitude); saved=feedback['cases'][key]
        samples=np.array(saved['plot']['exterior_time'])
        seed=lambda s:_exterior_quantum_seed(s,amplitude)
        record,values=nonlinear_pulse_solution(0.,loops=2,samples=samples,
                                             initial_chiral=seed,lambda_time=True)
        cases[key]=record; states[key]=values
        if amplitude==0:
            exact=np.zeros_like(values); exact[:,0]=samples; exact[:,4]=samples
            exact[:,[1,5]]=1; exact[:,[3,7]]=.5
            static_error=float(np.max(abs(values-exact)))
            assert static_error<2e-11 and record['lambda_foliation']['maximum_slice_energy_ledger_error']==0
        else:
            old,old_values=nonlinear_pulse_solution(0.,loops=2,samples=samples,initial_chiral=seed)
            assert old=={k:v for k,v in saved.items() if k not in ('diffeomorphism_amplitude','minimum_initial_diffeomorphism_derivative')}
            error=abs(values-old_values)
            scaled=float(np.max(error/(1+abs(values))))
            clock=float(np.max(error[:,:8])); energy=float(np.max(error[:,14:]))
            assert scaled<2e-8 and clock<2e-8 and energy<2e-8
            probe_error=0.
            for a,b in zip(record['launched_probes'],old['launched_probes']):
                assert a['exit_reached_in_interval']==b['exit_reached_in_interval']
                if a['complete_return_delay'] is not None:
                    probe_error=max(probe_error,abs(a['complete_return_delay']-b['complete_return_delay']))
            assert probe_error<2e-8
            overlap.append(dict(amplitude=amplitude,maximum_component_scaled_difference=scaled,
                                maximum_clock_state_difference=clock,maximum_integral_state_difference=energy,
                                maximum_probe_return_difference=probe_error))
    samples=np.array(feedback['cases']['0.05']['plot']['exterior_time'])
    fine,vector=nonlinear_pulse_solution(0.,.5,loops=2,samples=samples,
                    initial_chiral=lambda s:_exterior_quantum_seed(s,.05),lambda_time=True)
    diff=abs(vector-states['0.05']); scaled=float(np.max(diff/(1+abs(vector))))
    assert scaled<2e-8 and np.max(diff[:,:8])<2e-8 and np.max(diff[:,14:])<2e-8
    return dict(
        model='Coupled leading-JT/free-CFT evolution in lambda_i=s_i+(D/pi)*y_i(s_i), with the same synchronized initial quantum packets and fixed-D exterior as increment53',
        implementation='Optional lambda_time branch reuses the existing boundary/bulk RHS and quantum transport. Physical derivatives retain their exterior-time definitions. Separate outgoing flux primitives allow charge checks on tilted slices. Both clock and exterior-time inverse maps access completed lambda history only; steps shorten when necessary',
        energy='B + F_R(s_R)-F_R(s_L-D) + F_L(s_L)-F_L(s_R-D) + R_initial(s_R) = I_initial. Multiply primitive derivatives, not their energy values, by ds_i/dlambda',
        verification=dict(static_full_state_maximum_error=static_error,default_quantum_records_exact=True,
                          default_classical_path='Captured epsilon=.01/one-loop record and full state reproduced exactly after implementation',
                          overlap=overlap,step_refinement=dict(amplitude=.05,maximum_component_scaled_difference=scaled,
                          maximum_clock_state_difference=float(np.max(diff[:,:8])),maximum_integral_state_difference=float(np.max(diff[:,14:])),
                          refined_lambda_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error'])),
        cases=cases,
        decision='The new coupled implementation agrees with the old solver where both apply. The old exterior-one-way positivity restriction is absent from the lambda branch. These tested pulses still give no backward complete return; they do not themselves enter the extra kinematic domain established by the affine control',
        scope='No evolution through a closed-return zero. Positive history requirements are solver-domain conditions, not a chronology-protection result. Finite step/resolution limits are numerical failures, never physical obstruction claims. No4D exterior response, topology creation, new supporting material, matter transit or calendar mapping',
        sources=['https://arxiv.org/abs/1606.01857 (boundary-clock framework; lambda integration and tilted ledger are our implementation)',
                 'https://arxiv.org/abs/1807.04726 (leading throat and limitations of the parent construction)'])


def wide_packet_checks(previous):
    """Compare one energy-matched wider initial state and its reversed sequence."""
    from scipy.optimize import brentq
    b0=lambda d:math.pi*((2*math.pi/(math.pi+d))**2-1)/6
    energy=lambda a,d:quad(lambda s:_exterior_quantum_seed(s,a,d),0,d,epsabs=1e-11,epsrel=1e-11)[0]
    ratio=energy(.05,.1)/b0(.1); wide_delay=1.2
    amplitude=brentq(lambda a:energy(a,wide_delay)/b0(wide_delay)-ratio,.01,.2,xtol=5e-15)
    assert 1-amplitude*4.235640422134889>0
    reflection=max(abs(_exterior_quantum_seed(wide_delay*z,-amplitude,wide_delay)
                       -_exterior_quantum_seed(wide_delay*(1-z),amplitude,wide_delay)) for z in np.linspace(.01,.99,99))
    assert reflection<1e-10 and abs(energy(amplitude,wide_delay)-energy(-amplitude,wide_delay))<1e-12
    static,vector=nonlinear_pulse_solution(0.,loops=1,initial_chiral=lambda s:0.,
                                         lambda_time=True,exterior_delay=wide_delay)
    assert static['maximum_sampled_clock_rate_change']==0 and static['lambda_foliation']['maximum_slice_energy_ledger_error']==0
    cases={}; refinements=[]
    for name,delay,a in [('narrow',.1,.05),('wide',wide_delay,amplitude),('wide_reversed',wide_delay,-amplitude)]:
        length=math.pi+delay
        samples=(np.array(previous['cases']['0.05']['plot']['exterior_time']) if name=='narrow' else
                 np.unique(np.r_[np.linspace(0,2*length,241),np.linspace(0,delay,61),np.linspace(length,length+delay,61)]))
        seed=lambda s:_exterior_quantum_seed(s,a,delay)
        refinement=1. if name=='narrow' else .5
        record,values=nonlinear_pulse_solution(0.,refinement,loops=2,samples=samples,initial_chiral=seed,
                                lambda_time=True,exterior_delay=delay,causal_diagnostics=True)
        if name=='narrow':
            assert {k:v for k,v in record.items() if k!='causal_diagnostics'}==previous['cases']['0.05']
        else:
            fine,fine_values=nonlinear_pulse_solution(0.,.25,loops=2,samples=samples,initial_chiral=seed,
                                    lambda_time=True,exterior_delay=delay,causal_diagnostics=True)
            difference=abs(fine_values-values); scaled=float(np.max(difference/(1+abs(fine_values))))
            clock=float(np.max(difference[:,:8])); integral=float(np.max(difference[:,14:]))
            assert scaled<2e-8 and clock<2e-8 and integral<2e-8
            probe_error=0.
            for p,q in zip(record['causal_diagnostics']['phase_matched_probes'],fine['causal_diagnostics']['phase_matched_probes']):
                assert p['exit_reached_in_interval']==q['exit_reached_in_interval']
                assert p['complete_return_within_interval']==q['complete_return_within_interval']
                if p['complete_return_delay'] is not None:
                    probe_error=max(probe_error,abs(p['complete_return_delay']-q['complete_return_delay']))
            assert probe_error<2e-8
            refinements.append(dict(case=name,maximum_component_scaled_difference=scaled,
                maximum_clock_state_difference=clock,maximum_integral_state_difference=integral,
                maximum_probe_return_difference=probe_error,
                refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
        probes=record['causal_diagnostics']['phase_matched_probes']
        energy_ratio=record['initial_quantum_packet_energy_over_alpha']/b0(delay)
        assert abs(energy_ratio-ratio)<2e-12
        record['packet_parameters']=dict(amplitude=a,exterior_delay=delay,initial_loop_time=length,b0=b0(delay),
                                         initial_energy_over_b0=energy_ratio,base_refinement=refinement)
        record['probe_summary']=dict(total=len(probes),exits_reached=sum(p['exit_reached_in_interval'] for p in probes),
            full_returns_within_interval=sum(p['complete_return_within_interval'] for p in probes),
            first_right_launch_return_over_loop=probes[0]['complete_return_delay']/length,
            minimum_known_return_over_loop=min(p['complete_return_delay']/length for p in probes if p['complete_return_delay'] is not None))
        cases[name]=record
    common=[]
    for probes in zip(*(cases[name]['causal_diagnostics']['phase_matched_probes'] for name in ('narrow','wide','wide_reversed'))):
        assert len({(p['launch_phase'],p['launch_side']) for p in probes})==1
        if all(p['complete_return_within_interval'] for p in probes):
            common.append(dict(launch_phase=probes[0]['launch_phase'],launch_side=probes[0]['launch_side'],
                normalized_return_delays=[p['complete_return_delay']/cases[name]['packet_parameters']['initial_loop_time']
                    for p,name in zip(probes,('narrow','wide','wide_reversed'))]))
    assert common
    return dict(
        model='Same reduced free-CFT/JT equations with consistently changed fixed exterior delay. Synchronized initial mouths; matched I/b0 across one narrow and two wide packets',
        comparison='18 explicit nonnegative launches: both sides at s/L=0,.125,...,1. Each run covers two initial loop intervals. Compare identical phases/sides; report found exits separately from full returns inside the interval. Gate magnitude alone can shrink from clock slowing',
        state_control='V_negative(s)=D-V_positive(D-s), hence delta_negative(s)=delta_positive(D-s). The wide reversed state changes stress ordering at identical geometry and total initial energy, with positive map derivative',
        verification=dict(matched_initial_energy_over_b0=ratio,wide_amplitude=amplitude,
                          maximum_reflected_seed_difference=reflection,wide_static_clock_rate_error=static['maximum_sampled_clock_rate_change'],
                          refinements=refinements,default_narrow_record_exact=True,
                          coarse_failure='The initial wide refinement1 trial exceeded the unchanged3e-9 charge tolerance (~5.42e-9); refinement.5 resolved it. Tolerances were not relaxed'),
        common_complete_probes=common,cases=cases,
        decision='Wider packets improve finite-window passage at the matched normalized initial energy. Neither order produces a tested negative complete return or negative one-way clock gate; this is improved passage, not chronology formation',
        scope='Different D changes the assumed loop geometry and background support, not just a laboratory pulse setting. No4D exterior solution/preparation, invariant device-energy estimate, permanent stability, matter transit or calendar mapping. Missing exits are finite-window censoring; minima over different completed subsets are not matched performance comparisons',
        sources=['https://arxiv.org/pdf/math-ph/0412028 (stress transformation underlying the admissible initial state; ordering and energy matching are our calculations)',
                 'https://arxiv.org/abs/1807.04726 (parent throat construction, not this changing exterior solution)'])


def long_packet_checks(previous):
    """Check continued passage and a local weak-source diagnostic for four loops."""
    delay=1.2; length=math.pi+delay; alpha=1e-9
    control,_=nonlinear_pulse_solution(0.,loops=1,initial_chiral=lambda s:0.,lambda_time=True,
                                      exterior_delay=delay,source_alpha=alpha)
    diagnostic=control['source_diagnostics']; baseline=diagnostic['baseline_center_strength']
    static_error=max(abs(np.array(diagnostic['center_chiral_strength'])-baseline))
    assert static_error<1e-22
    cases={}; verification=[]
    for name in ('wide','wide_reversed'):
        old=previous['cases'][name]; amplitude=old['packet_parameters']['amplitude']
        old_times=np.array(old['plot']['exterior_time'])
        samples=np.unique(np.r_[old_times,np.linspace(2*length,4*length,241)])
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        record,values=nonlinear_pulse_solution(0.,.5,loops=4,samples=samples,initial_chiral=seed,
            lambda_time=True,exterior_delay=delay,causal_diagnostics=True,source_alpha=alpha)
        prefix_record=record; attempts=[]
        for refinement in (.25,.125):
            fine,vector=nonlinear_pulse_solution(0.,refinement,loops=4,samples=samples,initial_chiral=seed,
                lambda_time=True,exterior_delay=delay,causal_diagnostics=True,source_alpha=alpha)
            difference=abs(vector-values); scaled=float(np.max(difference/(1+abs(vector))))
            clock=float(np.max(difference[:,:8])); integral=float(np.max(difference[:,14:]))
            assert scaled<2e-8 and clock<2e-8 and integral<2e-8
            jp=np.array(record['source_diagnostics']['center_chiral_sources_over_alpha'])
            jq=np.array(fine['source_diagnostics']['center_chiral_sources_over_alpha'])
            stress_error=float(np.max(abs(jp-jq)/(1+abs(jq))))
            strength_error=float(np.max(abs(np.array(record['source_diagnostics']['center_chiral_strength'])-
                fine['source_diagnostics']['center_chiral_strength'])))
            attempts.append(dict(base_refinement=2*refinement,fine_refinement=refinement,
                maximum_scaled_center_chiral_stress_difference=stress_error,
                maximum_absolute_center_strength_difference=strength_error,passed=stress_error<2e-7))
            if stress_error<2e-7: break
            record,values=fine,vector
        assert stress_error<2e-7, ('center stress refinement',name,stress_error,scaled,clock,integral)
        probes=record['causal_diagnostics']['phase_matched_probes']; refined=fine['causal_diagnostics']['phase_matched_probes']
        probe_error=0.
        for p,q in zip(probes,refined):
            assert p['exit_reached_in_interval']==q['exit_reached_in_interval']
            assert p['complete_return_within_interval']==q['complete_return_within_interval']
            if p['complete_return_delay'] is not None:
                probe_error=max(probe_error,abs(p['complete_return_delay']-q['complete_return_delay']))
        assert probe_error<2e-8
        prefix=np.searchsorted(samples,old_times); prefix_error=0.
        for key in ('right_clock_rate','left_clock_rate','right_clock_shift','left_clock_shift'):
            prefix_error=max(prefix_error,float(np.max(abs(np.array(prefix_record['plot'][key])[prefix]-old['plot'][key]))))
        assert prefix_error<2e-11
        strength=np.array(record['source_diagnostics']['center_chiral_strength'])
        early=strength[samples<=2*length]; later=strength[samples>2*length]
        known=[p for p in probes if p['complete_return_delay'] is not None]
        fastest=min(known,key=lambda p:p['complete_return_delay'])
        initial=[p for p in probes if p['launch_phase']<=1]
        late=[p for p in probes if p['launch_phase']>=2]
        record['long_run_summary']=dict(total_launches=len(probes),exits_reached=len(known),
            complete_returns_within_interval=sum(p['complete_return_within_interval'] for p in probes),
            initial18_complete_returns=sum(p['complete_return_within_interval'] for p in initial),
            late18_complete_returns=sum(p['complete_return_within_interval'] for p in late),
            known_returns_faster_than_initial_loop=sum(p['complete_return_delay']<length for p in known),
            fastest_probe=fastest,fastest_return_over_initial_loop=fastest['complete_return_delay']/length,
            minimum_known_one_way_exterior_elapsed=min(p['exit_exterior_time']-p['launch_exterior_time'] for p in known),
            exterior_shortcut_roundtrip_threshold=2*delay/length,
            early_center_strength_maximum=float(early.max()),later_center_strength_maximum=float(later.max()),
            later_to_early_center_strength_peak_ratio=float(later.max()/early.max()),
            minimum_sampled_center_area_ratio=min(record['source_diagnostics']['center_area_ratio']))
        verification.append(dict(case=name,maximum_component_scaled_difference=scaled,maximum_clock_state_difference=clock,
            maximum_integral_state_difference=integral,maximum_probe_return_difference=probe_error,
            maximum_scaled_center_chiral_stress_difference=stress_error,old_two_loop_clock_prefix_difference=prefix_error,
            source_refinement_attempts=attempts,accepted_base_refinement=2*refinement,
            stress_comparison_scope='Same exterior samples select nearby center events in each run; this check includes source and event-time numerical error. No exact-center interpolation or tolerance relaxation',
            refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
        cases[name]=record
    return dict(model='Four initial loops of the same two wide initial states; no source or background change from increment56',
        parameters=dict(exterior_delay=delay,initial_loop_time=length,final_exterior_time=4*length,source_alpha=alpha),
        verification=dict(static_center_strength_error=float(static_error),refinements=verification),cases=cases,
        decision='Both pulse orders retain substantial finite-window passage. Some later returns are faster than the initial unperturbed loop, but remain positive and slower than a direct exterior transit followed by exterior return. No closed time loop or ambient shortcut is established',
        stress_scope='The center observer chiral-strength check is local and conditional on alpha=1e-9. Sampled early/later maxima refer to center events selected by the corresponding exterior-data windows, not proof of global/permanent stability. Trace-sector corrections, fluctuations, finite-mouth stress and full4D geometry remain untested',
        sources=['https://arxiv.org/abs/1807.04726 (weak-source leading throat; this finite evolution is our reduced-model calculation)',
                 'https://arxiv.org/pdf/math-ph/0412028 (CFT stress transformation, not a proof of stability)'])


def throat_strip_checks(previous, packets):
    """Check the finite interior on converging grids and unchanged dynamics."""
    tau,chi,t,b=sp.symbols('tau chi tangent b',real=True); h=sp.pi/2
    u,v=tau+chi,tau-chi
    A=(u-h)**3/6; B=(v-h)**4/24
    F=b*(1+chi*t)+(A-B)*t+sp.diff(A,tau)+sp.diff(B,tau)
    dx=lambda f:sp.diff(f,chi)+(1+t*t)*sp.diff(f,t)
    du=lambda f:(sp.diff(f,tau)+dx(f))/2
    dv=lambda f:(sp.diff(f,tau)-dx(f))/2
    jp=-b/2-1-(u-h)**2/2; jm=-b/2-(v-h)-(v-h)**3/6
    # Exact chain rule d(tan chi)/dchi=1+tan(chi)^2 keeps this a polynomial check.
    residuals=[sp.expand(e) for e in (du(du(F))-t*du(F)+jp,
        dv(dv(F))+t*dv(F)+jm,du(dv(F))+F*(1+t*t)/2)]
    assert residuals==[0,0,0]
    alpha=1e-9; delay=1.2; beta=((2*math.pi/(math.pi+delay))**2-1)/3
    static=_finite_strip_diagnostics(lambda u,side:[0.,0.,0.,-beta/2],[20.,20.],alpha,beta)
    exact=alpha*beta/(1+alpha*beta); cut=math.atan(4)
    static_error=abs(static['summary']['maximum_sampled_strip_strength']-exact)
    assert static_error<1e-22
    assert abs(static['summary']['maximum_sampled_abs_phi']-alpha*beta*(1+4*cut))<1e-22
    try: _finite_strip_diagnostics(lambda u,side:[0.,0.,0.,-beta/2],[1.,1.],alpha,beta)
    except AssertionError: pass
    else: raise AssertionError('Unavailable finite-strip history was accepted')
    cases={}; checks=[]
    for name,refinement in (('wide',.5),('wide_reversed',.25)):
        old=previous['cases'][name]; amplitude=packets['cases'][name]['packet_parameters']['amplitude']
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        kw=dict(loops=4,samples=np.array(old['plot']['exterior_time']),initial_chiral=seed,
            lambda_time=True,exterior_delay=delay,causal_diagnostics=True,source_alpha=alpha,strip_diagnostics=True)
        record,values=nonlinear_pulse_solution(0.,refinement,**kw)
        baseline=dict(record); baseline.pop('strip_diagnostics')
        saved=dict(old); saved.pop('long_run_summary')
        # Historical source-screen wording was corrected; compare all numerical
        # data exactly while allowing only that documented metadata correction.
        baseline['source_diagnostics']=dict(baseline['source_diagnostics'],scope=saved['source_diagnostics']['scope'])
        assert baseline==saved, ('diagnostic changed saved dynamics',name)
        fine,vector=nonlinear_pulse_solution(0.,refinement/2,**kw)
        difference=float(np.max(abs(values-vector)/(1+abs(vector))))
        assert difference<2e-8
        aa,bb=record['strip_diagnostics'],fine['strip_diagnostics']
        direct_peak_converged=aa['sampling_converged'] and bb['sampling_converged']
        envelope_error=None
        if not direct_peak_converged:
            # The direct-peak criterion remains failed. Test a separate,
            # conservative source envelope rather than mislabel that scan.
            assert aa['source_envelope']['numerically_converged'] and bb['source_envelope']['numerically_converged']
            ea=aa['source_envelope']['resolution_checks'][-1]['chiral_strength_envelope']
            eb=bb['source_envelope']['resolution_checks'][-1]['chiral_strength_envelope']
            envelope_error=abs(ea/eb-1)
            assert envelope_error<2e-7
        assert aa['summary']['half_space_intervals']==bb['summary']['half_space_intervals']
        assert np.array_equal(aa['profile']['global_time'],bb['profile']['global_time'])
        errors={}
        for key in ('maximum_strength','center_strength','maximum_abs_phi'):
            a=np.array(aa['profile'][key]); b=np.array(bb['profile'][key])
            errors[key]=float(np.max(abs(a-b)/(alpha+abs(b))))
        assert max(errors.values())<2e-7, ('finite-strip evolution refinement',name,errors)
        checks.append(dict(case=name,refinements=[refinement,refinement/2],old_record_exact=True,
            maximum_component_scaled_difference=difference,profile_normalized_errors=errors,
            direct_peak_sampling_converged=bool(direct_peak_converged),
            source_envelope_evolution_relative_difference=envelope_error,
            refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
        cases[name]=dict(strip_diagnostics=aa,
            passage_summary=old['long_run_summary'],refined_strip_summary=bb['summary'])
    return dict(model='Same four-loop states and boundary histories; same-event finite throat reconstruction',
        parameters=dict(alpha=alpha,exterior_delay=delay,initial_loop_time=math.pi+delay),
        verification=dict(manufactured_JT_constraint_residuals=[int(x) for x in residuals],
            static_source_absolute_error=static_error,unavailable_history_rejected=True,cases=checks,
            unresolved_direct_peak_checks=[c['case'] for c in checks if not c['direct_peak_sampling_converged']]),
        cases=cases,scope='Finite-grid leading mean-field validity test, not permanent stability, near-mouth completion, material preparation or chronology change',
        sources=['https://arxiv.org/html/1807.04726v3#S5.SS3.SSS2 (leading spherical throat and matter normalization; reconstruction here is our calculation)'])


def energy_sweep_checks(previous):
    """Test increased initial energy at fixed geometry and matched launches."""
    delay=1.2; length=math.pi+delay; alpha=1e-9
    b0=math.pi*((2*math.pi/length)**2-1)/6
    samples=np.array(previous['cases']['wide']['plot']['exterior_time'])
    rows=[]; cases={}; checks=[]
    def summarize(record,ratio,orientation,values=None):
        probes=record['causal_diagnostics']['phase_matched_probes']
        known=[p for p in probes if p['complete_return_delay'] is not None]
        common=probes[0]
        assert common['launch_phase']==0 and common['launch_side']==0 and common['complete_return_within_interval']
        strength=np.array(record['source_diagnostics']['center_chiral_strength'])
        times=np.array(record['plot']['exterior_time']); crossings=[]
        for threshold in (.001,.01,.1,1.):
            indices=np.flatnonzero(strength>=threshold)
            if len(indices):
                k=int(indices[0]); crossings.append(dict(threshold=threshold,
                    first_resolved_exterior_time=float(times[k]),exterior_phase=float(times[k]/length),
                    center_global_time=record['source_diagnostics']['center_global_time'][k],strength=float(strength[k])))
        order_one=next((p['exterior_phase'] for p in crossings if p['threshold']==1),None)
        row=dict(energy_ratio=ratio,orientation=orientation,total_launches=len(probes),exits_reached=len(known),
            complete_returns_within_interval=sum(p['complete_return_within_interval'] for p in probes),
            common_launch_return_over_loop=common['complete_return_delay']/length,
            fastest_resolved_return_over_loop=min(p['complete_return_delay'] for p in known)/length,
            minimum_sampled_normalized_complete_clock_gate=record['causal_diagnostics']['minimum_sampled_normalized_complete_clock_gate'],
            minimum_sampled_clock_rate=record['minimum_sampled_clock_rate'],
            maximum_sampled_center_strength=float(strength.max()),first_sampled_order_one_exterior_phase=order_one,
            first_resolved_source_crossings=crossings,
            common_return_precedes_first_resolved_order_one=bool(order_one is None or common['complete_return_delay']/length<order_one),
            late_source_screen='assumed weak-source screen fails' if order_one is not None else 'no sampled order-one center source',
            minimum_sampled_boundary_charge_sum=None if values is None else float(np.min(values[:,3]+values[:,7])))
        rows.append(row)
    for name,orientation in (('wide','wide'),('wide_reversed','reversed')):
        record=previous['cases'][name]
        summarize(record,record['initial_quantum_packet_energy_over_alpha']/b0,orientation)
    for ratio,amplitude in ((2.,.12412962953797861),(3.,.14360153622486665),(4.,.15752457643370363)):
        assert 1-amplitude*4.235640422134889>0
        for sign,orientation in ((1,'wide'),(-1,'reversed')):
            seed=lambda s:_exterior_quantum_seed(s,sign*amplitude,delay)
            energy=quad(seed,0,delay,epsabs=1e-12,epsrel=2e-12)[0]
            assert abs(energy/b0-ratio)<1e-12
            kw=dict(loops=4,samples=samples,initial_chiral=seed,lambda_time=True,exterior_delay=delay,
                causal_diagnostics=True,source_alpha=alpha)
            record,values=nonlinear_pulse_solution(0.,.25,**kw)
            fine,vector=nonlinear_pulse_solution(0.,.125,**kw)
            state_error=float(np.max(abs(values-vector)/(1+abs(vector))))
            a=np.array(record['source_diagnostics']['center_chiral_strength']);b=np.array(fine['source_diagnostics']['center_chiral_strength'])
            source_error=float(np.max(abs(a-b)/(alpha+abs(b))))
            assert state_error<2e-8 and source_error<2e-7, ('energy sweep refinement',ratio,orientation,state_error,source_error)
            probe_error=0.
            for p,q in zip(record['causal_diagnostics']['phase_matched_probes'],fine['causal_diagnostics']['phase_matched_probes']):
                assert p['exit_reached_in_interval']==q['exit_reached_in_interval']
                assert p['complete_return_within_interval']==q['complete_return_within_interval']
                if p['complete_return_delay'] is not None: probe_error=max(probe_error,abs(p['complete_return_delay']-q['complete_return_delay']))
            assert probe_error<2e-8
            for threshold in (.001,.01,.1,1.):
                assert np.array_equal(a>=threshold,b>=threshold), ('source crossing sampling changes',ratio,orientation,threshold)
            key=f'{ratio:g}_{orientation}'
            record['source_diagnostics']['scope']=WEAK_SOURCE_SCOPE
            record['packet_parameters']=dict(energy_ratio=ratio,amplitude=sign*amplitude,
                minimum_diffeomorphism_rate=1-amplitude*4.235640422134889)
            cases[key]=record; summarize(record,ratio,orientation,values)
            checks.append(dict(case=key,maximum_component_scaled_difference=state_error,
                maximum_normalized_center_strength_difference=source_error,maximum_return_time_difference=probe_error,
                all_exit_and_return_statuses_match=True,source_threshold_sample_masks_match=True,
                refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
    return dict(parameters=dict(exterior_delay=delay,initial_loop_time=length,final_exterior_time=4*length,alpha=alpha),
        rows=rows,cases=cases,verification=checks,
        decision='In this tested family, added initial energy suppresses finite-window returns and lengthens the same initial right launch. Smaller clock gates do not imply faster physical return. At ratio4 the assumed weak-source screen fails; late outputs are formal reduced-model continuations, not certified physical predictions',
        scope='No universal trapping threshold from I/b0 or the boundary-charge sum. Signed exterior excess energy must not be clipped. Missing exits are finite-window censoring; S>=1 is a warning against the adopted small-source assumption, not a horizon, collapse, or universal JT-validity theorem. No source preparation, complete4D spacetime or backward calendar transit',
        sources=['https://arxiv.org/html/1807.04726v3#S5.SS4 (parent construction discusses energy and stability; not a published threshold for this driven state)'],
        interpretation_update='S>=1 is an observer-specific weak-source warning, not a proved universal JT breakdown. This supersedes stronger historical wording; numerical history is unchanged.')


def center_residual_checks(previous):
    """Check the reconstructed spherical geometry, without changing its dynamics."""
    t,x,alpha,beta=sp.symbols('t x alpha beta', real=True)
    aj=sp.symbols('A0:4'); bj=sp.symbols('B0:4')
    aa=sum(aj[k]*(t+x)**k/sp.factorial(k) for k in range(4))
    bb=sum(bj[k]*(t-x)**k/sp.factorial(k) for k in range(4))
    phi=alpha*(beta*(1+x*sp.tan(x))+(aa-bb)*sp.tan(x)+sp.diff(aa,t)+sp.diff(bb,t))
    du=lambda f:(sp.diff(f,t)+sp.diff(f,x))/2
    dv=lambda f:(sp.diff(f,t)-sp.diff(f,x))/2
    for derivative,jets,sign in ((du,aj,1),(dv,bj,-1)):
        assert sp.expand(derivative(phi).subs({t:0,x:0})-alpha*(jets[2]+sign*(aj[0]-bj[0])/2))==0
        assert sp.expand(derivative(derivative(phi)).subs({t:0,x:0})-alpha*(jets[1]+jets[3]+beta/2))==0
    f,fu,fuu,j=sp.symbols('F Fu Fuu J', positive=True)
    radius=sp.sqrt(f)
    ru=sp.diff(radius,f)*fu
    ruu=sp.diff(ru,f)*fu+sp.diff(ru,fu)*fuu
    exact=-2*ruu/radius-alpha*j/f
    assert sp.simplify(exact.subs(fuu,-alpha*j)-fu**2/(2*f**2))==0
    delay=previous['parameters']['exterior_delay']; length=math.pi+delay
    alpha=previous['parameters']['alpha']; cases={}; rows=[]; checks=[]
    for key in ('3_wide','3_reversed','4_wide','4_reversed'):
        old=previous['cases'][key]; amplitude=old['packet_parameters']['amplitude']
        samples=np.array(old['plot']['exterior_time'])
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        kw=dict(loops=4,samples=samples,initial_chiral=seed,lambda_time=True,
                exterior_delay=delay,source_alpha=alpha,causal_diagnostics=True,null_residual=True)
        record,values=nonlinear_pulse_solution(0.,.25,**kw)
        fine,vector=nonlinear_pulse_solution(0.,.125,**kw)
        legacy=dict(record); legacy.pop('null_residual')
        expected=dict(old); expected.pop('packet_parameters')
        assert legacy==expected, ('new diagnostic changed legacy output',key)
        assert record['null_residual']['center_phi_u_v'][0]==[0.,0.]
        errors={}
        for name in ('center_phi','center_phi_u_v','absolute_observer_residual','background_plus_source_normalized_residual'):
            a=np.array(record['null_residual'][name]); b=np.array(fine['null_residual'][name])
            errors[name]=float(np.max(abs(a-b))/max(1e-30,float(np.max(abs(b)))))
            assert errors[name]<1e-6, ('center residual refinement',key,name,errors[name])
        source=record['source_diagnostics']; diagnostic=record['null_residual']
        absolute=np.array(diagnostic['absolute_observer_residual'])
        normalized=np.array(diagnostic['background_plus_source_normalized_residual'])
        gradients=np.array(diagnostic['center_phi_u_v']); area=np.array(source['center_area_ratio'])[:,None]
        assert np.all(absolute>=0) and np.all(normalized<=absolute)
        assert np.allclose(absolute,2*gradients**2/area**2,rtol=1e-14,atol=0)
        peak=int(np.argmax(np.max(absolute,axis=1)))
        rows.append(dict(case=key,energy_ratio=old['packet_parameters']['energy_ratio'],
            maximum_sampled_center_strength=max(source['center_chiral_strength']),
            maximum_absolute_radial_null_residual=float(absolute.max()),
            maximum_background_plus_source_normalized_residual=float(normalized.max()),
            maximum_absolute_phi=float(np.max(abs(np.array(diagnostic['center_phi'])))),
            absolute_residual_peak_exterior_phase=float(samples[peak]/length)))
        cases[key]=dict(exterior_phase=(samples/length).tolist(),center_global_time=source['center_global_time'],
            center_chiral_strength=source['center_chiral_strength'],**diagnostic)
        checks.append(dict(case=key,legacy_record_identical=True,relative_linf_refinement_errors=errors,
            maximum_component_scaled_state_difference=float(np.max(abs(values-vector)/(1+abs(vector)))),
            refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
    result=dict(parameters=dict(alpha=alpha,exterior_delay=delay,initial_loop_time=length,final_exterior_phase=4),
        rows=rows,cases=cases,verification=checks,
        algebra_checks=dict(center_derivatives=True,center_leading_constraints=True,exact_spherical_null_residual=True),
        scope='Sampled center events only. Exact residual of a specified reconstructed ansatz, not an evolved4D solution, global stability proof or physical source preparation. Passing radial-null equations does not test mixed/angular closure or state bandwidth.',
        sources=['https://arxiv.org/html/1807.04726v3#S5.SS3.SSS2 (published leading throat foundation; residual analysis here is our calculation)'])
    return center_polyakov_residuals(result)


def center_polyakov_residuals(result):
    """Evaluate the already specified radius-independent Polyakov closure."""
    p,p0,a,g=sp.symbols('p p0 a g', real=True); f=1+p; f0=1+p0
    cross=-2*(3*p*p+a*f)/(3*f*f); angular=(p*(3+2*p)+g)/(f*f)
    radius=sp.sqrt(f); ruv=-p/(4*radius)-g/(4*radius**3); ru_rv=g/(4*f)
    guv=2*ruv/radius+2*ru_rv/f+1/(2*f)
    assert sp.simplify(4*(guv-1/(2*f*f)-a/(6*f))-cross)==0
    assert sp.simplify(-4*ruv/radius+1-1/f**2-angular)==0
    dx=2*(p-p0)*(a*f*f0/3-p-p0-2*p*p0)/(f*f*f0*f0)
    da=(p-p0)*(f*f0+f+f0)/(f*f*f0*f0)+g/(f*f)
    assert sp.cancel(cross-cross.subs(p,p0)-dx)==0
    assert sp.cancel(angular-angular.subs({p:p0,g:0})-da)==0
    alpha=result['parameters']['alpha']; length=result['parameters']['initial_loop_time']
    phi0=alpha*(((2*math.pi/length)**2-1)/3); area0=1+phi0
    cross0=-2*(3*phi0**2+alpha*area0)/(3*area0**2)
    angular0=phi0*(3+2*phi0)/area0**2
    for row in result['rows']:
        case=result['cases'][row['case']]; phi=np.array(case['center_phi']); area=1+phi
        gradients=np.array(case['center_phi_u_v']); product=gradients[:,0]*gradients[:,1]
        cross=-2*(3*phi**2+alpha*area)/(3*area**2)
        angular=(phi*(3+2*phi)+product)/area**2
        cross_change=2*(phi-phi0)*(alpha*area*area0/3-phi-phi0-2*phi*phi0)/(area**2*area0**2)
        angular_change=(phi-phi0)*(area*area0+area+area0)/(area**2*area0**2)+product/area**2
        assert np.all(np.isfinite(cross+angular+cross_change+angular_change))
        assert abs(cross_change[0])<1e-25 and abs(angular_change[0])<1e-25
        case['conditional_polyakov']=dict(mixed_observer_residual=cross.tolist(),angular_residual=angular.tolist(),
            mixed_pulse_minus_static=cross_change.tolist(),angular_pulse_minus_static=angular_change.tolist())
        row.update(maximum_absolute_mixed_residual=float(np.max(abs(cross))),
            maximum_absolute_angular_residual=float(np.max(abs(angular))),
            maximum_absolute_mixed_pulse_change=float(np.max(abs(cross_change))),
            maximum_absolute_angular_pulse_change=float(np.max(abs(angular_change))))
    result['conditional_polyakov']=dict(model='Existing radius-independent Polyakov action plus fixed-charge spherical Maxwell source; p_perp=0 follows from this specified action',
        static_center_phi=phi0,static_mixed_residual=cross0,static_angular_residual=angular0,
        definition='Ecross=re^2*E(kplus,kminus)=-2*(3phi^2+alpha*F)/(3F^2); Eangular=re^2*E^theta_theta=[phi*(3+2phi)+phi_u*phi_v]/F^2 at center; E=G-8piG*Ttotal',
        calibration='The O(alpha) static mixed defect is the known trace/reference mismatch of the uncorrected leading metric. Report absolute and pulse-minus-static results separately; not a pulse instability',
        scope='Conditional mean-field residuals on the reconstructed metric at sampled center events. Not full4D renormalized fermion stress, an evolved corrected metric, off-center control, stability, state bandwidth, preparation or backward transit',
        verification='Mixed/angular difference identities checked symbolically; evaluated from the already step-refined phi and null gradients. No additional dynamics or independent mixed/angular mesh refinement')
    result['algebra_checks']['conditional_mixed_angular_difference_identities']=True
    result['algebra_checks']['conditional_mixed_angular_warped_geometry']=True
    result['scope']='Sampled center equations only. Radial-null residuals plus mixed/angular residuals for the specified Polyakov closure are evaluated, not solved away. No full4D fermion closure, global stability or physical source preparation.'
    if 'https://arxiv.org/pdf/hep-th/0006035' not in result['sources']:
        result['sources'].append('https://arxiv.org/pdf/hep-th/0006035')
    return result


def spatial_residual_checks(previous):
    """Check both strong-state orders over their full available finite interiors."""
    tau,chi,t,b=sp.symbols('tau chi tangent b',real=True)
    u=tau+chi; v=tau-chi
    A=u**3/6; B=v**4/24
    phi=b*(1+chi*t)+(A-B)*t+sp.diff(A,tau)+sp.diff(B,tau)
    dx=lambda f:sp.diff(f,chi)+(1+t*t)*sp.diff(f,t)
    du=lambda f:(sp.diff(f,tau)+dx(f))/2
    dv=lambda f:(sp.diff(f,tau)-dx(f))/2
    offset=(b*(t+chi*(1+t*t))+(A-B)*(1+t*t))/2
    assert sp.expand(du(phi)-sp.diff(A,tau,2)-sp.diff(A,tau)*t-offset)==0
    assert sp.expand(dv(phi)-sp.diff(B,tau,2)+sp.diff(B,tau)*t+offset)==0
    assert sp.expand(du(dv(phi))+phi*(1+t*t)/2)==0
    parameters=previous['parameters']; alpha=parameters['alpha']; delay=parameters['exterior_delay']
    length=math.pi+delay; beta=((2*math.pi/length)**2-1)/3
    # Manufactured static values test the same numeric evaluator independently.
    x=np.array([0.,.2,math.atan(4.)]); zeros=np.zeros((3,3))
    p,values=_spherical_residual_values(x,zeros,zeros,alpha,beta)
    p0=alpha*beta*(1+x*np.tan(x)); g0=alpha*beta*(np.tan(x)+x/np.cos(x)**2)/2
    assert np.allclose(p,p0,rtol=1e-14,atol=0)
    assert np.allclose(values[:,0],2*(np.cos(x)*g0/(1+p0))**2,rtol=1e-14,atol=0)
    assert np.allclose(values[:,0],values[:,1],rtol=1e-14,atol=0)
    cases={}; checks=[]
    for key in ('4_wide','4_reversed'):
        old=previous['cases'][key]; amplitude=old['packet_parameters']['amplitude']
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        kw=dict(loops=4,samples=np.array(old['plot']['exterior_time']),initial_chiral=seed,
            lambda_time=True,exterior_delay=delay,causal_diagnostics=True,source_alpha=alpha,
            null_residual=True,residual_envelope=True)
        record,values=nonlinear_pulse_solution(0.,.25,**kw)
        fine,vector=nonlinear_pulse_solution(0.,.125,**kw)
        legacy=dict(record); legacy.pop('residual_envelope'); legacy.pop('null_residual')
        expected=dict(old); expected.pop('packet_parameters')
        assert legacy==expected, ('spatial diagnostic changed legacy output',key)
        a,b=record['residual_envelope'],fine['residual_envelope']
        bounds=np.array(a['summary']['absolute_residual_envelope'])
        refined_bounds=np.array(b['summary']['absolute_residual_envelope'])
        state_error=float(np.max(abs(values-vector)/(1+abs(vector))))
        envelope_error=float(np.max(abs(bounds/refined_bounds-1)))
        assert state_error<2e-8 and envelope_error<1e-6
        assert a['numerical_envelope_converged'] and b['numerical_envelope_converged']
        center=np.array(record['null_residual']['absolute_observer_residual'])
        assert np.all(np.max(center,axis=0)<=bounds[:2])
        assert max(record['source_diagnostics']['center_global_time'])<=a['latest_center_global_time']+1e-13
        cases[key]=a
        checks.append(dict(case=key,legacy_record_identical=True,maximum_scaled_state_difference=state_error,
            maximum_relative_envelope_evolution_difference=envelope_error,
            coarse_and_fine_envelope_search_converged=True,
            sampled_peaks_converged_on_both_evolutions=bool(a['sampled_peaks_converged'] and b['sampled_peaks_converged']),
            exact_center_samples_below_envelope=True,all_center_times_covered=True,
            refined_absolute_residual_envelope=refined_bounds.tolist(),
            refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
    return dict(parameters=dict(alpha=alpha,exterior_delay=delay,initial_loop_time=length),cases=cases,
        verification=dict(manufactured_derivative_and_mixed_constraint_identities=True,static_numeric_evaluator=True,
            cases=checks,unresolved_sampled_peak_checks=[c['case'] for c in checks if not c['sampled_peaks_converged_on_both_evolutions']]),
        scope='Finite interior residual bounds conditional on numerically estimated six-jet maxima and specified Polyakov matter. Not certified continuum enclosures, an evolved corrected metric, near-mouth/exterior matching, stability or backward travel.',
        sources=['https://arxiv.org/html/1807.04726v3#A2','https://arxiv.org/pdf/hep-th/0006035'])


def timing_response_checks(previous):
    """Exact inverse-clock sensitivity; no assumption that residuals bound phase."""
    delay=previous['parameters']['exterior_delay']; alpha=previous['parameters']['alpha']
    cases={}; checks=[]
    for key in ('4_wide','4_reversed'):
        old=previous['cases'][key]; amplitude=old['packet_parameters']['amplitude']
        seed=lambda s:_exterior_quantum_seed(s,amplitude,delay)
        kw=dict(loops=4,samples=np.array(old['plot']['exterior_time']),initial_chiral=seed,
            lambda_time=True,exterior_delay=delay,causal_diagnostics=True,source_alpha=alpha,timing_response=True)
        record,values=nonlinear_pulse_solution(0.,.25,**kw)
        fine,vector=nonlinear_pulse_solution(0.,.125,**kw)
        legacy=dict(record);legacy.pop('timing_response');expected=dict(old);expected.pop('packet_parameters')
        assert legacy==expected
        a,b=record['timing_response'],fine['timing_response']; rate_error=0.; interval_error=0.; derivative_error=0.
        for p,q in zip(a['probes'],b['probes']):
            assert p['exit_reached_in_interval']==q['exit_reached_in_interval']
            assert p['return_margin_uses_initial_history']==q['return_margin_uses_initial_history']
            assert abs(p['endpoint_exit_margin']-q['endpoint_exit_margin'])<2e-8
            assert abs(p['complete_return_phase_margin']-q['complete_return_phase_margin'])<2e-8
            if not p['exit_reached_in_interval']: continue
            rate_error=max(rate_error,abs(p['linear_time_gain_per_phase']/q['linear_time_gain_per_phase']-1))
            assert p['phase_to_advance_one_percent_loop']>0 and p['phase_to_delay_one_percent_loop']>0
            for i,j in zip(p['conditional_response_intervals'],q['conditional_response_intervals']):
                assert i['assumed_combined_phase_error']==j['assumed_combined_phase_error']
                assert all(t is not None for t in i['exit_interval']+j['exit_interval'])
                lower,upper=i['exit_interval']; assert lower<p['exit_exterior_time']<upper
                interval_error=max(interval_error,max(abs(np.array(i['exit_interval'])-j['exit_interval'])))
                if i['assumed_combined_phase_error']==1e-6:
                    numerical=(upper-lower)/(2e-6)
                    derivative_error=max(derivative_error,abs(numerical/p['linear_time_gain_per_phase']-1))
        assert rate_error<2e-7 and interval_error<2e-8 and derivative_error<1e-6
        state_error=float(np.max(abs(values-vector)/(1+abs(vector))));assert state_error<2e-8
        cases[key]=a
        checks.append(dict(case=key,legacy_record_identical=True,maximum_relative_gain_refinement=rate_error,
            maximum_inverse_endpoint_refinement=interval_error,maximum_centered_derivative_relative_error=derivative_error,
            maximum_scaled_state_difference=state_error,
            refined_slice_ledger_error=fine['lambda_foliation']['maximum_slice_energy_ledger_error']))
    return dict(parameters=dict(alpha=alpha,exterior_delay=delay,initial_loop_time=math.pi+delay),
        cases=cases,verification=checks,
        scope='Conditional inverse-clock enclosures and exact launch-specific margins, not an estimate of actual metric/source phase errors, completed4D matching or chronology change')


def _magnetic_return_path(nu, s):
    # Prolate coordinates, mouth centers z/d=+-1/2. Rationalize u-1.
    h = 1-s*s
    root = math.hypot(h, nu*s)
    u = (h+root)/nu
    um = (2-nu)*h/(nu*(u+1)-2*h)
    rho = math.sqrt(um*(u+1)*h)/2
    up = s*(nu-2*u)/root
    rp = (u*up*h-s*um*(u+1))/(4*rho)
    zp = (u+s*up)/2
    return rho, u*s/2, math.hypot(rp, zp)


def _magnetic_return_length(nu, endpoint=1.):
    assert 0 < nu <= 2 and 0 < endpoint <= 1
    if nu == 2:
        return endpoint
    if nu < 1e-4:
        # Integrate g=nu*ds_path/ds. The original Cartesian derivatives
        # lose accuracy inside an endpoint layer of width O(nu).
        def scaled_speed(s):
            h=(1-s)*(1+s); root=math.hypot(h,nu*s); v=h+root
            # v=nu*u; f=2*v-nu**2, rationalized also near the axis nu=2.
            f=(2-nu)*(2+nu)*v/(v+2*s*s)
            return math.sqrt(v/2)*math.hypot(1.,s*math.sqrt(f)/root)
        points=sorted({s for s in (1-math.sqrt(nu),1-nu) if 0<s<endpoint})
        return 2*quad(scaled_speed,0,endpoint,points=points,
                      epsabs=1e-10*nu,epsrel=2e-12)[0]/nu
    return 2*quad(lambda s: _magnetic_return_path(nu, s)[2], 0, endpoint,
                  epsabs=1e-10, epsrel=2e-12)[0]


def _magnetic_return_length_precision_checks():
    """Independent Cartesian-derivative reference for the small-nu integral."""
    import mpmath as mp
    from polyakov_response import _magnetic_closure_paths
    failed_node=float(np.polynomial.legendre.leggauss(512)[0][0]+1)
    records=[]
    with mp.workdps(60):
        for nu in (failed_node,1e-6,8e-5):
            n=mp.mpf(nu)
            def speed(s):
                if s==1: return 1/n
                h=1-s*s; root=mp.sqrt(h*h+n*n*s*s); u=(h+root)/n
                rho=mp.sqrt((u*u-1)*h)/2; up=s*(n-2*u)/root
                rp=(u*up*h-s*(u*u-1))/(4*rho); zp=(u+s*up)/2
                return mp.sqrt(rp*rp+zp*zp)
            reference=2*mp.quad(speed,[0,mp.mpf('.5'),1-mp.sqrt(n),1-n,1-n/10,1])
            actual=_magnetic_return_length(nu)
            error=float(abs(mp.mpf(actual)-reference)/reference)
            assert error<2e-12
            records.append(dict(nu=nu,length=actual,reference_60_digits=mp.nstr(reference,60),
                relative_error=error))
    overlap=2*quad(lambda s:_magnetic_return_path(8e-5,s)[2],0,1,
                   epsabs=1e-10,epsrel=2e-12)[0]
    overlap_error=abs(records[-1]['length']-overlap)/records[-1]['length']
    assert overlap_error<2e-12
    saved=json.loads(Path(__file__).with_name('polyakov_closure_results.json').read_text(encoding='utf-8'))
    row=next(r for r in saved['quadrature_refinement'] if r['order']==128)['roots'][-1]
    distance=saved['parameters']['separation_over_re']
    _,_,lengths,pieces,_=_magnetic_closure_paths(128,distance,(100.,))
    piece=pieces['100.0']
    periods=row['throat_optical_time']+(distance*lengths-2*piece['isotropic']+piece['shapiro'])/row['clock_scale']
    preserved=bool(np.array_equal(periods,np.asarray(row['profile']['periods'])))
    assert preserved
    return dict(all_declared_checks_pass=True,relative_tolerance=2e-12,
        reference='60-digit integration of original Cartesian path derivatives; no scaled-speed formula reused',
        cases=records,overlap_relative_difference=overlap_error,
        legacy128_saved_periods_bitwise_equal=preserved,
        production_tolerances='epsrel=2e-12; epsabs=1e-10 for original speed, nu*1e-10 for scaled speed')


def magnetic_return_path_checks():
    """Finite flat-region flux paths and the leading MMP variational closure."""
    from scipy.optimize import brentq
    path, length = _magnetic_return_path, _magnetic_return_length

    field_lines, flux_errors = [], []
    for nu in (.01, .1, .5, 1., 1.5, 1.9, 2.):
        f = length(nu)
        assert f >= 2/nu-1e-11  # Broken straight route through the midplane.
        if nu < 2:
            ss = np.linspace(-.999, .999, 301)
            coordinates = [path(nu, float(s))[:2] for s in ss]
            for rho, z in coordinates:
                c1 = (z+.5)/math.hypot(rho, z+.5)
                c2 = (z-.5)/math.hypot(rho, z-.5)
                flux_errors.append(abs(c1-c2-nu))
            polygon = np.array([(0., -.5)]+coordinates+[(0., .5)])
            polygon_length = np.linalg.norm(np.diff(polygon, axis=0), axis=1).sum()
            assert polygon_length < f
            assert 1-polygon_length/f < 2e-4
        else:
            polygon = np.array([(0., -.5), (0., .5)])
        field_lines.append(dict(nu=nu, length_over_separation=f,
            rho_over_d=polygon[:, 0].tolist(), z_over_d=polygon[:, 1].tolist()))
    assert max(flux_errors) < 3e-12

    def quadrature(n):
        nodes, weights = np.polynomial.legendre.leggauss(n)
        return nodes+1, weights/2, np.array([length(float(v+1)) for v in nodes])

    nu, weights, f = quadrature(256)
    nu_coarse, weights_coarse, f_coarse = quadrature(128)
    assert np.all(f > 2/nu) and abs(sum(weights)-1) < 1e-14
    def omega2(delay, lengths=f, ww=weights):
        return float(ww@((2*math.pi/(math.pi+delay*lengths))**2))

    fixed, refinement = [], []
    for delay in (.1, 1., 10.):
        value = omega2(delay)
        error = abs(value-omega2(delay, f_coarse, weights_coarse))
        refinement.append(error)
        fastest = (2*math.pi/(math.pi+delay))**2
        assert 0 < value < fastest
        local_beta = ((2*math.pi/(math.pi+delay*f))**2-1)/3
        mean_beta = (value-1)/3
        rms = float(math.sqrt(weights@((local_beta-mean_beta)**2)))
        critical_nu = 2. if delay >= math.pi else brentq(
            lambda v: delay*length(v)-math.pi, 1e-7, 2., xtol=1e-13)
        fixed.append(dict(separation_over_ell=delay, mean_omega_squared=value,
            common_shortest_path_omega_squared=fastest,
            mean_beta_over_alpha=mean_beta,
            common_shortest_path_beta_over_alpha=(fastest-1)/3,
            angular_rms_beta_over_alpha=rms,
            angular_dipole_legendre_coefficient=float(3*weights@(local_beta*(nu-1))),
            nonnegative_traceless_null_source_flux_fraction=critical_nu/2,
            quadrature_128_256_difference=error))

    def matched(delta, lengths=f, ww=weights):
        def residual(z):
            return z*(omega2(delta/z, lengths, ww)-1)-3
        high = max(2., 2*delta)
        while residual(high) <= 0:
            high *= 2
        z = brentq(residual, 1., high, xtol=1e-13, rtol=1e-13)
        periods = math.pi*z+delta*lengths
        energy = 1/(16*z*z)+1/(24*z)-math.pi/6*float(ww@(1/periods))
        slope = -1/(8*z**3)-1/(24*z*z)+math.pi**2/6*float(ww@(1/periods**2))
        curvature = 3/(8*z**4)+1/(12*z**3)-math.pi**3/3*float(ww@(1/periods**3))
        assert z >= 1 and math.pi*z > delta and curvature > 0
        assert abs(residual(z)) < 2e-11 and abs(slope)*z*z < 1e-13
        return dict(separation_over_ell_v=delta, ell_over_ell_v=z,
            throat_time_over_fastest_exterior=math.pi*z/delta,
            mean_beta_over_alpha=1/z, energy_times_ell_v_over_q=energy,
            energy_first_derivative=slope, energy_second_derivative=curvature,
            matching_residual=residual(z))

    matching = [matched(delta) for delta in (.1, 1., 10., 100.)]
    for row in matching:
        coarse = matched(row['separation_over_ell_v'], f_coarse, weights_coarse)
        row['root_quadrature_128_256_difference'] = abs(
            row['ell_over_ell_v']-coarse['ell_over_ell_v'])
        refinement.append(row['root_quadrature_128_256_difference'])
    asymptotic_y = brentq(lambda y: float(weights@((2*y/(y+f))**2))-1,
                         1., 5., xtol=1e-13)
    assert 2.3 < asymptotic_y < 2.4
    assert max(refinement) < 2e-8
    q, delay = 100000, .1
    epsilon = 1/q
    return dict(model='Fixed nearly flat exterior of two opposite point magnetic charges; continuum flux-line NS state and leading MMP energy minimization. Reproduction and extension of a published approximation, not a new global spacetime',
        coordinates='All path coordinates divided by mouth separation d. Prolate s in [-1,1], u=(1-s²+sqrt((1-s²)²+nu²*s²))/nu, z=u*s/2, rho=sqrt((u²-1)*(1-s²))/2. Flux label nu=cos(theta_lower)-cos(theta_upper) is uniform on(0,2) with measure dnu/2',
        optical_length='Leading matched-region circumference L_nu=pi*ell+d*f(nu). Flat-region length equals optical time with c=1. This is not2*atan(4)+d*f/ell: the omitted throat and mouth overlap must also be included. The finite Polyakov cutoff is not a point-charge mouth boundary',
        source_formula='Omega2_bar=integral dnu/2*(2pi/(pi+(d/ell)*f(nu)))²; leading averaged beta/alpha=(Omega2_bar-1)/3. Channel dependence generates angular stress; its average does not solve nonspherical Einstein equations',
        path_bound='f(nu)>=2/nu>=1 by splitting the path at its midplane crossing. nu=0 is a divergent limiting path and is never treated as a finite closed channel',
        mean_path_warning='The continuum integral of f(nu) diverges logarithmically at nu=0. The inverse-period Casimir integrals converge. Replacing the path distribution by its continuum mean delay is not defined',
        field_line_samples=field_lines, max_flux_label_residual=max(flux_errors),
        fixed_scale_cases=fixed,
        matching_formula='ell_v=16*re³/(qG)=4re/(pi*alpha), z=ell/ell_v, delta=d/ell_v. E*ell_v/q=1/(16z²)+1/(24z)-pi/6*integral dnu/2/(pi*z+delta*f). Extremum:z*(Omega2_bar(delta/z)-1)=3',
        matched_cases=matching, asymptotic_throat_over_exterior=asymptotic_y,
        causality='Within this leading stationary approximation beta>0 and f>=1 imply pi*ell>d. Positive energy curvature only tests the chosen ell mode, not dynamical stability or a time offset',
        tail_diagnostic=dict(q=q, separation_over_ell=delay, illustrative_nu_cutoff=epsilon,
            continuum_flux_fraction=epsilon/2,
            omitted_omega_squared_integral_upper_bound=math.pi**2*epsilon**3/(6*delay**2),
            interpretation='Only a bound on the omitted formal Casimir integrand using f>=2/nu. nu*q near1 lacks controlled field-line localization;1/q is a diagnostic cutoff, not a sharp physical threshold. This is not a bound on leakage, changed anomaly, finite-q errors or geometry'),
        maximum_quadrature_change=max(refinement),
        scope='No exact finite-radius mouth matching, nonspherical gravitational response, field leakage, state preparation, formation or time machine. No replacement of the prior specified Polyakov source or corrected interior. Static slow-mouth, weak-gravity flat-region, large-q localization and leading matched-region assumptions require separate control. Positive variational curvature is not full stability',
        sources=['https://arxiv.org/pdf/1807.04726 (Eqs3.13-3.17,5.40,5.44-5.47; AppD-E). Prolate parametrization, numerical diagnostics and bound stated here are our calculations'])


def _magnetic_loop_source(position, throat_length, weak_mass=0., circumference=None):
    """Leading flux-channel tensor per q on a prescribed conformally flat static patch."""
    point = np.asarray(position, dtype=float)
    upper, lower = point-np.array([0., 0., .5]), point+np.array([0., 0., .5])
    rp, rm = np.linalg.norm(upper), np.linalg.norm(lower)
    identity = np.eye(3)
    magnetic = upper/rp**3-lower/rm**3
    derivative = ((identity/rp**3-3*np.outer(upper, upper)/rp**5)
                  -(identity/rm**3-3*np.outer(lower, lower)/rm**5))
    magnitude = np.linalg.norm(magnetic)
    tangent = magnetic/magnitude
    tangent_derivative = (identity-np.outer(tangent, tangent))@derivative/magnitude
    curvature = tangent_derivative@tangent
    potential = weak_mass*(1/rp+1/rm)
    assert 0 <= potential < .25
    grad_potential = -weak_mass*(upper/rp**3+lower/rm**3)
    hess_potential = weak_mass*(3*np.outer(upper, upper)/rp**5-identity/rp**3
                               +3*np.outer(lower, lower)/rm**5-identity/rm**3)
    lapse, spatial = math.sqrt(1-2*potential), math.sqrt(1+2*potential)
    grad_sigma = -grad_potential/(1-2*potential)
    hess_sigma = (-hess_potential/(1-2*potential)
                  -2*np.outer(grad_potential, grad_potential)/(1-2*potential)**2)
    grad_spatial = grad_potential/(1+2*potential)
    parallel = tangent@grad_sigma/spatial
    parallel_derivative = (tangent@hess_sigma@tangent+curvature@grad_sigma
                           -(tangent@grad_sigma)*(tangent@grad_spatial))/spatial**2
    nu = lower[2]/rm-upper[2]/rp
    assert 0 < nu < 2
    if circumference is None:
        circumference = throat_length+_magnetic_return_length(float(nu))
    assert math.isfinite(circumference) and circumference>0
    casimir = -math.pi/(6*lapse**2*circumference**2)
    density2 = casimir+(2*parallel_derivative+parallel**2)/(24*math.pi)
    pressure2 = casimir-parallel**2/(24*math.pi)
    # Physical B/q and n_B/q; all coordinates and lengths are in units of d.
    magnetic_per_q = magnetic/(2*spatial**2)
    channels_per_q = magnitude/(4*math.pi*spatial**2)
    density, pressure = channels_per_q*density2, channels_per_q*pressure2
    project = identity-np.outer(tangent, tangent)
    curvature_physical = (curvature-project@grad_spatial)/spatial
    acceleration_perp = project@grad_sigma/spatial
    force = pressure*curvature_physical+density*acceleration_perp
    current = np.cross(magnetic_per_q, force)/(magnetic_per_q@magnetic_per_q)
    return dict(nu=nu, circumference=circumference, lapse=lapse, spatial=spatial,
        grad_sigma=grad_sigma, grad_spatial=grad_spatial, tangent=tangent,
        coordinate_speed_direction=-tangent, curvature=curvature,
        magnetic=magnetic, magnetic_derivative=derivative, magnitude=magnitude,
        magnetic_per_q=magnetic_per_q, channels_per_q=channels_per_q,
        density=density, pressure=pressure, density2=density2, pressure2=pressure2,
        density2_anomaly=(2*parallel_derivative+parallel**2)/(24*math.pi),
        pressure2_anomaly=-parallel**2/(24*math.pi),
        curvature_physical=curvature_physical, acceleration_perp=acceleration_perp,
        stress_mixed=pressure*np.outer(tangent, tangent), force_physical=force,
        current_physical=current)


def magnetic_loop_source_checks(paths, optical):
    """Off-shell loop action/current checks; no coupled source or mouth solution."""
    throat_length = math.pi*optical['parameters']['ell_over_re']/optical['parameters']['separation_over_re']
    branch = paths['matched_cases'][0]
    assert math.isclose(throat_length, math.pi*branch['ell_over_ell_v']/branch['separation_over_ell_v'], rel_tol=2e-15)
    mass = 1/optical['parameters']['separation_over_re']
    charge = 100000
    source = lambda p, m=0.: _magnetic_loop_source(p, throat_length, m)
    profile, adiabatic = [], []
    for nu in (.5, 1., 1.5):
        samples = []
        for s in np.linspace(-.98, .98, 197):
            rho, z, speed = _magnetic_return_path(nu, float(s))
            row = source([rho, 0., z])
            assert abs(row['nu']-nu) < 1e-13
            b, db = row['magnetic'], row['magnetic_derivative']
            ell_b = math.sqrt(2/(charge*row['magnitude']))
            gradient_measure = ell_b*np.linalg.norm(db.T@b)/row['magnitude']**2
            curvature_measure = ell_b*np.linalg.norm(row['curvature'])
            # Orthonormal angular direction relative to the upper mouth.
            radial = np.array([rho, 0., z-.5]); radial /= np.linalg.norm(radial)
            theta = np.array([radial[2], 0., -radial[0]])
            radial_projection, angular_projection = row['tangent']@radial, row['tangent']@theta
            samples.append(dict(s=float(s), rho_over_d=rho, z_over_d=z,
                density_times_d4_over_q=row['density'], parallel_pressure_times_d4_over_q=row['pressure'],
                radial_pressure_times_d4_over_q=row['pressure']*radial_projection**2,
                polar_pressure_times_d4_over_q=row['pressure']*angular_projection**2,
                mixed_pressure_times_d4_over_q=row['pressure']*radial_projection*angular_projection,
                azimuthal_current_times_d3=float(row['current_physical'][1]),
                force_times_d5_over_q=row['force_physical'].tolist(),
                adiabatic_gradient=gradient_measure, adiabatic_curvature=curvature_measure))
        profile.append(dict(nu=nu, circumference_over_d=throat_length+_magnetic_return_length(nu), samples=samples))
        adiabatic.append(dict(nu=nu, sampled_max_gradient=max(v['adiabatic_gradient'] for v in samples),
            sampled_max_curvature=max(v['adiabatic_curvature'] for v in samples)))
    assert max(v['sampled_max_gradient'] for v in adiabatic) < .019
    assert max(v['sampled_max_curvature'] for v in adiabatic) < .019

    # Independent Cartesian differentiation of T^j_i, including the 4D connection.
    ward = []
    points = [(.5, 0.), (1., .4), (1.5, -.6)]
    for weak_mass in (0., mass, .01):
        for nu, s in points:
            rho, z, _ = _magnetic_return_path(nu, s)
            position = np.array([rho, 0., z]); base = source(position, weak_mass)
            errors = []
            for step in (2e-4, 1e-4):
                divergence = np.zeros(3)
                for axis in range(3):
                    unit = np.eye(3)[axis]*step
                    values = [source(position+j*unit, weak_mass)['stress_mixed'][axis]
                              for j in (-2, -1, 1, 2)]
                    divergence += (values[0]-8*values[1]+8*values[2]-values[3])/(12*step)
                divergence += ((base['grad_sigma']+3*base['grad_spatial'])@base['stress_mixed']
                               -base['pressure']*base['grad_spatial']+base['density']*base['grad_sigma'])
                divergence /= base['spatial']  # Covector components to the orthonormal frame.
                scale = max(np.linalg.norm(base['force_physical']), abs(base['pressure']), abs(base['density']))
                errors.append(float(np.linalg.norm(divergence-base['force_physical'])/scale))
                assert errors[-1] < 2e-8
            assert np.linalg.norm(np.cross(base['current_physical'], base['magnetic_per_q'])
                                  -base['force_physical']) < 1e-17
            ward.append(dict(nu=nu,s=s,weak_mass_over_d=weak_mass,relative_errors=errors))

    # Deform the actual saved return curve at fixed channel flux. The C3 bump
    # and its first derivatives vanish at both segment ends, so endpoint work is zero.
    def shape_check(nu, weak_mass, order):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        endpoint = .6
        nodes, weights = endpoint*nodes, endpoint*weights
        length0 = throat_length+_magnetic_return_length(nu)
        prepared = []
        predicted, predicted_anomaly = 0., 0.
        for s, weight in zip(nodes, weights):
            rho, z, speed = _magnetic_return_path(nu, float(s))
            position = np.array([rho, 0., z]); base = source(position, weak_mass)
            velocity = -speed*base['tangent']
            bump = (1-(s/endpoint)**2)**4
            bump_prime = -8*s/endpoint**2*(1-(s/endpoint)**2)**3
            displacement = np.array([bump, 0., 0.])
            factor = -weight*base['lapse']*base['spatial']**2*speed
            force2 = base['pressure2']*base['curvature_physical']+base['density2']*base['acceleration_perp']
            anomaly2 = (base['pressure2_anomaly']*base['curvature_physical']
                        +base['density2_anomaly']*base['acceleration_perp'])
            predicted += factor*(force2@displacement)
            predicted_anomaly += factor*(anomaly2@displacement)
            sigma_derivative = base['grad_sigma']@velocity
            anomaly_integrand = base['lapse']/base['spatial']*sigma_derivative**2/speed
            prepared.append((position, velocity, bump, bump_prime, weight,
                             base['spatial']/base['lapse']*speed, anomaly_integrand))
        def action_change(amplitude):
            delta_length, delta_anomaly = 0., 0.
            for position, velocity, bump, bump_prime, weight, optical0, anomaly0 in prepared:
                changed = position+np.array([amplitude*bump,0.,0.])
                derivative = velocity+np.array([amplitude*bump_prime,0.,0.])
                speed = np.linalg.norm(derivative)
                # Off-shell geometry alone: do not recompute a different channel label/state.
                upper, lower = changed-np.array([0.,0.,.5]), changed+np.array([0.,0.,.5])
                rp, rm = np.linalg.norm(upper), np.linalg.norm(lower)
                potential = weak_mass*(1/rp+1/rm)
                lapse, spatial = math.sqrt(1-2*potential), math.sqrt(1+2*potential)
                grad_sigma = weak_mass*(upper/rp**3+lower/rm**3)/(1-2*potential)
                delta_length += weight*(spatial/lapse*speed-optical0)
                delta_anomaly += weight*(lapse/spatial*(grad_sigma@derivative)**2/speed-anomaly0)/(24*math.pi)
            return -math.pi*delta_length/(6*length0*(length0+delta_length))+delta_anomaly, delta_anomaly
        derivatives = []
        for step in (2e-4, 1e-4, 5e-5):
            plus, minus = action_change(step), action_change(-step)
            derivatives.append(dict(step=step, full=(plus[0]-minus[0])/(2*step),
                                    anomaly=(plus[1]-minus[1])/(2*step)))
        relative = abs(derivatives[-1]['full']-predicted)/max(abs(predicted),1e-20)
        assert relative < 3e-7
        anomaly_relative = 0. if weak_mass==0 else abs(derivatives[-1]['anomaly']-predicted_anomaly)/abs(predicted_anomaly)
        assert anomaly_relative < 3e-7
        return dict(nu=nu,weak_mass_over_d=weak_mass,quadrature_order=order,
            action_derivative_predicted=predicted, anomaly_derivative_predicted=predicted_anomaly,
            finite_difference_derivatives=derivatives, relative_action_error=relative,
            relative_anomaly_error=anomaly_relative)
    shapes = [shape_check(nu, weak_mass, 64) for nu in (.5,1.,1.5) for weak_mass in (0.,mass,.01)]
    for case in shapes:
        refined = shape_check(case['nu'],case['weak_mass_over_d'],96)
        case['quadrature_prediction_change'] = abs(refined['action_derivative_predicted']-case['action_derivative_predicted'])
        assert case['quadrature_prediction_change'] < 1e-14
    return dict(model='Leading adiabatic flux-channel Dirac source on saved magnetic return curves. Flat exterior profiles; separate local prescribed weak-metric action tests. No solved semiclassical source or joined geometry',
        parameters=dict(q=charge,flux_interval=[.5,1.5],fraction_of_total_flux=.5,
            profile_prolate_interval=[-.98,.98],shape_support_prolate_interval=[-.6,.6],
            throat_optical_length_over_d=throat_length,physical_weak_mass_over_d=mass,
            amplified_off_shell_test_mass_over_d=.01),
        units='Natural hbar=c=1, electric fermion charge1 and flux2*pi*q. Coordinates in d. Report T*d^4/q, force*d^5/q and j*d^3. No SI current/material realization or rescaled human probe inferred',
        channel_action='Gamma_channel/time=pi/(6L)+(1/(24pi))*integral dl*N*(d_l logN)^2. Closed static antiperiodic channel; fixed local/reference terms. L=integral dl/N on a completed channel',
        state_scope='For flat profiles L/d=pi*ell/d+f(nu), the saved leading matched-region estimate. Nonconstant-lapse tests retain this baseline L and vary only a compact segment optical length; the unvaried complement is conditional state data, not a solved full metric. The .01 mass is an amplified off-shell identity test, not a physical correction prediction',
        prescribed_metric='On test patches N=sqrt(1-2U), A=sqrt(1+2U), gamma_ij=A^2 delta_ij, U=mu*(1/r_plus+1/r_minus). Exact as a prescribed test metric only; physical weak-gravity meaning through first order. The flux is the original flat point-monopole field, not a solved weak-metric Maxwell field',
        stress='rho=n_B[-pi/(6N²L²)+(2*d_l a_parallel+a_parallel²)/(24pi)], P=n_B[-pi/(6N²L²)-a_parallel²/(24pi)], a_parallel=d_l logN; T=rho*u*u+P*b*b',
        action_variation='deltaGamma/time=-integral dl*N*(p2*kappa+rho2*a_perp).eta per channel. At fixed flux a transverse displacement corresponds to deltaA_phi=-eta.grad A_phi. Both numerical action variation and independent Cartesian tensor differentiation are checked',
        current='f=P*kappa+rho*a_perp; gauge/shape variation predicts j_perp=B cross f/B². j is azimuthal for the poloidal axisymmetric source and has zero divergence. It is generally nonzero although the original flat point field is curl-free off the mouths; Maxwell backreaction is therefore required',
        sign_guard='With deltaGamma=+integral j^mu deltaA_mu and Maxwell action -F²/(4g²): nabla_mu F^(mu nu)=-g²j^nu, consistent with divT=F.j',
        flux_boundary='Band measure is q*dnu/2, not a normalized all-q replacement. No physical source taper is imposed at nu=.5 or1.5. Sampled local estimates do not certify all exterior points, mouth transfer, discarded channels or a global Hadamard state',
        endpoint_terms='General segment boundary work is [N*p2*b.eta+(a_parallel/(12pi))*delta_eta N]_ends. The radial deformation bump(1-(s/.6)^2)^4 and its first derivatives vanish at both support ends; these terms vanish in the tests. No physical attachment calculation is performed',
        profile=profile, adiabatic_samples=adiabatic,ward_checks=ward,shape_variation_checks=shapes,
        limits='Angular pressure from bent poloidal channels is negative semidefinite; azimuthal pressure is zero at this order. This cannot supply the isotropic positive taper pressure. No quantum correction to Maxwell/Einstein solved, no stable throat, creation mechanism, source preparation or13-to12Sep2027 route. Higher modes, spin transport, renormalization, flux boundaries and leakage remain',
        sources=['https://arxiv.org/pdf/1807.04726 (MMP leading flux channels and magnetic return paths)',
                 'https://arxiv.org/pdf/2606.12528 (fermion reduction on a given background, not our source closure)',
                 'QUANTUM_COMPLETION.md (our reviewed variational extension; current numerical tests are ours)'])


def magnetic_patch_response_checks(loop_source, optical):
    """Linear flat Maxwell response with diagnostic, explicitly artificial boundaries."""
    from scipy.sparse import diags, eye, kron
    from scipy.sparse.linalg import spsolve

    throat = loop_source['parameters']['throat_optical_length_over_d']
    charge = loop_source['parameters']['q']
    alpha = optical['parameters']['alpha']
    coupling = 8*math.pi**2*alpha  # 2*g^2/q; delta_u=coupling*w, u=2*psi/q.

    def solve_patch(a, b, extent, count):
        radial, axial = np.linspace(a, b, count), np.linspace(-extent, extent, count)
        hr, hz = radial[1]-radial[0], axial[1]-axial[0]
        ri, zi = radial[1:-1], axial[1:-1]
        nr, nz = len(ri), len(zi)
        face = 1/((radial[:-1]+radial[1:])/2)/hr**2
        minus, plus = face[:-1], face[1:]
        radial_operator = diags([minus[1:], -minus-plus, plus[:-1]], [-1, 0, 1])
        axial_operator = diags([np.ones(nz-1), -2*np.ones(nz), np.ones(nz-1)], [-1, 0, 1])/hz**2
        operator = (kron(radial_operator, eye(nz))
                    +kron(diags(1/ri), axial_operator)).tocsr()
        assert (operator-operator.T).nnz == 0
        current, magnetic, labels, adiabatic = [], [], [], []
        for rho in radial:
            for z in axial:
                source = _magnetic_loop_source([rho, 0., z], throat)
                independently_reduced = -np.cross(source['tangent'], source['curvature'])/(12*source['circumference']**2)
                assert np.linalg.norm(source['current_physical']-independently_reduced) < 2e-18
                current.append(float(source['current_physical'][1]))
                magnetic.append(source['magnetic'][[0, 2]]/2)
                labels.append(source['nu'])
                ell_b = math.sqrt(2/(charge*source['magnitude']))
                adiabatic.append(ell_b*max(np.linalg.norm(source['curvature']),
                    np.linalg.norm(source['magnetic_derivative'].T@source['magnetic'])/source['magnitude']**2))
        current = np.asarray(current).reshape(count, count)
        magnetic = np.asarray(magnetic).reshape(count, count, 2)
        assert min(labels) > .5 and max(labels) < 1.5
        rhs = -current[1:-1, 1:-1].ravel()  # div(grad(w)/rho)=-j_phi.
        inside = spsolve(operator, rhs)
        residual = float(np.max(np.abs(operator@inside-rhs))/np.max(np.abs(rhs)))
        assert residual < 2e-10
        w = np.zeros((count, count)); w[1:-1, 1:-1] = inside.reshape(nr, nz)
        dw_dr, dw_dz = np.gradient(w, hr, hz, edge_order=2)
        correction = coupling*np.stack([-dw_dz/radial[:, None], dw_dr/radial[:, None]], axis=-1)/2
        relative_field = np.linalg.norm(correction, axis=-1)/np.linalg.norm(magnetic, axis=-1)

        # An independent polynomial manufactured solution tests sign, boundaries and order.
        rr, zz = np.meshgrid(ri, zi, indexing='ij')
        polynomial = rr**2*(rr**2-a*a)*(rr**2-b*b)
        exact = polynomial*(extent**2-zz**2)
        forcing = ((24*rr**4-8*(a*a+b*b)*rr**2)*(extent**2-zz**2)-2*polynomial)/rr
        manufactured = spsolve(operator, forcing.ravel()).reshape(nr, nz)
        manufactured_error = float(np.max(np.abs(manufactured-exact))/np.max(np.abs(exact)))

        # Continuous bound: |Bunit|>=rmax^-3, ||DBunit||<=4*rmin^-3,
        # L>=throat+1. The z-only barrier also bounds the discrete solution.
        rmin, rmax = math.hypot(a, .5-extent), math.hypot(b, .5+extent)
        current_bound = 4*(rmax/rmin)**3/(12*(throat+1)**2)
        w_bound = b*current_bound*extent**2/2
        assert np.max(np.abs(w)) <= w_bound
        centre = count//2
        record = dict(radial_interval=[a,b], axial_interval=[-extent,extent], grid_count=count,
            spacing=[float(hr),float(hz)], flux_label_sample_range=[min(labels),max(labels)],
            sampled_max_adiabatic=max(adiabatic), max_current_times_d3=float(np.max(np.abs(current))),
            current_absolute_analytic_bound=current_bound, w_absolute_analytic_bound=w_bound,
            delta_u_absolute_analytic_bound=coupling*w_bound,
            max_abs_w=float(np.max(np.abs(w))), centre_w=float(w[centre,centre]),
            centre_delta_b_over_q=correction[centre,centre].tolist(),
            sampled_max_relative_field_interior=float(np.max(relative_field[2:-2,2:-2])),
            relative_linear_residual=residual, manufactured_relative_error=manufactured_error)
        return record, w, correction, current, magnetic

    base = [solve_patch(.7,1.1,.2,n) for n in (33,65,129)]
    enlarged = [solve_patch(.6,1.2,.3,n) for n in (97,193)]
    differences = [float(np.max(np.abs(base[i+1][1][::2,::2]-base[i][1]))) for i in (0,1)]
    assert 3.5 < differences[0]/differences[1] < 4.5
    assert 3.5 < base[0][0]['manufactured_relative_error']/base[1][0]['manufactured_relative_error'] < 4.5
    assert 3.5 < base[1][0]['manufactured_relative_error']/base[2][0]['manufactured_relative_error'] < 4.5
    # Identical spacing and overlapping nodes isolate the chosen boundary change.
    overlaps = [enlarged[0][1][16:81,16:81], enlarged[1][1][32:161,32:161]]
    boundary_changes = [float(np.max(np.abs(overlaps[i]-base[i+1][1]))) for i in (0,1)]
    boundary_refinement = abs(boundary_changes[1]-boundary_changes[0])
    assert boundary_refinement < .01*boundary_changes[-1]
    finest, w, correction, current, magnetic = base[-1]
    # ponytail: a finite patch diagnoses Maxwell response; global mouth/state matching is a separate unsolved problem.
    epsilon = finest['max_abs_w']/1.1**2
    centre_baseline = np.linalg.norm(correction[64,64])
    harmonic_field = coupling*epsilon  # w_h=epsilon*rho^2 => delta_Bz/q=coupling*epsilon.
    return dict(model='First quantum-current Maxwell response on a flat exterior rectangle with zero Dirichlet correction. Conditional finite-patch diagnostic, not a global solution',
        parameters=dict(q=charge,alpha=alpha,g_squared=4*math.pi**2*charge*alpha,
            two_g_squared_over_q=coupling,throat_optical_length_over_d=throat),
        units='Natural hbar=c=1 and unit electric fermion charge; rho,z,L in d. u=2*psi/q, w=delta_u/(2g²/q), jhat=j_phi*d³, reported delta_B*d²/q',
        equation='Delta_star w=-rho*jhat_phi; Delta_star=partial_rho²-rho^-1 partial_rho+partial_z². Symmetric discretization of div(grad w/rho)=-jhat_phi with reciprocal-radius face coefficients',
        source='Evaluate the existing flat complete-loop source at every grid point; j=-b cross kappa/(12L²), L/d=throat+f(nu). Current is frozen on the unperturbed geometry/state for this first response only',
        boundary='delta_u=0 on each artificial rectangle. These are not physical conducting walls, mouth data or flux-band edges. No current is cut off at the rectangle or at nu=.5,1.5; source outside the rectangle is unspecified by this solve',
        continuum_bound='rmin=hypot(a,.5-Z),rmax=hypot(b,.5+Z); |jhat|<=4(rmax/rmin)^3/[12(throat+1)^2]. Maximum principle with barrier b*j_bound*(Z²-z²)/2 gives |w|<=b*j_bound*Z²/2 for zero Dirichlet data; this bounds flux correction, not field derivatives',
        grid_results=[row[0] for row in base], enlarged_grid_results=[row[0] for row in enlarged],
        verification=dict(grid_restriction_max_changes=differences,refinement_ratio=differences[0]/differences[1],
            enlarged_boundary_max_changes_on_base=boundary_changes,boundary_change_refinement=boundary_refinement),
        homogeneous_boundary_example=dict(epsilon=epsilon,solution='w_h=epsilon*rho² solves Delta_star w_h=0 exactly; add its boundary values to the original solve',
            max_boundary_delta_u=coupling*epsilon*1.1**2,uniform_delta_bz_over_q=harmonic_field,
            field_change_over_centre_source_response=harmonic_field/centre_baseline,
            meaning='Different unspecified boundary data give a different field with the same interior current. This is boundary dependence, not an instability test'),
        plot=dict(radial=np.linspace(.7,1.1,129)[::4].tolist(),axial=np.linspace(-.2,.2,129)[::4].tolist(),
            w=w[::4,::4].tolist(),current_times_d3=current[::4,::4].tolist(),
            delta_b_over_q=correction[::4,::4].tolist(),background_b_over_q=magnetic[::4,::4].tolist()),
        scope='Only the flat quantum contribution, not the previously separate weak-metric forcing. No updated loop lengths/current susceptibility, nonlocal Newton solve, Einstein response, omitted-mode control, physical boundary matching, dynamical stability, material preparation or chronology change. Small local correction does not imply a completed medium or13-to12Sep2027 travel',
        sources=['https://arxiv.org/pdf/1807.04726 (leading magnetic flux channels and return geometry)',
                 'https://arxiv.org/pdf/gr-qc/0302106 (axisymmetric Maxwell/Grad-Shafranov framework; no plasma closure imported)',
                 'QUANTUM_COMPLETION.md section7 (our reviewed operator/current extension); discretization, bound and results here are our calculations'])


def sourced_rn_transfer_checks(loop_source, optical, exterior):
    """Casimir dipole forcing and affine RN transfer, with unknown physical end data."""
    import mpmath as mp

    t, r, theta, phi = sp.symbols('t r theta phi', real=True)
    mass, charge, amplitude = sp.symbols('M Q A', positive=True)
    lapse2 = 1-2*mass/r+charge**2/r**2
    lam = 6*mass/r-4*charge**2/r**2
    ht, hr, aa = [sp.Function(name)(r) for name in ('Ht','Hr','a')]
    metric = (-lapse2,1/lapse2,r*r,r*r*sp.sin(theta)**2)
    perturbation = (metric[0]*ht*sp.cos(theta),metric[1]*hr*sp.cos(theta),sp.Integer(0),sp.Integer(0))
    f0, f1 = sp.zeros(4), sp.zeros(4)
    f0[2,3] = 2*charge*sp.sin(theta); f0[3,2] = -f0[2,3]
    potential = (0,0,0,aa*sp.sin(theta)**2)
    coordinates = (t,r,theta,phi)
    for i in range(4):
        for j in range(4):
            f1[i,j] = sp.diff(potential[j],coordinates[i])-sp.diff(potential[i],coordinates[j])
    clean = lambda value: sp.factor(sp.trigsimp(sp.expand_trig(value)))
    reversed_einstein, maxwell, _ = _linear_diagonal_einstein_maxwell(
        coordinates,metric,perturbation,f0,f1,r*r*sp.sin(theta),clean)
    trace = clean(sum(reversed_einstein[i,i]/metric[i] for i in range(4)))
    einstein = {key:clean(value-(metric[key[0]]*trace/2 if key[0]==key[1] else 0))
                for key,value in reversed_einstein.items()}

    # The source is the radial channel's l=1 state term, not a manufactured master forcing.
    rho = amplitude/(r*r*lapse2)
    assert sp.factor(sp.diff(rho,r)+sp.diff(lapse2,r)*rho/lapse2+2*rho/r) == 0
    aval, bval, hval = sp.symbols('a_value a_prime Hr_value', real=True)
    ht_value = (4*charge*aval/r**2+4*charge*lapse2*bval/r
                -r*sp.diff(lapse2,r)*hval-2*r*r*rho)/lam
    hr_prime = r*rho/lapse2-(2*r*r-charge**2)*hval/(lapse2*r**3)+2*charge*aval/(lapse2*r**3)
    b_prime = ((2*aval+charge*(ht_value+hval))/r**2-sp.diff(lapse2,r)*bval)/lapse2
    derivative = lambda value: sp.factor(sp.diff(value,r)+sp.diff(value,aval)*bval
        +sp.diff(value,bval)*b_prime+sp.diff(value,hval)*hr_prime)
    ht_prime = derivative(ht_value)
    substitutions = {sp.diff(ht,r,2):derivative(ht_prime),sp.diff(ht,r):ht_prime,
        sp.diff(hr,r,2):derivative(hr_prime),sp.diff(hr,r):hr_prime,
        sp.diff(aa,r,2):b_prime,sp.diff(aa,r):bval,ht:ht_value,hr:hval,aa:aval}
    source = {(0,0):lapse2*rho*sp.cos(theta),(1,1):rho/lapse2*sp.cos(theta)}
    field_residuals = [clean((value-source.get(key,0)).subs(substitutions)) for key,value in einstein.items()]
    field_residuals += [clean(value.subs(substitutions)) for value in maxwell]
    assert all(value == 0 for value in field_residuals)

    potential_u = sp.sympify(exterior['general_polar_potential_U'],locals={'r':r,'M':mass,'Q':charge})
    master = 2*aval-4*charge*lapse2*hval/lam
    forcing = sp.factor(derivative(lapse2*derivative(master))-potential_u*master)
    expected_forcing = 8*amplitude*charge**3/(3*mass*r-2*charge**2)**2
    assert sp.factor(forcing-expected_forcing) == 0
    gauge = sp.symbols('V',real=True)
    assert sp.factor(master.subs({aval:2*charge*gauge,hval:lam*gauge/lapse2})) == 0

    throat = loop_source['parameters']['throat_optical_length_over_d']
    def moment_integrand(nu):
        return 1.5*(1-nu)/(throat+_magnetic_return_length(float(nu)))**2
    moment = quad(moment_integrand,.5,1.5,epsabs=1e-17,epsrel=2e-12)[0]
    split = sum(quad(moment_integrand,a,b,epsabs=1e-17,epsrel=2e-12)[0] for a,b in ((.5,1.),(1.,1.5)))
    nodes, weights = np.polynomial.legendre.leggauss(64)
    gauss = sum(w*moment_integrand(1+.5*x)*.5 for x,w in zip(nodes,weights))
    assert moment < 0 and max(abs(moment-split),abs(moment-gauss)) < 1e-16
    alpha, separation = [optical['parameters'][key] for key in ('alpha','separation_over_re')]
    amplitude_factor = -4*math.pi**2*alpha/(3*separation**2)
    retained_amplitude = amplitude_factor*moment

    # Reuse the verified physical homogeneous basis; do not evolve the exponential pure gauge inward.
    growing, decaying = [sp.sympify(exterior['exact_static_basis'][key],locals={'r':r}) for key in ('growing','decaying')]
    fext, uext = lapse2.subs({mass:1,charge:1}),potential_u.subs({mass:1,charge:1})
    basis = sp.Matrix([[growing,decaying],[fext*sp.diff(growing,r),fext*sp.diff(decaying,r)]])
    assert sp.factor(basis.det()+3) == 0
    phi_mp = sp.lambdify(r,basis,'mpmath')
    forcing_mp = sp.lambdify(r,expected_forcing.subs({mass:1,charge:1,amplitude:1}),'mpmath')
    particular = -4*r/(3*(3*r-2))
    assert sp.factor(sp.diff(fext*sp.diff(particular,r),r)-uext*particular-8/(3*r-2)**2) == 0
    particular_mp = sp.lambdify(r,sp.Matrix([particular,fext*sp.diff(particular,r)]),'mpmath')
    ff, uu = [sp.lambdify(r,value,'numpy') for value in (fext,uext)]
    rows = []
    with mp.workdps(60):
        def affine(lower,upper):
            lo, hi = mp.mpf(str(lower)),mp.mpf(str(upper))
            philo, phihi = phi_mp(lo),phi_mp(hi)
            transfer = philo*phihi**-1
            integrals = mp.matrix([mp.quad(lambda s:phi_mp(s)[0,1]*forcing_mp(s)/3,[hi,lo]),
                                   mp.quad(lambda s:-phi_mp(s)[0,0]*forcing_mp(s)/3,[hi,lo])])
            closed = particular_mp(lo)-transfer*particular_mp(hi)
            assert max(abs(v) for v in closed-philo*integrals) < mp.mpf('1e-44')
            return transfer,closed
        tf = lambda mat:[[float(mat[i,j]) for j in range(mat.cols)] for i in range(mat.rows)]
        for inner in (1.1,1.01):
            direct, response = affine(inner,4.)
            first, first_response = affine(inner,2.)
            second, second_response = affine(2.,4.)
            composition_error = max(abs(value) for value in direct-first*second)
            source_error = max(abs(value) for value in response-first_response-first*second_response)
            assert composition_error < mp.mpf('1e-44') and source_error < mp.mpf('1e-44')
            assert abs(mp.det(direct)-1) < mp.mpf('1e-40')
            assert response[0] > 0 and response[1] < 0
            mass_deficit = -mp.mpf(str(optical['parameters']['geometric_mass_minus_re_over_re']))
            # This monotone bound controls F only, not all approximation errors.
            relative_f_change = 2*mass_deficit*mp.mpf(str(inner))/(mp.mpf(str(inner))-1)**2
            rows.append(dict(inner_radius_over_re=inner,outer_radius_over_re=4.,transfer=tf(direct),
                source_response_per_A=[float(v) for v in response],
                retained_band_additive_contribution=[retained_amplitude*float(v) for v in response],
                transfer_composition_absolute_error=float(composition_error),
                source_composition_absolute_error=float(source_error),
                relative_extremal_F_approximation_bound=float(relative_f_change)))

    def ode(radius,flat):
        values = flat.reshape(2,3)
        result = np.array([[0.,1/ff(radius)],[uu(radius),0.]])@values
        result[1,2] += 8/(3*radius-2)**2
        return result.ravel()
    initial = np.array([[1.,0.,0.],[0.,1.,0.]]).ravel()
    solutions = [solve_ivp(ode,(4.,1.01),initial,method='DOP853',dense_output=True,
        rtol=2e-12/level,atol=2e-13/level,max_step=.02/level) for level in (1,2)]
    assert all(solution.success for solution in solutions)
    for row in rows:
        expected = np.column_stack([np.asarray(row['transfer']),row['source_response_per_A']])
        errors = [float(np.max(np.abs(solution.sol(row['inner_radius_over_re']).reshape(2,3)-expected)/np.maximum(1.,np.abs(expected))))
                  for solution in solutions]
        assert max(errors) < 2e-9
        row['independent_ODE_scaled_errors'] = errors
    radii = 1+np.geomspace(.01,3.,161)
    profile = solutions[-1].sol(radii).reshape(2,3,-1)[:,2,:]
    return dict(model='Leading radial-channel Casimir dipole drives the coupled static Einstein-Maxwell RN master. Affine transfer retains unknown end data and omitted-channel moment; no completed mouth matching',
        units='Natural hbar=c=1. Radius and RN time in r_e; geometric Q and M in r_e. Fcal=sqrt(16piG)*F_original/g_EM, background Fcal_theta_phi=2Q*sin(theta). S=8piG*T; computed dimensionless source is r_e²*S',
        source_convention='nabla_mu Fcal^(mu nu)=-Jcal^nu, Jcal=sqrt(16piG)*g_EM*j_original; div S=Fcal.Jcal/2. S_tt=F*rho*cos(theta), S_rr=pr*cos(theta)/F, S_rtheta=-tau*sin(theta), S_AB=r²*pt*Omega_AB*cos(theta)',
        ward='pr_prime+F_prime*(rho+pr)/(2F)+2*(pr-pt)/r-2*tau/r²=0; F*tau_prime+(F_prime+2F/r)*tau+pt=-Q*Jcal_phi_contravariant',
        source=dict(rho='rho=pr=A/(r²F), pt=tau=Jcal_phi=0; trace0; exact Ward identities at the stated leading radial background',
            origin='Radial channel state term rho_original=pr_original=-q/[24*r_physical²*F*L_RN(nu)²]. Near the positive mouth nu=1-cos(theta) at leading radial order. The radial geometric anomaly is angle independent here and has no l1 projection',
            retained_flux_interval=[.5,1.5],retained_fraction_of_flux=.5,retained_dipole_moment=moment,
            split_quadrature_difference=abs(moment-split),gauss_quadrature_difference=abs(moment-gauss),
            amplitude_per_total_moment_at_common_clock=amplitude_factor,retained_amplitude_at_common_clock=retained_amplitude,
            total_amplitude='A=-(4*pi²*alpha/[3*(d/re)²])*C_clock^-2*(I1_retained+I1_unknown). I1=(3/2)*integral dnu*(1-nu)/(L_global(nu)/d)²; t_RN=C_clock*t_global',
            missing='I1_unknown is not set to zero. The numerical retained contribution is additive known-band information, not the total source or a physical source taper. Full optical lengths, holonomies, higher modes, local terms and common clock completion remain conditional'),
        metric_gauge='delta_gtt=gtt*Ht*cos(theta), delta_grr=grr*Hr*cos(theta), other metric entries0; delta_Aphi=a*sin(theta)². Time-reflection-even static magnetic sector only',
        master=dict(definition='H=2a-4QF*Hr/Lambda, Lambda=6M/r-4Q²/r²; P_H=F*H_prime',
            equation='(F*H_prime)_prime-U*H=8*A*Q³/(3*M*r-2*Q²)²',potential_U=str(potential_u),
            forcing_per_A=str(sp.factor(forcing/amplitude)),
            gauge='Residual V_prime=-V/(rF), delta_a=2QV, delta_Hr=Lambda*V/F, delta_Ht=-rF_prime*V/F. H is invariant; raw inward gauge evolution can grow exponentially near extremality and is not used for the transfer'),
        reconstruction='Let b=a_prime. Hr_prime=r*rho/F-(2r²-Q²)Hr/(F*r³)+2Qa/(F*r³); Ht=[4Qa/r²+4QF*b/r-rF_prime*Hr-2r²*pr]/Lambda; b_prime=[(2a+Q(Ht+Hr))/r²-F_prime*b]/F. These sourced relations reconstruct all 4D equations; one coordinate-gauge integration constant is distinct from physical end data',
        transfer_equation='Y(inner)=T(inner,outer)*Y(outer)+A*B(inner,outer), Y=(H,FH_prime). The source column B is the integral of the existing homogeneous fundamental matrix against the derived forcing; outer physical amplitudes remain open',
        exact_extremal_particular='For M=Q=1, H_p/A=-4r/[3(3r-2)], (F H_p_prime)/A=8(r-1)²/[3r²(3r-2)²]. B(inner,outer)=Y_p(inner)-T(inner,outer)Y_p(outer), independently checked against the source integral. This algebraic choice of particular solution imposes no physical horizon/end condition',
        transfer_cases=rows,
        approximation=dict(numerical_background='M=Q=1 on finite r>=1.01 only, never a black-hole horizon boundary condition',
            extremal_F_bound='For the saved M=1-deltaM branch, abs(deltaF/F_ext)<=2*deltaM*r_in/(r_in-1)² for r>=r_in>1. This does not justify extending extremality into the throat or bound source/state errors',
            companion_point_field_fraction_bound=(4/(separation-4))**2,
            companion_scope='Leading point-field companion/self ratio bounded by r_out²/(d-r_out)²; radial RN source is a leading approximation, not the previously computed bent-loop current substituted at a mouth'),
        verification=dict(direct_four_dimensional_Einstein_and_Maxwell_residuals=[0]*len(field_residuals),
            radial_source_Ward_identity=True,source_master_elimination_identity=True,residual_coordinate_gauge_identity=True,
            exact_extremal_particular_identity=True,
            reused_extremal_weighted_Wronskian=-3,independent_ODE_evaluations=[int(s.nfev) for s in solutions]),
        plot=dict(radius_over_re=radii.tolist(),source_response_H_per_A=profile[0].tolist(),
            source_response_FHprime_per_A=profile[1].tolist(),retained_band_H_contribution=(retained_amplitude*profile[0]).tolist()),
        scope='A sourced dipole transfer and reconstruction, not fixed physical mouth amplitudes, a self-consistent full quantum state, exact finite matching, spherical-pressure repair, nonlinear/dynamical stability, matter preparation or chronology change. Unknown outer/inner matching and omitted channels can change the physical answer. No13-to12Sep2027 route',
        sources=['https://arxiv.org/pdf/1807.04726 (leading radial flux-channel/Casimir framework)',
                 'https://arxiv.org/pdf/1212.2557 (AppendixB electric-master definitions; magnetic duality mapping derived here)',
                 'https://arxiv.org/pdf/gr-qc/9507012 (Eq87 homogeneous polar potential)',
                 'Our derivation: the sourced4D equations, radial-Casimir forcing, retained-band projection and affine transfer; these forcing terms are not claimed as quoted published results'])


def companion_incident_checks(optical, exterior, sourced):
    """Classical incident normalization and a source-aware finite-cut relation."""
    import mpmath as mp

    r = sp.symbols('r', positive=True)
    cg, aa, ap = sp.symbols('C_g a a_prime', real=True)
    f, lam = (r-1)**2/r**2, 6/r-4/r**2
    g, d = [sp.sympify(exterior['exact_static_basis'][key], locals={'r':r})
            for key in ('growing', 'decaying')]
    hr = lam*(2*aa-cg*g)/(4*f)
    hr_prime = -(2*r*r-1)*hr/(f*r**3)+2*aa/(f*r**3)
    a_prime = sp.factor(sp.solve(sp.diff(hr,r)+sp.diff(hr,aa)*ap-hr_prime,ap)[0])
    assert sp.factor(a_prime-(3*cg*r/2-aa/(r*f))) == 0
    ht = sp.factor((4*aa/r**2+4*f*a_prime/r-r*sp.diff(f,r)*hr)/lam)
    a_asymptotic = cg*(r*r-r-1)/2
    residual = sp.factor(sp.diff(a_asymptotic,r)-a_prime.subs(aa,a_asymptotic))
    assert sp.limit(r*residual,r,sp.oo) == 0
    limits = [sp.limit(value,r,sp.oo) for value in
        (2*a_asymptotic/r**2, (1-1/r)*sp.diff(a_asymptotic,r)/r,
         hr.subs(aa,a_asymptotic), ht.subs(aa,a_asymptotic)/r)]
    assert limits == [cg,cg,cg/2,cg]
    gauge = sp.exp(1/(r-1))/(r-1)
    assert sp.factor(sp.diff(gauge,r)+gauge/(r*f)) == 0
    assert sp.limit(gauge/r**2,r,sp.oo) == 0

    # Lagrange identity propagates an inhomogeneous boundary condition without
    # choosing the unresolved decaying response or the exterior quantum data.
    impedance = sp.factor(-f*sp.diff(d,r)/d)
    particular = -4*r/(3*(3*r-2))
    source_wronskian = sp.factor(d*f*sp.diff(particular,r)-f*sp.diff(d,r)*particular)
    assert sp.factor(sp.diff(source_wronskian,r)-8*d/(3*r-2)**2) == 0
    assert sp.factor(f*sp.diff(g,r)+impedance*g-3/d) == 0
    phi = sp.Matrix([[g,d],[f*sp.diff(g,r),f*sp.diff(d,r)]])
    transfer = phi*phi.subs(r,4).inv()
    row = sp.Matrix([[impedance,1]])
    assert all(sp.factor(v)==0 for v in row*transfer-d.subs(r,4)/d*row.subs(r,4))
    yp = sp.Matrix([particular,f*sp.diff(particular,r)])
    column = yp-transfer*yp.subs(r,4)
    source_coefficient = sp.factor((source_wronskian-source_wronskian.subs(r,4))/d)
    assert sp.factor((row*column)[0]-source_coefficient) == 0

    separation = optical['parameters']['separation_over_re']
    ell = optical['parameters']['ell_over_re']
    incident = -2/separation**2
    retained = sourced['source']['retained_amplitude_at_common_clock']
    cuts = []
    for value in ('1.01','1.1','4'):
        cut = sp.Rational(value)
        coeff = float(source_coefficient.subs(r,cut))
        cuts.append(dict(radius_over_re=float(cut),
            decaying_impedance=float(impedance.subs(r,cut)),
            leading_classical_incident_forcing=incident*float((3/d).subs(r,cut)),
            exterior_quantum_residual_factor=float((d.subs(r,4)/d).subs(r,cut)),
            source_forcing_per_A=coeff,retained_band_source_forcing=retained*coeff,
            growing_classical_H_component=incident*float(g.subs(r,cut))))

    # Derive the leading radiation coefficient rather than silently importing
    # the differing prefactor printed in MMP Eq6.53.
    time, mass, distance, omega = sp.symbols('t m d omega', positive=True)
    position = sp.Matrix([distance*sp.cos(omega*time)/2,
                          distance*sp.sin(omega*time)/2,0])
    quadrupole = 2*mass*(position*position.T-sp.eye(3)*distance**2/12)
    third = quadrupole.diff(time,3)
    contracted = sp.trigsimp(sum(value**2 for value in third))
    assert sp.simplify(contracted-8*mass**2*distance**4*omega**6) == 0
    dd, correction = sp.symbols('D a', positive=True)
    primitive = dd**3/3-correction*dd**2/2+correction**2*dd-correction**3*sp.log(dd+correction)
    assert sp.factor(sp.diff(primitive,dd)-dd**3/(dd+correction)) == 0

    orbit_rows = []
    with mp.workdps(70):
        D, length = mp.mpf(str(separation)),mp.mpf(str(ell))
        duration = mp.pi*length
        omega0 = 2/D**mp.mpf('1.5')
        shrink_scale = 32*duration/(3*D**2)
        for label, ac in (('leading_EM_only',mp.mpf(0)),
                          ('EM_plus_standard_leading_GW',mp.mpf(48)/5),
                          ('MMP_printed_GW_prefactor_comparison',mp.mpf(4)/5)):
            # Deficit form avoids cancellation of two ~D^3 primitives.
            def elapsed(deficit):
                final = D-deficit
                difference = deficit*((D*D+D*final+final*final)/3
                    -ac*(D+final)/2+ac*ac)-ac**3*mp.log1p(deficit/(final+ac))
                return 3*difference/32
            deficit = mp.findroot(lambda delta:elapsed(delta)-duration,
                                  (shrink_scale,shrink_scale*(1+mp.mpf('1e-6'))))
            integral = 3*mp.quad(lambda z:(D-z)**3/(D-z+ac),[0,deficit])/32
            assert abs(integral/duration-1) < mp.mpf('1e-55')
            eps = shrink_scale/D
            def rhs(_,state):
                ratio = 1-float(eps)*state[0]
                return [(1+float(ac/D)/ratio)/ratio**2]
            solve = solve_ivp(rhs,(0.,1.),[0.],method='DOP853',rtol=3e-13,atol=3e-14,max_step=.05)
            assert solve.success
            numerical_deficit = float(shrink_scale)*solve.y[0,-1]
            error = abs(numerical_deficit/float(deficit)-1)
            assert error < 3e-12
            angle = 3*mp.quad(lambda z:(D-z)**mp.mpf('1.5')/(D-z+ac),[0,deficit])/16
            assert omega0*duration <= angle <= 2*duration/(D-deficit)**mp.mpf('1.5')
            orbit_rows.append(dict(model=label,relative_GW_correction_coefficient=float(ac),
                separation_deficit_over_re=float(deficit),fractional_separation_decrease=float(deficit/D),
                orbital_phase_radians=float(angle),orbital_phase_degrees=float(180*angle/mp.pi),
                independent_deficit_ODE_relative_error=error))
        orbital = dict(separation_over_re=separation,ell_over_re=ell,
            interval_over_re=float(duration),omega_times_re=float(omega0),
            omega_times_ell=float(omega0*length),each_speed_over_c=float(1/mp.sqrt(D)),
            period_over_re=float(2*mp.pi/omega0),
            leading_EM_power_times_G=float(32/(3*D**4)),
            standard_leading_GW_power_times_G=float(512/(5*D**5)),
            MMP_printed_GW_power_times_G=float(128/(15*D**5)),
            cases=orbit_rows)
    return dict(
        model='Leading classical companion incident data in an accelerating weak overlap, propagated through the finite static RN source model; orbital drift diagnostic, not a joined rotating wormhole',
        units='Natural c=hbar=1, radius/time/geometric mass in r_e; M=Q=1, Fcal=2*F_standard. D=d/re. Physical powers restore G*P/c^5, frequency omega*re/c and time t*c/re',
        incident=dict(upper_mouth_charge=1,companion_charge=-1,
            geometric_uniform_B=incident,leading_master_growing_amplitude=incident,
            reconstructed_equation='For H=Cg*g, a_prime+a/(rF)=3*Cg*r/2; a=Cg*(r²-r-1)/2+O(log(r)/r). Gauge a_h is proportional to exp[1/(r-1)]/(r-1)',
            physical_field='delta Bcal_hat_r=2a*cos(theta)/r², delta Bcal_hat_theta=-sqrt(F)*a_prime*sin(theta)/r; their weak-buffer limits give a uniform field Cg along +z',
            metric='Hr->Cg/2 and Ht/r->Cg; delta log N~Cg*z/2. Upper-mouth companion gravity contributes+1/D² to the lapse gradient, frame acceleration-2/D², sum-1/D²=Cg/2',
            scope='A leading classical matching coefficient in 1<<r<<D, not H(r=4)/16 identified with a measured field. Quantum exterior response can add its own incident term; no complete outer state is fixed'),
        finite_cut_relation=dict(equation='P(r)+K_d(r)*H(r)=3*Cg/d_mode(r)+[d_mode(4)/d_mode(r)]*Xi_q(4)+A*J(r)',
            definitions='P=F*H_prime; K_d=-F*d_mode_prime/d_mode; Xi_q(4)=P(4)+K_d(4)*H(4)-3*Cg/d_mode(4). Xi_q includes uncomputed quantum exterior data and corrections to the leading classical matching',
            source_coefficient=str(source_coefficient),source_Wronskian=str(source_wronskian),
            source_identity='W_p=d_mode*F*h_p_prime-F*d_mode_prime*h_p; W_p_prime=d_mode*S_per_A; J=[W_p(r)-W_p(4)]/d_mode(r)',
            remaining='The decaying homogeneous response, Xi_q(4), complete A, common-clock and inner/throat matching remain open. No missing datum is assigned zero',cuts=cuts),
        orbit=orbital,
        leading_common_clock=dict(
            relation='t_RN rate relative to global time on the comoving reference worldline: C=1-U_comp-v²/2=1-3/(2D)+O(D^-2), U_comp=1/D and v²=1/D',
            rate_shift_first_order=-3/(2*separation),
            source_amplitude_factor_first_order=1+3/separation,
            retained_band_amplitude_with_clock_correction_only=retained*(1+3/separation),
            relative_mouth_rate_difference_in_equal_symmetric_model=0.,
            transformation='T=t-integral(U_ext+v²/2)dt-v dot(x-z(t))+higher orders. The simultaneity term vanishes on the reference worldline, not on an entire finite mouth surface',
            scope='Matched local RN coordinate time, not central proper time; self-redshift sqrt(F) stays separate. Acceleration/surface lapse gradients and rotating shifts still require finite-surface matching. Identical synchronized exchange-symmetric mouths share the rate, so this motion creates no relative time offset. Applying the displayed source factor is only the isolated common-clock correction, not a full PN quantum-state update or a global chronology theorem'),
        orbit_equations='Omega*re=2/D^(3/2), Eorb*G/re=-1/D. Leading EM: GP=32/(3D^4). Standard mass quadrupole: GP_GW=512/(5D^5). Selected diagnostic dD/d(t/re)=-32/(3D²)*(1+a/D), a=48/5 for this combination',
        radiation_prefactor_audit='For point masses m at positions plus/minus(d/2)(cos Omega*t,sin Omega*t,0), sum(I_ij_tripleprime²)=8m²d^4Omega^6; standard P=G*sum/5 gives8/5 rather than MMP printed2/15 in Eq6.53. Factor12 difference is retained explicitly, not silently called the paper formula',
        orbit_scope='Fixed-mass/charge, circular point-orbit and leading far-zone radiation diagnostic over the initial pi*ell traversal scale. Leading GW plus leading EM omits other relative1/D conservative/EM corrections; displayed digits are numerical model precision, not physical precision. Radiation capture, rotating quantum state, induced fields/currents, mode mixing and throat response are not solved. Omega*ell is a scale-separation diagnostic, not an error or stability bound. Tiny radial drift does not imply fixed orientation',
        verification=dict(reconstructed_uniform_field_and_acceleration=True,residual_gauge_decay=True,
            source_Lagrange_identity=True,affine_boundary_propagation_identity=True,
            mass_quadrupole_contraction=True,drift_primitive_identity=True),
        scope='One leading classical incident combination is now determined. No physical source creation, completed rotating exterior/throat solution, nonlinear stability, payload passage, time offset or13-to12Sep2027 route is established',
        sources=['https://arxiv.org/html/1807.04726v3 (Sec6.1 orbit and EM power; printed GW coefficient distinguished explicitly)',
                 'https://arxiv.org/pdf/gr-qc/9808044 (Wells Eq3.28 weak-field acceleration/charge/field relation; no static pair inferred)',
                 'https://arxiv.org/pdf/gr-qc/0501041 (Flanagan-Hughes quadrupole wave and GW energy flux; our explicit two-point-mass evaluation)',
                 'https://astro.troja.mff.cuni.cz/davok/papers/iauresol.pdf (Soffel et al. Eq34 local/global time transformation; our leading charged-orbit application)',
                 'Our derivations: magnetic growing-amplitude normalization, source-aware finite-cut relation and dimensionless orbit checks'])


def rotating_source_checks(loop_source, optical, companion):
    """Rotating Ward obstruction and the published weak-buffer screening equation."""
    from scipy.integrate import solve_bvp

    t, r, theta, phi = sp.symbols('t r theta phi', real=True)
    amplitude, omega = sp.symbols('A omega', real=True)
    f = sp.Function('F')(r)
    coordinates = (t,r,theta,phi)
    metric = (-f,1/f,r*r,r*r*sp.sin(theta)**2)
    volume = r*r*sp.sin(theta)
    harmonic = sp.sin(theta)*sp.cos(phi-omega*t)
    rho = amplitude/(r*r*f)
    stress = sp.zeros(4)
    stress[0,0], stress[1,1] = f*rho*harmonic, rho*harmonic/f
    clean = lambda value: sp.factor(sp.trigsimp(sp.expand_trig(value)))
    # For symmetric covariant stress on this diagonal metric, this density
    # identity avoids constructing a second general-purpose geometry helper.
    def divergence(tensor):
        return [clean(sum(sp.diff(volume*tensor[a,b]/metric[a],coordinates[a])
                          for a in range(4))/volume-sum(tensor[a,a]*
                          sp.diff(metric[a],coordinates[b])/metric[a]**2
                          for a in range(4))/2) for b in range(4)]
    naive = divergence(stress)
    assert clean(naive[0]+rho*sp.diff(harmonic,t)) == 0
    assert naive[0] != 0 and naive[1:] == [0,0,0]
    for a in (2,3):
        stress[0,a] = stress[a,0] = -amplitude*sp.diff(harmonic,t,coordinates[a])/(2*f)
    carried = divergence(stress)
    expected = [0,0]+[-amplitude*omega**2*sp.diff(harmonic,coordinates[a])/(2*f*f)
                       for a in (2,3)]
    assert all(clean(x-y)==0 for x,y in zip(carried,expected))
    assert clean(sum(stress[a,a]/metric[a] for a in range(4))) == 0
    angular_laplacian = sp.diff(sp.sin(theta)*sp.diff(harmonic,theta),theta)/sp.sin(theta)
    angular_laplacian += sp.diff(harmonic,phi,2)/sp.sin(theta)**2
    assert clean(angular_laplacian+2*harmonic) == 0

    # MMP G.38 is a weak, corotating exterior equation. Keep its explicit
    # gauge coupling: large q alone does not justify deleting the Laplacian.
    charge = loop_source['parameters']['q']
    alpha = optical['parameters']['alpha']
    g_squared = 4*math.pi**2*alpha*charge
    kappa = g_squared*charge/(4*math.pi**2)
    k = sp.symbols('k', positive=True)
    particular = -1/(k+2)  # a/(q*Omega_global), for B=+q*rhat/(2r²).
    radial_operator = lambda u: -sp.diff(r*r*sp.diff(u,r),r)+(k+2)*u
    assert sp.simplify(radial_operator(particular)+1) == 0
    roots = ((-1+sp.sqrt(9+4*k))/2,(-1-sp.sqrt(9+4*k))/2)
    assert all(sp.simplify(s*(s+1)-(k+2))==0 for s in roots)
    # The finite-shell basis exposes arbitrary physical end data, never zeroed.
    # ponytail: monopole flat-buffer equation only; solve the pair and throat
    # matching before interpreting these diagnostic end responses physically.
    r_inner, r_outer = 100.,10000.
    width = math.log(r_outer/r_inner)
    x = np.linspace(0.,width,201)
    beta = math.sqrt(9+4*kappa)/2
    inner_weight = np.exp(-x/2)*np.sinh(beta*(width-x))/math.sinh(beta*width)
    outer_weight = np.exp((width-x)/2)*np.sinh(beta*x)/math.sinh(beta*width)
    assert np.all(inner_weight>=0) and np.all(outer_weight>=0)
    assert np.max(inner_weight+outer_weight) <= 1+1e-14
    def rhs(_,values):
        return np.array([values[1],(kappa+2)*values[0]-values[1],
                         values[3],(kappa+2)*values[2]-values[3]])
    def boundary(left,right):
        return np.array([left[0]-1,right[0],left[2],right[2]-1])
    errors, nodes = [], []
    for tolerance in (1e-8,1e-10):
        mesh = np.linspace(0.,width,41)
        guess = np.array([1-mesh/width,-np.ones_like(mesh)/width,
                          mesh/width,np.ones_like(mesh)/width])
        solution = solve_bvp(rhs,boundary,mesh,guess,tol=tolerance,max_nodes=12000)
        assert solution.success, solution.message
        actual = solution.sol(x)
        error = float(max(np.max(abs(actual[0]-inner_weight)),np.max(abs(actual[2]-outer_weight))))
        assert error < 2e-9
        errors.append(error); nodes.append(int(solution.x.size))
    rows = []
    for value in (0.,1.,float(kappa),100.):
        beta_i = math.sqrt(9+4*value)/2
        midpoint_common = math.sinh(beta_i*width/2)/math.sinh(beta_i*width)
        rows.append(dict(kappa=value,particular_over_q_omega=-1/(2+value),
            screened_to_unscreened_particular_ratio=2/(2+value),
            large_kappa_relative_particular_error=None if value==0 else 2/value,
            radial_exponents=[-.5+beta_i,-.5-beta_i],
            midpoint_inner_boundary_weight=math.exp(-width/4)*midpoint_common,
            midpoint_outer_boundary_weight=math.exp(width/4)*midpoint_common))
    omega_global = companion['orbit']['omega_times_re']
    separation = optical['parameters']['separation_over_re']
    omega_rn = omega_global*(1+3/(2*separation))
    cuts = [dict(radius_over_re=radius,
        angular_flux_amplitude_over_density_amplitude=omega_rn*radius/(2*(1-1/radius)))
        for radius in (1.01,1.1,4.)]
    return dict(
        model='Necessary rotating-source conservation test on fixed RN, plus a separate leading weak-buffer fermion screening calculation; no joined rotating quantum state',
        parameters=dict(q=charge,alpha=alpha,g_squared=g_squared,kappa=kappa,
            omega_global_times_re=omega_global,omega_RN_times_re_with_leading_clock_only=omega_rn,
            diagnostic_buffer_radii_over_re=[r_inner,r_outer]),
        ward=dict(harmonic='Y=sin(theta)*cos(phi-omega*t_RN), omega=Omega_global/C_clock',
            baseline='S_tt=F*rho*Y, S_rr=rho*Y/F, rho=A/(r²F); other covariant components zero',
            naive_divergence_covector=['-A*partial_tY/(r²F)','0','0','0'],
            angular_flux='S_tA=-A*D_A(partial_tY)/(2F); A as an index denotes an angular coordinate, while the scalar A is source amplitude',
            remaining_divergence_covector=['0','0','-A*omega²*partial_thetaY/(2F²)','-A*omega²*partial_phiY/(2F²)'],
            trace=0, angular_flux_amplitude_ratios=cuts,
            numerical_ratio_background='F=(1-1/r)^2 at r=1.01,1.1,4, in r_e units. Ratios use the absolute dipole density amplitude, not its angle-dependent zeros; symbolic identities allow general positive F',
            meaning='Rotating the old source alone fails conservation at first order. The displayed polar energy flux cancels that order and leaves a second-order angular force. It is a necessary-balance example, not a unique physical completion',
            assumptions='Fixed diagonal RN, isolated radial state term, zero Lorentz force. A coordinate rotation of both metric and tensor instead preserves covariance and introduces metric shifts. The full charged source obeys divT=F.j and requires its actual current; a dipole coefficient is not positive matter everywhere'),
        screening=dict(published_equation='(-Laplacian+g²*abs(B)/(2pi²))*A_t=-2*Omega_global dot B (MMP G.38)',
            convention='B=+q*rhat/(2r²); theta measured from the rotation axis. Reverse B to reverse potential/current signs, not the screening ratios. This is the original electromagnetic potential, not the geometric Fcal field',
            radial_equation='A_t=a(r)*cos(theta); -(r²*a_prime)_prime+(2+kappa)*a=-q*Omega_global, kappa=g²*q/(4pi²)=alpha*q²',
            particular='a_p=-q*Omega_global/(2+kappa); homogeneous a_h=C_plus*r^s_plus+C_minus*r^s_minus, s_plus/minus=(-1+/-sqrt(9+4*kappa))/2',
            current='j_paper^t=abs(B)*A_t/(2pi²), in MMP G.37 convention. The earlier action-variation current used nablaF=-g²*j_action, so j_action=-j_paper; no sign substitution between the two is made',
            actual_particular_over_q_omega=-1/(2+kappa),
            screened_to_unscreened_particular_ratio=2/(2+kappa),
            large_kappa_relative_particular_error=2/kappa,
            finite_shell='a(r)=a_p+W_inner(r)*(a_inner-a_p)+W_outer(r)*(a_outer-a_p). Both physical end potentials remain unknown; the basis is their linear response, not assigned end values',
            comparison_cases=rows,
            scope='Exact solution only of the stated leading monopole screening equation. No complete pair, finite RN metric correction, rotating Casimir stress, loop holonomy or throat boundary condition. The constant particular is not a full physical solution or a global boundary condition'),
        verification=dict(direct_four_dimensional_naive_and_flux_divergence=True,
            dipole_angular_laplacian=True,trace_free=True,particular_and_exponents_exact=True,
            independent_collocation_max_basis_errors=errors,collocation_mesh_nodes=nodes),
        plot=dict(radius_over_re=(r_inner*np.exp(x)).tolist(),inner_boundary_weight=inner_weight.tolist(),
                  outer_boundary_weight=outer_weight.tolist()),
        next='Derive the coupled rotating electromagnetic current, stress and quantum state with actual exterior/throat boundary data; screened electric potential alone supplies neither the required stress nor a relative time offset',
        scope='No stable medium, source preparation, payload passage or13-to12Sep2027 route. The forward branch is unchanged in this increment',
        sources=['https://arxiv.org/html/1807.04726v3 (AppendixG G.35-G.40; equations2.3 and6.50-6.51; finite-kappa monopole solution and Ward checks are our derivations)'])


def loop_charge_response_checks(loop_source, optical, rotating):
    """Closed-channel charge response, optical averaging and actual return-path parity."""
    length, potential_mean, total_charge = sp.symbols('L Abar Q', real=True)
    variance = sp.symbols('V', nonnegative=True)
    multiplier = potential_mean-sp.pi*total_charge/length
    routhian = variance/(2*sp.pi)+total_charge*potential_mean-sp.pi*total_charge**2/(2*length)
    assert sp.diff(routhian,total_charge) == multiplier
    # Hamiltonian density pi*Pi²/2+Pi*A; Q=-integral Pi. Work with the
    # paper's current j=-Pi and do not call the Routhian a positive ground energy.
    a, j, mu = sp.symbols('a j mu', real=True)
    density = sp.pi*j*j/2-a*j
    assert sp.simplify(sp.diff(density+mu*j,j).subs(j,(a-mu)/sp.pi)) == 0
    assert sp.simplify(length*(potential_mean-multiplier)/sp.pi-total_charge) == 0

    throat = loop_source['parameters']['throat_optical_length_over_d']
    kappa = rotating['parameters']['kappa']
    charge = loop_source['parameters']['q']
    scale = charge*rotating['parameters']['omega_global_times_re']*optical['parameters']['separation_over_re']
    def sample_path(label,count):
        nodes, weights = np.polynomial.legendre.leggauss(count)
        arcweights, bare, algebraic = [], [], []
        for s,w in zip(nodes,weights):
            radius,z,speed = _magnetic_return_path(label,float(s))
            rp,rm = math.hypot(radius,z-.5),math.hypot(radius,z+.5)
            source = _magnetic_loop_source([radius,0.,z],throat)
            b = source['magnetic']
            # Omega along y, azimuth pi/2. Rotate the saved x-z plane into y-z;
            # the radial magnetic component is unchanged by that azimuth change.
            bare.append(-radius*(1/rp-1/rm)/2)
            algebraic.append(-b[0]/(kappa*source['magnitude']))
            arcweights.append(w*speed)
        return nodes,np.asarray(arcweights),np.asarray(bare),np.asarray(algebraic)

    cases, plots = [], []
    for label in (.5,1.,1.5):
        coarse = sample_path(label,64)
        nodes,weights,bare,algebraic = sample_path(label,128)
        exterior = _magnetic_return_length(label)
        total = exterior+throat
        length_error = abs(float(np.sum(weights))-exterior)
        assert length_error < 3e-13
        assert max(np.max(abs(bare+bare[::-1])),np.max(abs(algebraic+algebraic[::-1]))) < 2e-14
        integrals = [float(weights@values) for values in (bare,algebraic)]
        assert max(abs(value) for value in integrals) < 3e-14
        changes = [abs(float(coarse[1]@values)-integral)
                   for values,integral in zip(coarse[2:],integrals)]

        # The last bin is the EXACT aggregate throat-mean degree of freedom
        # for charge balance. Its unresolved internal variance remains separate.
        # ponytail: fixed-geometry, no charge exchange between channels; a full
        # mixed-channel state requires its own charge-allocation susceptibility.
        w = np.r_[weights,throat]
        test = np.r_[bare,.17]  # An audit input, never a physical throat potential.
        norm = float(w.sum())
        projector = np.eye(len(w))-np.outer(np.ones_like(w),w)/norm
        hessian = np.diag(w)-np.outer(w,w)/norm
        mean = float(w@test/norm)
        predicted = test-mean
        bordered = np.zeros((len(w)+1,len(w)+1))
        bordered[:-1,:-1] = np.diag(w)
        bordered[:-1,-1] = w; bordered[-1,:-1] = w
        solution = np.linalg.solve(bordered,np.r_[w*test,0.])
        current,lagrange = solution[:-1],solution[-1]
        error = float(max(np.max(abs(current-predicted)),abs(lagrange-mean)))
        assert error < 1e-12 and abs(w@current) < 1e-12
        assert np.max(abs(projector@projector-projector)) < 3e-15
        shifted = test+.7
        gauge_error = float(np.max(abs(projector@shifted-projector@test)))
        assert gauge_error < 2e-14
        whitened = np.eye(len(w))-np.outer(np.sqrt(w),np.sqrt(w))/norm
        eigenvalues = np.linalg.eigvalsh(whitened)
        assert abs(eigenvalues[0]) < 2e-14 and np.max(abs(eigenvalues[1:]-1)) < 2e-14
        direction = np.r_[np.sin(nodes),-.2]
        functional = lambda values: float(values@hessian@values/2)
        step = 1e-4
        derivative = (functional(test+step*direction)-functional(test-step*direction))/(2*step)
        predicted_derivative = float(w*current@direction)
        derivative_error = abs(derivative-predicted_derivative)
        assert derivative_error < 1e-10
        cases.append(dict(nu=label,exterior_optical_length_over_d=exterior,
            throat_optical_length_over_d=throat,total_optical_length_over_d=total,
            throat_mean_weight=throat/total,exterior_mean_weight=exterior/total,
            exterior_trial_potential_integrals=integrals,quadrature_integral_changes=changes,
            quadrature_length_error=length_error,
            exterior_normalized_charge_per_unit_throat_mean=-exterior*throat/total,
            opposite_throat_normalized_charge_per_unit_throat_mean=exterior*throat/total,
            constrained_solve_max_error=error,charge_balance_residual=float(abs(w@current)),
            gauge_shift_current_error=gauge_error,response_variation_error=derivative_error,
            whitened_hessian_min_eigenvalue=float(eigenvalues[0]),
            nonzero_whitened_eigenvalue_max_error=float(np.max(abs(eigenvalues[1:]-1)))))
        plots.append(dict(nu=label,prolate_coordinate=nodes.tolist(),
            bare_potential_over_qOmega=bare.tolist(),
            algebraic_screening_potential_over_qOmega=algebraic.tolist()))

    radius,z,_ = _magnetic_return_path(1.,.4)
    point = np.array([radius,0.,z]); up=point-np.array([0.,0.,.5]); down=point+np.array([0.,0.,.5])
    rp,rm = np.linalg.norm(up),np.linalg.norm(down)
    ez=np.array([0.,0.,1.])
    grad_label=ez/rm-down[2]*down/rm**3-ez/rp+up[2]*up/rp**3
    tangent=_magnetic_loop_source(point,throat)['tangent']
    along=float(tangent@grad_label)
    assert abs(along) < 1e-13 and np.linalg.norm(grad_label) > .1
    return dict(
        model='Fixed-charge stationary response of complete independent fermion channels; actual exterior reflection checks retain an unknown throat-average potential',
        parameters=dict(q=charge,kappa=kappa,qOmega_times_d=scale,
            charge_sector='Each retained channel has Q=0; zero winding and compatible A_x=0 holonomy are declared state assumptions',
            geometry='Same prescribed leading paths and optical throat length; no new full-loop metric or quantum state solved'),
        canonical=dict(action='L=(1/2pi)*integral dx[(phi_dot-A_t)^2-(phi_x-A_x)^2]',
            momentum_and_charge='Pi=(phi_dot-A_t)/pi; Q=-integral Pi dx, j_paper=-Pi',
            fixed_charge_response='mu_phase=mean_loop(A_t)-pi*Q/L; j_paper=(A_t-mean_loop(A_t))/pi+Q/L',
            optical_average='For the static zeroth-order worldsheet dx=ds/N; mean_loop(A_t)=integral A_t ds/N / integral ds/N. Physical line charge=(A_t-mu_phase)/(pi*N)',
            response_routhian='Gamma_Q/T=integral(A_t-mean(A_t))^2 dx/(2pi)+Q*mean(A_t)-pi*Q^2/(2L), plus fixed winding/Casimir terms',
            sign='The canonical ground energy has the opposite polarization term. mu_phase is phi_dot, not a charge chemical potential unless its sign convention is supplied',
            susceptibility='At fixed Q and geometry: delta j=(I-P_loop)delta A_t/pi, with P_loop the complete optical average. Positive semidefinite Hessian; constants along a loop are polarization null directions',
            zero_mode='For zero-holonomy NS Dirac fermions, Q=N_R+N_L and J=N_R-N_L have the same parity; Q=J=0 is allowed. A smooth zero-mean potential displaces nonzero oscillators without changing that exact charge sector. phi=mu_phase*t denotes stationary derivatives, not a localized phase eigenstate or a demonstrated preparation protocol'),
        charge_scope='Fixed charge of each independent channel is stronger than fixed total charge. Continuous equilibration with exchange gives a shared multiplier weighted by all channel lengths; a finite compact system additionally needs allowed discrete charge allocations. Mode mixing is not solved',
        electromagnetic_response='In a flat patch, j_paper^t=|B|*(A_t-mu_phase)/(2pi²). For Q_loop=0 the local mass term becomes |B|*(I-P_loop)A_t. In curved space use the covariant current and optical measure, not this flat equation unchanged',
        reflection=dict(exterior='For equal mouths, rho(s) is even, z(s) is odd and optical speed is even. The bare rotating potential and the algebraic screening trial are odd, so their complete exterior integrals vanish',
            sufficient_full_loop_condition='Neutral state plus an actual reflection-symmetric throat/exterior attachment, optical measure and odd A_t give mean_loop(A_t)=0. Only then does the projected response reduce to the published local response on this symmetry sector',
            unresolved='Exterior parity alone does not prove full-loop parity. The symmetry must preserve each loop, not just exchange different loops. At fixed Omega the rotating metric shift is not z-reflection invariant; the present cancellation uses a symmetric zeroth-order optical measure and an odd first-order potential. Throat attachment, mixing, gauge/phase transport and boundary data remain requirements',
            throat_mean_formula='With zero exterior mean: mu_phase=(L_th/L_total)*mean_th(A_t). No throat mean is assigned zero unless the additional full-reflection condition is imposed'),
        cases=cases,
        numerical_units='Lengths in d, trial u=A_t/(q*Omega_global), normalized line current jhat=pi*j_paper/(q*Omega_global), normalized charge Qhat=pi*Q/(q*Omega_global*d). The .17 throat mean and .7 common shift only exercise the numerical audit',
        variational_audit='Nonuniform actual exterior quadrature plus one exact aggregate throat-mean variable. A bordered constrained Hamiltonian solve is compared to the projected response. Unresolved throat variations about its mean add a nonnegative variance; no constant physical throat profile is assumed',
        loop_constant_counterexample=dict(potential='A_t proportional to flux label nu, constant along each ideal loop',
            position_over_d=point.tolist(),gradient_nu_times_d=grad_label.tolist(),
            tangent_dot_gradient_nu_times_d=along,transverse_gradient_norm=float(np.linalg.norm(grad_label)),
            meaning='Zero fixed-channel polarization need not mean zero electromagnetic field: loop-dependent constants have transverse gradients. Only a common allowed constant shift is the tested gauge invariance'),
        plot=plots,
        scope='A conditional state-response constraint and parity test, not fixed physical mouth potentials, a coupled rotating Maxwell/Einstein solution, stability, creation or13-to12Sep2027 travel. The algebraic screened trial remains an approximation, as quantified in rotating_source',
        sources=['https://arxiv.org/html/1807.04726v3 (AppendixG bosonized action and local response; closed fixed-charge projector is our constrained application)',
                 'https://arxiv.org/pdf/cond-mat/9805275 (von Delft and Schoeller section7: finite-size charge sectors and nonzero oscillator basis)',
                 'https://arxiv.org/pdf/1902.06584 (Yao and Fukusumi Eqs28-31: spin/charge-sector distinctions in finite fermion-boson correspondence)',
                 'QUANTUM_COMPLETION.md (prescribed optical loops and unresolved state/holonomy context)'])


def loop_retarded_response_checks(charge_response):
    """Causal channel response and finite-duration control on prescribed optical loops."""
    t,x = sp.symbols('t x', real=True)
    phi,at,ax = (sp.Function(name)(t,x) for name in ('phi','At','Ax'))
    rho,current = (at-sp.diff(phi,t))/sp.pi,(sp.diff(phi,x)-ax)/sp.pi
    field = sp.diff(at,x)-sp.diff(ax,t)
    eom = sp.diff(phi,t,2)-sp.diff(phi,x,2)-sp.diff(at,t)+sp.diff(ax,x)
    assert sp.simplify(sp.pi*(sp.diff(rho,t)+sp.diff(current,x))+eom) == 0
    assert sp.simplify(sp.diff(current,t)+sp.diff(rho,x)-field/sp.pi) == 0
    assert sp.simplify(sp.pi*(sp.diff(rho,t,2)-sp.diff(rho,x,2))+sp.diff(field,x)+sp.diff(eom,t)) == 0
    omega,k,a,b = sp.symbols('omega k a b', nonzero=True, real=True)
    response = k*(k*a+omega*b)/(sp.pi*(k*k-omega*omega))
    assert sp.simplify(-omega*(omega/k)*response+k*response-(k*a+omega*b)/sp.pi) == 0
    assert sp.simplify(response.subs({a:-sp.I*omega,b:sp.I*k})) == 0

    # ponytail: prescribed fixed optical circle and external potential only;
    # coupled Maxwell/metric dynamics can shift the spectrum and control zeros.
    def compact(z):
        distance = (np.asarray(z)-.25+.5)%1-.5
        return np.maximum(1-(distance/.1)**2,0)**6
    def step_response(y,tau):
        left,right = compact(y-tau),compact(y+tau)
        return compact(y)-(left+right)/2,(right-left)/2
    fd_rows = []
    for count in (256,512,1024):
        y=np.arange(count)/count; potential=compact(y)
        final=.2; steps=round(final*count/.4); dt=final/steps
        laplacian=lambda v: count**2*(np.roll(v,1)-2*v+np.roll(v,-1))
        previous=np.zeros(count); present=-.5*dt*dt*laplacian(potential)
        for _ in range(1,steps):
            previous,present=present,2*present-previous+dt*dt*laplacian(present-potential)
        exact,_=step_response(y,final)
        error=float(np.max(abs(present-exact)))
        fd_rows.append(dict(nodes=count,courant=dt*count,max_density_error=error,
                            mean_charge_residual=float(abs(np.mean(present)))))
    assert fd_rows[-1]['max_density_error'] < 2e-4
    assert all(3.9 < fd_rows[i]['max_density_error']/fd_rows[i+1]['max_density_error'] < 4.1 for i in (0,1))
    assert max(row['mean_charge_residual'] for row in fd_rows) < 2e-13
    y=np.arange(4096)/4096; potential=compact(y); mean=float(np.mean(potential))
    initial_energy=float(np.mean((potential-mean)**2)/2)
    energy_errors=[]; charge_errors=[]
    for tau in (.07,.2,.35,.61,1.):
        density,flux=step_response(y,tau)
        energy_errors.append(abs(float(np.mean((density-(potential-mean))**2+flux*flux)/2)-initial_energy))
        charge_errors.append(abs(float(np.mean(density))))
    assert max(energy_errors+charge_errors) < 1e-13
    observer=.65; arrival=.3
    before=[float(step_response(observer,tau)[0]) for tau in (0.,.1,.2,.299)]
    assert before == [0.,0.,0.,0.]
    after=float(step_response(observer,.35)[0]); assert abs(after) > .01
    recurrence=float(np.max(abs(step_response(y,1.)[0]))); assert recurrence < 2e-14

    def raised_amplitude(z):
        if abs(abs(z)-math.pi) < 1e-10:
            return math.pi/4
        return abs(math.pi**2*math.cos(z/2)/(math.pi**2-z*z))
    ramps=[]
    for duration in (.25,.5,1.,2.3,10.3):
        frequency=2*math.pi
        drive=lambda time: (1-math.cos(math.pi*time/duration))/2
        sol=solve_ivp(lambda time,v:[v[1],frequency**2*(drive(time)-v[0])],
                      (0,duration),[0.,0.],method='DOP853',rtol=2e-11,atol=2e-13)
        assert sol.success
        measured=math.hypot(sol.y[0,-1]-1,sol.y[1,-1]/frequency)
        predicted=raised_amplitude(frequency*duration)
        assert abs(measured-predicted) < 1e-10
        ramps.append(dict(duration_over_L=duration,residual_amplitude=measured,
            analytic_amplitude=predicted,amplitude_error=abs(measured-predicted)))

    smoothing=.3
    def integrated_rise(time):
        if time <= 0: return 0.
        if time >= smoothing: return time-smoothing/2
        return time/2-smoothing*math.sin(math.pi*time/smoothing)/(2*math.pi)
    def matched_drive(time):
        return integrated_rise(time)-integrated_rise(time-1)
    def matched_modes(frequencies):
        frequencies=np.asarray(frequencies); count=len(frequencies)
        solution=solve_ivp(lambda time,v:np.r_[v[count:],frequencies**2*(matched_drive(time)-v[:count])],
            (0,1+smoothing),np.zeros(2*count),method='DOP853',rtol=2e-11,atol=2e-13,
            max_step=.02,dense_output=True)
        assert solution.success
        residual=np.hypot(solution.y[:count,-1]-1,solution.y[count:,-1]/frequencies)
        return residual,solution
    harmonics=np.arange(1,9)
    matched_residual,matched_solution=matched_modes(2*math.pi*harmonics)
    assert np.max(matched_residual) < 2e-10
    # Appendix A supplies this flat uniform 1+1 comparison, not the mass or
    # spectrum of the unresolved curved, inhomogeneous four-dimensional device.
    mass=sp.symbols('m',nonnegative=True)
    free=k*k/(sp.pi*(k*k-omega*omega))
    coupled=sp.simplify(free/(1+sp.pi*mass*mass*free/(k*k)))
    assert sp.simplify(coupled-k*k/(sp.pi*(k*k+mass*mass-omega*omega))) == 0
    maxwell=[]
    for mass_Lref in (0.,1.,2*math.pi):
        frequency=math.hypot(2*math.pi,mass_Lref)
        residual,_=matched_modes([frequency])
        predicted=abs(np.sinc(frequency/(2*math.pi)))*raised_amplitude(frequency*smoothing)
        assert abs(residual[0]-predicted) < 2e-10
        maxwell.append(dict(mass_times_Lref=mass_Lref,frequency_times_Lref=frequency,
            static_charge_over_free=4*math.pi**2/frequency**2,
            residual_relative_to_screened_static=float(residual[0]),
            predicted_residual=predicted,amplitude_error=abs(float(residual[0])-predicted)))
    reference=next(row['total_optical_length_over_d'] for row in charge_response['cases'] if row['nu']==1.)
    omega_d=charge_response['parameters']['qOmega_times_d']/charge_response['parameters']['q']
    cases=[]
    for row in charge_response['cases']:
        length=row['total_optical_length_over_d']; ratio=reference/length
        frequencies=2*math.pi*harmonics*ratio
        residual,_=matched_modes(frequencies)
        predicted=np.array([abs(np.sinc(value/(2*math.pi)))*raised_amplitude(value*smoothing) for value in frequencies])
        assert np.max(abs(residual-predicted)) < 2e-10
        drive_ratio=omega_d*length/(2*math.pi)
        cases.append(dict(nu=row['nu'],total_optical_length_over_d=length,
            common_ramp_reference_over_actual_length=ratio,
            residual_amplitudes=residual.tolist(),predicted_amplitudes=predicted.tolist(),
            max_amplitude_error=float(np.max(abs(residual-predicted))),
            omega_times_L=omega_d*length,omega_over_first_mode=drive_ratio,
            conditional_harmonic_response_relative_correction=drive_ratio**2/(1-drive_ratio**2)))
    plot_y=np.linspace(0,1,257)
    plot_t=np.linspace(0,3,301)
    matched_values=matched_solution.sol(np.minimum(plot_t,1+smoothing))
    after_t=np.maximum(plot_t-(1+smoothing),0)
    matched_first=np.where(plot_t<=1+smoothing,matched_values[0],
        1+(matched_solution.y[0,-1]-1)*np.cos(2*math.pi*after_t)
        +matched_solution.y[len(harmonics),-1]*np.sin(2*math.pi*after_t)/(2*math.pi))
    return dict(
        model='Retarded linear response of the same neutral independent massless channel on a prescribed static optical circle; external A_t drive, A_x=0, no damping',
        initial_state='Fixed Q=0 and zero winding/holonomy, with zero initial oscillator means and current. Interpreting matched means as the final displaced oscillator vacuum additionally assumes the initial oscillator vacuum and the isolated quadratic action; covariance, source preparation and the complete four-dimensional quantum state are not computed',
        equations=dict(current='rho=(A_t-phi_t)/pi, J=(phi_x-A_x)/pi; E=A_t,x-A_x,t',
            local='rho_t+J_x=0; J_t+rho_x=E/pi; (dt²-dx²)rho=-E_x/pi',
            fourier='exp(-i omega t+i k x), k_n=2pi n/L, n!=0: rho=k*(k*A_t+omega*A_x)/[pi*(k²-(omega+i0)²)]',
            susceptibility='chi_rho,At=k²/[pi*(k²-(omega+i0)²)]; retarded poles below real axis, not damping',
            zero_mode='Q=int rho dx stays fixed; delta rho_0=0. Uniform A_t is compensated by phi_t. A_x holonomy is held zero; its current zero mode is a different problem',
            step='rho=[f(x)-f(x-t)/2-f(x+t)/2]/pi; J=[f(x+t)-f(x-t)]/(2pi), periodic arguments, zero response before switch-on',
            static='The time average is (f-mean f)/pi. A sudden drive recurs at t=L and does not relax. Adiabatic switching or specially timed finite control can prepare the static response'),
        units='c=hbar=1. x is optical length dx=ds/N, t the same static clock. Numerical y=x/L, tau=t/L, u=A_t/A0, r=pi*rho/A0, h=pi*J/A0. Actual loop comparisons instead use Lref of nu1 as time unit; no seconds or calendar offset inferred',
        propagation=dict(profile='f(y)=max(1-((periodic_distance(y,.25))/.1)^2,0)^6',
            finite_difference=fd_rows,observer_y=observer,earliest_arrival_tau=arrival,
            response_before_arrival=before,response_at_tau_point35=after,
            erroneous_instant_static_response_at_observer=-mean,
            full_loop_density_recurrence_error=recurrence,
            max_charge_error=max(charge_errors),max_excitation_energy_drift=max(energy_errors),
            energy='H-H_static=A0²*L/pi times mean[(r-(u-mean u))²+h²]/2 after the fixed step. Positive conserved displacement energy is not a gravitational stability proof',
            causality='The exact characteristics, not finite-difference tails, establish the support bound. The abrupt switch is an ideal effective-model diagnostic, not a physically UV-complete drive'),
        generic_ramps=dict(profile='g(t)=[1-cos(pi*t/T)]/2 on0..T, then1',
            residual_formula='abs(pi²*cos(k*T/2)/(pi²-(k*T)²)), continuous value pi/4 at abs(k*T)=pi; mode amplitude relative to static nonzero mode',
            first_mode=ramps),
        matched_control=dict(reference_nu=1.,reference_length_over_d=reference,
            smoothing_over_reference_length=smoothing,total_duration_over_reference_length=1+smoothing,
            profile='g(t)=[B(t)-B(t-Lref)]/Lref; B=0 below0; B=t/2-delta*sin(pi*t/delta)/(2pi) on0..delta; B=t-delta/2 above delta. delta=.3Lref; g is C2',
            residual_formula='abs(sinc(k*Lref/2)*pi²*cos(k*delta/2)/(pi²-(k*delta)²)); sinc(z)=sin(z)/z, removable limits retained',
            exact_cancellation='The box factor vanishes at every nonzero k=2pi*n/Lref. All ideal channel oscillator means can reach static values after finite time; numerical check tests n1..8, analytic factor covers all',
            harmonics=harmonics.tolist(),matched_residual_amplitudes=matched_residual.tolist(),
            actual_loop_cases=cases,
            limits='Requires known complete optical spectrum and a prescribed spatial potential/drive. One common timing generally fails for unequal loops. No self-consistent source, quantum-state preparation apparatus, Maxwell solution, mixing, damping or evolving metric is constructed'),
        slow_drive_scope='Numbers use a hypothetical periodic potential at omega=Omega on a fixed loop with free transients absent. The MMP field is stationary in its rotating frame; Omega is not automatically a physical channel driving frequency. Three sampled lengths are not a bound on all channels',
        maxwell_comparison=dict(model='Separate flat uniform1+1 Maxwell theory with n identical massless charged channels, periodic zero-total-charge sector and a prescribed conserved external drive. Not the actual MMP geometry',
            definition='rho is per-flavor induced density. Nonzero-mode Gauss law gives -k²*A_ind=g2²*n*rho; A_total=A_external+A_ind; m²=g2²*n/pi. A_external is the bare imposed source potential, not a clamped total potential',
            susceptibility='rho=k²*A_external/[pi*(k²+m²-(omega+i0)²)]; static factor k²/(k²+m²), charged poles sqrt(k²+m²)',
            cases=maxwell,
            interpretation='Illustrative m*Lref values are not inferred material properties or throat masses. Only the charged collective combination acquires this mass; n-1 flavor-neutral combinations stay massless. Gauge feedback changes the charged spectrum and invalidates its generic free-channel exact-control claim. The normalized oscillator ODE checks residual relative to its screened static value; no damping assumed',
            remaining='Uniform1+1 electric feedback omits actual field-line geometry, transverse modes, variable coupling/lapse, metric response, external driver realization and matching. The zero-mode spatial holonomy is not driven'),
        plot=dict(y=plot_y.tolist(),potential=compact(plot_y).tolist(),
            step_profiles=[dict(tau=tau,density=step_response(plot_y,tau)[0].tolist()) for tau in (0.,.12,.24,.5,1.)],
            tau=plot_t.tolist(),step_first_mode=(1-np.cos(2*math.pi*plot_t)).tolist(),
            matched_first_mode=matched_first.tolist(),matched_drive=[matched_drive(float(time)) for time in plot_t]),
        scope='A causal response and ideal finite-control test, not medium formation, full wormhole stability, a new physics discovery or13-to12Sep2027 travel. Maxwell and throat/exterior matching remain unresolved',
        sources=['https://arxiv.org/html/1807.04726v3 (AppendixG action and AppendixA charged collective mass; retarded closed-loop and control calculations are our application)',
                 'https://arxiv.org/pdf/cond-mat/9805275 (finite-size bosonic oscillator and charge-sector foundation)'])


def throat_charge_response_checks(rotating, optical):
    """Curved spherical charged mode and finite-throat response with open exterior data."""
    import mpmath as mp
    from scipy.integrate import solve_bvp
    from scipy.linalg import eigh_tridiagonal
    from scipy.optimize import brentq

    q=rotating['parameters']['q']; coupling=rotating['parameters']['g_squared']
    kappa=rotating['parameters']['kappa']; ell=optical['parameters']['ell_over_re']
    assert math.isclose(kappa,coupling*q/(4*math.pi**2),rel_tol=2e-15)
    t,x=sp.symbols('t x',real=True)
    number=sp.symbols('q',positive=True); weight=sp.Function('W')(x)
    displacement=sp.Function('D')(t,x)
    density=sp.diff(displacement,x)/number; current=-sp.diff(displacement,t)/number
    wave=sp.diff(displacement,t,2)-sp.diff(displacement,x,2)+number*displacement/(sp.pi*weight)
    assert sp.simplify(sp.diff(density,t)+sp.diff(current,x)) == 0
    assert sp.simplify(number*(sp.diff(current,t)+sp.diff(density,x)-displacement/(sp.pi*weight))+wave) == 0
    energy=sp.pi/(2*number)*(sp.diff(displacement,t)**2+sp.diff(displacement,x)**2)+displacement**2/(2*weight)
    flux=-sp.pi/number*sp.diff(displacement,t)*sp.diff(displacement,x)
    assert sp.simplify(sp.diff(energy,t)+sp.diff(flux,x)-sp.pi/number*sp.diff(displacement,t)*wave) == 0
    z,h=sp.symbols('z h',real=True)
    for n in range(4):
        poly=sp.gegenbauer(n,h,z)
        assert sp.expand((1-z*z)*sp.diff(poly,z,2)-(2*h+1)*z*sp.diff(poly,z)+n*(n+2*h)*poly) == 0

    # ponytail: spherical charged collective mode on the prescribed leading
    # throat; exterior response stays unknown instead of imposing a physical wall.
    def basis(cut, rate_squared=0.):
        a=math.atan(cut)
        def rhs(chi,v):
            potential=kappa/math.cos(chi)**2+rate_squared
            return [v[1],potential*v[0],v[3],potential*v[2]]
        solution=solve_ivp(rhs,(0,a),[1.,0.,0.,1.],method='DOP853',
                           rtol=2e-12,atol=2e-14,dense_output=True)
        assert solution.success
        C,Cp,S,Sp=solution.y[:,-1]
        diagonal=(Cp/C+Sp/S)/2
        off=-1/(2*C*S)  # Exact unit Wronskian avoids destructive subtraction.
        matrix=np.array([[diagonal,off],[off,diagonal]])
        return solution,matrix

    rows=[]
    for cut in [4.]+list(optical['parameters']['inner_rho_values']):
        solution,matrix=basis(cut); a=math.atan(cut)
        C,Cp,S,Sp=solution.y[:,-1]
        with mp.workdps(40):
            hh=(1+mp.sqrt(1+4*mp.mpf(str(kappa))))/2
            def even(chi):
                return mp.cos(chi)**hh*mp.hyp2f1(hh/2,hh/2,mp.mpf('.5'),mp.sin(chi)**2)
            def odd(chi):
                return mp.cos(chi)**hh*mp.sin(chi)*mp.hyp2f1((hh+1)/2,(hh+1)/2,mp.mpf('1.5'),mp.sin(chi)**2)
            aa=mp.mpf(a)
            reference=np.array([float(even(aa)),float(mp.diff(even,aa)),float(odd(aa)),float(mp.diff(odd,aa))])
            exact_off=float(-1/(2*even(aa)*odd(aa)))
        error=float(np.max(abs(solution.y[:,-1]/reference-1)))
        assert error < 1e-10 and abs(matrix[0,1]/exact_off-1) < 2e-10
        assert np.linalg.eigvalsh(matrix)[0] > 0
        rows.append(dict(rho_cut=float(cut),chi_cut=a,optical_half_length_over_re=ell*a,
            boundary_matrix=matrix.tolist(),even_eigenvalue=Cp/C,odd_eigenvalue=Sp/S,
            common_endpoint_displacement_to_center=1/C,
            hyp2f1_max_relative_error=error,off_diagonal_relative_error=abs(matrix[0,1]/exact_off-1),
            endpoint_lapse_over_center=math.sqrt(1+cut*cut),
            leading_overlap_radius_fraction=cut/ell))

    solution,matrix=basis(4.); a=math.atan(4.)
    C,Cp,S,Sp=solution.y[:,-1]
    energy_errors=[]
    for parity,normalization in ((0,C),(2,S)):
        integral=2*quad(lambda chi:(solution.sol(chi)[parity+1]**2+
            kappa/math.cos(chi)**2*solution.sol(chi)[parity]**2)/normalization**2,
            0,a,epsabs=1e-10,epsrel=3e-11)[0]
        values=np.array([1.,1.]) if parity==0 else np.array([-1.,1.])
        energy_errors.append(abs(integral-float(values@matrix@values)))
    assert max(energy_errors) < 1e-9
    mesh=np.linspace(-a,a,81)
    bvp=solve_bvp(lambda chi,v:np.vstack([v[1],kappa/np.cos(chi)**2*v[0]]),
        lambda left,right:np.array([left[0],right[0]-1]),mesh,
        np.vstack([(mesh+a)/(2*a),np.full_like(mesh,1/(2*a))]),tol=1e-9,max_nodes=10000)
    assert bvp.success
    endpoint_error=float(np.max(abs(np.array([-bvp.y[1,0],bvp.y[1,-1]])-matrix[:,1])))
    assert endpoint_error < 1e-9
    sample=np.linspace(-a,a,241)
    fundamental=solution.sol(abs(sample))
    normalized=(fundamental[0]/C+np.sign(sample)*fundamental[2]/S)/2
    bvp_error=float(np.max(abs(bvp.sol(sample)[0]-normalized)))
    assert bvp_error < 1e-9
    laplace=[]
    for rate in (0.,.5,1.,2.):
        _,kernel=basis(4.,rate*rate)
        assert np.linalg.eigvalsh(kernel)[0] > 0
        laplace.append(dict(rate_times_ell=rate,boundary_matrix=kernel.tolist(),
                            parity_eigenvalues=np.linalg.eigvalsh(kernel).tolist()))

    hh=(1+math.sqrt(1+4*kappa))/2
    roots=[brentq(lambda frequency:basis(4.,-frequency*frequency)[0].y[2*(n%2),-1],
                   hh+n-.01,hh+n+.1,xtol=5e-13) for n in range(4)]
    grids=[]
    for count in (256,512,1024):
        step=2*a/(count+1); interior=np.linspace(-a+step,a-step,count)
        eigenvalues=eigh_tridiagonal(2/step**2+kappa/np.cos(interior)**2,
            np.full(count-1,-1/step**2),select='i',select_range=(0,3),eigvals_only=True)
        frequencies=np.sqrt(eigenvalues)
        grids.append(dict(nodes=count,frequencies_times_ell=frequencies.tolist(),
                           max_shooting_difference=float(np.max(abs(frequencies-roots)))))
    assert grids[-1]['max_shooting_difference'] < 3.2e-5
    assert all(3.9 < grids[i]['max_shooting_difference']/grids[i+1]['max_shooting_difference'] < 4.1 for i in (0,1))
    return dict(
        model='Gauge-invariant spherical electric displacement in the charged collective lowest-Landau-level sector on a prescribed leading curved throat. This angular l=0 mode does not solve the rotating l=1 source',
        parameters=dict(q=q,g_squared=coupling,kappa=kappa,ell_over_re=ell,ads_scaling_dimension=hh,
            cuts='rho4 is a retained interior diagnostic; rho10,100,1000 reuse saved overlap cuts. None is an established physical junction'),
        reduction=dict(metric='ds4²=N²(-dt²+dx²)+R²dOmega²',
            action='S=int dt dx[W*E²/2+q*((phi_t-A_t)²-(phi_x-A_x)²)/(2pi)], W=4piR²/(g²N²), E=A_t,x-A_x,t',
            displacement='D=W*E=4piR²*E_orthonormal/g², rho=(A_t-phi_t)/pi per flavor, J=(phi_x-A_x)/pi; D_x=q*rho,D_t=-q*J',
            wave='D_tt-D_xx+q*D/(pi*W)=0 on any fixed static spherical geometry in this truncation',
            throat='R=re, N=re/(ell*cos chi), x=ell*chi, tau=t/ell: D_tautau-D_chichi+kappa*sec²chi*D=0',
            charged_only='q-1 flavor-neutral combinations remain massless. The charged collective restoring term is fixed by the same saved q and g, not an adjustable new material mass'),
        energy=dict(global_time='H_t=pi/(2q)*int dx[D_t²+D_x²+q*D²/(pi*W)]',
            scaled='H_tau=ell*H_t=pi/(2q)*int dchi[D_tau²+D_chi²+kappa*sec²chi*D²]',
            flux='dH_tau/dtau=(pi/q)*[D_tau*D_chi] at right minus left. Unknown exterior data can exchange energy',
            static_identity='H_tau=pi/(2q)*b^T*K(0)*b, b=(D_left,D_right); applies to charged perturbation energy, excludes background Casimir and metric response',
            unit_amplitude_integral_errors=energy_errors,
            limit='Positive interior energy rules out growth for this fixed operator under homogeneous reflecting/passive conditions, not for an unspecified active exterior or the complete geometry'),
        matching=dict(convention='K maps b=(D(-a),D(a)) to outward chi derivatives (-D_chi(-a),D_chi(a)); physical optical derivative divides by ell, proper-normal derivative multiplies by cos(a)/re',
            formula='C(0)=1,C_prime(0)=0,S(0)=0,S_prime(0)=1. Kdiag=(C_prime/C+S_prime/S)/2, Koff=-1/(2*C*S). Wronskian C*S_prime-C_prime*S=1',
            electric_potentials='For a static solution, A_t-mu_phase=pi*D_x/q; in the common clock A_R-mu=pi*(K*b)_R/(q*ell), A_L-mu=-pi*(K*b)_L/(q*ell). Endpoint D values are not voltages',
            charge='Total throat-segment charge=D_right-D_left. Full-loop neutrality allows compensation in the exterior and does not set the odd throat response to zero',
            join='For a smooth source-free continuation within the same spherical reduced model, D and D_x match. With aligned field variables/outward normals, (K_throat+K_exterior)*b=f_external after common-clock conversion. K_exterior and f_external are not assigned',
            caveat='Actual non-spherical returning fields can mix angular/channels and require the full gauge/metric junction conditions. This2x2 map is the throat l=0 contribution, not a completed exterior relation'),
        static_cases=rows,
        independent_bvp=dict(rho_cut=4.,endpoint_basis=[0.,1.],mesh_nodes=int(len(bvp.x)),
                             normal_derivative_max_error=endpoint_error,profile_max_error=bvp_error),
        laplace_cases=laplace,
        laplace_scope='rate=s_global*ell>=0; solutions obey D_chichi=(kappa*sec²chi+rate²)D. These samples are not an inverse Laplace kernel, a transmission coefficient or a completed causal exterior',
        reflecting_benchmark=dict(scope='Mathematical homogeneous Dirichlet walls at rho4, solely to verify the interior operator; not physical mouth conditions or device frequencies',
            finite_cut_shooting_frequencies_times_ell=roots,finite_difference_grids=grids,
            full_ads='For the separate full interval (-pi/2,pi/2) with finite-energy reflecting endpoints, D_n=cos^h(chi)*Gegenbauer_n^h(sin chi), h(h-1)=kappa, omega_n*ell=n+h. Polynomial identity verified for n0..3',
            full_ads_frequencies_times_ell=[hh+n for n in range(4)]),
        plot=dict(chi=sample.tolist(),potential=(kappa/np.cos(sample)**2).tolist(),
                  same_end_displacement=(fundamental[0]/C).tolist(),
                  opposite_end_displacement=(np.sign(sample)*fundamental[2]/S).tolist()),
        scope='Actual leading throat redshift and fixed coupling are included, but creation, nonlinear/metric feedback, source apparatus, rotating multipoles, full-loop quantum state and exterior attachment remain unproved. No stable medium, novelty claim or13-to12Sep2027 travel',
        sources=['https://arxiv.org/html/1807.04726v3 (AppendixA spherical charged collective action; finite-cut response and verification are our application)',
                 'QUANTUM_COMPLETION.md (spherical mode reduction and need for physical junction response)'])


def charged_mouth_matching_checks(throat, optical):
    """Join the spherical throat response to finite RN tails without closing the exterior."""
    import mpmath as mp
    z,p=sp.symbols('z p',positive=True); radius=1+z
    F=z*z/radius**2
    growing=z**p*(radius+1/p); decaying=z**(-p-1)*(radius-1/(p+1))
    pg=(p+1)*z**(p+1)/radius; pd=-p*z**(-p)/radius
    for field,flux in ((growing,pg),(decaying,pd)):
        assert sp.simplify(F*sp.diff(field,z)-flux) == 0
        assert sp.simplify(sp.diff(flux,z)-p*(p+1)*field/radius**2) == 0
    assert sp.simplify(growing*pd-decaying*pg+2*p+1) == 0
    r,M,kk=sp.symbols('r M kappa',positive=True); f=sp.Function('f')(r)
    delta=r*r-2*M*r+1
    field=kk*r*f-delta*sp.diff(f,r); flux=kk*delta*sp.diff(f,r)/r
    second=(kk*f-sp.diff(delta,r)*sp.diff(f,r))/delta
    assert sp.simplify((delta/r**2*sp.diff(field,r)-flux).subs(sp.diff(f,r,2),second)) == 0
    assert sp.simplify((sp.diff(flux,r)-kk*field/r**2).subs(sp.diff(f,r,2),second)) == 0
    kappa=throat['parameters']['kappa']; ell=optical['parameters']['ell_over_re']
    b2=-2*optical['parameters']['geometric_mass_minus_re_over_re']
    assert math.isclose(b2*ell*ell,1.,rel_tol=3e-15)
    outer=optical['parameters']['outer_areal_radii_over_re'][0]
    clock=1-3/(2*optical['parameters']['separation_over_re'])
    # ponytail: finite spherical RN tails and leading overlap embedding only;
    # arbitrary outer data remain. This is not the nonspherical return region.
    def tail(cut, mass_term):
        ya=cut/ell; yb=outer-1; ua,ub=math.log(ya),math.log(yb)
        def rhs(u,v):
            y=math.exp(u); r=1+y
            A=r*r/(1+mass_term*r/y**2)  # y²/F, without subtracting near-equal masses.
            return [A*v[1],kappa*v[0]/r**2-v[1]]
        forward=solve_ivp(rhs,(ua,ub),[0.,1.],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        backward=solve_ivp(rhs,(ub,ua),[0.,-1.],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        assert forward.success and backward.success
        df,qf=forward.y[:,-1]; db,qb=backward.y[:,-1]
        off=-ya/df; reverse_off=-yb/db
        matrix=np.array([[-ya*qb/db,off],[off,yb*qf/df]])
        reciprocity=abs(off/reverse_off-1)
        assert reciprocity < 3e-10 and np.linalg.eigvalsh(matrix)[0] > 0
        return matrix,reciprocity,forward,backward

    def analytic_reference(cut, mass_term):
        with mp.workdps(70):
            h=(1+mp.sqrt(1+4*mp.mpf(str(kappa))))/2
            kk=mp.mpf(str(kappa)); mass=1-mp.mpf(str(mass_term))/2
            aa=1-mass*mass
            def fundamental(y):
                r=1+y
                delta=(r-mass)**2+aa; fields=[]; fluxes=[]
                for power in (h-1,-h):
                    def f(radial):
                        zz=radial-mass
                        return zz**power*mp.hyp2f1(-power/2,(1-power)/2,mp.mpf('.5')-power,-aa/zz**2)
                    value=f(r); derivative=mp.diff(f,r)
                    fields.append((kk*r*value-delta*derivative)/(kk-power))
                    fluxes.append(kk*delta*derivative/(r*(kk-power)))
                return (*fields,*fluxes)
            gi,di,pgi,pdi=fundamental(mp.mpf(str(cut))/mp.mpf(str(ell)))
            go,do,pgo,pdo=fundamental(mp.mpf(str(outer))-1)
            denominator=-go*di+do*gi; wronskian=1-2*h
            return np.array([[float((go*pdi-do*pgi)/denominator),float(-wronskian/denominator)],
                [float(-wronskian/denominator),float((-pgo*di+pdo*gi)/denominator)]])

    def joined(matrix,record,clock_factor=1.):
        kt=np.asarray(record['boundary_matrix'])*clock_factor/ell
        ka=matrix*clock_factor
        even=record['even_eigenvalue']*clock_factor/ell
        odd=record['odd_eigenvalue']*clock_factor/ell
        s_even,s_odd=ka[0,0]+even,ka[0,0]+odd
        inner_even=-ka[0,1]/s_even
        center=inner_even*record['common_endpoint_displacement_to_center']
        correction=ka[0,1]**2*(ka[0,0]+kt[0,0])/(s_even*s_odd)
        # Do not subtract almost equal parity eigenvalues: coupling is ~1e-75.
        off=ka[0,1]**2*kt[0,1]/(s_even*s_odd)
        result=np.array([[ka[1,1]-correction,off],[off,ka[1,1]-correction]])
        inner=np.linalg.solve(ka[0,0]*np.eye(2)+kt,-ka[0,1]*np.ones(2))
        assert np.max(abs(inner/inner_even-1)) < 2e-14
        return result,center,inner_even,correction

    records=[]; plot=None; reference_case=None
    for cut in optical['parameters']['inner_rho_values']:
        saved=next(row for row in throat['static_cases'] if row['rho_cut']==cut)
        corrected,reciprocity,forward,backward=tail(cut,b2)
        extreme,ext_reciprocity,_,_=tail(cut,0.)
        exact=analytic_reference(cut,0.)
        exact_error=float(np.max(abs(extreme/exact-1)))
        assert exact_error < 2e-10
        corrected_exact=analytic_reference(cut,b2)
        corrected_error=float(np.max(abs(corrected/corrected_exact-1)))
        assert corrected_error < 2e-10
        outer_matrix,center,inner,diagonal_correction=joined(corrected,saved)
        _,extreme_center,_,_=joined(extreme,saved)
        transformed,changed_center,_,_=joined(corrected,saved,clock)
        covariance=max(abs(changed_center/center-1),float(np.max(abs(transformed/(clock*outer_matrix)-1))))
        assert covariance < 2e-14
        ya=cut/ell; yb=outer-1
        # Normalized one-ended solutions independently check the radial energy identity.
        energy_errors=[]
        for solution,normalization,expected in ((forward,forward.y[0,-1],corrected[1,1]),
                                               (backward,backward.y[0,-1],corrected[0,0])):
            def integrand(u):
                y=math.exp(u); r=1+y; v=solution.sol(u)/normalization
                A=r*r/(1+b2*r/y**2)
                return y*(A*v[1]**2+kappa*v[0]**2/r**2)
            integral=quad(integrand,math.log(ya),math.log(yb),epsabs=1e-24,epsrel=1e-10,limit=120)[0]
            energy_errors.append(abs(integral/expected-1))
        assert max(energy_errors) < 2e-10
        records.append(dict(rho_cut=cut,inner_radius_minus_re_over_re=ya,
            relative_F_mass_correction=b2*(1+ya)/ya**2,
            first_order_sphere_profile_relative_to_linear_embedding=2*(1-cut*math.atan(1/cut))/(math.pi*cut),
            corrected_tail_matrix=corrected.tolist(),extremal_tail_matrix=extreme.tolist(),
            reciprocity_relative_error=reciprocity,extremal_reciprocity_relative_error=ext_reciprocity,
            extremal_closed_basis_max_relative_error=exact_error,
            corrected_hypergeometric_max_relative_error=corrected_error,
            radial_energy_identity_relative_errors=energy_errors,
            corrected_center_per_unit_outer_displacement=center,
            extremal_center_per_unit_outer_displacement=extreme_center,
            extremal_omission_fraction=extreme_center/center-1,
            inner_displacement_per_unit_outer_displacement=inner,
            composite_outer_matrix=outer_matrix.tolist(),outer_diagonal_reduction=diagonal_correction,
            clock_covariance_relative_error=covariance))
        if cut==100:
            reference_case=dict(rho_cut=cut,independent_method='70-digit hypergeometric endpoints versus DOP853',
                                max_matrix_relative_error=corrected_error)
            u=np.linspace(math.log(ya),math.log(yb),241)
            profile=inner*backward.sol(u)[0]/backward.y[0,-1]+forward.sol(u)[0]/forward.y[0,-1]
            assert math.isclose(profile[0],inner,rel_tol=2e-12) and math.isclose(profile[-1],1.,rel_tol=2e-12)
            plot=dict(radial_offset_over_re=np.exp(u).tolist(),unit_outer_static_displacement=profile.tolist(),
                      rho_cut=cut,inner_displacement=inner,center_displacement=center)
    centers=np.array([row['corrected_center_per_unit_outer_displacement'] for row in records])
    spread=float(np.max(centers)/np.min(centers)-1)
    assert spread < 2e-9
    return dict(
        model='Conditional static spherical charged-mode continuation through two finite RN regions, joined to the saved leading throat response. Arbitrary outer endpoint data remain; no global source/metric solution',
        parameters=dict(kappa=kappa,ell_RN_over_re=ell,mass_deficit_b_squared=b2,
            outer_areal_radius_over_re=outer,outer_radius_over_separation=outer/optical['parameters']['separation_over_re'],
            baseline_clock_factor=1.,test_clock_factor=clock,
            embedding='re=1; inner r=1+rho_cut/ell_RN is the leading overlap relation. It is not an exact physical finite junction'),
        equations=dict(metric='F=[(r-1)²+b²*r]/r², b²=-2*(M-1)=1/ell_RN². Evaluate this form, not1-2*M/r+1/r² in double precision',
            static='(F*D_r)_r-kappa*D/r²=0, P=F*D_r; radial state (D,P) has trace-zero transfer generator and exact transfer determinant1',
            stable='y=r-1,u=log y,Q=P/y: D_u=(y²/F)*Q, Q_u=kappa*D/r²-Q',
            exact_extremal='For b²=0: D_p=(r-1)^p*(r+1/p), P_p=(p+1)*(r-1)^(p+1)/r, p=h-1 or-h, h(h-1)=kappa. Ordered basis determinant=-(2h-1), checked symbolically',
            exact_corrected='M=1-b²/2, Delta=r²-2Mr+1, z=r-M, a²=1-M². f_p=z^p*2F1(-p/2,(1-p)/2;1/2-p;-a²/z²) solves (Delta f_p_prime)_prime=kappa*f_p. D_p=(kappa*r*f_p-Delta*f_p_prime)/(kappa-p), P_p=kappa*Delta*f_p_prime/[r*(kappa-p)], p=h-1,-h. These span both solutions; neither branch coefficient is discarded',
            energy='int_a^b [F*D_r²+kappa*D²/r²]dr=[D*P]_a^b. Multiply by pi/(2q) for static RN-clock charged energy'),
        matching=dict(tail_map='Ktail maps (D_inner,D_outer) to (-P_inner,P_outer). Two-ended integrations give its diagonals; reciprocity gives Koff=-y_inner/D_forward=-y_outer/D_backward',
            schur='Kth=Kchi/ell_RN; S=Ktail_aa*I+Kth; b_inner=-Ktail_ab*S^-1*b_outer; Kouter=Ktail_bb*I-Ktail_ab²*S^-1',
            orientation='All D values are aligned along one direction through the throat. Outward radial D at the left mouth has the opposite sign; at the right it has the same sign. Equal through-oriented outer values correspond to opposite outward radial fluxes',
            unknown_data='No outer endpoint value is selected physically. The displayed center/profile uses the unit basis test b_outer=(1,1), neutral for the finite composite; actual source and quantum state are not supplied',
            numerical='Keep the cross-mouth off-diagonal directly. It is lost by subtracting almost equal eigenvalues; no numerical determinant is inferred from subtraction of huge transfer-matrix products'),
        clock=dict(convention='Pure relabel t_RN=C*t_common with ell_common=ell_RN/C multiplies every throat/tail/composite derivative map by C. Endpoint and center displacements are invariant',
            limit='The saved geometry remains unchanged. The chosen C only tests covariance; it does not implement physical companion/rotation corrections, which also require new geometry and consistent embedding'),
        cases=records,independent_corrected_check=reference_case,
        corrected_center_relative_cut_spread=spread,
        plot=plot,
        interpretation='Retaining the saved below-extremal mass correction removes an artificial overlap dependence in this charged static response. The tiny center displacement is a normalized electrostatic screening response, NOT a particle-transmission probability, time-machine performance or physical field-strength prediction',
        embedding_limit='The first-order sphere shape 1+rho*atan(rho) differs from its asymptotic linear pi*rho/2 by the recorded relative amount. That correction is not applied alone: a consistent higher-order join also needs the lapse and radial metric corrections. Small response cutoff spread therefore does not establish physical accuracy',
        scope='Static angular l=0 only, fixed approximate geometry and leading finite-cut stitching. Cutoff agreement is not a bound on omitted metric, quantum source, angular mixing, finite separation or junction errors. The nonspherical return region, dynamic kernel, stress feedback and source preparation remain open; no stable medium or13-to12Sep2027 travel',
        sources=['https://arxiv.org/html/1807.04726v3 (AppendixA charged collective reduction; Eqs5.39-5.42 leading throat/RN overlap and negative mass correction. Transfer basis and composed response are our application)',
                 'throat_charge_response (saved curved spherical throat boundary response)'])


def optical_matching_checks(paths):
    """Cutoff sensitivity of leading throat/RN/weak-field optical pieces."""
    import mpmath as mp
    from scipy.optimize import brentq

    with mp.workdps(50):
        alpha = mp.mpf('1e-9')
        ell_v = 4/(mp.pi*alpha)  # re=1; geometric lengths, c=hbar=1.
        branch = paths['matched_cases'][0]
        assert branch['separation_over_ell_v'] == .1
        ell = ell_v*mp.mpf(str(branch['ell_over_ell_v']))
        distance = ell_v/10
        epsilon = -1/(2*ell**2)
        mass = 1+epsilon
        b2 = -epsilon*(2+epsilon)
        b = mp.sqrt(b2)

        def mouth(cutoff_minus_one, outer):
            v1, v2 = cutoff_minus_one-epsilon, outer-mass
            angle = mp.atan2(b*(v2-v1), b2+v2*v1)
            return (v2-v1+mass*mp.log((v2*v2+b2)/(v1*v1+b2))
                    +(2*mass*mass-1)*angle/b)

        exterior, junctions, primitive_errors, flux_root_errors = [], [], [], []
        for outer in map(mp.mpf, (100, 1000, 10000)):
            v = outer-mass
            isotropic = (v+mp.sqrt(v*v+b2))/2
            assert abs(isotropic+mass-b2/(4*isotropic)-outer) < mp.mpf('1e-44')
            xi = float(isotropic/distance)
            # Each ideal point-field cap is almost radial; bound the omitted excess.
            other_over_self = (isotropic/(distance-isotropic))**2
            slope_bound = other_over_self/(1-other_over_self)
            cap_excess_bound = 2*isotropic*slope_bound**2/(mp.sqrt(1+slope_bound**2)+1)
            for nu in (.5, 1., 1.5, 2.):
                root_fn = lambda t: t*(2-xi*t)*(1-xi*t+2*xi)-2*nu*(1-xi*t+xi)
                t = 2. if nu == 2. else brentq(root_fn, 0., 2., xtol=1e-14)
                flux_root_errors.append(abs(root_fn(t)))
                sc = 1-xi*t
                def integrand(s):
                    if nu == 2:
                        rho, z, speed = 0., s/2, .5
                    else:
                        rho, z, speed = _magnetic_return_path(nu, s)
                    return speed*(1/math.hypot(rho, z+.5)+1/math.hypot(rho, z-.5))
                first = 4*float(mass)*quad(integrand, 0., sc,
                    epsabs=1e-7, epsrel=1e-10, limit=200)[0]
                refined = 4*float(mass)*quad(integrand, 0., sc,
                    epsabs=1e-9, epsrel=2e-12, limit=200)[0]
                assert abs(first-refined) < 1e-7
                axis_error = None
                if nu == 2:
                    exact = 4*mass*mp.log((distance-isotropic)/isotropic)
                    axis_error = float(abs(mp.mpf(str(refined))-exact))
                    assert axis_error < 2e-8
                exterior.append(dict(outer_areal_radius_over_re=float(outer), nu=nu,
                    outer_isotropic_radius_over_re=float(isotropic), prolate_endpoint=sc,
                    fixed_path_shapiro_length_over_re=refined,
                    quadrature_refinement_change_over_re=abs(first-refined),
                    axis_analytic_difference_over_re=axis_error,
                    two_flat_caps_radial_length_over_re=float(2*isotropic),
                    two_flat_caps_excess_upper_bound_over_re=float(cap_excess_bound)))
            for rhoc in map(mp.mpf, (10, 100, 1000)):
                lower = rhoc/ell
                travel = mouth(lower, outer)
                # Independent logarithmic-coordinate integral spans both overlap scales.
                integral = mp.quad(lambda y: (1+mp.exp(y))**2*mp.exp(y)/
                    (mp.exp(2*y)+mp.exp(y)/ell**2+1/ell**2),
                    [mp.log(lower), 0, mp.log(outer-1)])
                primitive_errors.append(abs(travel/integral-1))
                assert primitive_errors[-1] < mp.mpf('1e-40')
                # Leading sphere relation only; no uncomputed sqrt(1+phi) resummation.
                alternative = 2*(1+rhoc*mp.atan(rhoc))/(mp.pi*ell)
                map_change = 2*(mouth(alternative, outer)-travel)
                residual = 2*ell*mp.atan(rhoc)+2*travel-mp.pi*ell-2*isotropic
                shapiro = next(row['fixed_path_shapiro_length_over_re'] for row in exterior
                               if row['nu'] == 1 and row['outer_areal_radius_over_re'] == float(outer))
                junctions.append(dict(inner_throat_rho=float(rhoc),
                    outer_areal_radius_over_re=float(outer),
                    leading_inner_radius_minus_re_over_re=float(lower),
                    leading_flat_length_residual_over_re=float(residual),
                    leading_flat_length_residual_over_ell=float(residual/ell),
                    fixed_path_weak_field_length_residual_over_re=float(residual)+shapiro,
                    fixed_path_weak_field_length_residual_over_ell=(float(residual)+shapiro)/float(ell),
                    finite_rho_sphere_map_change_over_ell=float(map_change/ell),
                    scope='nu1 exterior for the weak-field column; radial cap replacement has its separately bounded ideal-point-field error. The altered sphere map is sensitivity only, not an improved joined solution'))
        representative = [row for row in junctions if row['inner_throat_rho'] == 100]
        raw_span = np.ptp([row['leading_flat_length_residual_over_re'] for row in representative])
        corrected_span = np.ptp([row['fixed_path_weak_field_length_residual_over_re'] for row in representative])
        assert raw_span > 18 and 0 < corrected_span < .03
        rlo = mp.mpf(str(exterior[0]['outer_isotropic_radius_over_re']))
        rhi = mp.mpf(str(exterior[-1]['outer_isotropic_radius_over_re']))
        omitted_self = 2*mass**2*(1/rlo-1/rhi)
        omitted_companion = -4*mass*(rhi-rlo)/distance
        predicted_span = omitted_self+omitted_companion
        assert abs(float(predicted_span)-corrected_span) < 1e-7
        assert max(flux_root_errors) < 1e-13
        return dict(model='Optical overlap sensitivity: leading AdS throat, isolated RN mouth approximations with mass below the extremal value, and fixed flat magnetic paths, with and without a weak static optical-metric correction. Not a solution of the joined semiclassical equations',
            parameters=dict(alpha=float(alpha), ell_v_over_re=float(ell_v), ell_over_re=float(ell),
                separation_over_re=float(distance), geometric_mass_minus_re_over_re=float(epsilon),
                inner_rho_values=[10,100,1000], outer_areal_radii_over_re=[100,1000,10000],
                weak_field_path_labels=[.5,1.,1.5,2.]),
            rn_formula='F(R)=(1-1/R)^2+1/(ell^2*R) in re1 units; geometric m=1-1/(2ell^2). Primitive I=R+m*log[(R-m)^2+b^2]+(2m^2-1)/b*atan[(R-m)/b], b^2=1-m^2. Definite arctangent computed with atan2 and lower offset formed without subtracting near-equal radii',
            coordinate_matching='R_iso=[R_areal-m+sqrt((R_areal-m)^2+b^2)]/2; equivalently R_areal=R_iso+m-b^2/(4R_iso). Isolated RN transform only; companion corrections are not supplied',
            mass_scope='m is the local isolated RN matching parameter. At finite exterior length the full variational binding energy is not identified with2*(m-1)/G; exterior quantum energy and anomaly allocation matter',
            optical_pieces='L0=2ell*atan(rhoc)+2*integral_Rc^Rm dR/F+Lflat_cut. Rc-1=rhoc/ell. Report residual relative to pi*ell+d*f(nu); approximate Lflat_cut=d*f(nu)-2R_iso with a separately bounded flat-cap excess',
            cap_bound='Inside a matching sphere Riso<d/2, other-field/self-field<=e=[Riso/(d-Riso)]^2. The tangent transverse/radial ratio is at most e/(1-e). Two-cap excess<=2Riso*(sqrt(1+[e/(1-e)]^2)-1), evaluated without cancellation. This bounds only flat point-magnetic-field geometry',
            weak_metric='In static isotropic coordinates ds^2=-(1-2U)dt^2+(1+2U)dl^2, U=m/r1+m/r2. Along the fixed path dt=(1+2U)dl+O(U^2); using only1/N would miss the spatial-metric contribution. DeltaL=16m integral_0^sc u*speed/(u^2-s^2) ds',
            matching_log='Two isolated RN mouths contribute +4m*log(Riso) at large radius; exterior2U contributes -4m*log(Riso). This cancels the leading outer-cutoff logarithm; it does not fix finite missing terms',
            exterior_cases=exterior, matching_cases=junctions,
            outer_cutoff_comparison=dict(inner_rho=100, nu=1,
                raw_span_over_re=float(raw_span), fixed_path_weak_field_span_over_re=float(corrected_span),
                span_reduction_factor=float(raw_span/corrected_span),
                missing_exterior_self_second_order_span_over_re=float(omitted_self),
                missing_mouth_companion_potential_span_over_re=float(omitted_companion),
                predicted_residual_span_over_re=float(predicted_span),
                predicted_vs_computed_span_difference_over_re=abs(float(predicted_span)-float(corrected_span)),
                explanation='Near the matching boundary dL_total/dRiso approximately2m^2/Riso^2-4m/d. The terms identify omitted exterior self-field and companion-mouth contributions. This explains the residual; it neither constructs those corrections consistently nor bounds all missing terms'),
            inner_map_sensitivity='Keeping other pieces fixed, substitute Rc-1=2*(1+rhoc*atan(rhoc))/(pi*ell), the leading spherical area relation. The resulting change is approximately -4/(3pi*rhoc^4) in L/ell at large rhoc. This exposes omitted matching information; it is not a consistent correction with the unchanged throat metric. sqrt(1+phi) would import unjustified higher-order terms',
            verification=dict(maximum_rn_primitive_vs_log_quadrature_relative_error=mp.nstr(max(primitive_errors),8),
                maximum_stable_cutoff_root_residual=max(flux_root_errors)),
            scope='No exact junction, angular metric/source response, companion mouth redshift, perturbed magnetic paths, continued quantum stress, formation, stability or time offset. Cutoff dependence is not a physical error bound. The fixed curved field line is not a stationary free null path, so path deformation can change finite first-order contributions. No corrected Casimir state or past-travel claim follows',
            sources=['https://arxiv.org/pdf/1807.04726 (MMP5.33,5.39-5.41,AppD; leading matching pieces)',
                'https://arxiv.org/pdf/1409.3743 (Soffel-Han section2 Eqs1-3: weak metric and null propagation; local metric used on a fixed magnetic path, not the paper\'s straight-ray endpoint formula)'])


def _weak_magnetic_return_path(nu, s, epsilon=0.):
    h = 1-s*s
    root = np.sqrt((1-2*epsilon)**2*h*h+nu*nu*s*s+4*epsilon*nu*h)
    u = ((1-2*epsilon)*h+root)/nu
    # Rationalizing u-1 avoids cancellation next to the matching spheres.
    um = (2-nu)*h/(nu*(u+1)-2*(1-2*epsilon)*h)
    rho, z = np.sqrt(um*(u+1)*h)/2, u*s/2
    up = s*(nu-2*((1-2*epsilon)*u+2*epsilon))/root
    rp = (u*up*h-s*um*(u+1))/(4*rho)
    zp = (u+s*up)/2
    return rho, z, np.sqrt(rp*rp+zp*zp), u, up


def weak_magnetic_response_checks(optical):
    """First-order point-source Maxwell paths and constant-companion annulus matching."""
    import mpmath as mp
    from scipy.optimize import brentq

    parameters = optical['parameters']
    distance = parameters['separation_over_re']
    # Coordinates in the field calculation are divided by d; epsilon=m/d.
    def fields(rho, z):
        x1, x2 = np.array([rho, z-.5]), np.array([rho, z+.5])
        a, b = 1/np.sqrt(x1@x1), 1/np.sqrt(x2@x2)
        ga, gb = -x1*a**3, -x2*b**3
        return -ga+gb, 2*((a-1)*gb-(b-1)*ga), ga+gb

    def stream(rho, z):
        r1, r2 = np.sqrt(rho*rho+(z-.5)**2), np.sqrt(rho*rho+(z+.5)**2)
        base = (z+.5)/r2-(z-.5)/r1
        cosine = (rho*rho+z*z-.25)/(r1*r2)
        return base, 2*(1-cosine-base)

    curve=_weak_magnetic_return_path

    maxwell, stream_errors, contour_errors, baseline_errors, flux_errors = [], [], [], [], []
    step = 1e-25
    for rho, z in ((.2, -.3), (.7, .1), (1.3, .8), (.04, .54)):
        b0, b1, gu = fields(rho, z)
        dr = np.imag(fields(rho+1j*step, z)[1])/step
        dz = np.imag(fields(rho, z+1j*step)[1])/step
        maxwell.append(abs(dr[0]+b1[0]/rho+dz[1])/max(1., np.linalg.norm(b1)))
        curl = dz[0]-dr[1]
        expected = 2*(gu[1]*b0[0]-gu[0]*b0[1])
        maxwell.append(abs(curl-expected)/max(1., abs(expected)))
        dpr = np.imag(stream(rho+1j*step, z)[1])/step
        dpz = np.imag(stream(rho, z+1j*step)[1])/step
        stream_errors.append(float(np.max(np.abs(np.array([-dpz, dpr])/rho-b1))))
    nodes, weights = np.polynomial.legendre.leggauss(96)
    for center in (-.5, .5):
        for radius in (.05, .1, .2):
            integrands = []
            for cosine in nodes:
                sine = np.sqrt(1-cosine*cosine)
                b0, b1, _ = fields(radius*sine, center+radius*cosine)
                normal = np.array([sine, cosine])
                integrands.append((b0@normal, b1@normal))
            flux = 2*math.pi*radius**2*(weights@np.array(integrands))
            flux_errors.extend((abs(flux[0]/(4*math.pi)-np.sign(center)), abs(flux[1]/(4*math.pi))))
    for nu in (.1, .5, 1., 1.5, 1.9):
        for s in (0., .3, .8, .999):
            reference = _magnetic_return_path(nu, s)
            baseline_errors.append(float(np.max(np.abs(np.array(curve(nu, s)[:3])-reference))))
            for epsilon in (-1e-4, 0., 1e-4):
                rho, z, _, _, _ = curve(nu, s, epsilon)
                p0, p1 = stream(rho, z)
                contour_errors.append(abs(p0+epsilon*p1-nu))
    assert max(maxwell) < 2e-12 and max(stream_errors) < 2e-10
    assert max(flux_errors) < 2e-12 and max(contour_errors) < 2e-12
    assert max(baseline_errors) < 2e-9

    rows, independent_errors = [], []
    for radius in (99., 999., 9999.):
        xi = radius/distance
        for nu in (.1, .5, 1., 1.5, 1.9, 2.):
            if nu == 2:
                rows.append(dict(nu=nu, outer_isotropic_radius_over_re=radius,
                    field_path_length_change_over_m=0., integral_contribution_over_m=0.,
                    moving_endpoint_contribution_over_m=0., midplane_rho_change_over_m=0.))
                continue
            root_fn = lambda t: t*(2-xi*t)*(1+xi*(2-t))-2*nu*(1+xi*(1-t))
            t = brentq(root_fn, 0., 2., xtol=1e-14)
            sc = 1-xi*t
            # Implicit derivative of the cutoff cubic, not a frozen endpoint.
            ft = ((2-2*xi*t)*(1+xi*(2-t))-xi*t*(2-xi*t)+2*nu*xi)
            s1 = -2*xi*xi*t*(2-xi*t)*(2-t)/ft
            variation = lambda s: float(np.imag(curve(nu, s, 1j*step)[2])/step)
            integral = quad(variation, 0., sc, epsabs=2e-10, epsrel=2e-11, limit=150)[0]
            refined = quad(variation, 0., sc, epsabs=2e-11, epsrel=2e-12, limit=200)[0]
            endpoint = float(curve(nu, sc)[2]*s1)
            assert abs(integral-refined) < 2e-8
            value = 2*(refined+endpoint)
            midpoint = float(np.imag(curve(nu, 0., 1j*step)[0])/step)
            assert value < 0 and midpoint < 0
            rows.append(dict(nu=nu, outer_isotropic_radius_over_re=radius,
                field_path_length_change_over_m=value,
                integral_contribution_over_m=2*refined,
                moving_endpoint_contribution_over_m=2*endpoint,
                midplane_rho_change_over_m=midpoint,
                quadrature_refinement_difference_over_m=2*abs(integral-refined)))
            if radius == 999. and nu in (.5, 1., 1.9):
                # Independent high-precision coordinate derivatives; no complex step.
                with mp.workdps(40):
                    nn, cutoff = mp.mpf(str(nu)), mp.mpf(str(sc))
                    def coordinate(s, e, component):
                        h = 1-s*s
                        u = ((1-2*e)*h+mp.sqrt((1-2*e)**2*h*h+nn*nn*s*s+4*e*nn*h))/nn
                        return mp.sqrt((u*u-1)*h)/2 if component == 0 else u*s/2
                    def speed(s, e):
                        return mp.sqrt(sum(mp.diff(lambda ss: coordinate(ss,e,k),s)**2 for k in (0,1)))
                    independent = mp.quad(lambda s: mp.diff(lambda e: speed(s,e),0), [0,.5,.9,cutoff])
                    err = abs(float(2*(independent+mp.mpf(str(endpoint))))-value)
                    independent_errors.append(err)
                    assert err < 2e-8

    # This is an annulus pullback of an isolated local metric, not a two-body solution.
    r, m, beta, w = sp.symbols('r m beta w', positive=True)
    n = (r+m-beta/(4*r))**2/(r*r+beta/4)
    transformed = (1+w)/(1-w)*n.subs(r,(1+w)*r)
    response = 2*n+r*sp.diff(n,r)
    assert sp.simplify(sp.diff(transformed,w).subs(w,0)-response) == 0
    assert sp.simplify(response.subs(beta,0)-(2+2*m/r)) == 0
    with mp.workdps(40):
        dd = mp.mpf(str(distance)); uc = 1/dd
        lo, hi = mp.mpf(99), mp.mpf(9999)
        leading = 4*uc*(hi-lo)
        mixed = 4*uc*mp.log(hi/lo)
        annulus = dict(local_mass_over_re=1., constant_companion_potential=float(uc),
            global_isotropic_endpoints_over_re=[float(lo),float(hi)],
            universal_two_mouth_length_change_over_re=float(leading),
            previously_missing_companion_span_over_re=optical['outer_cutoff_comparison']['missing_mouth_companion_potential_span_over_re'],
            universal_cancellation_residual_over_re=abs(float(leading)+optical['outer_cutoff_comparison']['missing_mouth_companion_potential_span_over_re']),
            pullback_mixed_term_over_re=float(mixed))
        assert annulus['universal_cancellation_residual_over_re'] < 1e-16
    return dict(model='First-order Maxwell response on a prescribed static weak metric with fixed opposite point-monopole fluxes; separate constant-companion annulus coordinate diagnostic',
        parameters=dict(separation_over_re=distance, expansion_parameter_m_over_d=1/distance,
            local_mass_over_re=1., weak_exterior_cutoffs_over_re=[99.,999.,9999.]),
        metric_and_maxwell='ds2=-(1-2U)dt2+(1+2U)dx2; U=m/r1+m/r2. Coordinate magnetic flux B obeys divB=0,curl[(1-2U)B]=0 to first order, not the flat-space curlB=0 equation',
        scalar_solution='A=1/r1,B=1/r2,psi0=g_m(A-B),psi1=-U*psi0-(2m/d)*psi0. deltaB=2m*g_m[(A-1/d)gradB-(B-1/d)gradA]. Here g_m is magnetic charge amplitude, not the gauge coupling. Fixed charges, no extra point multipoles and no imposed field at infinity select the harmonic term',
        stream_function='Psi0/g_m=nu0=2u(1-s2)/(u2-s2); deltaPsi/g_m=2(m/d)(1/u-1)*nu0. The quadratic contour solves the truncated field exactly, not all orders of gravity',
        path_variation='At fixed flux nu and fixed global isotropic matching radius: deltaL/m=2[integral v_epsilon ds+v(sc)*sc_epsilon]. The moving endpoint on each cutoff sphere is included. Only the first derivative at m/d=0 is a physical perturbative coefficient',
        field_path_cases=rows,
        companion_annulus=annulus,
        companion_formula='t_local=(1-Uc)t_global,r_local=(1+Uc)r_global,with fixed LOCAL RN m,Q. t_local is RN reference time, not local observer proper time. n_global=(1+Uc)/(1-Uc)*n_RN((1+Uc)r). First variation2n+r*n_prime; extremal n=(1+m/r)^2 gives2+2m/r. The universal+4Uc two-mouth boundary derivative cancels the old-4m/d term. Its mixed coefficient is only that of this pullback',
        motion_scope='Static background only. In the MMP equal-mouth circular orbit, Eq6.50 gives each mouth v^2=re/d, so the kinematic clock term-v^2/2 is the same order as-Uc. Slow angular motion does not justify omitting it from full orbital matching',
        verification=dict(maximum_linear_maxwell_scaled_residual=float(max(maxwell)),
            maximum_stream_gradient_field_difference=float(max(stream_errors)),
            maximum_monopole_flux_error_in_4pi_units=float(max(flux_errors)),
            maximum_truncated_contour_flux_error=float(max(contour_errors)),
            maximum_zero_gravity_path_difference=float(max(baseline_errors)),
            maximum_independent_mp_path_variation_difference_over_m=max(independent_errors)),
        scope='Point-source harmonic boundary choice and weak exterior only. Actual finite mouths may supply additional harmonic multipoles. No throat/source/angular Einstein response, orbital motion, second-order metric, preparation or time offset solved. Path corrections and annulus matching are not substituted into a new Casimir source or stability calculation. Mixed-order pullback terms are not a complete two-body prediction; no all-order opposite-charge metric is assumed',
        sources=['https://arxiv.org/pdf/1807.04726 (leading opposite-charge field paths and matched regions; new perturbation is our calculation)',
            'https://arxiv.org/pdf/1409.3743 (weak metric and local optical factor)',
            'https://arxiv.org/pdf/1903.00223 (Gibbons-Werner Eqs34-42,52: coordinate magnetic flux density and static Maxwell constitutive factor; no physical material-medium claim)'])


def magnetic_boundary_response_checks(static):
    """Static Maxwell Dirichlet response on fixed throat metrics; not coupled polarizability."""
    x, mu = sp.symbols('x mu', real=True)
    even, odd = sp.Integer(1), x
    bases, symbolic = {}, []
    for degree in (1, 2, 3):
        previous_even, previous_odd = even, odd
        even = sp.simplify(sp.diff(previous_odd,x)+degree*sp.tan(x)*previous_odd)
        odd = sp.simplify((sp.diff(previous_even,x)+degree*sp.tan(x)*previous_even)/degree**2)
        potential = degree*(degree+1)/sp.cos(x)**2
        for mode in (even,odd):
            residual = sp.trigsimp(sp.diff(mode,x,2)-potential*mode)
            assert sp.simplify(sp.expand_trig(residual).rewrite(sp.sin)) == 0
        wronskian_error = sp.trigsimp(even*sp.diff(odd,x)-sp.diff(even,x)*odd-1)
        assert sp.simplify(sp.expand_trig(wronskian_error).rewrite(sp.sin)) == 0
        angular = (1-mu*mu)*sp.diff(sp.legendre(degree,mu),mu)
        assert sp.simplify((1-mu*mu)*sp.diff(angular,mu,2)+degree*(degree+1)*angular) == 0
        bases[degree] = sp.lambdify(x,(even,sp.diff(even,x),odd,sp.diff(odd,x)),'numpy')
        symbolic.append(dict(degree=degree, even=str(even), odd=str(odd)))

    def response(boundary):
        ee, ep, oo, op = boundary
        off = -1/(2*ee*oo)  # Exact Wronskian=1; avoid subtracting near-equal slopes.
        diagonal = (ep/ee+op/oo)/2
        return np.array([[diagonal,off],[off,diagonal]])

    def integrate_potential(potential, endpoint, refinement=1):
        def rhs(xx,state):
            ee,ep,oo,op = state
            v = potential(xx)
            return ep,v*ee,op,v*oo
        sol = solve_ivp(rhs,(0,endpoint),(1.,0.,0.,1.),method='DOP853',
            rtol=2e-12/refinement**2,atol=2e-13/refinement**2,
            max_step=.015/refinement,dense_output=True)
        assert sol.success
        return sol

    rows, mode_errors, energy_errors = [], [], []
    for degree in (1,2,3):
        exact = bases[degree]
        for rho in (4.,10.,100.,1000.):
            endpoint = math.atan(rho)
            boundary = np.asarray(exact(endpoint),dtype=float)
            dtn = response(boundary)
            numerical = integrate_potential(lambda xx:degree*(degree+1)/math.cos(xx)**2,endpoint)
            err = float(np.max(np.abs(numerical.y[:,-1]/boundary-1)))
            mode_errors.append(err)
            assert err<2e-9 and np.linalg.eigvalsh(dtn).min()>0
            assert dtn[0,1]<0
            boundary_energy = []
            for left,right in ((1.,0.),(1.,1.),(1.,-1.)):
                ce,co = (left+right)/(2*boundary[0]),(right-left)/(2*boundary[2])
                def integrand(xx):
                    ee,ep,oo,op = exact(xx)
                    vv = degree*(degree+1)/math.cos(xx)**2
                    return 2*(ce*ce*(ep*ep+vv*ee*ee)+co*co*(op*op+vv*oo*oo))
                energy = quad(integrand,0,endpoint,epsabs=2e-10,epsrel=2e-12,limit=150)[0]
                data = np.array([left,right])
                predicted = float(data@dtn@data)
                energy_errors.append(abs(energy/predicted-1))
                assert energy_errors[-1]<2e-10
                boundary_energy.append(dict(left_value=left,right_value=right,
                    radial_energy_integral=energy,boundary_quadratic_form=predicted))
            rows.append(dict(degree=degree, throat_cutoff_rho=rho,
                diagonal=float(dtn[0,0]), off_diagonal=float(dtn[0,1]),
                cross_to_diagonal_ratio=float(-dtn[0,1]/dtn[0,0]),
                even_eigenvalue=float(dtn[0,0]+dtn[0,1]),
                odd_eigenvalue=float(dtn[0,0]-dtn[0,1]),
                positive_parity_eigenvalue_split=float(-2*dtn[0,1]),
                right_unit_left_zero_center_value=float(1/(2*boundary[0])),
                right_unit_left_zero_center_derivative=float(1/(2*boundary[2])),
                ivp_relative_mode_difference=err, energy_checks=boundary_energy))

    endpoint = static['symmetric_x_interval'][1]
    frequency2 = static['Omega_squared']
    reference = bases[1]
    ref_boundary = np.asarray(reference(endpoint),dtype=float)
    reference_dtn = response(ref_boundary)
    corrected = []
    for alpha in (static['alpha'],1e-5,5e-6):
        background,_ = _polyakov_static_solution(alpha,frequency2,endpoint,refinement=2)
        def relative_change(xx):
            yy,_,zz,_ = background.sol(xx)
            return math.expm1(2*alpha*zz-math.log1p(alpha*yy))
        def delta_potential(xx):
            return 2/math.cos(xx)**2*relative_change(xx)
        def potential(xx):
            return 2/math.cos(xx)**2*(1+relative_change(xx))
        mode = integrate_potential(potential,endpoint,refinement=2)
        boundary = mode.y[:,-1]
        dtn = response(boundary)
        direct = dtn-reference_dtn
        # Exact two-operator Green identity, plus its first Born approximation.
        eigen_changes, born_changes = [], []
        for component in (0,2):
            def exact_integrand(xx):
                return (reference(xx)[component]/ref_boundary[component]*
                    mode.sol(xx)[component]/boundary[component]*delta_potential(xx))
            def born_integrand(xx):
                return (reference(xx)[component]/ref_boundary[component])**2*delta_potential(xx)
            eigen_changes.append(quad(exact_integrand,0,endpoint,epsabs=alpha*1e-10,epsrel=2e-12)[0])
            born_changes.append(quad(born_integrand,0,endpoint,epsabs=alpha*1e-10,epsrel=2e-12)[0])
        def matrix(changes):
            diagonal,off = (changes[0]+changes[1])/2,(changes[0]-changes[1])/2
            return np.array([[diagonal,off],[off,diagonal]])
        exact_change,born_change = matrix(eigen_changes),matrix(born_changes)
        identity_error = float(np.max(np.abs(direct-exact_change)))
        assert identity_error<2e-10 and np.linalg.eigvalsh(dtn).min()>0
        corrected.append(dict(alpha=alpha, diagonal=float(dtn[0,0]),off_diagonal=float(dtn[0,1]),
            exact_green_identity_change=exact_change.tolist(),
            exact_change_over_alpha=(exact_change/alpha).tolist(),
            first_born_change=born_change.tolist(),
            first_born_maximum_remainder=float(np.max(np.abs(exact_change-born_change))),
            direct_difference_vs_green_identity_error=identity_error,
            minimum_probe_potential=min(potential(xx) for xx in np.linspace(0,endpoint,101))))
    ratio = corrected[1]['first_born_maximum_remainder']/corrected[2]['first_born_maximum_remainder']
    assert 3.9<ratio<4.1
    radius,mass,charge = sp.symbols('R M Q',positive=True)
    rn_f = 1-2*mass/radius+charge**2/radius**2
    growing = radius**2-charge**2
    assert sp.simplify(sp.diff(rn_f*sp.diff(growing,radius),radius)-2*growing/radius**2) == 0
    extremal_f = rn_f.subs(mass,charge)
    decaying = (2*radius-charge)/(2*(radius-charge)**2)
    assert sp.simplify(sp.diff(extremal_f*sp.diff(decaying,radius),radius)-2*decaying/radius**2) == 0
    assert sp.simplify(extremal_f*(growing*sp.diff(decaying,radius)-sp.diff(growing,radius)*decaying)) == -3
    return dict(model='Static source-free Maxwell angular boundary operator on prescribed spherical throat backgrounds; a probe calculation, not the coupled Einstein-Maxwell mouth response',
        equation='For ds2=e2sigma(-dt2+dx2)+r(x)^2 dOmega2 and deltaAphi=a_l(x)(1-mu2)P_l_prime(mu), a_l_second-l(l+1)e2sigma/r2*a_l=0. This perturbation changes angular magnetic flux but not the monopole charge',
        boundary_convention='At x=+-xc use common time and coordinate x; outward canonical derivatives(-a_prime(-xc),a_prime(xc))=D*(a_left,a_right). These derivatives are not the locally measured normal field without the metric factors',
        leading_geometry='e2sigma/r2=sec2(x) on global AdS2 times S2; rho=tan(x). Each finite cutoff is a diagnostic boundary, not an already matched RN mouth',
        normalized_bases=symbolic, leading_cases=rows,
        numerical_scope='At large degree or cutoff, parity eigenvalues can round to the same float. Their positive split and the off-diagonal response are evaluated from the exact unit Wronskian, not by subtracting the rounded eigenvalues',
        uniqueness='For zero Dirichlet data, integration by parts gives integral(a_prime2+V*a2)=0. Positive V implies a=0. This proves uniqueness of this finite static probe boundary problem, not dynamical stability of the coupled spacetime',
        energy='For l1 and action-F2/(4g_EM2), static probe energy with respect to the chosen time is4pi/(3g_EM2)*integral(a_prime2+2e2sigma*a2/r2)dx=4pi/(3g_EM2)*boundary_values_transpose*D*boundary_values',
        corrected_background=dict(reference_alpha=static['alpha'],Omega_squared=frequency2,
            endpoint=endpoint,cases=corrected,first_born_remainder_halving_ratio=ratio,
            identity='Exact(D2-D1)_ij=integral f_i^(1)*f_j^(2)*(V2-V1)dx for the same unit Dirichlet data. Replacing the second profile by the first is a first variation with quadratic remainder, not a finite exact identity. Common coordinate cutoff and time normalization held fixed'),
        missing_linear_gravity_source='For background Ftheta_phi=P*sin(theta) and deltaAphi=a*sin2(theta), orthonormal deltaB_radial=2a*cos(theta)/r2 and deltaB_theta=-exp(-sigma)*a_prime*sin(theta)/r. Hence delta_rho=2P*a*cos(theta)/(g_EM2*r4), deltaT_radial_radial=-delta_rho, deltaT_angular_angular=delta_rho, deltaT_radial_theta=P*exp(-sigma)*a_prime*sin(theta)/(g_EM2*r3). The sphere average of delta_rho vanishes but its l1 Einstein source does not',
        rn_matching='RN probe equation(F*a_R)_R-l(l+1)a/R2=0; for l1 the growing branch isR2-Q2. At extremality the decaying branch is(2R-Q)/(2*(R-Q)^2), with F times their Wronskian=-3. With globalt=ell*throat_tau and localRNtime=C*globalt, no sheetcurrent requires outward(+/-a_chi/ell)=C*F*a_R and continuousa in a consistent angular convention, after the induced metric is matched. Blackhole horizon regularity cannot discard a branch at a traversable finite boundary',
        verification=dict(maximum_ivp_relative_mode_difference=max(mode_errors),
            maximum_direct_energy_relative_difference=max(energy_errors),
            angular_equations_and_ladder_odes_symbolically_verified=True,
            wronskians_equal_one_symbolically_verified=True,rn_growing_mode_verified=True,
            extremal_rn_decaying_mode_and_wronskian_verified=True),
        scope='No finite-mouth multipole amplitudes, coupled metric perturbation or fermion current response selected. Electromagnetic-gauge-invariant polar-cap flux on a prescribed background is not a complete spacetime-diffeomorphism invariant. No horizon-Meissner boundary is imposed on the traversable throat. Small cross-boundary static coupling does not prove material transport, chronology change or spacetime stability; prior Casimir state/equilibrium remain unchanged',
        sources=['https://arxiv.org/pdf/1807.04726 (Maxwell action2.1, global throat5.18-5.19; boundary operator is our derivation)',
            'https://arxiv.org/pdf/1603.07894 (RN magnetic test field and blackhole boundary selection; not a traversable-mouth condition)',
            'https://arxiv.org/pdf/gr-qc/0012006 (coupled electromagnetic/gravitational perturbations about a charged background; the present calculation is only a fixed-metric probe)'])


def _linear_diagonal_einstein_maxwell(coordinates,metric,perturbation,f0,f1,volume,clean):
    """Direct 4D linear equations for a diagonal metric in the R-F^2/4 convention."""
    connection,delta_connection = {},{}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                def diagonal_derivative(values):
                    return ((sp.diff(values[a],coordinates[b]) if a==c else 0)+
                            (sp.diff(values[a],coordinates[c]) if a==b else 0)-
                            (sp.diff(values[b],coordinates[a]) if b==c else 0))/(2*metric[a])
                connection[a,b,c] = clean(diagonal_derivative(metric))
                delta_connection[a,b,c] = clean(diagonal_derivative(perturbation)-
                    perturbation[a]/metric[a]*connection[a,b,c])
    f_squared = sum(f0[a,b]**2/(metric[a]*metric[b]) for a in range(4) for b in range(4))
    delta_f_squared = clean(sum((2*f0[a,b]*f1[a,b]-f0[a,b]**2*(
        perturbation[a]/metric[a]+perturbation[b]/metric[b]))/(metric[a]*metric[b])
        for a in range(4) for b in range(4)))
    einstein = {}
    for b in range(4):
        for c in range(b,4):
            delta_ricci = sum(sp.diff(delta_connection[a,b,c],coordinates[a])-
                sp.diff(delta_connection[a,b,a],coordinates[c])+sum(
                    delta_connection[a,a,d]*connection[d,b,c]+
                    connection[a,a,d]*delta_connection[d,b,c]-
                    delta_connection[a,c,d]*connection[d,b,a]-
                    connection[a,c,d]*delta_connection[d,b,a] for d in range(4))
                for a in range(4))
            delta_stress = sum((f1[b,a]*f0[c,a]+f0[b,a]*f1[c,a]-
                f0[b,a]*f0[c,a]*perturbation[a]/metric[a])/metric[a] for a in range(4))
            if b==c:
                delta_stress -= (perturbation[b]*f_squared+metric[b]*delta_f_squared)/4
            einstein[b,c] = clean(delta_ricci-delta_stress/2)
    volume_fraction = sum(perturbation[i]/metric[i] for i in range(4))/2
    maxwell = []
    for b in range(4):
        divergence = sum(sp.diff(volume/(metric[a]*metric[b])*(f1[a,b]+f0[a,b]*(
            volume_fraction-perturbation[a]/metric[a]-perturbation[b]/metric[b])),coordinates[a])
            for a in range(4))
        maxwell.append(clean(divergence))
    return einstein,maxwell,delta_f_squared


def coupled_dipole_checks(probe):
    """Leading coupled Einstein-Maxwell dipole, with explicitly fixed diagnostic cuts."""
    from scipy.optimize import brentq
    t,x,theta,phi = sp.symbols('t x theta phi', real=True)
    coordinates = (t,x,theta,phi)
    psi = sp.Function('Psi')(t,x)
    clean = lambda value: sp.simplify(sp.expand_trig(value).rewrite(sp.sin))
    metric = (-1/sp.cos(x)**2,1/sp.cos(x)**2,sp.Integer(1),sp.sin(theta)**2)
    perturbation = tuple(-psi*sp.cos(theta)*metric[i] if i<2 else sp.Integer(0)
                         for i in range(4))
    # P=2 and r_e=1 give the paper's R-F^2/4 convention.
    f0,f1 = sp.zeros(4),sp.zeros(4)
    f0[2,3] = 2*sp.sin(theta); f0[3,2] = -f0[2,3]
    potential = (0,0,0,-psi*sp.sin(theta)**2/2)
    for a in range(4):
        for b in range(4):
            f1[a,b] = sp.diff(potential[b],coordinates[a])-sp.diff(potential[a],coordinates[b])
    volume = sp.sin(theta)/sp.cos(x)**2
    einstein,maxwell,delta_f_squared = _linear_diagonal_einstein_maxwell(
        coordinates,metric,perturbation,f0,f1,volume,clean)
    master = sp.cos(x)**2*(sp.diff(psi,x,2)-sp.diff(psi,t,2))-6*psi
    einstein_residuals = []
    for (b,c),value in einstein.items():
        expected = metric[b]*sp.cos(theta)*master/2 if b==c and b<2 else 0
        residual = clean(value-expected)
        assert residual == 0, (b,c,residual)
        einstein_residuals.append(str(residual))
    assert clean(delta_f_squared+8*psi*sp.cos(theta)) == 0
    maxwell_residuals = []
    for b,value in enumerate(maxwell):
        expected = -sp.sin(theta)*master/(2*sp.cos(x)**2) if b==3 else 0
        residual = clean(value-expected)
        assert residual == 0
        maxwell_residuals.append(str(residual))

    # Reuse the already verified degree-2 radial basis; physical angular degree here is 1.
    base = next(row for row in probe['normalized_bases'] if row['degree']==2)
    even = sp.sympify(base['even'],locals={'x':x})
    odd = sp.sympify(base['odd'],locals={'x':x})
    for value in (even,odd):
        assert clean(sp.diff(value,x,2)-6*value/sp.cos(x)**2) == 0
    assert clean(even*sp.diff(odd,x)-sp.diff(even,x)*odd-1) == 0
    modes = sp.lambdify(x,(even,sp.diff(even,x),odd,sp.diff(odd,x)),'numpy')
    rows = []
    for previous in probe['leading_cases']:
        if previous['degree'] != 2:
            continue
        rho = previous['throat_cutoff_rho']; xc = math.atan(rho)
        ee,ep,oo,op = modes(xc)
        diagonal,cross = (ep/ee+op/oo)/2,-1/(2*ee*oo)
        assert abs(diagonal-previous['diagonal'])<2e-11
        assert abs(cross-previous['off_diagonal'])<1e-15
        dtn = np.array([[diagonal,cross],[cross,diagonal]])
        acceleration = (math.sin(xc)*np.eye(2)-math.cos(xc)*dtn)/2
        rows.append(dict(cutoff_rho=rho,master_dtn=dtn.tolist(),
            fixed_cut_lapse_log_response=(-np.eye(2)/2).tolist(),
            fixed_cut_re_times_acceleration_response=acceleration.tolist(),
            unit_right_zero_left_center_psi=float(1/(2*ee)),
            boundary_energy_quadratic_form_positive=bool(np.linalg.eigvalsh(dtn).min()>0)))

    spectrum = []
    for rho in (4.,10.):
        xc = math.atan(rho)
        for parity,lower in (('even',3.),('odd',4.)):
            initial = (1.,0.) if parity=='even' else (0.,1.)
            def integrate(frequency,refinement=1):
                sol = solve_ivp(lambda xx,state:(state[1],
                    (6/math.cos(xx)**2-frequency**2)*state[0]),(0,xc),initial,
                    method='DOP853',dense_output=True,rtol=2e-12/refinement**2,
                    atol=2e-13/refinement**2,max_step=.01/refinement)
                assert sol.success
                return sol
            frequencies = [brentq(lambda w:integrate(w,level).y[0,-1],lower,lower+.9,
                                  xtol=1e-13) for level in (1,2)]
            frequency = frequencies[-1]; sol = integrate(frequency,2)
            norm = 2*quad(lambda xx:float(sol.sol(xx)[0])**2,0,xc,epsabs=1e-12,epsrel=2e-12)[0]
            stiffness = 2*quad(lambda xx:float(sol.sol(xx)[1])**2+
                6*float(sol.sol(xx)[0])**2/math.cos(xx)**2,0,xc,epsabs=1e-12,epsrel=2e-12)[0]
            error = abs(stiffness/norm-frequency**2)
            assert abs(frequencies[1]-frequencies[0])<2e-10 and error<2e-10
            assert abs(sol.y[0,-1])<1e-10
            half_points = np.linspace(0,xc,151)
            half_profile = sol.sol(half_points)[0]
            assert np.min(half_profile[1:-1])>0
            points = np.r_[-half_points[:0:-1],half_points]
            profile = np.r_[half_profile[:0:-1]*(1 if parity=='even' else -1),half_profile]
            spectrum.append(dict(cutoff_rho=rho,parity=parity,frequency_in_tau_units=frequency,
                frequency_refinement_difference=abs(frequencies[1]-frequencies[0]),
                dirichlet_boundary_residual=float(abs(sol.y[0,-1])),
                rayleigh_identity_error=error,mode_norm=norm,
                positive_wave_energy_at_unit_displacement_phase=stiffness/2,
                plot=dict(chi=points.tolist(),mode_profile=profile.tolist()),
                evolution='Psi(tau,chi)=amplitude*mode_profile(chi)*cos(omega*tau), '
                    'with homogeneous Dirichlet values at the fixed diagnostic cuts'))

    rho = 4.; xc = math.atan(rho)
    points = np.linspace(-xc,xc,301)
    ee,_,oo,_ = modes(xc)
    profile = (modes(points)[0]/ee+modes(points)[2]/oo)/2
    assert np.all(np.diff(profile)>0) and np.max(np.abs(profile[[0,-1]]-[0,1]))<1e-14
    amplitude_symbol = sp.symbols('epsilon',real=True)
    multiplier = 1-amplitude_symbol*psi*sp.cos(theta)
    finite_metric = tuple(metric[i]*multiplier if i<2 else metric[i] for i in range(4))
    finite_divergence = sum(sp.diff(clean(volume*multiplier*(f0[a,3]+amplitude_symbol*f1[a,3])/
        (finite_metric[a]*finite_metric[3])),coordinates[a]) for a in range(4))
    finite_expected = (-amplitude_symbol*sp.sin(theta)*master/(2*sp.cos(x)**2)-
        2*amplitude_symbol**2*psi**2/sp.cos(x)**2*sp.sin(theta)*sp.cos(theta))
    assert clean(finite_divergence-finite_expected) == 0
    # This finite-amplitude continuation is deliberately tested as an ansatz, not an exact solution.
    # Its exact Maxwell density divergence on a mass-6 solution is quadratic in amplitude.
    nonlinear = []
    for amplitude in (1e-3,5e-4):
        residual = -2*amplitude**2*profile**2/np.cos(points)**2*.5
        nonlinear.append(dict(amplitude=amplitude,
            maximum_abs_exact_maxwell_density_residual_at_theta_pi_over_4=float(np.max(np.abs(residual))),
            metric_multiplier_lower_bound=1-amplitude))
    ratio = nonlinear[0]['maximum_abs_exact_maxwell_density_residual_at_theta_pi_over_4']/nonlinear[1][
        'maximum_abs_exact_maxwell_density_residual_at_theta_pi_over_4']
    assert abs(ratio-4)<1e-13
    return dict(
        model='Physical angular dipole scalar of coupled classical Einstein-Maxwell theory on leading global AdS2 times S2; finite cuts are diagnostic, not joined mouths',
        conventions='r_e=1 in computed geometry, P=2 in R-F^2/4 action. Y=cos(theta), pi=g_background^(mn)*delta_g_mn/Y; Psi=pi-4a/P. Restoring units gives mass_squared*r_e_squared=6 and physical time=ell*tau',
        reconstruction='Sphere-trace pi=0 gauge: delta_g_AB=-Psi*g_AB*cos(theta), delta_g_mn=delta_g_Am=0, delta_A_phi=-P*Psi*sin(theta)^2/4. The independent vector sector is set to zero',
        master_equation='(-partial_tau^2+partial_chi^2-6*sec(chi)^2)Psi=0',
        direct_four_dimensional_check=dict(einstein_components=einstein_residuals,
            maxwell_components=maxwell_residuals,
            method='Direct variation of the diagonal 4D metric connection, Ricci tensor, F tensor, inverse metric, magnetic stress and Maxwell density divergence for arbitrary Psi(tau,chi)',
            off_shell_trace_reversed_residual='Only AB block is nonzero: delta[R_AB-.5*(F_AC*F_B^C-.25*g_AB*F^2)]=.5*g_AB*cos(theta)*(Box_AdS2-6)Psi. Hence delta_R4=cos(theta)*(Box_AdS2-6)Psi; Einstein-tensor residual instead has only sphere block',
            magnetic_invariant='delta(F^2)/F0^2=-Psi*cos(theta). Constant background F^2 makes this scalar perturbation invariant under infinitesimal coordinate changes'),
        physical_boundary_observables=dict(
            cap_invariant='delta[Flux_cap/(2*pi*P)-Area_cap/(2*pi*r_e^2)]=-Psi*sin(theta_cap)^2/4, including displaced cap edges since the background difference vanishes',
            fixed_cut_lapse='delta(log N)=-Psi*cos(theta)/2 in pi=0 gauge and the common time normalization',
            fixed_cut_acceleration='r_e*delta(kappa)_+-=cos(chi_c)/2*[rho_c*Psi_+--d_+-]*cos(theta), d_+-=+-Psi_prime(+-chi_c). kappa=n dot grad logN equals the mixed tt extrinsic curvature for these static cuts',
            bending='For displaced boundary chi=c+eta1*cos(theta), delta(log N_boundary)=[-Psi(c)/2+eta1*tan(c)]*cos(theta). Fixed-cut lapse is not a boundary-embedding invariant',
            matching_gap='A joined exterior must match induced metric and extrinsic curvature (or supply a specified surface stress) along with flux/current. The old probe-Maxwell RN relation omits metric response at the same perturbation order'),
        static_master_boundary_cases=rows,
        static_profile=dict(cutoff_rho=rho,chi=points.tolist(),unit_right_zero_left_psi=profile.tolist(),
            illustrative_amplitude=1e-3,center_psi_per_unit_right_boundary=float(profile[len(points)//2])),
        clamped_scalar_modes=spectrum,
        conditional_linear_stability='For fixed homogeneous Dirichlet cuts, E=.5*integral(Psi_tau^2+Psi_chi^2+6*sec(chi)^2*Psi^2)dchi is positive and dE/dtau=[Psi_tau*Psi_chi]_left^right=0. It excludes growing finite-energy solutions in this linear scalar boundary problem. Its physical action prefactor is not derived here. It is not coupled exterior, vector-sector, quantum-source or nonlinear stability',
        finite_amplitude_maxwell_check=dict(cases=nonlinear,halved_amplitude_residual_ratio=ratio,
            exact_divergence='For g_AB=(1-epsilon*Psi*cos(theta))*g0_AB and A_phi=A0_phi-epsilon*P*Psi*sin(theta)^2/4 with P=2, partial_M(sqrt(-g)*F^(Mphi))=-epsilon*sin(theta)/(2*cos(chi)^2)*(Box-6)Psi-2*epsilon^2*Psi^2*sec(chi)^2*sin(theta)*cos(theta)',
            meaning='The linear solution does not solve finite-amplitude source-free Maxwell exactly. Higher-order corrections to fields and/or geometry are required; only this Maxwell residual is evaluated, not every nonlinear Einstein component'),
        scope='The angular geometry is now included in this leading classical scalar sector. A quantum-supported radius-varying MMP throat requires rederiving coupled equations with angular quantum stress/current response; replacing sec^2 by the earlier probe corrected potential is not justified. No exterior junction, formation, material generation, traversable payload route or past calendar map is established',
        sources=['https://arxiv.org/pdf/hep-th/9906056 (Einstein-Maxwell action and magnetic background2.1-2.6; scalar and metric equations3.28,3.36-3.39; reconstruction and finite-boundary calculations here are our derivation)',
            'https://arxiv.org/pdf/1212.2557 (section6, table1: invariant RN physical dipole has near-horizon mass coefficient6)',
            'https://arxiv.org/pdf/1807.04726 (leading throat and quantum support; not a published solution of the present angular boundary problem)'])


def coupled_exterior_checks(coupled):
    """Physical polar dipole response of a truncated extremal RN exterior."""
    r,m,q = sp.symbols('r M Q',positive=True)
    f_general = 1-2*m/r+q*q/r**2
    varpi = 3*m-2*q*q/r
    w = 3*m*(r*r-2*m*r+q*q)/(r*varpi**2)+m/varpi
    u_general = sp.factor((6*m*w+varpi-m)/r**3)
    expected = 2*(9*m*m*r**4-18*m*m*q*q*r*r+16*m*q**4*r-4*q**6)/(r**4*(3*m*r-2*q*q)**2)
    assert sp.factor(u_general-expected) == 0
    f = f_general.subs({m:1,q:1})
    u = u_general.subs({m:1,q:1})
    growing = 3*r*(r-1)**2/(3*r-2)
    decaying = r*(30*r*r-45*r+17)/(10*(r-1)**3*(3*r-2))
    for mode in (growing,decaying):
        assert sp.factor(sp.diff(f*sp.diff(mode,r),r)-u*mode) == 0
    assert sp.factor(f*(growing*sp.diff(decaying,r)-decaying*sp.diff(growing,r))+3) == 0
    assert sp.limit(growing/r**2,r,sp.oo) == 1
    assert sp.limit(r*decaying,r,sp.oo) == 1
    impedance = sp.factor(-f*sp.diff(decaying,r)/decaying)
    assert sp.factor(f*sp.diff(growing,r)+impedance*growing-3/decaying) == 0
    delta = sp.symbols('delta',positive=True)
    positive_u = sp.Poly(sp.fraction(sp.factor(u.subs(r,1+delta)))[0],delta)
    positive_k = sp.Poly(sp.fraction(sp.factor(impedance.subs(r,1+delta)))[0],delta)
    for polynomial in (positive_u,positive_k):
        assert all(v>=0 for v in polynomial.all_coeffs()) and any(v>0 for v in polynomial.all_coeffs())
    assert sp.limit(u,r,1,dir='+') == 6
    assert sp.limit(impedance/(r-1),r,1,dir='+') == 3
    assert sp.limit(r*impedance,r,sp.oo) == 1
    scale,rho = sp.symbols('epsilon rho',real=True)
    throat_subs = {m:q*(1-scale**2/2),r:q*(1+scale*rho)}
    assert sp.simplify(sp.limit(f_general.subs(throat_subs)/scale**2,scale,0)-(1+rho*rho)) == 0
    assert sp.simplify(sp.limit(q*q*u_general.subs(throat_subs),scale,0)-6) == 0
    assert 'mass_squared*r_e_squared=6' in coupled['conventions']
    fn = sp.lambdify(r,f,'numpy')
    un = sp.lambdify(r,u,'numpy')
    dn = sp.lambdify(r,(decaying,sp.diff(decaying,r)),'numpy')
    kn = sp.lambdify(r,impedance,'numpy')
    cuts = (1.01,1.1,1.5,2.,4.,10.)
    rows = []
    for cut in cuts:
        dc = float(dn(cut)[0])
        # Compactify infinity; integrate normalized fields, not huge near-horizon amplitudes.
        def energy_integrand(z):
            if z == 1:
                return 0.
            radius = cut+z/(1-z)
            dd,dp = dn(radius)
            return (fn(radius)*(dp/dc)**2+un(radius)*(dd/dc)**2)/(1-z)**2
        energy,error = quad(energy_integrand,0,1,epsabs=3e-12,epsrel=3e-12,limit=150)
        k = float(kn(cut))
        assert k>0 and abs(energy-k)<2e-10
        rows.append(dict(cut_radius_over_M=cut,decaying_master_impedance_M=k,
            unit_normal_log_derivative_M=k/math.sqrt(float(fn(cut))),
            normalized_exterior_static_energy=energy,energy_identity_error=abs(energy-k),
            quadrature_error_estimate=error))

    # Positive real Laplace s selects decay at infinity. The response depends on s;
    # its zero-frequency value is not a time-domain transparent boundary condition.
    laplace = []
    for rate in (.02,.1,.5):
        traces = []
        for outer,refinement in ((600.,1),(1200.,2)):
            def rhs(y,state):
                offset = math.exp(y); radius = 1+offset
                return [offset*((state[0]**2-rate**2)/fn(radius)-un(radius))]
            solution = solve_ivp(rhs,(math.log(outer-1),math.log(cuts[0]-1)),[rate],
                method='DOP853',rtol=2e-12/refinement**2,atol=2e-13/refinement**2,
                max_step=.02/refinement,dense_output=True)
            assert solution.success
            traces.append(solution)
        values = []
        for row in rows:
            cut = row['cut_radius_over_M']; y = math.log(cut-1)
            kvals = [float(solution.sol(y)[0]) for solution in traces]
            difference = abs(kvals[1]-kvals[0])
            assert difference<2e-9 and kvals[-1]>row['decaying_master_impedance_M']
            values.append(dict(cut_radius_over_M=cut,impedance_M=kvals[-1],
                combined_outer_boundary_and_step_refinement_difference=difference,
                excess_over_static=kvals[-1]-row['decaying_master_impedance_M']))
        laplace.append(dict(laplace_rate_times_M=rate,cases=values,
            outer_radii_over_M=[600.,1200.],outer_boundary='K(Rmax)=s; asymptotic truncation checked by doubling Rmax with finer integration'))

    # Reproduce the static logarithmic derivative via the ODE, independently of d(r).
    static_runs = []
    for outer in (100.,200.):
        def static_rhs(y,state):
            offset = math.exp(y); radius = 1+offset
            return [offset*(state[0]**2/fn(radius)-un(radius))]
        solution = solve_ivp(static_rhs,(math.log(outer-1),math.log(cuts[0]-1)),
            [float(kn(outer))],method='DOP853',rtol=3e-13,atol=3e-14,max_step=.01,
            dense_output=True)
        assert solution.success
        static_runs.append(solution)
    for row in rows:
        y = math.log(row['cut_radius_over_M']-1)
        error = max(abs(float(sol.sol(y)[0])-row['decaying_master_impedance_M']) for sol in static_runs)
        assert error<2e-10
        row['static_riccati_versus_exact_error'] = error
    radii = 1+np.geomspace(.01,30.,241)
    return dict(
        model='Coupled polar l1 Einstein-Maxwell Moncrief master on an extremal RN exterior r>=r_cut>M; no throat junction imposed',
        equation='(-partial_t^2+partial_rstar^2-V)H=0, drstar/dr=1/F, V=F*U. Static: (F H_prime)_prime-U H=0. F=1-2M/r+Q^2/r^2',
        general_polar_potential_U=str(expected),
        conventions='Computed M=Q=1; r means R/M, time means t_RN/M. Physical K has units1/M. H is the physical polar l1 master of electric RN, obtained from the magnetic scalar sector by electromagnetic duality; this is not the axial dipole or a test Maxwell field',
        gauge_invariant='Moncrief H=f1-(4Q/Lambda)*khat2, Lambda=6M/r-4Q^2/r^2, khat2=F*H2-F^(3/2)*partial_r(r*F^(-1/2)*Kmetric)+Kmetric-2F*h1/r. Definitions and harmonic normalization follow Lucietti et al. AppendixB Eq137; do not identify this H with the throat Psi by name',
        exact_static_basis=dict(growing=str(growing),decaying=str(decaying),weighted_wronskian=-3,
            normalization='g/r^2 ->1 and r*d ->1 at infinity; both are valid at a finite exterior cut. Infinity decay selects H=d/d(r_cut)',
            horizon_scope='The decaying mode diverges at r=M, which lies outside the truncated domain. Imposing black-hole horizon regularity would remove this mode and would change the boundary problem'),
        static_cases=rows,
        impedance='K0=-F H_prime/H at r_cut for decay at infinity. The outward normal of the exterior domain at its inner boundary points toward decreasing r. K0 is its positive canonical logarithmic derivative',
        normal_and_clock='The exterior domain outward normal at its inner boundary points toward decreasing r. The unit-normal logarithmic derivative is K0/sqrt(F), tending3/M near r=M, whereas canonical K0 tends0. If t_RN=C*t_global, the canonical impedance measured in global optical coordinates is C*K0',
        forced_static_relation='For H=A*g+B*d and p_ext=F*H_prime pointing into the exterior, p_ext+K0*H=3*A/d(r_cut). A is a growing master amplitude at infinity, not an already normalized physical magnetic field. A=0 selects decay; imposed external fields require the inhomogeneous term',
        leading_normalization='With action R-F_magnetic^2/4 and P=2Q, choose electric dual F_e=-star(F_magnetic), including the perturbed Hodge star. In the matching electric convention f1=2a. A regular near-horizon scaling of the sphere-trace-zero representative gives khat2->0, a->-Q*Psi/2 and H->-Q*Psi. This is only a leading overlap relation; finite-radius corrections, boundary motion and common clocks still require reconstruction',
        formal_global_throat_limit='epsilon=Q/ell, M=Q*(1-epsilon^2/2), r=Q*(1+epsilon*rho), t_RN=ell*tau gives F/epsilon^2->1+rho^2 and Q^2*U->6. With chi=atan(rho) this recovers the global mass6 master. This separate formal near-extremal scaling has M<Q; it is not the exactly extremal exterior used for the numeric impedance, a complete regular spacetime or a fitted MMP quantum mass/radius relation',
        static_energy='With H(r_cut)=1 and H decaying at infinity, integral_cut^infinity[F H_prime^2+U H^2]dr=K0. This is the reduced master quadratic form; its normalization to the throat physical action is not derived here',
        positivity=dict(U_numerator_in_delta=str(positive_u.as_expr()),
            K_numerator_in_delta=str(positive_k.as_expr()),
            proof='delta=r/M-1>0; numerator coefficients are nonnegative with positive terms, and all denominator factors are positive',
            near_horizon_U=6,near_horizon_K_over_delta=3,infinity_r_times_K=1),
        laplace_cases=laplace,
        laplace_equation='For exp(s*t), (F h_prime)_prime-(U+s^2/F)h=0; K=-F h_prime/h obeys dK/dr=(K^2-s^2)/F-U. Real s>0 and decay at infinity define the sampled response',
        laplace_scope='K(s) varies with s. Sampling positive real Laplace rates is not a complete causal convolution kernel or a driven throat evolution; no real-frequency transmission coefficient is inferred',
        matching_gap='Must reconstruct metric/field data in a common duality and harmonic convention and match induced metric, extrinsic curvature, flux/current and clocks at a moving finite mouth, including quantum stress response. No numerical identification H=Psi or sum of unrelated boundary operators is imposed',
        plot=dict(radius_over_M=radii.tolist(),static_impedance_M=[float(kn(rr)) for rr in radii],
            potential_U=[float(un(rr)) for rr in radii]),
        verification=dict(static_modes_and_wronskian_symbolically_verified=True,
            general_polar_potential_algebra_verified=True,
            maximum_static_energy_identity_error=max(row['energy_identity_error'] for row in rows),
            maximum_static_riccati_error=max(row['static_riccati_versus_exact_error'] for row in rows),
            maximum_laplace_refinement_difference=max(case['combined_outer_boundary_and_step_refinement_difference'] for trace in laplace for case in trace['cases'])),
        scope='Exterior response is a necessary new input for mouth matching, not a joined or created spacetime. Positivity in this reduced exterior sector does not establish full throat, quantum, nonlinear or payload stability. No past-calendar route follows',
        sources=['https://arxiv.org/pdf/gr-qc/9507012 (Burko, physical polar dipole equation87 and definitions of varpi/W; exact extremal basis and response calculations here are our derivation)',
            'https://arxiv.org/pdf/1212.2557 (Lucietti et al., Moncrief wave equation81/82, physical l1 polar master and AppendixB Eq137; near-horizon coefficient6)'])


def coupled_boundary_checks(exterior):
    """Finite-radius static RN fields and embedded boundary data; no junction imposed."""
    t,theta,phi = sp.symbols('t theta phi',real=True)
    r = sp.symbols('r',positive=True)
    h,integral = sp.Function('H')(r),sp.Function('I')(r)
    f,sqrtf,varpi = (r-1)**2/r**2,(r-1)/r,3-2/r
    u = sp.sympify(exterior['general_polar_potential_U'],locals={'r':r,'M':1,'Q':1})
    h_second = (u*h-sp.diff(f,r)*sp.diff(h,r))/f
    integral_prime = -h/(2*r*sqrtf*varpi)
    substitutions = {sp.diff(h,r,2):h_second,sp.diff(integral,r):integral_prime}
    substitutions[sp.diff(h,r,3)] = sp.diff(h_second,r).subs(substitutions)
    substitutions[sp.diff(integral,r,2)] = sp.diff(integral_prime,r).subs(substitutions)
    def clean(value):
        return sp.factor(sp.trigsimp(sp.expand(value).subs(substitutions,simultaneous=True)))
    n0 = sqrtf/r**2*(1-2/r)
    l0 = sqrtf*(3*r-4)/r**3
    n = n0*integral+f*sp.diff(h,r)/(2*varpi)+(varpi+2*f/r)*h/(2*r*varpi**2)
    l = l0*integral+h/(2*r**2)
    k = -sqrtf*varpi*integral/r**2
    a = h/2-2*sqrtf*integral/r
    metric = (-f,1/f,r*r,r*r*sp.sin(theta)**2)
    perturbation = tuple(2*value*sp.cos(theta)*g for value,g in zip((n,l,k,k),metric))
    f0,f1 = sp.zeros(4),sp.zeros(4)
    f0[2,3] = 2*sp.sin(theta); f0[3,2] = -f0[2,3]
    f1[1,3] = sp.diff(a,r)*sp.sin(theta)**2; f1[3,1] = -f1[1,3]
    f1[2,3] = 2*a*sp.sin(theta)*sp.cos(theta); f1[3,2] = -f1[2,3]
    einstein,maxwell,_ = _linear_diagonal_einstein_maxwell(
        (t,r,theta,phi),metric,perturbation,f0,f1,r*r*sp.sin(theta),clean)
    assert all(value==0 for value in einstein.values()),einstein
    assert all(value==0 for value in maxwell),maxwell
    constraint = clean(sp.diff(k,r)+(1/r-sp.diff(f,r)/(2*f))*k-l/r)
    khat2 = 2*f*l-f*sqrtf*sp.diff(2*r*k/sqrtf,r)+2*k
    moncrief = clean(2*a-4/(2*varpi/r)*khat2-h)
    assert constraint==0 and moncrief==0

    # Surface r=c+zeta*Y, theta=Theta+b*D_Theta Y. The same normal is used
    # on both sides of a proposed junction; domain-outward normals differ.
    zeta,b = sp.symbols('zeta b',real=True)
    lapse_gradient = sp.diff(f,r)/(2*f)
    def boundary(nn,ll,kk,aa,displacement,relabel):
        return tuple(clean(value) for value in (
            nn+displacement*lapse_gradient,
            kk+displacement/r-relabel,
            sqrtf*(sp.diff(nn,r)-ll*lapse_gradient)+displacement*sp.diff(sqrtf*lapse_gradient,r),
            sqrtf*(sp.diff(kk,r)-ll/r)+displacement*varpi/(r**3*sqrtf),
            aa/2-kk-displacement/r,
            sqrtf*(sp.diff(aa,r)+2*displacement/(f*r*r))))
    observables = boundary(n,l,k,a,zeta,b)
    names = ('lapse_log','angular_metric_half_fraction','mixed_tt_extrinsic_curvature',
             'mixed_angular_extrinsic_curvature','cap_flux_area_contrast','normal_F_phi')
    # A constant of integration is gauge, not another magnetic or gravitational mode.
    radial_gauge = sqrtf*(f-r*sp.diff(f,r)/2)
    angular_gauge = sqrtf/r
    gauge_fields = tuple(sp.diff(value,integral) for value in (n,l,k,a))
    lie_fields = (lapse_gradient*radial_gauge,
        sp.diff(radial_gauge,r)-lapse_gradient*radial_gauge,
        radial_gauge/r-angular_gauge,-2*angular_gauge)
    assert all(clean(x-y)==0 for x,y in zip(gauge_fields,lie_fields))
    assert clean(radial_gauge/f+r*r*sp.diff(angular_gauge,r))==0
    gauge_residuals = boundary(*gauge_fields,-radial_gauge,-angular_gauge)
    assert all(value==0 for value in gauge_residuals),gauge_residuals

    # Fixed geometric cut in the I(c)=0 representative; not a claim that
    # holding its coordinate fixed in another gauge preserves the surface.
    decaying = sp.sympify(exterior['exact_static_basis']['decaying'],locals={'r':r})
    logarithmic_derivative = sp.diff(decaying,r)/decaying
    fixed = tuple(clean(value).subs({sp.diff(h,r):logarithmic_derivative,h:1,
        integral:0,zeta:0,b:0},simultaneous=True) for value in observables)
    # Independent cut identities from the first-order radial constraints.
    assert clean(fixed[2]-(fixed[0]-1/(2*r*r))/(r*sqrtf))==0
    assert clean(fixed[5]-sqrtf*(logarithmic_derivative/2+1/(r*r*varpi)))==0
    values = sp.lambdify(r,fixed,'numpy')
    rows = []
    for case in exterior['static_cases']:
        cut = case['cut_radius_over_M']
        response = [float(value) for value in values(cut)]
        assert all(math.isfinite(value) for value in response)
        assert abs(response[1])<1e-15 and abs(response[4]-.25)<1e-15
        assert abs(response[3])<1e-15
        rows.append(dict(cut_radius_over_M=cut,**dict(zip(names,response))))
    return dict(
        model='Static physical polar l1 fields of the magnetic dual RN exterior, reconstructed at finite radius; no quantum throat joined',
        conventions='Computed M=Q=1, r>1, P=2 in R-F_m^2/4; g0=diag(-F,1/F,r^2,r^2*sin(theta)^2), delta_g_ii=2*(n,l,k,k)_i*Y*g0_ii, delta_A_phi=a*sin(theta)^2, Y=cos(theta). All fields are first order',
        reconstruction=dict(F=str(f),varpi=str(varpi),I_prime=str(integral_prime),
            n=str(n),l=str(l),k=str(k),a=str(a),
            master='(F H_prime)_prime-U H=0; I_prime=-Q*H/(2*r*sqrt(F)*varpi). I has one additive gauge constant',
            source_mapping='Burko Eqs88-99 with special functions from71-75 and Hbar=-H/2. Magnetic dual with P=2Q gives electric f1=2a, H2=2l, Kmetric=2k,h1=0 in Lucietti Eq137'),
        direct_four_dimensional_check=dict(einstein_components=[str(v) for v in einstein.values()],
            maxwell_components=[str(v) for v in maxwell],radial_constraint=str(constraint),
            moncrief_reextraction_residual=str(moncrief),
            method='Direct linear variation of all ten independent Ricci-minus-half-Maxwell-stress components and four Maxwell density divergences, for arbitrary H(r),I(r) subject to the master ODE and I quadrature. No field equations discarded'),
        embedded_boundary=dict(
            embedding='r=c+zeta*cos(Theta), theta=Theta+b*D_Theta cos(Theta). Common static time normalization. Formulas use the normal toward increasing r; reverse all extrinsic and normal-F amplitudes for the opposite normal',
            coefficients=dict(zip(names,map(str,observables))),
            angular_dependence='lapse, half fractional angular metric and both mixed extrinsic curvature amplitudes multiply cos(Theta); cap contrast multiplies sin(Theta_cap)^2; normal_F_phi multiplies sin(Theta)^2',
            cap_definition='delta[Flux_cap/(2*pi*P)-Area_cap/(2*pi*c^2)] at fixed background cut c, allowing displaced cap labels. Contrast coefficient=a/P-k-zeta/c; intrinsic relabel b cancels',
            maxwell_meaning='F_(normal phi) is normal contraction of the field tensor, not a current. Its jump under a consistently oriented normal would require a specified surface current',
            shape='A polar l1 conformal angular metric perturbation is an infinitesimal sphere relabeling. It does not by itself create an intrinsic shape or radius deformation'),
        residual_gauge=dict(radial_generator=str(radial_gauge),angular_generator=str(angular_gauge),
            convention='I->I+C has delta_g=Lie_xi g0, xi=C*(radial_generator*Y partial_r+angular_generator*D^A Y partial_A). The same physical surface has zeta->zeta-C*radial_generator(c), b->b-C*angular_generator(c)',
            all_six_boundary_residuals=[str(v) for v in gauge_residuals],mixed_r_theta_metric_residual='0'),
        static_boundary_cases=rows,
        case_normalization='At each cut independently H(c)=1 in M=1 units, infinity decay selects H_prime(c)=d_prime(c)/d(c), I(c)=0,zeta=b=0. These are linear response coefficients, not finite-amplitude geometries. Extrinsic curvature entries restore as1/M. The normal_F_phi entry retains a coordinate phi index; divide additionally by r*sin(Theta) for an orthonormal angular component',
        junction_obstruction='The leading exact product throat has background mixed angular curvature zero, whereas every regular finite RN cut has sqrt(F)/r>0. A first-order l1 displacement cannot cure that zeroth-order mismatch. The radius-varying Polyakov model instead has the separate previously verified quantum normal-pressure obstruction when abruptly joined to electrovac RN',
        next_required='A compatible radius-varying quantum-supported background and its angular quantum stress/current response, followed by one common embedded-surface junction. These exterior boundary data are necessary inputs only',
        scope='Static linear classical exterior only. No complete dynamic reconstruction, joined mouth, stability of the quantum source, formation mechanism, material recipe, payload transport or past-calendar map',
        sources=['https://arxiv.org/pdf/gr-qc/9507012 (Burko finite reconstruction71-75,88-99; magnetic dual and explicit finite-cut checks here)',
            'https://arxiv.org/pdf/1212.2557 (Lucietti et al. AppendixB Eq137 master normalization)',
            'https://arxiv.org/pdf/gr-qc/0507005 (Mars Proposition2 and Section8, perturbed hypersurface metric and extrinsic curvature; specialized boundary formulas here are our derivation)'])


def zero_delay_spectrum_checks():
    """Analytic characteristic-mode check for our D=0 linear return closure."""
    import mpmath as mp
    lam, tau, chi, x = sp.symbols('lam tau chi x', real=True)
    k = 4/(3*sp.pi)
    H = lam*(lam**2+1)*sp.sinh(sp.pi*lam/2)+k*(lam**2+4)*sp.cosh(sp.pi*lam/2)
    assert sp.simplify(H.subs(lam, sp.I)) == 0
    assert sp.simplify(sp.diff(H, lam).subs(lam, sp.I)) == 0
    assert sp.simplify(sp.diff(H, lam, 2).subs(lam, sp.I)/2+sp.Rational(13, 3)) == 0
    def c_mode(C):
        return ((C.subs(x, tau+chi)+C.subs(x, tau-chi))*sp.tan(chi)/2
                +(sp.diff(C, x).subs(x, tau+chi)-sp.diff(C, x).subs(x, tau-chi))/2)
    pure = c_mode(sp.exp(sp.I*x))
    generalized = c_mode(x*sp.exp(sp.I*x))
    assert sp.simplify(sp.expand_complex(pure)) == 0
    target = sp.I*sp.exp(sp.I*tau)*(sp.sin(chi)+chi/sp.cos(chi))
    assert sp.simplify(sp.expand_complex(generalized-target)) == 0
    # S=x produces a static background-shaped field, not a growing area.
    assert sp.simplify(((tau+chi)-(tau-chi))*sp.tan(chi)/2+1-(1+chi*sp.tan(chi))) == 0
    with mp.workdps(50):
        errors = []
        for ll in (mp.mpc('.4', '.7'), mp.mpc('.2', '2.3'), mp.mpc('.8', '4.1')):
            zz = ll*ll
            series = 4+3*zz*mp.nsum(lambda n: 1/(zz+(2*n+1)**2), [1, mp.inf])
            left = ll*(zz+1)*mp.sinh(mp.pi*ll/2)+4/(3*mp.pi)*(zz+4)*mp.cosh(mp.pi*ll/2)
            right = 4/(3*mp.pi)*(zz+1)*mp.cosh(mp.pi*ll/2)*series
            errors.append(abs(left-right))
            assert mp.im(series)*mp.im(zz)>0
        assert max(errors)<mp.mpf('1e-45')
        def characteristic(w):
            return -w*(1-w*w)*mp.sin(mp.pi*w/2)+4/(3*mp.pi)*(4-w*w)*mp.cos(mp.pi*w/2)
        frequencies = [mp.mpf(2)]+[mp.findroot(characteristic, (2*n+mp.mpf('.01'), 2*n+mp.mpf('.2'))) for n in range(2, 6)]
        assert all(abs(characteristic(w))<mp.mpf('1e-44') for w in frequencies)
    return dict(
        characteristic='H=lambda*(lambda^2+1)*sinh(pi*lambda/2)+(4/(3*pi))*(lambda^2+4)*cosh(pi*lambda/2). The full C determinant is2*exp(-pi*lambda/2)*H',
        factorization='z=lambda^2; H=(4/(3*pi))*(z+1)*cosh(pi*lambda/2)*F(z); F(z)=4+3*z*sum_(n=1..infinity)1/[z+(2*n+1)^2]',
        proof='Im F(z)=3*Im(z)*sum_(n=1..infinity)(2*n+1)^2/abs(z+(2*n+1)^2)^2. It cannot vanish for nonreal z. F is positive for real z>=0, so its zeros are negative real. Odd-integer poles n>=1 cancel cosine zeros without extra characteristic zeros. lambda=+/-i are separately double roots',
        conclusion='No exponentially growing characteristic modes in this specific D=0 linear closure; this is not a nonlinear or finite-delay stability theorem',
        double_root_check='H(i)=H_prime(i)=0; H_second(i)/2=-13/3',
        representation='C=exp(ix) reconstructs zero field. C=x*exp(ix) reconstructs i*exp(i*tau)*(sin(chi)+chi*sec(chi)), with secular x cancelling. It is an infinitesimal AdS2 frame change of the background. S=constant is redundant; S=x gives a static background-shaped area shift',
        sum_channel='Determinant=lambda*(lambda^2+1)*(1-exp(-pi*lambda)); zero is double, +/-i and nonzero even imaginary integers are simple',
        first_nonredundant_positive_frequencies=[mp.nstr(w, 25) for w in frequencies],
        numerical_factorization_errors=[mp.nstr(e, 5) for e in errors],
        scope='Characteristic analysis plus exact field-reconstruction identities for the zero-delay benchmark. No general nonlinear initial-history stability, damping, finite-D spectral theorem, quantum-state preparation, causality violation or physical apparatus follows')


def driven_clock_checks():
    """Retarded CFT on prescribed smooth mouth maps; test, do not assume, JT balance."""
    import mpmath as mp
    from scipy.optimize import brentq
    epsilon, width, delay = .05, 1., .1
    half, length = math.pi/2, math.pi+.1
    slope = 1+epsilon
    f0 = -math.pi/(12*length**2)  # Stress divided by q throughout this function.
    b0 = math.pi*((2*math.pi/length)**2-1)/6  # b0/alpha, including finite D.

    def step(x):
        if x <= 0: return (0., 0., 0., 0.)
        if x >= 1: return (1., 0., 0., 0.)
        logit = 1/(1-x)-1/x
        small = math.exp(-abs(logit))
        g = small/(1+small) if logit < 0 else 1/(1+small)
        variance = small/(1+small)**2
        a1 = 1/(1-x)**2+1/x**2
        a2 = 2/(1-x)**3-2/x**3
        a3 = 6/(1-x)**4+6/x**4
        return (g, variance*a1, variance*((1-2*g)*a1*a1+a2),
                variance*((1-6*g+6*g*g)*a1**3+3*(1-2*g)*a1*a2+a3))

    def clock(s):
        if s <= 0: return (s, 1., 0., 0., 0.)
        if s >= width: return (slope*s-epsilon*width/2, slope, 0., 0., 0.)
        g, g1, g2, g3 = step(s/width)
        integral = width*quad(lambda x: step(x)[0], 0, s/width,
                              epsabs=2e-14, epsrel=2e-13)[0]
        return (s+epsilon*integral, 1+epsilon*g, epsilon*g1/width,
                epsilon*g2/width**2, epsilon*g3/width**3)

    def inverse(y):
        if y <= 0: return y
        if y >= width*(1+epsilon/2): return (y+epsilon*width/2)/slope
        return brentq(lambda s: clock(s)[0]-y, 0, width, xtol=5e-15)

    exit_shortcut = (math.pi-delay)/epsilon+width/2
    exit_loop = length/epsilon+width/2
    entry_shortcut, entry_loop = exit_shortcut-delay, exit_loop+delay
    fixed = entry_loop+half
    with mp.workdps(40):
        fixed_precise = ((mp.pi+mp.mpf(str(delay)))/mp.mpf(str(epsilon))
                         +mp.mpf(str(width))/2+mp.mpf(str(delay))+mp.pi/2)
        slope_precise = 1+mp.mpf(str(epsilon))

    def ray(x, chirality):
        if chirality == 'plus':
            s = inverse(x-half)
            _, r1, r2, r3, _ = clock(s)
            return s-delay-half, 1/r1, -r2/r1**3, 3*r2*r2/r1**5-r3/r1**4
        if x >= half+delay+width:
            with mp.workdps(40):
                return fixed_precise-slope_precise*(fixed_precise-mp.mpf(x)), slope, 0., 0.
        r, r1, r2, r3, _ = clock(x-half-delay)
        return r-half, r1, r2, r3

    def transported(x, chirality):
        if chirality == 'minus' and x >= fixed:
            raise ValueError('Retarded initial state is not continued through the fixed point')
        initial = half if chirality == 'plus' else half+delay
        y, log_derivative, scaled_schwarzian, history = x, 0., 0., []
        if chirality == 'minus' and y >= half+delay+width:
            with mp.workdps(40):
                gap = float(fixed_precise-mp.mpf(y))
                affine_limit = float(fixed_precise-mp.pi/2-mp.mpf(str(delay))-width)
            # Keep delta itself between returns; reconstructing xstar-delta at
            # every step loses digits that the blueshift then amplifies.
            while gap <= affine_limit:
                gap *= slope
                history.append((slope, 0.))
                log_derivative += math.log(slope)
                assert len(history) < 10000
            y = float(fixed_precise)-gap
        while y > initial:
            old, first, second, third = ray(float(y), chirality)
            assert old < y and first > 0
            schwarzian = third/first-1.5*(second/first)**2
            history.append((first, schwarzian))
            log_derivative += math.log(first)
            scaled_schwarzian = (scaled_schwarzian+schwarzian)/first**2
            y = old
            assert len(history) < 10000
        derivative = math.exp(log_derivative)
        stress = derivative**2*(f0-scaled_schwarzian/(24*math.pi))
        direct = f0
        for first, schwarzian in reversed(history):
            direct = first**2*direct-schwarzian/(24*math.pi)
        error = abs(stress-direct)/max(1., abs(stress), abs(direct))
        assert error < 2e-12
        return dict(coordinate=float(x), uniformizer=float(y+len(history)*length),
                    uniformizer_derivative=derivative, backward_returns=len(history),
                    stress_over_q=stress, j_over_alpha=1/6+8*math.pi*stress,
                    direct_recurrence_scaled_error=error)

    # The Schwarzian boundary identity is nonlinear in the prescribed clock rate.
    s = sp.symbols('s', real=True)
    r = sp.Function('r')(s)
    v = sp.diff(r, s)
    K = sp.diff(r, s, 3)/v-sp.Rational(3, 2)*(sp.diff(r, s, 2)/v)**2+v**2/2
    expected = (sp.diff(r, s, 4)/v-4*sp.diff(r, s, 2)*sp.diff(r, s, 3)/v**2
                +3*sp.diff(r, s, 2)**3/v**3+v*sp.diff(r, s, 2))
    assert sp.simplify(sp.diff(K, s)-expected) == 0
    # b(r(s))/b0=r'(s); d/dtau=(1/r')d/ds. This checks the matching identity.
    d_tau = lambda f: sp.diff(f, s)/v
    assert sp.simplify(v**2*(d_tau(d_tau(d_tau(v)))+d_tau(v))-expected) == 0

    def balance(s):
        r, r1, r2, r3, r4 = clock(s)
        kd = r4/r1-4*r2*r3/r1**2+3*r2**3/r1**3+r1*r2
        rr = b0*kd+r1**2*(transported(r+half, 'plus')['j_over_alpha']
                          -transported(r-half, 'minus')['j_over_alpha'])
        ll = -(transported(s-half, 'plus')['j_over_alpha']
               -transported(s+half, 'minus')['j_over_alpha'])
        return dict(exterior_time=s, right_residual_over_alpha=rr,
                    left_residual_over_alpha=ll, right_clock_rate=r1,
                    imposed_right_clock_term_over_alpha=b0*kd)

    inverse_errors, ward_errors, past_errors = [], [], []
    for x in (-2., 0., .1, .5, 1., 2., 10., fixed-.01):
        inverse_errors.append(abs(inverse(clock(x)[0])-x))
        for ch in ('plus', 'minus'):
            p, first, second, third = ray(x, ch)
            forward = clock(p+half+delay)[0]+half if ch == 'plus' else inverse(p+half)+delay+half
            inverse_errors.append(abs(forward-x))
            state, old = transported(x, ch), transported(p, ch)
            ward_errors.append(abs(state['uniformizer']-old['uniformizer']-length))
            if x <= 0:
                past_errors.append(abs(state['stress_over_q']-f0))
    assert max(inverse_errors) < 5e-13 and max(ward_errors) < 5e-12, (max(inverse_errors), max(ward_errors))
    assert max(past_errors) == 0 and balance(-1.)['right_residual_over_alpha'] == 0
    assert abs(quad(lambda x: step(x)[0], 0, 1)[0]-.5) < 1e-14
    assert abs(inverse(entry_shortcut+math.pi)-entry_shortcut-delay) < 5e-13
    assert abs(inverse(entry_loop+math.pi)-entry_loop+delay) < 5e-13
    assert abs(ray(fixed, 'minus')[0]-fixed) < 5e-13

    # Independent high-precision stress propagation: differentiate the smooth step,
    # propagate raw stress backwards, and use no float uniformizer/jet recurrence.
    precision_rows = []
    with mp.workdps(60):
        ee, dd, ww = mp.mpf('.05'), mp.mpf('.1'), mp.mpf(1)
        hh, LL = mp.pi/2, mp.pi+dd
        def mp_step(z):
            return mp.exp(-1/z)/(mp.exp(-1/z)+mp.exp(-1/(1-z)))
        for coordinate in (half+delay+.25, half+delay+.7, fixed-.01, fixed-.00001):
            y, history = mp.mpf(coordinate), []
            while y > hh+dd:
                z = y-hh-dd
                if z >= ww:
                    old, first, second, third = (1+ee)*z-ee*ww/2-hh, 1+ee, mp.mpf(0), mp.mpf(0)
                else:
                    old = z+ee*ww*mp.quad(mp_step, [0, z/ww])-hh
                    first = 1+ee*mp_step(z/ww)
                    second = ee/ww*mp.diff(mp_step, z/ww)
                    third = ee/ww**2*mp.diff(mp_step, z/ww, 2)
                history.append((first, third/first-mp.mpf('1.5')*(second/first)**2))
                y = old
            value = -mp.pi/(12*LL**2)
            for first, schwarzian in reversed(history):
                value = first**2*value-schwarzian/(24*mp.pi)
            candidate = transported(coordinate, 'minus')['stress_over_q']
            error = abs(mp.mpf(candidate)-value)/max(mp.mpf(1), abs(value))
            assert error < mp.mpf('2e-9'), (coordinate, candidate, mp.nstr(value, 30), mp.nstr(error, 8))
            precision_rows.append(dict(coordinate=coordinate, backward_returns=len(history),
                                       stress_over_q_60_digits=mp.nstr(value, 50), scaled_error=float(error)))
        derivative_errors = []
        for point in (.15, .37, .5, .8):
            exact = [mp.diff(mp_step, mp.mpf(point), n) for n in range(4)]
            derivative_errors.append(max(float(abs(mp.mpf(a)-b)/max(1, abs(b)))
                                         for a, b in zip(step(point), exact)))
        assert max(derivative_errors) < 2e-13

    horizon_rows = []
    for gap in (.1, .01, .001, .0001, .00001):
        with mp.workdps(40):
            coordinate = fixed_precise-mp.mpf(str(gap))
        state = transported(coordinate, 'minus')
        state.update(gap_to_fixed_point=gap, gap_squared_stress_over_q=gap**2*state['stress_over_q'])
        horizon_rows.append(state)
    # Fixed log-phase geometric sequences must obey delta^2 T=constant exactly
    # in the affine region. Different phases need not share that constant.
    phase_rows = []
    for phase in (0., .2, .4, .6, .8):
        with mp.workdps(40):
            precise_gaps = [mp.mpf('.2')/slope_precise**(mp.mpf(str(phase))+n) for n in (0, 10, 40, 100)]
            coordinates = [fixed_precise-gap for gap in precise_gaps]
        gaps = [float(gap) for gap in precise_gaps]
        values = [gap**2*transported(x, 'minus')['stress_over_q'] for gap, x in zip(gaps, coordinates)]
        err = max(abs(value-values[0])/max(1., abs(values[0])) for value in values)
        assert err < 2e-9
        phase_rows.append(dict(log_phase=phase, gaps=gaps, scaled_stresses=values, relative_spread=err))
    times = np.unique(np.r_[np.linspace(0, exit_loop-.02, 241), .25, .5, .75,
                             entry_shortcut, exit_shortcut, 10., 30., 60.])
    balances = [balance(float(t)) for t in times]
    entries = np.linspace(0, entry_loop-.02, 321)
    arrivals = np.array([inverse(float(t)+math.pi) for t in entries])
    shortcut_state = transported(entry_shortcut+half, 'minus')
    assert shortcut_state['stress_over_q'] != f0
    assert max(abs(row['right_residual_over_alpha']) for row in balances) > 1
    return dict(
        model='Prescribed smooth physical sewing maps in fixed leading AdS2, common exterior delay and free NS CFT. This is an exact retarded state calculation within that ansatz, NOT a self-consistent driven wormhole solution',
        units='All times/null coordinates use one fixed dimensionless throat scale; stress is f/q, total JT source is j/alpha=1/6+8*pi*f/q. No seconds or calendar mapping is assigned',
        parameters=dict(epsilon=epsilon, transition_width=width, exterior_delay=delay,
                        fastest_exterior_delay_assumed=delay, initial_loop_length=length,
                        initial_stress_over_q=f0, b0_over_alpha=b0),
        drive='l(s)=s; r(s)=s+epsilon*H(s); H_prime=g(s/width). g=0 before0,1 after1, and exp(-1/x)/[exp(-1/x)+exp(-1/(1-x))] between. H=0 before0 and s-width/2 afterwidth. C-infinity, monotone, exactly static past and affine late branch',
        return_maps='p_plus(u)=l(r_inverse(u-pi/2)-D)-pi/2; p_minus(v)=r(l_inverse(v-pi/2)-D)-pi/2. One-way arrivals a_LR(s)=r_inverse(l(s)+pi), a_RL(s)=l_inverse(r(s)+pi)',
        transported_state='F(p(x))=F(x)-L0; F=x in static past. F=p^N(x)+N*L0, using enough backward returns to reach that past. f/q=-pi*F_prime^2/(12*L0^2)-Schwarzian(F,x)/(24*pi). This retains the state history; it does not select a new instantaneous vacuum',
        thresholds=dict(right_exit_time_first_shortcut=exit_shortcut, left_entry_time_first_shortcut=entry_shortcut,
                        right_exit_time_first_closed_null_loop=exit_loop, left_entry_time_first_closed_null_loop=entry_loop,
                        fixed_minus_coordinate=fixed, shortcut_to_fixed_null_coordinate_gap=fixed-entry_shortcut-half),
        source_at_shortcut=shortcut_state, formal_near_fixed_point=horizon_rows,
        affine_tail='p_minus=xstar+(1+epsilon)*(x-xstar). delta^2*f(xstar-delta) is periodic in log_(1+epsilon)(delta); the transition fixes its phase. Some negative divergent geometric sequence is forced for this transported NS vacuum, but a universal constant coefficient/sign at all phases is not assumed',
        affine_phase_checks=phase_rows,
        boundary_balance='K_r=Schwarzian(r,s)+r_prime^2/2. R_R=b0*K_r_prime+r_prime^2*(j_plus(r+pi/2)-j_minus(r-pi/2)); R_L=b0*K_l_prime-l_prime^2*(j_plus(l-pi/2)-j_minus(l+pi/2)). Residuals must vanish for unforced leading matching r_prime=bR(r)/b0,l_prime=bL(l)/b0 with the same static past',
        boundary_samples=[balance(t) for t in (-1., .25, .5, .75, 1., 10., 30., entry_shortcut, exit_shortcut, exit_loop-.02)],
        verification=dict(maximum_inverse_error=float(max(inverse_errors)), maximum_uniformizer_ward_error=max(ward_errors),
                          static_past_error=max(past_errors), boundary_symbolic_residuals=['0', '0'],
                          step_derivative_errors=derivative_errors, independent_60_digit_checks=precision_rows),
        plot=dict(entry_times=entries.tolist(), one_way_times=(arrivals-entries).tolist(),
                  complete_return_times=(arrivals+delay-entries).tolist(),
                  boundary_times=times.tolist(), right_residual_over_alpha=[r['right_residual_over_alpha'] for r in balances],
                  left_residual_over_alpha=[r['left_residual_over_alpha'] for r in balances]),
        conclusion='The imposed clock drive does not solve its own leading JT balance. Nonzero residual quantifies a missing boundary driving term; no realizable external source has been supplied. Near-fixed-point stress growth is a diagnostic on this prescribed background, not evolution through a horizon or proof of chronology protection',
        scope='No new topology, full4D exterior dynamics, finite-cutoff matching, interacting gauge-sector state, local stress fluctuations, material preparation, or stable transmission. Large repeated blueshift leaves the low-energy model; even earlier formal points are not a solved coupled geometry. Finite-cutoff corrections include moving-boundary velocity. Common channel delay need not equal fastest exterior route in an actual wormhole',
        sources=['https://arxiv.org/pdf/math-ph/0412028 (Fewster-Hollands Eq2.9: chiral stress transformation)',
                 'https://arxiv.org/pdf/1807.04726 (MMP Eq5.24 and AppendixF: vacuum normalization and geometric anomaly)',
                 'https://arxiv.org/pdf/1606.01857 (Maldacena-Stanford-Yang: fixed induced metric/dilaton cutoff data; residual and smooth-drive computation are our derivation)'])


def boundary_feedback_checks(alpha, epsilon, width):
    """Our reduced transparent-loop closure; exterior lapse/flight delay are fixed."""
    u, v, chi = sp.symbols('u v chi', real=True)
    A, B = sp.Function('A')(u), sp.Function('B')(v)
    du = lambda f: sp.diff(f, u)+sp.diff(f, chi)/2
    dv = lambda f: sp.diff(f, v)-sp.diff(f, chi)/2
    phi = (A-B)*sp.tan(chi)+sp.diff(A, u)+sp.diff(B, v)
    residuals = [sp.simplify(sp.trigsimp(f)) for f in (
        du(du(phi))-sp.tan(chi)*du(phi)-sp.diff(A, u, 3)-sp.diff(A, u),
        dv(dv(phi))+sp.tan(chi)*dv(phi)-sp.diff(B, v, 3)-sp.diff(B, v),
        du(dv(phi))+phi/(2*sp.cos(chi)**2))]
    assert residuals == [0, 0, 0]
    impulse = math.pi*epsilon/2  # J/alpha; energy-fraction calibration applies at D=0.
    def injection(x):
        return impulse*8/(3*width)*math.sin(math.pi*x/width)**4 if 0<x<width else 0.
    end = math.pi
    samples = np.unique(np.r_[np.linspace(0, end, 301), 1., math.pi/2+width/2])

    def solve(delay, max_step):
        omega2 = (2*math.pi/(math.pi+delay))**2
        coefficient = 2/(math.pi*(omega2-1))
        segments, ends = [], []
        def evaluate(x):
            if x <= 0:
                return np.zeros(6)
            assert segments and x <= ends[-1]+1e-12
            index = min(int(np.searchsorted(ends, x)), len(segments)-1)
            return segments[index].sol(min(x, ends[-1]))
        def rhs(x, y):
            # ponytail: first throat interval only; pi and full-loop histories are still zero.
            old = y if delay == 0 else evaluate(x-delay)
            qs, qc = y[2]+omega2*y[0], y[5]+omega2*y[3]
            old_qs, old_qc = old[2]+omega2*old[0], old[5]+omega2*old[3]
            return np.array((y[1], y[2], -y[1]-coefficient*(qs-old_qs)-injection(x),
                             y[4], y[5], -y[4]-coefficient*(qc+old_qc)-injection(x)))
        initial, start = np.zeros(6), 0.
        while start < end:
            stop = end if delay == 0 else min((len(segments)+1)*delay, end)
            segment = solve_ivp(rhs, (start, stop), initial, dense_output=True,
                                method='DOP853', rtol=2e-12, atol=2e-14, max_step=max_step)
            assert segment.success
            segments.append(segment); ends.append(stop)
            initial, start = segment.y[:, -1], stop
        return evaluate, rhs, omega2, coefficient

    records, solutions = {}, {}
    for delay in (0., .1):
        evaluate, rhs, omega2, coefficient = solve(delay, width/40)
        refined, _, _, _ = solve(delay, width/80)
        values = np.array([evaluate(x) for x in samples]).T
        refined_values = np.array([refined(x) for x in samples]).T
        refinement_error = float(np.max(np.abs(values-refined_values)))
        assert refinement_error < 1e-10
        beta_over_alpha = (omega2-1)/3
        b0_over_alpha = math.pi*beta_over_alpha/2
        def observables(x):
            y = evaluate(x); old = evaluate(x-delay)
            avec, bvec = (y[:3]+y[3:])/2, (y[:3]-y[3:])/2
            ad, bd = (old[:3]+old[3:])/2, (old[:3]-old[3:])/2
            qp = avec[2]+omega2*avec[0]-(bd[2]+omega2*bd[0])
            qm = bvec[2]+omega2*bvec[0]-(ad[2]+omega2*ad[0])
            jp, jm = injection(x)+coefficient*qp, coefficient*qm
            derivative = rhs(x, y)
            constraint_error = max(abs((derivative[2]+derivative[5]+y[1]+y[4])/2+jp),
                                   abs((derivative[2]-derivative[5]+y[1]-y[4])/2+jm))
            # eta' from inverse transparent ray maps; stress changes follow -alpha/3*Q eta'.
            eta_plus = (bd[0]-avec[0])/b0_over_alpha
            eta_minus = (ad[0]-bvec[0])/b0_over_alpha
            eta_plus_third = (bd[2]-avec[2])/b0_over_alpha
            eta_minus_third = (ad[2]-bvec[2])/b0_over_alpha
            ward_error = max(abs(jp-injection(x)+(eta_plus_third+omega2*eta_plus)/3),
                             abs(jm+(eta_minus_third+omega2*eta_minus)/3))
            return dict(x=x, exterior_time=x-math.pi/2,
                        S_prime_over_alpha=y[1], C_over_alpha=y[3],
                        j_plus_over_alpha=jp, j_minus_over_alpha=jm,
                        right_clock_rate=1+avec[0]/b0_over_alpha,
                        left_clock_rate=1+bvec[0]/b0_over_alpha,
                        inverse_ray_plus_rate=1+eta_plus, inverse_ray_minus_rate=1+eta_minus,
                        constraint_error=constraint_error, stress_return_error=ward_error)
        measured = [observables(x) for x in samples]
        identity_error = max(max(m['constraint_error'], m['stress_return_error']) for m in measured)
        assert identity_error < 1e-12
        rate_keys = ('right_clock_rate', 'left_clock_rate', 'inverse_ray_plus_rate', 'inverse_ray_minus_rate')
        minimum_rates = {key:min(m[key] for m in measured) for key in rate_keys}
        assert min(minimum_rates.values())>0
        if delay:
            before = np.array([evaluate(x) for x in np.linspace(0, delay, 21)]).T
            causal_error = float(np.max(np.abs(before[:3]-before[3:])))
            assert causal_error < 1e-12  # B and its first two derivatives vanish before exterior return.
        else:
            causal_error = None
            assert max(abs(m['j_plus_over_alpha']+m['j_minus_over_alpha']-injection(m['x']))
                       for m in measured)<1e-12
        records[str(delay)] = dict(exterior_delay=delay, loop_length=math.pi+delay,
            omega_squared=omega2, beta_over_alpha=beta_over_alpha, c0=coefficient,
            max_step_halving_error=refinement_error, maximum_identity_error=identity_error,
            early_opposite_chirality_error=causal_error, minimum_sampled_rates=minimum_rates,
            at_x_one=observables(1.), at_original_center_peak=observables(math.pi/2+width/2),
            plot=dict(x=samples.tolist(), center_delta_phi_over_alpha=values[1].tolist(),
                      j_plus_over_alpha=[m['j_plus_over_alpha'] for m in measured],
                      j_minus_over_alpha=[m['j_minus_over_alpha'] for m in measured],
                      right_clock_rate=[m['right_clock_rate'] for m in measured],
                      left_clock_rate=[m['left_clock_rate'] for m in measured]))
        solutions[delay] = evaluate

    k = 4/(3*math.pi)
    matrix_s = np.array(((0, 1, 0), (0, 0, 1), (0, -1, 0)), dtype=float)
    matrix_c = np.array(((0, 1, 0), (0, 0, 1), (-4*k, -1, -k)), dtype=float)
    convolution_errors = []
    for x in (.1, width, 1., math.pi/2+width/2, end):
        expected = np.array([-quad(lambda t: injection(t)*expm(matrix*(x-t))[i, 2],
                                   0, min(width, x), epsabs=1e-13, epsrel=1e-12)[0]
                             for matrix in (matrix_s, matrix_c) for i in range(3)])
        convolution_errors.append(float(np.max(np.abs(solutions[0.](x)-expected))))
    assert max(convolution_errors)<1e-10
    return dict(
        model='Our linear leading-JT mean-field closure with transparent CFT circulation, fixed exterior proper lapse and a common exterior delay D; first throat interval only. Not complete MMP four-dimensional evolution',
        assumptions='Retarded zero perturbation history; assumed initial antiperiodic cylinder vacuum; small pulse. All channels have the same fixed D, unlike actual channel-dependent MMP exterior paths. Physical cutoff curves have fixed dilaton and induced proper metric. Their overlap limit is used',
        bulk='delta_phi=[A(u)-B(v)]*tan(chi)+A_prime(u)+B_prime(v); P A=-j_plus, P B=-j_minus, P=d^3/dx^3+d/dx; S=A+B, C=A-B. Numerics integrate S,C divided by alpha',
        symbolic_bulk_residuals=[str(r) for r in residuals],
        clocks='tau_R,L(s)=s+xi_R,L(s); xi_R_prime=[A(s+pi/2)-B(s-pi/2)]/b0; xi_L_prime=[-A(s-pi/2)+B(s+pi/2)]/b0. Reported rates use s=x-pi/2, not s=x',
        state='Transform only the flat normal-ordered chiral state function f: f(x)=p_prime^2*f(p)-q*Schwarzian(p,x)/(24*pi). Add local geometric anomaly q/(48*pi) afterward. f0=-q*pi/[12*(pi+D)^2]',
        maps='p_plus(u)=tau_L(tau_R_inverse(u-pi/2)-D)-pi/2; p_minus swaps R,L. p=x-(pi+D)+eta. Linear delta_f(x)-delta_f(x-pi-D)=-q/(24*pi)*(eta_tripleprime+Omega^2*eta_prime)',
        background='Omega^2=[2*pi/(pi+D)]^2; beta=alpha*(Omega^2-1)/3; phi0=beta*(1+chi*tan(chi)); b0=pi*beta/2; c0=2/[pi*(Omega^2-1)]. beta>0 requires D<pi',
        first_interval_equations=['S_tripleprime=-S_prime-c0*(Q S-Q S_delayed_D)-i_sum',
                                  'C_tripleprime=-C_prime-c0*(Q C+Q C_delayed_D)-i_diff; Q=d^2/dx^2+Omega^2'],
        forcing='i_plus=j_packet, i_minus=0, fixed J/alpha=pi*.01/2 across D. epsilon=.01 is the D=0 binding-energy calibration only. New circulating injection and a separate one-pass carrier agree during this first interval; they differ after a full return',
        parameters=dict(alpha=alpha, epsilon_D0=epsilon, width=width, final_x=end),
        cases=records, D0_independent_matrix_convolution_errors=convolution_errors,
        D0_center_result='S_tripleprime+S_prime=-i: the first-interval center area matches the previous fixed-source calculation exactly, although the chiral stresses and mouth clocks change',
        scope='Sampled positive clock/ray rates are necessary linear admissibility checks, not a certified global quantum state. No complete return interval, long-time stability, quantum fluctuations, exterior gravity, formation, physical traversal or past-date mapping is established',
        sources=['https://arxiv.org/pdf/1807.04726 (MMP sections5.24-5.27 and AppendixF: background and anomaly)',
                 'https://arxiv.org/pdf/math-ph/0412028 (Fewster-Hollands equations2.7-2.9: CFT transformation law)',
                 'https://arxiv.org/pdf/1606.01857 (Maldacena-Stanford-Yang: fixed proper metric/dilaton JT cutoff data)'])


def loop_return_checks(alpha, epsilon, width, first_interval):
    """Retarded chiral-source recurrence through eight loops, with two pulse conventions."""
    def pulse(x):
        return math.pi*epsilon/2*8/(3*width)*math.sin(math.pi*x/width)**4 if 0<x<width else 0.

    def solve(delay, mode, step):
        length = math.pi+delay
        omega2 = (2*math.pi/length)**2
        c0 = 2/(math.pi*(omega2-1))
        ends, segments = [], []
        def history(x):
            if x <= 0:
                return np.zeros(8)
            assert segments and x <= ends[-1]+1e-11
            index = min(int(np.searchsorted(ends, x)), len(segments)-1)
            return segments[index].sol(min(x, ends[-1]))
        def q(y):
            return np.array((y[2]+omega2*y[0], y[5]+omega2*y[3]))
        def sources(x, y):
            # Sum the causal stress recurrence; no numerical third-derivative history is needed.
            local = q(y)
            past, other = np.zeros(2), np.zeros(2)
            driving = pulse(x)
            count = int(x/length)+1
            for k in range(count):
                old = x-(k+1)*length
                past += q(history(old))
                for offset in (delay, math.pi):
                    argument = x-k*length-offset
                    other += q(y if k == 0 and offset == 0 else history(argument))
                if mode == 'circulating' and k:
                    driving += pulse(x-k*length)
            flux = c0*(local+2*past-other[::-1])
            flux[0] += driving
            return flux
        def rhs(x, y):
            jp, jm = sources(x, y)
            return (y[1], y[2], -y[1]-jp, y[4], y[5], -y[4]-jm, y[0], y[3])
        interval = delay if delay else math.pi
        start, initial, final = 0., np.zeros(8), 8*length
        while start < final:
            stop = min((len(segments)+1)*interval, final)
            segment = solve_ivp(rhs, (start, stop), initial, dense_output=True,
                                method='DOP853', rtol=2e-11, atol=2e-13, max_step=step)
            assert segment.success
            segments.append(segment); ends.append(stop)
            initial, start = segment.y[:, -1], stop
        return history, sources, q, length, omega2

    records = {}
    for delay in (0., .1):
        for mode in ('one_pass', 'circulating'):
            history, sources, q, length, omega2 = solve(delay, mode, width/12)
            refined, _, _, _, _ = solve(delay, mode, width/24)
            samples = np.unique(np.r_[np.linspace(0, 8*length, 801),
                                      np.arange(9)*length, math.pi/2+width/2])
            values = np.array([history(x) for x in samples]).T
            fine = np.array([refined(x) for x in samples]).T
            refinement = float(np.max(np.abs(values-fine)))
            assert refinement < 2e-8
            b0_over_alpha = math.pi*(omega2-1)/6
            c0 = 2/(math.pi*(omega2-1))
            first = first_interval['cases'][str(delay)]['at_original_center_peak']
            at_peak = history(math.pi/2+width/2)
            first_error = abs(at_peak[1]+at_peak[4]-first['S_prime_over_alpha'])
            assert first_error < 2e-9
            clock_rates, ray_rates, ray_gaps, fluxes, residuals = [], [], [], [], []
            center_derivative = values[1]+values[4]
            for x, y in zip(samples, values.T):
                old_pi, old_d, old_l = history(x-math.pi), history(x-delay), history(x-length)
                clock = np.array((y[0]-old_pi[3], y[3]-old_pi[0]))/b0_over_alpha
                eta_prime = (old_d[[3, 0]]+old_pi[[3, 0]]-y[[0, 3]]-old_l[[0, 3]])/b0_over_alpha
                eta_third = (old_d[[5, 2]]+old_pi[[5, 2]]-y[[2, 5]]-old_l[[2, 5]])/b0_over_alpha
                eta = (old_d[[7, 6]]+old_pi[[7, 6]]-y[[6, 7]]-old_l[[6, 7]])/b0_over_alpha
                flux = sources(x, y)
                old_flux = sources(x-length, old_l) if x>length else np.zeros(2)
                new_injection = pulse(x)-(pulse(x-length) if mode=='one_pass' else 0.)
                residual = flux-old_flux+(eta_third+omega2*eta_prime)/3
                residual[0] -= new_injection
                residuals.append(float(np.max(np.abs(residual))))
                clock_rates.append(1+clock); ray_rates.append(1+eta_prime)
                ray_gaps.append(length-eta); fluxes.append(flux)
            clocks, rays, flux = np.array(clock_rates), np.array(ray_rates), np.array(fluxes)
            assert max(residuals)<1e-10
            if delay == 0:
                # At D=0, center S solves P S = -p for a one-pass carrier,
                # or -sum p(x-k*pi) for a circulating injection. Check by direct convolution.
                errors = []
                for x in (1., length+.7, 3*length+.7, 8*length):
                    expected = 0.
                    for k in range(1 if mode=='one_pass' else int(x/length)+1):
                        local_x = x-k*length
                        if local_x>0:
                            expected -= quad(lambda t: pulse(t)*math.sin(local_x-t),
                                             0, min(width, local_x), epsabs=1e-13)[0]
                    y = history(x)
                    errors.append(abs(y[1]+y[4]-expected))
                assert max(errors)<2e-8
            else:
                errors = None
            blocks = []
            for k in range(8):
                mask = (samples>=k*length)&(samples<=(k+1)*length)
                blocks.append(dict(loop=k+1, max_center_area_excess_fraction=float(np.max(np.abs(center_derivative[mask]))*3/(omega2-1)),
                                   max_clock_rate_change=float(np.max(np.abs(clocks[mask]-1))),
                                   max_inverse_ray_rate_change=float(np.max(np.abs(rays[mask]-1)))))
            # Reconstruct the field: raw A/B coefficient growth can cancel in delta_phi.
            interior_max, total_area_max = 0., 0.
            for tau in np.linspace(-.8, 8*length-.8, 301):
                for chi in np.linspace(-.8, .8, 17):
                    aa, bb = history(tau+chi), history(tau-chi)
                    delta = (aa[0]-bb[3])*math.tan(chi)+aa[1]+bb[4]
                    background = (omega2-1)/3*(1+chi*math.tan(chi))
                    interior_max = max(interior_max, abs(delta/background))
                    total_area_max = max(total_area_max, abs(alpha*delta/(1+alpha*background)))
            records[f'{delay}_{mode}'] = dict(exterior_delay=delay, forcing=mode,
                loops=8, final_x=8*length, maximum_step_halving_error=refinement,
                first_interval_center_error=first_error, maximum_stress_recurrence_error=max(residuals),
                D0_center_convolution_errors=errors, per_loop=blocks,
                minimum_sampled_clock_rate=float(np.min(clocks)), minimum_sampled_inverse_ray_rate=float(np.min(rays)),
                minimum_sampled_inverse_ray_time_gap=float(np.min(ray_gaps)),
                maximum_normalized_bulk_coefficient=float(np.max(np.abs(values[:6]))),
                interior_sample=dict(chi_min=-.8, chi_max=.8, chi_points=17, time_points=301,
                    maximum_change_relative_to_background_area_excess=interior_max,
                    maximum_fractional_total_sphere_area_change=total_area_max),
                plot=dict(x=samples.tolist(), loop_number=(samples/length).tolist(),
                          center_delta_phi_over_alpha=center_derivative.tolist(),
                          right_clock_rate=clocks[:, 0].tolist(), left_clock_rate=clocks[:, 1].tolist(),
                          j_plus_over_alpha=flux[:, 0].tolist(), j_minus_over_alpha=flux[:, 1].tolist()))
    return dict(
        model='Our linear transparent-loop JT/CFT closure extended through eight background loop lengths, fixed exterior lapse/common delay; not full 4D evolution',
        parameters=dict(alpha=alpha, epsilon_D0=epsilon, width=width),
        recurrence='j_plus-j_plus_L=i_plus+c0*Q(A+A_L-B_D-B_pi); j_minus-j_minus_L=i_minus+c0*Q(B+B_L-A_D-A_pi); P A=-j_plus, P B=-j_minus. Numerics use A,B divided by alpha',
        conventions='one_pass: i_plus=p(x)-p(x-L), separate carrier removed after one throat passage; circulating: i_plus=p(x), one injection subsequently circulates. i_minus=0. Both include induced circulating CFT response',
        clock_observables='Rates refer to exterior s=x-pi/2. Full delayed A,B contributions are retained after x=pi. The background is chronology-preserving; sampled perturbed rates and return gaps are checked in finite linear runs. A complete signal loop is not a closed timelike curve',
        cases=records,
        scope='Finite-time linear response only. Positive sampled rates and bounded center response alone do not prove nonlinear stability, global state positivity or a physical wormhole. Secular coordinate coefficients must be distinguished from invariant geometry. No formation, actual apparatus, backward calendar mapping or human/organic transit is demonstrated')


def _distributed_clock_response(delays,weights,step,end,pulse,scalar_reference=False,exterior_motion=None):
    """Retarded linear clock response; optional motion returns h,h',h'',h''' arrays."""
    assert not scalar_reference or exterior_motion is None
    delays=np.asarray(delays,float);weights=np.asarray(weights,float)
    assert np.all(delays>0) and np.all(weights>0) and abs(weights.sum()-1)<1e-13
    lengths=math.pi+delays;omega2=(2*math.pi/lengths)**2
    mean=float(weights@omega2);b0=math.pi*(mean-1)/6;c0=1/(3*b0)
    assert b0>0
    ends,segments=[],[]
    def history(arguments):
        arguments=np.asarray(arguments,float);flat=arguments.ravel()
        result=np.zeros((len(flat),8));active=flat>0
        if np.any(active):
            assert ends and np.max(flat[active])<=ends[-1]+2e-11
            times=np.minimum(flat[active],ends[-1]);indices=np.searchsorted(ends,times)
            positions=np.flatnonzero(active)
            for index in np.unique(indices):
                select=indices==index
                result[positions[select]]=segments[int(index)].sol(times[select]).T
        return result.reshape(arguments.shape+(8,))
    def channels(x,y):
        if x<=0:return np.zeros((len(delays),2))
        count=int(x/lengths.min())
        local=np.column_stack((y[2]+omega2*y[0],y[5]+omega2*y[3]))
        past=np.zeros_like(local)
        if count:
            old=history(x-np.arange(1,count+1)[:,None]*lengths)
            past=np.stack((old[...,2]+omega2*old[...,0],
                           old[...,5]+omega2*old[...,3]),axis=-1).sum(axis=0)
        offsets=np.arange(count+1)[:,None]*lengths
        old=history(np.stack((x-offsets-delays,x-offsets-math.pi)))
        other=np.stack((old[...,2]+omega2*old[...,0],
                        old[...,5]+omega2*old[...,3]),axis=-1).sum(axis=(0,1))
        result=c0*(local+2*past-other[:,::-1])
        if exterior_motion is not None:
            # Fixed A receives the delayed motion; moving B receives it locally.
            arguments=np.stack((x-offsets-math.pi/2-delays,x-offsets-math.pi/2))
            motion=exterior_motion(arguments)
            drive=(np.asarray(motion[3])+omega2*np.asarray(motion[1]))/3
            assert drive.shape==arguments.shape
            result+=drive.sum(axis=1).T
        return result
    def old_scalar_sources(x,y):
        # Literal one_pass source sum from loop_return_checks, including p.
        length=float(lengths[0]);frequency=float(omega2[0])
        def q(v):return np.array((v[2]+frequency*v[0],v[5]+frequency*v[3]))
        past,other=np.zeros(2),np.zeros(2)
        for k in range(int(x/length)+1):
            past+=q(history(x-(k+1)*length))
            for offset in (float(delays[0]),math.pi):
                other+=q(history(x-k*length-offset))
        result=c0*(q(y)+2*past-other[::-1]);result[0]+=pulse(x)
        return result
    def rhs(x,y):
        if scalar_reference:
            jp,jm=old_scalar_sources(x,y)
        else:
            jp,jm=weights@channels(x,y);jp+=pulse(x)
        return (y[1],y[2],-y[1]-jp,y[4],y[5],-y[4]-jm,y[0],y[3])
    interval=min(float(delays.min()),math.pi)
    start=0.;initial=np.zeros(8)
    while start<end:
        stop=min((len(segments)+1)*interval,end)
        solution=solve_ivp(rhs,(start,stop),initial,method='DOP853',
            rtol=2e-11,atol=2e-13,max_step=step,dense_output=True)
        assert solution.success
        ends.append(stop);segments.append(solution)
        start=stop;initial=solution.y[:,-1]
    return history,channels,dict(delays=delays,lengths=lengths,weights=weights,
        omega2=omega2,mean=mean,b0=b0,c0=c0,
        solver_evaluations=sum(s.nfev for s in segments),segments=len(segments))


def distributed_loop_response_checks(progress=None, orders=(64,128), delay_panels=False):
    """Finite-time leading-JT response with unequal, fixed channel delays."""
    import hashlib
    source_path=Path(__file__).with_name('polyakov_closure_results.json')
    saved=json.loads(source_path.read_text(encoding='utf-8'))
    pulse_epsilon,width,loops=1e-4,.4,4
    pulse_area=math.pi*pulse_epsilon/2
    reference=saved['quadrature_refinement'][-1]['roots'][-1]
    assert len(orders)==2 and all(isinstance(n,int) and n>0 for n in orders) and orders[0]<orders[1]
    shortest_delay=float(reference['shortest_exterior_time'])
    final=loops*(math.pi+shortest_delay)
    def pulse(x):
        return pulse_area*8/(3*width)*math.sin(math.pi*x/width)**4 if 0<x<width else 0.

    def solve(delays,weights,step,end,scalar_reference=False):
        return _distributed_clock_response(delays,weights,step,end,pulse,scalar_reference)

    samples=np.unique(np.r_[np.linspace(0,final,401),width,
                            np.arange(loops+1)*(math.pi+shortest_delay)])
    def record(order,step):
        state=None if delay_panels else next((v['roots'][-1] for v in saved['quadrature_refinement'] if v['order']==order),None)
        path_checks=None
        if state is not None:
            delays=np.asarray(state['profile']['periods'])-state['throat_optical_time']
            weights=np.asarray(state['profile']['weights'])
        else:
            from polyakov_response import _magnetic_closure_paths
            distance=saved['parameters']['separation_over_re']
            quadrature=None
            if delay_panels:
                # Saved interpolation selects integration panels only. Every new
                # path is integrated directly; no interpolated source is evolved.
                refnu=np.asarray(reference['profile']['nu'])
                refdelay=np.asarray(reference['profile']['periods'])-reference['throat_optical_time']
                assert np.all(np.diff(refdelay)<0)
                targets=np.arange(shortest_delay+width,final+2*width,width)
                boundaries=np.r_[0.,np.sort(np.interp(targets,refdelay[::-1],refnu[::-1])),2.]
                assert np.all(np.diff(boundaries)>0)
                nodes,gauss_weights=np.polynomial.legendre.leggauss(order)
                lengths_of_panels=np.diff(boundaries)
                nu=((boundaries[:-1]+boundaries[1:])[:,None]/2+lengths_of_panels[:,None]*nodes/2).ravel()
                weights=(lengths_of_panels[:,None]*gauss_weights/4).ravel()
                quadrature=(nu,weights)
            _,weights,lengths,pieces,path_checks=_magnetic_closure_paths(order,distance,(100.,),quadrature=quadrature)
            piece=pieces['100.0']
            delays=(distance*lengths-2*piece['isotropic']+piece['shapiro'])/reference['clock_scale']
        assert delays.min()>=shortest_delay-1e-12
        history,channels,p=solve(delays,weights,step,final)
        values=history(samples);oldpi=history(samples-math.pi)
        b0=p['b0']
        clocks=1+(values[:,[0,3]]-oldpi[:,[3,0]])/b0
        clock_displacements=(values[:,[6,7]]-oldpi[:,[7,6]])/b0
        stress=np.array([channels(float(x),y) for x,y in zip(samples,values)])
        aggregate=np.einsum('n,tnc->tc',weights,stress)
        baseline=(1-p['mean'])/6
        # Entry coordinate x=s+pi/2; evaluate arrivals only inside solved history.
        entries=np.linspace(0,final-math.pi,401)
        now=history(entries);past=history(entries-math.pi);future=history(entries+math.pi)
        crossing=math.pi+(2*now[:,[7,6]]-past[:,[6,7]]-future[:,[6,7]])/b0
        full_return=crossing+shortest_delay
        ray_min=math.inf;recurrence_error=0.
        for x,y in zip(samples,values):
            oldd=history(x-delays);oldl=history(x-p['lengths']);op=history(x-math.pi)
            eta_prime=(oldd[:,[3,0]]+op[[3,0]]-y[[0,3]]-oldl[:,[0,3]])/b0
            ray_min=min(ray_min,float(np.min(1+eta_prime)))
        # Independent channel-by-channel recurrence at selected times and labels.
        def single_flux(time,index):
            if time<=0:return np.zeros(2)
            om=p['omega2'][index];length=p['lengths'][index];delay=delays[index]
            def q(v):return np.array([v[2]+om*v[0],v[5]+om*v[3]])
            result=q(history(time))
            for k in range(int(time/length)+1):
                result+=2*q(history(time-(k+1)*length))
                result-=q(history(time-k*length-delay))[::-1]
                result-=q(history(time-k*length-math.pi))[::-1]
            return p['c0']*result
        for x in np.linspace(.07,final,13):
            y=history(x);all_flux=channels(float(x),y)
            for index in (0,len(delays)//2,len(delays)-1):
                om=p['omega2'][index];oldl=history(x-p['lengths'][index])
                oldd=history(x-delays[index]);oldp=history(x-math.pi)
                def q(v):return np.array([v[2]+om*v[0],v[5]+om*v[3]])
                expected=p['c0']*(q(y)+q(oldl)-q(oldd)[::-1]-q(oldp)[::-1])
                residual=all_flux[index]-single_flux(x-p['lengths'][index],index)-expected
                recurrence_error=max(recurrence_error,float(np.max(abs(residual))))
        center=values[:,1]+values[:,4]
        center_relative=3*center/(p['mean']-1)
        max_clock=float(np.max(abs(clocks-1)))
        max_transit=float(np.max(abs(crossing-math.pi))/math.pi)
        assert recurrence_error<1e-10 and max_clock<.02 and max_transit<.02
        result=dict(flux_order=order,maximum_step=step,final_coordinate=final,
            path_checks=path_checks,
            channel_count=len(delays),quadrature_scheme='arrival-time panels in flux; order per panel' if delay_panels else 'global flux Gauss',
            mean_omega_squared=p['mean'],b0_over_alpha=b0,c0=p['c0'],
            shortest_exterior_delay=shortest_delay,
            sampled_channel_delay_range=[float(delays.min()),float(delays.max())],
            background_mean_j_over_alpha=baseline,
            maximum_absolute_mean_stress_response=float(np.max(abs(aggregate))),
            maximum_absolute_channel_stress_response=float(np.max(abs(stress))),
            maximum_absolute_center_area_excess_response_fraction=float(np.max(abs(center_relative))),
            maximum_sampled_clock_rate_change=max_clock,
            minimum_sampled_clock_rate=float(clocks.min()),
            minimum_sampled_inverse_ray_rate=ray_min,
            maximum_sampled_linearized_throat_time_fractional_change=max_transit,
            minimum_sampled_linearized_throat_times=crossing.min(axis=0).tolist(),
            minimum_sampled_linearized_full_return_times=full_return.min(axis=0).tolist(),
            maximum_stress_recurrence_residual=recurrence_error,
            solver_evaluations=p['solver_evaluations'],method_of_steps_segments=p['segments'],
            plot=dict(x=samples.tolist(),shortest_loop_number=(samples/(math.pi+shortest_delay)).tolist(),
                center_delta_phi_over_alpha=center.tolist(),
                mean_j_plus_response_over_alpha=aggregate[:,0].tolist(),
                mean_j_minus_response_over_alpha=aggregate[:,1].tolist(),
                mean_j_plus_over_alpha=(baseline+aggregate[:,0]).tolist(),
                mean_j_minus_over_alpha=(baseline+aggregate[:,1]).tolist(),
                external_carrier=[pulse(float(x)) for x in samples],
                right_clock_rate=clocks[:,0].tolist(),left_clock_rate=clocks[:,1].tolist(),
                right_clock_displacement=clock_displacements[:,0].tolist(),
                left_clock_displacement=clock_displacements[:,1].tolist(),
                transit_entry_exterior_time=(entries-math.pi/2).tolist(),
                linearized_left_to_right_time=crossing[:,0].tolist(),
                linearized_right_to_left_time=crossing[:,1].tolist()))
        comparison=np.c_[values,aggregate,clocks]
        return result,comparison

    records=[];comparisons=[]
    for order,step in ((orders[0],width/12),(orders[1],width/12),(orders[1],width/24)):
        row,comparison=record(order,step);records.append(row);comparisons.append(comparison)
        if progress is not None:progress(row)
    quadrature_error=float(np.max(abs(comparisons[0]-comparisons[1])))
    time_error=float(np.max(abs(comparisons[1]-comparisons[2])))
    scale=float(np.max(abs(comparisons[2][:,:10])))
    observable_keys=('center_delta_phi_over_alpha','mean_j_plus_response_over_alpha',
        'mean_j_minus_response_over_alpha','right_clock_rate','left_clock_rate',
        'right_clock_displacement','left_clock_displacement',
        'linearized_left_to_right_time','linearized_right_to_left_time')
    observable_errors={name:{key:float(np.max(abs(np.asarray(records[a]['plot'][key])-
        np.asarray(records[b]['plot'][key])))) for key in observable_keys}
        for name,a,b in ((f'flux_{orders[0]}_{orders[1]}',0,1),('timestep_halving',1,2))}
    converged=quadrature_error<2e-7 and time_error<2e-8
    degeneracy_error=None
    if converged:
        degeneracy_end=loops*(math.pi+.1)
        general,_,_=solve([.1],[1.],width/24,degeneracy_end)
        scalar,_,_=solve([.1],[1.],width/24,degeneracy_end,scalar_reference=True)
        times=np.linspace(0,degeneracy_end,257)
        degeneracy_error=float(np.max(abs(general(times)-scalar(times))))
        assert degeneracy_error<2e-9
    area_error=abs(quad(pulse,0,width,epsabs=1e-14)[0]-pulse_area)
    assert area_error<1e-14
    return dict(model='Finite-time distributed-delay linear JT/CFT response on a leading throat',
        input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        parameters=dict(alpha=saved['parameters']['alpha'],pulse_epsilon=pulse_epsilon,
            pulse_area=pulse_area,pulse_width=width,shortest_loops=loops,
            shortest_exterior_delay=shortest_delay,final_coordinate=final,flux_orders=list(orders),
            new_order_clock_scale=reference['clock_scale'],delay_panels=delay_panels,
            panel_target_delay_width=width if delay_panels else None),
        refinement=records,
        equations=dict(background='D_nu=T67_nu−T67_throat>0; L_nu=pi+D_nu; Omega_nu²=(2pi/L_nu)². Positive flux weights sum to1. b0/alpha=pi*(meanOmega²−1)/6; c0=1/(3*b0/alpha).',
            recurrence='jnu_plus−jnu_plus_L=c0*Qnu(A+A_L−B_D−B_pi), and plus/minus with A/B exchanged; Qnu=d²/dx²+Omega_nu². These j are deviations from the channel vacuum stress.',
            evolution='(d³/dx³+d/dx)A=−sum(weights*jnu_plus)−p; (d³/dx³+d/dx)B=−sum(weights*jnu_minus). A,B and stress responses are divided by alpha; integrals of A/B are retained.',
            carrier='p is an external small test carrier with area pi*1e−4/2, sin^4 profile and width .4. It crosses once and is removed; it is not injected once into every channel or subsequently recirculated. No physical preparation or new matter model is supplied.',
            clocks='For x=s+pi/2: delta r=(integral A(x)−integral B(x−pi))/b0, delta l=(integral B(x)−integral A(x−pi))/b0, with normalized b0. The stored rates are1+their derivatives.',
            transit='T_LR(s)=pi+delta l(s)−delta r(s+pi), T_RL(s)=pi+delta r(s)−delta l(s+pi), to first order only. The shortest complete return adds the fixed shortest exterior delay. Sampled entries require both endpoints inside the solved history.'),
        verification=dict(all_declared_checks_pass=converged,
            status='CONVERGED' if converged else 'UNCONVERGED: supplementary numerical response only',
            acceptance_thresholds={f'flux_{orders[0]}_{orders[1]}_absolute':2e-7,'timestep_halving_absolute':2e-8},
            flux_convergence_pass=quadrature_error<2e-7,timestep_convergence_pass=time_error<2e-8,
            **{f'flux_{orders[0]}_{orders[1]}_maximum_absolute_difference':quadrature_error},
            timestep_halving_maximum_absolute_difference=time_error,comparison_response_scale=scale,
            observable_maximum_absolute_differences=observable_errors,
            flux_error_over_peak_mean_stress=quadrature_error/records[1]['maximum_absolute_mean_stress_response'],
            one_channel_D01_old_source_maximum_state_difference=degeneracy_error,
            one_channel_check_status='passed' if converged else 'not run after distributed convergence failure',
            pulse_area_absolute_error=area_error,
            comparison='Common time samples of A/B jets, primitives, aggregate stress and real mouth-clock rates. The separate literal-old-source one-channel solve is skipped if distributed convergence fails; no skipped check is recorded as passed.'),
        scope='Four shortest-loop times, small external pulse and a linear leading-pi-throat model only. The67 exterior delays/weights seed this new model; it is not full finite-collar dynamics and does not evolve angular geometry, the physical exterior, finite-q states or interface conditions. Long-delay channels may not complete one loop. Sampled positive clocks and linearized transit/return times are finite-run diagnostics, not exact nonlinear travel times, a global stability theorem, quantum-state preparation, formation or backward time travel.',
        sources=['https://arxiv.org/html/1807.04726v3: assumed magnetic-channel/leading-throat framework',
            'https://arxiv.org/abs/math-ph/0412028: chiral conformal stress transformation',
            'Our distributed extension of the existing retarded loop_return_checks recurrence; no new physical wormhole solution is asserted.'])


def moving_exterior_response_checks():
    """One prescribed small displacement coupled to the linear clock recurrence."""
    delay,duration,amplitude,loops=.7,4.,.001,4
    length=math.pi+delay;final=loops*length;half=math.pi/2
    frequency=(2*math.pi/length)**2
    def motion(s):
        z=np.asarray(s,dtype=float)/duration;w=1-z;active=(z>0)&(z<1)
        raw=(256*z**4*w**4,1024*z**3*w**3*(1-2*z),
             1024*z**2*w**2*(3-14*z+14*z*z),
             6144*z*w*(1-2*z)*(1-7*z+7*z*z))
        return tuple(np.where(active,amplitude*delay*value/duration**k,0.)
                     for k,value in enumerate(raw))
    def zero_pulse(x):return 0.
    def direct_drive(x):
        result=[]
        for shift in (half+delay,half):
            jet=motion(float(x)-shift)
            result.append(float((jet[3]+frequency*jet[1])/3))
        return np.asarray(result)
    fronts=np.array([half+k*length+d+t for k in range(loops+1)
                     for d in (0.,delay) for t in (0.,duration)])
    samples=np.unique(np.r_[np.linspace(0,final,601),fronts[(fronts>=0)&(fronts<=final)]])
    launch=np.unique(np.r_[np.linspace(0,final-3*half,501),0.,.2,.5])
    ledger_times=np.unique(np.r_[np.linspace(0,final-half,401),0.,duration,delay,duration+delay])
    def record(step):
        history,channels,p=_distributed_clock_response([delay],[1.],step,final,
            zero_pulse,exterior_motion=motion)
        b0=p['b0'];q0=-frequency/6
        def source(x):
            return channels(float(x),history(float(x)))[0]
        def clocks(s):
            now=history(np.asarray(s)+half);past=history(np.asarray(s)-half)
            return tuple((now[...,i]-past[...,j])/b0 for i,j in
                (([6,7],[7,6]),([0,3],[3,0]),([1,4],[4,1]),([2,5],[5,2])))
        values=history(samples);stress=np.array([source(t) for t in samples])
        clock_jets=clocks(samples-half);rates=1+clock_jets[1]
        initial=clocks(launch)[0];arrival=clocks(launch+math.pi)[0]
        direct_return=np.column_stack((motion(launch+math.pi)[0],
                                       motion(launch+math.pi+delay)[0]))
        clock_return=initial-arrival[:,::-1]
        return_changes=clock_return+direct_return
        returns=length+return_changes
        def qout(s):
            xi=clocks(float(s))
            return source(float(s)-half)[::-1]+2*q0*xi[1]-xi[3]/3
        def qout_primitive(s):
            old=history(np.asarray(s)-half);xi=clocks(s)
            # Integral j_other = -(A_other''+A_other), since the pulse is zero.
            return -(old[..., [5,2]]+old[..., [3,0]])+2*q0*xi[0]-xi[2]/3
        xi=clocks(ledger_times)
        boundary=-b0*np.sum(xi[3]+xi[1],axis=-1)
        exterior=np.sum(qout_primitive(ledger_times)-qout_primitive(ledger_times-delay),axis=-1)
        here,old=motion(ledger_times),motion(ledger_times-delay)
        effective=(here[2]+old[2]+frequency*(here[0]+old[0]))/3
        ledger_error=float(np.max(abs(boundary+exterior-effective)))
        # Check the independently integrated outgoing stresses at off-grid times.
        quadrature_error=0.
        for s in (0.371,0.933,1.417,2.133,2.831,4.317,7.193,11.731):
            numerical=quad(lambda u:float(np.sum(qout(u))),s-delay,s,
                epsabs=2e-11,epsrel=2e-9,limit=160)[0]
            primitive=float(np.sum(qout_primitive(s)-qout_primitive(s-delay)))
            quadrature_error=max(quadrature_error,abs(numerical-primitive))
        recurrence_error=0.;decomposition_error=0.
        for x in np.unique(np.r_[np.linspace(.031,final,47),fronts[(fronts>0)&(fronts<final)]+1e-7]):
            if x>final:continue
            def q(t):
                y=history(float(t));return np.array([y[2]+frequency*y[0],y[5]+frequency*y[3]])
            expected=p['c0']*(q(x)+q(x-length)-q(x-delay)[::-1]-q(x-math.pi)[::-1])+direct_drive(x)
            recurrence_error=max(recurrence_error,float(np.max(abs(source(x)-source(x-length)-expected))))
            literal=q(x);drive=np.zeros(2)
            for k in range(int(x/length)+1):
                literal+=2*q(x-(k+1)*length)
                literal-=q(x-k*length-delay)[::-1]+q(x-k*length-math.pi)[::-1]
                drive+=direct_drive(x-k*length)
            decomposition_error=max(decomposition_error,float(np.max(abs(source(x)-p['c0']*literal-drive))))
        pre_b=np.linspace(0,half,129);pre_a=np.linspace(0,half+delay,129)
        causal_state=np.array([np.max(abs(history(pre_a)[:,0:3])),
                               np.max(abs(history(pre_b)[:,3:6]))])
        causal_source=np.array([max(abs(source(x)[0]) for x in pre_a),
                                max(abs(source(x)[1]) for x in pre_b)])
        center=values[:,1]+values[:,4]
        relative_center=3*center/(p['mean']-1)
        row=dict(maximum_step=step,solver_evaluations=p['solver_evaluations'],
            method_of_steps_segments=p['segments'],
            maximum_absolute_mean_stress_response=float(np.max(abs(stress))),
            maximum_absolute_stress_response_over_background=float(np.max(abs(stress))/abs((1-frequency)/6)),
            maximum_sampled_clock_rate_change=float(np.max(abs(rates-1))),
            minimum_sampled_clock_rate=float(rates.min()),
            maximum_absolute_center_area_excess_response_fraction=float(np.max(abs(relative_center))),
            minimum_sampled_return_changes=return_changes.min(axis=0).tolist(),
            maximum_sampled_return_changes=return_changes.max(axis=0).tolist(),
            minimum_sampled_clock_return_changes=clock_return.min(axis=0).tolist(),
            maximum_sampled_clock_return_changes=clock_return.max(axis=0).tolist(),
            maximum_sampled_direct_return_changes=direct_return.max(axis=0).tolist(),
            selected_launches=[dict(launch_time=float(launch[k]),
                clock_return_change=clock_return[k].tolist(),
                direct_geometric_return_change=direct_return[k].tolist(),
                total_return_change=return_changes[k].tolist())
                for k in (int(np.searchsorted(launch,s)) for s in (0.,.2,.5))],
            minimum_sampled_full_return_times=returns.min(axis=0).tolist(),
            maximum_sampled_full_return_times=returns.max(axis=0).tolist(),
            recurrence_residual=recurrence_error,independent_source_decomposition_error=decomposition_error,
            maximum_source_before_first_arrival=causal_source.tolist(),
            maximum_state_jets_before_first_arrival=causal_state.tolist(),
            reduced_ledger_maximum_residual=ledger_error,
            outgoing_stress_quadrature_vs_primitive_error=quadrature_error,
            maximum_reduced_effective_drive=float(np.max(effective)),
            minimum_reduced_effective_drive=float(np.min(effective)),
            maximum_reduced_ledger_after_drive=float(np.max(abs((boundary+exterior)[ledger_times>=duration+delay]))),
            plot=dict(characteristic_time=samples.tolist(),exterior_time=(samples-half).tolist(),
                displacement=motion(samples-half)[0].tolist(),
                displacement_velocity=motion(samples-half)[1].tolist(),
                direct_motion_source=np.array([direct_drive(t) for t in samples]).tolist(),
                induced_plus_source=stress[:,0].tolist(),induced_minus_source=stress[:,1].tolist(),
                clock_displacement_A=clock_jets[0][:,0].tolist(),clock_displacement_B=clock_jets[0][:,1].tolist(),
                clock_rate_A=rates[:,0].tolist(),clock_rate_B=rates[:,1].tolist(),
                center_delta_phi_over_alpha=center.tolist(),return_launch_time=launch.tolist(),
                return_change_A=return_changes[:,0].tolist(),return_change_B=return_changes[:,1].tolist(),
                clock_return_change_A=clock_return[:,0].tolist(),clock_return_change_B=clock_return[:,1].tolist(),
                full_return_A=returns[:,0].tolist(),full_return_B=returns[:,1].tolist(),
                direct_return_displacement_A=direct_return[:,0].tolist(),
                direct_return_displacement_B=direct_return[:,1].tolist(),
                ledger_time=ledger_times.tolist(),reduced_boundary_energy=boundary.tolist(),
                reduced_exterior_energy=exterior.tolist(),reduced_effective_drive=effective.tolist()))
        comparison=np.r_[values.ravel(),stress.ravel(),rates.ravel(),return_changes.ravel(),
                         boundary.ravel(),exterior.ravel()]
        return row,comparison
    records=[];comparisons=[]
    for step in (duration/40,duration/80):
        row,comparison=record(step);records.append(row);comparisons.append(comparison)
    time_error=float(np.max(abs(comparisons[0]-comparisons[1])))
    keys=('induced_plus_source','induced_minus_source','clock_rate_A','clock_rate_B',
          'clock_displacement_A','clock_displacement_B','return_change_A','return_change_B',
          'center_delta_phi_over_alpha','reduced_boundary_energy','reduced_exterior_energy')
    differences={key:float(np.max(abs(np.asarray(records[0]['plot'][key])-
        np.asarray(records[1]['plot'][key])))) for key in keys}
    # Literal preexisting one-channel source formulation, with a nonzero pulse.
    old_width=.4;old_area=math.pi*1e-4/2;old_end=2*(math.pi+.1)
    def old_pulse(x):
        return old_area*8/(3*old_width)*math.sin(math.pi*x/old_width)**4 if 0<x<old_width else 0.
    general,_,_=_distributed_clock_response([.1],[1.],old_width/24,old_end,old_pulse)
    scalar,_,_=_distributed_clock_response([.1],[1.],old_width/24,old_end,old_pulse,scalar_reference=True)
    old_times=np.linspace(0,old_end,193)
    old_error=float(np.max(abs(general(old_times)-scalar(old_times))))
    polynomial=np.polynomial.Polynomial([0,0,0,0,256,-1024,1536,-1024,256])
    zs=np.linspace(.003,.997,173)
    jet_error=max(float(np.max(abs(motion(duration*zs)[k]-
        amplitude*delay*polynomial.deriv(k)(zs)/duration**k))) for k in range(4))
    support=np.asarray(motion(np.array([-1.,0.,duration,duration+1.])))
    direct_causality=max(float(np.max(abs(direct_drive(x)))) for x in (-1.,0.,half))
    direct_a_before=max(abs(direct_drive(x)[0]) for x in np.linspace(0,half+delay,101))
    threshold=2e-8
    assert time_error<threshold and old_error<2e-9
    assert jet_error<1e-12 and np.max(abs(support))==0 and direct_causality==0 and direct_a_before==0
    assert all(row['recurrence_residual']<1e-10 and row['independent_source_decomposition_error']<1e-10
               and row['reduced_ledger_maximum_residual']<threshold
               and row['outgoing_stress_quadrature_vs_primitive_error']<threshold
               and max(row['maximum_source_before_first_arrival'])<1e-10
               and max(row['maximum_state_jets_before_first_arrival'])<1e-10
               and row['maximum_sampled_clock_rate_change']<.02
               and row['maximum_absolute_center_area_excess_response_fraction']<.02
               for row in records)
    return dict(model='Linear moving-exterior response with one straight channel and prescribed displacement',
        parameters=dict(exterior_delay=delay,period=length,omega_squared=frequency,
            b0_over_alpha=math.pi*(frequency-1)/6,displacement_fraction=amplitude,
            peak_displacement=amplitude*delay,duration=duration,loops=loops,
            final_characteristic_time=final,external_pulse=0.,
            first_source_arrival_characteristic_times_A_B=[half+delay,half],
            first_source_arrival_exterior_times_A_B=[delay,0.]),
        equations=dict(displacement='h(s)=a*D*256*z^4*(1-z)^4 for z=s/T in (0,1), zero outside; h is C3 and already includes a.',
            motion_source='i=(h third derivative+Omega squared*h first derivative)/3. For each k>=0 the A drive uses x-kL-pi/2-D and B uses x-kL-pi/2.',
            recurrence='j_A-j_A(x-L)=c0*Q[A+A_L-B_D-B_pi]+i_A, with A/B interchanged for j_B; P A=-j_A and P B=-j_B, no injected carrier.',
            clocks='xi_A(s)=[integral A(s+pi/2)-integral B(s-pi/2)]/b0; xi_B exchanges A and B. Rates are 1+xi prime.',
            returns='C_A=L+xi_A(s)-xi_B(s+pi)+h(s+pi); C_B=L+xi_B(s)-xi_A(s+pi)+h(s+pi+D), only to first order.',
            outgoing='delta qout_i(s)=j_other(s-pi/2)+2Q0*xi_i prime-xi_i third derivative/3; Q0=-Omega squared/6.',
            reduced_ledger='B=-b0*sum(xi third derivative+xi prime); E0=integral from s-D to s of sum(delta qout). B+E0=[h second derivative(s)+h second derivative(s-D)+Omega squared*(h(s)+h(s-D))]/3.'),
        refinement=records,
        verification=dict(all_declared_checks_pass=True,timestep_threshold=threshold,
            timestep_halving_maximum_absolute_difference=time_error,
            observable_timestep_differences=differences,
            legacy_nonzero_pulse_D01_vector_vs_literal_scalar_state_error=old_error,
            displacement_jets_vs_polynomial_error=jet_error,
            exact_motion_source_before_first_arrival=direct_causality,
            exact_A_motion_source_before_delayed_arrival=direct_a_before),
        scope='Retarded finite-time, first-order leading-JT/CFT calculation with one externally prescribed straight-channel displacement. Background frequencies and induced histories use unshifted times; shifting an already first-order history would add second-order terms. Clock and return extrema are sampled; times are linearized, not exact finite-amplitude transit. Direct geometric and induced-clock return changes are separately recorded; T=4 allows early nonnegative launches to sample the moving exterior at exit. The reduced ledger omits moving-domain vacuum energy (including 2Q0*h), geometric terms and higher-order production. Its effective drive is not mechanical work, a full4D energy budget or net particle production. No physical mouth matching, preparation, nonlinear stability, apparatus or human passage is demonstrated.',
        independent_validation='moving_exterior_maps.py separately tests the exact moving-endpoint null maps, covariant stress transformation and amplitude-halving limit; this helper tests their coupled linear recurrence.')


def _finite_strip_diagnostics(evaluate, clock_limits, alpha, beta):
    """Sample an exact characteristic reconstruction, reusing each source query."""
    half,cut=math.pi/2,math.atan(4.)
    # A common characteristic lattice avoids repeating history inversions at
    # every spacetime point. All resolutions cover the identical strip.
    time_steps=int(12/(cut/40)); final=time_steps*cut/40
    assert final+cut-half<min(clock_limits), 'finite strip exceeds solved histories'
    levels=[]
    for n in (40,80,160,320,640):
        count=time_steps*(n//40); step=cut/n
        tau=np.arange(count+1)*step; chi=np.arange(-n,n+1)*step
        entries=np.arange(-n,count+n+1)*step-half
        right=np.array([evaluate(float(u),0) for u in entries])
        left=np.array([evaluate(float(u),1) for u in entries])
        ii=np.arange(count+1)[:,None]; jj=np.arange(2*n+1)[None,:]
        ri,li=ii+jj,ii+2*n-jj
        tangent=np.tan(chi); background=beta*(1+chi*tangent)
        phi=background+(right[ri,0]-left[li,0])*tangent+right[ri,1]+left[li,1]
        area=1+alpha*phi
        assert np.min(area)>0 and np.all(np.isfinite(phi))
        strength=alpha*np.cos(chi)**2*(abs(right[ri,3])+abs(left[li,3]))/area
        assert np.all(np.isfinite(strength))
        peak=np.unravel_index(np.argmax(strength),strength.shape)
        row=dict(half_space_intervals=n,time_samples=count+1,space_samples=2*n+1,
            maximum_sampled_strip_strength=float(strength.max()),
            maximum_sampled_center_strength=float(strength[:,n].max()),
            maximum_sampled_abs_phi=float(alpha*np.max(abs(phi))),
            minimum_sampled_area_ratio=float(area.min()),
            peak_global_time=float(tau[peak[0]]),peak_chi=float(chi[peak[1]]))
        if levels:
            row['peak_relative_change']=abs(row['maximum_sampled_strip_strength']/levels[-1]['maximum_sampled_strip_strength']-1)
            row['area_excursion_relative_change']=abs(row['maximum_sampled_abs_phi']/levels[-1]['maximum_sampled_abs_phi']-1)
        levels.append(row)
        if n>=160 and row['peak_relative_change']<.005 and row['area_excursion_relative_change']<.005: break
    summary=dict(row)
    summary['strip_to_center_peak_ratio']=row['maximum_sampled_strip_strength']/row['maximum_sampled_center_strength']
    summary['maximum_sampled_area_change_from_initial']=float(alpha*np.max(abs(phi-background)))
    stride=n//20; time_stride=max(1,count//180)
    result=dict(alpha=alpha,beta_over_alpha=beta,chi_interval=[-cut,cut],global_time_interval=[0.,final],
        clock_limits=list(map(float,clock_limits)),maximum_required_entry_clock=final+cut-half,
        minimum_history_margin=float(min(clock_limits)-(final+cut-half)),
        summary=summary,resolution_checks=levels,
        sampling_converged=bool(row['peak_relative_change']<.005 and row['area_excursion_relative_change']<.005),
        plot=dict(global_time=tau[::time_stride].tolist(),chi=chi[::stride].tolist(),
            chiral_strength=strength[::time_stride,::stride].tolist(),phi_over_alpha=phi[::time_stride,::stride].tolist()),
        profile=dict(global_time=tau.tolist(),maximum_strength=strength.max(axis=1).tolist(),center_strength=strength[:,n].tolist(),
            minimum_area_ratio=area.min(axis=1).tolist(),maximum_abs_phi=(alpha*abs(phi).max(axis=1)).tolist()),
        definition='At one event (tau,chi), entry clocks are tau+chi-pi/2 and tau-chi-pi/2. phi=alpha*[beta*(1+chi*tan chi)+(A-B)*tan chi+A_prime+B_prime]. S=alpha*cos(chi)^2*(abs(Jplus)+abs(Jminus))/(1+phi)',
        scope='Finite sampled strip |chi|<=atan4, common tau interval, fixed local observer. Display arrays downsampled; summaries use finest evaluated grid. Sampling convergence is not a certified continuum bound. Excludes trace/fluctuations, near-mouth region, full4D corrections and physical preparation')
    if not result['sampling_converged']:
        from scipy.optimize import minimize_scalar
        envelopes=[]
        for stride in (4,2,1):
            coordinates=entries[::stride]; bounds=[]
            for side,table in enumerate((right,left)):
                maxima=[]
                for column in (0,1,3):
                    sampled=abs(table[::stride,column])
                    peaks=np.flatnonzero((sampled[1:-1]>sampled[:-2])&(sampled[1:-1]>=sampled[2:]))+1
                    maximum=float(sampled.max())
                    for k in peaks:
                        opt=minimize_scalar(lambda z:-abs(evaluate(float(z),side)[column]),
                            bounds=(coordinates[k-1],coordinates[k+1]),method='bounded',options={'xatol':1e-11})
                        assert opt.success
                        maximum=max(maximum,float(-opt.fun))
                    maxima.append(maximum)
                bounds.append(maxima)
            area_excursion=alpha*(beta*(1+cut*math.tan(cut))+
                math.tan(cut)*(bounds[0][0]+bounds[1][0])+bounds[0][1]+bounds[1][1])
            assert area_excursion<1
            envelopes.append(dict(characteristic_samples=len(coordinates),
                optimized_absolute_maxima_A_Aprime_J=bounds,
                area_excursion_envelope=area_excursion,
                chiral_strength_envelope=alpha*(bounds[0][2]+bounds[1][2])/(1-area_excursion)))
        error=abs(envelopes[-1]['chiral_strength_envelope']/envelopes[-2]['chiral_strength_envelope']-1)
        result['source_envelope']=dict(resolution_checks=envelopes,last_relative_change=error,
            numerically_converged=bool(error<1e-6),
            derivation='If |A|<=a, |B|<=b, |A_prime|<=ap, |B_prime|<=bp and |J_i|<=ji over the characteristic interval, then |phi|<=alpha*[beta*(1+cut*tan(cut))+(a+b)*tan(cut)+ap+bp]=P and S<=alpha*(jplus+jminus)/(1-P). cos(chi)^2<=1; P<1 is required',
            scope='An exact conditional inequality supplied with numerically optimized one-dimensional maxima. Successive search grids and dynamics refinements are checked, but missed extrema are not excluded by interval arithmetic. This does not make the failed direct-grid peak criterion pass')
    return result


def _spherical_residual_values(chi, right, left, alpha, beta):
    """Same-event reconstruction; columns are A/B and their first two derivatives."""
    tangent=np.tan(chi); cosine=np.cos(chi); secant=1/cosine
    difference=right[...,0]-left[...,0]
    phi=alpha*(beta*(1+chi*tangent)+difference*tangent+right[...,1]+left[...,1])
    offset=(beta*(tangent+chi*secant**2)+difference*secant**2)/2
    pu=alpha*(right[...,2]+right[...,1]*tangent+offset)
    pv=alpha*(left[...,2]-left[...,1]*tangent-offset)
    area=1+phi; assert np.all(area>0)
    residual=np.stack((2*(cosine*pu/area)**2,2*(cosine*pv/area)**2,
        -2*(3*phi**2+alpha*area)/(3*area**2),
        (phi*(3+2*phi)+cosine**2*pu*pv)/area**2),axis=-1)
    return phi,residual


def _throat_residual_envelope(evaluate, final, alpha, beta):
    """Bound the available characteristic domain using exterior-time jet extrema."""
    from functools import lru_cache
    from scipy.optimize import minimize_scalar
    evaluate=lru_cache(maxsize=None)(evaluate)
    cut=math.atan(4.); half=math.pi/2; tangent=4.; secant=math.sqrt(17.)
    limits=[float(evaluate(final,side)[0]) for side in (0,1)]
    assert min(limits)>0 and alpha>0 and beta>0
    levels=[]; previous=None
    for count in (401,801,1601):
        coordinates=np.linspace(0,final,count)
        tables=[np.array([evaluate(float(s),side) for s in coordinates]) for side in (0,1)]
        maxima=[]; locations=[]
        for side,table in enumerate(tables):
            assert np.all(np.diff(table[:,0])>0)
            side_max=[]; side_locations=[]
            for column in (1,2,3):
                sampled=abs(table[:,column]); k=int(np.argmax(sampled))
                maximum=float(sampled[k]); location=float(coordinates[k])
                peaks=np.flatnonzero((sampled[1:-1]>sampled[:-2])&(sampled[1:-1]>=sampled[2:]))+1
                for k in peaks:
                    opt=minimize_scalar(lambda s:-abs(evaluate(float(s),side)[column]),
                        bounds=(coordinates[k-1],coordinates[k+1]),method='bounded',options={'xatol':1e-10})
                    assert opt.success
                    if -opt.fun>maximum: maximum,location=float(-opt.fun),float(opt.x)
                side_max.append(maximum); side_locations.append(dict(exterior_time=location,
                    entry_clock=float(evaluate(location,side)[0])))
            maxima.append(side_max); locations.append(side_locations)
        a,b=np.array(maxima); m0=a[0]+b[0]; m1=a[1]+b[1]
        low=alpha*(beta-tangent*m0-m1)
        high=alpha*(beta*(1+cut*tangent)+tangent*m0+m1)
        area_min=1+low; assert area_min>0, 'Conditional area envelope does not establish positivity'
        common=(beta*(tangent/secant+cut*secant)+secant*m0)/2
        up=alpha*(a[2]+tangent/secant*a[1]+common)
        vp=alpha*(b[2]+tangent/secant*b[1]+common)
        ends=np.array([low,high]); areas=1+ends
        bounds=np.array([2*up**2/area_min**2,2*vp**2/area_min**2,
            2*np.max((ends/areas)**2)+2*alpha/(3*area_min),
            np.max(abs(ends*(3+2*ends)/areas**2))+up*vp/area_min**2])
        # ponytail: stream pair grids to cap memory; use certified extrema only
        # if a continuum guarantee, rather than a numerical envelope, is needed.
        negative=np.linspace(-half-cut,0,(count-1)//4+1)[:-1]
        initial=np.column_stack((negative,np.zeros((len(negative),3))))
        full=[np.vstack((initial,table)) for table in tables]
        exterior=np.r_[negative,coordinates]; peak_values=np.zeros(4); peaks=[None]*4
        phi_range=[math.inf,-math.inf]; valid_count=0
        for start in range(0,len(exterior),64):
            right=full[0][start:start+64,None,:]; left=full[1][None,:,:]
            chi=(right[...,0]-left[...,0])/2; tau=(right[...,0]+left[...,0])/2+half
            mask=(tau>=0)&(abs(chi)<=cut)
            safe_chi=np.where(mask,chi,0.)
            phi,residual=_spherical_residual_values(safe_chi,right[...,1:],left[...,1:],alpha,beta)
            valid_count+=int(mask.sum())
            if not np.any(mask): continue
            phi_range=[min(phi_range[0],float(phi[mask].min())),max(phi_range[1],float(phi[mask].max()))]
            for component in range(4):
                values=np.where(mask,abs(residual[...,component]),-np.inf)
                i,j=np.unravel_index(np.argmax(values),values.shape)
                if values[i,j]>peak_values[component]:
                    peak_values[component]=float(values[i,j])
                    peaks[component]=dict(chi=float(chi[i,j]),global_time=float(tau[i,j]),
                        right_entry_exterior=float(exterior[start+i]),left_entry_exterior=float(exterior[j]))
        assert np.all(peak_values<=bounds*(1+1e-12))
        assert low-1e-22<=phi_range[0]<=phi_range[1]<=high+1e-22
        row=dict(exterior_samples_per_side=count,optimized_jet_absolute_maxima=maxima,jet_peak_locations=locations,
            phi_interval=[float(low),float(high)],minimum_area_ratio=float(area_min),
            weighted_gradient_bounds=[float(up),float(vp)],absolute_residual_envelope=bounds.tolist(),
            maximum_orthonormal_component_envelope=float(max((bounds[0]+bounds[1]+2*bounds[2])/4,bounds[3])),
            sampled_characteristic_pairs=valid_count,sampled_absolute_residual_maxima=peak_values.tolist(),
            sampled_peak_locations=peaks,sampled_phi_interval=phi_range)
        if previous is not None:
            row['envelope_relative_change']=float(np.max(abs(bounds/np.array(previous['absolute_residual_envelope'])-1)))
            row['sampled_peak_relative_change']=float(np.max(abs(peak_values/np.array(previous['sampled_absolute_residual_maxima'])-1)))
        levels.append(row); previous=row
    chi=np.linspace(-cut,cut,161)
    return dict(alpha=alpha,beta=beta,chi_interval=[-cut,cut],entry_clock_upper_limits=limits,
        characteristic_projection_intervals=[[-half-cut,min(limits[side],limits[1-side]+2*cut)] for side in (0,1)],
        exterior_interval=[0.,final],domain='tau>=0, abs(chi)<=atan4, tau+chi-pi/2<=rmax, tau-chi-pi/2<=lmax; all negative-entry jets are zero initial data',
        domain_plot=dict(chi=chi.tolist(),maximum_global_time=(half+np.minimum(limits[0]-chi,limits[1]+chi)).tolist()),
        latest_center_global_time=half+min(limits),
        full_width_rectangle_final_time=half+min(limits)-cut,
        components=['null_plus','null_minus','mixed','angular'],resolution_checks=levels,summary=levels[-1],
        numerical_envelope_converged=bool(levels[-1]['envelope_relative_change']<1e-6),
        sampled_peaks_converged=bool(levels[-1]['sampled_peak_relative_change']<.005),
        scope='Exact conditional inequalities for the specified Polyakov reconstruction, supplied with numerically optimized jet maxima, not certified interval enclosures. Available finite interior only; no near-mouth/exterior completion, solution-error bound, global stability or backward travel. Pair-grid peak convergence and envelope convergence are distinct checks')


def nonlinear_pulse_solution(epsilon, refinement=1., loops=4, samples=None, interior=False, initial_chiral=None, lambda_time=False, exterior_delay=.1, causal_diagnostics=False, source_alpha=None, strip_diagnostics=False, null_residual=False, residual_envelope=False, timing_response=False):
    """Solve retarded loop feedback, with an optional initially exterior quantum packet."""
    from functools import lru_cache
    from scipy.optimize import brentq
    delay, width = exterior_delay, .4
    assert 0<delay<math.pi, 'This supported background requires 0<D<pi'
    assert not causal_diagnostics or lambda_time
    assert source_alpha is None or (lambda_time and source_alpha>0)
    assert not strip_diagnostics or source_alpha is not None
    assert not null_residual or source_alpha is not None
    assert not residual_envelope or source_alpha is not None
    assert not timing_response or causal_diagnostics
    length = math.pi+delay
    omega2 = (2*math.pi/length)**2
    baseline = (1-omega2)/6
    b0 = math.pi*(omega2-1)/6
    final = loops*length
    tilt = delay/math.pi if lambda_time else 0.
    if lambda_time:
        assert initial_chiral is not None and epsilon == 0 and not interior
    if initial_chiral is not None:
        assert epsilon == 0 and not interior  # Existing interior bandwidth diagnostic is for the classical carrier.
    initial_energy = 0. if initial_chiral is None else quad(initial_chiral,0,delay,epsabs=1e-12,epsrel=2e-12)[0]
    def initial_remaining(s):
        if initial_chiral is None or s>=delay: return 0.
        return quad(initial_chiral,max(0.,s),delay,epsabs=1e-12,epsrel=2e-12)[0]
    if samples is None:
        samples = np.linspace(0, final, 100*loops+1)
    segments, ends, end_clocks = [], [], [[], []]
    end_exterior = [[], []]
    minimum_lambda_history_gap = math.inf
    class RetardedHistoryUnavailable(RuntimeError):
        pass
    def history(s):
        if s <= 0:
            clock = s/(1+tilt) if lambda_time else s
            return np.array([clock, 1., 0., .5, clock, 1., 0., .5]+[0.]*(10 if lambda_time else 9))
        assert ends and s <= ends[-1]+2e-12, ('history outside solved interval', s, ends[-1:])
        index = min(int(np.searchsorted(ends, s)), len(segments)-1)
        return segments[index].sol(min(s, ends[-1]))

    def inverse(target, side):
        if target <= 0: return (1+tilt)*target if lambda_time else target
        # A state-dependent throat delay must reach COMPLETED history. Never
        # extrapolate through the active integration segment or a causal failure.
        clocks = end_clocks[side]
        if lambda_time and (not clocks or target>clocks[-1]):
            raise RetardedHistoryUnavailable(('clock query beyond completed lambda history',target,clocks[-1:]))
        assert clocks and target <= clocks[-1], ('throat history unavailable', target, clocks[-1:])
        index = min(int(np.searchsorted(clocks, target)), len(segments)-1)
        start = 0. if index == 0 else ends[index-1]
        return brentq(lambda t: segments[index].sol(t)[4*side]-target,
                      start, ends[index], xtol=5e-14)

    def inverse_exterior(target, side):
        if not lambda_time: return target
        if target<=0: return (1+tilt)*target
        clocks=end_exterior[side]
        if not clocks or target>clocks[-1]:
            raise RetardedHistoryUnavailable(('exterior query beyond completed lambda history',target,clocks[-1:]))
        index=min(int(np.searchsorted(clocks,target)),len(segments)-1)
        start=0. if index==0 else ends[index-1]
        return brentq(lambda lam:lam-tilt*segments[index].sol(lam)[4*side]-target,
                      start,ends[index],xtol=5e-14)

    def power(s):
        return math.pi*epsilon/2*8/(3*width)*math.sin(math.pi*s/width)**4 if 0<s<width else 0.

    @lru_cache(maxsize=8192)
    def quantum_history(s, side):
        return incoming_quantum(s, history(s), side)

    def incoming_quantum(s, y, side):
        if s <= 0: return baseline
        i, other = 4*side, 4*(1-side)
        schwarzian = y[i+3]-.5*y[i+1]**2
        exterior_s=s-tilt*y[i] if lambda_time else s
        if initial_chiral is not None and exterior_s<delay:
            # Initial exterior data are functions of the exterior arrival clock,
            # not of the evolving throat clock. Later returns use this history.
            seed=initial_chiral(exterior_s) if side == 0 else 0.
            return 1/6+(baseline-1/6+seed+schwarzian/3)/y[i+1]**2
        old = history(inverse_exterior(exterior_s-delay,1-side))
        entry = inverse(old[other]-math.pi, side)
        past = quantum_history(entry, side)
        old_schwarzian = old[other+3]-.5*old[other+1]**2
        return (1/6+(old[other+1]/y[i+1])**2*(past-1/6)
                +(schwarzian-old_schwarzian)/(3*y[i+1]**2))

    def fluxes(s, y):
        nonlocal minimum_lambda_history_gap
        assert y[1] > 0 and y[5] > 0
        # Entry times of light which is leaving each mouth now.
        from_left = inverse(y[0]-math.pi, 1)
        from_right = inverse(y[4]-math.pi, 0)
        left_past, right_past = history(from_left), history(from_right)
        if lambda_time:
            minimum_lambda_history_gap=min(minimum_lambda_history_gap,s-from_left,s-from_right)
        carrier_in = power(s-tilt*y[0] if lambda_time else s)/y[1]**2
        carrier_out = power(from_right-tilt*right_past[0] if lambda_time else from_right)/right_past[1]**2
        incoming = np.array([incoming_quantum(s, y, 0)+carrier_in,
                             incoming_quantum(s, y, 1)])
        quantum_out = np.array([quantum_history(from_left, 1), quantum_history(from_right, 0)])
        outgoing = quantum_out+np.array([0., carrier_out])
        q_out = np.array([y[1]**2*(quantum_out[0]-1/6)-(y[3]-.5*y[1]**2)/3,
                          y[5]**2*(quantum_out[1]-1/6)-(y[7]-.5*y[5]**2)/3])
        return incoming, outgoing, q_out, y[5]**2*carrier_out, from_left, from_right

    def rhs(s, y):
        incoming, outgoing, q_out, absorbed, _, _ = fluxes(s, y)
        dy = np.zeros(len(y))
        for side in (0, 1):
            i, bulk = 4*side, 8+3*side
            v, a, k = y[i+1:i+4]
            dy[i:i+4] = (v, a, v*(k-.5*v*v)+1.5*a*a/v,
                          v*v*(outgoing[side]-incoming[side])/b0)
            dy[bulk:bulk+3] = v*np.array([y[bulk+1], y[bulk+2],
                                          -y[bulk+1]-(incoming[side]-baseline)])
            if lambda_time:
                dy[i:i+4]/=1+tilt*v
                dy[bulk:bulk+3]/=1+tilt*v
        if lambda_time:
            h=1/(1+tilt*y[[1,5]]); q=(q_out-(baseline-1/6))*h
            dy[14:]=(q.sum(),h[0]*power(s-tilt*y[0]),h[1]*absorbed,q[0])
        else:
            dy[14:] = (q_out.sum()-2*(baseline-1/6), power(s), absorbed)
        return dy

    start, initial = 0., history(0.)
    while (min(start-tilt*initial[[0,4]])<final if lambda_time else start<final):
        step=delay/2
        if lambda_time:
            assert len(segments)<10000, 'lambda segment budget exhausted; not a physical obstruction'
            ancestors=[inverse(initial[0]-math.pi,1),inverse(initial[4]-math.pi,0),
                       inverse_exterior(start-tilt*initial[0]-delay,1),
                       inverse_exterior(start-tilt*initial[4]-delay,0)]
            gap=start-max(ancestors)
            assert gap>0, 'lambda foliation no longer strictly retarded'
            minimum_lambda_history_gap=min(minimum_lambda_history_gap,gap)
            step=min(step,.4*gap)
        while True:
            stop=start+step if lambda_time else min(start+step,final)
            try:
                segment = solve_ivp(rhs, (start, stop), initial, dense_output=True,
                                    method='DOP853', rtol=3e-11*refinement**2,
                                    atol=3e-13*refinement**2,
                                    max_step=(width/20 if initial_chiral is None else delay/20)*refinement)
                break
            except RetardedHistoryUnavailable:
                # Shorten only a future-history step; never extrapolate a source
                # or reinterpret this numerical resolution limit as protection.
                step/=2
                if step<1e-10: raise
        assert segment.success
        segments.append(segment); ends.append(stop)
        initial, start = segment.y[:, -1], stop
        for side in (0, 1):
            assert not end_clocks[side] or initial[4*side] > end_clocks[side][-1]
            end_clocks[side].append(initial[4*side])
            if lambda_time:
                physical=start-tilt*initial[4*side]
                assert not end_exterior[side] or physical>end_exterior[side][-1]
                end_exterior[side].append(physical)
    if lambda_time:
        def at_exterior(t,side):
            lam=inverse_exterior(float(t),side)
            return lam,history(lam)
        def physical_entry(lam,side):
            return lam-tilt*history(lam)[4*side]
        def right_primitive(t): return at_exterior(t,0)[1][17]
        def left_primitive(t):
            state=at_exterior(t,1)[1]
            return state[14]-state[17]
        # Primary charge test uses the actual spacelike equal-lambda slices.
        lambda_grid=np.linspace(0,ends[-1],161); lambda_ledgers=[]
        for lam in lambda_grid:
            y=history(float(lam)); sr,sl=lam-tilt*y[[0,4]]
            outside=y[14]-right_primitive(sl-delay)-left_primitive(sr-delay)+initial_remaining(sr)
            lambda_ledgers.append(abs(-b0*(y[3]+y[7]-1)+outside-y[15]+y[16]-initial_energy))
        states=[]; matching=[]; ledgers=[]; transit=[]; returns=[]
        center=[]; center_times=[]; sources=[]; gate_rows=[]; center_strength=[]; center_sources=[]; center_areas=[]
        center_gradients=[]
        center_phi=[]
        for t in samples:
            lr,yr=at_exterior(t,0); ll,yl=at_exterior(t,1)
            fr,fl=fluxes(lr,yr),fluxes(ll,yl)
            left_entry,right_entry=fr[4],fl[5]
            predicted=np.array([1+(yr[8]-history(left_entry)[11])/b0,
                                1+(-history(right_entry)[8]+yl[11])/b0])
            matching.append(float(np.max(abs(np.array([yr[1],yl[5]])-predicted))))
            outside=right_primitive(t)-right_primitive(t-delay)+left_primitive(t)-left_primitive(t-delay)+initial_remaining(t)
            ledgers.append(abs(-b0*(yr[3]+yl[7]-1)+outside-yr[15]+yl[16]-initial_energy))
            transit.extend((t-physical_entry(left_entry,1),t-physical_entry(right_entry,0)))
            gates=[]
            for side in (0,1):
                old=at_exterior(t-delay,1-side)[1]
                entry=inverse(old[4*(1-side)]-math.pi,side)
                returns.append(t-physical_entry(entry,side))
                if causal_diagnostics:
                    gates.append(float(math.pi+(yr[0] if side==0 else yl[4])-old[4*(1-side)]))
            if causal_diagnostics: gate_rows.append(gates)
            target=min(yr[0],yl[4]); cr=inverse(float(target),0); cl=inverse(float(target),1)
            aa=history(cr); bb=history(cl)
            center.append((aa[9]+bb[12])*3/(omega2-1)); center_times.append(float(target+math.pi/2))
            if source_alpha is not None:
                # Both stresses refer to this SAME center event, not to the
                # different events represented by the two mouth-source columns.
                jp,jm=quantum_history(cr,0),quantum_history(cl,1)
                phi=source_alpha*((omega2-1)/3+aa[9]+bb[12]); area=1+phi
                assert area>0, ('center area expansion invalid',area)
                center_sources.append([float(jp),float(jm)]); center_areas.append(float(area))
                center_strength.append(float(source_alpha*(abs(jp)+abs(jm))/area))
                if null_residual:
                    center_phi.append(float(phi))
                    offset=(aa[8]-bb[11])/2
                    center_gradients.append([float(source_alpha*(aa[10]+offset)),
                                             float(source_alpha*(bb[13]-offset))])
            sources.append([fr[0][0],fl[0][1]])
            states.append(np.r_[yr[:4],yl[4:8],yr[8:11],yl[11:14],yr[17]+yl[14]-yl[17],yr[15],yl[16]])
        values=np.array(states)
        assert max(matching)<3e-8 and max(ledgers)<3e-9 and max(lambda_ledgers)<3e-9, (
            'lambda matching/common-time/slice charge residuals',max(matching),max(ledgers),max(lambda_ledgers))
        assert min(returns)>0 and np.min(values[:,[1,5]])>0 and minimum_lambda_history_gap>0
        record=dict(epsilon=epsilon,pulse_integral_over_alpha=math.pi*epsilon/2,
            normalized_launch_integral=float(values[-1,15]),normalized_absorption_integral=float(values[-1,16]),
            final_boundary_energy_change_over_alpha=-b0*float(values[-1,3]+values[-1,7]-1),
            final_exterior_energy_change_over_alpha=float(outside-initial_energy),
            maximum_bulk_boundary_rate_mismatch=max(matching),maximum_reduced_energy_ledger_error=float(max(ledgers)),
            minimum_sampled_one_way_time=float(min(transit)),minimum_sampled_complete_return_time=float(min(returns)),
            minimum_sampled_clock_rate=float(np.min(values[:,[1,5]])),
            maximum_sampled_clock_rate_change=float(np.max(abs(values[:,[1,5]]-1))),
            maximum_sampled_center_area_excess_change_fraction=float(np.max(abs(np.array(center)))),
            final_clock_offset=float(values[-1,0]-values[-1,4]),
            initial_quantum_packet_energy_over_alpha=initial_energy,
            maximum_sampled_clock_offset=float(np.max(abs(values[:,0]-values[:,4]))),
            plot=dict(exterior_time=samples.tolist(),right_clock_rate=values[:,1].tolist(),left_clock_rate=values[:,5].tolist(),
                      right_clock_shift=(values[:,0]-samples).tolist(),left_clock_shift=(values[:,4]-samples).tolist(),
                      center_global_time=center_times,center_area_excess_change_fraction=center,incoming_total_j_over_alpha=sources),
            lambda_foliation=dict(k=tilt,final_lambda=ends[-1],final_boundary_exterior_times=[row[-1] for row in end_exterior],
                                 maximum_slice_energy_ledger_error=float(max(lambda_ledgers)),
                                 minimum_evaluated_retarded_lambda_gap=float(minimum_lambda_history_gap),segments=len(segments),
                                 scope='Physical derivative variables retained; common-exterior outputs are reconstructed separately. No continuation through a closed-return zero or completed4D exterior'),
            launched_probes=[])
        def launched_probe(launch,side):
            target=at_exterior(launch,side)[1][4*side]+math.pi
            reached=bool(target<=at_exterior(final,1-side)[1][4*(1-side)])
            arrival=physical_entry(inverse(target,1-side),1-side) if reached else None
            return dict(launch_exterior_time=launch,launch_side=side,
                        exit_reached_in_interval=reached,exit_exterior_time=arrival,
                        complete_return_delay=None if arrival is None else arrival+delay-launch)
        record['launched_probes']=[launched_probe(launch,side) for launch in (0.,delay,1.) for side in (0,1)]
        if causal_diagnostics:
            gates=np.array(gate_rows)
            matched=[]
            for phase in np.linspace(0,loops-1,8*(loops-1)+1):
                for side in (0,1):
                    probe=launched_probe(float(phase*length),side)
                    probe.update(launch_phase=float(phase),complete_return_within_interval=bool(
                        probe['exit_reached_in_interval'] and probe['exit_exterior_time']+delay<=final))
                    matched.append(probe)
            record['causal_diagnostics']=dict(exterior_delay=delay,initial_loop_time=length,
                normalized_exterior_time=(samples/length).tolist(),normalized_complete_clock_gates=(gates/length).tolist(),
                minimum_sampled_normalized_complete_clock_gate=float(gates.min()/length),
                minimum_sampled_one_way_clock_gate=float(math.pi-np.max(abs(values[:,0]-values[:,4]))),
                phase_matched_probes=matched,
                scope='Clock-gate magnitude can shrink merely because the receiving clock slows. Compare complete physical return/L for identical launch phases and sides; missing exits are censored, not zero or permanent trapping. Legacy arrival-indexed minima include virtual initial ancestors')
            if timing_response:
                response=[]
                for probe in matched:
                    side=probe['launch_side']; other=1-side; launch=probe['launch_exterior_time']
                    target=at_exterior(launch,side)[1][4*side]+math.pi
                    end_clock=at_exterior(final,other)[1][4*other]
                    gate=target-at_exterior(launch-delay,other)[1][4*other]
                    row=dict(launch_side=side,launch_phase=probe['launch_phase'],matching_target=float(target),
                        endpoint_exit_margin=float(target-end_clock),complete_return_phase_margin=float(gate),
                        return_margin_uses_initial_history=bool(launch<delay),
                        exit_reached_in_interval=probe['exit_reached_in_interval'])
                    assert gate>0
                    if probe['exit_reached_in_interval']:
                        arrival=probe['exit_exterior_time']; rate=at_exterior(arrival,other)[1][4*other+1]
                        assert rate>0 and abs(at_exterior(arrival,other)[1][4*other]-target)<2e-12
                        row.update(exit_exterior_time=arrival,receiving_clock_rate=float(rate),
                            linear_time_gain_per_phase=float(1/rate),complete_return_delay=probe['complete_return_delay'])
                        step=.01*length
                        assert 0<arrival-step<arrival+step<final
                        row['phase_to_advance_one_percent_loop']=float(target-at_exterior(arrival-step,other)[1][4*other])
                        row['phase_to_delay_one_percent_loop']=float(at_exterior(arrival+step,other)[1][4*other]-target)
                        intervals=[]
                        for error in (1e-9,1e-6,1e-3):
                            endpoints=[]
                            for sign in (-1,1):
                                shifted=target+sign*error
                                exit_time=physical_entry(inverse(shifted,other),other) if shifted<=end_clock else None
                                endpoints.append(exit_time)
                            intervals.append(dict(assumed_combined_phase_error=error,
                                exit_interval=endpoints,return_delay_interval=[None if t is None else t+delay-launch for t in endpoints],
                                exit_shift_interval=[None if t is None else t-arrival for t in endpoints],
                                linear_shift_magnitude=float(error/rate)))
                        row['conditional_response_intervals']=intervals
                    else:
                        assert row['endpoint_exit_margin']>0
                    response.append(row)
                record['timing_response']=dict(probes=response,
                    minimum_censored_exit_phase_margin=min(p['endpoint_exit_margin'] for p in response if not p['exit_reached_in_interval']),
                    minimum_launched_complete_return_phase_margin=min(p['complete_return_phase_margin'] for p in response),
                    minimum_physical_return_phase_margin=min(p['complete_return_phase_margin'] for p in response if not p['return_margin_uses_initial_history']),
                    definition='For a fixed exterior return delay, l(e)=r(s)+pi. A uniformly bounded combined matching-phase error epsilon brackets e between inverse_l(target-epsilon) and inverse_l(target+epsilon), when both are available. Local de/dphase=1/l_prime(e)',
                    scope='Kinematic conditional sensitivity using actual completed exits and a fixed exterior delay, not a predicted phase correction. Einstein residual size is not a phase-error bound. A censored exit has an endpoint margin, not an invented exit rate. Early gates may reference virtual initial history; for launch<D, future exits alone imply C>=D-launch>0. Bounds for selected launch events are not a continuous-launch chronology theorem')
        if source_alpha is not None:
            beta=(omega2-1)/3
            record['source_diagnostics']=dict(alpha=source_alpha,
                baseline_center_strength=source_alpha*beta/(1+source_alpha*beta),
                center_global_time=center_times,center_chiral_sources_over_alpha=center_sources,
                center_area_ratio=center_areas,center_chiral_strength=center_strength,
                definition='alpha*(abs(Jplus_center)+abs(Jminus_center))/(1+phi_center). Same center event; equals8pi*G*re^2*max(abs(rho_chiral),abs(flux_chiral)) for the static throat-center observer',
                scope=WEAK_SOURCE_SCOPE)
            if null_residual:
                gradients=np.array(center_gradients); area=np.array(center_areas)[:,None]
                absolute=2*gradients**2/area**2
                source=4*source_alpha*abs(np.array(center_sources))/area
                record['null_residual']=dict(center_phi=center_phi,center_phi_u_v=gradients.tolist(),
                    absolute_observer_residual=absolute.tolist(),
                    background_plus_source_normalized_residual=(absolute/(1+source)).tolist(),
                    definition='At the same center event, F=1+phi, E_uu=phi_u^2/(2F^2), E_vv=phi_v^2/(2F^2). Absolute observer residual re^2*|E(k,k)|=2*phi_u/v^2/F^2; normalization divides by1+4*alpha*|Jplus/minus|/F',
                    scope='Exact radial-null residual of the reconstructed fixed2D-metric spherical ansatz and specified chiral source. Two center components only; not a full4D solution or stability proof. Background-plus-source normalization is not fractional source error')
        if strip_diagnostics:
            def evaluate(entry,side):
                lam=inverse(entry,side); bulk=8+3*side
                return [*history(lam)[bulk:bulk+3],quantum_history(lam,side)]
            record['strip_diagnostics']=_finite_strip_diagnostics(evaluate,
                [at_exterior(final,side)[1][4*side] for side in (0,1)],source_alpha,(omega2-1)/3)
        if residual_envelope:
            def exterior_jets(s,side):
                state=at_exterior(s,side)[1]; bulk=8+3*side
                return tuple(float(x) for x in (state[4*side],*state[bulk:bulk+3]))
            record['residual_envelope']=_throat_residual_envelope(exterior_jets,final,source_alpha,(omega2-1)/3)
        return record,values
    values = np.array([history(float(t)) for t in samples])
    matching_errors, ledger_errors, transit_gaps, return_gaps = [], [], [], []
    center, center_times, source_rows = [], [], []
    for t, y in zip(samples, values):
        incoming, outgoing, q_out, absorbed, left_entry, right_entry = fluxes(float(t), y)
        left_old, right_old = history(left_entry), history(right_entry)
        predicted = np.array([1+(y[8]-left_old[11])/b0,
                              1+(-right_old[8]+y[11])/b0])
        matching_errors.append(float(np.max(np.abs(y[[1, 5]]-predicted))))
        # Reduced boundary charges + energy currently in both exterior legs.
        outside = y[14]-history(float(t)-delay)[14]+initial_remaining(float(t))
        ledger = -b0*(y[3]+y[7]-1)+outside-y[15]+y[16]-initial_energy
        ledger_errors.append(abs(float(ledger)))
        transit_gaps.extend((float(t-left_entry), float(t-right_entry)))
        old = history(float(t)-delay)
        return_gaps.extend((float(t-inverse(old[4]-math.pi, 0)),
                            float(t-inverse(old[0]-math.pi, 1))))
        # Select a common global center time covered by BOTH incoming histories.
        target = min(y[0], y[4], end_clocks[0][-1], end_clocks[1][-1])
        aa = history(inverse(float(target), 0))
        bb = history(inverse(float(target), 1))
        center.append((aa[9]+bb[12])*3/(omega2-1))
        center_times.append(float(target+math.pi/2))
        source_rows.append(incoming.tolist())
    assert max(matching_errors)<3e-8 and max(ledger_errors)<3e-9
    assert min(transit_gaps)>delay/2 and min(return_gaps)>0
    assert np.min(values[:, [1, 5]])>0
    record = dict(epsilon=epsilon, pulse_integral_over_alpha=math.pi*epsilon/2,
                normalized_launch_integral=float(values[-1, 15]),
                normalized_absorption_integral=float(values[-1, 16]),
                final_boundary_energy_change_over_alpha=-b0*float(values[-1, 3]+values[-1, 7]-1),
                final_exterior_energy_change_over_alpha=float(values[-1, 14]-history(final-delay)[14]+initial_remaining(final)-initial_energy),
                maximum_bulk_boundary_rate_mismatch=max(matching_errors),
                maximum_reduced_energy_ledger_error=max(ledger_errors),
                minimum_sampled_one_way_time=min(transit_gaps),
                minimum_sampled_complete_return_time=min(return_gaps),
                minimum_sampled_clock_rate=float(np.min(values[:, [1, 5]])),
                maximum_sampled_clock_rate_change=float(np.max(np.abs(values[:, [1, 5]]-1))),
                maximum_sampled_center_area_excess_change_fraction=float(np.max(np.abs(center))),
                final_clock_offset=float(values[-1, 0]-values[-1, 4]),
                plot=dict(exterior_time=samples.tolist(), right_clock_rate=values[:, 1].tolist(),
                          left_clock_rate=values[:, 5].tolist(), right_clock_shift=(values[:, 0]-samples).tolist(),
                          left_clock_shift=(values[:, 4]-samples).tolist(),
                          center_global_time=center_times, center_area_excess_change_fraction=center,
                          incoming_total_j_over_alpha=source_rows))
    if initial_chiral is not None:
        record['initial_quantum_packet_energy_over_alpha']=initial_energy
        record['initial_data_scope']='Unshifted reduced throat and fixed-D exterior. The packet is supplied on the initial exterior segment; negative-time history represents initial characteristic data, not physical preparation or a static past containing that pulse. Quantum state recirculates, with no classical launch or absorber'
        probes=[]
        for launch in (0.,delay,1.):
            for side in (0,1):
                target=history(launch)[4*side]+math.pi
                reached=bool(target<=end_clocks[1-side][-1])
                arrival=inverse(target,1-side) if reached else None
                probes.append(dict(launch_exterior_time=launch,launch_side=side,
                                   exit_reached_in_interval=reached,exit_exterior_time=arrival,
                                   complete_return_delay=None if arrival is None else arrival+delay-launch))
        record['launched_probes']=probes
        record['maximum_sampled_clock_offset']=float(np.max(abs(values[:,0]-values[:,4])))
    if interior:
        alpha, beta_hat, half = 1e-9, (omega2-1)/3, math.pi/2
        def field(tau, chi, tangent=None, secant_squared=None):
            ra = history(inverse(tau+chi-half, 0))[8:11]
            lb = history(inverse(tau-chi-half, 1))[11:14]
            tangent = math.tan(chi) if tangent is None else tangent
            secant_squared = 1/math.cos(chi)**2 if secant_squared is None else secant_squared
            difference = ra[0]-lb[0]
            background = beta_hat*(1+chi*tangent)
            delta = difference*tangent+ra[1]+lb[1]
            g = beta_hat*(tangent+chi*secant_squared)+difference*secant_squared
            du = ra[2]+ra[1]*tangent+g/2
            dv = lb[2]-lb[1]*tangent-g/2
            return background+delta, du, dv, background, delta

        def sources(tau, chi):
            right, left = inverse(tau+chi-half, 0), inverse(tau-chi-half, 1)
            return np.array([quantum_history(right, 0)+power(right)/history(right)[1]**2,
                             quantum_history(left, 1)])

        times = np.unique(np.r_[np.linspace(0, final, 241), width/4, width/2, width])
        chi_grid = np.linspace(-.8, .8, 41)
        rows, min_full_area, max_perturbation, max_abs_phi = [], math.inf, 0., 0.
        trapped, anti_trapped, stationary_points = 0, 0, []
        local_source_strength = 0.
        for s in times:
            y = history(float(s))
            tau = min(y[0], y[4])  # This complete strip is covered by both histories.
            data = np.array([field(tau, chi) for chi in chi_grid])
            phi, du, dv, old, change = data.T
            min_full_area = min(min_full_area, float(np.min(1+alpha*phi)))
            max_abs_phi = max(max_abs_phi, float(np.max(np.abs(alpha*phi))))
            max_perturbation = max(max_perturbation, float(np.max(np.abs(alpha*change/(1+alpha*old)))))
            trapped += int(np.sum((du < -1e-8)&(dv < -1e-8)))
            anti_trapped += int(np.sum((du > 1e-8)&(dv > 1e-8)))
            def slope(chi):
                _, up, vp, _, _ = field(tau, chi)
                return up-vp
            if slope(-.8)<0<slope(.8):
                root = brentq(slope, -.8, .8, xtol=5e-14)
                value, up, vp, _, _ = field(tau, root)
                second = math.tan(root)*(up-vp)+value/math.cos(root)**2-sources(tau, root).sum()
                stationary_points.append(dict(exterior_time=float(s), global_time=float(tau), chi=root,
                    phi_over_alpha=value, time_derivative_over_alpha=up+vp,
                    spatial_second_derivative_over_alpha=float(second)))
            rows.append(dict(exterior_time=float(s), global_time=float(tau),
                             minimum_phi_over_alpha=float(np.min(phi)),
                             maximum_phi_over_alpha=float(np.max(phi))))
        # ponytail: a declared finite grid diagnoses source magnitude; continuum
        # bounds and stress fluctuations require a separate validity analysis.
        for s in times[::4]:
            y = history(float(s)); tau = min(y[0], y[4])
            for chi in chi_grid[::5]:
                phi = field(tau, chi)[0]
                jp, jm = sources(tau, chi)
                strength = alpha*math.cos(chi)**2*max(abs(jp+jm), abs(jp-jm))/(1+alpha*phi)
                local_source_strength = max(local_source_strength, float(strength))
        def bandwidth_integrands(s):
            v, a = history(s)[1:3]
            shape = math.sin(math.pi*s/width)**2
            derivative = math.pi/width*math.sin(2*math.pi*s/width)
            return (derivative*v-shape*a)**2/v**5, shape**2/v
        numerator = quad(lambda s: bandwidth_integrands(s)[0], 0, width, epsabs=1e-10)[0]
        denominator = quad(lambda s: bandwidth_integrands(s)[1], 0, width, epsabs=1e-12)[0]
        bandwidth = math.sqrt(numerator/denominator)
        if epsilon == 0:
            assert abs(bandwidth-2*math.pi/(math.sqrt(3)*width)) < 1e-10
        assert min_full_area > 0

        # This tests the finite-cutoff correction to the leading clock prescription,
        # not the omitted four-dimensional metric/exterior corrections.
        cutoffs = []
        for phi_cut in (1e-5, 1e-4):
            target = phi_cut/alpha
            delta0 = brentq(lambda d: d*beta_hat*(1+(half-d)/math.tan(d))-d*target,
                             1e-9, .01, xtol=5e-20)
            lapse_errors, speeds, deltas = [], [], []
            for s in times[::4]:
                y = history(float(s))
                for side, sign in ((0, 1), (1, -1)):
                    tau, rate = y[4*side], y[4*side+1]
                    def cutoff_field(delta):
                        return field(tau, sign*(half-delta), sign/math.tan(delta), 1/math.sin(delta)**2)
                    delta = brentq(lambda d: d*(cutoff_field(d)[0]-target), 1e-9, .01, xtol=5e-20)
                    _, up, vp, _, _ = cutoff_field(delta)
                    velocity = -(up+vp)/(up-vp)
                    assert abs(velocity) < 1
                    lapse = rate*math.sin(delta0)/math.sin(delta)*math.sqrt(1-velocity**2)
                    lapse_errors.append(abs(lapse-1)); speeds.append(abs(velocity)); deltas.append(delta)
            cutoffs.append(dict(phi_cut=phi_cut, initial_angular_distance=delta0,
                                minimum_angular_distance=min(deltas), maximum_angular_distance=max(deltas),
                                maximum_abs_dchi_dtau=max(speeds), maximum_relative_lapse_mismatch=max(lapse_errors)))
        record['interior'] = dict(
            alpha=alpha, beta=alpha*beta_hat, chi_interval=[-.8, .8], time_samples=len(times), chi_samples=41,
            minimum_sampled_full_area_over_extremal_area=min_full_area,
            maximum_sampled_abs_phi=max_abs_phi,
            maximum_sampled_fractional_total_area_change=max_perturbation,
            future_trapped_grid_points=trapped, past_trapped_grid_points=anti_trapped,
            sign_test='Both future null derivatives Phi_u,Phi_v below -1e-8 (above +1e-8) give future (past) trapped spheres in the reconstructed leading area. A single zero is not an event horizon; static center has both zero',
            spatial_stationary_points=stationary_points, slice_summary=rows,
            maximum_sampled_chiral_source_strength=local_source_strength,
            source_strength_definition='max(|Jplus+Jminus|,|Jplus-Jminus|)*alpha*cos(chi)^2/(1+phi) = max absolute modeled chiral energy/flux in a local orthonormal frame, multiplied by 8*pi*G*re^2. Excludes trace sector, metric corrections and fluctuations',
            source_grid=dict(time_samples=len(times[::4]), chi_samples=len(chi_grid[::5])),
            carrier_rms_u_frequency=bandwidth,
            carrier_bandwidth_scope='RMS frequency of the transported classical scalar derivative a(s)/r_prime(s), a=sin^2(pi*s/w). Amplitude-independent shape diagnostic, including at epsilon0. This is not a fermionic Landau-level occupation or a quantum-noise bound',
            finite_cutoff_checks=cutoffs,
            scope='Finite interior grid and implicit finite-cutoff curves in the reconstructed leading throat. Small area correction and timelike cutoffs do not prove stability, full4D matching or passage of a physical object. Trapped/anti-trapped spheres under NEC-violating matter need not form a permanent event horizon')
    return record, values


def _polyakov_static_solution(alpha, frequency2, endpoint, refinement=1):
    """Normalized static IVP shared by the interior and its matching diagnostics."""
    def rhs(x, y):
        Y, V, Z, W = y
        S = 1+alpha*Y
        sec2 = 1/math.cos(x)**2
        e2s = sec2*math.exp(2*alpha*Z)
        A = (e2s*(Y/S-1/(3*S*S))-alpha*V*V/(4*S))/(S-alpha/3)
        return (V, 2*S*A+alpha*V*V/(2*S), W,
                sec2*math.expm1(2*alpha*Z-2*math.log1p(alpha*Y))/alpha-A)
    initial = (frequency2/(3*(1-alpha*frequency2/3)), 0., 0., 0.)
    sol = solve_ivp(rhs, (0,endpoint), initial, method='DOP853', dense_output=True,
                    rtol=2e-12/refinement**2, atol=2e-13/refinement**2,
                    max_step=.015/refinement)
    assert sol.success
    return sol, rhs


def polyakov_junction_checks(static):
    """Israel matching requirements for abrupt joins to RN electrovac, not a shell design."""
    alpha, frequency2 = static['alpha'], static['Omega_squared']
    endpoint = static['symmetric_x_interval'][1]
    sol, _ = _polyakov_static_solution(alpha, frequency2, endpoint, refinement=2)
    w,v,R,A,K = sp.symbols('w v R A K', positive=True)
    sigma_surface = 2*(v-w)/R  # 8piG re times physical surface energy density.
    pressure_surface = (A/w+w-2*v)/(2*R)-K
    dp = sp.simplify(sp.diff(pressure_surface,w))
    assert sp.simplify(dp-(1-A/w**2)/(2*R)) == 0
    rn_radius,rn_mass = sp.symbols('RN_radius RN_mass', positive=True)
    rn_F = 1-2*rn_mass/rn_radius+1/rn_radius**2
    rn_R2 = sp.simplify(-sp.diff(rn_F,rn_radius,2))
    assert sp.simplify(rn_R2-(4*rn_mass/rn_radius**3-6/rn_radius**4)) == 0
    rows = []
    for x in (.5,1.,endpoint):
        Y,V,Z,W = sol.sol(x)
        S = 1+alpha*Y
        radius = math.sqrt(S)
        lapse = math.exp(alpha*Z)/math.cos(x)
        sx = math.tan(x)+alpha*W
        normal_r = alpha*V/(2*radius*lapse)
        inward_acceleration = sx/lapse
        radius_excess = alpha*Y/(radius+1)
        Avalue = alpha*Y/S  # 1-1/R², avoiding subtraction.
        # M-1 is O(alpha²); computing M first would erase this information.
        mass_excess = alpha**2/(2*radius)*(Y*Y/(radius+1)**2-V*V/(4*lapse*lapse))
        outward_acceleration = (radius_excess+radius*mass_excess)/(radius**3*normal_r)
        pressure_jump = outward_acceleration-inward_acceleration
        bulk_pressure = -alpha/(3*S*lapse*lapse)*(sx*sx+frequency2)
        independent_pressure = -radius*bulk_pressure/(2*normal_r)
        force_residual = abs(pressure_jump-independent_pressure)
        energy_upper = 2*normal_r/radius
        assert normal_r>0 and bulk_pressure<0 and Avalue>normal_r**2
        assert pressure_jump>energy_upper and force_residual<2e-10
        # Exterior branch is checked only on R>=Rcut; its formal horizons
        # belong to the discarded RN continuation, not the joined interior.
        horizon_excess = (mass_excess+math.sqrt(mass_excess*(2+mass_excess))) if mass_excess>=0 else None
        if horizon_excess is not None:
            assert radius_excess>horizon_excess
        eta = np.linspace(.05,1.,120)
        exterior_speed = eta*normal_r
        surface_energy = 2*(normal_r-exterior_speed)/radius
        surface_pressure = (Avalue/exterior_speed+exterior_speed-2*normal_r)/(2*radius)-inward_acceleration
        assert np.all(surface_pressure>surface_energy)
        rows.append(dict(x=x, sphere_radius_excess_over_re=radius_excess,
            outward_normal_radius_derivative=normal_r, interior_lapse_relative_to_center=lapse,
            interior_proper_log_lapse_derivative_times_re=inward_acceleration,
            quantum_normal_pressure_times_8piG_re2=bulk_pressure,
            derivative_matched_RN_mass_excess_over_re=mass_excess,
            mass_excess_divided_by_alpha_squared=mass_excess/alpha**2,
            derivative_matched_exterior_F=normal_r**2,
            derivative_matched_RN_horizon_excess_over_re=horizon_excess,
            derivative_matched_surface_energy_times_8piG_re=0.,
            derivative_matched_surface_pressure_times_8piG_re=pressure_jump,
            pressure_from_normal_force_identity=independent_pressure,
            normal_force_identity_residual=force_residual,
            all_nonnegative_energy_shells_pressure_lower_bound=pressure_jump,
            all_nonnegative_energy_shells_energy_upper_bound=energy_upper,
            dominant_energy_condition_gap_lower_bound=pressure_jump-energy_upper,
            exterior_to_interior_time_scale=lapse/normal_r,
            plot=dict(sqrt_F_over_interior_normal_radius_derivative=eta.tolist(),
                      surface_energy_times_8piG_re=surface_energy.tolist(),
                      surface_pressure_times_8piG_re=surface_pressure.tolist())))
    # Refinement is independent of the analytic sign argument and remains a
    # numerical check of the inherited interior data used in the tables.
    coarse, _ = _polyakov_static_solution(alpha,frequency2,endpoint)
    sampled = np.linspace(0,endpoint,161)
    refinement = float(np.max(np.abs(coarse.sol(sampled)-sol.sol(sampled))))
    assert refinement<2e-8
    return dict(model='Static finite-cutoff Israel-junction diagnostic: specified Polyakov throat joined abruptly to same-charge RN electrovac, outward increasing-radius exterior branch',
        units='R=physical_areal_radius/re, M=G*physical_mass/(c²*re); normal radius derivatives use proper radial distance. Natural units c=hbar=1: Sigma_hat=8piG re*surface_energy, p_hat=8piG re*surface_tangential_pressure; bulk_pressure_hat=8piG re²*T_quantum_nn. Restore1/c⁴ in all stress normalizations for SI; surface pressure is N/m and bulk pressure is Pa',
        definitions='v=exp(-sigma)*R_x>0, a_in=exp(-sigma)*sigma_x, w=sqrt(Fcut)>0, F=1-2M/R+1/R², A=1-1/R². A constant exterior time normalization matches the induced lapse but cannot change extrinsic curvature',
        junction='Sigma_hat=2(v-w)/R; p_hat=(A/w+w-2v)/(2R)-a_in. Same magnetic charge makes normal Maxwell stress continuous',
        shell_free_obstruction='A no-shell Darmois join requires continuous normal bulk pressure. In the normalized variables, quantum bulk_pressure_hat=-alpha*(sigma_x²+Omega²)/(3*R²*exp(2sigma))<0 for this NS source at every finite regular cutoff, while exact RN has zero quantum normal pressure. Thus abrupt source removal cannot be shell-free in this model',
        derivative_match='Setting w=v removes surface energy density but leaves p_hat=-R*bulk_pressure_hat/(2v)>0. This is a required surface stress, not an identified material or a stable shell',
        dominant_energy_condition_argument='For nonnegative surface energy,0<w<=v. If A>v², dp_hat/dw=(1-A/w²)/(2R)<0, so p_hat>=p_hat(v), while Sigma_hat<=2v/R. Every recorded cutoff has p_hat(v)>2v/R. At these recorded cutoffs no same-charge RN join in this outward branch satisfies the surface dominant energy condition Sigma>=|p|. This conditional obstruction does not prohibit other cutoffs, exteriors, continuing quantum matter or smooth transition regions',
        symbolic_pressure_derivative=str(dp), alpha=alpha, Omega_squared=frequency2,
        continued_CFT_on_RN_trace='Keeping the same CFT in exact RN instead of removing it also fails the coupled equations: R2=-F_second=4M/R³-6/R⁴ is nonzero on every open charged-RN radial region; t2_trace=q*R2/(24pi) is independent of the chiral state constants. RN may be an overlap approximation, not an exact quantum-source completion',
        RN_base_curvature=str(rn_R2),
        asymptotic_state_obstruction='Assume a globally spherical static extension of the same constant NS source to infinity with standard differentiable asymptotic flatness N->C>0,B->1. In dimensionless coordinates R=Rphys/re and time=tphys*c/re: rho_hat=8piG re²*rho_quantum/c⁴~-alpha*Omega²/(3*C²*R²), and d(m_geom/re)/d(Rphys/re)->-alpha*Omega²/(6*C²). This nonzero limit contradicts finite ADM mass for fixedalpha,Omega>0; it does not determine the actual alternative asymptotics',
        asymptotic_optical_length='The same asymptotic assumptions give integral_to_infinity B*dR/N=infinity, inconsistent with the prescribed finite optical circle. A finite returning magnetic path with nonspherical exterior geometry, or different source/state data, must replace this naive infinite radial continuation',
        cases=rows, maximum_normalized_interior_refinement_difference=refinement,
        quantum_state_scope='An Israel stress balance does not terminate the transparent fermion fields or supply their optical loop. The proposed abrupt matter removal and shell are diagnostic requirements, not a completed global quantum state. A real continuation must solve the fields, geometry and state together',
        scope='No full exterior solution with continued quantum stress, transition layer, shell material, stability, formation or time offset is constructed. RN formal horizons below the cutoff are discarded in this diagnostic and are not horizons of the actual unspecified joined interior. Surface DEC failure is not by itself a universal causality or time-travel theorem',
        sources=['https://arxiv.org/pdf/gr-qc/0310107 (Visser-Wiltshire Eq8: Israel convention; their surface tension is minus our tangential pressure)',
                 'https://arxiv.org/pdf/1807.04726 (MMP reference and fermion loop-state framework)',
                 'Stable mass-excess evaluation, finite-cutoff surface stresses and continuous DEC bound are our calculations'])


def quantum_taper_checks(static):
    """Conservation cost of tapering a fixed throat source, not a new solution."""
    from scipy.optimize import minimize_scalar
    x = sp.symbols('x',real=True)
    radius,sigma = sp.Function('R')(x),sp.Function('sigma')(x)
    tt,xx,perp = (sp.Function(name)(x) for name in ('t_tt','t_xx','p_perp'))
    theta = sp.symbols('theta',real=True)
    metric = (-sp.exp(2*sigma),sp.exp(2*sigma),radius**2,radius**2*sp.sin(theta)**2)
    mixed = (-sp.exp(-2*sigma)*tt/(4*sp.pi*radius**2),
        sp.exp(-2*sigma)*xx/(4*sp.pi*radius**2),perp,perp)
    divergence = sp.diff(mixed[1],x)+sum(sp.diff(g,x)/(2*g)*(mixed[1]-value)
        for g,value in zip(metric,mixed))
    ward = sp.diff(xx,x)+sp.diff(sigma,x)*(tt-xx)
    expected = sp.exp(-2*sigma)*ward/(4*sp.pi*radius**2)-2*sp.diff(radius,x)*perp/radius
    assert sp.simplify(divergence-expected)==0
    omega,w = sp.symbols('Omega_squared',positive=True),sp.Function('w')(x)
    s1,s2 = sp.diff(sigma,x),sp.diff(sigma,x,2)
    base_tt,base_xx = 2*s2-s1*s1-omega,-s1*s1-omega
    def reduced_ward(t0,t1):
        return sp.simplify(sp.diff(t1,x)+s1*(t0-t1))
    assert reduced_ward(base_tt,base_xx)==0
    state_tt,state_xx = 2*s2-s1*s1-omega*w,-s1*s1-omega*w
    state_force,full_force = -omega*sp.diff(w,x),-(s1*s1+omega)*sp.diff(w,x)
    assert sp.simplify(reduced_ward(state_tt,state_xx)-state_force)==0
    assert sp.simplify(reduced_ward(w*base_tt,w*base_xx)-full_force)==0
    assert sp.simplify(state_xx-state_tt+2*s2)==0
    assert sp.simplify(w*(base_xx-base_tt)+2*s2-2*(1-w)*s2)==0
    u = sp.symbols('u',real=True)
    window = 1-10*u**3+15*u**4-6*u**5
    assert window.subs(u,0)==1 and window.subs(u,1)==0
    assert sp.factor(sp.diff(window,u)+30*u*u*(1-u)**2)==0
    for order in (1,2):
        assert all(sp.diff(window,u,order).subs(u,end)==0 for end in (0,1))

    alpha,frequency2 = static['alpha'],static['Omega_squared']
    right = 1.2
    coarse,_ = _polyakov_static_solution(alpha,frequency2,right)
    fine,_ = _polyakov_static_solution(alpha,frequency2,right,2)
    def fields(solution,z):
        Y,V,Z,W = solution.sol(z)
        e2s = math.exp(2*alpha*Z)/math.cos(z)**2
        sp1 = math.tan(z)+alpha*W
        assert V>0
        return Y,V,e2s,sp1
    rows = []
    for end in (.8,1.,1.2):
        start=.5; width=end-start
        def data(z,solution=fine):
            Y,V,e2s,sp1 = fields(solution,z)
            uu=(z-start)/width
            weight=1-10*uu**3+15*uu**4-6*uu**5
            minus_derivative=30*uu*uu*(1-uu)**2/width
            denominator=3*V*e2s
            pressure=(frequency2*minus_derivative/denominator,
                      (sp1*sp1+frequency2)*minus_derivative/denominator)
            original_radial=-alpha*(sp1*sp1+frequency2)/(3*(1+alpha*Y)*e2s)
            return weight,pressure,original_radial,denominator,sp1*sp1+frequency2
        grid=np.linspace(start,end,401)
        values=np.array([data(z)[1] for z in grid]).T
        maxima=[]
        for index in (0,1):
            imax=int(np.argmax(values[index]))
            optimum=minimize_scalar(lambda z:-data(z)[1][index],
                bounds=(grid[imax-1],grid[imax+1]),method='bounded',options={'xatol':1e-14})
            assert optimum.success
            z=float(optimum.x);pressure=-float(optimum.fun)
            maxima.append(dict(x=z,required_angular_pressure_hat=pressure,
                ratio_to_absolute_original_radial_pressure=pressure/abs(data(z)[2])))
        state_bound=frequency2/quad(lambda z:data(z)[3],start,end,epsabs=1e-12,epsrel=1e-12)[0]
        full_bound=1/quad(lambda z:data(z)[3]/data(z)[4],start,end,epsabs=1e-12,epsrel=1e-12)[0]
        force_integrals=[quad(lambda z:data(z)[3]*data(z)[1][0],start,end,epsabs=1e-12)[0],
            quad(lambda z:data(z)[3]*data(z)[1][1]/data(z)[4],start,end,epsabs=1e-12)[0]]
        assert abs(force_integrals[0]-frequency2)<1e-11 and abs(force_integrals[1]-1)<1e-11
        assert maxima[0]['required_angular_pressure_hat']>=state_bound>0
        assert maxima[1]['required_angular_pressure_hat']>=full_bound>0
        ya,yb=fine.sol(start)[0],fine.sol(end)[0]
        delta_log_radius=.5*math.log1p(alpha*(yb-ya)/(1+alpha*ya))
        relative_bound=1/(2*delta_log_radius)
        # The logarithmic-radius integral avoids subtracting nearly equal radii.
        relative_integral=quad(lambda z:data(z)[1][1]/abs(data(z)[2])*
            alpha*fine.sol(z)[1]/(2*(1+alpha*fine.sol(z)[0])),start,end,epsabs=1e-12)[0]
        assert abs(relative_integral-.5)<1e-12
        refinement=max(abs(data(z,fine)[1][i]-data(z,coarse)[1][i]) for z in grid for i in (0,1))
        assert refinement<1e-8
        rows.append(dict(start_x=start,end_x=end,
            state_taper_peak=maxima[0],whole_tensor_taper_peak=maxima[1],
            any_monotone_state_taper_peak_pressure_lower_bound=state_bound,
            any_monotone_whole_taper_peak_pressure_lower_bound=full_bound,
            delta_log_areal_radius=delta_log_radius,
            any_monotone_whole_taper_max_pressure_ratio_lower_bound=relative_bound,
            weighted_force_integrals=force_integrals,log_radius_ratio_integral=relative_integral,
            maximum_background_refinement_difference=refinement,
            plot=dict(x=grid.tolist(),window=[data(z)[0] for z in grid],
                state_required_angular_pressure_hat=values[0].tolist(),
                whole_required_angular_pressure_hat=values[1].tolist(),
                original_radial_pressure_hat=[data(z)[2] for z in grid])))
    return dict(
        model='Necessary conservation and anomaly diagnostic for two imposed quantum-source tapers on the already solved static throat; geometry is held fixed',
        units='x is dimensionless conformal radius, R=Rphysical/re, R^2=1+alpha*Y,V=Y_prime. Natural units hbar=c=1 give alpha=G*q/(4*pi*re^2); in SI alpha=q*hbar*G/(4*pi*c^3*re^2). Pressure_hat=8*pi*G*re^2*p_perp/c^4 in SI',
        alpha=alpha,Omega_squared=frequency2,
        four_dimensional_Ward='D_a t^a_x=8*pi*R*R_prime*p_perp; equivalently t_xx_prime+sigma_prime*(t_tt-t_xx)=8*pi*R*R_prime*exp(2sigma)*p_perp, with t_ab=4*pi*R^2*T4_ab',
        direct_four_dimensional_divergence_residual='0',
        source_normalization='t_ab=q/(24*pi)*bar_t_ab; bar_t_tt=2*sigma_second-sigma_prime^2-Omega^2,bar_t_xx=-sigma_prime^2-Omega^2. Original action is radius-independent Polyakov with p_perp=0',
        effective_action_identity='(1/sqrt(-h))*delta W/delta R=8*pi*R*p_perp. A radius-independent W[h] gives p_perp=0; changing this requires deriving base stress, angular pressure and anomaly together. In the minimally coupled conformal theory, chiral state functions become constants under stationarity; an arbitrary static spatial Omega^2(x) is not such a state',
        state_taper=dict(prescription='Only replace Omega^2 by Omega^2*w(x)',
            scaled_Ward_force=str(state_force),
            required_pressure_hat='-Omega^2*w_prime/(3*V*exp(2sigma))',
            trace_anomaly_defect='0',
            limitation='A varying static state constant violates the original p_perp=0 Ward identity. Even w=0 leaves curvature-induced quantum stress, so it does not give exact electrovac RN'),
        whole_tensor_taper=dict(prescription='Replace the entire original t_ab by w(x)*t_ab',
            scaled_Ward_force=str(full_force),
            required_pressure_hat='-(sigma_prime^2+Omega^2)*w_prime/(3*V*exp(2sigma))',
            scaled_coordinate_trace_anomaly_defect='2*(1-w)*sigma_second; multiply by q*exp(-2sigma)/(24*pi) for the actual two-dimensional trace defect',
            limitation='Tapering all components changes the fixed-CFT trace anomaly. Adding a transverse pressure to restore 4D conservation does not repair this anomaly or derive a quantum effective action'),
        taper='w=1-10u^3+15u^4-6u^5 for u=(x-start)/(end-start) in[0,1], w=1 before and0 after. C2 profile; rows concern the right branch only',
        cases=rows,
        window_independent_bounds='For any monotone complete taper on the same interval/background: peak p_hat_state>=Omega^2/integral(3V e^(2sigma))dx; peak p_hat_whole>=1/integral[3V e^(2sigma)/(sigma_prime^2+Omega^2)]dx',
        relative_bound='p_perp/abs(p_radial_original)=-.5*dw/dlnR for a whole-tensor taper. Its integral over lnR is.5; hence max ratio>=1/[2*ln(Rend/Rstart)]. This concerns the original negative radial pressure, not the tapered pressure that vanishes at the endpoint',
        perturbation_order='R*R_prime=alpha*V/2 cancels alpha in the required pressure_hat. For fixed conformal layer and nonzero limiting V, this is order1 as alpha->0, while original radial quantum pressure is orderalpha. It is not a small correction to the held-fixed throat geometry',
        no_solution_claim='These tensors with compensating angular pressure satisfy the necessary Ward balance, but the original quantum action/state and fixed Einstein geometry are not solved by them. No material, stable shell, smooth join, formation or time-travel route is constructed',
        next_required='Derive source and angular pressure from a specified coupled quantum action/state and solve its backreaction, or use the actual nonspherical returning-field geometry. A fitted taper is not a substitute',
        sources=['https://arxiv.org/pdf/1807.04726 (MMP fermion loop and quantum support; not this taper)',
            'https://arxiv.org/pdf/hep-th/0309160 (Fabbri et al. Eqs37-40,42-43: dilaton-coupled scalar conservation, trace and state functions; without dilaton coupling the state functions are chiral. Not interchangeable with the present fermion source)',
            'The four-dimensional Ward specialization, two taper forces, bounds and fixed-background numerics are our derivations'])


def polyakov_static_checks(interior):
    """All static components for a specified reduced quantum source, on a finite strip."""
    import mpmath as mp
    from numpy.polynomial import chebyshev as cheb
    case = interior['cases']['0.1_16']['interior']
    alpha0 = case['alpha']
    beta_hat = case['beta']/alpha0
    omega2 = 1+3*beta_hat
    endpoint = math.atan(4.)
    samples = np.linspace(0, endpoint, 161)

    # A genuine constraint, not another evolution equation imposed at each step.
    r, rp, s, sp1, k, om = sp.symbols('r rp s sp1 k om', real=True)
    rpp = r/(r*r-k)*(sp.exp(2*s)*(1-1/r**2-k/r**4)-rp**2)
    spp = sp.exp(2*s)/r**4-rpp/r
    constraint = sp.exp(2*s)*(1-1/r**2)-rp**2-2*r*sp1*rp-k*sp1**2-k*om
    derivative = sum(sp.diff(constraint,v)*dv for v,dv in ((r,rp),(rp,rpp),(s,sp1),(sp1,spp)))
    assert sp.simplify(derivative) == 0
    z1,z2 = sp.symbols('sigma_prime sigma_second',real=True)
    ttt,txx = 2*z2-z1*z1-om, -z1*z1-om
    ward = sp.simplify(sp.diff(txx,z1)*z2+z1*(ttt-txx))
    trace = sp.simplify(txx-ttt+2*z2)
    assert ward == trace == 0

    def integrate(alpha, frequency2=omega2, refinement=1):
        return _polyakov_static_solution(alpha,frequency2,endpoint,refinement)

    def published_conformal(x):
        rho = math.tan(x)
        def gamma(z):
            return -(z*z+z*(3+z*z)*math.atan(z)-math.log1p(z*z))
        integral = quad(lambda z:gamma(z)/(1+z*z)**2, 0,rho,
                        epsabs=2e-13,epsrel=2e-13)[0]
        return rho*integral+gamma(rho)/(2*(1+rho*rho))
    lead_y = beta_hat*(1+samples*np.tan(samples))+1/3
    lead_z = np.asarray([beta_hat*published_conformal(x)-x*math.tan(x)/3 for x in samples])
    records = []
    for alpha in (alpha0, 1e-5, 5e-6, 2.5e-6):
        sol, rhs = integrate(alpha)
        Y,V,Z,W = sol.sol(samples)
        S = 1+alpha*Y
        sigmap = np.tan(samples)+alpha*W
        e2s = np.exp(2*alpha*Z)/np.cos(samples)**2
        H = e2s*Y/S-alpha*V*V/(4*S)-sigmap*V-(sigmap**2+omega2)/3
        fields = np.asarray([rhs(x,y) for x,y in zip(samples,sol.sol(samples).T)]).T
        A = (fields[1]-alpha*V*V/(2*S))/(2*S)
        curvature = -2/S**2+2*alpha*A/e2s
        uu = A-sigmap*V/S+(1/np.cos(samples)**2+alpha*fields[3]-sigmap**2-omega2)/(3*S)
        uv = -S*A-alpha*V*V/(4*S)+e2s*Y/S-(1/np.cos(samples)**2+alpha*fields[3])/3
        assert np.min(S-alpha/3)>0 and np.max(np.abs(H))<2e-9
        assert max(np.max(np.abs(uu)),np.max(np.abs(uv)))<2e-9
        records.append(dict(alpha=alpha, center_sphere_radius_squared_over_re2=float(S[0]),
            minimum_S_minus_alpha_over_3=float(np.min(S-alpha/3)),
            minimum_lapse_relative_to_center=float(np.min(np.sqrt(e2s))),
            maximum_H_over_alpha=float(np.max(np.abs(H))),
            maximum_null_constraint_residual_over_alpha=float(np.max(np.abs(uu))),
            maximum_mixed_residual_over_alpha=float(np.max(np.abs(uv))),
            maximum_scaled_area_difference_from_first_order=float(np.max(np.abs(Y-lead_y))),
            maximum_scaled_conformal_difference_from_first_order=float(np.max(np.abs(Z-lead_z))),
            plot=dict(x=samples.tolist(), area_excess_over_alpha=Y.tolist(),
                      conformal_correction_over_alpha=Z.tolist(),
                      lapse_relative_to_center=np.sqrt(e2s).tolist(),
                      re2_R2=curvature.tolist())))
    ratios = {}
    for field in ('maximum_scaled_area_difference_from_first_order',
                  'maximum_scaled_conformal_difference_from_first_order'):
        ratio = [records[i][field]/records[i+1][field] for i in (1,2)]
        assert all(1.98<v<2.02 for v in ratio)
        ratios[field] = ratio

    fine, _ = integrate(alpha0, refinement=2)
    main, _ = integrate(alpha0)
    refinement_error = float(np.max(np.abs(fine.sol(samples)-main.sol(samples))))
    assert refinement_error<2e-8
    # Differentiate a polynomial representation of solved fields, independently
    # of calling the evolution RHS to supply derivatives in the residual.
    spectral = []
    for degree in (40,60):
        nodes = np.cos(np.arange(degree+1)*math.pi/degree)
        xnodes = (nodes+1)*endpoint/2
        vals = fine.sol(xnodes)
        pts = np.linspace(.05,endpoint-.05,61)
        zpts = 2*pts/endpoint-1
        d = []
        for values in (vals[0],vals[2]):
            coef = cheb.chebfit(nodes,values,degree)
            d.append([cheb.chebval(zpts,cheb.chebder(coef,m))* (2/endpoint)**m for m in (0,1,2)])
        Y,Yp,Ypp = d[0]; Z,Zp,Zpp = d[1]
        S = 1+alpha0*Y
        sec2 = 1/np.cos(pts)**2
        e2s = sec2*np.exp(2*alpha0*Z)
        sp1 = np.tan(pts)+alpha0*Zp
        spp1 = sec2+alpha0*Zpp
        A = Ypp/(2*S)-alpha0*Yp*Yp/(4*S*S)
        enull = A-sp1*Yp/S+(spp1-sp1**2-omega2)/(3*S)
        emixed = -S*A-alpha0*Yp*Yp/(4*S)+e2s*Y/S-spp1/3
        eangular = Zpp+A-sec2*np.expm1(2*alpha0*Z-2*np.log1p(alpha0*Y))/alpha0
        spectral.append(dict(degree=degree,maximum_null_residual_over_alpha=float(np.max(np.abs(enull))),
            maximum_mixed_residual_over_alpha=float(np.max(np.abs(emixed))),
            maximum_angular_residual_over_alpha=float(np.max(np.abs(eangular)))))
    assert max(v for k,v in spectral[-1].items() if k!='degree')<2e-6

    # Exact anomaly-shifted product solution checks reference calibration.
    alpha = 1e-5
    Sref = (1+math.sqrt(1+4*alpha/3))/2
    ref_sol, _ = integrate(alpha, 1/Sref**2)
    with mp.workdps(45):
        am = mp.mpf(str(alpha)); sr = (1+mp.sqrt(1+4*am/3))/2
        exactY = float(1/(3*sr))
        exactZ = np.asarray([float((-mp.log(mp.cos(mp.mpf(float(x))/sr))+
                                  mp.log(mp.cos(mp.mpf(float(x)))))/am) for x in samples])
    ref_error_y = float(np.max(np.abs(ref_sol.sol(samples)[0]-exactY)))
    ref_error_z = float(np.max(np.abs(ref_sol.sol(samples)[2]-exactZ)))
    assert max(ref_error_y,ref_error_z)<2e-9
    return dict(model='Self-consistent static Einstein-Maxwell plus radius-independent Polyakov CFT on a finite symmetric strip; conditional reduced source model',
        action='I=I_spherical_Einstein_Maxwell - q/(96*pi)*integral sqrt(-h)*R2*Box_inverse*R2 + W_state[h], with radius-independent measure and all additional local dynamical source terms set to zero in this benchmark. Specified NS state. p_perp=0 follows from this action, not a fitted pressure',
        source='u=t+x,v=t-x; t_uu=t_vv=q/(48*pi)*(sigma_second-sigma_prime^2-Omega^2), t_uv=q*sigma_second/(48*pi). t_ab=4*pi*R^2*T4_ab; t2_trace=q*R2/(24*pi); p_perp=0',
        equations='Units re=1, k=alpha/3: R_second=R/(R^2-k)*[exp(2sigma)*(1-1/R^2-k/R^4)-R_prime^2]; sigma_second=exp(2sigma)/R^4-R_second/R. H=exp(2sigma)*(1-1/R^2)-R_prime^2-2R*sigma_prime*R_prime-k*sigma_prime^2-k*Omega^2=0',
        normalized_variables='R^2=1+alpha*Y, sigma=-log(cos x)+alpha*Z; solves both corrections without subtracting nearly equal unscaled fields',
        initial_data='Reflection symmetry R_prime=sigma_prime=0; sigma(0)=0; R(0)^2=1/(1-alpha*Omega^2/3), fixed by H=0',
        alpha=alpha0, beta_over_alpha=beta_hat, Omega_squared=omega2,
        symmetric_x_interval=[-endpoint,endpoint],
        assumed_total_conformal_NS_period=2*math.pi/math.sqrt(omega2),
        complementary_conformal_path_required=2*math.pi/math.sqrt(omega2)-2*endpoint,
        state_scope='Negative Casimir constant assumes a transparent antiperiodic fermion loop with total conformal length2pi/Omega. Only the throat strip is solved; positive complementary path length is necessary, not exterior matching or state preparation',
        reference='Exact product reference Sref=(1+sqrt(1+4alpha/3))/2, Omega_ref=1/Sref, sigma_ref=-log cos(x/Sref). At first order Yref=1/3,Zref=-x tan(x)/3. Subtract this reference shift and apply the radial gauge transformation to compare with MMP AppendixB',
        symbolic_constraint_derivative='0', symbolic_source_Ward_and_trace_residuals=['0','0'],
        residual_normalization='Reported enull,emixed,eangular are rescaled equations. For E=G-8piG(Tmag+TCFT), in re=1 units: E_uu/alpha=-enull/2; E_uv/alpha=emixed/(2S); E^theta_theta/alpha=exp(-2sigma)*eangular. The vv component equals uu by static symmetry. These are mean-field equation residuals, not quantum fluctuations or hardware errors',
        cases=records, scaled_first_order_error_halving_ratios=ratios,
        maximum_normalized_state_refinement_difference=refinement_error,
        independently_differentiated_residuals=spectral,
        exact_product_reference_errors=dict(area_over_alpha=ref_error_y,conformal_correction_over_alpha=ref_error_z),
        leading_comparison=dict(x=samples.tolist(),area_excess_over_alpha=lead_y.tolist(),conformal_correction_over_alpha=lead_z.tolist()),
        scope='All static components of the stated reduced-source equations are checked. This is not full4D renormalized fermion stress, evolving stability, corrected pulse/return clocks, exterior matching, topology formation, physical material production, or proof of past travel. Larger-alpha cases are numerical convergence tests, not validated weak-coupling physical designs',
        sources=['https://arxiv.org/pdf/1807.04726 (AppendicesB,C,F and reference-radius discussion; source and comparison framework)',
                 'https://arxiv.org/pdf/hep-th/0006035 (Eqs5.1-5.5: Polyakov action and state-dependent source conventions)',
                 'https://arxiv.org/pdf/hep-th/0111091 (spherical Einstein-Maxwell action)',
                 'Closed static IVP, normalized implementation and finite-strip checks are our conditional-model calculations'])


def angular_metric_checks(interior):
    """Published static metric correction and the conditional curvature response."""
    import mpmath as mp
    rho, b = sp.symbols('rho b', real=True)
    f0 = 1+rho**2
    profile = 1+rho*sp.atan(rho)
    correction = -(rho**2+rho*(3+rho**2)*sp.atan(rho)-sp.log(1+rho**2))
    static_dilaton = sp.simplify(f0*sp.diff(profile,rho,2)+2*rho*sp.diff(profile,rho)-2*profile)
    published_equation = sp.simplify(sp.diff(correction,rho,2)+6*profile)
    assert static_dilaton == published_equation == 0
    sphere, f = 1+b*profile, f0+b*correction
    # re^2*(G^theta_theta - 8*pi*G*Tmag^theta_theta), with p_perp=0.
    angular = (sp.diff(f,rho,2)/2+(f*sp.diff(sphere,rho,2)+sp.diff(f,rho)*sp.diff(sphere,rho))/(2*sphere)
               -f*sp.diff(sphere,rho)**2/(4*sphere**2)-1/sphere**2)
    assert sp.simplify(angular.subs(b,0)) == 0
    assert sp.simplify(sp.diff(angular,b).subs(b,0)) == 0
    mixed = correction*sp.diff(profile,rho,2)+sp.diff(correction,rho)*sp.diff(profile,rho)
    numerator2 = -16*profile**2+2*mixed-f0*sp.diff(profile,rho)**2
    numerator3 = -12*profile**3+2*profile*mixed-correction*sp.diff(profile,rho)**2
    factored = b*b*(numerator2+b*numerator3)/(4*sphere*sphere)
    assert sp.simplify(sp.expand((angular-factored)*4*sphere*sphere)) == 0
    # An independent integration of the angular ODE checks the closed solution.
    sample_rho = np.linspace(0,4,81)
    ode = solve_ivp(lambda r,y:(y[1],-6*(1+r*math.atan(r))), (0,4), (0,0),
                    t_eval=sample_rho, method='DOP853', rtol=2e-12, atol=2e-13, max_step=.1)
    assert ode.success
    expected = -(sample_rho**2+sample_rho*(3+sample_rho**2)*np.arctan(sample_rho)-np.log1p(sample_rho**2))
    ode_error = float(np.max(np.abs(ode.y[0]-expected)))
    assert ode_error<1e-9
    with mp.workdps(55):
        pp = sp.lambdify(rho,profile,'mpmath')
        gg = sp.lambdify(rho,correction,'mpmath')
        corrected_residual = sp.lambdify((rho,b),factored,'mpmath')
        old_residual = sp.lambdify((rho,b),
            b*(12*profile+b*(8*profile**2-f0*sp.diff(profile,rho)**2))/(4*sphere**2),'mpmath')
        beta = mp.mpf(str(interior['cases']['0.1_16']['interior']['beta']))
        def evaluate(r, amplitude):
            phi, gamma = amplitude*pp(r), amplitude*gg(r)
            base = 1+r*r
            assert 1+phi>0 and base+gamma>0
            return dict(rho=float(r), phi=float(phi), relative_metric_change=float(gamma/base),
                relative_same_coordinate_lapse_change=float(mp.sqrt(1+gamma/base)-1),
                required_dimensionless_2D_curvature_change=float(6*phi),
                uncorrected_angular_residual=float(old_residual(r,amplitude)),
                corrected_angular_residual=float(corrected_residual(r,amplitude)))
        points = [evaluate(mp.mpf(r),beta) for r in (0,1,4)]
        scaling = []
        for r in (0,1,4):
            values = [abs(corrected_residual(mp.mpf(r),amp)) for amp in map(mp.mpf,('0.0001','0.00005','0.000025'))]
            ratios = [float(values[k]/values[k+1]) for k in (0,1)]
            assert all(3.98<ratio<4.02 for ratio in ratios)
            scaling.append(dict(rho=r, residuals=list(map(float,values)), halving_ratios=ratios))
        cutoffs = []
        for phi_cut in map(mp.mpf,('0.00001','0.0001')):
            rc = mp.findroot(lambda r:beta*pp(r)-phi_cut, 2*phi_cut/(mp.pi*beta))
            row = evaluate(rc,beta)
            # Integrate in asinh(rho), resolving both center and overlap region.
            xmax = mp.asinh(rc)
            def time_integrand(x, resum):
                r = mp.sinh(x); base=1+r*r; gamma=beta*gg(r)
                return -2*mp.cosh(x)*gamma/(base*(base+gamma if resum else base))
            linear_time = mp.quad(lambda x:time_integrand(x,False),[0,1,xmax])
            resummed_time = mp.quad(lambda x:time_integrand(x,True),[0,1,xmax])
            assert 0<linear_time<resummed_time and abs(row['relative_metric_change'])<.001
            row.update(phi_cut=float(phi_cut),
                static_uncorrected_cutoff_to_cutoff_time=float(2*mp.atan(rc)),
                first_order_chart_transit_time_increase=float(linear_time),
                resummed_truncated_chart_transit_time_increase=float(resummed_time),
                resummation_difference=float(resummed_time-linear_time))
            cutoffs.append(row)
    pulse_curvature = {}
    for key,case in interior['cases'].items():
        data=case['plot']; opening=case['interior']['beta']
        change=6*opening*np.asarray(data['center_area_excess_change_fraction'])
        pulse_curvature[key]=dict(maximum_abs_conditional_delta_re2_R2=float(np.max(np.abs(change))),
            plot=dict(global_center_time=data['center_global_time'],conditional_delta_re2_R2=change.tolist()))
    return dict(
        model='Angular Einstein-Maxwell correction benchmark in the MMP reference-subtracted, zero-additional-angular-pressure sector; not a full4D quantum-source solve',
        exact_equation='Box_h R/R - R2[h]/2 - re^2/R^4 = 8*pi*G*p_perp, for ds4^2=h_ab dx^a dx^b+R^2 dOmega^2 and fixed magnetic charge',
        reduced_action='Igrav=(1/(4G))*integral sqrt(-h)[R^2*R2+2*(nabla R)^2+2-2*re^2/R^2], up to boundary terms; add a specified matter effective action W[h,R]',
        general_angular_residual='For h=exp(2psi)h0,Q=1+phi,H=re^2*(nabla phi)_0^2,D=re^2*Box0: Eangular=exp(-2psi)[1+Dpsi+Dphi/(2Q)-H/(4Q^2)]-1/Q^2-P. For frozen h0,Dphi=2phi: Eangular=(3phi+2phi^2-H/4)/Q^2-P',
        static_metric='ds4^2=re^2[-f(rho)dTau^2+drho^2/f(rho)+(1+phi)dOmega^2], f=1+rho^2+gamma. phi=beta*(1+rho*atan(rho)); gamma=-beta*[rho^2+rho*(3+rho^2)*atan(rho)-log(1+rho^2)], gamma(0)=gamma_prime(0)=0',
        angular_residual='re^2*(G^theta_theta-8*pi*G*Tmag^theta_theta) at p_perp0 equals f_second/2+(f*S_second+f_prime*S_prime)/(2S)-f*S_prime^2/(4S^2)-1/S^2, S=1+phi. This is one component, not the complete field-equation residual',
        symbolic_verification=dict(static_dilaton_equation=str(static_dilaton),published_gamma_equation=str(published_equation),
            zero_and_first_order_angular_residuals_vanish=True,exact_quadratic_factorization_checked=True),
        normalized_gamma_ODE_error=ode_error, beta=float(beta), compact_points=points,
        quadratic_residual_scaling=scaling, finite_overlap_cutoffs=cutoffs,
        dynamic_relation='For h=exp(2psi)h0 and D=re^2*Box0: (D-2)psi=-Dphi/2-2phi+P, P=8*pi*G*re^2*p_perp. If Dphi=2phi, delta(re^2 R2[h])=6phi-2P to first order. Relative to a calibrated static reference: Delta(re^2 R2)=6*Delta_phi-2*Delta_P',
        conditional_pulse_curvature=pulse_curvature,
        calibration='MMP removes metric-proportional, state-independent quantum stress into the reference radii. The absolute unsubtracted CFT trace cannot be added to this benchmark without recalibrating those radii. The static resummed metric is only a test of a first-order formula, not an all-orders solution',
        quantum_source_gap='The full angular pressure is not fixed by existing chiral mean stress alone. With t_ab=4*pi*R^2*T4_ab and D_R=(1/sqrt(-h))*deltaW/deltaR=8*pi*R*p_perp, conservation requires nabla_a t^a_b=D_R*partial_b R. A specified reduced matter effective action and its4D renormalization/matching are needed',
        source_ambiguity_example='Wlambda=lambda*integral sqrt(-h)*R2*log(R/R0) contributes t_ab=2lambda*(nabla_a nabla_b-h_ab*Box)logR and p_perp=lambda*R2/(8*pi*R^2). At constant radius its base stress vanishes while angular pressure does not. This demonstrates missing matching data, not a free physical parameter to tune for travel',
        reference_trace_offset='At a symmetric static center,1-1/S0=(4G/re^2)*(t_uv-t_uu). Full2D CFT gives alpha*Omega^2/3; trace-subtracted source gives alpha*(Omega^2-1)/3. Their alpha/3 difference is the reference-radius convention, not a clock-solver error',
        pulse_condition='The plotted pulse-minus-static curvature uses Delta_P=0 at leading order in opening and the omitted loop corrections, on a fixed compact interior. For smooth local radius terms with coefficients of orderq, Delta_P is orderalpha^2*epsilon there; large coefficients, angular excitation, enhanced derivatives or an alpha-dependent cutoff can invalidate that ordering. The general relation with Delta_P is retained',
        scope='The angular correction removes its first-order residual and the remaining angular residual scales quadratically on fixed compact regions. Other Einstein components, changed quantum state/return maps, moving-cutoff boundary data, corrected exteriors and formation remain unsolved. Static coordinate transit corrections do not establish a physical time offset or calendar return',
        sources=['https://arxiv.org/pdf/1807.04726 (Eq5.34, reference-radius convention afterEq5.38, AppB EqsB.4-B.5; published static gamma)',
                 'https://arxiv.org/pdf/hep-th/0111091 (Eq9: spherical Einstein-Maxwell reduction)',
                 'https://arxiv.org/pdf/hep-th/9901119 (section2: radius/dilaton-dependent anomaly and finite-renormalization conditions)',
                 'Exact warped-product reduction, factorized residual, finite-cutoff chart-time integrals and calibrated pulse-curvature diagnostic here are our calculations'])


def nonlinear_interior_checks():
    """Longer and stronger pulse checks on the same leading nonlinear closure."""
    records, states = {}, {}
    for epsilon, loops in ((0., 4), (.1, 16), (.3, 4), (.5, 4)):
        key = f'{epsilon}_{loops}'
        records[key], states[key] = nonlinear_pulse_solution(epsilon, loops=loops, interior=True)
        print(f'Interior completed: epsilon={epsilon}, loops={loops}', flush=True)
    refinements = {}
    for epsilon, loops in ((.1, 16), (.5, 4)):
        key = f'{epsilon}_{loops}'
        fine, fine_states = nonlinear_pulse_solution(epsilon, .5, loops=loops, interior=True)
        difference = float(np.max(np.abs(states[key]-fine_states)))
        assert difference < 3e-8
        coarse_cut = records[key]['interior']['finite_cutoff_checks']
        fine_cut = fine['interior']['finite_cutoff_checks']
        refinements[key] = dict(maximum_full_state_difference=difference,
            refined_bulk_matching_error=fine['maximum_bulk_boundary_rate_mismatch'],
            refined_energy_ledger_error=fine['maximum_reduced_energy_ledger_error'],
            maximum_cutoff_lapse_diagnostic_difference=max(abs(a['maximum_relative_lapse_mismatch']-
                b['maximum_relative_lapse_mismatch']) for a,b in zip(coarse_cut, fine_cut)),
            maximum_local_source_diagnostic_difference=abs(records[key]['interior']['maximum_sampled_chiral_source_strength']-
                fine['interior']['maximum_sampled_chiral_source_strength']))
        print(f'Interior refinement passed: epsilon={epsilon}, loops={loops}', flush=True)
    zero = records['0.0_4']['interior']
    assert zero['maximum_sampled_fractional_total_area_change'] < 1e-20
    assert zero['future_trapped_grid_points'] == zero['past_trapped_grid_points'] == 0
    assert all(p['spatial_second_derivative_over_alpha']>0 for p in zero['spatial_stationary_points'])
    return dict(model='Interior and finite-cutoff diagnostics for nonlinear_pulse; same initial assumed throat and ideal reservoirs',
        parameters=dict(exterior_delay=.1, pulse_width=.4, alpha=1e-9, magnetic_charge_q=100000),
        field='Phi=phi/alpha=betaHat*(1+chi*tan(chi))+(A(u)-B(v))*tan(chi)+A_u(u)+B_v(v), u=tau+chi,v=tau-chi. Both histories are restricted to completed causal data',
        expansions='Future null expansions have signs Phi_u,Phi_v; with kplus dot kminus=-1, theta_plus/minus=sqrt(2)*alpha*cos(chi)*Phi_u/v/[re*(1+alpha*Phi)]. Static center is doubly marginal. Phi=0 is extremal area, not sphere disappearance',
        cutoff='Solve phi(tau,chi)=phi_cut in delta=pi/2-|chi|; chi_dot=-Phi_tau/Phi_chi. Proper lapse relative to the exact static cutoff is clock_rate*sin(delta0)/sin(delta)*sqrt(1-chi_dot^2). This is a leading-JT finite-cutoff diagnostic, not a solved finite-cutoff boundary condition',
        cases=records, step_refinements=refinements,
        validity='Small phi on declared interior samples and beta<<delta<<1 at the tested cutoffs support the area/overlap hierarchy. The small cutoff lapse residual is not full4D accuracy: omitted metric corrections can be order phi_cut (1e-5 or1e-4). Initial Schwarzian scale re^3/(G ell0)=q*betaHat/16=5743.78 is not a bound on later fluctuations. Classical carrier bandwidth is not total pulse energy and not a quantum-state validity test',
        scope='Finite-time mean-field diagnostics, with a sampled interior and leading metric. No global stability, quantum fluctuations, topology formation, physical reservoir/particle transport, backward time mapping or actual apparatus is established. Contracting/expanding trapped spheres alone do not determine event horizons under NEC violation',
        sources=['https://arxiv.org/pdf/1807.04726 (MMP section5.1, section5.4 and appendixB: throat expansion, energy and matching limits)',
                 'https://arxiv.org/pdf/1606.01857 (MSY boundary proper-metric/dilaton conditions)',
                 'https://arxiv.org/abs/gr-qc/9805019 (Hayward: dynamic wormholes and timelike trapping horizons). This coupled solution and diagnostic calculations are ours'])



def nonlinear_pulse_checks(linear):
    """Dynamical leading JT clocks, retarded vacuum, and an independently launched carrier."""
    delay, width = .1, .4
    length = math.pi+delay
    omega2 = (2*math.pi/length)**2
    baseline = (1-omega2)/6
    b0 = math.pi*(omega2-1)/6  # All JT sources, area coefficients and charges / alpha.
    final = 4*length
    samples = np.asarray(linear['plot']['x'])
    samples = samples[samples <= final+1e-12]
    reference = np.array([linear['plot'][key][:len(samples)]
                          for key in ('right_clock_rate', 'left_clock_rate')]).T

    solve = lambda epsilon, refinement=1.: nonlinear_pulse_solution(epsilon, refinement, samples=samples)

    records, values = {}, {}
    for epsilon in (0., .01, .005, .0025, .1):
        records[str(epsilon)], values[epsilon] = solve(epsilon)
    refined, fine = solve(.01, .5)
    refinement_error = float(np.max(np.abs(values[.01]-fine)))
    assert refinement_error < 3e-8
    strong_refined, strong_fine = solve(.1, .5)
    strong_refinement_error = float(np.max(np.abs(values[.1]-strong_fine)))
    assert strong_refinement_error < 3e-8
    assert records['0.0']['maximum_sampled_clock_rate_change']<1e-12
    corrections = []
    for epsilon in (.01, .005, .0025):
        difference = float(np.max(np.abs(values[epsilon][:, [1, 5]]-1-epsilon/.01*(reference-1))))
        corrections.append(dict(epsilon=epsilon, maximum_clock_rate_difference_from_linear=difference,
                                difference_over_epsilon_squared=difference/epsilon**2))
    assert corrections[1]['maximum_clock_rate_difference_from_linear']<.4*corrections[0]['maximum_clock_rate_difference_from_linear']
    assert corrections[2]['maximum_clock_rate_difference_from_linear']<.4*corrections[1]['maximum_clock_rate_difference_from_linear']
    return dict(
        model='Nonlinear clocks and transported CFT mean state, self-consistent only within the stated leading JT transparent-loop ansatz. Separate classical one-pass scalar/null pulse with ideal launch/absorption reservoirs. Fixed common exterior delay and fixed leading AdS2 chart',
        parameters=dict(exterior_delay=delay, pulse_width=width, final_exterior_time=final,
                        background_loop_lengths=4, b0_over_alpha=b0, initial_j_over_alpha=baseline),
        drive='Independently specified right exterior input I_R/alpha=(pi*epsilon/2)*8/(3w)*sin^4(pi*s/w) for0<s<w. Integral=pi*epsilon/2. Left input0. Positive C3 flux; no prescribed mouth clock or retrospectively defined force',
        carrier='Classical massless scalar/null fluid: j_carrier(r(s)+pi/2)/alpha=I_R(s)/(alpha*r_prime(s)^2), constant along u, absorbed on first exit at left mouth. No classical Schwarzian term and no carrier in the circulating vacuum recurrence. This does not specify an arbitrary NS quantum excitation',
        clocks='Each side evolves (r,v=r_prime,a=r_second,K={r,s}+v^2/2): r_prime=v,v_prime=a,a_prime=v*(K-v^2/2)+1.5*a^2/v; K_prime=v^2*(J_out-J_in)/(b0/alpha), with J=j_total/alpha. Static past r=l=s,v1,a0,K.5',
        vacuum='For Rin, J_vac=1/6+[l_prime(s-D)/r_prime(s)]^2*(J_old_plus_out-1/6)+(S_r-S_l_old)/(3*r_prime^2), S=K-v^2/2; swap for Lin. Old outgoing source obtained from its earlier same-chirality entrance. Initial NS vacuum retained through retarded recurrence',
        bulk='Along incoming coordinates u=r(s)+h,v=l(s)+h integrate A,A_u,A_uu and B,B_v,B_vv, with P A=-(J_plus-J0),P B=-(J_minus-J0). Reconstructed boundary coefficients independently test r_prime=1+[A(r+h)-B(r-h)]/B0 and l_prime=1+[-A(l-h)+B(l+h)]/B0',
        energy='E_R,L/alpha=-B0*K_R,L; external CFT outgoing Q=v^2*(J_vac_out-1/6)-S/3. Reduced ledger: Delta(E_R+E_L)/alpha + integral_(s-D)^s delta(Q_Rout+Q_Lout)dt = integral_0^s(I_R-I_abs)/alpha dt. No4D ADM-energy construction. Multiply normalized ledger byq/(8pi) for dimensionless exterior physical energy',
        source_energy='Input dimensionless exterior energy E_s=q*epsilon/16. Physical natural-unit energy E_s/ell0. epsilon is not automatically the finite-D binding-energy fraction. Reservoir source and absorber are ideal boundary channels, not a built apparatus',
        verification=dict(full_state_step_refinement_error=refinement_error,
                          refined_bulk_matching_error=refined['maximum_bulk_boundary_rate_mismatch'],
                          refined_energy_ledger_error=refined['maximum_reduced_energy_ledger_error'],
                          stronger_pulse_full_state_step_refinement_error=strong_refinement_error,
                          stronger_refined_bulk_matching_error=strong_refined['maximum_bulk_boundary_rate_mismatch'],
                          stronger_refined_energy_ledger_error=strong_refined['maximum_reduced_energy_ledger_error'],
                          linear_limit_comparisons=corrections),
        cases=records,
        scope='Four background loops, small classical input, quantum mean-field only. Sampled positive clock rates and return times do not prove global nonlinear stability or a chronology theorem. No finite-cutoff/4D matching, metric corrections, channel-dependent delays, scalar exterior transmission, reservoir gravity, state fluctuations, formation, material preparation, capture or calendar mapping. No continuation when throat histories cease to lie in completed causal data',
        sources=['https://arxiv.org/pdf/1606.01857 (MSY Eq4.22: boundary energy-flux/redshift factors)',
                 'https://arxiv.org/pdf/1402.6334 (AP section3.2: classical chiral scalar sources and dilaton response)',
                 'https://arxiv.org/pdf/math-ph/0412028 (conformal stress transformation); coupled finite-delay implementation and ledger here are our derivation'])


def run_checks():
    u, chi = sp.symbols('u chi', real=True)
    # Chain rule for chi=(u-v)/2, keeping u and chi independent during algebra.
    du = lambda expr: sp.diff(expr, u)+sp.diff(expr, chi)/2
    dv = lambda expr: -sp.diff(expr, chi)/2
    A = sp.Function('A')(u)
    perturbation = A*sp.tan(chi)+sp.diff(A, u)
    expressions = (
        du(du(perturbation))-sp.tan(chi)*du(perturbation)
        -sp.diff(A, u, 3)-sp.diff(A, u),
        dv(dv(perturbation))+sp.tan(chi)*dv(perturbation),
        du(dv(perturbation))+perturbation/(2*sp.cos(chi)**2))
    residuals = [sp.simplify(sp.trigsimp(expr)) for expr in expressions]
    assert all(expr == 0 for expr in residuals)
    a0 = sp.symbols('alpha', positive=True)
    background = a0*(1+chi*sp.tan(chi))
    background_residuals = [sp.simplify(sp.trigsimp(expr)) for expr in (
        du(du(background))-sp.tan(chi)*du(background)-a0/2,
        dv(dv(background))+sp.tan(chi)*dv(background)-a0/2,
        du(dv(background))+background/(2*sp.cos(chi)**2))]
    assert all(expr == 0 for expr in background_residuals)
    assert sp.simplify(8*sp.cos(chi)**2*du(dv(-sp.log(sp.cos(chi))))+2)==0
    re, G, q = sp.symbols('re G q', positive=True)
    ell, binding, alpha_mmp = 16*re**3/(G*q), G*q**2/(256*re**3), q*G/(4*sp.pi*re**2)
    assert sp.simplify(2*G*ell*binding/(re**2*alpha_mmp)-sp.pi/2)==0
    alpha, epsilon, width, charge = 1e-9, .01, .4, 100000
    coupling = math.sqrt(4*math.pi**2*charge*alpha)
    impulse = math.pi*alpha*epsilon/2
    def packet(s):
        return 8/(3*width)*math.sin(math.pi*s/width)**4 if 0<s<width else 0.
    assert abs(quad(packet, 0, width, epsabs=1e-13)[0]-1)<1e-12
    # ponytail: integrate A/J, not tiny absolute amplitudes; full boundary/state dynamics are separate.
    end = math.pi
    samples = np.unique(np.r_[np.linspace(0, end, 101), width, math.pi/2+width/2])
    solution = solve_ivp(lambda s, y: (y[1], y[2], -y[1]-packet(s)),
                         (0, end), (0, 0, 0), t_eval=samples, dense_output=True,
                         method='DOP853', rtol=2e-12, atol=2e-13, max_step=width/25)
    assert solution.success
    errors = []
    for time in (0, width/4, width, 1., math.pi/2+width/2, end):
        kernels = (lambda x: 1-math.cos(x), math.sin, math.cos)
        expected = [-quad(lambda s: packet(s)*kernel(time-s), 0, min(time, width),
                          epsabs=1e-13, epsrel=1e-12)[0] for kernel in kernels]
        errors.append(float(np.max(np.abs(solution.sol(time)-expected))))
    fourier_factor = (math.sin(width/2)/(width/2)*64*math.pi**4
                      /((4*math.pi**2-width**2)*(16*math.pi**2-width**2)))
    after = samples >= width
    phase = samples[after]-width/2
    closed = np.array([-1+fourier_factor*np.cos(phase),
                       -fourier_factor*np.sin(phase), -fourier_factor*np.cos(phase)])
    closed_error = float(np.max(np.abs(solution.y[:, after]-closed)))
    assert max(errors+[closed_error])<1e-10
    peak_time = math.pi/2+width/2
    center_change = -impulse*fourier_factor
    assert abs(impulse*solution.sol(peak_time)[1]-center_change)<impulse*1e-10
    feedback = boundary_feedback_checks(alpha, epsilon, width)
    loops = loop_return_checks(alpha, epsilon, width, feedback)
    result = dict(
        model='Linear spherical MMP/JT near-throat area response, fixed leading AdS2 conformal metric and prescribed background Casimir stress plus a classical chiral pulse. Not full four-dimensional semiclassical evolution',
        conventions='chi=atan(rho), u=tau+chi, v=tau-chi; ds2=re^2*sec(chi)^2*(-dtau^2+dchi^2); sphere area=4*pi*re^2*(1+phi), phi0=alpha*(1+chi*tan(chi))',
        source='j(u)=8*pi*G4*deltaT4_uu=(2*G4/re^2)*F(u), F=t2_uu. Only uu source is added; it is conserved and traceless in the leading 2D metric. No prepared quantum pulse state is claimed',
        solution='delta_phi=A(u)*tan(chi)+A_prime(u); A_tripleprime+A_prime=-j; retarded A=A_prime=A_second=0 beforeu0. A=-integral_0^u j(s)*(1-cos(u-s))ds',
        equations=['delta_phi_uu-tan(chi)*delta_phi_u=-j(u)',
                   'delta_phi_vv+tan(chi)*delta_phi_v=0',
                   'delta_phi_uv+delta_phi/(2*cos(chi)^2)=0'],
        symbolic_residuals=[str(expr) for expr in residuals],
        background_symbolic_residuals=[str(expr) for expr in background_residuals],
        energy_normalization='J=integral jdu=(pi*alpha/2)*epsilon, epsilon=E_infinity/abs(E_bind), E_tau=ell*E_infinity. Uses ell=16*re^3/(G*q) and alpha=4*re/(pi*ell)',
        parameters=dict(alpha=alpha, epsilon=epsilon, width=width, magnetic_charge_q=charge,
                        gauge_coupling_g=coupling, J=impulse, pulse_global_energy=epsilon*charge/16),
        parameter_scope='alpha=g^2/(4*pi^2*q); q=100000 and g~.0628 avoid choosing a small-alpha example incompatible with the weak-coupling/large-q hierarchy. This does not establish full approximation control, quantum fluctuations, formation or material availability',
        packet='j=J*8/(3*w)*sin(pi*u/w)^4 for0<u<w, else0. Nonnegative compact C3 pulse with integralJ; not a C-infinity quantum-state construction',
        post_pulse='C=sinc(w/2)*64*pi^4/[(4*pi^2-w^2)*(16*pi^2-w^2)], sinc(x)=sin(x)/x; A/J=-1+C*cos(u-w/2), A_prime/J=-C*sin(u-w/2)',
        fourier_factor=fourier_factor,
        peak_center=dict(global_time=peak_time, delta_phi=center_change,
                         delta_phi_over_alpha=center_change/alpha,
                         fractional_sphere_area_change=center_change/(1+alpha)),
        boundary='Formal throat-matching coefficients bR,L=lim_coschi->0+ cos(chi)*phi; baseline b0=pi*alpha/2. delta_bR(tau)=A(tau+pi/2), delta_bL(tau)=-A(tau-pi/2). A_prime contributes no leading coefficient. This asymptotic extraction is not a uniform small-phi four-dimensional limit; physical matching needs an overlap region and moving cutoffs',
        boundary_result='Clamping either leading coefficient to its old value at all global times while retaining this retarded solution forces A=0, hence j=0. A nonzero pulse needs boundary/exterior response. This is not an impossibility theorem for moving-cutoff JT boundary conditions or actual wormholes',
        boundary_sample_at_global_time_half=dict(
            right_coefficient_over_background=1+epsilon*float(solution.sol(.5+math.pi/2)[0]),
            left_coefficient_over_background=1.,
            scope='First passage with the left mouth not yet reached. Boundary/state response has not been fed back into the source'),
        independent_convolution_errors=errors, independent_post_pulse_error=closed_error,
        sampled_response=dict(u=samples.tolist(), A_over_J=solution.y[0].tolist(),
                              A_prime_over_J=solution.y[1].tolist()),
        remaining='Boundary trajectories, matching to both exteriors, Casimir-state change, corrected 2D metric and higher-order four-dimensional backreaction are not solved. Small interior area change is not a stability or transmission proof; no closure threshold or calendar mapping',
        boundary_feedback=feedback,
        loop_return=loops,
        zero_delay_spectrum=zero_delay_spectrum_checks(),
        time_offset=time_offset_checks(),
        stationary_balance=stationary_balance_checks(),
        driven_clock=driven_clock_checks(),
        nonlinear_pulse=nonlinear_pulse_checks(loops['cases']['0.1_one_pass']),
        nonlinear_interior=nonlinear_interior_checks(),
        sources=['https://arxiv.org/pdf/1807.04726 (MMP, equations5.31 and5.34-5.40, background/source normalization)',
                 'https://arxiv.org/pdf/1402.6334 (Almheiri-Polchinski, dilaton-gravity constraints; not a proof of this MMP boundary closure)'])
    result['angular_metric'] = angular_metric_checks(result['nonlinear_interior'])
    result['polyakov_static'] = polyakov_static_checks(result['nonlinear_interior'])
    result['polyakov_junction'] = polyakov_junction_checks(result['polyakov_static'])
    result['magnetic_return_paths'] = magnetic_return_path_checks()
    result['optical_matching'] = optical_matching_checks(result['magnetic_return_paths'])
    result['weak_magnetic_response'] = weak_magnetic_response_checks(result['optical_matching'])
    result['magnetic_boundary_response'] = magnetic_boundary_response_checks(result['polyakov_static'])
    result['coupled_dipole'] = coupled_dipole_checks(result['magnetic_boundary_response'])
    result['coupled_exterior'] = coupled_exterior_checks(result['coupled_dipole'])
    result['coupled_boundary'] = coupled_boundary_checks(result['coupled_exterior'])
    result['quantum_taper'] = quantum_taper_checks(result['polyakov_static'])
    result['magnetic_loop_source'] = magnetic_loop_source_checks(result['magnetic_return_paths'], result['optical_matching'])
    result['magnetic_patch_response'] = magnetic_patch_response_checks(result['magnetic_loop_source'], result['optical_matching'])
    result['sourced_rn_transfer'] = sourced_rn_transfer_checks(result['magnetic_loop_source'], result['optical_matching'], result['coupled_exterior'])
    result['companion_incident'] = companion_incident_checks(result['optical_matching'], result['coupled_exterior'], result['sourced_rn_transfer'])
    result['rotating_source'] = rotating_source_checks(result['magnetic_loop_source'], result['optical_matching'], result['companion_incident'])
    result['loop_charge_response'] = loop_charge_response_checks(result['magnetic_loop_source'], result['optical_matching'], result['rotating_source'])
    result['loop_retarded_response'] = loop_retarded_response_checks(result['loop_charge_response'])
    result['throat_charge_response'] = throat_charge_response_checks(result['rotating_source'], result['optical_matching'])
    result['charged_mouth_matching'] = charged_mouth_matching_checks(result['throat_charge_response'], result['optical_matching'])
    result['nonstationary_support'] = nonstationary_support_checks(result['stationary_balance'])
    result['plateau_gravity_response'] = plateau_gravity_response_checks(result['nonstationary_support'], result['parameters'])
    result['quantum_seed_feedback'] = quantum_seed_feedback_checks()
    result['causal_foliation'] = causal_foliation_checks(result['quantum_seed_feedback'])
    result['lambda_feedback'] = lambda_feedback_checks(result['quantum_seed_feedback'])
    result['wide_packet'] = wide_packet_checks(result['lambda_feedback'])
    result['long_packet'] = long_packet_checks(result['wide_packet'])
    result['throat_strip'] = throat_strip_checks(result['long_packet'],result['wide_packet'])
    result['energy_sweep'] = energy_sweep_checks(result['long_packet'])
    result['center_residual'] = center_residual_checks(result['energy_sweep'])
    result['spatial_residual'] = spatial_residual_checks(result['energy_sweep'])
    result['timing_response'] = timing_response_checks(result['energy_sweep'])
    return result


if __name__ == '__main__':
    result = run_checks()
    Path(__file__).with_name('mmp_pulse_results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: {k:v for k,v in value.items() if k!='plot'}
                      for key,value in result['loop_return']['cases'].items()}, indent=2))
