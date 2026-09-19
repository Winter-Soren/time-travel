"""Local ideal-cavity EM vacuum checks; run directly to reproduce the JSON.

One electromagnetic field, zero temperature, six perfectly conducting faces.
Units are hbar=c=1; each check states its cavity side lengths. Subtract only
the Minkowski bulk singularity at the center; do not use global energy/V.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json
import math

import mpmath as mp
import numpy as np


@lru_cache(maxsize=None)
def epstein_zeta(s, cutoff):
    """Z_3(s)=sum_{n!=0}|n|^(-2s), via the theta split at t=1.

    Return a truncated value and an analytic omitted-image bound. Neither
    value is interval arithmetic; cached special functions use 60 digits.
    """
    assert s >= 2 and cutoff >= 1
    with mp.workdps(60):
        s = mp.mpf(s)
        counts = Counter(i*i+j*j+k*k for i in range(-cutoff, cutoff+1)
                         for j in range(-cutoff, cutoff+1)
                         for k in range(-cutoff, cutoff+1) if i or j or k)
        def integral(K):
            return (K**(-s)*mp.gammainc(s, K, mp.inf)
                    +K**(s-mp.mpf('1.5'))*mp.gammainc(mp.mpf('1.5')-s, K, mp.inf))
        prefactor = mp.pi**s/mp.gamma(s)
        value = prefactor*(1/(s-mp.mpf('1.5'))-1/s
                          +mp.fsum(mult*integral(mp.pi*r2) for r2, mult in counts.items()))
        # theta(t)^3-theta_cutoff(t)^3 <= 8 exp[-pi(cutoff+1)^2 t], t>=1.
        omitted_bound = 8*prefactor*integral(mp.pi*(cutoff+1)**2)
        return value, omitted_bound


def center_images(tau=0, ewald_cutoff=4, odd_cutoff=13, series_order=20):
    """Minkowski-subtracted center density at imaginary time separation tau."""
    tau = mp.mpf(str(tau))
    assert 0 <= tau < 2 and odd_cutoff >= 1 and odd_cutoff % 2 == 1
    u = tau*tau/4
    order = series_order if tau else 0
    terms, bounds = [], []
    for k in range(order+1):
        zeta, zeta_bound = epstein_zeta(k+2, ewald_cutoff)
        coefficient = (k+1)*(2*k+1)*u**k/(16*mp.pi**2)
        terms.append(-(-1)**k*coefficient*zeta)
        bounds.append(coefficient*zeta_bound)
    translation = mp.fsum(terms)
    edge = -3*mp.pi/8*mp.fsum(
        mp.coth(mp.pi*mp.sqrt(p*p+q*q+tau*tau)/2)
        /(mp.sqrt(p*p+q*q+tau*tau)*mp.sinh(mp.pi*mp.sqrt(p*p+q*q+tau*tau)/2)**2)
        for p in range(-odd_cutoff, odd_cutoff+1, 2)
        for q in range(-odd_cutoff, odd_cutoff+1, 2))
    edge_bound = 15*mp.pi*mp.exp(-mp.pi*(odd_cutoff+1))/(1-mp.exp(-mp.pi))
    # Z_3(k+2)<=Z_3(2)<=26*zeta(2)=13*pi^2/3 by max-norm shells.
    M = order+1
    series_bound = mp.mpf(13)/48*u**M*((M+1)*(2*M+1)/(1-u)
        +(4*M+3)*u/(1-u)**2+2*u*(1+u)/(1-u)**3)
    return dict(value=translation+edge, translation=translation, edge=edge,
                ewald_bound=mp.fsum(bounds), edge_bound=edge_bound,
                series_bound=series_bound,
                analytic_truncation_bound=mp.fsum(bounds)+edge_bound+series_bound)


def transverse_modes(tau):
    """Independent physical transverse-mode sum, including one-zero modes.

    All-positive modes have normalization sqrt(8/V); one-zero modes sqrt(4/V).
    The polarization projector and center parity give W=2 omega(1-n_j^2/n^2)
    times 1/2 for a zero index. j is the unique odd or unique even index.
    Modes with three equal parities vanish at the center.
    """
    assert tau > 0
    cutoff = math.ceil(60/(math.pi*tau))
    ns = np.arange(cutoff+1, dtype=float)
    yy, zz = np.meshgrid(ns, ns, indexing='ij')
    yz2, odd_yz = yy*yy+zz*zz, yy % 2+zz % 2
    odd_squares = (yy % 2)*yy*yy+(zz % 2)*zz*zz
    zero_yz = (yy == 0).astype(int)+(zz == 0).astype(int)
    partial = []
    for ix in range(cutoff+1):
        r2, odd, zeros = yz2+ix*ix, odd_yz+ix % 2, zero_yz+(ix == 0)
        active = ((odd == 1) | (odd == 2)) & (zeros <= 1)
        denominator = np.where(active, r2, 1.)
        special = np.where(odd == 1, odd_squares+(ix % 2)*ix*ix,
                           r2-odd_squares-(ix % 2)*ix*ix)
        omega = math.pi*np.sqrt(denominator)
        weight = np.where(active, 2*omega*(1-special/denominator)
                          *np.where(zeros == 1, .5, 1.), 0.)
        partial.append(float(np.sum(weight*np.exp(-tau*omega))))
    raw = math.fsum(partial)
    minkowski = 3/(math.pi**2*tau**4)
    # Shell m has <=7m^2 nonnegative-index modes, each W<=2*pi*sqrt(3)*m.
    M, r = cutoff+1, mp.exp(-mp.pi*tau)
    tail = 14*mp.pi*mp.sqrt(3)*r**M*(M**3/(1-r)+3*M*M*r/(1-r)**2
        +3*M*r*(1+r)/(1-r)**3+r*(1+4*r+r*r)/(1-r)**4)
    return dict(cutoff=cutoff, raw=float(raw), minkowski=float(minkowski),
                value=float(raw-minkowski), analytic_mode_tail_bound=float(tail))


def exact_normalization_checks():
    """Rational identities; the factors hbar*c/pi^2 are kept outside."""
    # Parallel plates: -2*zeta(4)/(pi^2*(2L)^4), zeta(4)=pi^4/90.
    assert Fraction(-2, 16*90) == Fraction(-1, 720)
    # Two perpendicular reflecting planes: d=(2x,2y,0).
    assert -Fraction((2*2-1)*(2*2+11), 720) == -Fraction(1, 16)
    bilinears = [Fraction(0) for _ in range(3)]
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                if not (i or j or k):
                    continue
                m = (i, j, k)
                signs = tuple(1 if n % 2 == 0 else -1 for n in m)
                S, r2 = math.prod(signs), sum(n*n for n in m)
                for axis in range(3):
                    bilinears[axis] += ((1+S)*signs[axis]
                        *(-Fraction(1, r2*r2)+Fraction(2*m[axis]**2, r2**3)))
    rho = sum(bilinears)/2
    pressures = [rho-bilinear for bilinear in bilinears]
    assert pressures[0] == pressures[1] == pressures[2] == rho/3
    assert -rho+sum(pressures) == 0
    # At D=-rho=4A, four equal beams have n_z=-1/3 and transverse |n|^2=8/9.
    assert Fraction(8, 9)+Fraction(1, 9) == 1
    radiation_moments = [4*Fraction(3, 4), 2*Fraction(3, 4)*Fraction(8, 9),
                         2*Fraction(3, 4)*Fraction(8, 9), 4*Fraction(3, 4)*Fraction(1, 9)]
    assert radiation_moments == [3, Fraction(4, 3), Fraction(4, 3), Fraction(1, 3)]
    assert 4*Fraction(3, 4)*Fraction(-1, 3) == -1
    return dict(parallel_plates='rho=-pi^2*hbar*c/(720*L^4): exact rational coefficient checked',
        right_angle_wedge='rho=-hbar*c/(16*pi^2*r^4): exact match to published wedge formula at angle pi/2',
        single_reflection='E^2+B^2 cancels identically because 1+S=0',
        center_tensor='Static vacuum reflections kill flux and off-diagonal stress; cubic symmetry and zero bulk Maxwell trace give p_i=rho/3',
        rational_image_partial_sum_isotropic=True, rational_trace_zero=True,
        isotropy_scope='Finite symmetric image sums are exactly isotropic; absolute convergence preserves the identities in their limit')


def coherent_center_checks():
    """Exact one-event state in the same free cavity field; no persistent source."""
    import sympy as sp

    xyz = sp.symbols('x y z', real=True)
    center = dict.fromkeys(xyz, sp.Rational(1, 2))
    modes, values, curls = [], [], []
    for magnetic in (False, True):
        for axis in range(3):
            n, a = [2 if magnetic else 1]*3, sp.zeros(3, 1)
            n[axis] = 1 if magnetic else 0
            if magnetic:
                a[(axis+1) % 3], a[(axis+2) % 3] = -2, 2
            else:
                a[axis] = 2
            f = sp.Matrix([a[i]*sp.cos(sp.pi*n[i]*xyz[i])
                *sp.prod(sp.sin(sp.pi*n[j]*xyz[j]) for j in range(3) if j != i)
                for i in range(3)])
            curl = sp.Matrix([sp.diff(f[(i+2) % 3], xyz[(i+1) % 3])
                -sp.diff(f[(i+1) % 3], xyz[(i+2) % 3]) for i in range(3)])
            assert sp.simplify(sum(sp.diff(f[i], xyz[i]) for i in range(3))) == 0
            assert (sp.Matrix([sum(sp.diff(f[i], x, 2) for x in xyz)
                for i in range(3)])+sp.pi**2*sum(j*j for j in n)*f).applyfunc(sp.simplify) == sp.zeros(3, 1)
            assert sp.integrate(f.dot(f), *((x, 0, 1) for x in xyz)) == 1
            for face in range(3):
                for side in (0, 1):
                    assert all(sp.simplify(f[j].subs(xyz[face], side)) == 0
                        for j in range(3) if j != face)
                    assert sp.simplify(curl[face].subs(xyz[face], side)) == 0
            v, b = f.subs(center), curl.subs(center)
            unit = sp.eye(3)[:, axis]
            assert v == (sp.zeros(3, 1) if magnetic else 2*unit)
            assert b == (8*sp.pi*unit if magnetic else sp.zeros(3, 1))
            values.append(v)
            curls.append(b)
            modes.append(dict(family='magnetic' if magnetic else 'electric', indices=n,
                vector_amplitudes=list(map(str, a)), omega=str(sp.pi*sp.sqrt(sum(j*j for j in n))),
                center_f=list(map(str, v)), center_curl=list(map(str, b)),
                divergence_helmholtz_norm_and_all_six_faces_exact=True))

    A = sp.symbols('A', positive=True)
    e, b = sp.symbols('electric_factor magnetic_factor', real=True)
    directions = [sp.Matrix([sg*2*sp.sqrt(2)/3, 0, -sp.Rational(1, 3)]) for sg in (1, -1)]
    directions += [sp.Matrix([0, sg*2*sp.sqrt(2)/3, -sp.Rational(1, 3)]) for sg in (1, -1)]
    radiation, component_energies = sp.zeros(4), []
    for j, direction in enumerate(directions):
        E = sp.sqrt(3*A)*sp.eye(3)[:, 1 if j < 2 else 0]
        B = direction.cross(E)
        assert sp.Matrix.hstack(*values[:3])*(-E/2) == -E
        assert (sp.Matrix.hstack(*curls[3:])*(B/(8*sp.pi))-B).applyfunc(sp.simplify) == sp.zeros(3, 1)
        assert sp.simplify(direction.dot(direction)-1) == 0
        assert sp.simplify(E.dot(B)) == sp.simplify(E.dot(E)-B.dot(B)) == 0
        classical_energy = sp.simplify((E/2).dot(E/2)/2+(3*sp.pi)**2*(B/(8*sp.pi)).dot(B/(8*sp.pi))/2)
        assert classical_energy == 75*A/128
        component_energies.append(str(classical_energy/A))
        rho, flux = sp.simplify((e*e*E.dot(E)+b*b*B.dot(B))/2), e*b*E.cross(B)
        stress = sp.zeros(4)
        stress[0, 0], stress[0, 1:4], stress[1:4, 0] = rho, flux.T, flux
        stress[1:4, 1:4] = rho*sp.eye(3)-e*e*E*E.T-b*b*B*B.T
        k = sp.Matrix([1, *direction])
        assert (stress.subs({e: 1, b: 1})-3*A*k*k.T).applyfunc(sp.simplify) == sp.zeros(4)
        radiation += stress/4
    vacuum = sp.diag(-4*A, -4*A/3, -4*A/3, -4*A/3)
    total = (vacuum+radiation).applyfunc(sp.simplify)
    ell = sp.Matrix([1, 0, 0, 1])
    assert (total.subs({e: 1, b: 1})+A*ell*ell.T).applyfunc(sp.simplify) == sp.zeros(4)
    assert sp.simplify(total[0, 0]/A+4-sp.Rational(3, 2)*(e*e+b*b)) == 0
    assert sp.simplify(total[0, 3]/A+e*b) == 0
    assert (sp.sqrt(2)/3).is_irrational is True
    samples = []
    for tau in (0., .05, .1, .2):
        ce, cb = math.cos(math.pi*math.sqrt(2)*tau), math.cos(3*math.pi*tau)
        samples.append(dict(tau=tau, flux_ratio=ce*cb, total_rho_over_A=-4+1.5*(ce*ce+cb*cb)))

    # Parity is necessary for a nonzero center amplitude, even when a polarization vanishes.
    integers = sp.symbols('k0:3', integer=True)
    parity = []
    for bits in np.ndindex(2, 2, 2):
        electric = [i for i in range(3) if bits[i] == 0 and all(bits[j] == 1 for j in range(3) if j != i)]
        magnetic = [i for i in range(3) if bits[i] == 1 and all(bits[j] == 0 for j in range(3) if j != i)]
        assert sp.expand(sum((2*k+p)**2 for k, p in zip(integers, bits))-sum(bits)
            -4*sum(k*(k+p) for k, p in zip(integers, bits))) == 0
        assert not (electric and magnetic)
        assert not electric or sum(bits) == 2
        assert not magnetic or sum(bits) == 1
        parity.append(dict(bits=list(bits), squared_frequency_over_pi2_mod4=sum(bits) % 4,
            possible_electric_axes=electric, possible_magnetic_axes=magnetic))
    return dict(
        convention='Contravariant orthonormal stress, signature -+++; hbar=c=cube side=1; center at x=y=z=1/2, reference time t0',
        six_modes=modes,
        mode_definition='f_i=a_i*cos(pi*n_i*x_i)*product_{j!=i}sin(pi*n_j*x_j); A_field=sum Q_lambda*f_lambda, E=-sum P_lambda*f_lambda, B=sum Q_lambda*curl(f_lambda)',
        exact_checks='Six divergence, vector Helmholtz, unit spatial norm, tangential f and normal curl on all six faces, center values, four null tensors and mixture equality passed',
        classical_initial_data='For any prescribed center E,B: electric-mode P_i=-E_i/2, Q_i=0; magnetic-mode Q_i=B_i/(8*pi), P_i=0. Finite smooth modes solve source-free Maxwell equations and the conducting boundary conditions',
        quantum_state='Use alpha_lambda=(sqrt(omega_lambda)*Q_lambda+i*P_lambda/sqrt(omega_lambda))/sqrt(2) to displace the same cavity vacuum. D^dagger*F*D=F+F_classical gives <T>_alpha,ren=<T>_vac,ren+T[F_classical]; the finite displacement leaves the ultraviolet singularity unchanged',
        mixed_state='density=(1/4)*sum_j |alpha_j><alpha_j|. Choose E_j=sqrt(3*A)*ey for the first two directions, sqrt(3*A)*ex for the other two; B_j=n_j cross E_j. Linearity of a probabilistic mixture removes inter-state interference; this is not a coherent superposition of the four states',
        directions=[list(map(str, n)) for n in directions], component_excitation_energy_over_A=component_energies,
        radiation_tensor_over_A=[[str(v) for v in row] for row in (radiation.subs({e: 1, b: 1})/A).tolist()],
        total_tensor_over_A=[[str(v) for v in row] for row in (-ell*ell.T).tolist()],
        match='At A=-rho_vac/4=-q/4>0, vacuum=diag(-4*A,-4*A/3,-4*A/3,-4*A/3); mixture gives exactly -A*ell^a*ell^b for ell=(1,0,0,1), at one event',
        witness_evolution=dict(time='Real evolution time tau=c*(t-t0)/L, distinct from the imaginary-time regulator', fields='E_center=E0*cos(pi*sqrt(2)*tau); B_center=B0*cos(3*pi*tau)',
            flux_ratio='T0z(tau)/T0z(0)=cos(pi*sqrt(2)*tau)*cos(3*pi*tau)',
            total_rho_over_A='-4+(3/2)*(cos(pi*sqrt(2)*tau)^2+cos(3*pi*tau)^2)', samples=samples,
            no_exact_nonzero_recurrence='rho/A=-1 requires both cosines squared=1, hence sqrt(2)*tau and 3*tau integers. Irrational sqrt(2)/3 permits only tau=0. This applies to this six-mode witness and this fixed target amplitude, not every cavity state or an approximate tolerance'),
        stationary_center_audit=dict(parity=parity,
            identity='n_i=2*k_i+b_i implies |n|^2=sum(b_i)+4*sum(k_i*(k_i+b_i)); all eight bit patterns checked exactly',
            disjoint_spectra='Nonzero center E requires two odd indices, so |n|^2=2 mod4; nonzero center B requires one odd index, so |n|^2=1 mod4. Their frequencies cannot coincide',
            argument='For a stationary density operator commuting with the free cavity Hamiltonian, expectation(a_lambda^dagger*a_mu)=0 unless omega_lambda=omega_mu, while expectation(a_lambda*a_mu)=0 because every physical frequency is positive. The aa^dagger commutator term is diagonal in the mode and no mode contributes both center E and B. Thus every symmetrized center E-B correlation is zero, including temporal point splitting; its Minkowski subtraction is also zero, so all T0i(center)=0',
            conclusion='No stationary state with defined renormalized stress in this free empty lossless closed ideal cube can supply the target nonzero center flux, even if the state breaks spatial symmetry. The nonstationary one-event mixed coherent state is allowed',
            scope='External drives, ports, changed shape, realistic walls or curved geometry change the hypotheses. This supplies no universal finite-duration or beat-time bound: higher-frequency gaps can approach zero'),
        sources=['https://doi.org/10.1103/PhysRev.131.2766 (Glauber: coherent radiation states and factorization; specific cavity construction and parity audit derived here)',
                 'https://arxiv.org/pdf/0804.3916 (Lim-Teo: free physical conducting-cavity modes)'],
        scope='Same electromagnetic field and boundary-compatible quantum state at one event; no independent-field additivity assumption. Null-ray tensors describe only that event, not global plane beams. No persistent finite region, stress-fluctuation control, wall/material validation, self-consistent gravity, stability or traversal has been established')


def rectangular_resonance_checks():
    """Stationary positive center radiation in a resonant rectangular cavity."""
    import sympy as sp

    xyz = sp.symbols('x y z', real=True)
    Lx, Ly, Lz = sp.symbols('Lx Ly Lz', positive=True)
    nE, nB = (2, 1, 3), (2, 3, 2)
    frequency_difference = sp.pi**2*sum((nE[i]**2-nB[i]**2)/L**2
                                       for i, L in enumerate((Lx, Ly, Lz)))
    assert sp.simplify(frequency_difference-sp.pi**2*(5/Lz**2-8/Ly**2)) == 0
    assert sp.simplify(frequency_difference.subs(Lz, sp.sqrt(sp.Rational(5, 8))*Ly)) == 0
    lengths = (sp.Integer(1), sp.Integer(1), sp.sqrt(10)/4)
    volume = sp.prod(lengths)
    normalization = sp.sqrt(8/volume)
    center = dict(zip(xyz, (L/2 for L in lengths)))
    polarizations = (sp.Matrix([1, -2, 0])/sp.sqrt(5),
                     sp.Matrix([-2*sp.sqrt(2), 0, sp.sqrt(5)])/sp.sqrt(13))
    values, curls, modes = [], [], []
    for indices, polarization in zip((nE, nB), polarizations):
        wavevector = sp.Matrix([sp.pi*indices[i]/lengths[i] for i in range(3)])
        assert sp.simplify(polarization.dot(polarization)-1) == 0
        assert sp.simplify(wavevector.dot(polarization)) == 0
        omega2 = sp.simplify(wavevector.dot(wavevector))
        assert omega2 == 97*sp.pi**2/5
        f = sp.Matrix([normalization*polarization[i]*sp.cos(wavevector[i]*xyz[i])
            *sp.prod(sp.sin(wavevector[j]*xyz[j]) for j in range(3) if j != i)
            for i in range(3)])
        curl = sp.Matrix([sp.diff(f[(i+2) % 3], xyz[(i+1) % 3])
            -sp.diff(f[(i+1) % 3], xyz[(i+2) % 3]) for i in range(3)])
        assert sp.simplify(sum(sp.diff(f[i], xyz[i]) for i in range(3))) == 0
        assert all(sp.simplify(sum(sp.diff(f[i], coordinate, 2) for coordinate in xyz)
                               +omega2*f[i]) == 0 for i in range(3))
        norm_squared = 0
        for i in range(3):
            factors = [sp.cos(wavevector[j]*xyz[j]) if j == i
                       else sp.sin(wavevector[j]*xyz[j]) for j in range(3)]
            norm_squared += normalization**2*polarization[i]**2*sp.prod(
                sp.integrate(factors[j]**2, (xyz[j], 0, lengths[j])) for j in range(3))
        assert sp.simplify(norm_squared-1) == 0
        for face in range(3):
            for side in (0, lengths[face]):
                assert all(sp.simplify(f[i].subs(xyz[face], side)) == 0 for i in range(3) if i != face)
                assert sp.simplify(curl[face].subs(xyz[face], side)) == 0
        values.append(f)
        curls.append(curl)
        modes.append(dict(indices=list(indices), polarization=list(map(str, polarization)),
            omega_squared=str(omega2), center_f=list(map(str, f.subs(center).applyfunc(sp.simplify))),
            center_curl=list(map(str, curl.subs(center).applyfunc(sp.simplify))),
            transverse_unit_norm_helmholtz_and_six_PEC_faces_exact=True))
    assert sp.integrate(values[0].dot(values[1]),
                        *((xyz[i], 0, lengths[i]) for i in range(3))) == 0
    t, theta, E0, B0 = sp.symbols('t theta E0 B0', real=True)
    omega, K = sp.pi*sp.sqrt(sp.Rational(97, 5)), sp.pi*sp.sqrt(sp.Rational(52, 5))
    phase = omega*t+theta
    qE = -E0*sp.sqrt(5)/(omega*normalization)*sp.sin(phase)
    qB = B0/(normalization*K)*sp.cos(phase)
    E = -sp.diff(qE, t)*values[0]-sp.diff(qB, t)*values[1]
    B = qE*curls[0]+qB*curls[1]
    Ec, Bc = E.subs(center).applyfunc(sp.simplify), B.subs(center).applyfunc(sp.simplify)
    assert (Ec-sp.Matrix([E0*sp.cos(phase), 0, 0])).applyfunc(sp.simplify) == sp.zeros(3, 1)
    assert (Bc-sp.Matrix([0, B0*sp.cos(phase), 0])).applyfunc(sp.simplify) == sp.zeros(3, 1)
    for face in range(3):
        for side in (0, lengths[face]):
            Ef = E.subs(xyz[face], side).applyfunc(sp.simplify)
            Bf = B.subs(xyz[face], side).applyfunc(sp.simplify)
            assert sp.simplify(Ef.cross(Bf)[face]) == 0
    rho, flux = (Ec.dot(Ec)+Bc.dot(Bc))/2, Ec.cross(Bc)
    stress = sp.zeros(4)
    stress[0, 0], stress[0, 1:4], stress[1:4, 0] = rho, flux.T, flux
    stress[1:4, 1:4] = rho*sp.eye(3)-Ec*Ec.T-Bc*Bc.T
    averaged = stress.subs(t, 0).applyfunc(lambda value:
        sp.simplify(sp.integrate(value, (theta, 0, 2*sp.pi))/(2*sp.pi)))
    expected = sp.Matrix([[(E0**2+B0**2)/4, 0, 0, E0*B0/2],
        [0, (B0**2-E0**2)/4, 0, 0], [0, 0, (E0**2-B0**2)/4, 0],
        [E0*B0/2, 0, 0, (E0**2+B0**2)/4]])
    assert averaged == expected and sp.trace(sp.diag(-1, 1, 1, 1)*averaged) == 0
    for direction in (-1, 1):
        ell = sp.Matrix([1, 0, 0, direction])
        assert averaged.subs(B0, direction*E0) == E0**2/2*ell*ell.T
    excitation = sp.simplify((sp.diff(qE, t)**2+sp.diff(qB, t)**2
                              +omega**2*(qE*qE+qB*qB))/2)
    assert sp.simplify(excitation-volume*(5*E0**2+sp.Rational(97, 52)*B0**2)/16) == 0

    # Conditional full center tensor; D and p are not numerical vacuum inputs.
    D, A = sp.symbols('D A', positive=True)
    p = sp.symbols('p', negative=True)
    theta_z = sp.symbols('theta_z', real=True)
    fz = sp.Matrix([0, 0, 2/sp.sqrt(volume)*sp.sin(sp.pi*xyz[0])*sp.sin(sp.pi*xyz[1])])
    curl_z = sp.Matrix([sp.diff(fz[(i+2) % 3], xyz[(i+1) % 3])
                       -sp.diff(fz[(i+1) % 3], xyz[(i+2) % 3]) for i in range(3)])
    omega_z = sp.pi*sp.sqrt(2)
    assert sum(sp.diff(fz[i], xyz[i]) for i in range(3)) == 0
    assert all(sp.simplify(sum(sp.diff(fz[i], coordinate, 2) for coordinate in xyz)
                           +omega_z**2*fz[i]) == 0 for i in range(3))
    assert sp.integrate(fz.dot(fz), *((xyz[i], 0, lengths[i]) for i in range(3))) == 1
    for previous_mode in values:
        assert sp.integrate(fz.dot(previous_mode),
                            *((xyz[i], 0, lengths[i]) for i in range(3))) == 0
    for face in range(3):
        for side in (0, lengths[face]):
            assert all(sp.simplify(fz[i].subs(xyz[face], side)) == 0 for i in range(3) if i != face)
            assert sp.simplify(curl_z[face].subs(xyz[face], side)) == 0
    qz = -sp.sqrt(volume)*sp.sqrt(-4*p)*sp.sin(omega_z*t+theta_z)/(2*omega_z)
    Ez = (-sp.diff(qz, t)*fz).subs(center).applyfunc(sp.simplify)
    assert (Ez-sp.Matrix([0, 0, sp.sqrt(-4*p)*sp.cos(omega_z*t+theta_z)]))\
        .applyfunc(sp.simplify) == sp.zeros(3, 1)
    assert curl_z.subs(center) == sp.zeros(3, 1)
    mean_Ez2 = sp.simplify(sp.integrate(Ez[2].subs(t, 0)**2, (theta_z, 0, 2*sp.pi))/(2*sp.pi))
    assert mean_Ez2 == -2*p
    z_stress = mean_Ez2/2*sp.diag(1, 1, 1, -1)
    assert z_stress == sp.diag(-p, -p, -p, p)
    z_energy = sp.simplify((sp.diff(qz, t)**2+omega_z**2*qz*qz)/2)
    assert z_energy == -volume*p/2
    # Independent uniform phases have zero first moments, so cross terms vanish
    # when the coherent amplitudes are added in the SAME three-mode space.
    assert sp.integrate(sp.cos(theta_z), (theta_z, 0, 2*sp.pi)) == 0
    maximum_A = (D+p)/2
    w_minus, w_plus = maximum_A, maximum_A-A
    ell_minus, ell_plus = sp.Matrix([1, 0, 0, -1]), sp.Matrix([1, 0, 0, 1])
    conditional_vacuum = sp.diag(-D, p, p, -D-2*p)
    added = w_minus*ell_minus*ell_minus.T+w_plus*ell_plus*ell_plus.T+z_stress
    total = (conditional_vacuum+added).applyfunc(sp.simplify)
    assert total == -A*ell_plus*ell_plus.T
    Cxy, Cz = D-A+p, -2*p
    assert sp.simplify(added[0, 0]-(2*Cxy+Cz)/2) == 0
    assert added[0, 3] == -A
    assert sp.expand(Cxy**2-A**2-(D+p)*(D+p-2*A)) == 0
    assert sp.simplify(Cxy.subs(A, maximum_A)-maximum_A) == 0
    w = sp.symbols('w', nonnegative=True)
    ray_energy = sp.simplify(excitation.subs({E0: sp.sqrt(2*w), B0: sp.sqrt(2*w)}))
    assert ray_energy == 357*volume*w/416
    total_energy = sp.simplify(ray_energy.subs(w, w_minus+w_plus)+z_energy)
    assert sp.simplify(total_energy-volume*(357*(D-A)+149*p)/416) == 0
    return dict(convention='Contravariant orthonormal added EM stress, signature -+++; hbar=c=1; physical side lengths (1,1,sqrt(10)/4)',
        general_resonance='For nE=(2,1,3), nB=(2,3,2), omegaE^2-omegaB^2=pi^2*(5/Lz^2-8/Ly^2); Lz^2=5*Ly^2/8 removes the difference for any Lx>0',
        lengths_exact=list(map(str, lengths)), volume_exact=str(volume), modes=modes,
        mode_definition='u_i=N*alpha_i*cos(k_i*x_i)*product_{j!=i}sin(k_j*x_j), N=sqrt(8/V), k_i=pi*n_i/L_i',
        two_modes_orthonormal_exact=True,
        center_amplitudes='uE(center)=(N/sqrt(5),0,0), curl(uE)=0; uB(center)=0, curl(uB)=(0,N*K,0); K=pi*sqrt(52/5)',
        classical_solution='A_field=qE*uE+qB*uB; qE=-E0*sqrt(5)*sin(omega*t+theta)/(omega*N), qB=B0*cos(omega*t+theta)/(N*K), omega=pi*sqrt(97/5)',
        center_fields='E=(E0*cos(omega*t+theta),0,0); B=(0,B0*cos(omega*t+theta),0)',
        stationary_quantum_state='Displace the SAME rectangular-cavity vacuum with these two finite mode coordinates and momenta, alpha_j=(sqrt(omega)*q_j+i*p_j/sqrt(omega))/sqrt(2). Uniformly mix theta over0..2*pi. Equal frequencies make free evolution a common theta shift, leaving this density operator invariant; no external drive is introduced. Normal coherence between the degenerate modes survives phase averaging',
        renormalization='For each coherent component the renormalized stress is rectangular-vacuum stress plus its classical stress. The finite displacement leaves the ultraviolet singularity unchanged; averaging the density operator averages the classical tensor',
        phase_averaged_added_tensor=[[str(value) for value in row] for row in averaged.tolist()],
        positive_null_case='E0=B0 gives added tensor=(E0^2/2)*ell_plus*ell_plus, ell_plus=(1,0,0,1). B0=-E0 reverses flux while preserving positive energy; both identities checked exactly',
        excitation_energy_exact=str(excitation), excitation_energy_formula='V*(5*E0^2+(97/52)*B0^2)/16, finite for finite E0,B0; mean added photon number=excitation_energy/omega',
        wall_flux='Normal E cross B is identically zero on all six PEC faces, checked for the full two-mode solution. Nonzero center flux does not imply net transport through the closed walls',
        exact_checks='Unit polarizations, transversality, both vector Helmholtz equations, unit norms and mutual orthogonality, six-face PEC conditions, center fields, all16 averaged tensor components, zero trace, both null directions, conserved finite excitation and zero wall-normal flux',
        sources=['https://arxiv.org/pdf/0804.3916 (Lim-Teo: physical common cavity modes; this resonance is our derivation)',
                 'https://doi.org/10.1103/PhysRev.131.2766 (Glauber: coherent radiation states; this common-phase mixture is our construction)'],
        conditional_compensation=dict(
            assumptions='Conditional vacuum tensor diag(-D,p,p,-D-2p), D>0, p<0, D+p>0; no numerical rectangular vacuum value is asserted or borrowed from the cube',
            target='-A*ell_plus^a*ell_plus^b, ell_plus=(1,0,0,1), 0<A<=(D+p)/2',
            added_weights='w_minus=(D+p)/2, w_plus=(D+p)/2-A; both nonnegative on the stated interval. For each resonant family set E0=sqrt(2*w), B0=+/-sqrt(2*w) for the corresponding flux direction',
            extra_mode=dict(indices=[1, 1, 0], definition='fz=(0,0,2*sin(pi*x)*sin(pi*y)/sqrt(V)), omega_z=pi*sqrt(2)',
                coordinate='qz=-sqrt(V)*sqrt(-4*p)*sin(omega_z*t+theta_z)/(2*omega_z), pz=dqz/dt',
                center_fields='Ez=sqrt(-4*p)*cos(omega_z*t+theta_z), B=0',
                mean_Ez_squared=str(mean_Ez2), phase_averaged_tensor=[[str(v) for v in row] for row in z_stress.tolist()],
                excitation_energy=str(z_energy), exact_checks='Unit norm, orthogonality to both resonant modes, divergence, Helmholtz, all six PEC faces, center fields, phase-averaged stress and conserved finite excitation'),
            added_tensor=[[str(v) for v in row] for row in added.applyfunc(sp.simplify).tolist()],
            total_tensor=[[str(v) for v in row] for row in total.tolist()], full_target_identity_exact=True,
            state='Integrate |alpha_minus(theta_minus)+alpha_plus(theta_plus)+alpha_z(theta_z)><same| against the product of three uniform phase measures. All amplitudes belong to the SAME cavity mode space; the +/- families share the two resonant modes. This is a mixture of coherent states of summed amplitudes, explicitly not a Hilbert-space tensor product for those families',
            stationarity='Each coherent family evolves by its phase shift; the joint uniform measure is invariant under free evolution. Zero first phase moments remove cross terms. The full center stress EXPECTATION is stationary at all times; individual coherent realizations need not have a stationary tensor. No drive or independent electromagnetic field is introduced',
            finite_excitation='Each ray family has mean energy357*V*w/416 and the extra mode has energy -V*p/2. Their sum is finite and positive under the stated assumptions; all three physical mode frequencies are finite and positive',
            total_mean_excitation_energy=str(total_energy),
            necessity=dict(class_of_states='Arbitrary positive mixtures of coherent displacements of this fixed vacuum, with finite classical second moments',
                required_moments='Cx=Cy=D-A+p, Cz=-2p, where Ci=<Ei^2+Bi^2>_classical>=0',
                inequality='Let X=(Ex,Bx), Y=(By,-Ey). Cauchy-Schwarz in the classical ensemble gives A^2=<X dot Y>^2<=Cx*Cy=(D-A+p)^2. Since Cx>=0, 2*A<=D+p',
                factorization='(D-A+p)^2-A^2=(D+p)*(D+p-2*A), checked exactly; the construction attains A=(D+p)/2',
                exclusions='Not a bound for arbitrary squeezed or other nonclassical states: their vacuum-subtracted additions need not be positive classical second moments'),
            previous_bound='D/4 was the null-beam-only compensation bound for the isotropic cube vacuum. The extra positive-energy electric mode supplies longitudinal Maxwell tension, so that bound does not restrict this broader coherent-mixture class. Formal isotropic input p=-D/3 gives D/3 here, but the actual cube retains its separate stationary parity obstruction',
            scope='Conditional stationary full tensor at the center only. Neither its vacuum hypotheses nor a finite-region stress match, wall/material realization, self-consistent gravitational source, fluctuations, stability or traversal are established by this check'),
        scope='Stationary positive added center radiation and a conditional full center-tensor construction in a flat ideal rectangle. Vacuum parameters D,p remain conditional and are not replaced by the cube tensor. No finite-region source, material validation, self-consistent gravitational solution, stability or traversal is established')


def rectangle_vacuum_checks():
    """Independent local vacuum integrals and physical modes; ideal static rectangle."""
    with mp.workdps(40):
        squared_sides = (mp.mpf(1), mp.mpf(1), mp.mpf(5)/8)
        theta_cutoff = 8

        @lru_cache(maxsize=None)
        def theta_moment(u, side_squared, bit):
            v = u*side_squared
            if v >= 1:
                terms = [(-1)**(bit*n)*mp.exp(-mp.pi*v*n*n)
                         for n in range(1, theta_cutoff+1)]
                delta = 2*mp.fsum(terms)
                moment = 2*mp.fsum(n*n*t for n, t in enumerate(terms, 1))
                return 1+delta, delta, moment
            shift = mp.mpf('.5') if bit else mp.mpf(0)
            indices = range(theta_cutoff+1) if bit else range(1, theta_cutoff+1)
            terms = [(n+shift, mp.exp(-mp.pi*(n+shift)**2/v)) for n in indices]
            S0 = (0 if bit else 1)+2*mp.fsum(t for n, t in terms)
            S2 = 2*mp.fsum(n*n*t for n, t in terms)
            theta = S0/mp.sqrt(v)
            moment = S0/(2*mp.pi*v**mp.mpf('1.5'))-S2/v**mp.mpf('2.5')
            return theta, theta-1, moment

        @lru_cache(maxsize=None)
        def image_sums(u, axis, geometry):
            total, moment = mp.mpf(0), mp.mpf(0)
            for complement in (False, True):
                bits = [int((i == axis) != complement) for i in range(3)]
                rows = [theta_moment(u, geometry[i], bits[i]) for i in range(3)]
                d = [row[1] for row in rows]
                # Product(theta)-1 without cancellation when all theta approach 1.
                total += sum(d)+d[0]*d[1]+d[0]*d[2]+d[1]*d[2]+mp.fprod(d)
                moment += rows[axis][2]*mp.fprod(rows[i][0] for i in range(3) if i != axis)
            return total, moment

        def image_stress(tau, geometry=squared_sides):
            C, estimated_errors = [], []
            for axis in (0, 2):
                def integrand(u):
                    if not u:
                        return mp.mpf(0)
                    S, Q = image_sums(u, axis, geometry)
                    v = mp.pi*u*tau*tau
                    return mp.exp(-v)*(u*(v-1)*S+mp.pi*u*u*geometry[axis]*Q)
                value, estimate = mp.quad(integrand, [0, mp.mpf('.25'), mp.mpf('.5'), 1, 2, 4, 8, 16, mp.inf], error=True)
                C.append(value)
                estimated_errors.append(mp.nstr(estimate, 8))
            rho, px, pz = C[0]+C[1]/2, C[1]/2, C[0]-C[1]/2
            assert abs(-rho+2*px+pz) < mp.mpf('1e-30')
            return rho, px, pz, estimated_errors

        def physical_modes(tau):
            sides = (1., 1., math.sqrt(5/8))
            volume, smallest = math.prod(sides), min(sides)
            cutoff = math.ceil(60/(math.pi*tau))
            ns = np.arange(cutoff+1, dtype=float)
            y, z = np.meshgrid(ns, ns, indexing='ij')
            y2, z2 = (y/sides[1])**2, (z/sides[2])**2
            by, bz = y % 2, z % 2
            zero_yz = (y == 0).astype(int)+(z == 0).astype(int)
            parts = [[], [], []]
            for ix in range(cutoff+1):
                bx, x2 = ix % 2, (ix/sides[0])**2
                r2, odd, zeros = x2+y2+z2, bx+by+bz, zero_yz+(ix == 0)
                active = ((odd == 1) | (odd == 2)) & (zeros <= 1)
                denominator = np.where(active, r2, 1.)
                odd_square = bx*x2+by*y2+bz*z2
                special = np.where(odd == 1, odd_square, r2-odd_square)
                omega = math.pi*np.sqrt(denominator)
                C = np.where(active, 4*omega/volume*(1-special/denominator)
                             *np.where(zeros == 1, .5, 1.), 0.)*np.exp(-tau*omega)
                parity_needed = np.where(odd == 1, 1, 0)
                for i, bit in enumerate((bx, by, bz)):
                    parts[i].append(float(np.sum(np.where(bit == parity_needed, C, 0.))))
            raw = np.array([math.fsum(part) for part in parts])
            minkowski_C = 2/(math.pi**2*tau**4)
            bilinear = raw-minkowski_C
            rho = float(np.sum(bilinear)/2)
            M, r = cutoff+1, mp.exp(-mp.pi*tau)
            tail = 14*mp.pi*mp.sqrt(3)/(volume*smallest)*r**M*(
                M**3/(1-r)+3*M*M*r/(1-r)**2+3*M*r*(1+r)/(1-r)**3+r*(1+4*r+r*r)/(1-r)**4)
            return dict(cutoff=cutoff, rho=rho, pressures=[float(rho-c) for c in bilinear],
                        raw_rho=float(np.sum(raw)/2), minkowski_rho=1.5*minkowski_C,
                        analytic_mode_tail_bound=float(tail))

        rho, p, pz, estimates = image_stress(mp.mpf(0))
        D = -rho
        assert rho < 0 and p < 0 and pz < 0 and D+p > 0
        comparisons = []
        for tau in (.4, .15, .08):
            reference = image_stress(mp.mpf(str(tau)))
            modes = physical_modes(tau)
            expected = [float(reference[0]), float(reference[1]), float(reference[1]), float(reference[2])]
            observed = [modes['rho'], *modes['pressures']]
            error = max(abs(a-b) for a, b in zip(expected, observed))
            acceptance = max(1e-11, 8*np.finfo(float).eps*(modes['raw_rho']+modes['minkowski_rho']))
            assert math.isfinite(error) and error <= acceptance
            comparisons.append(dict(tau=tau, **modes, image_rho=expected[0],
                image_pressures=expected[1:], maximum_stress_disagreement=error,
                floating_acceptance=acceptance, quadrature_error_estimates=reference[3]))
        cube_checks = []
        for tau in map(mp.mpf, ('0', '.15')):
            crho, cpx, cpz, errors = image_stress(tau, (mp.mpf(1),)*3)
            reference = center_images(tau, ewald_cutoff=5, odd_cutoff=17)
            difference = abs(crho-reference['value'])
            isotropy = max(abs(cpx-cpz), abs(crho/3-cpx))
            assert difference <= reference['analytic_truncation_bound']+mp.mpf('1e-30')
            assert isotropy < mp.mpf('1e-30')
            cube_checks.append(dict(tau=float(tau), rho=mp.nstr(crho, 30),
                pressure=mp.nstr(cpx, 30), existing_cube_image_rho=mp.nstr(reference['value'], 30),
                disagreement=float(difference), existing_cube_analytic_tail=float(reference['analytic_truncation_bound']),
                isotropy_residual=float(isotropy), quadrature_error_estimates=errors))
        audited_cutoff = 2501
        c = mp.sqrt(squared_sides[2])
        edge_bound = mp.mpf('1.001')*(1/c+2/c**3)/(mp.pi*(audited_cutoff+1))
        return dict(
            model='Zero-temperature vacuum of one free EM field in a flat static conducting rectangle; local center only',
            normalization='Lx=Ly=L, Lz=sqrt(5/8)*L; reported stress coefficients multiply hbar*c/L^4. tau is imaginary-time separation in units L/c',
            rho_coefficient=float(rho), rho_decimal=mp.nstr(rho, 30),
            D_coefficient=float(D), D_decimal=mp.nstr(D, 30),
            transverse_pressure_coefficient=float(p), transverse_pressure_decimal=mp.nstr(p, 30),
            axial_pressure_coefficient=float(pz), axial_pressure_decimal=mp.nstr(pz, 30),
            sign_checks=dict(negative_rho=True, negative_transverse_pressure=True, negative_axial_pressure=True, D_plus_p_positive=True),
            trace_residual=float(-rho+2*p+pz),
            conditional_compensation_numeric_link=dict(A_max_coefficient=float((D+p)/2), A_max_decimal=mp.nstr((D+p)/2, 30),
                relation='Amax=(D+p)/2 for the separately checked rectangular_resonance.conditional_compensation coherent-mixture theorem; coefficient multiplies hbar*c/L^4 and applies at the center only'),
            unit_cube_limit_checks=cube_checks,
            image_kernel='d_i=m_i*Li, s_i=(-1)^m_i, S=product(s_i), U=sum(d_i^2)+tau^2: C_i=sum_{m!=0}(1+S)*s_i/pi^2*(-1/U^2+2*(d_i^2+tau^2)/U^3); rho=sum(C_i)/2, p_i=rho-C_i',
            theta_definition='theta_b(v)=sum_n(-1)^(b*n)*exp(-pi*v*n^2); Theta_b(u)=product_j theta_bj(u*Lj^2); M_i,b=sum_m m_i^2*(-1)^(b dot m)*exp(-pi*u*sum_j Lj^2*m_j^2)',
            theta_integral='S_i=Theta_ei+Theta_(1-ei)-2, Q_i=M_i,ei+M_i,(1-ei). C_i(tau)=integral_0^infinity exp(-pi*u*tau^2)*[u*(pi*u*tau^2-1)*S_i+pi*u^2*Li^2*Q_i] du. The -2 removes only the zero Minkowski image',
            modular_identity='For v<1, theta_b(v)=v^(-1/2)*sum_n exp[-pi*(n+b/2)^2/v]; its n^2 moment is S0/(2*pi*v^(3/2))-S2/v^(5/2). Direct signed Gaussian sums are used for v>=1',
            physical_mode_weights='R^2=sum_i(n_i/Li)^2, omega=pi*sqrt(R^2), V=product(Li), n_i>=0 with at most one zero. For mixed parity j is the unique odd or unique even index; W_rho=2*omega/V*(1-(n_j/Lj)^2/R^2)*2^(-zero_count). Only C_j=2*W_rho is nonzero. Multiply by exp(-tau*omega); other parities vanish',
            mode_subtraction='Subtract 3/(pi^2*tau^4) from rho and 2/(pi^2*tau^4) from each C_i; compute each pressure as rho-C_i',
            finite_tau_checks=comparisons,
            maximum_mode_image_disagreement=max(row['maximum_stress_disagreement'] for row in comparisons),
            precision=dict(decimal_working_digits=40, reported_digits=30, theta_cutoff=theta_cutoff,
                zero_tau_quadrature_error_estimates=estimates,
                mode_tail='Max-index shell m has <=7*m^2 modes and W_rho<=2*pi*sqrt(3)*m/(V*Lmin). Each mode has only C_j=2*W_rho, so every pressure is +/-W_rho; the positive density tail also bounds each absolute pressure tail. This fixed aspect has Lmax=1',
                scope='Quadrature estimates and refinement evidence are not interval certification. Mode sums are float64; analytic mode tails exclude arithmetic rounding. The direct/modular theta split avoids the algebraic edge tail'),
            separate_edge_sign_audit=dict(recomputed_in_default_checker=False, odd_cutoff=audited_cutoff,
                reported_partial_transverse_pressure=-.001280973317475672,
                complete_hyperbolic_bilinear_tail_bound=float(edge_bound),
                reported_partial_plus_tail=float(mp.mpf('-.001280973317475672')+edge_bound),
                derivation='For an odd pair, |C_i term|<=F(x)/(pi*h*a^3), F=coth(x)+x*csch(x)^2+x^2*coth(x)*csch(x)^2. F decreases for x>=1 and F(10)<1.001; omitted x>3108. Odd max-norm shells have 4*k points, giving sum_tail r^-3<=2/(M+1). Sum all three pair bounds; omitted rho is negative, so p_tail<=|C_i tail|',
                scope='Partial pressure is from the separately reported M2501 audit; its large finite sum is not rerun here. This analytic tail upper bound does not certify arithmetic rounding'),
            sources=['https://arxiv.org/pdf/0804.3916 (Lim-Teo: physical EM conducting-cavity mode framework)',
                     'https://arxiv.org/pdf/hep-th/0202089 (Nesterenko et al.: local wedge/plate stress framework)'],
            scope='Our derived ideal-rectangle center stress and independent numerical verification. No originality claim, realistic-conductor result, spatially extended source, self-consistent curved geometry or traversal conclusion')


def temporal_density_noise_checks(D, p):
    """Time-smeared center density variance of the threshold coherent mixture."""
    assert math.isfinite(D) and math.isfinite(p) and D > 0 and p < 0 and D+p > 0
    A, volume = (D+p)/2, math.sqrt(5/8)
    Omega, Omega_z = math.pi*math.sqrt(97/5), math.pi*math.sqrt(2)

    def spectrum(cutoff):
        n = np.indices((cutoff+1,)*3).reshape(3, -1)
        bits = n % 2
        odd, zeros = bits.sum(axis=0), (n == 0).sum(axis=0)
        lam = (np.array([5, 5, 8])[:, None]*n*n).sum(axis=0)
        active = ((odd == 1) | (odd == 2)) & (zeros <= 1)
        special_axis = np.argmax(bits == np.where(odd == 1, 1, 0), axis=0)
        special_square = (np.array([5, 5, 8])[special_axis]
                          *n[special_axis, np.arange(n.shape[1])]**2)
        denominator = np.where(active, lam, 1)
        omega = math.pi*np.sqrt(denominator/5)
        weight = np.where(active, 4*omega/volume*(1-special_square/denominator)
                          *np.where(zeros == 1, .5, 1.), 0.)
        channel = special_axis+np.where(odd == 1, 3, 0)
        W = np.array([np.bincount(lam[active & (channel == c)],
            weights=weight[active & (channel == c)], minlength=int(lam.max())+1)
            for c in range(6)])
        keep = np.any(W > 0, axis=0)
        return np.flatnonzero(keep), math.pi*np.sqrt(np.flatnonzero(keep)/5), W[:, keep]

    def omitted_bound(sigma, cutoff):
        # Max-index shell has <=7*m^2 triples; total channel weight <=K*m^3.
        with mp.workdps(50):
            a = (mp.pi*mp.mpf(str(sigma)))**2
            beta = mp.sqrt(mp.mpf(97)/5)
            N = mp.mpf(cutoff)
            K = 28*mp.pi*mp.sqrt(3)/mp.mpf(str(volume*math.sqrt(5/8)))
            assert N > beta and 2*a*(N-beta)*N > 3
            tail = K*mp.exp(-a*N*N)*(a*N*N+1)/(2*a*a)
            # Integral plus maximum bounds the sum of a nonnegative unimodal function.
            all_weight = K*(1/(2*a*a)+(3/(2*a*mp.e))**mp.mpf('1.5'))
            u, e = N-beta, mp.exp(-a*(N-beta)**2)
            I0 = mp.sqrt(mp.pi)*mp.erfc(mp.sqrt(a)*u)/(2*mp.sqrt(a))
            I1 = e/(2*a)
            I2 = u*e/(2*a)+I0/(2*a)
            I3 = (a*u*u+1)*e/(2*a*a)
            shifted_tail = K*(I3+3*beta*I2+3*beta*beta*I1+beta**3*I0)
            return tail*all_weight+2*mp.mpf(str(A-p))*shifted_tail

    def evaluate(sigma, cutoff):
        lam, omega, W = spectrum(cutoff)
        vacuum_parts = []
        for start in range(0, len(omega), 96):
            block = slice(start, start+96)
            kernel = np.exp(-(sigma*(omega[block, None]+omega[None, :]))**2)
            vacuum_parts.append(float(np.sum(kernel*(W[:, block].T @ W))))
        vacuum = math.fsum(vacuum_parts)/2
        pair_kernel = (np.exp(-(sigma*(omega-Omega))**2)
                       +np.exp(-(sigma*(omega+Omega))**2))
        z_kernel = (np.exp(-(sigma*(omega-Omega_z))**2)
                    +np.exp(-(sigma*(omega+Omega_z))**2))
        linear = float(A/2*np.dot(W[0]+W[4], pair_kernel)-p*np.dot(W[2], z_kernel))
        phase = A*A/2*math.exp(-4*(sigma*Omega)**2)+p*p/2*math.exp(-4*(sigma*Omega_z)**2)
        total = vacuum+linear+phase
        assert all(math.isfinite(v) and v >= 0 for v in (vacuum, linear, phase, total))
        bound = omitted_bound(sigma, cutoff)
        return dict(sigma=sigma, cutoff=cutoff, spectral_bins=len(lam),
            vacuum_variance=vacuum, displaced_linear_variance=linear,
            classical_phase_variance=phase, total_variance=total,
            mean_density=-A, relative_standard_deviation=math.sqrt(total)/A,
            analytic_omitted_mode_bound=mp.nstr(bound, 12))

    lam, omega, W = spectrum(5)
    expected_weights = [308*Omega/(97*volume), 208*Omega/(97*volume), 2*Omega_z/volume]
    measured_weights = [float(W[0, lam == 97][0]), float(W[4, lam == 97][0]),
                        float(W[2, lam == 10][0])]
    assert np.allclose(measured_weights, expected_weights, rtol=2e-15, atol=0)
    floor = A/2*(expected_weights[0]+expected_weights[1])-p*expected_weights[2]
    rows = []
    for sigma in (.25, .5, 1., 2., 4., 8., 32.):
        cutoff = math.ceil(Omega/math.pi+10/(math.pi*sigma))
        row = evaluate(sigma, cutoff)
        assert row['total_variance'] >= floor-2e-13
        rows.append(row)
    refined = evaluate(.25, rows[0]['cutoff']+2)
    disagreement = max(abs(rows[0][key]-refined[key]) for key in
        ('vacuum_variance', 'displaced_linear_variance', 'total_variance'))
    floating_tolerance = 5e-12*max(1., rows[0]['total_variance'])
    assert disagreement <= floating_tolerance
    assert abs(rows[-1]['total_variance']-floor) < 2e-13
    # Six channels, Wick factor1/2, two omega^3/(6*pi^2) spectra, beta and radial integrals.
    assert Fraction(3)*Fraction(1, 36)*Fraction(36, 5040)*3 == Fraction(1, 560)
    return dict(
        model='Final threshold state of the ideal rectangle, A=(D+p)/2; one common phase for Ex/By and an independent Ez phase',
        observable='rho_f = integral dt f_sigma(t)*T00(t,center), f_sigma=exp(-t^2/(2*sigma^2))/(sqrt(2*pi)*sigma); sigma in L/c, rho in hbar*c/L^4, variance in (hbar*c/L^4)^2',
        renormalization='Subtract the same cavity-vacuum c-number mean; this does not change variance. No replacement by photon-number variance or normal-ordering away of vacuum noise',
        mode_convention='F=(Ex,Ey,Ez,Bx,By,Bz). Each physical center mode has one nonzero channel; W_cc(dt)=sum_omega w_c(omega)*exp(-i*omega*dt), all off-diagonal vacuum correlations vanish',
        spectrum='lambda=5*nx^2+5*ny^2+8*nz^2, omega=pi*sqrt(lambda/5); w=4*omega/V*(1-kj^2/k^2)*2^(-zero_count). Unique even parity gives an E channel; unique odd gives a B channel. Polarizations are summed by the transverse projector',
        vacuum_formula='Var0=1/2*sum_c,omega,nu w_c(omega)*w_c(nu)*exp[-sigma^2*(omega+nu)^2]. Ordered pairs include equal modes. Wick contraction is W_ab^2, not absolute(W_ab)^2',
        general_mode_formula='1/2*sum_mn absolute(v_m^T*v_n)^2*absolute(f_tilde(omega_m+omega_n))^2; ordinary transpose inside the modulus, not a dagger. v=(i*sqrt(omega/2)*u,curl(u)/sqrt(2*omega))',
        displaced_linear_formula='A/2*sum(w_Ex+w_By)*(exp[-sigma^2*(omega-Omega)^2]+exp[-sigma^2*(omega+Omega)^2]) - p*sum w_Ez*(exp[-sigma^2*(omega-Omega_z)^2]+exp[-sigma^2*(omega+Omega_z)^2])',
        phase_formula='A^2/2*exp(-4*sigma^2*Omega^2)+p^2/2*exp(-4*sigma^2*Omega_z^2); law of total variance, with independent uniform phases',
        total_variance='vacuum + displaced_linear + classical_phase; all contributions nonnegative',
        shell_weights=dict(lambda_97_Ex=expected_weights[0], lambda_97_By=expected_weights[1], lambda_10_Ez=expected_weights[2],
            exact='w_Ex(97)=308*Omega/(97V), w_By(97)=208*Omega/(97V), w_Ez(10)=2*Omega_z/V, independently enumerated; unused degenerate polarizations remain in the vacuum noise'),
        gaussian_long_time_limit=dict(variance=floor, relative_standard_deviation=math.sqrt(floor)/A,
            formula='A/2*(w_Ex(Omega)+w_By(Omega))-p*w_Ez(Omega_z)',
            reason='The resonant difference-frequency terms remain1 for every sigma; all other terms are nonnegative. This is both a Gaussian-width-independent lower bound and the sigma->infinity limit'),
        general_normalized_sampling_bound=dict(variance_lower_bound=floor,
            relative_standard_deviation_lower_bound=math.sqrt(floor)/A,
            assumptions='Any real normalized temporal sampling g with integral(g)=1 and a well-defined finite second moment, for this exact stationary state and center. Normalization alone does not ensure finite variance',
            argument='Replace each Gaussian kernel by absolute(g_tilde(frequency))^2. Independent phases eliminate cross terms. All coefficients are nonnegative; the resonant difference-frequency terms contain absolute(g_tilde(0))^2=1 and supply the stated bound, including signed, asymmetric or time-shifted samplings',
            scope='This state and center-density observable only; not every quantum source or a gravitational stability theorem'),
        samples=rows, refinement=dict(sigma=.25, second_cutoff=refined['cutoff'],
            maximum_component_disagreement=disagreement, floating_tolerance=floating_tolerance),
        mode_tail=dict(shell_bound='K*m^3, K=28*pi*sqrt(3)/(V*min_side), omega>=pi*m',
            vacuum_bound='(sum_tail w*exp(-sigma^2*omega^2))*(sum_all w*exp(-sigma^2*omega^2)); integral/maximum bounds on m^3*exp(-a*m^2)',
            linear_bound='2*(A-p)*K*integral_N^infinity m^3*exp[-a*(m-sqrt(97/5))^2] dm; integrand decreasing beyond N',
            scope='Bounds cover the omitted physical modes, not floating-point roundoff or uncertainty in the ideal-wall model. Finite sums use float64; positive-tail expressions use50-digit mpmath. At wide sampling widths, positive exponentially small vacuum/phase terms can underflow to numerical zero; they are not exactly zero at finite width'),
        minkowski_benchmark='Var_vac=1/(560*pi^4*sigma^8), the analytic free-space short-sampling limit, not an exact finite-width cavity formula',
        sources=['https://arxiv.org/pdf/1204.3570 (Fewster-Ford-Roman, Eq16 and Eqs37,41-42: time-averaged quadratic stress moments; our Gaussian cavity substitution)',
                 'https://arxiv.org/pdf/0804.3916 (Lim-Teo: physical conducting-cavity modes)',
                 'https://arxiv.org/pdf/gr-qc/0209075 (Anderson et al.: gravitational linear-response validity requires more than a local variance ratio)'],
        scope='Actual time-smeared center-density second moment for this specified ideal quantum state. No material-wall noise, spatially extended source, detector dynamics or gravitational backreaction has been calculated; a large ratio alone is not a general semiclassical-gravity no-go theorem')


def spatial_source_checks(D,p):
    """Full spatial vacuum and fixed threshold state; constant-center diagnostic only."""
    import sympy as sp

    assert math.isfinite(D) and math.isfinite(p) and D>0 and p<0 and D+p>0
    with mp.workdps(28):
        SIDES=(mp.mpf(1),mp.mpf(1),mp.sqrt(mp.mpf(5)/8))
        GROUPS=((1,1,1),(-1,-1,1),(-1,1,-1),(1,-1,-1))
        THETA_N=6
        
        def image_stress_at(offset,tau=0):
            position=tuple(L/2+mp.mpf(str(d)) for L,d in zip(SIDES,offset))
            tau=mp.mpf(str(tau))
            @lru_cache(maxsize=None)
            def theta(v,shift):
                symmetric=shift==0 or shift==mp.mpf('.5')
                if v>=1:
                    if shift==0:
                        terms=[(mp.mpf(n),mp.exp(-mp.pi*v*n*n)) for n in range(1,THETA_N+1)]
                        delta=2*mp.fsum(t for n,t in terms)
                        return 1+delta,mp.mpf(0),2*mp.fsum(n*n*t for n,t in terms),delta
                    terms=[(mp.mpf(n)-shift,mp.exp(-mp.pi*v*(n-shift)**2)) for n in range(-THETA_N,THETA_N+1)]
                    T=mp.fsum(t for n,t in terms)
                    M1=mp.mpf(0) if symmetric else mp.fsum(n*t for n,t in terms)
                    M2=mp.fsum(n*n*t for n,t in terms)
                    return T,M1,M2,T-1
                terms=[(mp.mpf(n),mp.exp(-mp.pi*n*n/v)) for n in range(1,THETA_N+1)]
                H=1+2*mp.fsum(t*mp.cos(2*mp.pi*n*shift) for n,t in terms)
                M1=mp.mpf(0) if symmetric else -2/v**mp.mpf('1.5')*mp.fsum(n*t*mp.sin(2*mp.pi*n*shift) for n,t in terms)
                M2=H/(2*mp.pi*v**mp.mpf('1.5'))-2/v**mp.mpf('2.5')*mp.fsum(n*n*t*mp.cos(2*mp.pi*n*shift) for n,t in terms)
                T=H/mp.sqrt(v)
                return T,M1,M2,T-1
            @lru_cache(maxsize=None)
            def integrands(u):
                if not u: return (mp.mpf(0),)*6
                entries=[mp.mpf(0) for i in range(6)]
                z=mp.pi*u*tau*tau
                for eps in GROUPS:
                    rows=[theta(4*L*L*u,mp.mpf(0) if eps[i]==1 else position[i]/L) for i,L in enumerate(SIDES)]
                    if eps==(1,1,1):
                        d=[row[3] for row in rows]
                        K=sum(d)+d[0]*d[1]+d[0]*d[2]+d[1]*d[2]+mp.fprod(d)
                    else:
                        K=mp.fprod(row[0] for row in rows)
                    for i in range(3):
                        M2=4*SIDES[i]**2*rows[i][2]*mp.fprod(rows[j][0] for j in range(3) if j!=i)
                        entries[i]+=2*eps[i]*(u*(z-1)*K+mp.pi*u*u*M2)
                    for h,(i,j) in enumerate(((0,1),(0,2),(1,2)),3):
                        if eps[i]+eps[j]:
                            Mij=4*SIDES[i]*SIDES[j]*rows[i][1]*rows[j][1]*rows[3-i-j][0]
                            entries[h]+=mp.pi*(eps[i]+eps[j])*u*u*Mij
                return tuple(mp.exp(-z)*v for v in entries)
            values=[]; errors=[]
            for h in range(6):
                if h>=3:
                    i,j=((0,1),(0,2),(1,2))[h-3]
                    if offset[i]==0 or offset[j]==0:
                        values.append(mp.mpf(0)); errors.append(mp.mpf(0)); continue
                value,error=mp.quad(lambda u:integrands(u)[h],[0,mp.mpf('.1'),mp.mpf('.25'),mp.mpf('.5'),1,2,4,8,16,mp.inf],error=True)
                values.append(value); errors.append(error)
            C=mp.matrix([[values[0],values[3],values[4]],[values[3],values[1],values[5]],[values[4],values[5],values[2]]])
            rho=(C[0,0]+C[1,1]+C[2,2])/2
            T=rho*mp.eye(3)-C
            return dict(offset=offset,tau=float(tau),rho=mp.nstr(rho,25),stress=[[mp.nstr(T[i,j],25) for j in range(3)] for i in range(3)],quad_error_estimate=mp.nstr(max(errors),8))
        
        def physical_stress_at(offset,tau):
            L=np.array([1.,1.,math.sqrt(5/8)])
            V=float(np.prod(L)); x=L/2+np.array(offset,dtype=float)
            N=math.ceil(60/(math.pi*tau))
            ns=np.arange(N+1,dtype=float)
            ny,nz=np.meshgrid(ns,ns,indexing='ij')
            ky,kz=ny/L[1],nz/L[2]
            sy,cy=np.sin(math.pi*ky*x[1]),np.cos(math.pi*ky*x[1])
            sz,cz=np.sin(math.pi*kz*x[2]),np.cos(math.pi*kz*x[2])
            zero_yz=(ny==0).astype(int)+(nz==0).astype(int)
            parts=[[] for _ in range(6)]
            pairs=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
            for nx in range(N+1):
                kx=nx/L[0]
                sx,cx=math.sin(math.pi*kx*x[0]),math.cos(math.pi*kx*x[0])
                q=(kx,ky,kz); R2=kx*kx+ky*ky+kz*kz
                zeros=zero_yz+(nx==0); active=zeros<=1
                den=np.where(active,R2,1.)
                omega=math.pi*np.sqrt(den)
                pref=np.where(active,4*omega/V*np.where(zeros==1,.5,1.),0.)*np.exp(-tau*omega)
                v=(cx*sy*sz,sx*cy*sz,sx*sy*cz)
                w=(sx*cy*cz,cx*sy*cz,cx*cy*sz)
                for h,(i,j) in enumerate(pairs):
                    P=float(i==j)-q[i]*q[j]/den
                    parts[h].append(float(np.sum(pref*P*(v[i]*v[j]+w[i]*w[j]))))
            vals=np.array([math.fsum(p) for p in parts])
            vals[:3]-=2/(math.pi**2*tau**4)
            C=np.array([[vals[0],vals[3],vals[4]],[vals[3],vals[1],vals[5]],[vals[4],vals[5],vals[2]]])
            rho=float(np.trace(C)/2)
            return dict(cutoff=N,rho=rho,stress=(rho*np.eye(3)-C).tolist())
        
        def added_stress(xyz,D,p):
            xyz=np.asarray(xyz,dtype=float)
            lengths=np.array([1.,1.,np.sqrt(5/8)])
            omega=np.pi*np.sqrt(97/5); K=np.pi*np.sqrt(52/5)
            def mode(n,alpha):
                k=np.pi*np.asarray(n)/lengths
                sn,cs=np.sin(k*xyz),np.cos(k*xyz)
                f=np.array([alpha[i]*cs[i]*sn[(i+1)%3]*sn[(i+2)%3] for i in range(3)])
                b=np.array([(alpha[(i+2)%3]*k[(i+1)%3]-alpha[(i+1)%3]*k[(i+2)%3])*sn[i]*cs[(i+1)%3]*cs[(i+2)%3] for i in range(3)])
                return f,b
            f,cf=mode((2,1,3),np.array([1.,-2.,0.])/np.sqrt(5))
            g,cg=mode((2,3,2),np.array([-2*np.sqrt(2),0.,np.sqrt(5)])/np.sqrt(13))
            ec,es,bc,bs=np.sqrt(5)*f,-omega*g/K,-cg/K,-np.sqrt(5)*cf/omega
            sx,sy=np.sin(np.pi*xyz[:2]); cx,cy=np.cos(np.pi*xyz[:2])
            ez=np.array([0.,0.,sx*sy]); bz=np.array([-sx*cy,cx*sy,0.])/np.sqrt(2)
            A=(D+p)/2
            C=A*sum(np.outer(v,v) for v in (ec,es,bc,bs))-2*p*(np.outer(ez,ez)+np.outer(bz,bz))
            rho=np.trace(C)/2; flux=A*(np.cross(ec,bc)+np.cross(es,bs))
            T=np.zeros((4,4)); T[0,0]=rho; T[0,1:]=T[1:,0]=flux
            T[1:,1:]=rho*np.eye(3)-C
            return T
        
        
        center=np.array([float(L/2) for L in SIDES])
        A=(D+p)/2
        ell=np.array([1.,0.,0.,1.])
        target=-A*np.outer(ell,ell)
        
        @lru_cache(maxsize=None)
        def vacuum_tensor(offset,tau=0.):
            if offset[0]==0 and offset[1]!=0 and offset[2]==0:
                matrix,estimate=vacuum_tensor((offset[1],0.,0.),tau)
                perm=[0,2,1,3]
                return matrix[np.ix_(perm,perm)],estimate
            result=image_stress_at(offset,tau)
            matrix=np.zeros((4,4))
            matrix[0,0]=float(result['rho'])
            matrix[1:,1:]=np.array(result['stress'],dtype=float)
            return matrix,result['quad_error_estimate']
        
        def sample(offset):
            vacuum,estimate=vacuum_tensor(tuple(offset))
            added=added_stress(center+np.array(offset),D,p)
            total=vacuum+added
            delta=total-target
            index=np.unravel_index(np.argmax(np.abs(delta)),delta.shape)
            assert np.max(np.abs(total-total.T))<1e-12
            assert abs(-total[0,0]+np.trace(total[1:,1:]))<1e-12
            return dict(offset_L=list(offset),vacuum_tensor=vacuum.tolist(),
                added_tensor=added.tolist(),total_tensor=total.tolist(),
                maximum_componentwise_error_over_A=float(np.max(np.abs(delta))/A),
                worst_component=[int(i) for i in index],
                vacuum_quadrature_error_estimate=estimate)
        
        at_center=sample((0.,0.,0.))
        expected_vacuum=np.diag([-D,p,p,-D-2*p])
        expected_added=np.diag([A-p,-p,-p,A+p])
        expected_added[0,3]=expected_added[3,0]=-A
        assert np.max(np.abs(np.array(at_center['vacuum_tensor'])-expected_vacuum))<1e-12
        assert np.max(np.abs(np.array(at_center['added_tensor'])-expected_added))<1e-12
        assert at_center['maximum_componentwise_error_over_A']<1e-12
        offsets=(.001,.005,.01,.02,.05)
        axes={}
        for axis,name in enumerate(('x','y','z')):
            rows=[]
            for delta in offsets:
                offset=[0.,0.,0.];offset[axis]=delta
                rows.append(sample(tuple(offset)))
            axes[name]=rows
        generic=[sample(offset) for offset in ((.01,.01,.01),(.01,.015,.02))]
        envelope=[]
        for i,delta in enumerate(offsets):
            worst=max(axes,key=lambda axis:axes[axis][i]['maximum_componentwise_error_over_A'])
            envelope.append(dict(offset_L=delta,worst_axis=worst,
                maximum_componentwise_error_over_A=axes[worst][i]['maximum_componentwise_error_over_A']))
        brackets={}
        for axis,rows in (*axes.items(),('sampled_axis_envelope',envelope)):
            crossings=[]
            for i in range(len(rows)-1):
                if rows[i]['maximum_componentwise_error_over_A']<=.01<rows[i+1]['maximum_componentwise_error_over_A']:
                    crossings.append([offsets[i],offsets[i+1]])
            brackets[axis]=crossings
        
        probe=(.01,.015,.02);tau=.15
        image_matrix,estimate=vacuum_tensor(probe,tau)
        physical=physical_stress_at(probe,tau)
        mode_matrix=np.zeros((4,4));mode_matrix[0,0]=physical['rho']
        mode_matrix[1:,1:]=physical['stress']
        mode_error=float(np.max(np.abs(mode_matrix-image_matrix)))
        minkowski_rho=3/(math.pi**2*tau**4)
        acceptance=max(1e-11,8*np.finfo(float).eps*(2*minkowski_rho+abs(physical['rho'])))
        assert math.isfinite(mode_error) and mode_error<=acceptance
        M=physical['cutoff']+1;r=mp.exp(-mp.pi*tau)
        volume=mp.fprod(SIDES)
        mode_tail=14*mp.pi*mp.sqrt(3)/(volume*min(SIDES))*r**M*(
            M**3/(1-r)+3*M*M*r/(1-r)**2+3*M*r*(1+r)/(1-r)**3+r*(1+4*r+r*r)/(1-r)**4)
        
        # Each double reflection is conserved separately; physical reflected derivatives add a common factor2.
        su,sv,sw,st=sp.symbols('du dv dw imaginary_tau',real=True)
        R=su*su+sv*sv+sw*sw+st*st
        image_rho=(-su*su-sv*sv+3*sw*sw-st*st)/R**3
        image_spatial=sp.Matrix([
            [(su*su-3*sv*sv+sw*sw+st*st)/R**3,4*su*sv/R**3,0],
            [4*su*sv/R**3,(-3*su*su+sv*sv+sw*sw+st*st)/R**3,0],
            [0,0,(su*su+sv*sv+sw*sw-3*st*st)/R**3]])
        assert sp.simplify(sp.trace(image_spatial)-image_rho)==0
        assert sp.simplify(sp.diff(image_spatial[0,0],su)+sp.diff(image_spatial[0,1],sv))==0
        assert sp.simplify(sp.diff(image_spatial[1,0],su)+sp.diff(image_spatial[1,1],sv))==0
        step=1e-5;point=center+np.array(probe)
        divergence=sum((added_stress(point+np.eye(3)[i]*step,D,p)[i+1]
            -added_stress(point-np.eye(3)[i]*step,D,p)[i+1])/(2*step) for i in range(3))
        assert np.max(np.abs(divergence))<1e-6
        wall_flux=[]
        for axis in range(3):
            for side in (0.,float(SIDES[axis])):
                point=np.array([.37,.43,.29]);point[axis]=side
                wall_flux.append(abs(added_stress(point,D,p)[0,axis+1]))
        assert max(wall_flux)<1e-12
        
        return dict(
            model='Actual ideal-rectangle vacuum plus the fixed stationary threshold coherent-mixture state; amplitudes fixed by center D,p and not retuned at other positions',
            convention='Contravariant Cartesian orthonormal tensor ordered t,x,y,z, signature -+++; Lx=Ly=L, Lz=sqrt(5/8)*L. xyz and offsets are in units L even along z; tensor coefficients multiply hbar*c/L^4',
            center_position_L=center.tolist(),D=float(D),p=float(p),A=float(A),
            reference_tensor=target.tolist(),
            diagnostic='max_ab |T_total^{ab}(x)-T_reference^{ab}|/A, with T_reference=-A*ell^a*ell^b constant in these flat Cartesian components. This is a declared center-reference diagnostic, not the actual spatially varying curved Ori target, an invariant tensor norm, or a gravitational solution',
            center=at_center,positive_axis_samples=axes,generic_samples=generic,
            sampled_axis_envelope=envelope,
            one_percent_sample_brackets=brackets,
            bracket_scope='Each listed adjacent sample pair brackets at least one crossing of the 1% componentwise diagnostic by continuity. No interpolated location, uniqueness, monotonicity, or three-dimensional ball is inferred from axis samples',
            image_geometry='eps_i=+/-1, S=product(eps_i), d_i=(1-eps_i)*x_i-2*Li*n_i, U=sum(d_i^2)+tau^2. Exclude only the zero translation image. Single/triple reflections cancel; keep translations and the three double-reflection families',
            full_image_bilinear='Cii=sum(1+S)*eps_i/pi^2*(-1/U^2+2*(d_i^2+tau^2)/U^3); Cij=sum(1+S)*(eps_i+eps_j)*d_i*d_j/(pi^2*U^3) for i!=j. rho=trace(C)/2, Tij=rho*deltaij-Cij; physical static-vacuum T0i=0',
            shifted_theta='theta(v,beta)=sum_n exp[-pi*v*(n-beta)^2], v=4*Li^2*u and beta=x_i/Li for a reflected coordinate, otherwise0. Direct shifted Gaussians for v>=1; Poisson form v^(-1/2)*[1+2*sum_n>=1 exp(-pi*n^2/v)*cos(2*pi*n*beta)] for v<1. Differentiate to obtain first/second moments before the Laplace integrals',
            physical_modes='No center parity reduction: vi=cos(ki*xi)*prod_{j!=i}sin(kj*xj), wi=sin(ki*xi)*prod_{j!=i}cos(kj*xj), Pij=deltaij-ki*kj/omega^2. Cij_mode=(4*omega/V)*2^(-zero_count)*Pij*(vi*vj+wi*wj)*exp(-tau*omega), for physical all-positive and one-zero modes',
            noncenter_finite_tau_check=dict(offset_L=list(probe),tau=tau,cutoff=physical['cutoff'],
                image_tensor=image_matrix.tolist(),physical_mode_tensor=mode_matrix.tolist(),
                maximum_energy_or_spatial_stress_disagreement=mode_error,floating_acceptance=acceptance,
                analytic_mode_tail_bound=float(mode_tail),image_quadrature_error_estimate=estimate,
                scope='Imaginary-time split for checking energy and every spatial stress component; not real-time evolution. Subtract Minkowski Cii=2/(pi^2*tau^4). The positive mode-energy tail bounds every spatial stress entry; it excludes arithmetic rounding'),
            checks=dict(center_match=True,symmetric_and_traceless_samples=True,
                image_spatial_conservation='Exact double-reflection algebra checked; translations are constant. Each double-reflection family is symmetric, traceless and spatially conserved',
                axis_vacuum_symmetry='Off-diagonal vacuum stresses vanish on a single coordinate axis by the remaining midplane reflections. The y-axis vacuum is obtained from x by exchanging x/y tensor indices; the actual classical state is evaluated separately on every axis',
                classical_divergence_finite_difference_step_L=step,
                classical_divergence_residual=divergence.tolist(),
                maximum_classical_divergence_residual=float(np.max(np.abs(divergence))),
                maximum_added_normal_wall_flux=float(max(wall_flux)),
                classical_exact_framework='The same full classical mode construction was independently checked for exact conservation, PEC conditions and phase averaging; this callable is also checked numerically here'),
            precision=dict(decimal_working_digits=28,theta_cutoff=THETA_N,
                scope='Displayed tensors are float64 after high-precision image quadrature. Quadrature estimates are not interval certification and do not include every arithmetic/truncation error',
                prior_read_only_audit='The full tensor at center+(.01,.015,.02) also matched physical modes at imaginary tau=.08 within3.86e-13. At offsets(.02,0,0) and(.01,.01,.01), 32-digit/cutoff7 and44-digit/cutoff9 runs with different quadrature partitions agreed in25 displayed digits. Those refinements are not repeated by this default check'),
            sources=['https://arxiv.org/pdf/0804.3916 (Lim-Teo: physical EM conducting-cavity mode framework)',
                     'https://doi.org/10.1103/PhysRev.131.2766 (Glauber: coherent-state expectation framework; the specific spatial construction is derived here)'],
            scope='A finite set of spatial diagnostics for this fixed state and ideal flat cavity. No spatial retuning, three-dimensional error guarantee, realistic walls, self-consistent curved quantum state, stability or traversal result')


def run_checks():
    with mp.workdps(50):
        normalization = exact_normalization_checks()
        refinement = []
        for N, M in ((2, 5), (3, 9), (4, 13), (5, 17)):
            row = center_images(ewald_cutoff=N, odd_cutoff=M)
            refinement.append(dict(ewald_cutoff=N, odd_cutoff=M,
                q=mp.nstr(row['value'], 36),
                analytic_truncation_bound=float(row['analytic_truncation_bound'])))
        for first, second in zip(refinement, refinement[1:]):
            assert abs(mp.mpf(first['q'])-mp.mpf(second['q'])) <= (
                first['analytic_truncation_bound']+second['analytic_truncation_bound'])
        final = center_images(ewald_cutoff=5, odd_cutoff=17)
        q = float(final['value'])
        assert final['translation'] < 0 and final['edge'] < 0
        rows = []
        for tau in (.4, .25, .15, .1, .06):
            observed, reference = transverse_modes(tau), center_images(tau)
            error = abs(observed['value']-float(reference['value']))
            # Practical floating-point acceptance, not a certified rounding bound.
            acceptance = max(1e-11, 8*np.finfo(float).eps*(observed['raw']+observed['minkowski']))
            assert math.isfinite(observed['value']) and error <= acceptance
            rows.append(dict(tau=tau, **observed, image_reference=float(reference['value']),
                image_analytic_truncation_bound=float(reference['analytic_truncation_bound']),
                absolute_disagreement=error, floating_acceptance=acceptance))
        # Stated constants; the following is a formal perfect-conductor extrapolation.
        light_speed, planck_h, newton_G = 299792458., 6.62607015e-34, 6.67430e-11
        hbar = planck_h/(2*math.pi)
        mu, geometry_length = -11251/11250, 1.
        target_A = abs(mu)*light_speed**4/(16*math.pi*newton_G*geometry_length**2)
        formal_side = (abs(q)*hbar*light_speed/(4*target_A))**.25
        rectangle = rectangle_vacuum_checks()
        return dict(
            model='Single electromagnetic vacuum in a flat static cube with six perfectly conducting faces, zero temperature; local center only',
            normalization='rho_center=q*hbar*c/L^4; L is the cube side; imaginary time regulator tau is measured in L/c',
            subtraction='Point-split local field-strength correlations minus the Minkowski bulk singularity 3/(pi^2*tau^4); no global-energy/V substitution or subtraction of finite edge images',
            q=q, q_decimal=mp.nstr(final['value'], 36), pressure_coefficient=q/3,
            image_kernel='For m in Z^3 except zero, s_i=(-1)^m_i, S=product(s_i), r=|m|: f=(1+S)/(2*pi^2)*[-sum(s_i)/r^4+2*sum(s_i*m_i^2)/r^6]',
            image_basis='A0 Dirichlet on every face; Ai Neumann on its normal faces and Dirichlet on the others. Covariant image propagator D00=-S*G, Dii=S*s_i*G, G=[4*pi^2*(r^2-(dt-i0)^2)]^-1; differentiate to E and B correlations',
            resummation='q=-Z3(2)/(16*pi^2)-(3*pi/8)*sum_{p,q odd}[coth(pi*r/2)*csch(pi*r/2)^2/r], r=sqrt(p^2+q^2); both contributions strictly negative',
            zeta_definition='Z3(s)=sum_{n!=0}|n|^(-2*s)',
            ewald_identity='Z3(s)=pi^s/Gamma(s)*[1/(s-3/2)-1/s+sum_{n!=0}(K^-s*Gamma(s,K)+K^(s-3/2)*Gamma(3/2-s,K))], K=pi*|n|^2',
            bounds=dict(ewald='theta^3-theta_N^3<=8*exp[-pi*(N+1)^2*t] for t>=1; integrate the positive theta-split kernel',
                edge='For odd-grid cutoff M, omitted |q_edge|<=15*pi*exp[-pi*(M+1)]/(1-exp[-pi])',
                regulator_series='With u=tau^2/4<1, translation=-sum_k (-1)^k*(k+1)*(2*k+1)*u^k*Z3(k+2)/(16*pi^2); polynomial-geometric tail uses Z3(k+2)<=13*pi^2/3',
                modes='For cutoff N, omitted positive mode energy<=14*pi*sqrt(3)*sum_{m=N+1}^infinity m^3*exp[-pi*tau*m]',
                final_analytic_truncation=float(final['analytic_truncation_bound']),
                rounding='Analytic tail bounds exclude arithmetic rounding. Ewald special functions use60 decimal digits, other image arithmetic50, mode sums float64; no interval certification'),
            refinement=refinement, exact_checks=normalization, coherent_center=coherent_center_checks(),
            rectangular_resonance=rectangular_resonance_checks(),
            rectangle_vacuum=rectangle,
            temporal_density_noise=temporal_density_noise_checks(rectangle['D_coefficient'], rectangle['transverse_pressure_coefficient']),
            spatial_source=spatial_source_checks(rectangle['D_coefficient'], rectangle['transverse_pressure_coefficient']),
            mode_basis='Lim-Teo physical transverse modes, omega=pi*|n|; two polarizations for all indices positive, one for exactly one zero, none for two zero. Center parity and polarization projector determine local weights',
            mode_weights='For odd count1 or2, W=2*omega*(1-n_j^2/|n|^2)*2^(-zero_count), where j is the unique odd or unique even index; W=0 for odd count0 or3',
            finite_regulator_checks=rows, maximum_mode_image_disagreement=max(row['absolute_disagreement'] for row in rows),
            local_compensation=dict(convention='Contravariant orthonormal tensor, signature -+++; target=-A*ell^a*ell^b, ell=(1,0,0,1), A>0; D=-rho_center>0',
                criterion='A<=D/4, equivalently L<=(abs(q)*hbar*c/(4*A))^(1/4)',
                necessity='Residual radiation moment minor: (D-A)*(D/3-A)-A^2=D*(D-4*A)/3>=0',
                exact_threshold_witness='At D=4A, four beams each of weight3A/4 with directions (+/-2sqrt(2)/3,0,-1/3) and (0,+/-2sqrt(2)/3,-1/3). Null norms and diagonal/flux moments checked with rational coefficients',
                below_threshold='At fixed D, add (D/4-A)*ell^a*ell^b to the threshold witness',
                limitation='A cavity-compatible mixed coherent state is supplied only at the center at one event. A persistent finite spatial region, supporting walls and their stresses, and a self-consistent gravitational solution are not supplied'),
            formal_SI_example=dict(c_m_per_s=light_speed, h_J_s=planck_h, G_m3_per_kg_s2=newton_G,
                constants='c and h are exact SI defining values; G is the stated central value6.67430e-11, with its uncertainty not propagated',
                rho_at_one_micrometer_J_per_m3=q*hbar*light_speed/1e-24,
                mu=mu, geometry_length_m=geometry_length,
                central_observer='At Ori N, g_Tz=-1 and g_zz=0; U=(partial_T+partial_z)/sqrt(2), so U^z squared=1/2',
                target_formula='A=abs(mu)*c^4/(16*pi*G*ell^2)', target_A_J_per_m3=target_A,
                formal_max_cube_side_m=formal_side, scaling='formal Lmax scales as sqrt(ell)',
                limitation='ell is the geometry unit, not a traveler aperture. Perfect conductivity extrapolated to Lmax is not a reliable material model or fabrication recipe; no parallel-plate conductivity correction has been applied to the cube'),
            global_energy_limit='The ideal sharp-edge local stress has nonintegrable edge behavior. A globally regularized total energy divided by volume need not equal this interior point value',
            sources=['https://arxiv.org/pdf/0804.3916 (Lim-Teo, section5.C.1: common conducting-cavity physical modes)',
                     'https://arxiv.org/pdf/hep-th/0202089 (Nesterenko et al., Eqs3.38-3.41: local wedge/plate normalization and nonintegrable edge behavior)',
                     'https://www.quantumfields.com/old/physreva_cavity.pdf (Maclay, spatially averaged/global cube energy; not the local input)',
                     'https://arxiv.org/pdf/quant-ph/9907105 (Lambrecht-Reynaud: finite reflectivity; no cube correction inferred)'],
            scope='Our local flat-ideal-cube derivation and numerical check; no originality claim, no realistic-conductor validation, no curved periodic Ori quantum state, and no stability or traversal conclusion')


if __name__ == '__main__':
    result = run_checks()
    output = Path(__file__).with_name('casimir_cavity_results.json')
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(q=result['q'], analytic_tail=result['bounds']['final_analytic_truncation'],
        maximum_mode_image_disagreement=result['maximum_mode_image_disagreement'], output=str(output))))
