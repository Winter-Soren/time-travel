"""Pointwise regular-chart coefficients of the unchanged physical Dirac field.

Coordinates are (tau,u,y1,y2), with x=pole*(1-|y|^2)/(1+|y|^2).
E[a,mu] is the positive coframe; I[a,mu] is its inverse frame (both
matrices are symmetric here). In canonical coordinate L2,
  D = gamma[a] I[a,mu] (partial_mu-i A_mu) + gamma[a] divI[a]/2.
The triple-Clifford connection vanishes for this block coframe. The
physical half-density is v=2*P*R/(1+|y|^2); physical kernels retain both
endpoint factors. Curvature and its derivatives below are PER UNIT q.

The paired probe changes only E00 by exp(epsilon*f) and E11 by
exp(-epsilon*f). The angular coframe and gauge are held fixed, so v is
unchanged. This is a point evaluator, not an interval enclosure or a heat
kernel/source calculation. Geometry conventions: PHYSICAL_EUCLIDEAN_OPERATOR
checkpoint179, and Emparan hep-th/9906160 Eqs8--16.
"""
from mpmath import mp

from dihole_channel_bounds import calibration, field_function


class PointwiseDirac:
    """Compile once; optionally pass radial=(R,Rprime,Rsecond) to point()."""

    def __init__(self, context=mp):
        self.ctx = context
        self.field = field_function(context)
        self._parameters = {}
        self._radial = {}

    def _params(self):
        precision = self.ctx.dps
        if precision not in self._parameters:
            self._parameters[precision] = calibration(self.ctx)
        return self._parameters[precision]

    def step(self, t):
        c = self.ctx
        if t <= 0:
            return c.mpf(0), c.mpf(0), c.mpf(0)
        if t >= 1:
            return c.mpf(1), c.mpf(0), c.mpf(0)
        h = -1/t+1/(1-t)
        s = 1/(1+c.exp(-h))
        hp = 1/t**2+1/(1-t)**2
        hpp = -2/t**3+2/(1-t)**3
        first = s*(1-s)*hp
        second = s*(1-s)*((1-2*s)*hp**2+hpp)
        return s, first, second

    def radial(self, u):
        """Exact smooth definition: constant left, original exterior R=u right."""
        c = self.ctx
        u = c.mpf(u)
        L = c.mpf(1)/20
        if u <= 0:
            return c.mpf(1)/40, c.mpf(0), c.mpf(0)
        if u >= L:
            return u, c.mpf(1), c.mpf(0)
        key = (c.dps, u._mpf_)
        if key not in self._radial:
            t = u/L
            # Reflection keeps quadrature on the non-growing endpoint half.
            upper = min(t, 1-t)
            integral = c.quad(lambda z: self.step(z)[0], [0, upper])
            if t > c.mpf(1)/2:
                integral += t-c.mpf(1)/2
            s, first, _ = self.step(t)
            self._radial[key] = (c.mpf(1)/40+L*integral, s, first/L)
        return self._radial[key]

    def point(self, u, y1, y2, pole=1, radial=None):
        c = self.ctx
        if pole not in (-1, 1):
            raise ValueError('pole must be +1 or -1')
        u, y1, y2 = map(c.mpf, (u, y1, y2))
        t = y1*y1+y2*y2
        if t > 2:
            raise ValueError('The certified chart disk has |y|^2 <= 2')
        R, rate, second = self.radial(u) if radial is None else tuple(map(c.mpf, radial))
        if R <= 0:
            raise ValueError('Positive radius required')
        # Guard removable quotients; linked balanced parameters are recomputed
        # at this precision, rather than dividing a rounded polar residual by t.
        extra = 12 if not t else max(12, 12-int(c.floor(c.log10(t))))
        with c.workdps(c.dps+extra):
            x = pole*(1-t)/(1+t)
            raw = self.field(R, x, *self._params())
            F, W, potential, FR, WR, aR, Fx, Wx, ax, _, _, aRR, _, _, aRx, _, _, axx = raw
            N, P = c.sqrt(F/W), c.sqrt(F*W)
            J = 1/P
            Nx = N*(Fx/F-Wx/W)/2
            Jx = -J*(Fx/F+Wx/W)/2
            JR = -J*(FR/F+WR/W)/2
            Px = P*(Fx/F+Wx/W)/2
            inverseNx = -Nx/N**2
            y = c.matrix([y1, y2])
            if t:
                dI = (J-N)/t
                dE = (P-1/N)/t
                ar_over_t, arr_over_t = aR/t, aRR/t
            else:
                dI = -2*pole*(Jx-Nx)
                dE = -2*pole*(Px-inverseNx)
                ar_over_t = -2*pole*aRx
                arr_over_t = c.mpf(0)  # Multiplied by y=0 in d_u F_{u,y}.
            E, I = c.zeros(4), c.zeros(4)
            E[0,0], E[1,1] = N, P
            I[0,0], I[1,1] = 1/N, J
            for a in range(2):
                for b in range(2):
                    delta = int(a == b)
                    E[a+2,b+2] = 2*R/(1+t)*(delta/N+dE*y[a]*y[b])
                    I[a+2,b+2] = (1+t)/(2*R)*(N*delta+dI*y[a]*y[b])
            angular_div = (N-2*pole*Jx/(1+t)+(1+3*t)*dI/2)/R
            divI = c.matrix([0, rate*JR, angular_div*y1, angular_div*y2])
            curvature = c.zeros(4)
            dF = [c.zeros(4) for _ in range(4)]
            def skew(matrix, a, b, value):
                matrix[a,b], matrix[b,a] = value, -value
            turn = c.matrix([-y2, y1])
            rot = c.matrix([[0,-1],[1,0]])
            radial_curvature = pole*rate*ar_over_t
            for j in range(2):
                skew(curvature, 1, j+2, radial_curvature*turn[j])
                skew(dF[1], 1, j+2,
                     pole*(second*ar_over_t+rate*rate*arr_over_t)*turn[j])
                for k in range(2):
                    # y_k*(J y)_j/t is bounded: this form avoids constructing
                    # a singular derivative of the removable quotient itself.
                    correction = (2*pole*rate*y[k]*turn[j]/t
                                  *(-2*pole*aRx/(1+t)**2-ar_over_t)) if t else 0
                    skew(dF[k+2], 1, j+2, radial_curvature*rot[j,k]+correction)
            skew(curvature, 2, 3, -4*ax/(1+t)**2)
            skew(dF[1], 2, 3, -4*rate*aRx/(1+t)**2)
            for k in range(2):
                skew(dF[k+2], 2, 3,
                     16*y[k]*ax/(1+t)**3+16*pole*y[k]*axx/(1+t)**4)
            _, step_first, step_second = self.step(u/(c.mpf(1)/20))
            f, fu = step_first/2, 10*step_second
            return dict(E=E, I=I, divI=divI, d_u_inverse=rate*JR,
                        f=f, f_u=fu, F=curvature, dF=dF,
                        halfdensity=2*P*R/(1+t), radius=R,
                        Rprime=rate, Rsecond=second, N=N, P=P,
                        x=x, unit_A_phi=potential)


def build_point_evaluator(context=mp):
    return PointwiseDirac(context)


def self_check(context=mp):
    """Focused actual-field checks; no source integral or geometry solve."""
    c = context
    evaluator = PointwiseDirac(c)
    from dihole_uv_bound import gamma_matrices
    gamma = [c.matrix(g.tolist()) for g in gamma_matrices()]
    max_frame = max_div = max_spin = max_curvature = max_variation = c.mpf(0)
    max_potential = max_mask = c.mpf(0)
    points = [(c.mpf('.019'), c.mpf('.23'), c.mpf('-.17'), 1),
              (c.mpf('.037'), c.mpf('-.31'), c.mpf('.12'), -1)]
    for u, y1, y2, pole in points:
        coords = [c.mpf(0), u, y1, y2]
        def at(values): return evaluator.point(*values[1:], pole=pole)
        point = at(coords)
        E, I = point['E'], point['I']
        max_frame = max(max_frame, c.norm(E*I-c.eye(4)))
        derivatives = []
        for k in range(4):
            def changed(z):
                values = list(coords); values[k] = z
                return at(values)
            derivatives.append(c.zeros(4) if k == 0 else c.matrix([
                [c.diff(lambda z: changed(z)['E'][a,b], coords[k]) for b in range(4)] for a in range(4)]))
        dI = [-I*d*I for d in derivatives]
        div = c.matrix([sum(dI[k][a,k] for k in range(4)) for a in range(4)])
        max_div = max(max_div, c.norm(div-point['divI']))
        metric, inverse = E.T*E, I*I.T
        dg = [d.T*E+E.T*d for d in derivatives]
        connection = [[[sum(inverse[a,h]*(dg[b][h,d]+dg[d][h,b]-dg[h][b,d])/2
                            for h in range(4)) for d in range(4)] for b in range(4)] for a in range(4)]
        zeroth = c.zeros(4)
        for mu in range(4):
            omega = c.zeros(4)
            for a in range(4):
                for b in range(4):
                    value = sum(E[a,j]*(dI[mu][j,b]+sum(connection[j][mu][k]*I[k,b]
                                for k in range(4))) for j in range(4))
                    omega += value*(gamma[a]*gamma[b])/4
            density_log = sum((I*derivatives[mu])[j,j] for j in range(4))/2
            cmu = sum((I[a,mu]*gamma[a] for a in range(4)), c.zeros(4))
            zeroth += cmu*(omega-density_log*c.eye(4))
        divergence_term = sum((point['divI'][a]*gamma[a]/2 for a in range(4)), c.zeros(4))
        max_spin = max(max_spin, c.norm(zeroth-divergence_term))
        for k in range(1,4):
            for a,b in ((1,2),(1,3),(2,3)):
                def value(z):
                    values = list(coords); values[k] = z
                    return at(values)['F'][a,b]
                max_curvature = max(max_curvature, abs(c.diff(value,coords[k])-point['dF'][k][a,b]))
        def gauge(values):
            yy1, yy2 = values[2:]
            t = yy1*yy1+yy2*yy2
            scalar = pole*(at(values)['unit_A_phi']-(1-pole)/2)/t
            return c.matrix([0,0,-scalar*yy2,scalar*yy1])
        def potential_derivative(k, component):
            def value(z):
                values = list(coords); values[k] = z
                return gauge(values)[component]
            return c.diff(value,coords[k])
        for a,b in ((1,2),(1,3),(2,3)):
            observed = potential_derivative(a,b)-potential_derivative(b,a)
            max_potential = max(max_potential,abs(observed-point['F'][a,b]))
        max_mask = max(max_mask,abs(c.diff(lambda z:evaluator.point(z,y1,y2,pole)['f'],u)-point['f_u']))
        # Actual finite coframe variation: angular entries and volume fixed.
        def varied(e):
            ans = c.matrix(E)
            ans[0,0] *= c.exp(e*point['f']); ans[1,1] *= c.exp(-e*point['f'])
            return ans
        max_variation = max(max_variation, abs(c.diff(lambda e:c.det(varied(e)),0)))
        max_variation = max(max_variation,
                            abs(c.diff(lambda e:(varied(e)**-1)[1,1],0)-point['f']*I[1,1]))
    # Axis identities include the curvature's nonzero first radial-angular jet.
    pole_error = c.mpf(0)
    for pole in (-1,1):
        point = evaluator.point(c.mpf('.021'),0,0,pole)
        for k in (2,3):
            for a,b in ((1,2),(1,3),(2,3)):
                def value(z):
                    yy = [c.mpf(0),c.mpf(0)]; yy[k-2] = z
                    return evaluator.point(c.mpf('.021'),*yy,pole)['F'][a,b]
                pole_error = max(pole_error,abs(c.diff(value,0)-point['dF'][k][a,b]))
    extension = all(evaluator.point(u,0,0)['f'] == evaluator.point(u,0,0)['f_u'] == 0
                    for u in map(c.mpf,('-.001','0','.05','.051')))
    result = dict(frame_inverse_error=max_frame,analytic_divergence_error=max_div,
                  Cartan_halfdensity_error=max_spin,curvature_derivative_error=max_curvature,
                  curvature_from_potential_error=max_potential,mask_derivative_error=max_mask,
                  polar_curvature_derivative_error=pole_error,paired_variation_error=max_variation,
                  exact_zero_mask_extension=extension)
    threshold = c.power(10, -min(30,c.dps//2))
    result['passed'] = extension and all(v < threshold for k,v in result.items() if k.endswith('_error'))
    if not result['passed']:
        raise ArithmeticError({k:str(v) for k,v in result.items()})
    return result


if __name__ == '__main__':
    mp.dps = 50
    import json
    print(json.dumps({k:(v if isinstance(v,bool) else str(v)) for k,v in self_check().items()},indent=2))
