"""Validate the two integration premises of the frozen126 Schur certificate.

Exact balanced geometry, Dirichlet A4/N18 central sector, composite validated
Gauss rule. The final conclusion concerns this box, not the natural exterior.
"""
import ast
from fractions import Fraction
import hashlib
import inspect
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

import dihole_interval_arithmetic as ia
from dihole_interval_quadrature import _certified_gauss
from dihole_quadrature_remainder import dyadic, outward
from dihole_continuum_heat_tail import target
from dihole_optical_heat import optical_fields
from project_paths import project_path


def scalar_upper(a):
    return float(np.asarray(a[1]))


def norm_bound(a):
    total=ia.sum_interval(ia.square(a))
    # A sum of squares is nonnegative; outward additions of exact zeros may
    # otherwise produce a negative subnormal lower endpoint.
    return ia.sqrt((np.maximum(total[0],0.),total[1]))


def gamma(count):
    u=ia.point(2.**-53)
    n=ia.mul(ia.point(count),u)
    return ia.div(n,ia.sub(ia.point(1.),n))


def midpoint_radius(a):
    mid=a[0]/2+a[1]/2
    # Radius computed outward around the ACTUAL dyadic midpoint, even if it rounded.
    low=ia.sub(ia.point(mid),ia.point(a[0]))[1]
    high=ia.sub(ia.point(a[1]),ia.point(mid))[1]
    return mid,np.maximum(low,high)


def compile_geometry(iv):
    """Lower a restricted internally generated rational/sqrt expression to intervals."""
    rho,z=sp.symbols('rho z',real=True)
    adapter=SimpleNamespace(mpf=sp.Rational,sqrt=sp.sqrt,log=sp.log)
    a,b,g=optical_fields(adapter)(rho,z)
    raw=target(iv)
    rational=lambda x:sp.Rational(*float(x.a).as_integer_ratio())
    L,vs=(rational(raw[k]) for k in ('L','vstar'))
    E=vs/L
    w=sp.exp(-a)/vs
    mass=sp.exp(-b)*(sp.Rational(797711,2)-512000*g)/E
    wr,wz=sp.diff(w,rho),sp.diff(w,z)
    mr,mz=sp.diff(mass,rho),sp.diff(mass,z)
    expressions=[w,L*(wr+wz)/sp.sqrt(2),L*(wr-wz)/sp.sqrt(2),
        L**2*(sp.diff(wr,rho)+sp.diff(wz,z)),mass,
        L*(mr+mz)/sp.sqrt(2),L*(mr-mz)/sp.sqrt(2)]
    template=sp.lambdify((rho,z),expressions,'mpmath',cse=True)
    source=inspect.getsource(template)
    tree=ast.parse(source)
    def exact(node):
        if isinstance(node,ast.Constant) and isinstance(node.value,int):return Fraction(node.value)
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='mpf':return exact(node.args[0])
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -exact(node.operand)
        if isinstance(node,ast.BinOp):
            left,right=exact(node.left),exact(node.right)
            if isinstance(node.op,ast.Div):return left/right
            if isinstance(node.op,ast.Mult):return left*right
        raise ValueError('Nonconstant or unsupported exponent.')
    class Lower(ast.NodeTransformer):
        def visit_Constant(self,node):
            if not isinstance(node.value,int):raise ValueError('Rounded literal in exact expression.')
            return ast.Call(ast.Name('_number',ast.Load()),[node],[])
        def visit_UnaryOp(self,node):
            if not isinstance(node.op,ast.USub):raise ValueError('Unsupported unary operation.')
            return ast.Call(ast.Name('sub',ast.Load()),[ast.Call(ast.Name('_number',ast.Load()),[ast.Constant(0)],[]),self.visit(node.operand)],[])
        def visit_BinOp(self,node):
            if isinstance(node.op,ast.Pow):
                exponent=exact(node.right)
                return ast.Call(ast.Name('_power',ast.Load()),[self.visit(node.left),ast.Constant(exponent.numerator),ast.Constant(exponent.denominator)],[])
            name={ast.Add:'add',ast.Sub:'sub',ast.Mult:'mul',ast.Div:'div'}.get(type(node.op))
            if name is None:raise ValueError('Unsupported generated operation.')
            return ast.Call(ast.Name(name,ast.Load()),[self.visit(node.left),self.visit(node.right)],[])
        def visit_Call(self,node):
            if not isinstance(node.func,ast.Name):raise ValueError('Unsupported generated call.')
            if node.func.id=='mpf':return self.visit(node.args[0])
            if node.func.id=='sqrt':return ast.Call(ast.Name('sqrt',ast.Load()),[self.visit(node.args[0])],[])
            raise ValueError('Unexpected transcendental call: '+node.func.id)
    numbers={}
    def number(n):
        if n not in numbers:
            d=outward(iv.mpf(n));numbers[n]=(np.asarray(d['lower']),np.asarray(d['upper']))
        return numbers[n]
    def power(a,n,d):
        if d==2:a=ia.sqrt(a)
        elif d!=1:raise ValueError('Unsupported rational power.')
        result=ia.point(1.);base=a;exponent=abs(n)
        while exponent:
            if exponent&1:result=ia.mul(result,base)
            exponent//=2
            if exponent:base=ia.square(base)
        return ia.div(ia.point(1.),result) if n<0 else result
    lowered=ast.fix_missing_locations(Lower().visit(tree))
    namespace={name:getattr(ia,name) for name in ('add','sub','mul','div','sqrt')}
    namespace.update(_number=number,_power=power)
    exec(compile(lowered,'<exact_interval_geometry>','exec'),namespace)
    return namespace[template.__name__],dict(symbolic_source=source,interval_source=ast.unparse(lowered))


def sine_grid(order):
    iv,roots,weights,diagnostics=_certified_gauss(order,70)
    nodes=[];ws=[];basis=[];derivative=[];lap=[]
    for center in (-3,-1,1,3):
        for root,weight in zip(roots,weights):
            x=center+root;nodes.append(outward(x));ws.append(outward(weight))
            bs=[];ds=[];ls=[]
            for n in range(1,19):
                frequency=n*iv.pi/8;phase=frequency*(x+4)
                value=iv.sin(phase)/2
                bs.append(outward(value));ds.append(outward(frequency*iv.cos(phase)/2));ls.append(outward(-frequency**2*value))
            basis.append(bs);derivative.append(ds);lap.append(ls)
    def pack(rows):
        if isinstance(rows[0],dict):return tuple(np.array([r[e] for r in rows]) for e in ('lower','upper'))
        return tuple(np.array([[r[e] for r in row] for row in rows]) for e in ('lower','upper'))
    return dict(nodes=pack(nodes),weights=pack(ws),basis=pack(basis),derivative=pack(derivative),lap=pack(lap),
        diagnostics=diagnostics,order=order)


def lift(a):
    return tuple(np.einsum('pk,ab->pakb',x,np.eye(2)).reshape(len(x),2,648) for x in a)


def rotate(a,axis,sqrt2):
    x,y=(tuple(v[:,s,:] for v in a) for s in (0,1))
    first=ia.add(x,y) if axis=='u' else ia.sub(y,x)
    second=ia.sub(x,y) if axis=='u' else ia.add(x,y)
    first,second=ia.div(first,sqrt2),ia.div(second,sqrt2)
    return tuple(np.stack((first[i],second[i]),axis=1) for i in (0,1))


def j_action(a):
    lo,hi=a
    return np.stack((-hi[:,1,:],lo[:,0,:]),axis=1),np.stack((-lo[:,1,:],hi[:,0,:]),axis=1)


def scale(c,a):
    return ia.mul(tuple(np.asarray(x)[:,None,None] for x in c),a)


def fields(grid,indices,function,raw):
    n=len(grid['nodes'][0]);i,j=indices//n,indices%n
    take=lambda a,idx:tuple(x[idx] for x in a)
    product=lambda a,b:tuple(x.reshape(len(i),324) for x in ia.mul(tuple(x[:,:,None] for x in a),tuple(x[:,None,:] for x in b)))
    Bu,Bv=take(grid['basis'],i),take(grid['basis'],j)
    value=lift(product(Bu,Bv))
    du=lift(product(take(grid['derivative'],i),Bv));dv=lift(product(Bu,take(grid['derivative'],j)))
    lap=lift(ia.add(product(take(grid['lap'],i),Bv),product(Bu,take(grid['lap'],j))))
    u,v=take(grid['nodes'],i),take(grid['nodes'],j)
    rho=ia.add(raw['rho_star'],ia.div(ia.mul(raw['L'],ia.add(u,v)),raw['sqrt2']))
    z=ia.div(ia.mul(raw['L'],ia.sub(u,v)),raw['sqrt2'])
    w,wu,wv,lapw,m,mu,mv=function(rho,z)
    Du=ia.add(scale(w,du),scale(ia.div(wu,ia.point(2.)),value))
    Dv=ia.add(scale(w,dv),scale(ia.div(wv,ia.point(2.)),value))
    first=ia.add(ia.sub(ia.point(0.),ia.add(rotate(Du,'u',raw['sqrt2']),rotate(Dv,'v',raw['sqrt2']))),scale(m,j_action(value)))
    potential=ia.sub(ia.sub(ia.square(m),ia.div(ia.mul(w,lapw),ia.point(2.))),ia.div(ia.add(ia.square(wu),ia.square(wv)),ia.point(4.)))
    second=ia.add(ia.sub(scale(potential,value),scale(ia.square(w),lap)),
        ia.sub(scale(w,j_action(ia.sub(scale(wu,dv),scale(wv,du)))),scale(ia.mul(ia.point(2.),w),ia.add(scale(wu,du),scale(wv,dv)))))
    second=ia.add(second,scale(w,ia.sub(scale(mv,rotate(value,'u',raw['sqrt2'])),scale(mu,rotate(value,'v',raw['sqrt2'])))))
    root_weight=ia.sqrt(ia.mul(take(grid['weights'],i),take(grid['weights'],j)))
    return tuple(tuple(x.reshape(2*len(i),648) for x in scale(root_weight,a)) for a in (value,first,second))


def dot_underflow(inner,rows,columns):
    # Conservative absolute dot-product allowance, including gradual underflow.
    count=2*inner+2
    numerator=ia.mul(ia.point(count),ia.point(np.nextafter(0.,1.)))
    denominator=ia.sub(ia.point(1.),ia.mul(ia.point(count),ia.point(2.**-53)))
    return ia.mul(ia.div(numerator,denominator),ia.sqrt(ia.point(rows*columns)))


def push_sum(levels,matrix):
    """Streaming balanced outward sum of exact dyadic Gram blocks."""
    current=ia.point(matrix);level=0
    while level<len(levels) and levels[level] is not None:
        current=ia.add(levels[level],current);levels[level]=None;level+=1
    if level==len(levels):levels.append(current)
    else:levels[level]=current


def finish_sum(levels):
    value=ia.point(np.zeros((648,648)))
    for level in levels:
        if level is not None:value=ia.add(value,level)
    return value


def gram_block(Z):
    C=Z.T@Z
    if not np.isfinite(C).all():raise FloatingPointError('Nonfinite Gram product.')
    i,j=np.triu_indices(len(C),1);C[j,i]=C[i,j]
    return C


def accumulate(state,Z,delta):
    square_norm=ia.sum_interval(ia.square(ia.point(Z)))
    state['S']=ia.add(state['S'],square_norm)
    state['V']=ia.add(state['V'],ia.square(delta))
    state['rounding']=ia.add(state['rounding'],ia.add(ia.mul(gamma(2*len(Z)+2),square_norm),
        dot_underflow(len(Z),648,648)))
    push_sum(state['levels'],gram_block(Z))


def integration(order,grid,function,raw,M):
    def fresh():return dict(S=ia.point(0.),V=ia.point(0.),rounding=ia.point(0.),levels=[])
    first,second=fresh(),fresh();matrix_norm=norm_bound(ia.point(M))
    n=len(grid['nodes'][0]);started=time.perf_counter()
    for begin in range(0,n*n,128):
        indices=np.arange(begin,min(begin+128,n*n))
        phi,H,A=fields(grid,indices,function,raw)
        p,pr=midpoint_radius(phi);h,hr=midpoint_radius(H);a,ar=midpoint_radius(A)
        product=p@M
        R=a-product
        exact_sub=ia.sub(ia.point(a),ia.point(product))
        sub_error=norm_bound(ia.sub(exact_sub,ia.point(R)))
        product_error=ia.add(ia.mul(gamma(2*648+2),ia.mul(norm_bound(ia.point(p)),matrix_norm)),dot_underflow(648,len(p),648))
        delta=ia.add(ia.add(norm_bound(ia.point(ar)),ia.mul(norm_bound(ia.point(pr)),matrix_norm)),ia.add(product_error,sub_error))
        accumulate(first,h,norm_bound(ia.point(hr)));accumulate(second,R,delta)
        if begin%(128*32)==0:
            print(json.dumps(dict(stage='validated_composite_integrals',completed=len(indices)+begin,total=n*n,
                seconds=time.perf_counter()-started)),flush=True)
    summaries=[];values=[]
    for state in (first,second):
        enclosure=finish_sum(state['levels']);H,radius=midpoint_radius(enclosure)
        sum_error=ia.sum_interval(ia.point(radius),axis=1)[1].max()
        field_error=ia.add(ia.mul(ia.point(2.),ia.sqrt(ia.mul(state['S'],state['V']))),state['V'])
        total=ia.add(ia.add(field_error,state['rounding']),ia.point(sum_error))
        summaries.append(dict(midpoint_squared_frobenius_upper=scalar_upper(state['S']),
            field_radius_squared_frobenius_upper=scalar_upper(state['V']),field_gram_error_upper=scalar_upper(field_error),
            block_dot_rounding_upper=scalar_upper(state['rounding']),sum_rounding_upper=float(sum_error),
            quadrature_sum_operator_error_upper=scalar_upper(total)))
        values.append(H)
    return values,summaries


def controls(iv,function,raw):
    from dihole_spatial_complement import evaluators,coefficients,basis_grid,square_action
    # Independent high-precision differentiation of the original optical fields.
    mp=MPContext();mp.dps=90;source=optical_fields(mp)
    constants=target(iv)
    L,vs=(mp.mpf(float(constants[k].a)) for k in ('L','vstar'));E=vs/L
    errors=[]
    for r0,z0 in ((.6672354660415234,0.),(.61,.03),(.73,-.04)):
        r,z=mp.mpf(r0),mp.mpf(z0)
        w=lambda x,y:mp.exp(-source(x,y)[0])/vs
        mass=lambda x,y:mp.exp(-source(x,y)[1])*(mp.mpf(797711)/2-512000*source(x,y)[2])/E
        wr,wz=mp.diff(w,(r,z),(1,0)),mp.diff(w,(r,z),(0,1))
        mr,mz=mp.diff(mass,(r,z),(1,0)),mp.diff(mass,(r,z),(0,1))
        truth=[w(r,z),L*(wr+wz)/mp.sqrt(2),L*(wr-wz)/mp.sqrt(2),
            L*L*(mp.diff(w,(r,z),(2,0))+mp.diff(w,(r,z),(0,2))),mass(r,z),L*(mr+mz)/mp.sqrt(2),L*(mr-mz)/mp.sqrt(2)]
        got=function(ia.point(np.array([r0])),ia.point(np.array([z0])))
        assert all(mp.mpf(float(a[0][0]))<=t<=mp.mpf(float(a[1][0])) for a,t in zip(got,truth))
        errors.append(max(float(a[1][0]-a[0][0]) for a in got))
    grid=sine_grid(8);indices=np.arange(0,len(grid['nodes'][0])**2,31)
    phi,H,A=fields(grid,indices,function,raw)
    assert all(np.all(x[0]<=x[1]) for x in (phi,H,A))
    # Check against a separately implemented second-order action on the same point trials.
    _,floating,_=evaluators(iv);coord={k:float(v.a) for k,v in constants.items()}
    n=len(grid['nodes'][0]);i,j=indices//n,indices%n
    from dihole_finite_heat import sine_table
    nodes=(grid['nodes'][0]+grid['nodes'][1])/2;u,v=nodes[i],nodes[j]
    bu,du=sine_table(u,18,4);bv,dv=sine_table(v,18,4)
    product=lambda x,y:np.einsum('pi,pj->pij',x,y).reshape(len(x),324)
    spin=lambda x:np.einsum('pk,ab->pakb',x,np.eye(2)).reshape(len(x),2,648)
    value=product(bu,bv);freq=np.arange(1,19)*np.pi/8
    baseline=square_action(spin(value),spin(product(du,bv)),spin(product(bu,dv)),
        spin(-value*(freq[:,None]**2+freq[None,:]**2).ravel()),coefficients(floating,coord,u,v,0))
    weight=np.sqrt(((grid['weights'][0]+grid['weights'][1])/2)[i]*((grid['weights'][0]+grid['weights'][1])/2)[j])
    baseline=(baseline*weight[:,None,None]).reshape(-1,648)
    difference=float(np.max(np.maximum(A[0]-baseline,baseline-A[1])))
    assert difference<2e-9,difference  # Independent floating control, not an enclosure premise.
    # Check the midpoint Gram bound against exact rational arithmetic on a tiny case.
    Z=np.array([[.1,.25],[-.3,.7],[.125,-.2]])
    got=Z.T@Z
    error=max(abs(Fraction(float(got[i,j]))-sum((Fraction(float(row[i]))*Fraction(float(row[j])) for row in Z),Fraction())) for i in range(2) for j in range(2))
    bound=scalar_upper(ia.mul(gamma(2*len(Z)+2),ia.sum_interval(ia.square(ia.point(Z)))))
    assert float(error)<bound
    # Runtime probe for gradual underflow in the selected binary64 dot path.
    tiny=np.float64(2.**-1022)
    assert (np.array([[tiny]])@np.array([[.5]]))[0,0]==tiny/2>0
    return dict(passed=True,high_precision_geometry_max_widths=errors,
        independent_square_max_outside_interval=difference,exact_rational_dot_error=float(error),
        exact_rational_dot_allowance=bound,gradual_underflow_probe=True)


def run():
    from scipy.linalg import eigh
    from dihole_sine_remainder import remainder_bounds
    started=time.perf_counter();iv=MPIntervalContext();iv.dps=70
    previous_path=project_path('dihole_schur_spectrum_results.json')
    assert hashlib.sha256(previous_path.read_bytes()).hexdigest()=='7002975a714c32d50a028dbee17dc8892aaab397b6d161d3d2c90d43ccee303c'
    previous=json.loads(previous_path.read_text())
    assert previous['conditional_matrix_certificate']['passed']
    for name,sha in previous['input_sha256'].items():assert hashlib.sha256(project_path(name).read_bytes()).hexdigest()==sha,name
    assert hashlib.sha256(project_path('dihole_schur_spectrum.py').read_bytes()).hexdigest()==previous['source_sha256']
    assert hashlib.sha256(project_path('dihole_schur_spectrum_state.npz').read_bytes()).hexdigest()==previous['state_sha256']
    with np.load(project_path('dihole_heat_angular_state.npz')) as s:M=s['M0'].copy()
    with np.load(project_path('dihole_schur_spectrum_state.npz')) as s:G=s['j0_gram'].copy()
    assert np.array_equal(M,M.T) and np.array_equal(G,G.T)
    function,generated=compile_geometry(iv);raw={}
    for key,value in dict(target(iv),sqrt2=iv.sqrt(2)).items():
        o=outward(value);raw[key]=(np.asarray(o['lower']),np.asarray(o['upper']))
    checks=controls(iv,function,raw)
    remainder=remainder_bounds();order=remainder['quadrature']['order_per_panel']
    grid=sine_grid(order)
    values,summaries=integration(order,grid,function,raw,M)
    comparisons=[]
    for label,value,frozen,summary,rem,budget in zip(('M','residual_Gram'),values,(M,G),summaries,
        (remainder['operator_remainders']['M']['upper'],remainder['operator_remainders']['residual_Gram']['upper']),
        (previous['conditional_matrix_certificate']['required_M_operator_error'],previous['conditional_matrix_certificate']['required_G_upper_error'])):
        difference=norm_bound(ia.sub(ia.point(value),ia.point(frozen)))
        total=ia.add(ia.add(difference,ia.point(summary['quadrature_sum_operator_error_upper'])),ia.point(rem))
        comparisons.append(dict(quantity=label,midpoint_to_frozen_frobenius_upper=scalar_upper(difference),
            analytic_quadrature_remainder_upper=rem,total_operator_error_upper=scalar_upper(total),
            required_error_budget=budget,budget_met=bool(total[1]<budget),**summary))
    passed=all(c['budget_met'] for c in comparisons)
    # A numerical eigenvector is just a dyadic trial for the rigorous Rayleigh upper bound.
    trial=eigh(M,subset_by_index=(0,0),check_finite=False)[1][:,0]
    vt=ia.point(trial[None,:]);vc=ia.point(trial[:,None])
    numerator=ia.matmul(vt,ia.matmul(ia.point(M),vc))
    denominator=ia.matmul(vt,vc)
    upper=ia.add(ia.div(numerator,denominator),ia.point(comparisons[0]['total_operator_error_upper']))
    spectral=dict(lower=previous['conditional_matrix_certificate']['candidate_lower'] if passed else None,
        rayleigh_upper=float(upper[1][0,0]),trial_norm_squared_interval=[float(denominator[0][0,0]),float(denominator[1][0,0])],
        certified=passed,meaning='Lowest dimensionless squared-energy eigenvalue of the exact central angular Dirichlet-box Friedrichs square, subject to the documented arithmetic model.')
    state=project_path('dihole_sine_certificate_state.npz')
    np.savez_compressed(state,M_midpoint=values[0],residual_Gram_midpoint=values[1],rayleigh_trial=trial)
    blas=np.__config__.CONFIG['Build Dependencies']['blas']
    inputs=['dihole_sine_remainder.py','dihole_schur_spectrum.py','dihole_schur_spectrum_results.json',
        'dihole_schur_spectrum_state.npz','dihole_heat_angular_state.npz','dihole_interval_arithmetic.py',
        'dihole_interval_quadrature.py','dihole_quadrature_remainder.py','dihole_continuum_heat_tail.py',
        'dihole_optical_heat.py','dihole_spatial_complement.py','dihole_finite_heat.py',
        'dihole_probe_certificate_results.json','project_paths.py']
    result=dict(passed=passed,controls=checks,quadrature=dict(order_per_panel=order,panel_edges=[-4,-2,0,2,4],
        points_per_axis=len(grid['nodes'][0]),total_spatial_points=len(grid['nodes'][0])**2,
        scalar_rows_per_block=256,validated_rule=grid['diagnostics']),
        analytic_remainder=remainder,error_budgets=comparisons,spectrum=spectral,
        arithmetic_assumptions=['The existing outward IEEE754 interval-ufunc contract, correct mpmath interval transcendental evaluation and no overflow.',
            'The identified OpenBLAS binary64 DGEMM uses conventional double product/add or FMA dot arithmetic, without reduced precision or nonclassical matrix algorithms. Standard gamma_(2k+2) error allowances cover these operations.',
            'Gradual underflow is assumed and probed at runtime; explicit absolute underflow allowances accompany dot bounds. The implementation is not a formal verification of BLAS, mpmath or hardware.'],
        blas={k:blas[k] for k in ('name','version','openblas configuration')},
        method='Validated fields and positive square-root weights; Frobenius field-radius bounds; 256-row midpoint Gram blocks with explicit dot-rounding allowance; outward balanced Gram summation; analytic matrix quadrature remainder.',
        scope=['A single angular sector of the exact A4 Dirichlet box, with the entire spatial complement covered by the Schur theorem.',
            'No finite-basis convergence difference is used as an error certificate.',
            'Other spectral levels/sectors, pointwise heat, exterior boundary effects, renormalized stress, supporting material and stable backward passage remain unresolved.'],
        generated_geometry=generated,generated_sha256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in generated.items()},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),state_sha256=hashlib.sha256(state.read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in inputs},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_sine_certificate_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=passed,seconds=result['runtime_seconds'],error_budgets=comparisons,spectrum=spectral),indent=2),flush=True)
    if not passed:raise ArithmeticError('An integration premise failed its required budget; saved diagnostic is not a spectral certificate.')
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.5,4.5))
    x=np.arange(2);rows=result['error_budgets']
    ax.bar(x-.18,[r['total_operator_error_upper'] for r in rows],.36,label='Enclosed total error')
    ax.bar(x+.18,[r['required_error_budget'] for r in rows],.36,label='Allowed error')
    ax.set_yscale('log');ax.set_xticks(x,['Retained matrix','Full residual integral'])
    ax.set(ylabel='Operator-norm error bound',title='Checking the two integration premises')
    ax.legend();ax.grid(axis='y',alpha=.2)
    if result['passed']:
        text=f"Central Dirichlet squared-energy level: {result['spectrum']['lower']:.5f} to {result['spectrum']['rayleigh_upper']:.5f}"
    else:text='At least one premise remains unproved: inspect the saved error budget.'
    fig.text(.5,.045,text,ha='center',fontsize=10)
    fig.text(.5,.008,'Arithmetic assumptions documented. This is not a physical stability or time-travel result.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,1))
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_sine_certificate_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
