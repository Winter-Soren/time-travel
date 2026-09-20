"""Positive quadrature on a triangle mesh restricted outside round mouths.

Horizontal slices split at vertices and edge/circle intersections. Curved
radial limits remove the true disks; nodes are never accepted by a mask.
Floating geometry and finite quadrature order are not interval certified.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from two_mouth_mesh import mesh


def rule(vertices, triangles, order=6, b=.2, d=1.):
    """Return points, rho*d_rho*d_z weights, and triangle indices."""
    if (isinstance(order,bool) or not isinstance(order,int) or order<2
            or not all(math.isfinite(v) for v in (b,d)) or not 0<2*b<d):
        raise ValueError('Require integer order>=2 and disjoint positive mouths')
    vertices=np.asarray(vertices,float);triangles=np.asarray(triangles,int)
    if (vertices.ndim!=2 or vertices.shape[1]!=2 or not np.all(np.isfinite(vertices))
            or np.min(vertices[:,0]) < -1e-12 or triangles.ndim!=2 or triangles.shape[1]!=3
            or np.min(triangles)<0 or np.max(triangles)>=len(vertices)):
        raise ValueError('Finite nonnegative-rho vertices and valid triangle indices required')
    nodes,weights=np.polynomial.legendre.leggauss(order)
    nodes=(nodes+1)/2;weights=weights/2
    points=[];volume=[];owners=[]
    for index,tri in enumerate(triangles):
        xy=vertices[tri]
        cuts=list(xy[:,1]);zmin,zmax=min(cuts),max(cuts)
        if zmin==zmax:raise ValueError('Degenerate triangle')
        edges=[(xy[i],xy[j]) for i,j in ((0,1),(1,2),(2,0))]
        for center in (-d/2,d/2):
            cuts.extend(z for z in (center-b,center+b) if zmin<z<zmax)
            for p0,p1 in edges:
                delta=p1-p0;offset=p0-np.array([0.,center])
                aa=float(delta@delta)
                if aa==0:raise ValueError('Degenerate edge')
                bb=2*float(offset@delta);cc=float(offset@offset)-b*b
                disc=bb*bb-4*aa*cc
                if disc<0:continue
                root=math.sqrt(disc)
                for t in ((-bb-root)/(2*aa),(-bb+root)/(2*aa)):
                    if 0<t<1:cuts.append(float(p0[1]+t*delta[1]))
        cuts=sorted(set(cuts))
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            if hi-lo<1e-14:continue
            mid=(lo+hi)/2
            crossing=[]
            for p0,p1 in edges:
                if min(p0[1],p1[1])<mid<max(p0[1],p1[1]):
                    slope=(p1[0]-p0[0])/(p1[1]-p0[1])
                    intercept=p0[0]-slope*p0[1]
                    crossing.append((slope*mid+intercept,slope,intercept))
            if len(crossing)!=2:raise ArithmeticError('Ambiguous triangle slice')
            crossing.sort();left,right=crossing
            center=next((c for c in (-d/2,d/2) if abs(mid-c)<b),None)
            curved=center is not None and math.sqrt(max(0.,b*b-(mid-center)**2))>left[0]
            if curved and math.sqrt(max(0.,b*b-(mid-center)**2))>=right[0]:continue
            zz=lo+(hi-lo)*nodes
            lower=(np.sqrt(np.maximum(0.,b*b-(zz-center)**2)) if curved else left[1]*zz+left[2])
            upper=right[1]*zz+right[2]
            if np.min(upper-lower)<-1e-12:raise ArithmeticError('Circle crossing missing from split')
            width=np.maximum(0.,upper-lower)
            rr=lower[:,None]+width[:,None]*nodes[None,:]
            ww=(hi-lo)*weights[:,None]*weights[None,:]*width[:,None]*rr
            points.append(np.column_stack((rr.ravel(),np.repeat(zz,order))))
            volume.append(ww.ravel());owners.append(np.full(order*order,index,int))
    if not points:raise ValueError('Empty physical integration region')
    return np.vstack(points),np.concatenate(volume),np.concatenate(owners)


def controls():
    geometry=mesh(h=.12)
    vertices,triangles=geometry['vertices'],geometry['triangles']
    # Divergence theorem: F=(rho^2*z^k/2,0), over the outer polygon only.
    nodes,weights=np.polynomial.legendre.leggauss(6)
    nodes=(nodes+1)/2;weights=weights/2
    exact=np.zeros(3)
    for edge in geometry['boundary_edges']:
        if edge['kind']!='outer':continue
        p0,p1=vertices[edge['vertices']]
        at=p0+nodes[:,None]*(p1-p0)
        for k in range(3):exact[k]+=(p1[1]-p0[1])/2*np.sum(weights*at[:,0]**2*at[:,1]**k)
    b,d=.2,1.
    exact-=np.array([4*b**3/3,0.,4*b**3/3*(d*d/4+b*b/5)])
    rows=[]
    for order in (3,6):
        points,weights,owners=rule(vertices,triangles,order)
        radius=points[:,0];z=points[:,1]
        moments=np.array([np.sum(weights*z**k) for k in range(3)])
        distance=np.minimum(np.hypot(radius,z-.5),np.hypot(radius,z+.5))
        tri=vertices[triangles[owners]]
        edge1,edge2=tri[:,1]-tri[:,0],tri[:,2]-tri[:,0];offset=points-tri[:,0]
        det=edge1[:,0]*edge2[:,1]-edge1[:,1]*edge2[:,0]
        alpha=(offset[:,0]*edge2[:,1]-offset[:,1]*edge2[:,0])/det
        beta=(edge1[:,0]*offset[:,1]-edge1[:,1]*offset[:,0])/det
        rows.append(dict(order=order,points=len(points),moments=moments.tolist(),
            maximum_moment_error=float(np.max(abs(moments-exact))),
            minimum_weight=float(min(weights)),minimum_mouth_radius=float(min(distance)),
            minimum_barycentric=float(min(np.min(alpha),np.min(beta),np.min(1-alpha-beta)))))
    tests=dict(physical_volume_moments=max(r['maximum_moment_error'] for r in rows)<1e-11,
        positive_weights=min(r['minimum_weight'] for r in rows)>=0,
        actual_round_holes=min(r['minimum_mouth_radius'] for r in rows)>=b-1e-13,
        original_triangle_ownership=min(r['minimum_barycentric'] for r in rows)>-1e-11)
    return dict(passed=all(tests.values()),tests=tests,geometry=geometry['parameters'],
        exact_physical_moments=exact.tolist(),controls=rows,
        scope='Numerical quadrature on the round-hole exterior inside the outer polygon. No interval bound or certification of a field residual.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={'two_mouth_mesh.py':hashlib.sha256(project_path('two_mouth_mesh.py').read_bytes()).hexdigest()})


if __name__=='__main__':
    result=controls()
    project_path('physical_exterior_quadrature_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2));assert result['passed'],result['tests']
