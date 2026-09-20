"""Conforming polygonal meridional exterior of two finite spherical mouths.

Coordinates are (rho,z). Circular cuts are approximated by chords and refined
with the interior mesh; geometry errors are reported, not hidden in physics.
"""
import math

import numpy as np
from scipy.spatial import ConvexHull, Delaunay


def mesh(b=.2, d=1., R=2., h=.1):
    """Return vertices, CCW triangles and oriented, classified boundary edges."""
    if not np.all(np.isfinite([b, d, R, h])) or not (0 < 2*b < d and R > d/2+b and 0 < h < b):
        raise ValueError('Require 0<2b<d, R>d/2+b, and 0<h<b')
    scale = d/2

    def positive_grid(limit):
        values = [0.]
        while values[-1] < limit:
            values.append(min(limit, values[-1]+min(4*h, h*math.hypot(1., values[-1]/scale))))
        return np.array(values)

    positive = positive_grid(R)
    zs = np.r_[-positive[:0:-1], positive]
    interior = []
    for rho in positive[1:-1]:
        for z in zs[1:-1]:
            r = math.hypot(rho, z)
            # Leave room for the independently sampled boundary curves.
            if r >= R-.23*min(4*h, h*math.hypot(1., r/scale)):
                continue
            if min(math.hypot(rho, z-d/2), math.hypot(rho, z+d/2)) <= b+.23*h:
                continue
            interior.append((rho, z))
    points = list(interior)
    circles = {'upper': (b, d/2), 'lower': (b, -d/2), 'outer': (R, 0.)}
    circle_points = {}
    for name, (radius, center) in circles.items():
        spacing = h if name != 'outer' else min(4*h, h*math.hypot(1., R/scale))
        count = max(5, math.ceil(math.pi*radius/spacing))
        angles = np.linspace(0., math.pi, count+1)
        curve = [(radius*math.sin(t), center+radius*math.cos(t)) for t in angles]
        curve[0] = (0., center+radius)
        curve[-1] = (0., center-radius)
        circle_points[name] = curve
        points.extend(curve)
    intervals = [(-R, -d/2-b), (-d/2+b, d/2-b), (d/2+b, R)]
    for lower, upper in intervals:
        candidates = zs[(zs > lower) & (zs < upper)]
        points.extend((0., z) for z in np.r_[lower, candidates, upper])
    # Shared axis/arc endpoints must be the same mesh vertex.
    vertices = np.unique(np.round(np.asarray(points), 13), axis=0)
    raw = Delaunay(vertices).simplices
    centers = vertices[raw].mean(axis=1)
    outside_holes = np.ones(len(raw), dtype=bool)
    for name in ('lower', 'upper'):
        hull = ConvexHull(np.round(np.asarray(circle_points[name]), 13))
        # Test the polygon actually represented by the boundary vertices.
        # A true-circle test would discard thin but valid exterior triangles.
        signed = centers @ hull.equations[:, :2].T+hull.equations[:, 2]
        outside_holes &= np.max(signed, axis=1) > 1e-12*max(1., R)
    triangles = raw[outside_holes].copy()
    p = vertices[triangles]
    cross = (p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1])-(p[:, 1, 1]-p[:, 0, 1])*(p[:, 2, 0]-p[:, 0, 0])
    triangles[cross < 0] = triangles[cross < 0][:, [0, 2, 1]]
    triangles = triangles[abs(cross) > 1e-13*h*h]
    edge_owners = {}
    for element, triangle in enumerate(triangles):
        for i, j in ((0, 1), (1, 2), (2, 0)):
            key = tuple(sorted((int(triangle[i]), int(triangle[j]))))
            edge_owners.setdefault(key, []).append((element, i, j))
    boundary = []
    for key, owners in edge_owners.items():
        if len(owners) == 2:
            continue
        if len(owners) != 1:
            raise RuntimeError('Nonmanifold triangulation')
        element, i, j = owners[0]
        ids = triangles[element, [i, j]]
        ends = vertices[ids]
        tangent = ends[1]-ends[0]
        length = np.linalg.norm(tangent)
        normal = np.array([tangent[1], -tangent[0]])/length
        if np.max(abs(ends[:, 0])) < 2e-12*max(1., R):
            kind = 'axis'
        else:
            matches = [name for name, (radius, center) in circles.items()
                       if np.max(abs(np.hypot(ends[:, 0], ends[:, 1]-center)-radius)) < 2e-12*max(1., R)]
            if len(matches) != 1:
                raise RuntimeError(f'Unexpected cut edge {ends.tolist()}; refine or repair the conforming mesh')
            kind = matches[0]
            radial = ends.mean(axis=0)-np.array([0., circles[kind][1]])
            expected_sign = 1 if kind == 'outer' else -1
            if expected_sign*(radial @ normal) <= 0:
                raise RuntimeError('Incorrect exterior-domain normal')
        boundary.append(dict(vertices=ids.tolist(), triangle=element, local_edge=[i,j],
                             normal=normal.tolist(), kind=kind))
    p = vertices[triangles]
    areas = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1])-(p[:, 1, 1]-p[:, 0, 1])*(p[:, 2, 0]-p[:, 0, 0]))/2
    exact_area = math.pi*R*R/2-math.pi*b*b
    # The weighted volume is the full 3D volume divided by 2pi.
    exact_weighted_volume = 2*R**3/3-4*b**3/3
    return dict(vertices=vertices, triangles=triangles, boundary_edges=boundary,
                parameters=dict(b=b,d=d,R=R,h=h),
                geometry=dict(elements=len(triangles),vertices=len(vertices),minimum_area=float(areas.min()),
                    area=float(areas.sum()), exact_area=exact_area,
                    area_relative_error=float(abs(areas.sum()/exact_area-1)),
                    weighted_volume=float(np.sum(areas*p[:, :, 0].mean(axis=1))),
                    exact_weighted_volume=exact_weighted_volume,
                    weighted_volume_relative_error=float(abs(np.sum(areas*p[:, :, 0].mean(axis=1))/exact_weighted_volume-1))))


def self_check():
    """Check topology, orientation and refinement of the actual curved cuts."""
    records=[]
    for spacing in (.12,.08,.05):
        data=mesh(h=spacing)
        vertices,triangles=data['vertices'],data['triangles']
        used=np.unique(triangles)
        edges={tuple(sorted((int(t[i]),int(t[j])))) for t in triangles for i,j in ((0,1),(1,2),(2,0))}
        # Mouth half-disks attach to the axis: the meridional domain has one boundary.
        assert len(used)-len(edges)+len(triangles)==1
        assert {e['kind'] for e in data['boundary_edges']}=={'axis','upper','lower','outer'}
        divergence=np.zeros(2)
        for edge in data['boundary_edges']:
            ends=vertices[edge['vertices']]
            divergence+=np.linalg.norm(ends[1]-ends[0])*np.asarray(edge['normal'])
        assert np.linalg.norm(divergence)<1e-12
        records.append(data['geometry'])
    assert records[-1]['area_relative_error']<records[0]['area_relative_error']
    assert records[-1]['weighted_volume_relative_error']<records[0]['weighted_volume_relative_error']
    return records


if __name__=='__main__':
    import json
    print(json.dumps(self_check(),indent=2))
