"""Plane-stress annular-rim compliance for an internal ring tooth.

The model represents the material from the internal tooth-root circle to the
housing-integral outer rim.  Through holes are explicitly removed from the
annular mesh; the outer circular boundary is fixed.  The reported output is
the azimuthal mean normal compliance of all 84 tooth-root load locations.
It is an engineering 2-D rim model, not a substitute for a 3-D CAD/FE model
with an exact tooth-root fillet.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


@dataclass(frozen=True)
class RingRimResult:
    compliance_mean: float
    compliance_min: float
    compliance_max: float
    stiffness_mean: float
    n_nodes: int
    n_triangles: int
    nr: int
    ntheta: int


def _triangle_stiffness(xy: np.ndarray, young: float, poisson: float, thickness: float) -> np.ndarray:
    x1, y1 = xy[0]; x2, y2 = xy[1]; x3, y3 = xy[2]
    twice_area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    area = abs(twice_area)/2
    if area <= 0:
        raise ValueError("Degenerate annular-rim triangle")
    b1, b2, b3 = y2-y3, y3-y1, y1-y2
    c1, c2, c3 = x3-x2, x1-x3, x2-x1
    bmat = np.array([[b1, 0, b2, 0, b3, 0], [0, c1, 0, c2, 0, c3],
                     [c1, b1, c2, b2, c3, b3]], dtype=float) / twice_area
    dmat = young/(1-poisson**2) * np.array([[1, poisson, 0], [poisson, 1, 0], [0, 0, (1-poisson)/2]])
    return thickness*area*(bmat.T @ dmat @ bmat)


def annular_rim_normal_compliance(*, root_radius: float, outer_radius: float,
                                  face_width: float, young: float, poisson: float,
                                  hole_radius: float, hole_diameter: float,
                                  hole_count: int, tooth_count: int,
                                  root_tooth_thickness: float, pressure_angle: float,
                                  nr: int = 20, ntheta: int = 336) -> RingRimResult:
    """Return mean tooth-root normal compliance from a perforated annular FE mesh."""
    if ntheta % tooth_count:
        raise ValueError("ntheta must be divisible by tooth_count")
    radii = np.linspace(root_radius, outer_radius, nr+1)
    theta = 2*math.pi*np.arange(ntheta)/ntheta
    rr, tt = np.meshgrid(radii, theta, indexing="ij")
    xy_all = np.column_stack(((rr*np.cos(tt)).ravel(), (rr*np.sin(tt)).ravel()))
    hole_angles = math.radians(22.5) + 2*math.pi*np.arange(hole_count)/hole_count
    hole_centres = np.column_stack((hole_radius*np.cos(hole_angles), hole_radius*np.sin(hole_angles)))
    hole_r2 = (hole_diameter/2)**2
    triangles: list[tuple[int, int, int]] = []
    def node(ir: int, it: int) -> int:
        return ir*ntheta + (it % ntheta)
    def in_hole(point: np.ndarray) -> bool:
        return bool(np.any(np.sum((hole_centres-point)**2, axis=1) < hole_r2))
    for ir in range(nr):
        for it in range(ntheta):
            q00, q01, q10, q11 = node(ir, it), node(ir, it+1), node(ir+1, it), node(ir+1, it+1)
            for tri in ((q00, q10, q11), (q00, q11, q01)):
                if not in_hole(xy_all[list(tri)].mean(axis=0)):
                    triangles.append(tri)
    used = np.unique(np.asarray(triangles, dtype=int).ravel())
    new_index = -np.ones(len(xy_all), dtype=int)
    new_index[used] = np.arange(len(used))
    xy = xy_all[used]
    tri = np.asarray([[new_index[n] for n in item] for item in triangles], dtype=int)
    rows: list[int] = []; cols: list[int] = []; data: list[float] = []
    for element in tri:
        ke = _triangle_stiffness(xy[element], young, poisson, face_width)
        dof = np.ravel(np.column_stack((2*element, 2*element+1)))
        rows.extend(np.repeat(dof, 6)); cols.extend(np.tile(dof, 6)); data.extend(ke.ravel())
    ndof = 2*len(xy)
    stiffness = coo_matrix((data, (rows, cols)), shape=(ndof, ndof)).tocsc()
    r_node = np.hypot(xy[:, 0], xy[:, 1])
    fixed_nodes = np.flatnonzero(np.isclose(r_node, outer_radius, rtol=0, atol=1e-11))
    fixed = np.ravel(np.column_stack((2*fixed_nodes, 2*fixed_nodes+1)))
    free_mask = np.ones(ndof, dtype=bool); free_mask[fixed] = False
    free = np.flatnonzero(free_mask)
    lu = splu(stiffness[free][:, free])
    root_nodes = np.flatnonzero(np.isclose(r_node, root_radius, rtol=0, atol=1e-11))
    root_theta = np.mod(np.arctan2(xy[root_nodes, 1], xy[root_nodes, 0]), 2*math.pi)
    tooth_centres = 2*math.pi*np.arange(tooth_count)/tooth_count
    rhs = np.zeros((len(free), tooth_count))
    half_width = root_tooth_thickness/(2*root_radius)
    dtheta = 2*math.pi/ntheta
    for tooth, centre in enumerate(tooth_centres):
        delta = np.abs(np.angle(np.exp(1j*(root_theta-centre))))
        loaded = root_nodes[delta <= half_width + dtheta/2]
        if len(loaded) == 0:
            raise RuntimeError("No inner-boundary nodes selected for tooth load")
        weights = np.full(len(loaded), 1/len(loaded))
        angles = np.arctan2(xy[loaded, 1], xy[loaded, 0])
        # Unit normal load: radial and tangential components at the tooth root.
        fx = weights*(math.cos(pressure_angle)*np.cos(angles) - math.sin(pressure_angle)*np.sin(angles))
        fy = weights*(math.cos(pressure_angle)*np.sin(angles) + math.sin(pressure_angle)*np.cos(angles))
        global_dof = np.ravel(np.column_stack((2*loaded, 2*loaded+1)))
        force = np.empty(2*len(loaded)); force[0::2] = fx; force[1::2] = fy
        positions = np.searchsorted(free, global_dof)
        keep = (positions < len(free)) & (free[np.minimum(positions, len(free)-1)] == global_dof)
        rhs[positions[keep], tooth] = force[keep]
    displacement = lu.solve(rhs)
    compliances = np.empty(tooth_count)
    for tooth, centre in enumerate(tooth_centres):
        delta = np.abs(np.angle(np.exp(1j*(root_theta-centre))))
        loaded = root_nodes[delta <= half_width + dtheta/2]
        weights = np.full(len(loaded), 1/len(loaded))
        angles = np.arctan2(xy[loaded, 1], xy[loaded, 0])
        fx = weights*(math.cos(pressure_angle)*np.cos(angles) - math.sin(pressure_angle)*np.sin(angles))
        fy = weights*(math.cos(pressure_angle)*np.sin(angles) + math.sin(pressure_angle)*np.cos(angles))
        dof = np.ravel(np.column_stack((2*loaded, 2*loaded+1)))
        pos = np.searchsorted(free, dof)
        u = np.zeros_like(dof, dtype=float)
        ok = (pos < len(free)) & (free[np.minimum(pos, len(free)-1)] == dof)
        u[ok] = displacement[pos[ok], tooth]
        force = np.empty(2*len(loaded)); force[0::2] = fx; force[1::2] = fy
        compliances[tooth] = float(force @ u)
    if not np.all(np.isfinite(compliances)) or np.any(compliances <= 0):
        raise RuntimeError("Invalid annular-rim compliance")
    return RingRimResult(float(compliances.mean()), float(compliances.min()), float(compliances.max()),
                         float(1/compliances.mean()), len(xy), len(tri), nr, ntheta)
