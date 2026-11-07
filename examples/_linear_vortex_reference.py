"""Review-only independent linear-vortex boundary element reference.

No imports from aeroopt.solvers, no source panels and no constant circulation
unknown. N+1 nodal vortex strengths per sharp closed contour, continuous at
ordinary nodes and independent at the two TE copies. Endpoint-strength Kutta:
gamma_first + gamma_last = 0. Complex logarithmic integrals are derived directly
from the CCW point vortex. This is verification code, not a new public adapter.
"""

import numpy as np


def linear_influence(targets, start, end, *, self_index=None, orientation=1):
    """Complex u-iv coefficients for the two linear shape functions.

    With q=(z-z_start)/t and L=|z_end-z_start|,
    I0=log(q/(q-L))/t, I1=(q*log(q/(q-L))/L-1)/t.
    Velocity=(gamma_start*(I0-I1)+gamma_end*I1)/(2*pi*i).
    The explicit self limit uses log(q/(q-L))=orientation*i*pi.
    """
    start, end = complex(start), complex(end)
    length = abs(end-start)
    if length <= 0:
        raise ValueError("linear reference needs positive panel lengths")
    tangent = (end-start)/length
    q = (np.asarray(targets, dtype=complex) - start) / tangent
    # Do not evaluate a logarithm on the self branch; insert its exterior limit.
    q_safe = q.copy()
    if self_index is not None:
        q_safe[self_index] = .5*length + 1j*length
    logarithm = -np.log1p(-length/q_safe)
    if self_index is not None:
        logarithm[self_index] = orientation * 1j*np.pi
    # q/L*log(q/(q-L))-1 also cancels in the far field. Integrating the
    # geometric series directly gives sum_(k>=1) (L/q)^k/(k+1). The series
    # for I0 also avoids cancellation in complex log1p implementations.
    ratio = length/q_safe
    far = np.abs(ratio) < .05
    logarithm[far] = sum(ratio[far]**k/k for k in range(1,13))
    i0 = logarithm/tangent
    first_moment = q/length*logarithm-1
    first_moment[far] = sum(ratio[far]**k/(k+1) for k in range(1,13))
    i1 = first_moment/tangent
    return (i0-i1)/(2j*np.pi), i1/(2j*np.pi)


def solve_linear_vortex(contours, alpha_deg, reference=(.25, 0.)):
    """Solve sharp closed simple contours in their given global coordinates.

    Float64 dense reference, no resampling and no geometry normalization. The
    input benchmark is generated analytically and independently validated by
    exact potential solutions and off-panel quadrature in tests.
    """
    contours = tuple(np.asarray(points, dtype=float) for points in contours)
    if not contours or any(len(p) < 5 or not np.array_equal(p[0], p[-1]) for p in contours):
        raise ValueError("linear reference supports exactly closed sharp contours only")
    nodes = [p[:,0] + 1j*p[:,1] for p in contours]
    counts = [len(p)-1 for p in nodes]
    starts = np.concatenate([p[:-1] for p in nodes])
    ends = np.concatenate([p[1:] for p in nodes])
    length = np.abs(ends-starts)
    if np.any(length <= 0):
        raise ValueError("linear reference needs positive panel lengths")
    tangent = (ends-starts)/length
    centers = (starts+ends)/2
    signs = [1 if np.sum(p[:-1].real*p[1:].imag-p[1:].real*p[:-1].imag) > 0 else -1 for p in nodes]
    normal = np.concatenate([-1j*s*tangent[sum(counts[:e]):sum(counts[:e+1])] for e,s in enumerate(signs)])
    size = sum(counts)
    influence = np.zeros((size, size+len(nodes)), dtype=complex)
    panel_offset = node_offset = 0
    for e, points in enumerate(nodes):
        for j in range(counts[e]):
            left, right = linear_influence(centers, points[j], points[j+1],
                                           self_index=panel_offset+j, orientation=signs[e])
            influence[:,node_offset+j] += left
            influence[:,node_offset+j+1] += right
        panel_offset += counts[e]
        node_offset += counts[e]+1
    # For w=u-iv, V.n=real(w*n_complex).
    normal_matrix = np.real(influence*normal[:,None])
    tangent_matrix = np.real(influence*tangent[:,None])
    free = np.exp(-1j*np.deg2rad(alpha_deg))
    kutta = np.zeros((len(nodes), size+len(nodes)))
    node_offset = 0
    for e, count in enumerate(counts):
        kutta[e,node_offset] = kutta[e,node_offset+count] = 1
        node_offset += count+1
    matrix = np.vstack((normal_matrix, kutta))
    rhs = np.concatenate((-np.real(free*normal), np.zeros(len(nodes))))
    condition = float(np.linalg.cond(matrix))
    if not np.isfinite(condition) or condition > 1e10:
        raise ArithmeticError(f"independent linear-vortex reference cond_2={condition:.17g} exceeds 1e10")
    strengths = np.linalg.solve(matrix,rhs)
    vt = tangent_matrix@strengths + np.real(free*tangent)
    vn = normal_matrix@strengths + np.real(free*normal)
    if (not np.isfinite(strengths).all() or not np.isfinite(vt).all()
            or np.max(np.abs(vn)) > 1e-9 or np.max(np.abs(kutta@strengths)) > 1e-9):
        raise ArithmeticError("independent linear-vortex reference failed its finite/residual checks")
    cp = 1-vt**2
    df = -cp*normal*length
    force = np.sum(df)
    arm = centers-complex(*reference)
    cm = float(np.imag(np.sum(np.conj(arm)*df)))
    circulation, elements = [], []
    panel_offset = node_offset = 0
    for count in counts:
        sl = slice(panel_offset,panel_offset+count)
        gamma = strengths[node_offset:node_offset+count+1]
        total = float(np.sum(.5*(gamma[:-1]+gamma[1:])*length[sl]))
        circulation.append(total)
        elements.append({"cp": cp[sl].tolist(), "tangential_velocity": vt[sl].tolist(),
                         "nodal_vorticity": gamma.tolist(), "circulation": total,
                         "centers": np.column_stack((centers[sl].real,centers[sl].imag)).tolist()})
        panel_offset += count
        node_offset += count+1
    return {"accepted": True,"cl": float(np.imag(force*free)), "cm_ccw": cm, "circulations": circulation,
            "condition_number_2": condition,
            "cl_circulation": -2*sum(circulation), "normal_residual": float(np.max(np.abs(vn))),
            "kutta_residual": float(np.max(np.abs(kutta@strengths))), "elements": elements}
