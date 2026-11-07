"""Explicit research variants for v0.3.1; never automatic input modifications.

The dense-vector assembly below is the pre-review strategy, used only for
equivalence and cost measurements. It is NOT an independent flow reference.
The independent formulation lives in _linear_vortex_reference.py.
"""

import numpy as np

from aeroopt.solvers import _panel_method as pm


def assemble_dense_vectors(meshes, alpha_deg):
    """Original v0.3.0 tensor strategy with the current analytical influences."""
    source, vortex = pm.influences(meshes)
    normal = np.vstack([m.normals for m in meshes])
    tangent = np.vstack([m.tangents for m in meshes])
    count, elements = len(normal),len(meshes)
    velocity = np.empty((count,count+elements,2))
    velocity[:,:count] = source
    offset = 0
    pairs = []
    for e,mesh in enumerate(meshes):
        velocity[:,count+e] = np.sum(vortex[:,offset:offset+mesh.size],axis=1)
        pairs.append((offset,offset+mesh.surface_count-1))
        offset += mesh.size
    vn = np.einsum("ijk,ik->ij",velocity,normal)
    vt = np.einsum("ijk,ik->ij",velocity,tangent)
    angle = np.deg2rad(alpha_deg)
    free = np.array([np.cos(angle),np.sin(angle)])
    matrix = np.vstack([vn,*[vt[a]+vt[b] for a,b in pairs]])
    rhs = np.concatenate((-normal@free,[-((tangent[a]+tangent[b])@free) for a,b in pairs]))
    return pm.PanelSystem(meshes,matrix,rhs,vn,vt,free,tuple(pairs))


def joukowski_sample(count, skew=0., epsilon=.08):
    """New analytic samples theta=beta+skew*sin(beta), |skew|<1.

    Returns points and exact Cp at parameter midpoints. No redistribution of
    imported points. skew<0 explicitly clusters the cusp more strongly.
    """
    beta = np.linspace(0,2*np.pi,count+1)
    theta = beta+skew*np.sin(beta)
    radius = 1+epsilon
    z = -epsilon+radius*np.exp(1j*theta)
    w = z+1/z
    chord = 2+(1+2*epsilon)+1/(1+2*epsilon)
    points = np.column_stack(((w.real-(2-chord))/chord,w.imag/chord))
    points[-1] = points[0]
    midtheta = (theta[:-1]+theta[1:])/2
    relative = radius*np.exp(1j*midtheta)
    zm = -epsilon+relative
    alpha = np.deg2rad(4.)
    gamma = -4*np.pi*radius*np.sin(alpha)
    speed = (np.exp(-1j*alpha)-radius**2*np.exp(1j*alpha)/relative**2
             +gamma/(2j*np.pi*relative))/(1-1/zm**2)
    return points,1-np.abs(speed)**2


def joukowski_blasius(alpha_deg, epsilon=.08, reference=(.25,0.)):
    """Closed-form force/moment from the far-field Laurent coefficients.

    W'(w)=U+A/w+B/w^2+..., A=Gamma/(2*pi*i),
    B=U-R^2*conj(U)+i*Gamma*epsilon/(2*pi).
    M_origin=pi*rho*imag(A^2+2*U*B), CCW positive. Independent of
    the surface-pressure quadrature used in _panel_reference.py.
    """
    alpha = np.deg2rad(alpha_deg)
    radius = 1+epsilon
    chord = 2+(1+2*epsilon)+1/(1+2*epsilon)
    gamma = -4*np.pi*radius*np.sin(alpha)
    u = np.exp(-1j*alpha)
    a = gamma/(2j*np.pi)
    b = u-radius**2*np.conj(u)+1j*gamma*epsilon/(2*np.pi)
    x = 2-chord+reference[0]*chord
    y = reference[1]*chord
    cm = 2*(np.pi*np.imag(a*a+2*u*b)+gamma*(x*np.cos(alpha)+y*np.sin(alpha)))/chord**2
    return {"cl": float(-2*gamma/chord), "cm_ccw": float(cm)}


def integrate_pressure(meshes, cp, alpha_deg, reference=(.25,0.)):
    """Independent complex force arithmetic for the pressure-integration audit."""
    centers = np.vstack([m.centers for m in meshes])
    normals = np.vstack([m.normals for m in meshes])
    lengths = np.concatenate([m.lengths for m in meshes])
    n = normals[:,0]+1j*normals[:,1]
    z = centers[:,0]+1j*centers[:,1]
    df = -np.asarray(cp)*n*lengths
    force = np.sum(df)
    return {"cl": float(np.imag(force*np.exp(-1j*np.deg2rad(alpha_deg)))),
            "cm_ccw": float(np.imag(np.sum(np.conj(z-complex(*reference))*df)))}


def evaluate_strengths(system, strengths):
    """Research diagnostic; deliberately does not declare numerical acceptance."""
    normals = np.vstack([m.normals for m in system.meshes])
    tangents = np.vstack([m.tangents for m in system.meshes])
    vt = system.tangent@strengths+tangents@system.freestream
    vn = system.normal@strengths+normals@system.freestream
    alpha = np.rad2deg(np.arctan2(system.freestream[1],system.freestream[0]))
    result = integrate_pressure(system.meshes,1-vt**2,alpha)
    result.update(normal_residual=float(np.max(np.abs(vn))),
                  terminal_kutta_residual=float(max(abs(vt[a]+vt[b]) for a,b in system.kutta_pairs)),
                  acceptance="diagnostic_only_not_production_solution")
    return result


def extrapolated_kutta(system):
    """Research alternative: linear arclength extrapolation to each TE corner.

    The normal equations are unchanged. The terminal velocity extrapolants,
    rather than terminal collocation values, have opposite signs. This variant
    is compared and NOT adopted by the production adapter.
    """
    a,b = system.matrix.copy(),system.rhs.copy()
    count = sum(m.size for m in system.meshes)
    offset = 0
    for e,m in enumerate(system.meshes):
        first,last = offset,offset+m.surface_count-1
        row = np.zeros_like(b)
        rhs = 0.
        for i,j in ((first,first+1),(last,last-1)):
            li,lj = m.lengths[i-offset],m.lengths[j-offset]
            d1,d2 = li/2,li+lj/2
            w1,w2 = d2/(d2-d1),-d1/(d2-d1)
            row += w1*system.tangent[i]+w2*system.tangent[j]
            rhs -= (w1*m.tangents[i-offset]+w2*m.tangents[j-offset])@system.freestream
        a[count+e],b[count+e] = row,rhs
        offset += m.size
    z = np.linalg.solve(a,b)
    result = evaluate_strengths(system,z)
    result["variant_equation_residual"] = float(np.max(np.abs(a@z-b)))
    return result


def subdivide_cap(mesh, parts):
    """Explicit study variant: identical polygon, only the base segmentation changes."""
    if mesh.size == mesh.surface_count or parts == 1:
        return mesh
    cap = np.linspace(mesh.starts[-1],mesh.ends[-1],parts+1)
    starts = np.vstack((mesh.starts[:-1],cap[:-1]))
    ends = np.vstack((mesh.ends[:-1],cap[1:]))
    d = ends-starts
    length = np.linalg.norm(d,axis=1)
    t = d/length[:,None]
    n = mesh.orientation*np.column_stack((t[:,1],-t[:,0]))
    return pm.PanelMesh(starts,ends,(starts+ends)/2,t,n,length,mesh.orientation,mesh.surface_count)
