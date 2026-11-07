"""Research-only higher-order integration of an already solved panel field.

The no-penetration equations are not re-solved at quadrature points. Their
normal velocity is recorded explicitly. Neither variant is a production
coefficient or an automatically applied numerical correction.
"""

import numpy as np

from aeroopt.solvers._panel_method import source_velocity


def pressure_quadrature(meshes,solution,alpha_deg,order=2,reference=(.25,0.)):
    nodes,weights = np.polynomial.legendre.leggauss(order)
    fractions,weights = (nodes+1)/2,weights/2
    starts = np.vstack([m.starts for m in meshes])
    ends = np.vstack([m.ends for m in meshes])
    normals = np.vstack([m.normals for m in meshes])
    tangents = np.vstack([m.tangents for m in meshes])
    lengths = np.concatenate([m.lengths for m in meshes])
    sources = np.concatenate([e["source_strengths"] for e in solution["elements"]])
    vortex = np.concatenate([np.full(m.size,e["vortex_sheet_strength"]) for m,e in zip(meshes,solution["elements"],strict=True)])
    targets = (starts[:,None,:]+fractions[None,:,None]*(ends-starts)[:,None,:]).reshape(-1,2)
    alpha = np.deg2rad(alpha_deg)
    velocity = np.tile([np.cos(alpha),np.sin(alpha)],(len(targets),1))
    target_panel = np.repeat(np.arange(len(starts)),order)
    for j in range(len(starts)):
        other = target_panel != j
        source = np.empty_like(targets)
        source[other] = source_velocity(targets[other],starts[j:j+1],ends[j:j+1])[:,0]
        # Exact exterior limit at a general interior point of a straight panel.
        radial = np.log(fractions/(1-fractions))/(2*np.pi)
        source[~other] = radial[:,None]*tangents[j]+.5*normals[j]
        rotated = np.column_stack((-source[:,1],source[:,0]))
        velocity += sources[j]*source+vortex[j]*rotated
    n,t = np.repeat(normals,order,axis=0),np.repeat(tangents,order,axis=0)
    vn,vt = np.sum(velocity*n,axis=1),np.sum(velocity*t,axis=1)
    ds = np.repeat(lengths,order)*np.tile(weights,len(lengths))
    arm = targets-reference

    def integrate(cp):
        df = -cp[:,None]*n*ds[:,None]
        force = np.sum(df,axis=0)
        return {"cl": float(force@[-np.sin(alpha),np.cos(alpha)]),
                "cm_ccw": float(np.sum(arm[:,0]*df[:,1]-arm[:,1]*df[:,0]))}

    return {"gauss_order": order,"normal_velocity_max": float(np.max(np.abs(vn))),
            "tangential_cp": integrate(1-vt**2),"total_speed_cp": integrate(1-np.sum(velocity**2,axis=1)),
            "limitation": "research_only_no_impermeability_solution_at_quadrature_points"}
