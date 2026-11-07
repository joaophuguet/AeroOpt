"""Independent exact Joukowski potential, not a panel-method implementation.

Analytic circle potential plus conformal mapping w=z+1/z. Symmetric generating
circle: center=-epsilon, radius=1+epsilon, cusp at z=1. All outputs normalized
by the mapped chord. No external coordinates or experimental data.
"""

import numpy as np


def joukowski_points(panel_count: int, epsilon: float = .08) -> np.ndarray:
    theta = np.linspace(0, 2 * np.pi, panel_count + 1)
    z = -epsilon + (1 + epsilon) * np.exp(1j * theta)
    w = z + 1 / z
    chord = 2 + (1 + 2 * epsilon) + 1 / (1 + 2 * epsilon)
    points = np.column_stack(((w.real - (2 - chord)) / chord, w.imag / chord))
    points[-1] = points[0]  # analytical cusp, explicitly generated
    return points


def joukowski_exact(alpha_deg: float, sample_count: int = 32768,
                    epsilon: float = .08) -> dict:
    """Exact CL and midpoint quadrature of exact surface pressure/moment.

    Independent of production influence coefficients and linear system. A
    midpoint theta grid avoids the removable cusp singularity at theta=0.
    cm is CCW about (.25,0), normalized by q_inf*c_main^2.
    """
    alpha = np.deg2rad(alpha_deg)
    radius = 1 + epsilon
    chord = 2 + (1 + 2 * epsilon) + 1 / (1 + 2 * epsilon)
    theta = (np.arange(sample_count) + .5) * 2 * np.pi / sample_count
    relative = radius * np.exp(1j * theta)
    z = -epsilon + relative
    w = z + 1 / z
    circulation = -4 * np.pi * radius * np.sin(alpha)
    derivative = 1 - 1 / z**2
    complex_velocity = (np.exp(-1j * alpha) - radius**2 * np.exp(1j * alpha) / relative**2
                        + circulation / (2j * np.pi * relative)) / derivative
    cp = 1 - np.abs(complex_velocity)**2
    dw = derivative * 1j * relative * 2 * np.pi / sample_count / chord
    # CCW contour: n*ds=(dy,-dx), so dF=(-Cp*dy,Cp*dx).
    dfx, dfy = -cp * dw.imag, cp * dw.real
    x, y = (w.real - (2 - chord)) / chord, w.imag / chord
    return {"cl": float(-2 * circulation / chord),
            "cm_ccw_quarter_chord": float(np.sum((x - .25) * dfy - y * dfx)),
            "cl_pressure_quadrature": float(np.sum(-np.sin(alpha) * dfx + np.cos(alpha) * dfy)),
            "circulation": float(circulation / chord), "cp": cp}
