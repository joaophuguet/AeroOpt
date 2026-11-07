"""Independent two-circle image solution and conformal-map benchmark research.

Derived from Williams, ARC R&M 3717, sections 2--3 and appendices I--IV.
No import from the AeroOpt solver kernels. Gamma is CCW-positive here.
"""

from dataclasses import dataclass

import numpy as np

SOURCE_URL = "https://reports.aerade.cranfield.ac.uk/bitstream/handle/1826.2/2993/arc-rm-3717.pdf?isAllowed=y&sequence=1"
SOURCE_SHA256 = "db8ae04e423bd17becf7a1615db477246a576d4dd8300bc6f1dcfdd3c5b9a369"


@dataclass(frozen=True)
class TwoCircleFlow:
    centers: tuple[complex, complex] = (0j,1.9-.4j)
    radii: tuple[float, float] = (1.096,.5)
    reflections: int = 64

    def __post_init__(self):
        if (len(self.centers) != 2 or len(self.radii) != 2
                or not np.isfinite(self.centers).all() or not np.isfinite(self.radii).all()):
            raise ValueError("two finite circle centers/radii required")
        if abs(self.centers[1]-self.centers[0]) <= sum(self.radii):
            raise ValueError("reference circles must be disjoint")
        if min(self.radii) <= 0 or type(self.reflections) is not int or self.reflections < 2:
            raise ValueError("positive radii and at least two reflections required")

    def inverse(self,z,body):
        c,r = self.centers[body],self.radii[body]
        return c+r*r/np.conj(z-c)

    def stream_velocity(self,z,alpha_rad):
        z = np.asarray(z,dtype=complex)
        free = np.exp(-1j*alpha_rad)
        velocity = np.full(z.shape,free,dtype=complex)
        for body in (0,1):
            pole = self.centers[body]
            strength = self.radii[body]**2*np.conj(free)
            velocity -= strength/(z-pole)**2
            for step in range(self.reflections):
                target = 1-body
                c,r = self.centers[target],self.radii[target]
                strength = -r*r*np.conj(strength)/np.conj(pole-c)**2
                pole = self.inverse(pole,target)
                velocity -= strength/(z-pole)**2
                body = target
        return velocity

    def vortex_velocity(self,z,body):
        z = np.asarray(z,dtype=complex)
        # Unit CCW circulation, with zero net circulation around the other body.
        velocity = 1/(2j*np.pi*(z-self.centers[body]))
        target = 1-body
        positive = self.centers[target]
        negative = self.inverse(self.centers[body],target)
        for step in range(self.reflections):
            # Group opposite images before summation (Appendix IV), avoiding
            # cancellation of individual non-vanishing vortex terms.
            velocity += (positive-negative)/(2j*np.pi*(z-positive)*(z-negative))
            target = 1-target
            positive,negative = self.inverse(negative,target),self.inverse(positive,target)
        return velocity

    def circulations(self,trailing_points,alpha_rad):
        te = np.asarray(trailing_points,dtype=complex)
        if te.shape != (2,) or not np.allclose(abs(te-self.centers),self.radii,rtol=1e-12,atol=0):
            raise ValueError("each prescribed stagnation point must lie on its circle")
        tangent = 1j*(te-np.asarray(self.centers))/np.asarray(self.radii)
        rhs = -np.real(self.stream_velocity(te,alpha_rad)*tangent)
        matrix = np.column_stack([np.real(self.vortex_velocity(te,e)*tangent) for e in (0,1)])
        return np.linalg.solve(matrix,rhs)

    def velocity(self,z,alpha_rad,circulations):
        result = self.stream_velocity(z,alpha_rad)
        for body,gamma in enumerate(circulations):
            result += gamma*self.vortex_velocity(z,body)
        return result


def karman_trefftz(z,c,n):
    """Principal logarithm, cut between -c and +c; derivative at both criticals zero."""
    z = np.asarray(z,dtype=complex)
    special = (z == c) | (z == -c)
    safe = np.where(special,2*c+1j*c,z)
    power = np.exp(n*np.log((safe-c)/(safe+c)))
    value = n*c*(1+power)/(1-power)
    derivative = 4*n*n*c*c*power/((safe*safe-c*c)*(1-power)**2)
    value = np.where(z == c,n*c,value)
    value = np.where(z == -c,-n*c,value)
    derivative = np.where(special,0j,derivative)
    return value,derivative


@dataclass(frozen=True)
class PrintedParameterCase:
    """Fully defined two-body analytic case inspired by the printed parameters.

    The report does not explicitly prescribe the intermediate distorted-circle
    radius/centre construction. We use |K1(T2)-K1(C2)| for that radius, retain
    the declared second-map direction, and place its critical point at K1(T2).
    This choice is NEVER fitted to published coefficients or coordinates.
    It closely matches, but is not certified as, the published airfoil A/B.
    The two-circle part does reproduce Table 1. See VALIDATION_V032.md.
    """
    deflection_deg: float = 30.
    alpha_raw_deg: float = 0.
    reflections: int = 64

    def __post_init__(self):
        if not np.isfinite([self.deflection_deg,self.alpha_raw_deg]).all() or not 10 <= self.deflection_deg <= 40:
            raise ValueError("finite inputs and mapping direction in [10,40] degrees required")
        self.flow  # Validate the circle/image-series configuration.

    @property
    def n(self):
        return 2-.17279/np.pi

    @property
    def centers(self):
        first = 1-1.096*np.exp(-.001j)
        return first,first+1.9-.4j

    @property
    def trailing_points(self):
        eta = np.deg2rad(self.deflection_deg)
        return 1+0j,self.centers[1]+.5*np.exp(-1j*(.08725+eta))

    @property
    def second_map(self):
        eta = np.deg2rad(self.deflection_deg)
        te = karman_trefftz(self.trailing_points[1],1.,self.n)[0]
        center = karman_trefftz(self.centers[1],1.,self.n)[0]
        scale = float(abs(te-center)/1.15)
        return complex(te-scale*np.exp(-1j*eta)),scale

    @property
    def flow(self):
        return TwoCircleFlow(self.centers,(1.096,.5),self.reflections)

    @property
    def gamma_raw(self):
        return self.flow.circulations(self.trailing_points,np.deg2rad(self.alpha_raw_deg))

    def map_raw(self,z):
        z = np.asarray(z,dtype=complex)
        first,derivative1 = karman_trefftz(z,1.,self.n)
        origin,c = self.second_map
        rotation = np.exp(1j*np.deg2rad(self.deflection_deg))
        argument = (first-origin)*rotation
        argument = np.where(z == self.trailing_points[1],c+0j,argument)
        second,derivative2 = karman_trefftz(argument,c,self.n)
        return origin+second/rotation,derivative1*derivative2

    @property
    def normalization(self):
        # Nominal KT chord reference: the left intersection of the original
        # main circle with the first-map real axis. Not a sampled minimum x.
        le_z = self.centers[0]-1.096*np.exp(.001j)
        le = complex(self.map_raw(le_z)[0])
        te = complex(self.map_raw(self.trailing_points[0])[0])
        return le,te-le

    @property
    def alpha_deg(self):
        return float(self.alpha_raw_deg-np.rad2deg(np.angle(self.normalization[1])))

    def surface(self,body,theta):
        """theta is displacement from the declared TE, counterclockwise."""
        theta = np.asarray(theta,dtype=float)
        c,r = self.centers[body],(1.096,.5)[body]
        te = self.trailing_points[body]
        z = c+(te-c)*np.exp(1j*theta)
        endpoint = (theta == 0) | (theta == 2*np.pi)
        z = np.where(endpoint,te,z)
        if body == 0:
            z = np.where(theta == np.pi+.002,self.centers[0]-1.096*np.exp(.001j),z)
        mapped,derivative = self.map_raw(z)
        le,chord = self.normalization
        points = (mapped-le)/chord
        cp = np.ones(theta.shape)
        regular = ~endpoint
        velocity = self.flow.velocity(z[regular],np.deg2rad(self.alpha_raw_deg),self.gamma_raw)
        cp[regular] = 1-np.abs(velocity/derivative[regular])**2
        return points,cp,derivative*1j*(z-c)/chord

    def parameter_nodes(self,body,panels_per_element):
        if type(panels_per_element) is not int or panels_per_element < 8 or panels_per_element % 2:
            raise ValueError("even panel count >=8 required")
        if body == 0:
            # Include the explicit nominal leading reference at index N/2.
            split = np.pi+.002
            return np.concatenate((np.linspace(0,split,panels_per_element//2+1),
                                   np.linspace(split,2*np.pi,panels_per_element//2+1)[1:]))
        return np.linspace(0,2*np.pi,panels_per_element+1)

    def contours(self,panels_per_element):
        result = []
        for body in (0,1):
            theta = self.parameter_nodes(body,panels_per_element)
            z,_,_ = self.surface(body,theta)
            points = np.column_stack((z.real,z.imag))
            points[-1] = points[0]
            result.append(points)
        return tuple(result)

    def published_frame(self,z):
        """Explicit x-chord frame for checking the printed coordinate tables.

        Unlike the AeroOpt frame it retains the small nonzero main TE y. No
        scale/rotation is fitted to table values. The disagreement is reported.
        """
        le,chord = self.normalization
        return (self.map_raw(z)[0]-le)/chord.real

    def reference(self,quadrature_points=16384):
        """Midpoint integration of exact analytic Cp plus independent Blasius.

        Coefficients use the explicit main chord and CCW quarter-chord moment.
        These are exact-flow approximations for THIS defined mapping, not
        inferred published values. Truncation/quadrature must be checked.
        """
        if type(quadrature_points) is not int or quadrature_points < 32:
            raise ValueError("integer quadrature count >=32 required")
        theta = (np.arange(quadrature_points)+.5)*2*np.pi/quadrature_points
        forces,moments = [],[]
        for body in (0,1):
            points,cp,tangent = self.surface(body,theta)
            df = -cp*(-1j*tangent)*2*np.pi/quadrature_points
            forces.append(complex(np.sum(df)))
            moments.append(float(np.imag(np.sum(np.conj(points-.25)*df))))
        force = sum(forces)
        alpha = np.deg2rad(self.alpha_deg)
        le,chord = self.normalization
        gamma = self.gamma_raw/abs(chord)
        # A far contour in the original circle plane. Blasius is independent
        # of the surface pressure quadrature and is insensitive to TE sampling.
        t = np.arange(1024)*2*np.pi/1024
        z = 8*np.exp(1j*t)
        mapped,derivative = self.map_raw(z)
        velocity = self.flow.velocity(z,np.deg2rad(self.alpha_raw_deg),self.gamma_raw)
        integrand = velocity**2/derivative*1j*z*2*np.pi/1024
        complex_force = 1j*np.sum(integrand)/abs(chord)
        cm_blasius = -float(np.real(np.sum((mapped-(le+.25*chord))*integrand)))/abs(chord)**2
        return {"cl": float(np.imag(force*np.exp(-1j*alpha))),"cm_ccw": float(sum(moments)),
                "cl_circulation": float(-2*sum(gamma)),"cm_blasius": cm_blasius,
                "cl_blasius": -float(np.imag(complex_force*np.exp(1j*np.deg2rad(self.alpha_raw_deg)))),
                "element_circulations_ccw": gamma.tolist(),
                "element_cl": [float(np.imag(f*np.exp(-1j*alpha))) for f in forces],
                "alpha_deg": self.alpha_deg,"quadrature_points": quadrature_points,
                "reference_kind": "analytic_two_body_declared_mapping",
                "published_airfoil_reproduction_certified": False}

    def flow_checks(self):
        """Independent image truncation, circle boundary and stagnation checks."""
        theta = (np.arange(512)+.5)*2*np.pi/512
        alpha = np.deg2rad(self.alpha_raw_deg)
        normal_max = 0.
        for body in (0,1):
            radial = np.exp(1j*theta)
            z = self.centers[body]+self.flow.radii[body]*radial
            normal_max = max(normal_max,float(np.max(abs(np.real(
                self.flow.velocity(z,alpha,self.gamma_raw)*radial)))))
        half = TwoCircleFlow(self.centers,self.flow.radii,self.reflections//2)
        return {"circle_normal_velocity_max": normal_max,
                "stagnation_velocity_max": float(np.max(abs(self.flow.velocity(
                    np.asarray(self.trailing_points),alpha,self.gamma_raw)))),
                "circulation_half_series_difference_max": float(np.max(abs(
                    self.gamma_raw-half.circulations(self.trailing_points,alpha)))),
                "reflections": self.reflections,
                "half_series_reflections": self.reflections//2}

    def provenance(self):
        le,chord = self.normalization
        return {"source_url": SOURCE_URL,"source_sha256": SOURCE_SHA256,
                "circle_radii": [1.096,.5],"circle_center_offset": [1.9,-.4],
                "beta_rad": [.001,.08725],"k": [1.096,1.15],"tau_rad": [.17279,.17279],
                "mapping_direction_deg": self.deflection_deg,"alpha_raw_deg": self.alpha_raw_deg,
                "intermediate_radius_rule": "abs(K1(T2)-K1(C2)); no fitting",
                "reflections": self.reflections,"leading_reference_raw": [le.real,le.imag],
                "chord_vector_raw": [chord.real,chord.imag],"alpha_normalized_deg": self.alpha_deg,
                "published_airfoil_reproduction_certified": False}
