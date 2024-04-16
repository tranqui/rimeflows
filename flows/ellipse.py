import numpy as np

try: from .planar import PlanarFlowFieldInertial
except: from planar import PlanarFlowFieldInertial

class FlowField(PlanarFlowFieldInertial):
    """Implementation of inviscid flow around an ellipse.

    It may seem odd we do not base the implementation on our potential flow
    around a cylinder which uses polar coordinates. Because the ellipse may
    be rotated, it is more natural to base our implementation on the flow
    around a chord for which we used Cartesian coordinates and so we use Cartesian
    here.
    """

    def __init__(self, a, b, alpha=0, U=1):
        """
        Args:
            a: first semiaxis
            b: second semiaxis
            alpha: angle of rotation of ellipse relative to its standard form where
                   it's aligned with the axes with equation:
                       $$\frac{x^2}{a^2} + \frac{y^2}{b^2} = 1$$.
            U: flow velocity, which can be normalised to 1.
        """

        self.a, self.b = a, b
        self.alpha = alpha
        self.U = U

    @property
    def semi_major_axis(self):
        return max(self.a, self.b)

    @property
    def semi_minor_axis(self):
        return min(self.a, self.b)

    @property
    def eccentricity(self):
        return np.sqrt(1 - (self.semi_minor_axis/self.semi_major_axis)**2)

    @property
    def a_chord(self):
        return 0.5 * (self.a + self.b)

    @property
    def c_conformal(self):
        return np.sqrt(0.5*self.a_chord * (self.semi_major_axis - self.semi_minor_axis))

    @property
    def matplotlib_draw_angle(self):
        return -180*alpha/np.pi

    def complex_potential(self, x, y):
        """Representation of flow field via its complex potential.

        Returns:
            Complex number u - i v.
        """

        z = x + y*1j

        exp_fac = np.exp(-1j*(0.5*np.pi - self.alpha))
        exp_fac_sq = exp_fac**2
        if np.real(z*exp_fac) > 0:
            bigZ = 0.5*(z*exp_fac + np.sqrt(z**2*exp_fac_sq - 4*self.c_conformal**2))
        else:
            bigZ = 0.5*(z*exp_fac - np.sqrt(z**2*exp_fac_sq - 4*self.c_conformal**2))

        numerator = self.U*(1/exp_fac - exp_fac*self.a_chord**2/bigZ**2)
        denominator = (1 - self.c_conformal**2/bigZ**2) / exp_fac

        return numerator / denominator

    def u(self, x, y):
        return np.real(self.complex_potential(x, y))

    def v(self, x, y):
        return -np.imag(self.complex_potential(x, y))

    def collision(self, x, y):
        """Indicator function for detecting collisions with the ellipse.

        Args:
            x: x-coordinate.
            y: y-coordinate.
    
        Returns:
            The value -1 or 1. It changes sign when crossing to the interior of the ellipse.
        """

        th = np.arctan2(y, x) + self.alpha # angle relative to a-axis of ellipse
        r2_border = (self.a * self.b)**2 / (self.a**2*np.sin(th)**2 + self.b**2*np.cos(th)**2)
        r2 = x**2 + y**2
        return r2 > r2_border

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    a, b = 0.1, 1
    alpha = np.pi/4
    flow = FlowField(a, b, alpha)
    ax = plt.gca()

    L = 1.5
    eps = 1e-6
    x0 = -L
    for y0 in np.linspace(-L, L, 11):
        x, y = flow.streamline([x0, y0], ymax=(L+eps))
        plt.plot(x, y, 'k-', lw=0.5)

    y0 = 0.25
    for St in [0.1, 0.9, 1, 10]:
        t, (x, y, _, _), collides = flow.trajectory([x0, y0], St=St, max_step=1e-2, ymax=(L+eps), return_collision=True)
        pl, = plt.plot(x, y, 'b-', lw=1)
        if collides: plt.plot(x[-1], y[-1], 'o', c=pl.get_color())

    ellipse = Ellipse((0,0), 2*a, 2*b, ec='None', fc='red', angle=flow.matplotlib_draw_angle)
    ax.add_patch(ellipse)

    ax.set_xlim([-L,L])
    ax.set_ylim([-L,L])
    plt.show()