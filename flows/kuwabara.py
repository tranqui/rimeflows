#!/usr/bin/env python3

import numpy as np

try: from .cylinder import CylindricalFlowField
except: from cylinder import CylindricalFlowField

class FlowField(CylindricalFlowField):
    """Kuwabara flow field around an idealised cylindrical grain in porous media."""

    def __init__(self, alpha):
        self.alpha = alpha
        self.hydrodynamic_factor = -0.5*np.log(alpha) - 0.75 + alpha - 0.25*alpha**2
        self.outer_boundary = 1/np.sqrt(alpha)

        self.A = 0.5 * (1 - 0.5*alpha)
        self.B = 0.5 * (alpha - 1)
        self.C = -0.25 * alpha
        self.D = 1

    @property
    def l(self):
        """Short-hand for the size of the outer boundary (in units of fibre radius)."""
        return self.outer_boundary

    @property
    def default_starting_distance(self):
        """Default initial condition for stagnation point flows (in units of fibre radius)."""
        return self.outer_boundary

    def f(self, r):
        """Radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return (self.A/r + self.B*r + self.C*r**3 + self.D*r*np.log(r)) / self.hydrodynamic_factor

    def fp(self, r):
        """First (radial) derivative of radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return (-self.A/r**2 + self.B + 3*self.C*r**2 + self.D*(1 + np.log(r))) / self.hydrodynamic_factor

    def fpp(self, r):
        """Second (radial) derivative of radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return (2*self.A/r**3 + 6*self.C*r + self.D/r) / self.hydrodynamic_factor

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    alpha = 0.125
    flow = FlowField(alpha)

    L = 3
    x0 = -L
    for y0 in np.linspace(-L, L, 25):
        x, y = flow.streamline_cartesian([x0, y0])
        plt.plot(x, y, 'k-', lw=0.5)

    from matplotlib.patches import Circle
    plt.gca().add_patch(Circle((0,0), 1, fc='r', ec='None', zorder=10))
    plt.gca().add_patch(Circle((0,0), flow.l, fc='None', ec='k', ls='--', zorder=10))

    plt.xlim([-L, L])
    plt.ylim([-L, L])

    plt.show()
