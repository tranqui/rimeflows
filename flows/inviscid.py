#!/usr/bin/env python3

import numpy as np

try: from .cylinder import CylindricalFlowField
except: from cylinder import CylindricalFlowField

class FlowField(CylindricalFlowField):
    """Potential flow around a cylinder."""

    def f(self, r):
        """Radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return r - 1/r

    def fp(self, r):
        """First (radial) derivative of radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return 1 + 1/r**2

    def fpp(self, r):
        """Second (radial) derivative of radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            return -2/r**3

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    flow = FlowField()

    L = 3
    x0 = -L
    for y0 in np.linspace(-L, L, 25):
        x, y = flow.streamline_cartesian([x0, y0])
        plt.plot(x, y, 'k-', lw=0.5)

    from matplotlib.patches import Circle
    plt.gca().add_patch(Circle((0,0), 1, fc='r', ec='None', zorder=10))

    plt.xlim([-L, L])
    plt.ylim([-L, L])

    plt.show()
