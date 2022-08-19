#!/usr/bin/env python3

import numpy as np
from scipy import special
from scipy.special import iv, kn
def knp(m, x): return special.kvp(m, x, 1)
def ivp(m, x): return special.ivp(m, x, 1)

try: from .cylinder import CylindricalFlowField
except: from cylinder import CylindricalFlowField

def phi(m, n, x):
    return (kn(m+1, x) + kn(m-1, x)) * (iv(m-n, x) + iv(m+n, x)) + \
        kn(m, x) * (iv(m-n-1, x) - iv(m-n+1, x) - iv(m+n-1, x) + iv(m+n+1, x))

def phip(m, n, x):
    return (knp(m+1, x) + knp(m-1, x)) * (iv(m-n, x) + iv(m+n, x)) + \
        (kn(m+1, x) + kn(m-1, x)) * (ivp(m-n, x) + ivp(m+n, x)) + \
        knp(m, x) * (iv(m-n-1, x) - iv(m-n+1, x) - iv(m+n-1, x) + iv(m+n+1, x)) + \
        kn(m, x) * (ivp(m-n-1, x) - ivp(m-n+1, x) - ivp(m+n-1, x) + ivp(m+n+1, x))

class FlowField(CylindricalFlowField):
    """Flow field obtained for low Reynolds number via the Renormalization Group (RG) approach
    to finding uniformly valid asymptotic approximations.

    Reference:
        Veysey, J. and Goldenfeld, N, Rev. Mod. Phys. 79, 883 (2007).
    """

    def __init__(self, Re):
        self.Re = Re
        self.A1 = 0
        self.B1 = -Re**2 * phip(0, 1, 0.5*Re) / (4*phi(0, 1, 0.5*Re) + Re*phip(0, 1, 0.5*Re))
        self.X0 = -4/Re / (4*phi(0, 1, 0.5*Re) + Re*phip(0, 1, 0.5*Re))

    def f(self, r):
        """Radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            p = self.Re * r
            Psi = p + self.Re * (self.A1*p + self.B1/p + self.X0*p*phi(0, 1, 0.5*p))
            return Psi / self.Re**2

    def fp(self, r):
        """First (radial) derivative of radial part of streamfunction."""
        with np.errstate(divide='ignore'):
            p = self.Re * r
            dPsi = self.Re * ( 1 + self.Re * (self.A1 - self.B1/p**2 +
                                              self.X0*(phi(0, 1, 0.5*p) + p*phip(0, 1, 0.5*p))) )
            return dPsi / self.Re**2

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    Re = 2
    flow = FlowField(Re)

    L = 3
    x0 = -L
    for y0 in np.linspace(-L, L, 31):
        if abs(y0) < 1e-2:
            plt.axhline(y=y0, lw=0.5)
            continue

        print(y0)
        x, y = flow.streamline_cartesian([x0, y0])
        plt.plot(x, y, 'k-', lw=0.5)

    from matplotlib.patches import Circle
    plt.gca().add_patch(Circle((0,0), 1, fc='r', ec='None', zorder=10))
    plt.gca().add_patch(Circle((0,0), 1 + 1/np.sqrt(Re), fc='None', ec='k', ls='--', zorder=10))

    plt.xlim([-L, L])
    plt.ylim([-L, L])

    plt.show()
