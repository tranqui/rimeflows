#!/usr/bin/env python3
# Copyright (C) 2022 Joshua Robinson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from scipy import special
from scipy.differentiate import derivative
from scipy.special import iv, kn
def knp(m, x): return special.kvp(m, x, 1)
def ivp(m, x): return special.ivp(m, x, 1)

try: from .cylinder import CylindricalFlowField
except: from cylinder import CylindricalFlowField

def phi(m, n, x):
    return (kn(m+1, x) + kn(m-1, x)) * (iv(m-n, x) + iv(m+n, x)) + \
        kn(m, x) * (iv(m-n-1, x) + iv(m-n+1, x) + iv(m+n-1, x) + iv(m+n+1, x))

def phip(m, n, x):
    return (knp(m+1, x) + knp(m-1, x)) * (iv(m-n, x) + iv(m+n, x)) + \
        (kn(m+1, x) + kn(m-1, x)) * (ivp(m-n, x) + ivp(m+n, x)) + \
        knp(m, x) * (iv(m-n-1, x) + iv(m-n+1, x) + iv(m+n-1, x) + iv(m+n+1, x)) + \
        kn(m, x) * (ivp(m-n-1, x) + ivp(m-n+1, x) + ivp(m+n-1, x) + ivp(m+n+1, x))

class FlowField(CylindricalFlowField):
    """Flow field obtained for low Reynolds number via the Renormalization Group (RG) approach
    to finding uniformly valid asymptotic approximations.

    Reference:
        Veysey, J. and Goldenfeld, N., Rev. Mod. Phys. 79, 883 (2007).
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
            dpdr = self.Re
            dPsi_dp = 1 + self.Re * ( self.A1 - self.B1/p**2 +
                                      self.X0*(phi(0, 1, 0.5*p) + 0.5*p*phip(0, 1, 0.5*p)) )
            return (dPsi_dp / self.Re**2) * dpdr

    def fpp(self, r, dx=1e-4):
        """Second (radial) derivative of radial part of streamfunction."""
        return derivative(self.fp, r, dx=dx)

    def fppp(self, r, dx=1e-4):
        """Third (radial) derivative of radial part of streamfunction."""
        return derivative(self.fp, r, dx=dx, n=2)

    @property
    def drag_coefficient(self):
        """Drag coefficient for flow field."""
        return -np.pi * self.fppp(1)

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

    Re = np.linspace(0, 1, 1001)
    C = np.empty(Re.size)
    for i, R in enumerate(Re):
        flow = FlowField(R)
        C[i] = R * flow.drag_coefficient / (4*np.pi)

    plt.figure()
    plt.plot(Re, C)
    plt.xlabel(r'$\mathrm{Re}$')
    plt.ylabel(r'$\mathrm{Re} \, C_D / 4\pi$')

    plt.show()
