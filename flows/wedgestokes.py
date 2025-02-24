#!/usr/bin/env python3
# Copyright (C) 2025 Joshua Robinson
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
from scipy.optimize import brentq, minimize_scalar

try: from .wedgebase import WedgeFlowField
except: from wedgebase import WedgeFlowField

def wedge_exponent(alpha, N=10000):
    """Determine the exponent of the limiting stagnation point flow at
    zero Reynolds number into a wedge at a specified angle.

    This requires a numerical procedure to solve the underlying equation.

    Parameters:
        alpha (float): Number in [0, 1] parameterising angle of wedge ($\theta = \alpha \pi$).
        N (int): Number of points to discretise the interval to initially guess the root.

    Returns:
        The value of the exponent (float).
    """

    # Solution is the zero of this residual function:
    def residual(m):
        term1 = (m+1) * np.cos((m+1) * (np.pi - alpha*np.pi/2.0)) * \
                        np.sin((m-1) * (np.pi - alpha*np.pi/2.0))
        term2 = (m-1) * np.sin((m+1) * (np.pi - alpha*np.pi/2.0)) * \
                        np.cos((m-1) * (np.pi - alpha*np.pi/2.0))
        return term1 - term2

    m_array = np.linspace(0, 10, N)[1:]
    f_array = residual(m_array)

    # Check for sign changes
    sign_changes = np.where(np.diff(np.sign(f_array)))[0]

    # Find root in each interval with a sign change
    roots = []
    for index in sign_changes:
        root = brentq(residual, m_array[index], m_array[index + 1])
        roots.append(root)

    eps = 1e-10
    try:
        assert len(roots) <= 2
        roots = np.array(roots)

        if abs(roots[0] - 1.0) < eps:
            root = roots[1]
        if abs(roots[1] - 1.0) < eps:
            root = roots[0]

        return root
    except:
        # As m -> 1 the two roots converge on each other, so we may need
        # to look for a single root around m=1. Also, at small and large
        # angles we do not seem to find two roots so need to modify the
        # search there.
        bounds = [0, 2]
        if alpha < 0.4: bounds = [0, 0.75]
        elif alpha > 0.8:  bounds = [1.25, 2]
        return minimize_scalar(lambda m: residual(m)**2, bounds=bounds).x

def critical_alpha(guess=0.5):
    """Critical wedge angle where the limiting stagnation point flow at zero Reynolds number
    reduces to simple harmonic motion."""
    return brentq(lambda alpha: wedge_exponent(alpha) - 1, 0.2, 0.8)

class FlowField(WedgeFlowField):
    def __init__(self, alpha, B=1.):
        self.alpha = alpha
        self.m = wedge_exponent(alpha)
        self.B = B

        angle1 = (self.m+1) * (np.pi - 0.5*self.alpha*np.pi)
        angle2 = (self.m-1) * (np.pi - 0.5*self.alpha*np.pi)
        self.D = -np.sin(angle1)/np.sin(angle2) * self.B

    def stream_function(self, x, y):
        """Stream-function for flow field.

        Args:
            x (float): x-position (Cartesian).
            y (float): y-position (Cartesian).
        Returns:
            Stream function. The derivatives of this stream function yield the flow field.
        """
        # 1st convert from Cartesians to circular polars
        r = np.sqrt(x**2 + y**2)
        th = np.arctan2(y, x)
        if th < 0: th += 2*np.pi

        angle1 = (self.m+1) * (np.pi-th)
        angle2 = (self.m-1) * (np.pi-th)
        f_th = (self.B/(self.m+1)) * np.sin(angle1) + (self.D/(self.m+1)) * np.sin(angle2)
        psi = r**(self.m+1) * f_th
        if th < 0.5*self.alpha*np.pi: psi = 0.

        return psi

    def u_polar(self, x, y):
        """Flow field in polar coordinates.

        Args:
            x: x-position (Cartesian).
            y: y-position (Cartesian).
        Returns:
            r and theta polar components of flow field.
        """
        # 1st convert from Cartesians to circular polars
        r = np.sqrt(x**2 + y**2)
        th = np.arctan2(y, x)
        if th < 0: th += 2*np.pi

        angle1 = (self.m+1) * (np.pi-th)
        angle2 = (self.m-1) * (np.pi-th)
        square_brac_r = self.B * np.cos(angle1) + self.D * ((self.m-1)/(self.m+1)) * np.cos(angle2)
        u_r = -(r**self.m) * square_brac_r
        square_brac_theta = self.B * np.sin(angle1) + self.D * np.sin(angle2)
        u_th = -(r**self.m) * square_brac_theta

        return u_r, u_th

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print('critical (i.e. m=1) alpha is alpha={:.4f}'.format(critical_alpha()))

    for alpha in np.linspace(0, 1, 5):
        print('    alpha={:.4f}: m={:.4f}'.format(alpha, wedge_exponent(alpha)))

    f_m = np.vectorize(wedge_exponent)

    # Plot in finer detail:
    # alpha = np.linspace(0.5, 0.6, 1000)
    # m = f_m(alpha)
    # plt.plot(alpha, m-1)
    # plt.axhline(y=0)
    # plt.xlabel(r'$\alpha$')
    # plt.ylabel(r'$m$')
    # plt.show()

    flow = FlowField(alpha=0.4)

    # Draw streamlines.
    L = 0.25
    x0 = -L
    for y0 in np.linspace(0, L**0.5, 7)**2:
        x, y = flow.streamline([x0, y0], max_step=0.1)
        plt.plot(x, y, 'k-', lw=0.5)
        plt.plot(x, -y, 'k-', lw=0.5)

    # Draw inertial trajectories that collide with the wedge.
    for y0 in np.linspace(0, (0.25*L)**0.5, 5)[1:]**2:
        St = 1
        t, (x, y, _, _), collides = flow.trajectory([x0, y0], St, max_step=0.1, return_collision=True)
        assert collides
        pl, = plt.plot(x, y, 'b-', lw=0.5)
        plt.plot(x[-1], y[-1], 'o', c=pl.get_color())

    # Draw obstacle.
    x = np.linspace(0, 10*L, 100)
    plt.fill_between(x, -flow.y_wedge(x), flow.y_wedge(x), color='red', ec='None')

    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim([-L, L])
    plt.ylim([-L, L])

    plt.show()
