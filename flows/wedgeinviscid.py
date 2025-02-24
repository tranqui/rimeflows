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

try: from .wedgebase import WedgeFlowField
except: from wedgebase import WedgeFlowField

def wedge_exponent(alpha):
    """Determine the exponent of the limiting stagnation point flow at
    infinite Reynolds number into a wedge at a specified angle.

    Parameters:
        alpha (float): Number in [0, 1] parameterising angle of wedge ($\theta = \alpha \pi$).

    Returns:
        The value of the exponent (float).
    """
    return alpha / (2 - alpha)

class FlowField(WedgeFlowField):
    def __init__(self, alpha, B=1.):
        self.alpha = alpha
        self.m = wedge_exponent(alpha)
        self.B = B

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

        u_r = -self.B * r**self.m * np.cos((self.m+1) * (th - np.pi))
        u_th = -self.B * r**self.m * np.sin((self.m+1) * (th - np.pi))

        return u_r, u_th

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print('exponents for different wedges:')
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
