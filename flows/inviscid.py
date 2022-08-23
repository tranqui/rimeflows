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
