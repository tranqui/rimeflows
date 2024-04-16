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
from scipy.integrate import solve_ivp

try: from .planar import PlanarFlowField
except: from planar import PlanarFlowField

class FlowField(PlanarFlowField):
    def __init__(self, epsilon, k=1):
        """
        Args:
            epsilon: number in the range [0, 1]. epsilon=0 corresponds to Stokes case, whereas
                   epsilon=1 is SHM.
        """
        self.epsilon = epsilon
        self.k = k

    def u(self, x, y):
        return -self.k * np.abs(x)**(2-self.epsilon)

    def v(self, x, y):
        return 2*self.k * x * y

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    flow = FlowField()

    L = 3
    x0 = 1.5*L
    for y0 in np.linspace(0, L**0.5, 11)**2:
        x, y = flow.streamline([x0, y0], max_step=0.1)
        plt.plot(x, y, 'k-', lw=0.5)
        plt.plot(x, -y, 'k-', lw=0.5)

    plt.fill_between([-0.5, 0], -L, L, fc='r', zorder=10)

    plt.xlim([-0.5, 1.5*L])
    plt.ylim([-L, L])

    plt.show()
