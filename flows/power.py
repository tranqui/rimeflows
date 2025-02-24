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

try: from .planar import PlanarFlowField
except: from planar import PlanarFlowField

class FlowField(PlanarFlowField):
    def __init__(self, m, k=1):
        self.m = m
        self.k = k

    def u(self, x, y):
        return -self.k * np.abs(x)**self.m

    def v(self, x, y):
        return self.m*self.k * np.abs(x)**(self.m-1) * y

    def critical_stokes(self, *args, **kwargs):
        if self.m == 1: return 0.25/self.k
        else: return PlanarFlowField.critical_stokes(self, *args, **kwargs)

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    m = [1, 1.5, 2]
    fig, axes = plt.subplots(nrows=len(m), sharex=True, sharey=True)

    for ax, m in zip(axes, m):
        flow = FlowField(m=m)

        L = 3
        x0 = 1.5*L
        for y0 in np.linspace(0, L**0.5, 11)**2:
            x, y = flow.streamline([x0, y0], max_step=0.1)
            ax.plot(x, y, 'k-', lw=0.5)
            ax.plot(x, -y, 'k-', lw=0.5)

        ax.fill_between([-0.5, 0], -L, L, fc='r', zorder=10)
        ax.set_ylabel('$y$')

    plt.xlabel('$x$')
    plt.xlim([-0.5, 1.5*L])
    plt.ylim([-L, L])

    plt.show()
