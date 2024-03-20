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

try: from .planar import PlanarFlowFieldInertial
except: from planar import PlanarFlowFieldInertial

class FlowField(PlanarFlowFieldInertial):
    def __init__(self, chord_length, alpha=0.5*np.pi, U=1):
        """
        Args:
            chord_length
            alpha: angle of chord to the incoming stream, e.g. alpha=0 for parallel
                   or alpha=pi/2 for perpendicular.
            U: flow velocity, which can be normalised to 1.
        """

        self.chord_length = chord_length
        self.alpha = alpha
        self.U = U

    @property
    def a(self):
        return self.chord_length / 4

    def complex_potential(self, x, y):
        """Representation of flow field via its complex potential.

        Returns:
            Complex number u - i v.
        """

        z = y + x*1j
        if y < 0: bigZ = 0.5 * (z - np.sqrt(z**2 - 4*self.a**2))
        else: bigZ = 0.5 * (z + np.sqrt(z**2 - 4*self.a**2))

        numerator = self.U*(np.exp(-self.alpha*1j) - np.exp(self.alpha*1j)*self.a**2/bigZ**2)
        denominator = (1.0 - self.a**2/bigZ**2)

        return numerator / denominator

    def u(self, x, y):
        return -np.imag(self.complex_potential(x, y))

    def v(self, x, y):
        return np.real(self.complex_potential(x, y))

    def collision(self, x, y):
        """Indicator function for detecting collisions with the chord.

        We have to detect crossings of the y=0 line within the range x \in [-0.5, 0.5].

        Args:
            x: x-coordinate.
            y: y-coordinate.
    
        Returns:
            The value -1 or 1. It changes sign when crossing y=0 within the specified x-range.
        """

        if -0.5*self.chord_length <= y <= 0.5*self.chord_length:
            return np.sign(x)
        else:
            return np.nan

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    chord_length = 1
    alpha = np.pi/2
    flow = FlowField(chord_length, alpha)

    L = 1.5
    eps = 1e-6
    x0 = -L
    for y0 in np.linspace(-L, L, 21):
        x, y = flow.streamline([x0, y0], max_step=0.1, ymax=(L+eps))
        plt.plot(x, y, 'k-', lw=0.5)

    plt.plot([0, 0], [-0.5*chord_length, 0.5*chord_length], 'r-', lw=1)

    Stc = 0.12583899
    x0, y0, St = -1, -0.078, Stc + 1e-3
    t, (x, y, _, _), collides = flow.trajectory([x0, y0], St=St, max_step=1e-2, tmax=1e2, ymax=(L+eps), return_collision=True)
    pl, = plt.plot(x, y, 'b-', lw=1)
    if collides: plt.plot(x[-1], y[-1], 'o', c=pl.get_color())

    x0, y0, St = -1, 0.077, Stc + 1e-3
    t, (x, y, _, _), collides = flow.trajectory([x0, y0], St=St, max_step=1e-2, tmax=1e2, ymax=(L+eps), return_collision=True)
    pl, = plt.plot(x, y, 'b-', lw=1)
    if collides: plt.plot(x[-1], y[-1], 'o', c=pl.get_color())

    plt.ylabel('$y$')

    plt.xlabel('$x$')
    plt.xlim([-L, L])
    plt.ylim([-L, L])

    plt.show()
