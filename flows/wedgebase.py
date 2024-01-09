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
from scipy import integrate

try: from .planar import PlanarFlowField
except: from planar import PlanarFlowField

class WedgeFlowField(PlanarFlowField):
    def u_polar(self, x, y):
        """Flow field in polar coordinates.

        Args:
            x: x-position (Cartesian).
            y: y-position (Cartesian).
        Returns:
            r and theta polar components of flow field.
        """

        raise NotImplementedError

    def u(self, x, y):
        """x-component of flow field."""
        th = np.arctan2(y, x)
        ur, ut = self.u_polar(x, y)
        return ur*np.cos(th) - ut*np.sin(th)

    def v(self, x, y):
        """y-component of flow field."""
        th = np.arctan2(y, x)
        ur, ut = self.u_polar(x, y)
        return ur*np.sin(th) + ut*np.cos(th)

    def y_wedge(self, x):
        """y-position of wedge boundary at a given value of x >= 0."""
        return np.tan(0.5*self.alpha*np.pi) * x

    def collision(self, x, y):
        """Function that changes sign to indicate when a collision occurs."""
        return self.y_wedge(x) - y

    def trajectory(self, r0, St, xmax=None, ymax=5, return_collision=False,
                   tmax=1e2, max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field.

        Args:
            r0: initial condition (n-component vector, n=2 or 4).
                  When n=2 only r0 = (x, y) are specified, and the initial particle velocity will
                    be set to be aligned with the flow field.
                  For St > 0, n=4 one can optionally set also the initial velocities, i.e.
                    r0 = (x, y, u, v).
            St: Stokes number
            xmax: max value of x before particle is considered to have escaped the plane.
                    If None, will default to initial x coordinate.
            ymax: max value of abs(y) before particle is considered to have escaped the plane.
            return_collision: If True, will return additional boolean stating whether a collision
                    event occurred (i.e. the particle impacted on the cylinder to terminate).
            tmax: max time before aborting integration.
            max_step: maximum step during integration (sent to scipy.integrate.solve_ivp).
            rtol: relative error tolerance during integration (sent to scipy.integrate.solve_ivp).
            atol: absolute error tolerance during integration (sent to scipy.integrate.solve_ivp).
            kwargs: additional arguments to send to scipy.integrate.solve_ivp.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            x: coordinates at each frame (with dimensionality 2 for St = 0 or 4 for St > 0).
        """

        if St == 0:
            assert len(r0) == 2
            eom = lambda t, x: self.streamline_eom(*x)
        else:
            if len(r0) == 2: r0 = np.concatenate([r0, [self.u(*r0), self.v(*r0)]])
            else: assert len(r0) == 4
            eom = lambda t, x: self.eom(*x, St=St)

        collide = lambda t, x: self.collision(*x[:2])
        collide.terminal = True

        out_of_bounds = lambda t, x: np.abs(x[1]) - ymax
        out_of_bounds.terminal = True

        events = (collide, out_of_bounds)

        with np.errstate(divide='ignore', invalid='ignore'):
            trajectory = integrate.solve_ivp(eom, [0, tmax], r0, events=events,
                                             vectorized=True, dense_output=True,
                                             max_step=max_step, rtol=rtol, atol=atol,
                                             **kwargs)

        if return_collision: return trajectory.t, trajectory.y, trajectory.t_events[0].size > 0
        else: return trajectory.t, trajectory.y
