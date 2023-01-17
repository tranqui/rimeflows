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

try: from .refine import logical_refine
except: from refine import logical_refine

class OnAxisFlow:
    def __init__(self, flow):
        self.full_flow = flow

    def u(self, x):
        return self.full_flow.u(x, 0)

    def trajectory(self, r0, St, step=1e-12, atol=1e-12, rtol=1e-12, tmax=1e2, max_step=np.inf, xmin=-5, xmax=5, vmin=-5, vmax=5):
        terminate = lambda t,r: r[0] < xmax and r[0] > xmin and r[1] < vmax and r[1] > vmin
        terminate.terminal = True
        terminate.direction = -1

        xp = lambda x,v: v
        vp = lambda x,v: (self.u(x) - v)/St
        rp = lambda t,r: [xp(*r), vp(*r)]
        return integrate.solve_ivp(rp, (0, tmax), r0, first_step=step, atol=atol, rtol=rtol, max_step=max_step, dense_output=True, events=terminate)

    def streamline(self, x0, v0, St=1, tmax=1e2, N=2000, **kwargs):
        r0 = np.array([x0, v0])
        sol1 = self.trajectory(r0, St, tmax=-tmax, **kwargs)
        sol2 = self.trajectory(r0, St, tmax=tmax, **kwargs)

        t1 = np.linspace(sol1.t[-1], 0, N//2)
        t2 = np.linspace(0, sol2.t[-1], N-N//2)
        x1, v1 = sol1.sol(t1)
        x2, v2 = sol2.sol(t2)
        x = np.concatenate([x1, x2])
        v = np.concatenate([v1, v2])

        return x, v

    def separatrix(self, *args, St=1, step=1e-12, tmax=1e2, N=2000, return_critical_point=False, **kwargs):
        r0 = step * np.array([1, -1])
        sol1 = self.trajectory(r0, St, *args, tmax=-tmax, **kwargs)
        sol2 = self.trajectory(-r0, St, *args, tmax=-tmax, **kwargs)

        # Find location of critical point where separatrix crosses v' = 0:
        umin = lambda t: sol1.sol(t)[1]
        from scipy.optimize import minimize_scalar
        imin = np.argmin(sol1.y[1])
        tguess1 = sol1.t[imin-1]
        tguess2 = sol1.t[imin+1]

        t1 = np.linspace(sol1.t[-1], 0, N//2)
        t2 = np.linspace(0, sol2.t[-1], N-N//2)
        x1, v1 = sol1.sol(t1)
        x2, v2 = sol2.sol(t2)
        x = np.concatenate([x1, x2])
        v = np.concatenate([v1, v2])

        if return_critical_point:
            critical_point = sol1.sol(minimize_scalar(umin, (tguess1, tguess2)).x)
            return x, v, critical_point

        return x, v

    def nullcline(self, x):
        return self.u(x)

    def limit_colliding_streamline(self, St=1, N=2000, step=1e-12, **kwargs):
        r0 = np.array([-step, 0])
        sol = self.trajectory(r0, St, step=step, **kwargs)

        t = np.linspace(0, sol.t[-1], N-1)
        x, v = sol.sol(t)
        x = np.concatenate([[0], x])
        v = np.concatenate([[0], v])
        return x, v

class PlanarFlowField:
    @property
    def on_axis(self):
        return OnAxisFlow(self)

    @property
    def default_starting_distance(self):
        """Default initial condition for stagnation point flows."""
        return 1

    def u(self, x, y):
        """x-component of flow field."""
        raise NotImplementedError

    def v(self, x, y):
        """y-component of flow field."""
        raise NotImplementedError

    def eom(self, x, y, ux, uy, St):
        """Equation of motion for massive particle immersed in flow field."""
        return ux, uy, uy**2 - (ux - self.u(x,y))/St, -2*ux*uy - (uy - self.v(x,y))/St

    def streamline_eom(self, x, y):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(x,y), self.v(x,y)

    def trajectory(self, r0, St, xmax=None, ymax=5, return_collision=False,
                   tmax=1e2, max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field (in polar coordinates).

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

        collide = lambda t, x: -x[0]
        collide.terminal = True

        if xmax is None:
            eps = 1e-6
            xmax = r0[0] + eps

        out_of_bounds1 = lambda t, x: x[0] - xmax
        out_of_bounds2 = lambda t, x: np.abs(x[1]) - ymax
        out_of_bounds1.terminal = True
        out_of_bounds2.terminal = True

        events = (collide, out_of_bounds1, out_of_bounds2)

        with np.errstate(divide='ignore', invalid='ignore'):
            trajectory = integrate.solve_ivp(eom, [0, tmax], r0, events=events,
                                             vectorized=True, dense_output=True,
                                             max_step=max_step, rtol=rtol, atol=atol,
                                             **kwargs)

        if return_collision: return trajectory.t, trajectory.y, trajectory.t_events[0].size > 0
        else: return trajectory.t, trajectory.y

    def streamline(self, r0, tmax=1e2, **kwargs):
        """Obtain the streamline for the flow field.

        Args:
            r0: reference point [2d r0 = (x, y)] along streamline. Streamline will be obtained
                  by following flow in both directions from this point.
            tmax: maximum time to follow streamline for.
        Returns:
            x, y: coordinates along streamline.
        """

        t1, (x1, y1) = self.trajectory(r0, St=0, tmax=-tmax, **kwargs)
        t2, (x2, y2) = self.trajectory(r0, St=0, tmax=tmax, **kwargs)

        x = np.concatenate([np.flipud(x1), x2])
        y = np.concatenate([np.flipud(y1), y2])

        return x, y

    def does_collide(self, *args, **kwargs):
        """Determine whether particle collides with plane under specified conditions."""
        _, _, collides = self.trajectory(*args, return_collision=True, **kwargs)
        return collides

    def capture_efficiency(self, x, St, niters=25, quiet=True, *args, return_bounds=False, message_header='', **kwargs):
        """Estimate efficiency of point particle capture.

        Args:
            x: initial condition for x.
            St: value of Stokes number we are evaluating efficiency at.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
            message_header: preface to iteration updates if not quiet.
        Returns:
            Efficiency of particle capture efficiency.
        """

        collides = lambda y: self.does_collide([x, y], St, *args, **kwargs)

        if not collides(0):
            if not quiet: print('on-axis does not lead to collision, so efficiency is zero!')
            return 0

        ylow, yhigh = logical_refine(collides, niters=niters, quiet=quiet, message_header=message_header)
        if return_bounds: return ylow, yhigh
        else:
            yestimate = 0.5 * (ylow + yhigh)
            return 2*yestimate # capture efficiency covers both +ve and -ve y, so x2

    def on_axis_eom(self, x, u, St):
        """Equation of motion for massive particle immersed in flow field approaching on-axis."""
        return u, -(u - self.u(x, 0))/St

    def on_axis_streamline_eom(self, x):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(x, 0)

    def on_axis_trajectory(self, r0, St, tmax=1e2, return_collision=False,
                           max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field that's approaching the
            stagnation point on-axis (y = 0).

        Args:
            r0: initial condition (n-component vector, n=1 or 2).
                  When n=1 only r0 = x are specified, and the initial particle velocity will
                    be set to be aligned with the flow field.
                  For St > 0, n=2 one can optionally set also the initial velocity, i.e.
                    r0 = (x, u).
            St: Stokes number
            tmax: max time before aborting integration.
            return_collision: If True, will return additional boolean stating whether a collision
                  event occurred (i.e. the particle impacted on the plane to terminate).
            max_step: maximum step during integration (sent to scipy.integrate.solve_ivp).
            rtol: relative error tolerance during integration (sent to scipy.integrate.solve_ivp).
            atol: absolute error tolerance during integration (sent to scipy.integrate.solve_ivp).
            kwargs: additional arguments to send to scipy.integrate.solve_ivp.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            x: coordinates at each frame (with dimensionality 1 for St = 0 or 2 for St > 0).
            collides (optional): True/False depending on whether particle terminates on collision
                                   with plane. Only returned if return_collision set to True.
        """

        try: len(r0)
        except: r0 = [r0]

        if St == 0:
            assert len(r0) == 1
            eom = lambda t, x: self.on_axis_streamline_eom(*x)
        else:
            if len(r0) == 1: r0 = np.concatenate([r0, [self.u(*r0, 0)]])
            else: assert len(r0) == 2
            eom = lambda t, x: self.on_axis_eom(*x, St=St)

        collide = lambda t, x: -x[0]
        collide.terminal = True
        events = (collide,)

        with np.errstate(divide='ignore', invalid='ignore'):
            trajectory = integrate.solve_ivp(eom, [0, tmax], r0, events=events,
                                             vectorized=True, dense_output=True,
                                             max_step=max_step, rtol=rtol, atol=atol,
                                             **kwargs)

        if len(r0) == 1: trajectory.y = trajectory.y.reshape(-1)

        if return_collision: return trajectory.t, trajectory.y, trajectory.t_events[0].size > 0
        else: return trajectory.t, trajectory.y

    def on_axis_does_collide(self, *args, **kwargs):
        """Determine whether particle collides with plane under specified conditions."""
        _, _, collides = self.on_axis_trajectory(*args, return_collision=True, **kwargs)
        return collides

    def critical_stokes(self, x, niters=50, quiet=True, *args, message_header='', **kwargs):
        """Estimate critical Stokes number above which point particle capture occurs.

        This requires solving the on-axis problem only.

        Args:
            x: initial condition for x.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
            message_header: preface to iteration updates if not quiet.
        """

        collides = lambda St: not self.on_axis_does_collide(x, St, *args, **kwargs)
        Stc_low, Stc_high = logical_refine(collides, niters=niters, quiet=quiet, message_header=message_header)
        return 0.5 * (Stc_low + Stc_high)
