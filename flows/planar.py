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

class PlanarOnAxisFlow:
    """Planar flow under constraint of being confined to the y-axis."""

    def __init__(self, flow):
        self.full_flow = flow

    def u(self, x):
        """Only the x-component of the velocity field is retained on-axis."""
        return self.full_flow.u(x, 0)

    def eom(self, x, u, St):
        """Equation of motion for massive particle immersed in flow field."""
        partial_eom = lambda full: [full[0], full[2]]
        return partial_eom(self.full_flow.eom(x, 0, u, 0, St))

    def streamline_eom(self, x):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(x)

    def trajectory(self, r0, St, xmin=None, xmax=None, umin=None, umax=None,
                   test_collision=True, return_collision=False,
                   tmax=1e2, max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field.

        Args:
            r0: initial condition (n-component vector, n=1 or 2).
                  When n=1 only x0 = x is specified, and the initial particle velocity will
                    be set to be aligned with the flow field.
                  For St > 0, n=2 one can optionally set also the initial velocity, i.e.
                    r0 = (x, u).
            St: Stokes number
            xmin: min value of x before particle is considered to have escaped the plane.
            xmax: max value of x before particle is considered to have escaped the plane.
            umin: min value of u before particle is considered to have escaped the plane.
            umax: max value of u before particle is considered to have escaped the plane.
            test_collision: if True, will test for and terminate on collision with object.
            return_collision: If True, will return additional boolean stating whether a collision
                    event occurred (i.e. the particle impacted on the cylinder to terminate).
            tmax: max time before aborting integration.
            max_step: maximum step during integration (sent to scipy.integrate.solve_ivp).
            rtol: relative error tolerance during integration (sent to scipy.integrate.solve_ivp).
            atol: absolute error tolerance during integration (sent to scipy.integrate.solve_ivp).
            kwargs: additional arguments to send to scipy.integrate.solve_ivp.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            x: coordinates at each frame (with dimensionality 1 for St = 0 or 2 for St > 0).
        """

        try: len(r0)
        except: r0 = [r0]

        if St == 0:
            assert len(r0) == 1
            eom = lambda t, x: self.streamline_eom(*x)
        else:
            if len(r0) == 1: r0 = np.concatenate([r0, [self.u(*r0)]])
            else: assert len(r0) == 2
            eom = lambda t, x: self.eom(*x, St=St)

        events = []
        if test_collision:
            collide = lambda t, x: self.full_flow.collision(x[0], 0)
            collide.terminal = True
            events += [collide]

        if xmin is not None:
            out_of_bounds_x1 = lambda t, x: xmin - x[0]
            out_of_bounds_x1.terminal = True
            events += [out_of_bounds_x1]

        if xmax is not None:
            out_of_bounds_x2 = lambda t, x: x[0] - xmax
            out_of_bounds_x2.terminal = True
            events += [out_of_bounds_x2]

        if umin is not None:
            out_of_bounds_u1 = lambda t, x: umin - x[1]
            out_of_bounds_u1.terminal = True
            events += [out_of_bounds_u1]

        if umax is not None:
            out_of_bounds_u2 = lambda t, x: x[1] - umax
            out_of_bounds_u2.terminal = True
            events += [out_of_bounds_u2]

        with np.errstate(divide='ignore', invalid='ignore'):
            trajectory = integrate.solve_ivp(eom, [0, tmax], r0, events=events,
                                             vectorized=True, dense_output=True,
                                             max_step=max_step, rtol=rtol, atol=atol,
                                             **kwargs)

        if return_collision: return trajectory.t, trajectory.y, trajectory.t_events[0].size > 0
        else: return trajectory.t, trajectory.y

    def does_collide(self, *args, **kwargs):
        """Determine whether particle collides with obstacle under specified conditions."""
        _, _, collides = self.trajectory(*args, return_collision=True, **kwargs)
        return collides

    def separatrix(self, *args, St=1, step=1e-12, tmax=1e2, N=2000,
                   return_critical_point=False, **kwargs):
        r0 = step * np.array([1, -1])
        t1, (x1, u1) = self.trajectory( r0, St, *args, tmax=-tmax, test_collision=False, **kwargs)
        t2, (x2, u2) = self.trajectory(-r0, St, *args, tmax=-tmax, test_collision=False, **kwargs)
        t1, x1, u1 = [np.flipud(e) for e in [t1, x1, u1]]
        t2, x2, u2 = [np.flipud(e) for e in [t2, x2, u2]]

        from scipy.interpolate import CubicSpline
        sol_x1, sol_u1 = CubicSpline(t1, x1), CubicSpline(t1, u1)
        sol_x2, sol_u2 = CubicSpline(t2, x2), CubicSpline(t2, u2)

        t1_interp = np.linspace(t1[0], 0, N//2)
        t2_interp = np.linspace(0, t2[0], N-N//2)
        x1_interp, u1_interp = sol_x1(t1_interp), sol_u1(t1_interp)
        x2_interp, u2_interp = sol_x2(t2_interp), sol_u2(t2_interp)
        x_interp = np.concatenate([x1_interp, x2_interp])
        u_interp = np.concatenate([u1_interp, u2_interp])

        # Optionally find location of critical point where separatrix crosses u' = 0:
        if return_critical_point:
            # Find guess for location of intersection with nulclline.
            delta = u1 - self.u(x1)                       # delta = 0 on the nullcline
            iguess = np.where(np.diff(np.sign(delta)))[0] # sign change corresponds to intersection

            if len(iguess) == 0:
                # No critical point found!
                critical_point = (np.nan, np.nan)
            else:
                assert len(iguess) == 1
                tguess1 = t1[iguess[0]-1]
                tguess2 = t1[iguess[0]+1]

                # Optimise guess.
                from scipy.optimize import minimize_scalar
                residual = lambda t: (sol_u1(t) - self.u(sol_x1(t)))**2
                crit_t = minimize_scalar(residual, bounds=[tguess1, tguess2]).x
                critical_point = (sol_x1(crit_t), sol_u1(crit_t))

            return x_interp, u_interp, critical_point

        return x_interp, u_interp

    def limit_colliding_trajectory(self, St=1, N=1000, step=1e-12, **kwargs):
        r0 = np.array([-step, 0])
        t, (x, u) = self.trajectory(r0, St, **kwargs)

        from scipy.interpolate import CubicSpline
        sol_x, sol_u = CubicSpline(t, x), CubicSpline(t, u)
        t_interp = np.linspace(0, t[-0], N-1)
        x, u = sol_x(t_interp), sol_u(t_interp)
        return np.concatenate([[0], x]), np.concatenate([[0], u])


class BasePlanarFlowField:
    @property
    def on_axis(self):
        return PlanarOnAxisFlow(self)

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
        raise NotImplementedError

    def streamline_eom(self, x, y):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(x,y), self.v(x,y)

    def collision(self, x, y):
        """Function that changes sign to indicate when a collision occurs."""
        return -x

    def trajectory(self, r0, St, xmax=None, ymax=None,
                   test_collision=True, return_collision=False,
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
            test_collision: if True, will test for and terminate on collision with object.
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

        events = []

        if test_collision:
            collide = lambda t, x: self.collision(*x[:2])
            collide.terminal = True
            events += [collide]

        if xmax is not None:
            out_of_bounds1 = lambda t, x: x[0] - xmax
            out_of_bounds1.terminal = True
            events += [out_of_bounds1]

        if ymax is not None:
            out_of_bounds2 = lambda t, x: np.abs(x[1]) - ymax
            out_of_bounds2.terminal = True
            events += [out_of_bounds2]

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
        """Determine whether particle collides with obstacle under specified conditions."""
        _, _, collides = self.trajectory(*args, return_collision=True, **kwargs)
        return collides

    def capture_efficiency(self, x, St, niters=25, quiet=True, *args, return_bounds=False, yguess=1, message_header='', **kwargs):
        """Estimate efficiency of point particle capture.

        Args:
            x: initial condition for x.
            St: value of Stokes number we are evaluating efficiency at.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
            yguess: initial guess for collision value of y for iterative algorithm.
            message_header: preface to iteration updates if not quiet.
        Returns:
            Efficiency of particle capture efficiency.
        """

        collides = lambda y: self.does_collide([x, y], St, *args, **kwargs)

        if not collides(0):
            if not quiet: print('on-axis does not lead to collision, so efficiency is zero!')
            return 0

        ylow, yhigh = logical_refine(collides, niters=niters, quiet=quiet, xguess=yguess, message_header=message_header)
        if return_bounds: return ylow, yhigh
        else:
            yestimate = 0.5 * (ylow + yhigh)
            return 2*yestimate # capture efficiency covers both +ve and -ve y, so x2


    def critical_stokes(self, x, niters=50, quiet=True, *args, message_header='', **kwargs):
        """Estimate critical Stokes number above which point particle capture occurs.

        This requires solving the on-axis problem only.

        Args:
            x: initial condition for x.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
            message_header: preface to iteration updates if not quiet.
        """

        collides = lambda St: not self.on_axis.does_collide(x, St, *args, **kwargs)
        Stc_low, Stc_high = logical_refine(collides, niters=niters, quiet=quiet, message_header=message_header)
        return 0.5 * (Stc_low + Stc_high)


# Specialisations for inertial and non-inertial reference frames.

class PlanarFlowFieldInertial(BasePlanarFlowField):
    def eom(self, x, y, ux, uy, St):
        """Equation of motion for massive particle immersed in flow field with
        no inertial forces."""
        return ux, uy, -(ux - self.u(x,y))/St, - (uy - self.v(x,y))/St

class PlanarFlowFieldNonInertial(BasePlanarFlowField):
    def eom(self, x, y, ux, uy, St):
        """Equation of motion for massive particle immersed in flow field with
        limiting inertial forces."""
        return ux, uy, uy**2 - (ux - self.u(x,y))/St, -2*ux*uy - (uy - self.v(x,y))/St

    def trajectory(self, r0, *args, xmax=None, ymax=5, **kwargs):
        if xmax is None:
            eps = 1e-6
            xmax = r0[0] + eps

        return super().trajectory(r0, *args, xmax=xmax, ymax=ymax, **kwargs)

PlanarFlowField = PlanarFlowFieldNonInertial
