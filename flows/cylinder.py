#!/usr/bin/env python3

import numpy as np
from scipy import integrate

try: from .refine import logical_refine
except: from refine import logical_refine

class CylindricalFlowField:
    @property
    def default_starting_distance(self):
        """Default initial condition for stagnation point flows (in units of fibre radius)."""
        return 2

    def f(self, r):
        """Radial part of streamfunction."""
        raise NotImplementedError

    def fp(self, r):
        """First (radial) derivative of radial part of streamfunction."""
        raise NotImplementedError

    def u(self, r, th):
        """Radial component of flow field."""
        return np.cos(th) * self.f(r)/r

    def v(self, r, th):
        """Angular component of flow field."""
        return -np.sin(th) * self.fp(r)

    def u_cartesian(self, r, th):
        """Cartesian components of flow field."""
        u = self.u(r, th)
        v = self.v(r, th)
        ux = u*np.cos(th) - v*np.sin(th)
        uy = u*np.sin(th) + v*np.cos(th)
        return ux, uy

    def u_cartesian2(self, x, y):
        """Cartesian components of flow field from cartesian position."""
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return self.u_cartesian(r, theta)

    def eom(self, r, th, ur, ut, St):
        """Equation of motion for massive particle immersed in flow field."""
        return ur, ut/r, -(ur - self.u(r,th))/St + ut**2/r, -(ut - self.v(r,th))/St - ur*ut/r

    def streamline_eom(self, r, th):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(r,th), self.v(r,th)/r

    def trajectory(self, r0, St, tmax=1e2, xmax=None, ymax=None, return_collision=False,
                   max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field (in polar coordinates).

        Args:
            r0: initial condition (n-component vector, n=2 or 4).
                  When n=2 only r0 = (r, th) are specified, and the initial particle velocity will
                    be set to be aligned with the flow field.
                  For St > 0, n=4 one can optionally set also the initial velocities, i.e.
                    r0 = (r, th, u, v).
            St: Stokes number
            tmax: max time before aborting integration.
            xmax: max value of x before terminating integration.
                  If None, will be set to the negative of the initial x.
            ymax: max value of y before terminating integration.
            return_collision: If True, will return additional boolean stating whether a collision
                  event occurred (i.e. the particle impacted on the cylinder to terminate).
            max_step: maximum step during integration (sent to scipy.integrate.solve_ivp).
            rtol: relative error tolerance during integration (sent to scipy.integrate.solve_ivp).
            atol: absolute error tolerance during integration (sent to scipy.integrate.solve_ivp).
            kwargs: additional arguments to send to scipy.integrate.solve_ivp.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            x: coordinates at each frame (with dimensionality 2 for St = 0 or 4 for St > 0).
            collides (optional): True/False depending on whether particle terminates on collision
                                   with cylinder. Only returned if return_collision set to True.
        """

        if St == 0:
            assert len(r0) == 2
            eom = lambda t, x: self.streamline_eom(*x)
        else:
            if len(r0) == 2: r0 = np.concatenate([r0, [self.u(*r0), self.v(*r0)]])
            else: assert len(r0) == 4
            eom = lambda t, x: self.eom(*x, St=St)

        r = lambda x: x[0]
        theta = lambda x: x[1]

        collide = lambda t, x: 1 - r(x)
        collide.terminal = True

        if xmax is None:
            x0 = r0[0] * np.cos(r0[1])
            xmax = -x0
        out_of_bounds_x2 = lambda r, th: r*np.cos(th) - xmax
        out_of_bounds_x = lambda t, x: out_of_bounds_x2(r(x), theta(x))
        out_of_bounds_x.terminal = True

        events = (collide, out_of_bounds_x)

        if ymax is not None:
            out_of_bounds_y2 = lambda r, th: np.abs(r*np.sin(th)) - ymax
            out_of_bounds_y = lambda t, x: out_of_bounds_y2(r(x), theta(x))
            out_of_bounds_y.terminal = True
            events += (out_of_bounds_y,)

        with np.errstate(divide='ignore', invalid='ignore'):
            trajectory = integrate.solve_ivp(eom, [0, tmax], r0, events=events,
                                             vectorized=True, dense_output=True,
                                             max_step=max_step, rtol=rtol, atol=atol,
                                             **kwargs)

        if return_collision: return trajectory.t, trajectory.y, trajectory.t_events[0].size > 0
        else: return trajectory.t, trajectory.y

    def trajectory_cartesian(self, r0, *args, **kwargs):
        """Obtain trajectory for particle immersed in flow field (in Cartesian coordinates).

        Args:
            r0: reference point [2d r0 = (x, y)] along streamline. Streamline will be obtained
                  by following flow in both directions from this point.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            coords: coordinates of particle at each time frame.
            collides (optional): True/False depending on whether particle terminates on collision
                                   with cylinder. Only returned if return_collision set to True.
        """

        if len(r0) == 2:
            x, y = r0
            r = np.sqrt(x**2 + y**2)
            th = np.arctan2(y, x)

            ret = self.trajectory([r, th], *args, **kwargs)
            (t, coords), extra = ret[:2], ret[2:]
            r, th, ur, uth = coords

            try:
                r, th = coords
                x, y = r * np.cos(th), r * np.sin(th)
                coords = np.array((x, y))
            except:
                r, th, ur, uth = coords
                x, y = r * np.cos(th), r * np.sin(th)
                ux = ur*np.cos(th) - uth*np.sin(th)
                uy = ur*np.sin(th) + uth*np.cos(th)
                coords = np.array((x, y, ux, uy))

        else:
            assert len(r0) == 4
            x, y, ux, uy = r0
            r = np.sqrt(x**2 + y**2)
            th = np.arctan2(y, x)

            ret = self.trajectory([r, th], *args, **kwargs)
            (t, coords), extra = ret[:2], ret[2:]
            r, th, ur, uth = coords

            x, y = r * np.cos(th), r * np.sin(th)
            ux = ur*np.cos(th) - uth*np.sin(th)
            uy = ur*np.sin(th) + uth*np.cos(th)
            coords = np.array((x, y, ux, uy))

        return t, coords, *extra

    def streamline(self, r0, tmax=1e2, **kwargs):
        """Obtain the streamline for the flow field.

        Args:
            r0: reference point [2d r0 = (r, th)] along streamline. Streamline will be obtained
                  by following flow in both directions from this point.
            tmax: maximum time to follow streamline for.
        Returns:
            r: radii along streamline.
            th: corresponding theta along streamline.
        """

        t1, (r1, th1) = self.trajectory(r0, St=0, tmax=-tmax, **kwargs)
        t2, (r2, th2) = self.trajectory(r0, St=0, tmax=tmax, **kwargs)

        r = np.concatenate([np.flipud(r1), r2])
        th = np.concatenate([np.flipud(th1), th2])

        return r, th

    def streamline_cartesian(self, r0, *args, **kwargs):
        """Obtain the streamline for the flow field in a Cartesian basis.

        Args:
            r0: reference point [2d r0 = (x, y)] along streamline. Streamline will be obtained
                  by following flow in both directions from this point.
        Returns:
            x: x-coordinate along streamline.
            y: corresponding y-coordinates along streamline.
        """

        x, y = r0
        r = np.sqrt(x**2 + y**2)
        th = np.arctan2(y, x)

        r, th = self.streamline([r, th], *args, **kwargs)
        x, y = r * np.cos(th), r * np.sin(th)
        return np.array((x, y))

    def does_collide(self, *args, **kwargs):
        """Determine whether particle collides with cylinder under specified conditions."""
        _, _, collides = self.trajectory(*args, return_collision=True, **kwargs)
        return collides

    def does_collide_cartesian(self, *args, **kwargs):
        """Determine whether particle collides with cylinder under specified conditions."""
        _, _, collides = self.trajectory_cartesian(*args, return_collision=True, **kwargs)
        return collides

    def capture_efficiency(self, x, St, niters=25, cartesian=False, quiet=True, *args, **kwargs):
        """Estimate efficiency of point particle capture.

        Args:
            x: initial condition for r or x depending on whether we are using Cartesian coordinates or not.
            St: value of Stokes number we are evaluating efficiency at.
            cartesian: if True, will assume initial conditions are specified in Cartesian coordinates, otherwise polars.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
        Returns:
            Efficiency of particle capture efficiency.
        """

        if cartesian:
            r = lambda y: [x, y]
            does_collide = self.does_collide_cartesian
        else:
            r = lambda y: [x, y + np.pi]
            does_collide = self.does_collide
        collides = lambda y: does_collide(r(y), St, *args, **kwargs)

        if not collides(0):
            if not quiet: print('on-axis does not lead to collision, so efficiency is zero!')
            return 0

        ylow, yhigh = logical_refine(collides, niters=niters, quiet=quiet)
        return 0.5 * (ylow + yhigh)

    def on_axis_eom(self, r, ur, St):
        """Equation of motion for massive particle immersed in flow field approaching on-axis."""
        return ur, -(ur - self.u(r, np.pi))/St

    def on_axis_streamline_eom(self, r):
        """Equation of motion for inertialess particle immersed in flow field."""
        return self.u(r, np.pi)

    def on_axis_trajectory(self, r0, St, tmax=1e2, return_collision=False,
                           max_step=np.inf, rtol=1e-12, atol=1e-12, **kwargs):
        """Obtain trajectory for particle immersed in flow field (in polar coordinates)
            that's approaching the stagnation point on-axis (theta = \pi).

        Args:
            r0: initial condition (n-component vector, n=1 or 2).
                  When n=1 only r0 = r are specified, and the initial particle velocity will
                    be set to be aligned with the flow field.
                  For St > 0, n=2 one can optionally set also the initial velocity, i.e.
                    r0 = (r, u).
            St: Stokes number
            tmax: max time before aborting integration.
            return_collision: If True, will return additional boolean stating whether a collision
                  event occurred (i.e. the particle impacted on the cylinder to terminate).
            max_step: maximum step during integration (sent to scipy.integrate.solve_ivp).
            rtol: relative error tolerance during integration (sent to scipy.integrate.solve_ivp).
            atol: absolute error tolerance during integration (sent to scipy.integrate.solve_ivp).
            kwargs: additional arguments to send to scipy.integrate.solve_ivp.
        Returns:
            t: time sequence for each sampled frame along trajectory.
            x: coordinates at each frame (with dimensionality 1 for St = 0 or 2 for St > 0).
            collides (optional): True/False depending on whether particle terminates on collision
                                   with cylinder. Only returned if return_collision set to True.
        """

        try: len(r0)
        except: r0 = [r0]

        if St == 0:
            assert len(r0) == 1
            eom = lambda t, x: self.on_axis_streamline_eom(*x)
        else:
            if len(r0) == 1: r0 = np.concatenate([r0, [self.u(*r0, np.pi)]])
            else: assert len(r0) == 2
            eom = lambda t, x: self.on_axis_eom(*x, St=St)

        collide = lambda t, x: 1 - x[0]
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
        """Determine whether particle collides with cylinder under specified conditions."""
        _, _, collides = self.on_axis_trajectory(*args, return_collision=True, **kwargs)
        return collides

    def critical_stokes(self, x, niters=25, quiet=True, *args, **kwargs):
        """Estimate critical Stokes number above which point particle capture occurs.

        This requires solving the on-axis problem only.

        Args:
            x: initial condition for x.
            niters: number of refinement iterations to estimate Stokes number.
            quiet: if True, will suppress iteration updates.
        """

        collides = lambda St: not self.on_axis_does_collide(x, St, *args, **kwargs)
        Stc_low, Stc_high = logical_refine(collides, niters=niters, quiet=quiet)
        return 0.5 * (Stc_low + Stc_high)
