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
        return -self.k * x**(2-self.epsilon)

    def v(self, x, y):
        raise NotImplementedError
        return 2*self.k * x * y

class OnAxis:
    def __init__(self, epsilon, k=1):
        self.epsilon = epsilon
        self.k = k

    def u(self, x):
        #return -self.k*x**(2-self.epsilon)
        return -self.k * ( (1-self.epsilon) * x**2 + self.epsilon * x )

    def trajectory(self, r0, step=1e-12, rtol=1e-12, tmax=1e2, xmin=-5, xmax=5, vmin=-5, vmax=5):
        terminate = lambda t,r: r[0] < xmax and r[0] > xmin and r[1] < vmax and r[1] > vmin
        terminate.terminal = True
        terminate.direction = -1

        xp = lambda x,v: v
        vp = lambda x,v: self.u(x) - v
        rp = lambda t,r: [xp(*r), vp(*r)]

        return solve_ivp(rp, (0, tmax), r0, first_step=step, rtol=rtol, dense_output=True, events=terminate)

    def separatrix(self, *args, step=1e-12, tmax=1e2, N=2000, return_critical_point=False, **kwargs):
        r0 = step * np.array([1, -1])
        sol1 = self.trajectory(r0, *args, tmax=-tmax, **kwargs)
        sol2 = self.trajectory(-r0, *args, tmax=-tmax, **kwargs)

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

    def nullcline(self, xmax=5):
        x = np.linspace(0, xmax, 1000)
        v = self.u(x)
        return x, v

    def streamline(self, x0, v0, *args, tmax=1e2, N=2000, **kwargs):
        r0 = np.array([x0, v0])
        sol1 = self.trajectory(r0, *args, tmax=-tmax, **kwargs)
        sol2 = self.trajectory(r0, *args, tmax=tmax, **kwargs)

        t1 = np.linspace(sol1.t[-1], 0, N//2)
        t2 = np.linspace(0, sol2.t[-1], N-N//2)
        x1, v1 = sol1.sol(t1)
        x2, v2 = sol2.sol(t2)
        x = np.concatenate([x1, x2])
        v = np.concatenate([v1, v2])

        return x, v

    def limit_unstable_streamline(self, *args, N=2000, step=1e-12, **kwargs):
        r0 = np.array([-step, 0])
        sol = self.trajectory(r0, *args, step=step, **kwargs)

        t = np.linspace(0, sol.t[-1], N-1)
        x, v = sol.sol(t)
        x = np.concatenate([[0], x])
        v = np.concatenate([[0], v])
        return x, v

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
