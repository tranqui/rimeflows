#!/usr/bin/env python3

# Copyright (c) 2025 Patrick B Warren <patrick.warren@stfc.ac.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

try: from .planar import PlanarFlowField
except: from planar import PlanarFlowField

class Hiemenz:
    """Hiemenz flow has (vx, vy) = (−\sqrt(nu*k) phi(eta), k*y*phi'(eta)) where
    eta = \sqrt(k/nu)*x is a scaled distance from the surface, nu is viscosity, and the
    function solves phi''' + phi phi'' − phi'^2 + 1 = 0 with phi(0) = phi'(0) = 0 and
    phi'(\infty) = 1 as boundary conditions; see H. Schlichting, "Boundary
    layer theory" (McGraw-Hill, New York, 1979).

    For eta -> 0 the function must be phi(eta) = beta*eta^2 to leading order because
    of the boundary conditions.  Here 2*beta = phi''(0) \approx 1.2326 is found by a
    shooting method (below); see also Table 5.1 in Schlichting (1979).
    For eta -> \infty the function becomes phi(eta) = eta - eta* to leading order, where
    eta* \approx 0.6479 is the displacement thickness in units of \sqrt(k/nu).

    The accurate values are 2*beta = 1.2325876568, eta* = 0.6479004744, for
    the default eta cutoff = 10, and tol = 1e-12.

    For the present purposes it's convenient to identify delta = \sqrt(nu/k) as a
    proxy for the boundary layer thickness, so that the scaled variable
    eta = x/delta.  With this (vx, vy) = (−k*phi*delta, k*phi'*y)

    In the limit eta -> \infty therefore, (vx, vy) = (-k*x, k*y) as expected.
    Conversely in the limit eta -> 0, (vx, vy) = (-beta*k*x^2/delta, 2*beta*k*x*y/delta).

    To match to the potential flow around a cylinder, k = 2*U/R, where U
    is the flow in far-field and R is the cylinder radius.  Then delta/R =
    \sqrt(nu/2UR) = 1 / \sqrt(2*Re) where Re = U*R/nu = (R/delta)^2/2 is the Reynolds
    number.  We can work in units where R = 2U = 1, so that k = 1, but
    then the Stokes number St = (2*U)*m/(R*xi) should be divided by 2.  This
    does not affect the Reynolds number calculation.

    The class here instantiates a solution for the Hiemenz function phi(eta).
    It provides a method for evaluating phi, phi', and phi'' as a function of
    eta, and the above-mentioned attributes.
    """

    def ode(self, eta, y): # encodes the above ODE
        phi0, phi1, phi2 = y
        phi3 = phi1**2 - phi0 * phi2 - 1
        return phi1, phi2, phi3

    def root(self, phi2init, cutoff, tol): # return phi'(\infty) - 1 where '\infty' = cutoff (eg 10.0)
        y0 = np.array([0.0, 0.0, phi2init]) # the ODE is solved by the default method 'RK45'
        soln = solve_ivp(self.ode, [0, cutoff], y0, t_eval=[0, cutoff], atol=tol, rtol=tol)
        return soln.y[1][1] - 1 # endpoint value of phi', minus one

    def __init__(self, cutoff=10.0, init_bracket=[1.22, 1.24], tol=1e-12):

        soln = root_scalar(self.root, args=(cutoff, tol), bracket=init_bracket, xtol=tol, rtol=tol)
        self.two_beta = soln.root # extract the solved value of phi''(0) = 2*beta

        y0 = np.array([0.0, 0.0, self.two_beta]) # construct an interpolating function (dense_output) [re-use soln]
        soln = solve_ivp(self.ode, [0, cutoff], y0, t_eval=[0, cutoff], dense_output=True, atol=tol, rtol=tol)

        self.tol = tol # keep this as an attribute for reference
        self.cutoff = soln.t[1] # ditto; this was the maximum value of the argument
        self.eta_star = soln.t[1] - soln.y[0][1] # the displacement thickness estimated from the function value here
        self.sol = soln.sol # the interpolating function used in the function evaluations

    def __call__(self, eta):
        """Value of Hiemenz function and first two derivatives."""

        phi, phip, phipp = self.sol(eta) # this extrapolates ..
        if np.isscalar(eta):
            if eta > self.cutoff: # .. but we should use the known asymptotic form phi, phi', phi'' = eta - eta*, 1, 0
                phi, phip, phipp = eta - self.eta_star, 1.0, 0.0
        else: # vectorised situation: deal with the case where eta is a numpy array
            phi[eta > self.cutoff] = eta[eta > self.cutoff] - self.eta_star
            phip[eta > self.cutoff] = 1.0
            phipp[eta > self.cutoff] = 0.0
        return phi, phip, phipp

class FlowField(PlanarFlowField):
    def __init__(self, reynolds, *args, k=1, **kwargs):
        self.reynolds = reynolds
        if reynolds == np.inf: self.delta = 0
        elif reynolds == 0: self.delta = np.inf
        else: self.delta = np.sqrt(1/(2*reynolds))

        self.hiemenz = Hiemenz(*args, **kwargs)
        self.k = k

    def u(self, x, y):
        phi, _, _ = self.hiemenz(x / self.delta)
        return -self.k * phi * self.delta

    def v(self, x, y):
        _, phip, _ = self.hiemenz(x / self.delta)
        return self.k * phip * y

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    L = 3
    x0 = 1.5*L

    for Re in [1e-4, 1e4]:
        flow = FlowField(Re)

        plt.figure()
        for y0 in np.linspace(0, L**0.5, 11)**2:
            x, y = flow.streamline([x0, y0], tmax=1e6)
            plt.plot(x, y, 'k-', lw=0.5)
            plt.plot(x, -y, 'k-', lw=0.5)

        plt.fill_between([-0.5, 0], -L, L, fc='r', zorder=10)

        plt.xlim([-0.5, 1.5*L])
        plt.ylim([-L, L])

    plt.show()
