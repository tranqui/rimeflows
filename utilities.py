#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# utilities.py -- supporting code for filtration calculations.

# Copyright (c) 2022 Patrick B Warren <patrick.warren@stfc.ac.uk>

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

# Hiemenz flow has vx, vy = − √(νk) ϕ(η), k y ϕ'(η) where η = √(k/ν) x
# is a scaled distance from the surface, ν is viscosity, and the
# function solves ϕ‴ + ϕ ϕ" − ϕ'² + 1 = 0 with ϕ(0) = ϕ'(0) = 0 and
# ϕ'(∞) = 1 as boundary conditions; see H. Schlichting, "Boundary
# layer theory" (McGraw-Hill, New York, 1979).

# For η → 0 the function must be ϕ(η) = βη² to leading order because
# of the boundary conditions.  Here 2β = ϕ"(0) ≈ 1.2326 is found by a
# shooting method (below); see also Table 5.1 in Schlichting (1979).
# For η → ∞ the function becomes ϕ(η) = η - η* to leading order, where
# η* ≈ 0.6479 is the displacement thickness in units of √(k/ν).

# The accurate values are 2β = 1.2325876568, η* = 0.6479004744, for
# the default eta cutoff = 10, and tol = 1e-12.

# The class here instantiates a solution for the Hiemenz function ϕ(η).
# It provides a method for evaluating ϕ, ϕ', and ϕ" as a function of
# η, and the above-mentioned attributes.

class Hiemenz:

    def __init__(self, cutoff=10.0, init_bracket=[1.22, 1.24], tol=1e-12):

        def hiemenz_ode(eta, y): # encodes the above ODE
            phi0, phi1, phi2 = y
            phi3 = phi1**2 - phi0 * phi2 - 1
            return phi1, phi2, phi3

        def hiemenz_root(phi2init, cutoff, tol): # return ϕ'(∞) - 1 where '∞' = cutoff (eg 10.0)
            y0 = np.array([0.0, 0.0, phi2init]) # the ODE is solved by the default method 'RK45'
            soln = solve_ivp(hiemenz_ode, [0, cutoff], y0, t_eval=[0, cutoff], atol=tol, rtol=tol)
            return soln.y[1][1] - 1 # endpoint value of ϕ', minus one

        soln = root_scalar(hiemenz_root, args=(cutoff, tol), bracket=init_bracket, xtol=tol, rtol=tol)
        self.two_beta = soln.root # extract the solved value of ϕ"(0) = 2β

        y0 = np.array([0.0, 0.0, self.two_beta]) # construct an interpolating function (dense_output) [re-use soln]
        soln = solve_ivp(hiemenz_ode, [0, cutoff], y0, t_eval=[0, cutoff], dense_output=True, atol=tol, rtol=tol)

        self.tol = tol # keep this as an attribute for reference
        self.cutoff = soln.t[1] # ditto; this was the maximum value of the argument
        self.eta_star = soln.t[1] - soln.y[0][1] # the displacement thickness estimated from the function value here
        self.sol = soln.sol # the interpolating function used in the function evaluations

    def func(self, eta):
        '''method returns the Hiemenz function and first two derivatives, for a float or numpy array'''
        phi, phip, phipp = self.sol(eta) # this extrapolates ..
        if np.isscalar(eta):
            if eta > self.cutoff: # .. but we should use the known asymptotic form ϕ, ϕ', ϕ" = η - η*, 1, 0
                phi, phip, phipp = eta - self.eta_star, 1.0, 0.0
        else: # vectorised situation: deal with the case where η is a numpy array
            phi[eta > self.cutoff] = eta[eta > self.cutoff] - self.eta_star
            phip[eta > self.cutoff] = 1.0
            phipp[eta > self.cutoff] = 0.0
        return phi, phip, phipp

# Utility function used for plotting.

def make_ticks(x_min, x_max, Δx, δx):
    '''Make a set of major and minor ticks with the indicated spacing'''
    major = list(np.linspace(x_min, x_max, 1 + int((x_max - x_min)/Δx)))
    minor = list(set(np.linspace(x_min, x_max, 1 + int((x_max - x_min)/δx))) - set(major))
    return major, minor

# Evaluate a string argument to a number, or np.inf

def string_evaluate(arg):
    s = arg.lstrip('0') # avoid issues with strings that begin '00...' which might be evaluated to octal
    ans = 0 if len(s) == 0 else np.inf if s.lower().startswith('inf') else eval(s) if s[0].isdigit() else eval('0' + s)
    return ans

# Convert a number to a string, handling the np.inf case

def str_or_infinity(v):
    return '∞' if v == np.inf else str(v)

if __name__ == "__main__":

    hiemenz = Hiemenz() # report default results
    res = (hiemenz.two_beta, hiemenz.eta_star, hiemenz.cutoff, hiemenz.tol)
    print('Hiemenz defaults: 2β = %0.10f, η* = %0.10f, for cutoff = %g, tol = %g' % res)

