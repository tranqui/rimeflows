#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# hiemenz.py -- calculate capture efficiency in Hiemenz flow.

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

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from utilities import Hiemenz, string_evaluate, str_or_infinity

parser = argparse.ArgumentParser(description='calculate collection efficiency')
parser.add_argument('--yout', default=3.0, type=float, help='escape distance, default 3.0')
parser.add_argument('--ymax', default=1.0, type=float, help='notional maximum in y, default 1.0')
parser.add_argument('--tmax', default=30.0, type=float, help='max time, default 30.0')
parser.add_argument('--x0', default=1.0, type=float, help='starting position in x, default 1.0')
parser.add_argument('--dy', default=1e-4, type=float, help='y-spacing, default 1e-4')
parser.add_argument('--k', default=1.0, type=float, help='value of k in the flow field, default 1.0')
parser.add_argument('--delta', default='0.2', help='value of δ/R in the flow field, default 0.2')
parser.add_argument('--reynolds', help='Reynolds number to compute δ/R in the flow field, if required')
parser.add_argument('--stokes', default='0.2,1.0,0.01', help='Stokes number range, default 0.5,1.0,0.01')
parser.add_argument('--method', default='LSODA', help='solver method, default LSODA')
parser.add_argument('--tol', default=1.0e-12, type=float, help='solver tolerance, default 1.0e-12')
parser.add_argument('--squared', action='store_true', help='include y^2')
parser.add_argument('--only', action='store_true', help='(with --include) only plot y^2')
parser.add_argument('-o', '--output', help='write results as tab-separated values to a file')
args = parser.parse_args()

k = args.k
b = args.yout

if args.reynolds is not None:
    re = string_evaluate(args.reynolds)
    δ = 0 if re == np.inf else np.inf if re == 0 else np.sqrt(1/(2*re))
else:
    δ = string_evaluate(args.delta)
    re = 0 if δ == np.inf else np.inf if δ == 0 else 1/(2*δ**2)

# print('Computing for δ, Re = ', δ, re)

# Hiemenz flow has vx, vy = − √(νk) ϕ(η), k y ϕ'(η) where η = √(k/ν) x
# is a scaled distance from the surface, ν is viscosity, and the
# function solves ϕ‴ + ϕ ϕ" − ϕ'² + 1 = 0 with ϕ(0) = ϕ'(0) = 0 and
# ϕ'(∞) = 1 as boundary conditions; see H. Schlichting, "Boundary
# layer theory" (McGraw-Hill, New York, 1979).

# For the present purposes it's convenient to identify δ = √(ν/k) as a
# proxy for the boundary layer thickness, so that the scaled variable
# η = x/δ.  With this (vx, vy) = (−kϕδ, kϕ'y)

# For η → 0 the function must be ϕ(η) = βη² to leading order because
# of the boundary conditions.  Here 2β = ϕ"(0) ≈ 1.2326 is found by a
# shooting method; see also Table 5.1 in Schlichting (1979).
# Conversely for η → ∞ one has ϕ(η) = η - η* to leading order, where
# η* ≈ 0.6479 is the displacement thickness in units of δ.

# The accurate values are 2β = 1.2325876568, η* = 0.6479004744, for
# the default numerics.

# In the limit η → ∞ therefore, (vx, vy) = (-kx, ky) as expected.
# Conversely in the limit η → 0, (vx, vy) = (-βkx²/δ, 2βkxy/δ).

# To match to the potential flow around a cylinder, k = 2U/R, where U
# is the flow in far-field and R is the cylinder radius.  Then δ/R =
# √(ν/2UR) = 1 / √(2Re) where Re = UR/ν = (R/δ)²/2 is the Reynolds
# number.  We can work in units where R = 2U = 1, so that k = 1, but
# then the Stokes number St = (2U)m/Rξ should be divided by 2.  This
# does not affect the Reynolds number calculation.

# More precisely, let Newton's EoM be mx" + ξ(x' - u) = 0, where m is
# the mass of the particle and ξ is the drag coefficient.  Defining a
# transit time τ = R/γU, with γ an arbitrary numerical prefactor, and
# Stokes number St = mU/Rξ (which is the conventional definition),
# this can be written as γSt τ²x" + τx' − Rũ = 0 where ũ = u/γU.  Now,
# Hiemenz has u = (−k ϕδ, kϕ'y) where k = 2U/R from the potential
# flow. So Rũ = (−(kR/γU) ϕδ, (kR/γU) ϕ'y).  By choosing γ = 2 we can
# make kR/γU go away so that finally γSt τ²x" + τx' − Rũ = 0 where Rũ
# = (−ϕδ, ϕ'y).  This shows that the Stokes number in the case where
# 'k = 1' is assumed is actually γ (= 2) times the true Stokes number.
# Finally, dividing through by R gives γSt τ²x"/R + τx'/R − ũ = 0
# where ũ = (−ϕδ/R, ϕ'y/R).  The argument of the Hiemenz function is
# unaffected by this since η = x/δ = (x/R) / (δ/R).  Similarly, the
# Reynolds number remains as Re = (R/δ)²/2.

# Instantiate the solution with the default numerics.

hiemenz = Hiemenz()

β = 0.5 * hiemenz.two_beta # capture the limiting behavior as η → 0

# Solve directly St [x" - (y')²] + x' + k δ ϕ = 0,
#             && St [y" + 2x'y'] + y' = k y ϕ'.

def hiemenz_flow(x, y):
    if δ == np.inf: # limiting behaviour as δ → ∞, rescaled by δ
        ux, uy = -β*k*x**2, 2*β*k*x*y
    elif δ == 0: # limiting behaviour as δ → 0
        ux, uy = -k*x, k*y
    else:
        η = x / δ
        φ, dφdη, _ = hiemenz.func(η)
        ux, uy = -k*ϕ*δ, k*dφdη*y
    return ux, uy

def eom(t, ξ):
    x, y, vx, vy = ξ
    ux, uy = hiemenz_flow(x, y)
    return vx, vy, -(vx - ux)/St + vy**2, -(vy - uy)/St - 2*vx*vy

def collide(t, ξ):
    x, y, vx, vy = ξ
    return x

def escape(t, ξ):
    x, y, vx, vy = ξ
    return y**2 - b**2

collide.terminal = True
escape.terminal = True

res = []
x0 = args.x0
y_min = 0.0

for St in np.arange(*eval('(' + args.stokes + ')')):
    for y0 in np.arange(y_min, args.ymax, args.dy):
        ux0, uy0 = hiemenz_flow(x0, y0)
        soln = solve_ivp(eom, [0, args.tmax], np.array([x0, y0, ux0, uy0]),
                         events=[collide, escape], method=args.method, atol=args.tol, rtol=args.tol)
        collision = soln.t_events[0].size > 0
        if not collision: # the trajectory escaped
            break
    if y0 > 0: # the previous trajectory hit
        y_prev = y0 - args.dy # previous trajectory starting height
        y_min = max(y_min, y_prev) # new initial starting height for arange above
        y_mean = 0.5*(y_prev + y0) # consensus for crossover for this value of Stokes
        res.append((St, y_mean))

df = pd.DataFrame(res, columns=['St', 'y'])

df.set_index('St', inplace=True)

dbyr = str_or_infinity(δ)
reynolds = str_or_infinity(re)
st_crit = df.index[0] # the first entry

if args.output:

    with open(args.output, 'w') as f:
        print('## δ/R, Re, St_crit =\t', '\t'.join([dbyr, reynolds, str(st_crit)]), file=f)
        print('#' + df.index.name + '\t' + '\t'.join(df.columns), file=f)
        print(df.to_csv(sep='\t', header=False), end='', file=f)

    print(f'Data for δ/R = {dbyr}, Re = {reynolds} written to {args.output}')

else:

    print(f'Plotting data for δ/R = {dbyr}, Re = {reynolds}, St_crit = {str(st_crit)}')

    if args.squared:
        df['y^2'] = df['y']**2
        if args.only:
            df[['y^2']].plot(style='x')
        else:
            df[['y', 'y^2']].plot(style='x')
    else:
        df[['y']].plot(style='x')

    plt.show()
