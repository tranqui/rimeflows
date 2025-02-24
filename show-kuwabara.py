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

# Note this script largely repeats the content of show-flow.py, but specialised
# to the Kuwabara use so it can be called inside a Jupyter notebook rather
# than from the command-line. This is for the purpose of tinkering with a
# figure for our manuscript.

import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from palettes import PurpleGold as palette
#from palettes import BlueRed as palette

streamline_line_width = 0.5
trajectory_line_width = 0.5

from flows import kuwabara

efficiency_data_path = 'data/efficiency_kuwabara_alpha=0.15.csv'

alpha = 0.15
flow = kuwabara.FlowField(alpha)

xmin = -flow.l
xmax = flow.l
ymin = -flow.l
ymax = flow.l

delta = 0.5
calc_xmin = xmin - delta
calc_xmax = xmax + delta
calc_ymin = ymin - delta
calc_ymax = ymax + delta


figsize = (3.375, 4.25)
# plt.figure(figsize=figsize)
# ax_eff = plt.gca()
fig, (ax_flow, ax_eff) = plt.subplots(nrows=2, figsize=figsize,
                                height_ratios=[1, 1])

St, dSt, eff = np.genfromtxt(efficiency_data_path, skip_header=1).T
Stc = St[0] - dSt[0]
ax_eff.loglog(dSt, eff, 'o', mfc='w')
ax_eff.set_xlabel(r'$\mathrm{St} - \mathrm{St}_\mathrm{c}$')
ax_eff.set_ylabel(r'$\lambda$')

#figsize = tuple(0.1*s for s in figsize)
figsize = (0.35, 0.35)
# lin = ax_eff.inset_axes([0.125, 0.975 - figsize[1], figsize[0], figsize[1]])
lin = ax_eff.inset_axes([0.96 - figsize[0], 0.175, figsize[0], figsize[1]])
lin.plot(St, eff, 'o', mfc='w')
lin.set_xlabel(r'$\mathrm{St}$', labelpad=-5)
lin.set_ylabel(r'$\lambda$', labelpad=-5)

lin.set_xlim([0., 2])
lin.set_xticks([0, 2], minor=False)
lin.set_xticks([1], minor=True)

ax_eff_xlim = ax_eff.get_xlim()
ax_eff_ylim = ax_eff.get_ylim()
lin_xlim = lin.get_xlim()
lin_ylim = lin.get_ylim()

# Fit the straight line of gradient 1/2 to obtain the intercept, and plot this best fit.
from scipy.optimize import curve_fit
f = lambda x,c: c*x**0.5
c, cerr = curve_fit(f, dSt[:5], eff[:5], p0=0.1)
c = c[0]
x = np.geomspace(1e-8, 1e8, 100)
pl1, = ax_eff.plot(x, f(x, c), '-k', zorder=-1)
pl2, = lin.plot(Stc + x, f(x, c), '-k', zorder=-1)

ax_eff.set_xlim(ax_eff_xlim)
ax_eff.set_ylim(ax_eff_ylim)
lin.set_xlim(lin_xlim)
lin.set_ylim(lin_ylim)

# Label location of critical Stokes
lin.annotate(r'$\mathrm{St}_\mathrm{c}$',
                xytext=(Stc, 0.75*lin_ylim[-1]),
                xy=(Stc, 0),
                ha='center', va='bottom',
                arrowprops=dict(lw=0.5, arrowstyle='-', ls='--',
                                ec=palette.critical_trajectory_colour),
                zorder=-1)

# A triangle to show gradient of straight line fit (in log space).
x = 3e-3
from mpltools import annotation
annotation.slope_marker((x, f(x, 0.75*c)), (1, 2),
                        poly_kwargs=dict(ec='black', fill=False, lw=0.5),
                        text_kwargs=dict(fontsize=8))

ax_flow.set_aspect(1)

St = 2
all_y0 = np.linspace(calc_ymin, calc_ymax, 21)
for y0 in all_y0:
    r0 = (xmin, y0)

    if abs(y0) < 1e-2:
        pl, = ax_flow.plot(np.array([xmin, xmax]), np.array([y0, y0]),
                           lw=streamline_line_width,
                           c=palette.streamline_colour, zorder=-1)
    else:
        x, y = flow.streamline_cartesian(r0, xmax=calc_xmax, ymax=calc_ymax)

        r2 = x**2 + y**2
        select = r2 <= flow.l**2
        x = x[select]
        y = y[select]

        pl, = ax_flow.plot(x, y, '-', lw=streamline_line_width,
                           c=palette.streamline_colour, zorder=-1)

    add_arrow(pl, x=-1.5, zorder=-1)
    add_arrow(pl, x=1.5, zorder=-1)

eff = flow.capture_efficiency(flow.l, St, max_step=0.1)

thc = 0.5*eff
yc = flow.l * np.sin(thc)
yc -= 1e-6 # Ensure that we are slightly below transition so trajectory will actually terminate on surface.

x = -1.15*flow.l
ax_flow.text(x, 0, r'$\lambda$', ha='right', va='center', fontsize=12)
x = -1.1*flow.l
y = 1.2*yc
ax_flow.annotate('', (x, y), (x, -y), annotation_clip=False,
                 weight=0.1, fontsize=8,
                 arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<->'))


N = 11
th0_arr = np.pi - np.linspace(0, eff, N)
y0_arr = flow.l * np.sin(th0_arr)

# Draw critical trajectory.
y0 = yc
x0 = -np.sqrt(flow.l**2 - y0**2)
r0 = (x0, y0)
t, (x, y, _, _), collides = flow.trajectory_cartesian(r0, St, max_step=0.1, return_collision=True)

r2 = x**2 + y**2
select = r2 <= flow.l**2
x = x[select]
y = y[select]

ax_flow.plot(-1, 0, 'o', mfc='w', ms=2, mew=0.25, zorder=100) # forward stagnation point
ax_flow.plot(x, y, '--', lw=trajectory_line_width,
             c=palette.critical_trajectory_colour, zorder=1)
ax_flow.plot(x, -y, '--', lw=trajectory_line_width,
             c=palette.critical_trajectory_colour, zorder=1)

# Ignore critical trajectory from ax_eff render
mask = np.ones(y0_arr.size, dtype=bool)
mask[N//2] = False
y0_arr = y0_arr[mask]

for y0 in y0_arr:
    x0 = -np.sqrt(flow.l**2 - y0**2)
    r0 = (x0, y0)
    t, (x, y, _, _), collides = \
        flow.trajectory_cartesian(r0, St, max_step=0.1, return_collision=True)

    r2 = x**2 + y**2
    select = r2 <= flow.l**2
    x = x[select]
    y = y[select]

    c = palette.colliding_trajectory_colour if collides else \
        palette.noncolliding_trajectory_colour
    ax_flow.plot(x, y, lw=trajectory_line_width, c=c, zorder=1)
    ax_flow.plot(x, -y, lw=trajectory_line_width, c=c, zorder=1)

ax_flow.add_patch(Circle((0,0), flow.l, fc='None', ec='k', ls='dashed', lw=1, zorder=10))


eps = 5e-2
background_colour = palette.flow_colour
ax_flow.fill_between([xmin-eps, xmax+eps], ymin-eps, ymax+eps,
                      facecolor=background_colour, zorder=-10)
ax_flow.set_xlim([xmin-eps, xmax+eps])
ax_flow.set_ylim([ymin-eps, ymax+eps])
ax_flow.add_patch(Circle((0,0), 1, color=palette.cylinder_colour, lw=0))

ax_flow.axis('off')

axes = [ax_flow, ax_eff]
for ax, y, letter in zip(axes, [0.9, 0.45], 'ab'):
    label = ax.text(0.05, y, r'\textbf{%s}' % letter, transform=fig.transFigure, zorder=20,
                    fontsize=18, horizontalalignment='left', verticalalignment='bottom')
    label.set_in_layout(False)

plt.savefig('kuwabara-efficiency.pdf')

plt.show()
