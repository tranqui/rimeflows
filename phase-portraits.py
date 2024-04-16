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

import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate, optimize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')
plt.rcParams['figure.constrained_layout.use'] = False

from palettes import PurpleGold as palette
#from palettes import BlueRed as palette

import argparse
parser = argparse.ArgumentParser(description='show phase portraits of toy models')
parser.add_argument('-o', '--output', type=str, help='output plot to file (e.g. pdf or png)')
args = parser.parse_args()

from flows import stokes, shm, power
stokes = stokes.FlowField().on_axis
shm = shm.FlowField().on_axis

figsize = 3.375 # (inches)
fig, big_axes = plt.subplots(ncols=3, figsize=(2*figsize, figsize), gridspec_kw={'hspace': 0})

axes = []
for ax in big_axes:
    div = make_axes_locatable(ax)
    new_ax = div.append_axes("bottom", size="100%", pad=0)#, sharex=ax)
    axes += [(ax, new_ax)]

axes = np.array(axes).T
ax1, ax2, ax3, ax4, ax5, ax6 = axes.T.flatten()

# div = make_axes_locatable(ax1)
# ax4 = div.append_axes("top", size="100%", pad=0.05)
# axes = [ax1, ax2, ax3, ax4]

#bbox = dict(pad=0.1, fc='None', ec='None')
bbox = dict(boxstyle='round', pad=0.1, fc='white', ec='none', alpha=0.75)

N = 20000
xmin = -0.5
xmax = 1.5
umin = -1
umax = 1
m = 5

u0 = umax
max_step = np.inf
tmax = 1e2
#max_step = 1e-2
#tmax = 1e2

for St,ax in zip([1, 1.4], (ax1, ax2)):
    # Find the location of the separatrix at our initial u:
    x, u = stokes.separatrix(St=St, max_step=max_step, tmax=tmax, N=N)
    x, u = x[:len(x)//2], u[:len(x)//2] # only want initial crossing on right side of manifold
    guess = np.argmin((u-u0)**2)
    fx = interpolate.interp1d(np.arange(len(x)), x, kind='quadratic')
    fu = interpolate.interp1d(np.arange(len(u)), u, kind='quadratic')
    separatrix_x = fx( optimize.minimize(lambda i: (fu(i) - u0)**2, guess).x )

    # Sample points on either side of the separatrix across the full visible range so we
    # fill space with stream lines
    eps = 2e-3
    dx = 0.15
    delta_x = np.arange(0, m*(xmax - xmin), dx)[1:]
    central_x = separatrix_x + eps
    x_range = np.hstack((central_x - delta_x, central_x, central_x + delta_x))
    for x0 in x_range:
        t, (x, u) = stokes.trajectory([x0, u0], St=St, max_step=max_step, tmax=tmax, test_collision=False,
                                       xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))

        if x0 < separatrix_x and x[-1] > 0: c = palette.noncolliding_trajectory_colour
        else: c = palette.colliding_trajectory_colour

        pl, = ax.plot(x, u, lw=0.5, c=c, zorder=-2)
        add_arrow(pl, y=0.8, zorder=-1)
        add_arrow(pl, y=0.1, zorder=-1)
        #add_arrow(pl, y=-0.41, zorder=-10)
        add_arrow(pl, y=-0.81, zorder=-1)

    x, u = stokes.limit_colliding_trajectory(St=St, tmax=tmax, max_step=max_step, N=N)
    #ax.plot(x, u, lw=0.5, c=palette.colliding_trajectory_colour, zorder=-8)

    x = np.linspace(0, xmax, 100)
    u = stokes.u(x)
    ax.plot(x, u, '--', c=palette.nullcline_colour)

    x, u, (xc, uc) = stokes.separatrix(St=St, max_step=max_step, tmax=tmax, N=N, return_critical_point=True)
    index = np.where(np.logical_and(x > 0, u < 0))[0][-1]
    pl, = ax.plot(x[index:]-1e-2, u[index:], c=palette.colliding_trajectory_colour, lw=0.5)
    add_arrow(pl, y=0.1, zorder=-1, direction='backward')
    #ax1.plot(x, u, c=palette.separatrix_colour)
    ax.plot(0, 0, 'o', c=palette.separatrix_colour, mfc='w', zorder=20)
    ax.plot(xc, uc, 'o', c=palette.separatrix_colour, mfc=palette.critical_trajectory_colour, zorder=20)
    #ax1.plot(0, 0, 'o', c=palette.separatrix_colour, mfc=palette.separatrix_colour)

    # xs = np.linspace(0, 1/8, 100)
    # us = (1 - np.sqrt(1 - 8*xs) - 12*xs + 8*xs*np.sqrt(1 - 8*xs)) / 24
    # ax1.plot(xs, us, 'b-.', lw=1, label='slow manifold', zorder=5)

    ax.fill_between([xmin, xmax], umin, umax, facecolor=palette.colliding_manifold_colour, zorder=-10)
    ax.fill_between(x[x > 0], 0, u[x > 0], facecolor=palette.noncolliding_manifold_colour, zorder=-5)
    ax.fill_between(x[x < 0.1], 0, u[x < 0.1], facecolor=palette.noncolliding_manifold_colour, zorder=-5)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([umin, umax])

eps = 1e-6

xmin = -0.5
xmax = 1.5
umin = -4
umax = 1
m = 3

St = 0.2
u0 = 0

c1 = -(1 + np.sqrt(1 - 4*St)) / (2*St)
c2 = -(1 - np.sqrt(1 - 4*St)) / (2*St)

for x0 in np.linspace(xmin, xmax, 11):
    if np.abs(x0) < eps: continue

    for u0 in [umin, umax]:
        stable = u0 > 0 or x0 > u0/c1
        if stable: c = palette.noncolliding_trajectory_colour
        else: c = palette.colliding_trajectory_colour

        _, (x, u) = shm.trajectory([x0, u0], St=St, test_collision=False,
                                    xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
        pl, = ax3.plot(x, u, lw=0.5, c=c, zorder=-2)
        if u0 > 0: add_arrow(pl, y=0.5, zorder=-1)
        else:
            add_arrow(pl, y=-1.5, zorder=-1)
            add_arrow(pl, y=-3, zorder=-1)


x0 = xmax
for u0 in np.linspace(-3.25, -1.5, 4):
    stable = u0 > 0 or x0 > u0/c1
    if stable: c = palette.noncolliding_trajectory_colour
    else: c = palette.colliding_trajectory_colour

    _, (x, u) = shm.trajectory([x0, u0], St=St, test_collision=False,
                                xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
    pl, = ax3.plot(x, u, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=-3, zorder=-1)
    if u0 < c2*x0:
        add_arrow(pl, y=-1.5, zorder=-1)

#ax3.fill_between([xmin, xmax], umin, umax, facecolor=palette.colliding_manifold_colour, zorder=-10)
#ax3.fill_between([, R], -R, R, facecolor=palette.noncolliding_manifold_colour, zorder=-5)
ax3.set_xlim([xmin, xmax])
ax3.set_ylim([umin, umax])

c1 = -(1 + np.sqrt(1 - 4*St)) / (2*St)
c2 = -(1 - np.sqrt(1 - 4*St)) / (2*St)
x1 = np.linspace(xmin, 0, 10)
x2 = np.linspace(0, xmax, 10)
# Relevant separatrices
#ax3.plot(x2, c1*x2, '-', c=palette.separatrix_colour)
#ax3.plot(x1, c2*x1, '-', c=palette.separatrix_colour)
# Irrelevant separatrices
#ax3.plot(x1, c1*x1, '-', c=palette.separatrix_colour)
#ax3.plot(x2, c2*x2, '-', c=palette.separatrix_colour)
ax3.fill_between(x1, umin, c2*x1, facecolor=palette.colliding_manifold_colour, zorder=-10)
ax3.fill_between(x1, umax, c2*x1, facecolor=palette.noncolliding_manifold_colour, zorder=-10)
ax3.fill_between(x2, umin, c1*x2, facecolor=palette.colliding_manifold_colour, zorder=-10)
ax3.fill_between(x2, umax, c1*x2, facecolor=palette.noncolliding_manifold_colour, zorder=-10)

R = 2
th0 = np.pi/3
St = 0.5
for RR in np.linspace(-R, R, 11):
    if np.abs(RR) < eps: continue
    x0, u0 = RR*np.cos(th0), RR*np.sin(th0)

    c = palette.colliding_trajectory_colour
    _, (x1, u1) = shm.trajectory([x0, u0], St=St, tmax=tmax, test_collision=False,
                                  xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
    _, (x2, u2) = shm.trajectory([x0, u0], St=St, tmax=-tmax, test_collision=False,
                                  xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
    x = np.concatenate([np.flipud(x2), x1])
    u = np.concatenate([np.flipud(u2), u1])
    pl, = ax4.plot(x, u, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=0.5, zorder=-1)
    add_arrow(pl, y=-0.5, zorder=-1)
    add_arrow(pl, y=-2.5, zorder=-1)

ax4.fill_between([xmin, xmax], umin, umax, facecolor=palette.colliding_manifold_colour, zorder=-5)
# ax4.set_xlim([-R/2, R/2])
# ax4.set_ylim([-R/2, R/2])
ax4.set_xlim([xmin, xmax])
ax4.set_ylim([umin, umax])


St = 1
xmin = -0.75
xmax = 5
umin = -2
umax = 2
m = 10
u0 = 0

for m,ax in zip([1.5, 0.5], [ax5,ax6]):
    flow = power.FlowField(m).on_axis

    x, u, (xc, uc) = flow.separatrix(St=St, max_step=max_step, tmax=tmax, N=N, return_critical_point=True)
    if np.abs(xc) < 1e-8: xc, uc = 0, 0
    separatrix_x = np.max(x)

    ax.plot(0, 0, 'o', c=palette.separatrix_colour, mfc='w', zorder=20)
    ax.fill_between([xmin, xmax], umin, umax, facecolor=palette.colliding_manifold_colour, zorder=-10)

    x_range = np.linspace(xmin, xmax, 60)

    if xc > 0:
        #ax.plot(x, u, c='k', lw=0.5)
        ax.plot(xc, uc, 'o', c=palette.separatrix_colour, mfc=palette.critical_trajectory_colour, zorder=20)
        ax.fill_between(x[x > 0], 0, u[x > 0], facecolor=palette.noncolliding_manifold_colour, zorder=-5)
        ax.fill_between(x[x < 0.1], 0, u[x < 0.1], facecolor=palette.noncolliding_manifold_colour, zorder=-5)

        eps = 1e-3
        dx = x_range[1] - x_range[0]
        delta_x = np.arange(0, m*(xmax - xmin), dx)[1:]
        central_x = separatrix_x + eps
        x_range = np.hstack((central_x - delta_x, central_x, central_x + delta_x))

    for x0 in x_range:
        _, (x1, u1) = flow.trajectory([x0, u0], St=St, tmax=tmax, test_collision=False,
                                       xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
        _, (x2, u2) = flow.trajectory([x0, u0], St=St, tmax=-tmax, test_collision=False,
                                       xmin=(m*xmin), xmax=(m*xmax), umin=(m*umin), umax=(m*umax))
        x = np.concatenate([np.flipud(x2), x1])
        u = np.concatenate([np.flipud(u2), u1])

        if xc == 0: c = palette.colliding_trajectory_colour
        else:
            if x0 < 0 or x0 > separatrix_x: c = palette.colliding_trajectory_colour
            else: c = palette.noncolliding_trajectory_colour

        pl, = ax.plot(x, u, lw=0.5, c=c)

        add_arrow(pl, y=0.4)
        add_arrow(pl, y=0.1, eps=1e-6)
        add_arrow(pl, y=-0.35)

    x = np.linspace(0, xmax, 1000)
    u = flow.u(x)
    ax.plot(x, u, '--', c=palette.nullcline_colour)

    ax.set_xlim([-0.25, 1])
    ax.set_ylim([-0.5, 0.5])
    ax.set_yticks(np.arange(-0.4, 0.41, 0.2))


for ax in [ax3, ax4]:
    ax.plot(0, 0, 'o', mfc='w', c=palette.separatrix_colour, zorder=20)

    ax.set_ylabel('$\dot{x}$')

    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, -x, '--', c=palette.nullcline_colour)

for ax in axes[-1]:
    ax.set_xlabel('$x$')

for ax in axes.flatten():
    ax.fill_between([-1e3, 0], -1e3, 1e3, facecolor='black', alpha=0.125, zorder=10)
    ax.axvline(x=0, lw=0.5)

label = ax1.text( 0.35, 0.8, 'noncolliding\ntrajectories', fontsize=6,
                  bbox=bbox, ha='center', va='center', transform=ax1.transAxes)
label.set_in_layout(False)
label = ax1.text( 0.8, 0.8, 'colliding\ntrajectories', fontsize=6,
                  bbox=bbox, ha='center', va='center', transform=ax1.transAxes)
label.set_in_layout(False)
#ax1.text(-0.06, -0.45, 'collision at\nfinite $t$', fontsize=8, rotation=90,
#         ha='center', va='center')
# label = ax1.text(0., -0.5, 'collision within\nfinite $t$', fontsize=8, rotation=90,
#                  ha='center', va='center')
# label.set_in_layout(False)

# label = ax2.text( 0.325, 0.5, r'$(x,\dot{x},t)^\top \to (\mathrm{St} x, \mathrm{St}^2 \dot{x},t/\mathrm{St})^\top$', fontsize=8,
#                   bbox=bbox, ha='center', va='center')
# label.set_in_layout(False)

# # # Annotation for direction of increasing Stokes number.
# # delta = 0.05
# # shift = 0.1
# # x = np.linspace(delta, xc - delta, 100)
# # u = stokes.u(x) - shift
# # pl, = ax1.plot(x[:-1], u[:-1], '-', c='black', lw=0.5)
# # ax1.annotate('',
# #              xytext=(x[-2], u[-2]),
# #              xy=(x[-1], u[-1]),
# #              arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='->'))
# #              #arrowprops=dict(arrowstyle="->", color=pl.get_colour())#,
# #              #size=10
# # ax1.text(0.5*xc, stokes.u(0.5*xc) - shift - 0.05, r'$\big\uparrow \mathrm{St}$', ha='center', va='top', fontsize=8)

# ax1.annotate(r'$x \sim 1/t$', (0, 0), (0.25, 0.1), #weight=0.1,
#              fontsize=10, ha='center', va='bottom',
#              bbox=bbox, arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='-'))
# ax3.annotate(r'$x \sim e^{-t}$', (0, 0), (0.5, 0.75), #weight=0.1,
#              fontsize=10, ha='center', va='top',
#              bbox=bbox, arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='-'))

label = plt.text(0.85, 0.12, r'$\mathrm{St} = 1$', transform=ax1.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)
label = plt.text(0.85, 0.12, r'$\mathrm{St} > 1$', transform=ax2.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)

label = plt.text(0.85, 0.12, r'$\mathrm{St} < \frac{1}{4}$', transform=ax3.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)
label = plt.text(0.85, 0.12, r'$\mathrm{St} > \frac{1}{4}$', transform=ax4.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)

label = plt.text(0.85, 0.12, r'$m = \frac{3}{2}$', transform=ax5.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)
label = plt.text(0.85, 0.12, r'$m = \frac{1}{2}$', transform=ax6.transAxes,
                 fontsize=10, bbox=bbox, ha='center', va='center')
label.set_in_layout(False)

for ax,l in zip(axes.flatten(), 'acebdf'):
    label = plt.text(-0.325, 0.8, (r'\textbf{%s}' % l), transform=ax.transAxes,
                     fontsize=18, ha='left', va='bottom')
    label.set_in_layout(False)

for ax in axes.flatten(): ax.set_ylabel('$\dot{x}$')
for ax in axes[0]: ax.set_xticklabels([])
for ax in axes[-1]: ax.set_label('$x$')

plt.subplots_adjust(bottom=0.125, top=0.925, left=0.1, right=0.975, wspace=0.4)

# Prevent overlap between ticks in top and bottom panels.
ax2.set_yticks([-1, -0.5, 0, 0.5])

ax1.set_title('$u = -x^2$')
ax3.set_title('$u = -x$')
ax5.set_title('$u = -|x|^m$, $\mathrm{St} = 1$')

if args.output: plt.savefig(args.output, pad_inches=0.1)

plt.show()
