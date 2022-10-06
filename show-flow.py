#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from palettes import PurpleGold as palette
#from palettes import BlueRed as palette

from flows import kuwabara, inviscid, rg, stokes, shm
flow_fields = {'kuwabara': kuwabara, 'inviscid': inviscid, 'rg': rg, 'stokes': stokes, 'shm': shm}

parser = argparse.ArgumentParser(description='make Kuwabara trajectory plot')
parser.add_argument('-f', '--flow', default='kuwabara', choices=flow_fields.keys(),
                    help='flow field to use (default=kuwabara)')
parser.add_argument('-e', '--efficiency', type=str,
                    help='path to data file containing efficiency vs Stokes number. The flow will be plotted as an inset within a plot of the efficiency')
parser.add_argument('-o', '--output', type=str, help='output plot to (e.g. pdf or png)')
parser.add_argument('-st', '--stokes', default=2., type=float, help='Stokes number (default=2.0)')
parser.add_argument('-c', '--crop', action='store_true', help='crop flow plot inside outer circle (of radius rout)')
parser.add_argument('--rout', default=3.0, type=float, help='outer radius if not automatically set by flow field (default=3.0)')
parser.add_argument('--tmax', default=10.0, type=float, help='max time (default=10.0)')
parser.add_argument('--y', default=None, help='launch values, e.g. 0.2,0.3')

parser.add_argument('-fc', '--flow_colour', default=palette.flow_colour, type=str, help='colour of background flow')
parser.add_argument('-sc', '--streamline_colour', default=palette.streamline_colour, type=str, help='colour of streamlines')
parser.add_argument('-xtc', '--critical_trajectory_colour', default=palette.critical_trajectory_colour, type=str, help='colour of critical trajectory')
parser.add_argument('-ctc', '--colliding_trajectory_colour', default=palette.colliding_trajectory_colour, type=str, help='colour of colliding trajectories')
parser.add_argument('-ntc', '--noncolliding_trajectory_colour', default=palette.noncolliding_trajectory_colour, type=str, help='colour of noncolliding trajectories')
parser.add_argument('-slw', '--streamline_line_width', default=0.5, type=str, help='thickness of streamlines')
parser.add_argument('-tlw', '--trajectory_line_width', default=0.5, type=str, help='thickness of particle trajectories')

parser.add_argument('-cc', '--cylinder_colour', default=palette.cylinder_colour, type=str, help='colour of cylindrical obstacle')

parser.add_argument('-cm', '--colliding_manifold', action='store_true', help='show colliding manifold')
parser.add_argument('-nm', '--noncolliding_manifold', action='store_true', help='show noncolliding manifold')
parser.add_argument('-cmc', '--colliding_manifold_colour', default=palette.colliding_manifold_colour, type=str, help='colour of colliding trajectories background')
parser.add_argument('-nmc', '--noncolliding_manifold_colour', default=palette.noncolliding_manifold_colour, type=str, help='colour of noncolliding trajectories background')

parser.add_argument('flowparams', nargs='*',
                    help='parameters to pass to flow field')
args = parser.parse_args()

args.flowparams = [eval(p) for p in args.flowparams]
flow_module = flow_fields[args.flow]
flow = flow_module.FlowField(*args.flowparams)

if args.flow == 'kuwabara':
    args.rout = flow.l

xmin = -args.rout
xmax = args.rout
ymin = -args.rout
ymax = args.rout

delta = 0.5
calc_xmin = xmin - delta
calc_xmax = xmax + delta
calc_ymin = ymin - delta
calc_ymax = ymax + delta

if args.efficiency:
    figsize = (3.375, 3)
    plt.figure(figsize=figsize)
    main = plt.gca()

    St, dSt, eff = np.genfromtxt(args.efficiency, skip_header=1).T
    Stc = St[0] - dSt[0]
    main.loglog(dSt, eff, 'o', mfc='w')
    main.set_xlabel('$\mathrm{St} - \mathrm{St}_c$')
    main.set_ylabel('capture efficiency $\lambda$')

    figsize = tuple(0.1*s for s in figsize)
    lin = main.inset_axes([0.125, 0.975 - figsize[1], figsize[0], figsize[1]])
    lin.plot(St, eff, 'o', mfc='w')
    lin.set_xlabel('$\mathrm{St}$', labelpad=-5)
    lin.set_ylabel('$\lambda$', labelpad=-5)

    lin.set_xlim([0., 2])
    lin.set_xticks([0, 2], minor=False)
    lin.set_xticks([1], minor=True)

    main_xlim = main.get_xlim()
    main_ylim = main.get_ylim()
    lin_xlim = lin.get_xlim()
    lin_ylim = lin.get_ylim()

    # Fit the straight line of gradient 1/2 to obtain the intercept, and plot this best fit.
    from scipy.optimize import curve_fit
    f = lambda x,c: c*x**0.5
    c, cerr = curve_fit(f, dSt[:5], eff[:5], p0=0.1)
    c = c[0]
    x = np.geomspace(1e-8, 1e8, 100)
    pl1, = main.plot(x, f(x, c), '-k', zorder=-1)
    pl2, = lin.plot(Stc + x, f(x, c), '-k', zorder=-1)

    main.set_xlim(main_xlim)
    main.set_ylim(main_ylim)
    lin.set_xlim(lin_xlim)
    lin.set_ylim(lin_ylim)

    # Label location of critical Stokes
    lin.annotate('$\mathrm{St}_c$',
                 xytext=(Stc, 0.75*lin_ylim[-1]),
                 xy=(Stc, 0),
                 ha='center', va='bottom',
                 arrowprops=dict(lw=0.5, arrowstyle='-', ec=args.critical_trajectory_colour, ls='--'),
                 zorder=-1)

    # A triangle to show gradient of straight line fit (in log space).
    x = 3e-3
    from mpltools import annotation
    annotation.slope_marker((x, f(x, 0.75*c)), (1, 2),
                            poly_kwargs=dict(ec='black', fill=False, lw=0.5),
                            text_kwargs=dict(fontsize=8))

    ax = main.inset_axes([0.55, 0.05, 0.45, 0.45])
    ax.set_aspect(1)

else:
    figsize = 3.375 / 1.5
    plt.figure(figsize=(figsize, figsize))
    ax = plt.gca()

#ax = plt.axes(projection='polar')
St = args.stokes

all_y0 = np.linspace(calc_ymin, calc_ymax, 21)
for y0 in all_y0:
    r0 = (xmin, y0)

    if abs(y0) < 1e-2:
        pl, = ax.plot(np.array([xmin, xmax]), np.array([y0, y0]),
                      lw=args.streamline_line_width, c=args.streamline_colour, zorder=-1)
    else:
        x, y = flow.streamline_cartesian(r0, xmax=calc_xmax, ymax=calc_ymax)

        if args.crop:
            r2 = x**2 + y**2
            select = r2 <= args.rout**2
            x = x[select]
            y = y[select]

        pl, = ax.plot(x, y, '-', lw=args.streamline_line_width, c=args.streamline_colour, zorder=-1)
        #pl, = ax.plot(x, y, '-', lw=streamline_lw, c=args.streamline_colour, zorder=-1)

    add_arrow(pl, x=-1.5, zorder=-1)
    add_arrow(pl, x=1.5, zorder=-1)

eff = flow.capture_efficiency(args.rout, St, max_step=0.1)

if args.y is None:
    thc = 0.5*eff
    yc = args.rout * np.sin(thc)
    yc -= 1e-6 # Ensure that we are slightly below transition so trajectory will actually terminate on surface.

    if args.colliding_manifold and args.efficiency:
        x = -1.15*args.rout
        ax.text(x, 0, r'$\lambda$', ha='right', va='center', fontsize=8)
        x = -1.1*args.rout
        y = 1.4*yc
        ax.annotate('', (x, y), (x, -y), annotation_clip=False,
                    weight=0.1, fontsize=8,
                    arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<->'))

    N = 11
    th0_arr = np.pi - np.linspace(0, eff, N)
    y0_arr = args.rout * np.sin(th0_arr)

    # Draw critical trajectory.
    y0 = yc
    x0 = -np.sqrt(args.rout**2 - y0**2)
    r0 = (x0, y0)
    t, (x, y, _, _), collides = flow.trajectory_cartesian(r0, St, max_step=0.1, return_collision=True)

    if args.crop:
        r2 = x**2 + y**2
        select = r2 <= args.rout**2
        x = x[select]
        y = y[select]

    ax.plot(x, y, '--', lw=args.trajectory_line_width, c=args.critical_trajectory_colour, zorder=1)
    ax.plot(x, -y, '--', lw=args.trajectory_line_width, c=args.critical_trajectory_colour, zorder=1)
    if args.colliding_manifold:
        ax.fill_between(x, -y, y, facecolor=args.colliding_manifold_colour, zorder=-9)

    # Ignore critical trajectory from main render
    mask = np.ones(y0_arr.size, dtype=bool)
    mask[N//2] = False
    y0_arr = y0_arr[mask]
else:
    y0_arr = np.array(eval('[' + args.y + ']'))

for y0 in y0_arr:
    x0 = -np.sqrt(args.rout**2 - y0**2)
    r0 = (x0, y0)
    t, (x, y, _, _), collides = flow.trajectory_cartesian(r0, St, max_step=0.1, return_collision=True)

    # t0 = np.arctan2(y0, -x0)
    # r0 = (args.rout, t0)
    # t, (x, y, _, _), collides = flow.trajectory(r0, St, max_step=0.1, return_collision=True)

    if args.crop:
        r2 = x**2 + y**2
        select = r2 <= args.rout**2
        x = x[select]
        y = y[select]

    c = args.colliding_trajectory_colour if collides else args.noncolliding_trajectory_colour
    ax.plot(x, y, lw=args.trajectory_line_width, c=c, zorder=1)
    ax.plot(x, -y, lw=args.trajectory_line_width, c=c, zorder=1)

ax.add_patch(Circle((0,0), args.rout, fc='None', ec='k', ls='dashed', lw=1, zorder=10))


eps = 5e-2
background_colour = args.flow_colour
if args.noncolliding_manifold: background_colour = args.noncolliding_manifold_colour
if args.crop:
    ax.add_patch(Circle((0,0), args.rout, color=background_colour, lw=0, zorder=-10))
else:
    ax.fill_between([xmin-eps, xmax+eps], ymin-eps, ymax+eps, facecolor=background_colour, zorder=-10)
ax.set_xlim([xmin-eps, xmax+eps])
ax.set_ylim([ymin-eps, ymax+eps])
ax.add_patch(Circle((0,0), 1, color=args.cylinder_colour, lw=0))
#ax.plot(1, 0, 'ko', mfc='w')

ax.axis('off')

if args.output:
    if args.output[-4:] == '.png' and not args.efficiency:
        plt.savefig(args.output, bbox_inches='tight', pad_inches=0, dpi=1000)
    else:
        plt.savefig(args.output)

plt.show()
