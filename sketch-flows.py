#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpltools import annotation
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from palettes import PurpleGold as palette
obstacle_colour = palette.cylinder_colour
# obstacle_colour = palette.noncolliding_manifold_colour


from flows import kuwabara, rg, stokes, shm, hiemenz, wedge

figsize = (3.375, 1.5*3.375)
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=figsize, constrained_layout=True)


# Left column showing various wedge flows:

alpha_crit = wedge.stokes.critical_alpha()

for alpha, ax in zip([0.9, alpha_crit, 0.2], axes[:,0]):
    flow = wedge.stokes.FlowField(alpha)
    # print(alpha, flow.m)

    # Draw streamlines.
    L = 0.25
    offset = (alpha - 0.2)**2 * L
    x0 = -L - offset
    for y0 in np.linspace(0, L**0.5, 7)**2:
        x, y = flow.streamline([x0, y0], max_step=0.1)
        pl, = ax.plot(x, y, '-', c=palette.streamline_colour, lw=0.5)
        ax.plot(x, -y, '-', c=pl.get_color(), lw=0.5)

    # Draw inertial trajectories that collide with the wedge.
    for y0 in np.linspace(0, (0.25*L)**0.5, 5)[1:]**2:
        St = 1
        t, (x, y, _, _), collides = flow.trajectory([x0, y0], St, max_step=0.1, return_collision=True)
        # assert collides
        c = palette.colliding_trajectory_colour if collides else palette.noncolliding_trajectory_colour
        pl, = ax.plot(x, y, '-', c=c, lw=0.5)
        if collides:
            ax.plot(x[-1], y[-1], 'o', c=c)

    # Draw obstacle.
    x = np.linspace(0, 10*L, 100)
    ax.fill_between(x, -flow.y_wedge(x), flow.y_wedge(x), color=obstacle_colour, ec='None')

    ax.set_xlim([x0, x0 + 2*L])
    ax.set_ylim([-L, L])


# Right column varying Reynolds number in Hiemenz flow:

for Re, ax in zip([1e-4, 1e8], axes[:,1]):
    # print(Re)
    flow = hiemenz.FlowField(Re)

    L = 0.5
    x0 = 1.5*L
    for y0 in np.linspace(0, L**0.5, 11)**2:
        x, y = flow.streamline([x0, y0], tmax=1e6)
        pl, = ax.plot(-x, y, '-', c=palette.streamline_colour, lw=0.5)
        ax.plot(-x, -y, '-', c=pl.get_color(), lw=0.5)

    # Draw inertial trajectories that collide with the wedge.
    # for y0 in np.linspace(0, (0.05*L)**0.5, 5)[1:]**2:
    for y0 in np.linspace(0, (0.25*L)**0.5, 5)[1:]**2:
        St = 1e4
        t, (x, y, _, _), collides = flow.trajectory([x0, y0], St, tmax=1e6, return_collision=True)
        c = palette.colliding_trajectory_colour #if collides else palette.noncolliding_trajectory_colour
        pl, = ax.plot(-x, y, '-', c=c, lw=0.5)
        if collides:
            ax.plot(-x[-1], y[-1], 'o', c=c, zorder=100)

    ax.fill_between([0, 0.5], -L, L, fc=obstacle_colour, zorder=10)

    ax.set_xlim([-x0, -x0 + 2*L])
    ax.set_ylim([-L, L])

ax_empty = axes[-1,1]
for spine in ax_empty.spines.values():
    spine.set_visible(False)
ax_empty.text(0.5, 0.5, 'no equivalent\nplanar flow',
              ha='center', va='center', transform=ax_empty.transAxes)

for ax in axes.ravel():
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.set_aspect('equal')
    # ax.set_xlabel('$x$')
    # ax.set_ylabel('$y$')

for ax, letter in zip(axes.T.ravel(), 'abcde'):
    label = ax.text(2e-2, 0, r'\textbf{%s}' % letter, transform=ax.transAxes, zorder=100,
                    fontsize=18, ha='left', verticalalignment='bottom')
    label.set_in_layout(False)

for ax, text1, text2 in zip(axes[:,0],
                            [r'blunt wedge', r'critical wedge', r'sharp wedge'],
                            [r'$m > 1$', r'$m = 1$', r'$m < 1$']):
    label = ax.text(0.15, 0.06, text1, transform=ax.transAxes, zorder=100,
                    fontsize=10, ha='left', verticalalignment='center')
    label.set_in_layout(False)
    label = ax.text(0.15, 0.94, text2, transform=ax.transAxes, zorder=100,
                    fontsize=10, ha='left', verticalalignment='center')
    label.set_in_layout(False)

for ax, text in zip(axes[:-1,1],
                      [r'$\mathrm{Re} \ll 1$',
                       r'$\mathrm{Re} \gg 1$']):
    label = ax.text(0.15, 0.06, text, transform=ax.transAxes, zorder=100,
                    fontsize=10, ha='left', verticalalignment='center')
    label.set_in_layout(False)

plt.savefig('wedge-flows.pdf')
plt.savefig('wedge-flows.png')
plt.show()
