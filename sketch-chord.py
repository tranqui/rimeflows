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
from matplotlib.patches import Circle
from flows import kuwabara, chord

from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

fig, axes = plt.subplots(nrows=2)
ax1, ax2 = axes
# ax2 = plt.gca()
# axes = np.array([ax2])


## First plot: flow around a circle.

volume_fraction = 0.05
circle_flow = kuwabara.FlowField(volume_fraction)

L = circle_flow.l
x0 = -L
for y0 in np.linspace(-0.2*L, 0.2*L, 4):
    x, y = circle_flow.streamline_cartesian([x0, y0], xmax=1.2*L, ymax=1.2*L)
    select = np.abs(x) < L
    x, y = x[select], y[select]
    pl, = ax1.plot(x, y, 'k-', lw=0.5)

    for x in [-0.5*L, 0.5*L]:
        add_arrow(pl, x=x, zorder=-1)

ax1.add_patch(Circle((0,0), 1, fc='r', ec='None', zorder=10))

ax1.set_xlim([-1.2*L, 1.2*L])
ax1.set_ylim([-0.5*L, 0.5*L])
ax1.set_aspect(1)


## Second plot: flow around a chord.

chord_length = 1
alpha = np.pi/2 - 0.25*np.pi
chord_flow = chord.FlowField(chord_length, alpha)

# By default chord.FlowField places the chord along y-axis, and alpha
# measures the angle of the flow field. To get the flow going from
# left to right we need to rotate.
theta = alpha - np.pi/2
R = np.matrix([[np.cos(theta), -np.sin(theta)],
               [np.sin(theta),  np.cos(theta)]])

x_chord, y_chord = [0, 0], [-0.5*chord_length, 0.5*chord_length]
x_chord, y_chord = (R * np.array((x_chord, y_chord))).A

L = 1.5
eps = 1e-6
x0 = -L
for y0 in np.linspace(-0.5*chord_length, 0.5*chord_length, 4):
    x, y = chord_flow.streamline(([x0, y0] * R).tolist()[0], max_step=0.1, ymax=(3*L+eps))
    x, y = (R * np.array((x, y))).A
    pl, = ax2.plot(x, y, 'k-', lw=0.5)

    for x in [-0.5*L, 0.5*L]:
        add_arrow(pl, x=x, zorder=-1)

# Render the chord itself.
pl, = ax2.plot(x_chord, y_chord, 'r-', lw=1)
ax2.plot([0, x_chord[-1]], [0, 0], '--', c=pl.get_color())
# Annotate the angle \alpha.
th_arc = np.linspace(0, alpha, 1000)
r_arc = 0.25 * chord_length
ax2.plot(r_arc * np.cos(th_arc), r_arc * np.sin(th_arc), '-', c=pl.get_color())
r_arc *= 1.4
ax2.text(r_arc * np.cos(alpha/2), r_arc * np.sin(alpha/2), r'$\alpha$',
         fontsize=14, c=pl.get_color(), ha='center', va='center')

# Stc = 0.12583899
# x0, y0, St = -1, -0.078, Stc + 1e-3
# t, (x, y, _, _), collides = chord_flow.trajectory([x0, y0], St=St, max_step=1e-2, tmax=1e2, ymax=(L+eps), return_collision=True)
# pl, = ax2.plot(x, y, 'b-', lw=1)
# if collides: ax2.plot(x[-1], y[-1], 'o', c=pl.get_color())

# x0, y0, St = -1, 0.077, Stc + 1e-3
# t, (x, y, _, _), collides = chord_flow.trajectory([x0, y0], St=St, max_step=1e-2, tmax=1e2, ymax=(L+eps), return_collision=True)
# pl, = ax2.plot(x, y, 'b-', lw=1)
# if collides: ax2.plot(x[-1], y[-1], 'o', c=pl.get_color())

ax2.set_xlim([-L,L])
# ax2.set_ylim([-L,L])
ax2.set_ylim([-0.5*L,0.5*L])
ax2.set_aspect(1)


for ax in axes.ravel():
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

for ax, y, letter in zip(axes, [0.9, 0.45], 'ab'):
    label = ax.text(0.05, y, r'\textbf{%s}' % letter, transform=fig.transFigure, zorder=20,
                    fontsize=18, horizontalalignment='left', verticalalignment='bottom')
    label.set_in_layout(False)

plt.savefig('sketch-obstacles.pdf')
plt.show()
