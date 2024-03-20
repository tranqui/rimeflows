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
from flows import chord

from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

# fig, axes = plt.subplots(nrows=2)
# ax1, ax2 = axes
ax2 = plt.gca()
axes = np.array([ax2])

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
for y0 in np.linspace(-3*L, 3*L, 21):
    x, y = chord_flow.streamline([x0, y0], max_step=0.1, ymax=(3*L+eps))
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
r_arc *= 1.3
ax2.text(r_arc * np.cos(alpha/2), r_arc * np.sin(alpha/2), r'$\alpha$',
         fontsize=10, c=pl.get_color(), ha='center', va='center')

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
ax2.set_ylim([-L,L])
ax2.set_aspect(1)

for ax in axes.ravel():
    ax.set_xticks([])
    ax.set_yticks([])

plt.show()
