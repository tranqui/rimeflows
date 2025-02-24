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
from flows import kuwabara, chord, ellipse

from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

fig, axes = plt.subplots(nrows=2)
ax1, ax2 = axes


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

ax1.add_patch(Circle((0,0), 1, fc='r', alpha=0.75, ec='None', zorder=10))

ax1.set_xlim([-1.2*L, 1.2*L])
ax1.set_ylim([-0.5*L, 0.5*L])
ax1.set_aspect(1)


## Second plot: flow around an ellipse.

a, b = 0.3, 0.5
alpha = 0 #np.pi/2
flow = ellipse.FlowField(a, b, alpha)
ax = plt.gca()

L = 1.5
eps = 1e-6
x0 = -L
for y0 in np.linspace(-L/4, L/4, 4):
    x, y = flow.streamline([x0, y0], ymax=(L+eps))
    pl, = plt.plot(x, y, 'k-', lw=0.5)

    for x in [-0.5*L, 0.5*L]:
        add_arrow(pl, x=x, zorder=-1)

arrowprops = dict(facecolor='black', lw=0.5, arrowstyle='<->')
ax2.annotate('', [0, 0], [a, 0], arrowprops=arrowprops)
ax2.annotate('', [0, 0], [0, b], arrowprops=arrowprops)
offset = 0.04
ax2.text(0.5*a, -offset, '$\hat\ell_x$', ha='center', va='top')
ax2.text(-offset, 0.5*b, '$\hat\ell_y$', ha='right', va='center')

from matplotlib.patches import Ellipse
patch = Ellipse((0,0), 2*a, 2*b, ec='None', fc='r', alpha=0.75, angle=flow.matplotlib_draw_angle)
ax.add_patch(patch)

ax2.set_xlim([-L, L])
ax2.set_ylim([-0.5*L, 0.5*L])
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
