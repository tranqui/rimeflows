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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from flows import stokes, shm
stokes = stokes.OnAxis()
shm = shm.OnAxis()

figsize = 3.375 # (inches)
aspect_ratio = 2
fig = plt.figure(figsize=(figsize, figsize/aspect_ratio))

ax = plt.gca()
ax.axis('off')

palette = ['#4d0e5a', '#871c90', '#cc67c2', '#dabbe3', '#efda4d', '#ffffc0']

stable_manifold_color = palette[-1]
unstable_manifold_color = palette[3]

stable_streamline_color = palette[-2]
unstable_streamline_color = palette[2]

separatrix_color = 'k'
nullcline_color = separatrix_color

#bbox = dict(pad=0.1, fc='None', ec='None')
bbox = dict(boxstyle='round', pad=0.1, fc='white', ec='none', alpha=0.75)

N = 20000
xmin = -0.5
xmax = 1.5
vmin = -1
vmax = 1
m = 3

v0 = vmax

#for x0 in np.linspace(m*xmin, m*xmax, 40):
delta = -0.105
for x0 in np.linspace(m*xmin + delta, m*xmax + delta, 40):
    x, v = stokes.streamline(x0, v0, tmax=1e3, N=N,
                             xmin=(m*xmin), xmax=(m*xmax), vmin=(m*vmin), vmax=(m*vmax))

    if abs(x[-1]) < 1e-2: c = stable_streamline_color
    else: c = unstable_streamline_color

    pl, = ax.plot(x, v, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=-0.075, zorder=-1)
    add_arrow(pl, y=-0.4, zorder=-1)

x, v = stokes.limit_unstable_streamline()
ax.plot(x, v, lw=0.5, c=unstable_streamline_color, zorder=-8)

x, v = stokes.nullcline()
ax.plot(x, v, '--', c=nullcline_color)

x, v, (xc, vc) = stokes.separatrix(return_critical_point=True)
select = np.logical_and(x < 0, v < 1)
pl, = ax.plot(x[select]-1e-2, v[select], c=unstable_streamline_color, lw=0.5)
add_arrow(pl, y=0.1, zorder=-1, direction='backward')
#ax.plot(x, v, c=separatrix_color)
#ax.plot(0, 0, 'o', c=separatrix_color, mfc='w', zorder=20)
#ax.plot(xc, vc, 'o', c=separatrix_color, mfc=unstable_streamline_color, zorder=20)
#ax.plot(0, 0, 'o', c=separatrix_color, mfc=separatrix_color)

ax.fill_between([xmin, xmax], vmin, vmax, facecolor=unstable_manifold_color, zorder=-10)
ax.fill_between(x[x > 0], 0, v[x > 0], facecolor=stable_manifold_color, zorder=-5)
ax.fill_between(x[x < 0.1], 0, v[x < 0.1], facecolor=stable_manifold_color, zorder=-5)

#ax.set_xlim([xmin, xmax])
#ax.set_ylim([vmin, vmax])

ax.set_xlim([0, 1])
ax.set_ylim([-0.5, 0])

plt.savefig('data/backdrop.png')
