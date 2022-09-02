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

import argparse
parser = argparse.ArgumentParser(description='show streamlines for hybrid field')
parser.add_argument('epsilon', type=float, nargs='?', default=0,
                    help='epsilon parameter for hybrid flow field (default=0 for Stokes flow field)')
args = parser.parse_args()

from flows import hybrid
flow = hybrid.OnAxis(args.epsilon)

figsize = 3.375 # (inches)
aspect_ratio = 2
fig = plt.figure(figsize=(figsize, figsize/aspect_ratio))

ax = plt.gca()

palette = ['#4d0e5a', '#871c90', '#cc67c2', '#dabbe3', '#efda4d', '#ffffc0']

stable_manifold_color = palette[-1]
unstable_manifold_color = palette[3]

stable_streamline_color = palette[-2]
unstable_streamline_color = palette[2]

separatrix_color = 'k'
nullcline_color = separatrix_color

#bbox = dict(pad=0.1, fc='None', ec='None')
bbox = dict(boxstyle='round', pad=0.1, fc='white', ec='none', alpha=0.75)

xmin = -0.5
xmax = 1.5
vmin = -1
vmax = 1
m = 3

v0 = vmax

def draw_streamline(x0, tmax=1e3, N=50000):
    x, v = flow.streamline(x0, v0, tmax=1e3, N=N,
                           xmin=(m*xmin), xmax=(m*xmax), vmin=(m*vmin), vmax=(m*vmax))

    if x0 < separatrix_x and x[-1] > 0: c = stable_streamline_color
    else: c = unstable_streamline_color

    pl, = ax.plot(x, v, lw=0.5, c=c, zorder=-2)
    arrows = []
    arrows += add_arrow(pl, y=-0.075, zorder=-1)
    arrows += add_arrow(pl, y=-0.4, zorder=-1)
    return pl, arrows

# Find the location of the separatrix at our initial v:
x, v = flow.separatrix()
guess = np.where((v-v0)**2 < 0.1)[0]
guess = guess[np.argmax(x[guess])]
ax.plot(x[guess], v[guess], 'ko', mfc='None', mew=0.5, zorder=100)
fx = interpolate.interp1d(np.arange(len(x)), x, kind='quadratic')
fv = interpolate.interp1d(np.arange(len(v)), v, kind='quadratic')
separatrix_x = fx( optimize.minimize(lambda i: (fv(i) - v0)**2, guess).x )[0]

# Sample points on either side of the separatrix across the full visible range so we
# fill space with stream lines
eps = 2e-3
dx = 0.15
delta_x = np.arange(0, m*(xmax - xmin), dx)[1:]
central_x = separatrix_x + eps
x_range = np.hstack((central_x - delta_x, central_x + delta_x))
for x0 in x_range: draw_streamline(x0)
sep_pl, sep_arrows = draw_streamline(central_x)

x, v = flow.limit_unstable_streamline()
ax.plot(x, v, lw=0.5, c=unstable_streamline_color, zorder=-8)

x, v = flow.nullcline()
ax.plot(x, v, '--', c=nullcline_color)

x, v, (xc, vc) = flow.separatrix(return_critical_point=True)
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

ax.set_xlim([m*xmin, m*xmax])
ax.set_ylim([m*vmin, m*vmax])

ax.set_xlim([0, 1])
ax.set_ylim([-0.5, 0])

ax.set_xlabel(r'$\mathrm{St}^{\frac{1}{1 - \epsilon}} \, x$')
ax.set_ylabel(r'$\mathrm{St}^{\frac{2 - \epsilon}{1 - \epsilon}} \, \dot{x}$')

label = ax.text(0.99, 0.06, r'$\epsilon={:.3f}$'.format(args.epsilon), transform=ax.transAxes,
                fontsize=10, bbox=bbox, horizontalalignment='right', verticalalignment='bottom')

plt.savefig('data/hybrid_streamlines_eps={:.3f}.png'.format(args.epsilon), dpi=1000)

ax.axis('off')

label.remove()

sep_pl.set_visible(False)
for arrow in sep_arrows: arrow.set_visible(False)

plt.savefig('data/hybrid_streamlines_noaxis_eps={:.3f}.png'.format(args.epsilon), bbox_inches='tight', pad_inches=0, dpi=1000)
#plt.show()
