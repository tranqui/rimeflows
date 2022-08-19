#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from flows import stokes, shm
stokes = stokes.OnAxis()
shm = shm.OnAxis()

figsize = 3.375 # (inches)
fig, (ax1, ax3) = plt.subplots(nrows=2, figsize=(figsize, 1.5*figsize))

div = make_axes_locatable(ax3)
ax2 = div.append_axes("top", size="100%", pad=0.05)
axes = [ax1, ax2, ax3]

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

    pl, = ax1.plot(x, v, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=0.8, zorder=-1)
    add_arrow(pl, y=0.1, zorder=-1)
    #add_arrow(pl, y=-0.41, zorder=-10)
    add_arrow(pl, y=-0.81, zorder=-1)

x, v = stokes.limit_unstable_streamline()
ax1.plot(x, v, lw=0.5, c=unstable_streamline_color, zorder=-8)

x, v = stokes.nullcline()
ax1.plot(x, v, '--', c=nullcline_color)

x, v, (xc, vc) = stokes.separatrix(return_critical_point=True)
select = np.logical_and(x < 0, v < 1)
pl, = ax1.plot(x[select]-1e-2, v[select], c=unstable_streamline_color, lw=0.5)
add_arrow(pl, y=0.1, zorder=-1, direction='backward')
#ax1.plot(x, v, c=separatrix_color)
ax1.plot(0, 0, 'o', c=separatrix_color, mfc='w', zorder=20)
ax1.plot(xc, vc, 'o', c=separatrix_color, mfc=unstable_streamline_color, zorder=20)
#ax1.plot(0, 0, 'o', c=separatrix_color, mfc=separatrix_color)

# xs = np.linspace(0, 1/8, 100)
# vs = (1 - np.sqrt(1 - 8*xs) - 12*xs + 8*xs*np.sqrt(1 - 8*xs)) / 24
# ax1.plot(xs, vs, 'b-.', lw=1, label='slow manifold', zorder=5)


ax1.fill_between([xmin, xmax], vmin, vmax, facecolor=unstable_manifold_color, zorder=-10)
ax1.fill_between(x[x > 0], 0, v[x > 0], facecolor=stable_manifold_color, zorder=-5)
ax1.fill_between(x[x < 0.1], 0, v[x < 0.1], facecolor=stable_manifold_color, zorder=-5)

ax1.set_xlim([xmin, xmax])
ax1.set_ylim([vmin, vmax])

ax1.text( 0.325,  0.5, 'noncolliding\ntrajectories', fontsize=8,
          bbox=bbox, horizontalalignment='center', verticalalignment='center')
ax1.text( 1.15, -0.3, 'colliding\ntrajectories', fontsize=8,
          bbox=bbox, horizontalalignment='center', verticalalignment='center')
#ax1.text(-0.06, -0.45, 'collision at\nfinite $t$', fontsize=8, rotation=90,
#         horizontalalignment='center', verticalalignment='center')
ax1.text(0., -0.5, 'collision at\nfinite $t$', fontsize=8, rotation=90,
         horizontalalignment='center', verticalalignment='center')

ax1.set_xlabel('$\mathrm{St} \, x$')
ax1.set_ylabel('$\mathrm{St}^2 \, \dot{x}$')

eps = 1e-6

xmin = -0.5
xmax = 1.5
vmin = -4
vmax = 1
m = 3

St = 0.2
v0 = 0

c1 = -(1 + np.sqrt(1 - 4*St)) / (2*St)
c2 = -(1 - np.sqrt(1 - 4*St)) / (2*St)

for x0 in np.linspace(xmin, xmax, 11):
    if np.abs(x0) < eps: continue

    for v0 in [vmin, vmax]:
        stable = v0 > 0 or x0 > v0/c1
        if stable: c = stable_streamline_color
        else: c = unstable_streamline_color

        x, v = shm.streamline(x0, v0, St=St, N=N,
                              xmin=(m*xmin), xmax=(m*xmax), vmin=(m*vmin), vmax=(m*vmax))
        pl, = ax2.plot(x, v, lw=0.5, c=c, zorder=-2)
        if v0 > 0: add_arrow(pl, y=0.5, zorder=-1)
        else:
            add_arrow(pl, y=-1.5, zorder=-1)
            add_arrow(pl, y=-3, zorder=-1)


x0 = xmax
for v0 in np.linspace(-3.25, -1.5, 4):
    stable = v0 > 0 or x0 > v0/c1
    if stable: c = stable_streamline_color
    else: c = unstable_streamline_color

    x, v = shm.streamline(x0, v0, St=St, N=N,
                          xmin=(m*xmin), xmax=(m*xmax), vmin=(m*vmin), vmax=(m*vmax))
    pl, = ax2.plot(x, v, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=-3, zorder=-1)
    if v0 < c2*x0:
        add_arrow(pl, y=-1.5, zorder=-1)

#ax2.fill_between([xmin, xmax], vmin, vmax, facecolor=unstable_manifold_color, zorder=-10)
#ax2.fill_between([, R], -R, R, facecolor=stable_manifold_color, zorder=-5)
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([vmin, vmax])

c1 = -(1 + np.sqrt(1 - 4*St)) / (2*St)
c2 = -(1 - np.sqrt(1 - 4*St)) / (2*St)
x1 = np.linspace(xmin, 0, 10)
x2 = np.linspace(0, xmax, 10)
# Relevant separatrices
#ax2.plot(x2, c1*x2, '-', c=separatrix_color)
#ax2.plot(x1, c2*x1, '-', c=separatrix_color)
# Irrelevant separatrices
#ax2.plot(x1, c1*x1, '-', c=separatrix_color)
#ax2.plot(x2, c2*x2, '-', c=separatrix_color)
ax2.fill_between(x1, vmin, c2*x1, facecolor=unstable_manifold_color, zorder=-10)
ax2.fill_between(x1, vmax, c2*x1, facecolor=stable_manifold_color, zorder=-10)
ax2.fill_between(x2, vmin, c1*x2, facecolor=unstable_manifold_color, zorder=-10)
ax2.fill_between(x2, vmax, c1*x2, facecolor=stable_manifold_color, zorder=-10)

R = 2
th0 = np.pi/3
St = 0.5
for RR in np.linspace(-R, R, 11):
    if np.abs(RR) < eps: continue
    x0, v0 = RR*np.cos(th0), RR*np.sin(th0)

    c = unstable_streamline_color
    #x, v = shm.streamline(x0, v0, St=0.26, xmin=(-R-eps), xmax=(R+eps), vmin=(-R-eps), vmax=(R+eps))
    x, v = shm.streamline(x0, v0, St=St, N=N,
                          xmin=(m*xmin), xmax=(m*xmax), vmin=(m*vmin), vmax=(m*vmax))
    pl, = ax3.plot(x, v, lw=0.5, c=c, zorder=-2)
    add_arrow(pl, y=0.5, zorder=-1)
    add_arrow(pl, y=-0.5, zorder=-1)
    add_arrow(pl, y=-2.5, zorder=-1)

ax3.fill_between([xmin, xmax], vmin, vmax, facecolor=unstable_manifold_color, zorder=-5)
# ax3.set_xlim([-R/2, R/2])
# ax3.set_ylim([-R/2, R/2])
ax3.set_xlim([xmin, xmax])
ax3.set_ylim([vmin, vmax])

for ax in [ax2, ax3]:
    ax.plot(0, 0, 'o', mfc='w', c=separatrix_color, zorder=20)

    ax.set_ylabel('$\dot{x}$')

    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, -x, '--', c=nullcline_color)

for ax in axes:
    ax.fill_between([-1e3, 0], -1e3, 1e3, facecolor='black', alpha=0.125, zorder=10)
    ax.axvline(x=0, lw=0.5)

ax3.set_xlabel('$x$')
ax2.set_xticklabels([])

ax1.annotate(r'$x \sim 1/t$', (0, 0), (0.25, 0.1), #weight=0.1,
             fontsize=10, horizontalalignment='center', verticalalignment='bottom',
             bbox=bbox, arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='-'))
ax2.annotate(r'$x \sim e^{-t}$', (0, 0), (0.5, 0.75), #weight=0.1,
             fontsize=10, horizontalalignment='center', verticalalignment='top',
             bbox=bbox, arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='-'))

plt.text(0.9, 0.12, r'$\mathrm{St} < \frac{1}{4}$', transform=ax2.transAxes,
         fontsize=10, bbox=bbox, horizontalalignment='center', verticalalignment='center')
plt.text(0.9, 0.12, r'$\mathrm{St} > \frac{1}{4}$', transform=ax3.transAxes,
         fontsize=10, bbox=bbox, horizontalalignment='center', verticalalignment='center')

for ax,l in zip(axes, 'abc'):
    label = plt.text(-0.2, 1, (r'\textbf{%s}' % l), transform=ax.transAxes,
                     fontsize=18, horizontalalignment='left', verticalalignment='top')
    label.set_in_layout(False)

plt.savefig('phase-portraits.pdf')
plt.savefig('phase-portraits.png')
