#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpltools import annotation
from arrows import add_arrow
plt.style.use('figstyle.mplstyle')

from flows import kuwabara, rg, stokes, shm

palette = ['#4d0e5a', '#871c90', '#cc67c2', '#dabbe3', '#efda4d', '#ffffc0']

flow_color = palette[-1]
flow_streamline_color = palette[-2]
cylinder_color = palette[3]
noncolliding_particle_trajectory_color = palette[2]
colliding_particle_trajectory_color = palette[1]

ny = 21
streamline_lw = 0.5

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(2*3.375, 3.375), constrained_layout=True)
#                         gridspec_kw={'hspace': 0, 'wspace': 0})
#plt.figure(figsize=(3.375, 0.5*3.375))

# Axes showing different flow fields.

diagram_axes = axes[:,:2]
for ax in diagram_axes.ravel():
    #ax.axis('off')
    ax.set_xticks([], minor=[])
    ax.set_yticks([], minor=[])
    ax.set_aspect(1)

ax_kuwabara, ax_rg = axes[0,:2]
ax_stokes, ax_inviscid = axes[1,:2]

# Main axes showing efficiency scalings for each flow field.

main_axes = axes[:,-2:]
gs = axes[0,0].get_gridspec()
for ax in main_axes.ravel(): ax.remove()
axbig = fig.add_subplot(gs[:, -2:])
inset = axbig.inset_axes([0.175, 0.675, 0.3, 0.3])

for ls, path in zip(['-', '--', '-.', ':'],
                    ['data/efficiency_kuwabara_alpha=0.1.csv',
                     'data/efficiency_stokes.csv',
                     'data/efficiency_inviscid_gamma=1.2.csv',
                     'data/efficiency_shm.csv']):
    st, delta, eff = np.genfromtxt(path).T
    axbig.loglog(delta, eff, ls)
    inset.plot(st, eff, ls)

eps = 0.1
axbig.set_xlim([10**(-3-eps), 10**eps])
axbig.set_ylim([3e-3, 0.2])
inset.set_xlim([0, 1])
inset.set_ylim([0, 0.2])

x = 1e-2
f = lambda x: 0.1*x**0.5
annotation.slope_marker((x, f(x)), (1, 2),
                        poly_kwargs=dict(ec='black', fill=False, lw=0.5),
                        text_kwargs=dict(fontsize=8))
axbig.set_xlabel(r'$\mathrm{St} - \mathrm{St}_c$', labelpad=0)
axbig.set_ylabel(r'capture window $\lambda$', labelpad=0)
inset.set_xlabel(r'inertia $\mathrm{St}$', labelpad=-5)
inset.set_ylabel(r'$\lambda$', labelpad=5)

# Kuwabara subplot.

alpha = 0.175
kuwabara_flow = kuwabara.FlowField(alpha)

xmin = -3
xmax = 3
ymin = -3
ymax = 3

delta = 0.5
calc_xmin = xmin - delta
calc_xmax = xmax + delta
calc_ymin = ymin - delta
calc_ymax = ymax + delta

#x0 = xmin
x0 = -1/alpha**0.5

all_y0 = np.linspace(calc_ymin, calc_ymax, ny)
for y0 in all_y0:
    if abs(y0) < 1e-2:
        pl, = ax_kuwabara.plot(np.array([xmax, xmin]), np.array([y0, y0]),
                               lw=streamline_lw, c=flow_streamline_color, zorder=-1)
    else:
        x, y = kuwabara_flow.streamline_cartesian([x0, y0], xmax=calc_xmax, ymax=calc_ymax)
        pl, = ax_kuwabara.plot(-x, y, '-', lw=streamline_lw, c=flow_streamline_color, zorder=-1)

    add_arrow(pl, x=-1.5, zorder=-1)
    add_arrow(pl, x=1.5, zorder=-1)

ax_kuwabara.add_patch(Circle((0,0), 1/alpha**0.5, fc='None', ec='k', ls='dashed', lw=1, zorder=10))

ax_kuwabara.annotate('', (0, 0), (0, 1/alpha**0.5), weight=0.1,
                     arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<->'))
ax_kuwabara.annotate(r'$\approx\! \frac{1}{\sqrt{\alpha}}$', (0, 1 + 0.5*(1/alpha**0.5 - 1)), (1, 0),
                     textcoords='offset points', ha='left', va='center', fontsize=6)

# RG subplot.

Re = 1/(1/alpha**0.5 - 1)**2
rg_flow = rg.FlowField(Re)

for y0 in np.linspace(calc_ymin, calc_ymax, ny):
    if abs(y0) < 1e-2:
        pl, = ax_rg.plot(np.array([xmax, xmin]), np.array([y0, y0]),
                         lw=streamline_lw, c=flow_streamline_color, zorder=-1)
    else:
        x, y = rg_flow.streamline_cartesian([x0, y0], xmax=calc_xmax, ymax=calc_ymax)
        pl, = ax_rg.plot(-x, y, lw=streamline_lw, c=flow_streamline_color, zorder=-1)

    add_arrow(pl, x=-1.5, zorder=-1)
    add_arrow(pl, x=1.5, zorder=-1)

ax_rg.add_patch(Circle((0,0), 1 + 1/Re**0.5, fc='None', ec='k', ls='--', lw=1, zorder=10))

ax_rg.annotate('', (0, 1), (0, 1 + 1/Re**0.5), weight=0.1,
               arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<->'))
ax_rg.annotate(r'$\approx\! \frac{1}{\sqrt{\mathrm{Re}}}$', (0, 1 + 0.5/Re**0.5), (1, 0),
               textcoords='offset points', ha='left', va='center', fontsize=6)

# Both cylinder subplots

for ax in [ax_kuwabara, ax_rg]:
    ax.fill_between([xmin, xmax], ymin, ymax, facecolor=flow_color, zorder=-10)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.add_patch(Circle((0,0), 1, color=cylinder_color, lw=0))
    ax.plot(1, 0, 'ko', mfc='w')

# Stokes subplot.

L = 3
x0 = 1.75*L

stokes_flow = stokes.FlowField()

for y0 in np.linspace(-L**0.5, L**0.5, 11):
    y0 = np.sign(y0) * y0**2
    x, y = stokes_flow.streamline([x0, y0], max_step=0.1)
    pl, = ax_stokes.plot(x, y, lw=streamline_lw, c=flow_streamline_color, zorder=-1)

    if y0 == 0:
        add_arrow(pl, x=1.5, zorder=-1)
        add_arrow(pl, x=3, zorder=-1)
        add_arrow(pl, x=4.5, zorder=-1)

    for yy in [1.5, 2.5]:
        add_arrow(pl, y=-yy, zorder=-1)
        add_arrow(pl, y=yy, zorder=-1)

ax_stokes.fill_between([-2*L, 2*L], -L, L, fc=flow_color, zorder=-10)
ax_stokes.fill_between([x0 - 2*L, 0], -L, L, fc=cylinder_color)

ax_stokes.set_xlim([x0 -2*L, x0])
ax_stokes.set_ylim([-L, L])

# Inviscid subplot.

inviscid_flow = shm.FlowField()

for y0 in np.linspace(-L**0.5, L**0.5, 11):
    y0 = np.sign(y0) * y0**2
    x, y = inviscid_flow.streamline([x0, y0], max_step=0.1)
    pl, = ax_inviscid.plot(x, y, lw=streamline_lw, c=flow_streamline_color, zorder=-1)

    if y0 == 0:
        add_arrow(pl, x=0.7, zorder=-1)
        add_arrow(pl, x=1.7, zorder=-1)
        add_arrow(pl, x=3.7, zorder=-1)

    for yy in [1.5, 2.5]:
        add_arrow(pl, y=-yy, zorder=-1)
        add_arrow(pl, y=yy, zorder=-1)

St = 500
y0_list = np.linspace(-0.4*L, 0.4*L, 6)
for y0 in y0_list:
    r0 = [x0, y0, inviscid_flow.u(x0, y0), 0]
    t, (x, y, _, _), collides = inviscid_flow.trajectory(r0, St, max_step=0.1, return_collision=True)

    # Chop off starts of innermost trajectories so we can insert an annotation there.
    if abs(y0) < 0.5:
        select = x < x0 - 0.4*L
        x = x[select]
        y = y[select]

    c = colliding_particle_trajectory_color if collides else noncolliding_particle_trajectory_color
    pl, = ax_inviscid.plot(x, y, lw=streamline_lw, c=c, zorder=1)
    add_arrow(pl, x=2.5)

y = y0_list[-2]
x = x0 - 0.25*L
ax_inviscid.annotate('', (x, -y), (x, y), weight=0.1,
               arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<->'))
ax_inviscid.annotate(r'$\lambda$', (x, 0), (4, 0),
                     textcoords='offset points', ha='left', va='center', fontsize=8)

ax_inviscid.fill_between([-2*L, 2*L], -L, L, fc=flow_color, zorder=-10)
ax_inviscid.fill_between([x0 - 2*L, 0], -L, L, fc=cylinder_color)

ax_inviscid.set_xlim([x0 -2*L, x0])
ax_inviscid.set_ylim([-L, L])

# Both toy models.

for ax in [ax_stokes, ax_inviscid]:
    ax.plot(0, 0, 'ko', mfc='w')

    offset = 0.25 * (xmax - xmin)
    ax.annotate('$x$', (0, 0), (offset, 0), fontsize=6, weight=0.1,
                horizontalalignment='center', verticalalignment='center',
                arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<-'))
    ax.annotate('$y$', (0, 0), (0, offset), fontsize=6, weight=0.1,
                horizontalalignment='center', verticalalignment='center',
                arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<-'))

# y0 = all_y0[ny//2+1]
# colors = list(reversed(palette[:3])) + ['k']
# for St,c in zip([0.2, 0.8, 1.5], colors):
#     t, (x, y, _, _) = flow.trajectory_cartesian([x0, y0], St)
#     t, (x2, y2) = flow.trajectory_cartesian([x0, y0], St=0, tmax=-1e2)
#     x = np.concatenate((np.flipud(x2), x))
#     y = np.concatenate((np.flipud(y2), y))

#     pl, = axk.plot(x, y, '-', lw=streamline_lw, c=c)
#     pl2, = ax_kuwabara.plot(x, y, '-', lw=streamline_lw, c=c)
#     r_final = np.sqrt(x[-1]**2 + y[-1]**2)
#     if (r_final - 1) < 1e-2:
#         axk.plot(x[-1], y[-1], 'o', c=c, mfc=c)
#         ax_kuwabara.plot(x[-1], y[-1], 'o', c=c, mfc=c)

#     add_arrow(pl, y=0.3)
#     add_arrow(pl, y=0.8)
#     if St > 1: add_arrow(pl2, x=-2)
#     add_arrow(pl2, x=2)

bbox = dict(pad=0.1, fc='None', ec='None')
#bbox = dict(boxstyle='round', pad=0.1, fc='white', ec='none', alpha=0.75)
fontsize = 8
ha, va = 'center', 'top'
x, y = 0.5, 0.99
ax_kuwabara.text(x, y, r'cylindrical grain',
                 transform=ax_kuwabara.transAxes, fontsize=fontsize,
                 bbox=bbox, horizontalalignment=ha, verticalalignment=va)
ax_rg.text(x, y, r'free cylinder',
           transform=ax_rg.transAxes, fontsize=fontsize,
           bbox=bbox, horizontalalignment=ha, verticalalignment=va)

#x, y = 0.98, 0.975
#ha, va = 'right', 'top'
x, y = 0.9, 1.025
ha, va = 'right', 'bottom'
label = ax_stokes.text(x, y, r'$\mathrm{Re} \to 0^+$',
                       transform=ax_stokes.transAxes, fontsize=fontsize,
                       bbox=bbox, horizontalalignment=ha, verticalalignment=va)
label.set_in_layout(False)
label = ax_inviscid.text(x, y, r'$\mathrm{Re} \to \infty$',
                         transform=ax_inviscid.transAxes, fontsize=fontsize,
                         bbox=bbox, horizontalalignment=ha, verticalalignment=va)
label.set_in_layout(False)

x, y = 0.5, 0.025
x, y = 0.5, -0.025
ha, va = 'center', 'top'
label = ax_stokes.text(x, y, r'$\mathbf{u} \propto (-x^2, 2 x y)^\top$ as $\epsilon \to 0^+$',
                       transform=ax_stokes.transAxes, fontsize=fontsize,
                       bbox=bbox, horizontalalignment=ha, verticalalignment=va)
label.set_in_layout(False)
label = ax_inviscid.text(x, y, r'$\mathbf{u} \propto (-x, y)^\top$ as $\epsilon \to 0^+$',
                         transform=ax_inviscid.transAxes, fontsize=fontsize,
                         bbox=bbox, horizontalalignment=ha, verticalalignment=va)
label.set_in_layout(False)

# Give each subplot a letter.

x, y, fontsize, ha, va = 0.025, 0.825, 18, 'left', 'bottom'
for letter,ax in zip('abcde', np.concatenate([diagram_axes.ravel(), [axbig]])):
    label = ax.text(x, y, r'\textbf{%s}' % letter, transform=ax.transAxes, zorder=20,
                    fontsize=fontsize, horizontalalignment=ha, verticalalignment=va)
    label.set_in_layout(False)

    if letter == 'b':
        x, y, ha, va = 0.975, 0.01, 'right', 'bottom'

# Smaller letters on main axes showing which efficiency lines belong to which flow field.
fontsize, ha, va = 14, 'center', 'center'
for letter, (x, y) in zip('acdb', [[3e-3, 5.25e-3], [2e-3, 1.1e-2], [7e-2, 7e-3], [4e-1, 7e-3]]):
    label = axbig.text(x, y, r'\textbf{%s}' % letter, zorder=20,
                       fontsize=fontsize, horizontalalignment=ha, verticalalignment=va)
    label.set_in_layout(False)

# Draw fancy lines between toy model subplot and region around stagnation point in
# cylindrical flow plot, so it looks like the toy model is the zoomed in region there.

fig.canvas.draw() # Need to do a preliminary draw to position the axes with constrained layout.
transFigure = fig.transFigure.inverted()

xleft, xright = 0.7, 1.5
yleft, yright = -0.4, 0.4
box_coords = [[xleft, yleft], [xright, yleft], [xright, yright], [xleft, yright], [xleft, yleft]]
box_coords = [transFigure.transform(ax_kuwabara.transData.transform(coord)) for coord in box_coords]

lw = 0.5
lines = []
for coord1, coord2 in zip(box_coords[:-1], box_coords[1:]):
    lines += [plt.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]), transform=fig.transFigure, lw=lw, ls='--')]

coord1 = transFigure.transform(ax_stokes.transAxes.transform([0, 1]))
coord2 = box_coords[0]
lines += [plt.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]), transform=fig.transFigure, lw=lw)]

coord1 = transFigure.transform(ax_stokes.transAxes.transform([1, 1]))
coord2 = box_coords[1]
lines += [plt.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]), transform=fig.transFigure, lw=lw)]

box_coords = [[xleft, yleft], [xright, yleft], [xright, yright], [xleft, yright], [xleft, yleft]]
box_coords = [transFigure.transform(ax_rg.transData.transform(coord)) for coord in box_coords]
for coord1, coord2 in zip(box_coords[:-1], box_coords[1:]):
    lines += [plt.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]), transform=fig.transFigure, lw=lw, ls='--')]

coord1 = box_coords[0]
coord2 = 0.5 * ( transFigure.transform(ax_kuwabara.transAxes.transform([1, 0])) +
                 transFigure.transform(ax_stokes.transAxes.transform([1, 1])) )


# Arrow connecting stagnation point in inviscid plot to stokes flow plot.

offset1 = np.array((-0.011, -0.005)) # annotate annoyingly isn't quite aligned with fig.transFigure for some reason...
offset2 = np.array((-0.015, -0.01))
arrow = plt.annotate('', coord1 + offset1, coord2 + offset2,
                     xycoords='figure fraction', textcoords='figure fraction',
                     arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<-'))
arrow.set_in_layout(False)

# Arrow connecting inviscid flow plot to shm plot (showing limit Re -> inf).

x = 0.9
coord1 = transFigure.transform(ax_rg.transAxes.transform([x, 0.5]))
coord2 = transFigure.transform(ax_inviscid.transAxes.transform([x, 0.975]))
arrow = plt.annotate('', coord1, coord2,
                     xycoords='figure fraction', textcoords='figure fraction',
                     arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='<-'))
arrow.set_in_layout(False)

# Line connecting efficiency definition on inviscid flow plot to y-axis of big plot.

coord1 = transFigure.transform(ax_inviscid.transAxes.transform([0.98, 0.5]))
coord2 = transFigure.transform(axbig.transAxes.transform([-0.175, 0.3]))
lines += [plt.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]), transform=fig.transFigure, lw=lw, ls='-')]

fig.lines.extend(lines)

plt.savefig('flow-fields.pdf')
plt.savefig('flow-fields.png')
