#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate#, fsolve
from mpltools import annotation

from glob import glob
from natsort import natsorted

def adjust_lightness(color, amount):
    """Source: https://stackoverflow.com/a/49601444"""
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

cmap = plt.get_cmap('turbo_r')
paths = natsorted(glob('data/efficiency_hiemenz*.csv'))

Stc = {}

figsize1 = (3.375, 2.1)
figsize2 = (3.375, 2.5)
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(figsize1[0], 2*figsize1[1]), constrained_layout=True)
plt.figure(figsize=figsize2)
ax3 = plt.gca()

St_cross= []
eff_cross = []

first = True
for i, path in enumerate(paths):

    Re_str = path.split('.csv')[0].split('Re=')[-1]
    Re = float(Re_str) if Re_str != 'inft' else np.inf
    if Re >= 1e4: Re_str = '$10^{:.0g}$'.format(np.log10(Re))
    else: Re_str = '{:.0f}'.format(Re)

    St, dSt, eff = np.genfromtxt(path, skip_header=1).T

    Stc[Re] = np.average(St - dSt)

    boundary_layer_thickness = 1 / np.sqrt(Re)
    f = interpolate.interp1d(St, eff)
    guess = np.argmin((eff - boundary_layer_thickness)**2)
    if dSt[guess] < 1e-2:
        St_cross += [dSt[guess]]
        eff_cross += [eff[guess]]
        #print('{:.1g} {:.4g}'.format(Re, boundary_layer_thickness), St_cross[-1], eff_cross[-1])

    c = cmap(i/(len(paths)-1))
     # default colour scheme is very bright which obscures some of the lighter colours  (esp. greens)
    c = adjust_lightness(c, 0.75)

    label = Re_str
    if first:
        #label = '$\mathrm{Re} = %s$' % label
        first = False

    for ax in [ax2, ax3]: ax.plot(dSt, eff, c=c, lw=0.75, label=label)

    St = np.concatenate([[Stc[Re]], St])
    eff = np.concatenate([[0], eff])
    ax1.plot(St, eff, c=c, lw=0.75, label=label)

St, dSt, eff = np.genfromtxt('data/efficiency_shm.csv', skip_header=1).T
for ax in [ax2, ax3]:
    ax.plot(dSt, eff, c='k', lw=0.75, label='$\infty$')
    ax.plot(St_cross, eff_cross, 'k--')
St = np.concatenate([[0.25], St])
eff = np.concatenate([[0], eff])
ax1.plot(St, eff, c='k', lw=0.75, label='$\infty$')

# for ax in [ax1, ax2]:
#     #ax.legend(loc='upper left', title='Re', title_fontsize=8, alignment='left', ncol=2, fontsize=8, borderpad=0)
#     ax.legend(loc='lower right', title='Re', title_fontsize=8, alignment='left', ncol=2, fontsize=6, borderpad=0)
#     ax.set_xlabel('$\mathrm{St}$')
#     #ax2.set_ylabel('squared efficiency $\lambda^2$')
#     ax.set_ylabel('efficiency $\lambda$')
#     #ax.set_xlim([0, 2])
#     #ax.set_ylim([0, 0.35])
#     ax.set_xscale('log')
#     ax.set_yscale('log')

for ax in [ax1, ax2, ax3]:
    ax.set_ylabel('efficiency $\lambda = 2y_c$')

ax1.set_xlabel('$\mathrm{St}$')
#ax2.set_ylabel('squared efficiency $\lambda^2$')

ax1.legend(loc='lower right', title='Re', title_fontsize=8, alignment='center', ncol=3, fontsize=8, borderpad=0, labelspacing=0., handlelength=1.)
ax1.set_xlim([0, 2])
ax1.set_ylim([0, 0.35])

for ax in [ax2, ax3]:
    ax.set_xlabel('$\mathrm{St} - \mathrm{St}_c$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-5, 0.5])
    ax.set_xlim([1e-5, 1])

ax3.legend(loc='upper left', title='Re', title_fontsize=8, alignment='center', ncol=3, fontsize=8, borderpad=0, labelspacing=0., handlelength=1.)

f = lambda x,c: c*x**0.5
for ax in [ax2, ax3]:
    x, c = 2e-5, 5e-3
    #x, c = 1e-4, 0.6
    annotation.slope_marker((x, f(x, c)), (1, 2), invert=False, ax=ax,
                            poly_kwargs=dict(ec='black', fill=False, lw=0.5),
                            text_kwargs=dict(fontsize=8))

    x, y = np.average(St_cross[-2:]), np.average(eff_cross[-2:])
    ax.annotate('$\lambda = \mathrm{Re}^{-1/2}$', size=8,
                 xytext=(1.5*x, 0.2*y), xy=(x, y),
                 arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='-'))

    x, y = 5e-2, 1e-2
    ax.annotate('inviscid\nscaling', size=8, va='center', ha='center',
                 xytext=(7*x, y), xy=(x, y),
                 arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='->'))

for ax, figsize, gap, inset_width in zip([ax2, ax3], [figsize1, figsize2],
                                         [0.05, 0.035], [0.3, 0.225]):
    #inset_width = 0.3
    aspect = figsize[1] / figsize[0]
    #inset = ax.inset_axes([gap*aspect, 1-gap - inset_width, inset_width, inset_width])
    inset = ax.inset_axes([1 - inset_width - gap*aspect, gap, inset_width, inset_width])

    Re,Stc1 = np.array([(k,v) for k,v in Stc.items()]).T
    inset.semilogx(Re, Stc1)
    inset.set_ylim([0, 1])
    inset.axhline(y=0.25, ls='dotted')

    xticks = np.geomspace(1, 1e8, 9)
    inset.set_xticks(xticks)
    labels = inset.get_xticklabels()
    labels[1:-1] = ['' for _ in labels[1:-1]]
    inset.set_xticklabels(labels, fontsize=8)

    yticks = np.linspace(0, 1, 5)
    inset.set_yticks(yticks)
    labels = ['%.1f' % y for y in yticks]
    labels[1:-1] = ['' for _ in labels[1:-1]]
    inset.set_yticklabels(labels, fontsize=8)

    inset.set_xlabel('$\mathrm{Re}$', labelpad=-3, fontsize=8)
    inset.set_ylabel('$\mathrm{St}_\mathrm{c}$', labelpad=-3, fontsize=8)

    inset.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    inset.xaxis.set_label_position('top')

x, y, fontsize, ha, va = 0.025, 0.825, 18, 'left', 'bottom'
for ax, letter in zip([ax1, ax2], 'ab'):
    label = ax.text(x, y, r'\textbf{%s}' % letter, transform=ax.transAxes, zorder=20,
                     fontsize=fontsize, horizontalalignment=ha, verticalalignment=va)
    label.set_in_layout(False)

plt.show()
