#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate#, fsolve
from mpltools import annotation
plt.style.use('figstyle.mplstyle')

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

lw = 0.75
cmap = plt.get_cmap('turbo_r')
paths_hiemenz = natsorted(glob('data/efficiency_hiemenz*.csv'))
paths_power = glob('data/efficiency_power*.csv')
indices = np.argsort([float(p.split('m=')[-1].split('.csv')[0]) for p in paths_power])
paths_power = [paths_power[i] for i in indices]
for p in paths_power: print(p)
# paths_power = natsorted(glob('data/efficiency_power*.csv'))
# for p in paths_power: print(p)

figsize = (3.375, 2.5)
fig1, (ax3, ax1) = plt.subplots(nrows=2, figsize=(figsize[0], 2*figsize[1]), constrained_layout=True)
fig2, (ax4, ax2) = plt.subplots(nrows=2, figsize=(figsize[0], 2*figsize[1]), constrained_layout=True)

# Show collection efficiencies for Hiemenz model (ax1):

St_cross= []
eff_cross = []
Stc_hiemenz = {}
for i, path in enumerate(paths_hiemenz):

    Re_str = path.split('.csv')[0].split('Re=')[-1]
    Re = float(Re_str) if Re_str != 'inft' else np.inf
    if Re >= 1e4: Re_str = '$10^{:.0g}$'.format(np.log10(Re))
    else: Re_str = '{:.0f}'.format(Re)

    St, dSt, eff = np.genfromtxt(path, skip_header=1).T

    Stc_hiemenz[Re] = np.average(St - dSt)
    if Re in [5, 50, 500, 2000, 5000]: continue

    boundary_layer_thickness = 1 / np.sqrt(Re)
    f = interpolate.interp1d(St, eff)
    guess = np.argmin((eff - boundary_layer_thickness)**2)
    if dSt[guess] < 1e-2:
        St_cross += [dSt[guess]]
        eff_cross += [eff[guess]]
        #print('{:.1g} {:.4g}'.format(Re, boundary_layer_thickness), St_cross[-1], eff_cross[-1])

    c = cmap(i/(len(paths_hiemenz)-1))
     # default colour scheme is very bright which obscures some of the lighter colours  (esp. greens)
    c = adjust_lightness(c, 0.75)

    ax1.plot(dSt, eff, c=c, lw=lw, label=Re_str)

# Limiting inviscid efficiencies for Hiemenz (ax1) and power law (ax2) models:

St, dSt, eff = np.genfromtxt('data/efficiency_shm.csv', skip_header=1).T
ax1.plot(dSt, eff, c='k', lw=lw, label='$\infty$')
cross_plot, = ax1.plot(St_cross, eff_cross, 'k--')
ax2.plot(dSt, eff, c='k', lw=lw, label='1.0')

# Show collection efficiencies for power law model (ax2):

Stc_power = {}
Stc_power[1] = 0.25
first = True
for i, path in enumerate(paths_power):
    m_str = path.split('.csv')[0].split('m=')[-1]
    m = float(m_str) if m_str != 'inft' else np.inf

    St, dSt, eff = np.genfromtxt(path, skip_header=1).T

    Stc_power[m] = np.average(St - dSt)

    c = cmap(i/(len(paths_power)-1))
    # default colour scheme is very bright which obscures some of the lighter colours  (esp. greens)
    c = adjust_lightness(c, 0.75)

    ax2.plot(dSt, eff, c=c, lw=lw, label=m_str)

# Formatting/legends of efficiencies plots:

for ax in [ax1, ax2]:
    ax.set_ylabel('efficiency $\lambda = 2y_c$')
    ax.set_xlabel('$\mathrm{St} - \mathrm{St}_c$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-5, 0.5])
    ax.set_xlim([1e-5, 1])

leg1 = ax1.legend(loc='upper left', title='Re', title_fontsize=8, alignment='center', ncol=3, fontsize=8, borderpad=0, labelspacing=0., handlelength=1.)
ax2.legend(loc='upper left', title='$m$', title_fontsize=8, alignment='center', ncol=3, fontsize=8, borderpad=0, labelspacing=0., handlelength=1.)

St, dSt, eff = np.genfromtxt('data/efficiency_kuwabara_alpha=0.15_2.csv', skip_header=1).T
Stc_kuwabara = St[0] - dSt[0]
kuwabara_plot, = ax1.plot(dSt, eff, '-.', lw=lw)
ax1.legend([kuwabara_plot, cross_plot],
           [r'$\mathrm{Re}=0$ (Kuwabara $\alpha=0.15$)', r'$\lambda = \mathrm{Re}^{-1/2}$'],
           loc='lower left', fontsize=6, borderpad=-0.25, labelspacing=0., handlelength=2.1)
ax1.add_artist(leg1)

# Annotate the square-root slopes of the small inertia regions, and the inviscid limits:

f = lambda x,c: c*x**0.5
for ax in [ax1, ax2]:
    x, c = 1e-4, 0.6
    annotation.slope_marker((x, f(x, c)), (1, 2), invert=True, ax=ax,
                            poly_kwargs=dict(ec='black', fill=False, lw=0.5),
                            text_kwargs=dict(fontsize=8))

    x, y = 5e-2, 1e-2
    ax.annotate('inviscid\nscaling', size=8, va='center', ha='center',
                 xytext=(7*x, y), xy=(x, y),
                 arrowprops=dict(facecolor='black', lw=0.5, arrowstyle='->'))

# Show the critical Stokes numbers for the Hiemenz flow (ax3) and the power law model (ax4)

Re,Stc1 = np.genfromtxt('data/Stc_hiemenz.csv').T
pl1, = ax3.semilogx(Re, Stc1, lw=lw)
Re,Stc1 = np.genfromtxt('data/Stc_hiemenz_x05.csv').T
pl2, = ax3.semilogx(Re, Stc1, '--', lw=lw)
Re,Stc1 = np.genfromtxt('data/Stc_hiemenz_x2.csv').T
pl3, = ax3.semilogx(Re, Stc1, ':', lw=lw)

ax3.set_ylim([0, 1])
ax3.set_xlim([1, 1e8])
ax3.axhline(y=0.25, ls='dotted')
ax3.text(1e2, 0.25, 'inviscid limit', ha='center', va='bottom', fontsize=8)
ax3.axhline(y=Stc_kuwabara, ls='-.')
ax3.text(1e4, Stc_kuwabara, r'Kuwabara $\alpha=0.15$ ($\mathrm{Re} = 0$)', ha='center', va='bottom', fontsize=8)
ax3.legend([pl2, pl1, pl3], ['$x(t=0) = 1/2$', '1', '2'], loc='best')

xticks = np.geomspace(1, 1e8, 9)
ax3.set_xticks(xticks)
labels = ax3.get_xticklabels()
labels[1::2] = ['' for _ in labels[1::2]]
ax3.set_xticklabels(labels)

m,Stc1 = np.concatenate([[(1,0)], np.genfromtxt('data/Stc_power_m.csv')]).T
pl1, = ax4.plot(m, Stc1, '-', lw=lw, zorder=10)
m,Stc1 = np.concatenate([[(1,0)], np.genfromtxt('data/Stc_power_m_x05.csv')]).T
pl2, = ax4.plot(m, Stc1, '--', lw=lw)
m,Stc1 = np.concatenate([[(1,0)], np.genfromtxt('data/Stc_power_m_x2.csv')]).T
pl3, = ax4.plot(m, Stc1, ':', lw=lw)

ax4.set_xlim([0, 2])
ax4.set_ylim([0, 0.6])
ax4.axhline(y=0.25, ls='dotted')
ax4.text(0.5, 0.25, 'inviscid limit', ha='center', va='bottom', fontsize=8)
ax4.legend([pl2, pl1, pl3], ['$x(t=0) = 1/2$', '1', '2'], loc='best')

# ax4.legend(loc='best')
# ax4.axhline(y=0.25, ls='-.', label='inviscid limit')

for ax in [ax3, ax4]: ax.set_ylabel('$\mathrm{St}_\mathrm{c}$')
ax3.set_xlabel('$\mathrm{Re}$')
ax4.set_xlabel('$m$')

# Letter each subpanel:

x, y, fontsize, ha, va = -0.175, 0.9, 18, 'left', 'bottom'
for axes in [[ax3, ax1], [ax4, ax2]]:
    for ax, letter in zip(axes, 'ab'):
        label = ax.text(x, y, r'\textbf{%s}' % letter, transform=ax.transAxes, zorder=20,
                        fontsize=fontsize, horizontalalignment=ha, verticalalignment=va)
        label.set_in_layout(False)

fig1.savefig('singular-onset-hiemenz.pdf')
fig1.savefig('singular-onset-hiemenz.png')
fig2.savefig('singular-onset-power.pdf')
fig2.savefig('singular-onset-power.png')
# plt.show()
