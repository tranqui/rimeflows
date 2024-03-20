#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt

x0, Stc = np.genfromtxt('data/Stc_chord_x0.csv').T
plt.plot(x0, Stc)
plt.xlabel('$x_0$')
plt.ylabel('$\mathrm{St}_\mathrm{c}$')

plt.figure()

first = True
for x0 in np.linspace(0.1, 1, 10):
    St, dSt, eff = np.genfromtxt(f'data/efficiency_chord_x0=-{x0:.1f}.csv').T
    label = f'{x0:.1f}'
    if first:
        label = f'$x_0 = {label}$'
        first = False
    plt.semilogx(dSt, eff, label=label)

plt.xlabel(r'$\mathrm{St} - \mathrm{St}_\mathrm{c}$')
plt.ylabel('efficiency $\lambda$')
plt.legend(loc='best')

plt.show()

