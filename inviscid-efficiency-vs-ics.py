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
# along with this program.  If not, see <https://www.gnu.org/licenses/>

import sys
from natsort import natsorted

import matplotlib.pyplot as plt, numpy as np
plt.style.use('figstyle.mplstyle')

first = True
for path in natsorted(sys.argv[1:]):
    gamma = float(path.split('.csv')[0].split('=')[-1])
    print(gamma)

    delta_St, St, y = np.genfromtxt(path).T

    label = '{.1f}'.format(gamma)
    if first:
        label = '$\gamma={}$'.format(label)
        first = False

    plt.plot(delta_St, 2*y, lw=0.5, label=label)

plt.xscale('log')
plt.yscale('log')

x = np.geomspace(1e-3, 1, 100)
f = lambda x: 10*x**0.5
plt.plot(x, f(x), 'k--', lw=0.5)
from mpltools import annotation
x = 1e-2
annotation.slope_marker((x, f(x)), 0.5, text_kwargs={'fontsize': 8.})

plt.ylim([1e-3, 10])
plt.legend(loc='upper left', ncol=4, borderpad=0)
plt.xlabel('$\mathrm{St} - \mathrm{St}_c$')
plt.ylabel('efficiency $2y_0$')
plt.show()
