#!/usr/bin/env python3
# Copyright (C) 2025 Joshua Robinson
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

import numpy as np

try:
    from . import wedgestokes as stokes
    from . import wedgeinviscid as inviscid
except:
    import wedgestokes as stokes
    import wedgeinviscid as inviscid

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print('critical (i.e. m=1) alpha for Stokes flow is alpha={:.4f}'.format(stokes.critical_alpha()))

    print('Stokes:')
    for alpha in np.linspace(0, 1, 5):
        print('    alpha={:.4f}: m={:.4f}'.format(alpha, stokes.wedge_exponent(alpha)))

    print('inviscid:')
    for alpha in np.linspace(0, 1, 5):
        print('    alpha={:.4f}: m={:.4f}'.format(alpha, inviscid.wedge_exponent(alpha)))

    f_stokes = np.vectorize(stokes.wedge_exponent)
    f_inviscid = np.vectorize(inviscid.wedge_exponent)

    # Plot in finer detail:
    # alpha = np.linspace(0.5, 0.6, 1000)
    # m = f_stokes(alpha)
    # plt.plot(alpha, m-1)
    # plt.axhline(y=0)
    # plt.xlabel(r'$\alpha$')
    # plt.ylabel(r'$m$')
    # plt.show()

    plt.figure(figsize=(3.375, 2.5))
    alpha = np.linspace(0, 1, 100)
    m_stokes = f_stokes(alpha)
    m_inviscid = f_inviscid(alpha)
    plt.plot(alpha * 180, m_stokes, lw=1, label='Stokes')
    plt.plot(alpha * 180, m_inviscid, 'b--', lw=1, label='inviscid')
    plt.legend(loc='best')
    plt.axhline(y=1, ls='dotted')
    plt.xlabel(r'$\Phi \, (\si{\degree})$')
    plt.ylabel('$m$')
    plt.xlim([0, 180])
    plt.xticks(np.linspace(0, 180, 7))
    plt.ylim([0, 2])
    plt.savefig('mofPhi.pdf')
    plt.savefig('mofPhi.png')
    plt.show()

    # m = [1, 1.5, 2]
    # fig, axes = plt.subplots(nrows=len(m), sharex=True, sharey=True)

    # for ax, m in zip(axes, m):
    #     flow = FlowField(m=m)

    #     L = 3
    #     x0 = 1.5*L
    #     for y0 in np.linspace(0, L**0.5, 11)**2:
    #         x, y = flow.streamline([x0, y0], max_step=0.1)
    #         ax.plot(x, y, 'k-', lw=0.5)
    #         ax.plot(x, -y, 'k-', lw=0.5)

    #     ax.fill_between([-0.5, 0], -L, L, fc='r', zorder=10)
    #     ax.set_ylabel('$y$')

    # plt.xlabel('$x$')
    # plt.xlim([-0.5, 1.5*L])
    # plt.ylim([-L, L])

    # plt.show()
