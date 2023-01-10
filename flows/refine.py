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

import sys
import numpy as np

def logical_refine(f, niters, xmin=0, xmax=np.inf, quiet=True):
    """Refine location of transition of boolean from True to False.

    Args:
        f: objective function we are determining sign change for.
        niters: number of refinement iterations to estimate change over.
        xmin: initial lower bound for sign change.
        xmax: initial upper bound for sign change.
        quiet: if True, will suppress iteration updates.
    Returns:
        Lower and upper bounds for where transition occurs.
    """

    xlower = xmin
    xupper = xmax

    assert f(xlower)

    xguess = 1
    while xguess < xupper:
        if f(xguess): xlower = xguess
        else: xupper = xguess
        if not quiet: sys.stderr.write('initial refinement: [%.4f %.4f] guess=%.4f\n' % (xlower, xupper, xguess))
        xguess *= 2

    for i in range(niters):
        xguess = 0.5 * (xlower + xupper)
        if f(xguess): xlower = xguess
        else: xupper = xguess
        if not quiet: sys.stderr.write('refining %d/%d: [%.4f %.4f] guess=%.4f\n' % (i+1, niters, xlower, xupper, xguess))

    if not quiet: sys.stderr.write('best guess: %.8f\n' % (0.5*(xlower + xupper)))

    return xlower, xupper
