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

def logical_refine(f, niters, xmin=0, xmax=np.inf, xguess=1, quiet=True, message_header=''):
    """Refine location of transition of boolean from True to False.

    Args:
        f: objective function we are determining sign change for.
        niters: number of refinement iterations to estimate change over.
        xmin: initial lower bound for sign change.
        xmax: initial upper bound for sign change.
        xguess: initial guess to establish starting bounds.
        quiet: if True, will suppress iteration updates.
        message_header: preface to iteration updates if not quiet.
    Returns:
        Lower and upper bounds for where transition occurs.
    """

    xlower = xmin
    xupper = xmax

    assert f(xlower)
    assert xguess < xmax

    while xguess < xupper:
        if f(xguess): xlower = xguess
        else: xupper = xguess
        if not quiet: sys.stderr.write('{} initial refinement: [{:.4g} {:.4g}] guess={:.4g}\n'.format(message_header, xlower, xupper, xguess))
        xguess *= 2

    for i in range(niters):
        xguess = 0.5 * (xlower + xupper)
        if f(xguess): xlower = xguess
        else: xupper = xguess
        if not quiet: sys.stderr.write('{} refining {}/{}: [{:.4g} {:.4g}] guess={:.4g}\n'.format(message_header, i+1, niters, xlower, xupper, xguess))

    if not quiet: sys.stderr.write('{} best guess: {:.8g}\n'.format(message_header, 0.5*(xlower + xupper)))

    return xlower, xupper
