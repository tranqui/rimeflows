#!/usr/bin/env python3

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
        if not quiet: print('initial refinement: [%.4f %.4f] guess=%.4f' % (xlower, xupper, xguess))
        xguess *= 2

    for i in range(niters):
        xguess = 0.5 * (xlower + xupper)
        if f(xguess): xlower = xguess
        else: xupper = xguess
        if not quiet: print('refining %d/%d: [%.4f %.4f] guess=%.4f' % (i+1, niters, xlower, xupper, xguess))

    return xlower, xupper
