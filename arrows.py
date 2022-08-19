#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt

def add_arrow(line, x=None, y=None, direction='forward', size=7, color=None, eps=1e-4, **kwargs):
    """
    Add an arrow to a line, for showing the direction of lines in arrow plots (for streamlines etc).

    line:       Line2D object
    x,y:   x and/or y-position of the arrow. If both are None, mean of data is taken
    direction:  'forward' or 'backward'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    # find closest index
    if x is None and y is None:
        dr = np.linalg.norm([xdata - xdata.mean(), ydata - ydata.mean()], axis=0)
        start_ind = np.argmin(dr)

        if direction == 'forward':
            end_ind = start_ind + 1
        else:
            end_ind = start_ind - 1

        start = np.array((xdata[start_ind], ydata[start_ind]))
        end = np.array((xdata[end_ind], ydata[end_ind]))
        direct = (end - start)
        direct /= np.linalg.norm(direct)
        points = [(start, direct)]

    elif y is None:
        delta = np.sign(xdata - x)
        sign_change = ((np.roll(delta, 1) - delta) != 0).astype(int)
        sign_change[0] = 0

        points = []
        for index in np.where(sign_change)[0]:
            if direction == 'forward':
                start = index - 1
                end = index
            else:
                start = index
                end = index - 1

            start = np.array((xdata[start], ydata[start]))
            end = np.array((xdata[end], ydata[end]))
            direct = end - start
            if np.linalg.norm(direct) < eps: continue
            direct /= np.linalg.norm(direct)

            if start[0] < end[0]:
                y = np.interp(x, [start[0], end[0]], [start[1], end[1]])
            else:
                y = np.interp(x, [end[0], start[0]], [end[1], start[1]])
            position = np.array([x, y])
            points += [(position, direct)]

    elif x is None:
        delta = np.sign(ydata - y)
        sign_change = ((np.roll(delta, 1) - delta) != 0).astype(int)
        sign_change[0] = 0

        points = []
        for index in np.where(sign_change)[0]:
            if direction == 'forward':
                start = index - 1
                end = index
            else:
                start = index
                end = index - 1

            start = np.array((xdata[start], ydata[start]))
            end = np.array((xdata[end], ydata[end]))
            direct = end - start
            if np.linalg.norm(direct) < eps: continue
            direct /= np.linalg.norm(direct)

            if start[1] < end[1]:
                x = np.interp(y, [start[1], end[1]], [start[0], end[0]])
            else:
                x = np.interp(y, [end[1], start[1]], [end[0], start[0]])
            position = np.array([x, y])
            points += [(position, direct)]

    else:
        dr = np.linalg.norm([xdata - x, ydata - y], axis=0)
        start_ind = np.argmin(dr)

        if direction == 'forward':
            end_ind = start_ind + 1
        else:
            end_ind = start_ind - 1

        start = np.array((xdata[start_ind], ydata[start_ind]))
        end = np.array((xdata[end_ind], ydata[end_ind]))
        direct = (end - start)
        direct /= np.linalg.norm(direct)
        points = [(start, direct)]

    size *= line.get_linewidth()
    for pos, direct in points:
        # size=size*line.get_linewidth()
        # line.axes.arrow(start[0], start[1], direction[0], direction[1],
        #                 length_includes_head=True,
        #                 head_width=0.1, head_length=0.3,
        #                 fc=color, ec=color, lw=0.)
        #                 #head_starts_at_zero=True)

        start = pos
        end = pos + eps*direct
        arrow = line.axes.annotate('',
            xytext=start, xy=end,
            #arrowprops=dict(arrowstyle='simple', color=color),
            arrowprops=dict(headwidth=size, headlength=size, lw=0., fc=color),
            **kwargs
        )
        arrow.set_in_layout(False)
        arrow.arrow_patch.set_clip_box(line.axes.bbox)
        #line.axes.plot(pos[0], pos[1], 'o', c=pl.get_color())
