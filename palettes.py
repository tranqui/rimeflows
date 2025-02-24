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

class PurpleGold:
    """Colours generated from:
    https://coolors.co/palette/4d0e5a-871c90-cc67c2-dabbe3-efda4d
    """

    palette = ['#4d0e5a', '#871c90', '#cc67c2', '#dabbe3', '#efda4d', '#ffffc0']

    cylinder_colour = 'grey'

    flow_colour = 'white'
    streamline_colour = 'steelblue'

    critical_trajectory_colour = 'red'
    separatrix_colour = 'black'
    nullcline_colour = 'black'

    colliding_trajectory_colour = palette[2]
    colliding_manifold_colour = palette[3]

    noncolliding_trajectory_colour = palette[-2]
    noncolliding_manifold_colour = palette[-1]

class BlueRed:
    cylinder_colour = 'grey'

    flow_colour = '#e6e6ff'
    streamline_colour = 'steelblue'

    critical_trajectory_colour = 'crimson'
    separatrix_colour = 'black'
    nullcline_colour = 'black'

    colliding_trajectory_colour = 'red'
    colliding_manifold_colour = 'lightsalmon'

    noncolliding_trajectory_colour = 'blue'
    noncolliding_manifold_colour = '#e6e6ff'
