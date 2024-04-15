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

import argparse
import numpy as np, matplotlib.pyplot as plt

from flows import kuwabara, inviscid, rg, stokes, shm, hiemenz, power, chord, eliipse
flow_fields = {'kuwabara': kuwabara, 'inviscid': inviscid, 'rg': rg, 'stokes': stokes, 'shm': shm, 'hiemenz': hiemenz, 'power': power, 'chord': chord, 'ellipse': ellipse}

parser = argparse.ArgumentParser(description='calculate critical inertia for point particles to be captured at stagnation point flows')
parser.add_argument('-f', '--flow', default='kuwabara', choices=flow_fields.keys(),
                    help='flow field to use (default=kuwabara)')
parser.add_argument('-x', '--x', type=float,
                    help='starting position in x. If not set will default to default value for chosen flow field.')
parser.add_argument('--maxstep', default=np.inf, type=float,
                    help='max step size during integrations (default=infinite)')
parser.add_argument('--tmax', default=1e4, type=float,
                    help='max time before aborting integrations (default=1e4)')
parser.add_argument('--niters', default=50, type=int,
                    help='number of refinement iterations to estimate efficiency (default=50)')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='show iteration steps')
parser.add_argument('flowparams', nargs='*',
                    help='parameters to pass to flow field')
args = parser.parse_args()
quiet = not args.verbose

args.flowparams = [eval(p) for p in args.flowparams]
flow_module = flow_fields[args.flow]
flow = flow_module.FlowField(*args.flowparams)

if args.x is None: args.x = flow.default_starting_distance

Stc = flow.critical_stokes(args.x, tmax=args.tmax, max_step=args.maxstep, niters=args.niters,
                           quiet=quiet, message_header='critical stokes param={}'.format(repr(args.flowparams)))

flowparams = ' '.join([repr(element) for element in args.flowparams])
print('{} {:^16.12f}'.format(flowparams, Stc))
