#!/usr/bin/env python3

import numpy as np, matplotlib.pyplot as plt
from flows import shm, power

import argparse
parser = argparse.ArgumentParser(description='calculate contact angle for critical trajectories for particles being transported by flow (u,v) = (x**m, m*x**(m-1)*y)')
parser.add_argument('-x', '--x', type=float,
                    help='starting position in x. If not set will default to default value for chosen flow field.')
parser.add_argument('-s', '--stokes', nargs=3, default=[1e-3, 1, 4],
                    help='distribution of Stokes numbers *above* critical value as arguments to np.geomspace(...)')
parser.add_argument('--maxstep', default=1e-2, type=float,
                    help='max step size during integrations (default=1e-2)')
parser.add_argument('--tmax', default=1e2, type=float,
                    help='max time before aborting integrations (default=1e2)')
parser.add_argument('--niters', default=50, type=int,
                    help='number of refinement iterations to estimate efficiency (default=50)')
parser.add_argument('-H', '--half', action='store_true',
                    help='evaluate trajectories at half the critical value of y (instead of critical trajectories)')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='show iteration steps')
parser.add_argument('-li', '--left-inset', action='store_true',
                    help='if true, move inset to (upper) left')
parser.add_argument('-m0', '--m0', default=1, type=float,
                    help='reference value of m')
parser.add_argument('-m', '--m', nargs=3, default=[1e-2, 1, 21],
                    help='distribution of values of dm (m = m0+dm) to take in range [0, 1] as arguments to np.geomspace(...)')
args = parser.parse_args()
quiet = not args.verbose

figsize = (3.375, 3)
plt.figure(figsize=figsize)
ax = plt.gca()

import ast
dSt = np.geomspace(*[ast.literal_eval(str(s)) for s in args.stokes])
dM = np.geomspace(*[ast.literal_eval(str(m)) for m in args.m])
M = args.m0 + dM

def critical_stokes(m):
    if m == 1: flow = shm.FlowField()
    else: flow = power.FlowField(m=m)

    x0 = args.x if args.x is not None else flow.default_starting_distance
    x0 = np.longdouble(x0)
    Stc = flow.critical_stokes(x0, tmax=args.tmax, max_step=args.maxstep, niters=args.niters, quiet=quiet, message_header='critical stokes dm={:.4f}'.format(m - args.m0))

    contact_angles = np.empty(dSt.size)
    for i, delta in enumerate(dSt):
        St = Stc + delta

        y0, y1 = flow.capture_efficiency(x0, St, return_bounds=True, tmax=args.tmax, max_step=args.maxstep, niters=args.niters, quiet=quiet, message_header='capture efficiency dm={:.4f}'.format(m - args.m0))
        y = 0.25*(y0+y1) if args.half else y0
        _, traj = flow.trajectory([x0, y], St, tmax=args.tmax, max_step=args.maxstep)
        u, v = traj[-2:,-1]
        assert u <= 0
        contact_angles[i] = np.arctan2(-u, np.abs(v))

    return Stc, contact_angles

from multiprocessing import Pool
with Pool(M.size) as pool:
    Stc, contact_angles = list(zip(*pool.map(critical_stokes, M)))
    Stc = np.array(Stc)
    contact_angles = np.array(contact_angles)

print(np.hstack([np.array(Stc).reshape(-1,1), contact_angles]))

#import sys; sys.exit(0)

first = True
for delta, angle in zip(dSt, contact_angles.T):
    label = '{:.2g}'.format(delta)
    if first: label = r'$\mathrm{St} - \mathrm{St}_c = %s$' % label
    ax.plot(dM, angle/np.pi, label=label)
    first = False

ax.legend(loc='upper right', fontsize=8)
ax.set_ylabel(r'$\theta / \pi$')
ax.set_xlabel('$m-{}$'.format(args.m0))
ax.set_xscale('log')
ax.set_yscale('log')

#ax.set_xlim([1e-2, 1])

inset_width = 0.3
gap = 0.05
aspect = figsize[1] / figsize[0]
if args.left_inset:
    inset = ax.inset_axes([gap*aspect, 1-gap - inset_width, inset_width, inset_width])
else:
    inset = ax.inset_axes([1 - inset_width - gap*aspect, gap, inset_width, inset_width])

M = np.concatenate([[1], M])
Stc = np.concatenate([[0.25], Stc])
inset.plot(M, Stc)
inset.set_ylim([0, 0.6])

labels = inset.get_xticklabels()
inset.set_xticklabels(labels, fontsize=8)
labels = inset.get_yticklabels()
inset.set_yticklabels(labels, fontsize=8)
inset.set_xlabel('$m$', labelpad=-3, fontsize=8)
inset.set_ylabel(r'$\mathrm{St}_c$', labelpad=-3, fontsize=8)

if args.left_inset:
    inset.tick_params(right=True, labelright=True, left=False, labelleft=False)
    inset.yaxis.set_label_position('right')
else:
    inset.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    inset.xaxis.set_label_position('top')

plt.show()
