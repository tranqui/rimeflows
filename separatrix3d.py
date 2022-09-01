#!/usr/bin/env python3

# separatrix3d.py -- calculate 3d separatrix for toy stagnation point flow.

# Copyright (c) 2022 Patrick B Warren <patrick.warren@stfc.ac.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import numpy as np
import pandas as pd
from numpy import pi as pi
from numpy import sqrt, cos, abs
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from povray import PovrayLine, PovrayMesh2, triangulate_grid

parser = argparse.ArgumentParser(description='separatrix construction in 3d')
parser.add_argument('--tmax', default=20.0, type=float, help='max time, default 20.0')
parser.add_argument('--nsvals', default=50, type=int, help='number of arc points, default 50')
parser.add_argument('--nyvals', default=100, type=int, help='number of y0 points, default 100')
parser.add_argument('--ymax', default=1.5, type=float, help='max y0, default 1.5')
parser.add_argument('--eps', default=0.01, type=float, help='start offset, default 0.01')
parser.add_argument('--method', default='LSODA', help='integration method, default LSODA')
parser.add_argument('--tol', default=1e-12, type=float, help='integrator tolerance, default 1e-12')
parser.add_argument('--mayavi', action='store_true', help='use mayavi to render')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--povray', action='store_true', help='produce a povray input script for the separatrix triangulation')
group.add_argument('--surface', action='store_true', help='use a regular surface')
group.add_argument('--triangulate', action='store_true', help='use a triangulated surface')
group.add_argument('--wireframe', action='store_true', help='use a wireframe')
group.add_argument('--mesh', action='store_true', help='show the underlying (y0, t) mesh')
parser.add_argument('--zero', action='store_true', help='include the zero acceleration surface')
parser.add_argument('--nxzero', default=21, type=int, help='number of x points in above, default 21')
parser.add_argument('--nyzero', default=21, type=int, help='number of y points in above, default 21')
parser.add_argument('-o', '--output', help='save plot to, eg, pdf')
args = parser.parse_args()

args.nsvals = 2*(args.nsvals//2) # this has to be even for the slice function to work properly

def flatten_everything(xl): # used below
    '''flatten a list of numpy arrays'''
    return [x.flatten() for x in xl]

# The original equations are St [x" - (y')^2] = - x' - k x^2 and 0 = - y' + 2 k x y.

# If we set x → x k St ; y → y sqrt(k St) ; t → t / St then the parameters disappear and
# the equations become x" - (y')^2 = - x' + vx and y' = vy where (vx, vy) = (-x^2, 2xy).
# To convert this to first order ODEs set u = x' so that
# dx/dt = u , du/dt = vy^2 - u + vx , dy/dt = vy.

# For the backwards-in-time calculation, the right-hand sides are negated.

# For meshing purposes we define the metric ds^2 = dx^2 + du^2 + dy^2 and
# integrate the arc length s of the trajectory alongside the ODEs.

# The arc length slicing function below is a cosine-based smooth
# function of s which vanishes at evenly-spaced intermediate points
# determined by dL.  The start (s = 0) and end (s = L) are omitted by
# the conditional evaluation.

def eom_backwards(t, xi): # as described above, with negated RHS
    x, u, y, s = xi
    vx, vy = -x**2, 2*x*y
    #vx, vy = -x**2 - 0.25*x, 2*x*y + 0.25*y
    dx, du, dy = u, vy**2 - u + vx, vy
    ds = np.sqrt(dx**2 + du**2 + dy**2)
    return -dx, -du, -dy, ds

def acceleration(t, xi): # acceleration in the x-direction
    x, u, y, s = xi
    return u + x**2

def arcslicer(t, xi): # slice condition on backwards trajectory
    x, u, y, s = xi
    return cos(pi*(s/dL - 0.5)) if abs(s - 0.5*L) < 0.5*(L-dL) else 1.0

def stop(t, xi): # stopping condition on backwards trajectory
    x, u, y, l = xi
    return u # change of sign; other stopping conditions are possible

stop.terminal = True

# In the below, we first integrate the equation of motion backwards in
# time without keeping any intermediate values (t_eval is set to an
# empty list), past the zero acceleration condition, until we reach
# the stopping condition.  We record the times and positions of these
# events, and the total arc length L of the trajectory, then use the
# stopping time as the end point for a second integration, now
# capturing the position (x, u, y, s) at a fixed number of
# intermediate arc length positions, with spacing dL computed from L.

# This is repeated, incrementing the starting position y0 along the
# y-axis.  This generates a set of trajectories each of which consists
# of (x, u, y, s) values at a set of discrete arc length positions,
# with the same number of points for each trajectory.  This
# discretises the mapping (s, y0) -> (x, u, y) for the separatrix
# surface, or in other words is equivalent to parametrising the
# surface as x(s, y0), u(s, y0) and y(s, y0) on a (s, y0) mesh.

# For wireframe and simple surface plots, the underlying (s, y0) mesh
# is not required.  For triangulation, and of course for displaying
# the mesh itself (--mesh option), it is required.

solns = []

for y0 in np.linspace(0, args.ymax, args.nyvals):

    start = np.array([args.eps, -args.eps, y0, 0.0])

    soln = solve_ivp(eom_backwards, [0, args.tmax], start, t_eval=[],
                     events=[acceleration, stop],
                     method=args.method, atol=args.tol, rtol=args.tol)

    intersect = soln.y_events[0][0] # keep the intersection with the zero-acceleration condition
    tend = soln.t_events[1][0] # the end time
    L = soln.y_events[1][0][3] # the total arc length
    dL = L / (args.nsvals - 1) # discrete spacing used in the arc length slicing function

    soln = solve_ivp(eom_backwards, [0, tend], start, t_eval=[0, tend], events=[arcslicer],
                     method=args.method, atol=args.tol, rtol=args.tol)

    start = soln.y.transpose()[0] # this should actually be the same as above 'start' above
    middle = [x for x in soln.y_events[0].transpose()] # the points captured by the arc length slicer
    end = soln.y.transpose()[1] # and the end point

    # Build the trajectory by adding the start and end points to the
    # points captured by the slice function.

    traj = [np.append(np.insert(middle[i], 0, start[i]), end[i]) for i in range(4)]

    # The data structure (solns) is a list of lists containing the
    # zero-acceleration intersection points and trajectories.

    solns.append([intersect, traj])

# We now wrangle this data into a three-dimensional separatrix curve
# in (x, u, y) for the intersection of the backwards-in-time
# trajectories with the zero-acceleration initial condition.

intersect = np.stack([soln[0] for soln in solns]).transpose()
xintersect, uintersect, yintersect = intersect[0:3] # drop the arc length

# Similarly, two-dimensional arrays representing the separatrix
# surface in (x, u, y) are built to match the two-dimensional mesh in
# (s, y0), where s is the set of arc length positions along each
# trajectory and y0 is the set of initial starting points.

xsurf = np.stack([soln[1][0] for soln in solns])
usurf = np.stack([soln[1][1] for soln in solns])
ysurf = np.stack([soln[1][2] for soln in solns])

smesh = np.stack([soln[1][3] for soln in solns])

# Extract the starting positions and expand them out into a
# two-dimensional array to match the above.

y0vals = [soln[1][2][0] for soln in solns] # this is a list of values
y0mesh = np.vstack([y0vals] * smesh.shape[1]).transpose() # replicated

# At this point (xsurf, usurf, ysurf) are two-dimensional arrays which
# represent, parametrically, the surface x(s, y0), u(s, y0), y(s, y0)
# in phase space, and (smesh, y0mesh) is the underlying mesh in the
# (s, y0) parameter space.

fig = plt.figure()

if args.mesh: # display the parameter mesh (2d)

    ax = fig.add_subplot()

    smesh, y0mesh = flatten_everything([smesh, y0mesh])
    ax.plot(smesh, y0mesh, '.')

    ax.set_xlabel('s')
    ax.set_ylabel('y0')

elif args.mayavi: # 3d plots with mayavi

    from mayavi import mlab

    dpi = 300
    figsize = (dpi*3.375, dpi*3)
    mlab.figure(size=figsize)

    # Show the separatrix as a 3d surface.
    representation = 'surface'
    if args.wireframe: representation = 'wireframe'
    surface = mlab.mesh(-xsurf, 5*ysurf, usurf, representation=representation, colormap='viridis')

    # Show zero-acceleration condition.
    intersection = mlab.plot3d(-xintersect, 5*yintersect, uintersect, color=(1,0,0), tube_radius=None)
    intersection.actor.property.line_width = 2

    # Draw stagnation point.
    mlab.points3d(0, 0, 0, color=(1,1,1), resolution=100, scale_factor=0.1)

    # Orthographic projection is cleaner for showing aspect ratios.
    mlab.gcf().scene.parallel_projection = True

    # Label axes.
    ax = mlab.axes(surface, color=(0,0,0), xlabel='$x$', ylabel='$y$', zlabel='$u$')
    ax.axes.label_format = ''

else: # 3d plots with matplotlib

    if args.povray:
        coordinates = np.array((xintersect, 10*yintersect, uintersect)).T
        print(PovrayLine(coordinates, name='zeroAcceleration'))
        coordinates, triangles = triangulate_grid(xsurf, 10*ysurf, usurf)
        print(PovrayMesh2(coordinates, triangles, name='separatrix'))

        for row in range(len(xsurf)):
            coordinates = np.array((xsurf[row], 10*ysurf[row], usurf[row])).T
            #import sys; sys.stderr.write(repr(coordinates.shape))
            print(PovrayLine(coordinates, name=('line{}'.format(row))))
        import sys; sys.exit(0)

    ax = fig.add_subplot(projection='3d')

    if args.triangulate:

        coordinates, triangles = triangulate_grid(xsurf, usurf, ysurf)
        ax.plot_trisurf(*coordinates.T, triangles=triangles)

    elif args.wireframe: # this is much easier !!

        ax.plot_wireframe(xsurf, usurf, ysurf)

    elif args.surface: # and this !!!

        ax.plot_surface(xsurf, usurf, ysurf, color='w')

    if args.zero:

        # Add a second surface representing the zero-acceleration
        # initial condition, being u + x^2 = 0.  Since there is no
        # y-dependence there only really needs to be two points in the
        # y-direction at the minimum and maximum values, but the
        # number of points in this direction can of course be
        # increased for aesthetic reasons.  Unfortunately adding this
        # second surface doesn't work very well because matplotlib
        # does not do true 3d rendering ('z-order' problem).

        xinit = np.linspace(0, np.sqrt(2.0), args.nxzero)
        yinit = np.linspace(-1.5, 1.5, args.nyzero)
        xinit, yinit = np.meshgrid(xinit, yinit)
        uinit = - xinit**2

        if args.wireframe:
            ax.plot_wireframe(xinit, uinit, yinit, color='g')
        else:
            ax.plot_surface(xinit, uinit, yinit, color='g')

    # Finally, add the intersection curve where the zero-acceleration
    # initial condition meets the separatrix surface foliated by the
    # backwards-in-time trajectories.  This works reasonably well but
    # there is still a problem of the 'z-order' of the objects.

    ax.plot3D(xintersect, uintersect, yintersect, 'r', linewidth=2)

    # Set the plotting limits appropriately

    ax.set_xlim(0, 2) # x-axis
    ax.set_ylim(-2, 0) # u-axis
    ax.set_zlim(-1, 1) # y-axis

    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_zlabel('y')

if args.output:
    plt.savefig(args.output)
    print('created', args.output)
elif args.mayavi:
    mlab.show()
else:
    plt.show()
