#!/usr/bin/env python3

# Helper utilities to export numpy data into povray for 3d rendering.
#
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
import io

def triangulate_grid(*args):
    """Triangulate a grid of coordinates to obtain a triangle mesh.

    This only works on gridded coordinates because they can be easily
    turned into quadstrips.

    Args:
        Meshgrid coordinates (e.g. from output of numpy.meshgrid).
    Returns:
        coordinates: vertex coordinates referenced by triangulation.
        triangulation: (n x 3) int array giving vertex indices for each triangle.
    """
    nrows, ncols = args[0].shape
    triangles_per_strip = 2*(ncols - 1)
    nstrips = nrows - 1
    ntriangles = nstrips * triangles_per_strip

    triangulation = np.empty((ntriangles, 3), dtype=int)
    # Two basic triangles in a strip:
    triangulation[::2] = [0, ncols, 1]
    triangulation[1::2] = [1, ncols, ncols+1]
    # Shift vertex indices along the appropriate number of columns within a strip.
    triangulation += np.tile(np.repeat(np.arange(ncols - 1), 2), nstrips)[:, np.newaxis]
    # Shift vertex indices along the rows to get the correct indices for each strip.
    triangulation += np.repeat(ncols*np.arange(nstrips), triangles_per_strip)[:,np.newaxis]

    coordinates = np.array([x.flatten() for x in args]).T
    return coordinates, triangulation

def pov_vector(v):
    return str(v.tolist()).replace('[', '<').replace(']', '>')

class PovrayLine:
    """Generates a line in 3d as the union of cylinders to be rendered with ray-tracing."""

    def __init__(self, coordinates, name='line'):
        self.coordinates = coordinates
        self.name = name

    def write(self, f=sys.stdout):
        """Generate the povray script.

        Args:
            f: filestream to output to.
        """
        f.write('#macro %s(lw)\n  merge\n  {' % self.name)
        for v1, v2 in zip(self.coordinates, self.coordinates[1:]):
            f.write('\n    cylinder { %s, %s, lw }' % (pov_vector(v1), pov_vector(v2)))
        for v in self.coordinates:
            f.write('\n    sphere { %s, lw }' % pov_vector(v))
        f.write('\n  }\n#end')

    def __repr__(self):
        f = io.StringIO()
        self.write(f)
        return f.getvalue()

class PovrayMesh2:
    """Generates a povray "mesh2" object ready to be rendered with ray-tracing."""

    def __init__(self, coordinates, triangulation, name='triangulation'):
        """Create the mesh from raw data.

        Args:
            coordinates: vertex coordinates.
            triangulation: vertex indices.
            name: name of mesh to be used in resulting povray script (i.e. the name of
                    the variable).
        """
        ntriangles, d = triangulation.shape
        assert d == 3

        self.coordinates = np.array(coordinates).reshape(-1,d)
        self.triangulation = triangulation
        self.name = name

        A, B, C = self.coordinates[self.triangulation.T]
        self.face_normals = np.cross(B - A, C - A)
        self.vertex_normals = np.zeros(self.coordinates.shape)

        for vertices in self.triangulation.T:
            self.vertex_normals[vertices] += self.face_normals

        self.vertex_normals = (self.vertex_normals.T / np.linalg.norm(self.vertex_normals, axis=1)).T

    def write(self, f=sys.stdout):
        """Generate the povray script.

        Args:
            f: filestream to output to.
        """
        f.write('#declare %s = mesh2\n{\n  vertex_vectors\n  {\n    %d' % (self.name, len(self.coordinates)))
        for v in self.coordinates: f.write(',\n    %s' % pov_vector(v))
        f.write('\n  }\n  normal_vectors\n  {\n    %d' % len(self.vertex_normals))
        for v in self.vertex_normals: f.write(',\n    %s' % pov_vector(v))
        f.write('\n  }\n  face_indices\n  {\n    %d' % len(self.triangulation))
        for v in self.triangulation: f.write(',\n    %s' % pov_vector(v))
        f.write('\n  }\n}')

    def __repr__(self):
        f = io.StringIO()
        self.write(f)
        return f.getvalue()

if __name__ == '__main__':
    N = 25
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)
    Z = X*(1-X) * Y*(1-Y)

    coordinates, triangles = triangulate_grid(X, Y, Z)
    print(PovrayMesh2(coordinates, triangles))

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(*coordinates.T, triangles=triangles)
    plt.show()

