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
import re
from itertools import chain

def triangulate_grid(*args):
    """Triangulate a 2d grid of coordinates to obtain a triangle mesh.

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

def grid_lines(*args):
    """Coordinates of lines along a 2d grid of coordinates to obtain a triangle mesh.

    Args:
        Meshgrid coordinates (e.g. from output of numpy.meshgrid).
    Returns:
        xcoords: sequence of coordinates for first set of grid lines.
        ycoords: sequence of coordinates for second set of grid lines.
    """

    nrows, ncols = args[0].shape

    xcoords, ycoords = [], []
    for row in range(nrows): xcoords += [np.array([x[row] for x in args]).T]
    for col in range(ncols): ycoords += [np.array([x[:,col] for x in args]).T]
    return xcoords, ycoords

def pov_vector(v):
    """Convert numpy vector into form parseable by povray.

    Args:
        v: numpy array.
    Returns:
        String representation of vector ready for povray.
    """
    if type(v) is str: return v
    else: return str(v.tolist()).replace('[', '<').replace(']', '>')

def to_snake_case(s):
    """Convert PascalCase or camelCase (AKA mixedCase in PEP 8) string to snake_case.

    snake_case is the povray object type naming convention, whereas we make use of PascalCase for
    class definitions. We can scrape the class names and convert to snake_case to create
    a 1:1 correspondence between python class and povray objects.

    >>> to_snake_case('myStringInitiallyInCamelCase')
    'my_string_initially_in_camel_case'

    Args:
        s: string in PascalCase or camelCase.
    Returns:
        String in snake_case.
    """

    return re.sub(r'(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z])', '_', s).lower()

class Primitive:
    def __init__(self, *args, single_line=False, **kwargs):
        self.single_line = single_line
        self.children = []
        for arg in args:
            try: self.children += arg
            except: self.children += [arg]
        self.children = list(chain(self.children)) # flatten list of lists

    @property
    def header(self):
        return to_snake_case(self.__class__.__name__)

    @property
    def body_open(self):
        return '{'

    @property
    def body_close(self):
        return '}'

    @property
    def body(self):
        raise NotImplementedError

    def write(self, f=sys.stdout, nest=0):
        """Generate the povray script.

        Args:
            f: filestream to output to.
        """

        indent = '  ' * nest
        f.write('{}{}'.format(indent, self.header))

        if len(self.children):
            if len(self.body_open): f.write('\n{}{}'.format(indent, self.body_open))
            for child in self.children:
                f.write('\n')
                child.write(f, nest+1)
            f.write('\n{}{}'.format(indent, self.body_close))

        else:
            f.write(' {} {} {}'.format(self.body_open, self.body, self.body_close))

    def __repr__(self):
        f = io.StringIO()
        self.write(f)
        return f.getvalue()

# Basic combinations of objects:
class Union(Primitive): pass
class Merge(Primitive): pass
class Intersection(Primitive): pass
class Difference(Primitive): pass

class Sphere(Primitive):
    def __init__(self, position, radius, *args, **kwargs):
        self.position = pov_vector(position)
        self.radius = str(radius)
        super().__init__(*args, **kwargs)

    @property
    def body(self):
        return '{}, {}'.format(self.position, self.radius)

class Cylinder(Primitive):
    def __init__(self, position1, position2, radius, *args, **kwargs):
        self.position1 = pov_vector(position1)
        self.position2 = pov_vector(position2)
        self.radius = str(radius)
        super().__init__(*args, **kwargs)

    @property
    def body(self):
        return '{}, {}, {}'.format(self.position1, self.position2, self.radius)

# Derivative structures.

def line(coordinates, line_width, *args, **kwargs):
    """Generates a line in 3d as the union of cylinders that can be rendered with ray-tracing."""

    objects = []

    # Create the line as a series of cylinders joining adjacent points.
    for v1, v2 in zip(coordinates, coordinates[1:]):
        objects += [Cylinder(pov_vector(v1), pov_vector(v2), line_width)]

    # Round the edges of the cylinders where the lines join to make it smooth.
    for v in coordinates:
        objects += [Sphere(pov_vector(v), line_width)]

    return Merge(*objects, *args, **kwargs)

class VectorBundle(Primitive):
    def __init__(self, coordinates, *args, **kwargs):
        self.coordinates = coordinates
        self.vectors = [pov_vector(x) for x in self.coordinates]
        super().__init__(*args, **kwargs)

    @property
    def body(self):
        return '{}, {}'.format(len(self.vectors), ', '.join(self.vectors))

class VertexVectors(VectorBundle): pass
class NormalVectors(VectorBundle): pass
class FaceIndices(VectorBundle): pass

class Mesh2(Primitive):
    """Generates a povray "mesh2" object ready to be rendered with ray-tracing."""

    def __init__(self, coordinates, triangulation, *args, **kwargs):
        """Create the mesh from raw data.

        Args:
            coordinates: vertex coordinates.
            triangulation: vertex indices.
        """
        triangulation = np.atleast_2d(triangulation)
        ntriangles, d = triangulation.shape
        assert d == 3

        coordinates = np.array(coordinates).reshape(-1,d)
        A, B, C = coordinates[triangulation.T]
        self.face_normals = np.cross(B - A, C - A)
        vertex_normals = np.zeros(coordinates.shape)

        for vertices in triangulation.T:
            vertex_normals[vertices] += self.face_normals

        vertex_normals = (vertex_normals.T / np.linalg.norm(vertex_normals, axis=1)).T

        super().__init__(*args, **kwargs)
        self.vertex_vectors = VertexVectors(coordinates)
        self.normal_vectors = NormalVectors(vertex_normals)
        self.face_indices = FaceIndices(triangulation)
        self.children += [self.vertex_vectors, self.normal_vectors, self.face_indices]

# Variable declarations.

class Macro(Primitive):
    def __init__(self, name, *args, arguments=[], **kwargs):
        self.name = name
        self.arguments = arguments
        super().__init__(*args, **kwargs)

    @property
    def header(self):
        return '#{directive} {name}({arguments})'.format(directive=super().header, name=self.name,
                                                         arguments=','.join(self.arguments))

    @property
    def body_open(self):
        return ''

    @property
    def body_close(self):
        return '#end'

class Declare(Primitive):
    def __init__(self, name, *args, **kwargs):
        self.name = name
        self.arguments = []
        super().__init__(*args, **kwargs)

    @property
    def header(self):
        return '#{directive} {name}'.format(directive=super().header, name=self.name)

    @property
    def body_open(self):
        return ' = '

    @property
    def body_close(self):
        return ';'

def Local(Declare): pass

if __name__ == '__main__':
    spheres1 = Union([Sphere(np.random.random(3), 0.5) for i in range(3)])
    print(spheres1)
    spheres2 = Union(*[Sphere(np.random.random(3), 0.5) for i in range(3)])
    print(spheres2)
    print(Merge(spheres1, spheres2))

    print(Macro('abc', Sphere(np.random.random(3), 'line_width'), arguments=['line_width']))
    print(Declare('abc', Sphere(np.random.random(3), 0.5)))

    N = 25
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)
    Z = X*(1-X) * Y*(1-Y)

    coordinates, triangles = triangulate_grid(X, Y, Z)
    print(Mesh2(coordinates, triangles))

    # import matplotlib.pyplot as plt
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.plot_trisurf(*coordinates.T, triangles=triangles)
    # plt.show()

