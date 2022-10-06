#version 3.7;

#include "shapes.inc"
#include "colors.inc"
#include "functions.inc"

global_settings
{
  assumed_gamma 2.2
  max_trace_level 256
  ambient_light White
}

#declare noncollidingManifold = srgb <1, 1, 0.753>;         // #ffffc0
#declare noncollidingParticle = srgb <0.937, 0.855, 0.302>; // #efda4d

//#declare intersectionManifold = srgb <0.573, 0.365, 0.227>; // #925d3a
#declare intersectionManifold = srgb <1, 0, 0>;

#declare collidingManifold = srgb <0.855, 0.733, 0.89>; // #dabbe3
#declare collidingParticle = srgb <0.8, 0.404, 0.761>;  // #cc67c2

//background { colour White }

#declare shinyFinish = finish
{
  ambient 0.1
  diffuse 0.8
  phong 0.1
  phong_size 40
}

#declare matteFinish = finish
{
  ambient 1
  diffuse 0
}

#declare glassFinish = finish
{
  ambient 0.1
  diffuse 0.6
  phong 0.1
}

// Axes and bounding box vectors.
#declare v0 = <0 0 0>;
#declare vx = <1 0 0>;
#declare vy = <0 0.5 0>;
#declare vz = <0, 0, -0.5>;

#declare focusPos = 0.5 * (vx + vy + vz);
#declare offsetDirection = vnormalize(<0, 1, 0.9>);
#declare cameraPos = focusPos + 1.125*offsetDirection;

// // Uncomment to zoom in on one edge of box.
// #declare focusPos = vx+vy;
// #declare cameraPos = focusPos + <0.025, 0.025, 0.025>;

// #declare focusPos = v0;
// #declare cameraPos = focusPos + 0.1*<1 1 1>;

camera
{
  perspective
  //orthographic
  location cameraPos
  look_at  focusPos
  sky      <0, 0, 1>
  up       <0, 0.75, 0>
  right    <4/3, 0, 0>
}

#macro beam(position, intensity)
  light_source
  {
    position
    colour intensity*White
    parallel
    point_at focusPos
    shadowless
  }
#end

#macro illuminatedArea(position, v1, v2, n1, n2, intensity)
  light_source
  {
    position
    colour intensity*White
    area_light v1, v2, n1, n2
    area_illumination on
    adaptive 1
    shadowless
  }
#end

// Light up box from a mixture of angles to add some more shading to the 3d surface.
#if (AreaLighting)
  beam(cameraPos, 1.0)
  // Light up the front xz plane and right yz planes to cast soft shadows:
  illuminatedArea(vy, vx, vz, 10, 5, 0.5)
  illuminatedArea(vx, vy, vz, 5, 5, 0.5)
#else
  beam(cameraPos, 1.0)
  // Cheaper but leads to hard shadows: add two more parallel beams at oblique angles:
  beam(focusPos + <1.0, 0.25, 0.25>, 0.25)
  beam(focusPos + <0.25, 1.0, 0.25>, 0.25)
#end

// Frames for the domain boundaries.

#declare box_edge_radius = 0.0025;
#macro edge(v0, v1)
  cylinder { v0, v1, box_edge_radius }
#end
#macro node(v1)
  sphere { v1, box_edge_radius }
#end

// Black lines for labelled axis edges
merge
{
  edge(v0, vx)
  edge(v0, vy)
  node(v0)
  node(vx)
  node(vy)

  edge(vy, vy+vz)
  node(vy+vz)

  no_shadow
  pigment { colour Black }
  finish { matteFinish }
}

// Glass lines for unlabelled non-axis edges.
merge
{
  // edge(v0, vx)
  // edge(v0, vy)
  edge(v0, vz)

  edge(vx, vx+vy)
  edge(vy, vx+vy)
  edge(vx, vx+vz)
  edge(vz, vx+vz)
  //edge(vy, vy+vz)
  edge(vz, vy+vz)

  edge(vx+vy, vx+vy+vz)
  edge(vx+vz, vx+vy+vz)
  edge(vy+vz, vx+vy+vz)

  // Cap the ends with spheres so the edges join smoothly
  // node(v0)
  // node(vx)
  // node(vy)
  node(vz)
  node(vx+vy)
  node(vx+vz)
  //node(vy+vz)
  node(vx+vy+vz)

  // We apply some non-physical properties to these edges so the viewer unconsciously disregards
  // the edges: their presence acts as a guide to the eye to help get some 3d perspective when
  // looking at the image projected in 2d.
  // We first make the non-axis edges highly transparent, then prevent them from casting any
  // shadows on the scene to prevent optical artifacts from distracting the viewer.
  pigment { colour White transmit 0.9 }
  interior { ior 1.5 } // index of refraction of glass
  //interior { ior 2.5 } // index of refraction of diamond gives an ethereal-look
  no_shadow
  finish { glassFinish }
}

// Draw axes labels/directions

#macro axis_label(label, v0, v1)
  #local rawLabel = text
  {
    ttf "latinmodern-math.ttf" label
    0.001, 0
    scale 0.05
    pigment { colour Black }
  }
  #local com = 0.5*(max_extent(rawLabel) + min_extent(rawLabel));

  object
  {
    rawLabel
    translate -com
    rotate <90, 0, 0>  // by default place text in xz-plane rather than xy-plane
  }
#end

object
{
  #local axis = axis_label("x", v0, vx);
  #local height = (max_extent(axis) - min_extent(axis)).z;
  axis
  translate 0.5*vx
  translate height * vnormalize(-vz)
}

object
{
  #local axis = axis_label("y", v0, vy);
  #local height = (max_extent(axis) - min_extent(axis)).z;
  axis
  rotate <0, 0, -90>
  rotate <0, 0, 30>
  rotate <30, 0, 0>
  translate v0 + 0.5*vy
  translate height * vnormalize(-vz)
}

object
{
  #local axis = axis_label("-x'", v0, vz);
  #local length = (max_extent(axis) - min_extent(axis)).x;
  axis
  rotate <0, 0, -90>
  rotate <0, 0, 60>
  rotate <60, 0, 0>
  translate v0 + vy + 0.5*vz
  translate length * vnormalize(vy)
  translate 0.25*length * vnormalize(-vx)
}

// Let's now render the actual main objects in the scene.

// Retrieve mesh.
#include "separatrix3d.inc"
#declare nullcline_line_width = 0.001;
#declare intersection_line_width = 2*nullcline_line_width;
#declare streamline_width = 0.00125;
#declare arrow_width = 0.0075;

// Stagnation point.
#declare stagnation_point_radius = 0.015;
sphere {
  v0, stagnation_point_radius
  no_shadow
  pigment { colour White }
  finish { shinyFinish }
}

// Zero acceleration surface.
isosurface
{
  function { x*x + z }
  clipped_by { box { v0, vx+vy+vz } }
  no_shadow
  pigment { colour Grey transmit 0.7 }
  finish { glassFinish }
}

// Outline where zero-acceleration surface intersects separatrix and box boundaries

intersection
{
  union
  {
    onAxisNullcline(nullcline_line_width)
    object
    {
      onAxisNullcline(nullcline_line_width)
      translate vy
    }
  }
  separatrix

  // Only show part contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour Black }
  finish { matteFinish }
}

object
{
  zeroAccelerationIntersection(intersection_line_width)

  // Only show part contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour intersectionManifold }
  finish { matteFinish }
}

// xz-plane shows projection of streamlines onto the on-axis problem.
#declare backScreen = mesh2
{
  vertex_vectors { 4, v0, vx, vx+vz, vz }
  uv_vectors { 4, <0,1>, <1,1>, <1,0>, <0,0> }
  face_indices { 2, <0,1,2>, <2,3,0> }
  inside_vector  <0, -1, 0>
  uv_mapping

  no_shadow
  pigment { image_map { png "hybrid_streamlines_noaxis_eps=0.000" } }
  finish { matteFinish }
}

backScreen

// 3d separatrix surface and interior.

#declare collidingTrajectories = difference
{
  box { v0, vx+vy+vz }
  separatrix
}

difference
{
  object { collidingTrajectories }
  #if (DrawStreamLines)
    sGridLines(streamline_width, arrow_width)
  #end

  // Only show part contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour collidingManifold }
  finish { shinyFinish }
}

#if (DrawStreamLines)
  intersection
  {
    object { collidingTrajectories }
    sGridLines(streamline_width, arrow_width)

    // Only show part contained within our axes.
    bounded_by { box { v0, vx+vy+vz } }
    clipped_by { bounded_by }

    pigment { colour collidingParticle }
    finish { shinyFinish }
  }
#end
