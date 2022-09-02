#version 3.7;
global_settings { assumed_gamma 2.2 }

#include "shapes.inc"
#include "colors.inc"
#include "functions.inc"

global_settings
{
  max_trace_level 256
  ambient_light White
}

#declare noncollidingManifold = srgb <1, 1, 0.753>;         // #ffffc0
#declare noncollidingParticle = srgb <0.937, 0.855, 0.302>; // #efda4d

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
  sky      <0.0, 0.0, 1.0>
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
  beam(focusPos + <1.0, 0.25, 0.25>, 0.5)
  beam(focusPos + <0.25, 1.0, 0.25>, 0.5)
#end

// Frames for the domain boundaries.

#declare box_edge_radius = 0.005;
#macro edge(v0, v1)
  cylinder { v0, v1, box_edge_radius }
#end
#macro node(v1)
  sphere { v1, box_edge_radius }
#end

merge
{
  edge(v0, vx)
  edge(v0, vy)
  edge(v0, vz)

  edge(vx, vx+vy)
  edge(vy, vx+vy)
  edge(vx, vx+vz)
  edge(vz, vx+vz)
  edge(vy, vy+vz)
  edge(vz, vy+vz)

  edge(vx+vy, vx+vy+vz)
  edge(vx+vz, vx+vy+vz)
  edge(vy+vz, vx+vy+vz)

  // Cap the ends with spheres so the edges join smoothly
  node(v0)
  node(vx)
  node(vy)
  node(vz)
  node(vx+vy)
  node(vx+vz)
  node(vy+vz)
  node(vx+vy+vz)

  // We apply some non-physical properties to these edges so the viewer unconsciously disregards
  // the edges: their presence acts as a guide to the eye to help get some 3d perspective when
  // looking at the image projected in 2d.
  // We first make the non-axis edges highly transparent, then prevent them from casting any
  // shadows on the scene to prevent optical artifacts from distracting the viewer.
  pigment { colour White transmit 0.9 }
  interior { ior 2.5 } // index of refraction of diamond gives an ethereal-look
  no_shadow
  finish { glassFinish }
}

// Let's now render the actual main objects in the scene.

// Retrieve mesh.
#include "separatrix3d.inc"
#declare line_width = 0.00125;

// Stagnation point.
#declare stagnation_point_radius = 0.015;
sphere {
  v0, stagnation_point_radius
  no_shadow
  texture
  {
    pigment { colour White }
    finish { ambient 1 }
  }
}

// Zero acceleration surface.
isosurface
{
  function { x*x + z }
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour Grey transmit 0.7 }
  no_shadow
  finish { glassFinish }
}

// Outline where zero-acceleration surface intersects separatrix and box boundaries
#declare zeroAccelerationLine = merge
{
  onAxisNullcline(line_width)
  zeroAccelerationIntersection(line_width)
  object
  {
    onAxisNullcline(line_width)
    translate vy
  }

  // Only show part contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { color Black }
  finish { matteFinish }
}

zeroAccelerationLine

// xz-plane shows projection of streamlines onto the on-axis problem.
#declare backScreen = mesh2
{
  vertex_vectors { 4, v0, vx, vx+vz, vz }
  uv_vectors { 4, <0,1>, <1,1>, <1,0>, <0,0> }
  face_indices { 2, <0,1,2>, <2,3,0> }

  no_shadow
  uv_mapping
  pigment { image_map { png "hybrid_streamlines_noaxis_eps=0.000" } }
  finish { matteFinish }
}

backScreen

// 3d separatrix surface.
difference
{
  object { separatrix }
  #if (DrawStreamLines)
    sGridLines(line_width)
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
    object { separatrix }
    sGridLines(line_width)

    // Only show part contained within our axes.
    bounded_by { box { v0, vx+vy+vz } }
    clipped_by { bounded_by }

    pigment { colour collidingParticle }
    finish { shinyFinish }
  }
#end
