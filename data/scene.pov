#include "shapes.inc"
#include "colors.inc"

global_settings
{
  max_trace_level 256
  ambient_light White
}

#declare noncollidingManifold = rgb <1, 1, 0.753>;
#declare noncollidingParticle = rgb <0.937, 0.855, 0.302>;

#declare collidingManifold = rgb <0.855, 0.733, 0.89>;
#declare collidingParticle = rgb <0.8, 0.404, 0.761>;

background { colour Blue }

#declare shinyFinish = finish
{
  ambient 0.1
  diffuse 0.6
  phong 0.1
  phong_size 40
}

#declare matteFinish = finish
{
  ambient 0.1
  diffuse 0.6
  phong 0.
}

#declare v0 = <0 0 0>;
#declare vx = <1 0 0>;
#declare vy = <0 0.5 0>;
#declare vz = <0, 0, -0.5>;

#declare focusPos = 0.5 * (vx + vy + vz);
#declare cameraPos = focusPos + 0.9*<0 1 0.8>;

// // Uncomment to zoom in on one edge of box.
//#declare focusPos = vx+vy;
//#declare cameraPos = focusPos + <0.025, 0.025, 0.025>;

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
    //shadowless
  }
#end

//beam(focusPos + <0.0, 1.0, 0.0>, 1.0)
//beam(focusPos + <0.0, 0.0, 1.0>, 1.0)

//beam(focusPos + <1.0, 1.0, 1.0>, 1.0)
beam(cameraPos, 1.0)
beam(focusPos + <1.0, 0.25, 0.25>, 0.5)
beam(focusPos + <0.25, 1.0, 0.25>, 0.5)

// Frames for the domain boundaries.

#declare axis_radius = 0.005;
#macro axis(v0, v1)
  cylinder{v0, v1, axis_radius}
#end

union
{
  axis(v0, vx)
  axis(v0, vy)
  axis(v0, vz)

  // Axes should be in matte black.
  pigment { colour Black }
  finish { matteFinish }
}

#declare box_edge_radius = 0.005;
#macro edge(v0, v1)
  union
  {
    cylinder{v0, v1, box_edge_radius}
    sphere{v0, box_edge_radius}
    sphere{v1, box_edge_radius}
  }
#end

union
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

  // We apply some non-physical properties to these edges so the viewer unconsciously disregards
  // the edges: their presence acts as a guide to the eye to help get some 3d perspective when
  // looking at the image projected in 2d.
  // We first make the non-axis edges highly transparent, then prevent them from casting any
  // shadows on the scene to prevent optical artifacts from distracting the viewer.
  pigment { colour White transmit 0.9 }
  interior { ior 2.5 } // index of refraction of diamond gives an ethereal-look
  no_shadow

  finish { matteFinish }
}

// Stagnation point.

#declare stagnation_point_radius = 0.02;
sphere {
  v0, stagnation_point_radius
  no_shadow
  texture
  {
    pigment { colour White }
    finish { matteFinish }
  }
}

// xz-plane shows projection of streamlines there.

mesh2 {
  vertex_vectors {
    4, v0, vx, vx+vz, vz
  }
  face_indices {
    2, <0,1,2>, <2,3,0>
  }
  no_shadow
  pigment { colour noncollidingParticle }
  finish { matteFinish }
}

// 3d separatrix surface.

#include "separatrix3d.inc"

mesh2 {
  separatrix
  // Only show part of separatrix contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour collidingParticle }
  finish { shinyFinish }
}
