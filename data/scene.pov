#include "shapes.inc"
#include "colors.inc"
#include "functions.inc"

global_settings
{
  max_trace_level 256
  ambient_light White
}

#declare noncollidingManifold = rgb <1, 1, 0.753>;
#declare noncollidingParticle = rgb <0.937, 0.855, 0.302>;

#declare collidingManifold = rgb <0.855, 0.733, 0.89>;
#declare collidingParticle = rgb <0.8, 0.404, 0.761>;

//background { colour White }

#declare shinyFinish = finish
{
  ambient 0.1
  diffuse 0.6
  phong 0.1
  phong_size 40
  //brilliance 0.5
}

#declare matteFinish = finish
{
  ambient 0.1
  diffuse 0.6
  phong 0.
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
  }
#end

//beam(focusPos + <0.0, 1.0, 0.0>, 1.0)
//beam(focusPos + <0.0, 0.0, 1.0>, 1.0)

//beam(focusPos + <1.0, 1.0, 1.0>, 1.0)
beam(cameraPos, 1.0)
beam(focusPos + <1.0, 0.25, 0.25>, 0.5)
beam(focusPos + <0.25, 1.0, 0.25>, 0.5)

// Frames for the domain boundaries.

// #declare axis_radius = 0.005;
// #macro axis(v0, v1)
//   cylinder{v0, v1, axis_radius}
// #end

#declare box_edge_radius = 0.005;
#macro edge(v0, v1)
  cylinder { v0, v1, box_edge_radius }
#end
#macro node(v1)
  sphere { v1, box_edge_radius }
#end

// merge
// {
//   edge(v0, vx)
//   edge(v0, vy)
//   edge(v0, vz)
//   node(vx)
//   node(vy)
//   node(vz)

//   // Axes should be in matte black.
//   pigment { colour Black }
//   finish { matteFinish }
// }

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

// Stagnation point.

#declare stagnation_point_radius = 0.015;
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
  vertex_vectors { 4, v0, vx, vx+vz, vz }
  uv_vectors { 4, <0,1>, <1,1>, <1,0>, <0,0> }
  face_indices { 2, <0,1,2>, <2,3,0> }

  uv_mapping
  no_shadow
  pigment { image_map { png "backdrop" } }
  finish { shinyFinish }
}

// Let's now render the actual main objects in the scene.

// Retrieve mesh.
#include "separatrix3d.inc"
#declare line_width = 0.0025;

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

// Projection of zero-acceleration surface onto y=0.
object
{
  zeroAcceleration(line_width)
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour Black }
  finish { matteFinish }
}

// Streamlines
union
{
  line0(line_width)
  line1(line_width)
  line2(line_width)
  line3(line_width)
  line4(line_width)
  line5(line_width)
  line6(line_width)
  line7(line_width)
  line8(line_width)
  line9(line_width)
  line10(line_width)
  line11(line_width)
  line12(line_width)
  line13(line_width)
  line14(line_width)
  line15(line_width)
  line16(line_width)
  line17(line_width)
  line18(line_width)
  line19(line_width)
  line20(line_width)
  line21(line_width)
  line22(line_width)
  line23(line_width)
  line24(line_width)
  line25(line_width)
  line26(line_width)
  line27(line_width)
  line28(line_width)
  line29(line_width)
  line30(line_width)
  line31(line_width)
  line32(line_width)
  line33(line_width)
  line34(line_width)
  line35(line_width)
  line36(line_width)
  line37(line_width)
  line38(line_width)
  line39(line_width)
  line40(line_width)
  line41(line_width)
  line42(line_width)
  line43(line_width)
  line44(line_width)
  line45(line_width)
  line46(line_width)
  line47(line_width)
  line48(line_width)
  line49(line_width)
  line50(line_width)
  line51(line_width)
  line52(line_width)
  line53(line_width)
  line54(line_width)
  line55(line_width)
  line56(line_width)
  line57(line_width)
  line58(line_width)
  line59(line_width)
  line60(line_width)
  line61(line_width)
  line62(line_width)
  line63(line_width)
  line64(line_width)
  line65(line_width)
  line66(line_width)
  line67(line_width)
  line68(line_width)
  line69(line_width)
  line70(line_width)
  line71(line_width)
  line72(line_width)
  line73(line_width)
  line74(line_width)
  line75(line_width)
  line76(line_width)
  line77(line_width)
  line78(line_width)
  line79(line_width)
  line80(line_width)
  line81(line_width)
  line82(line_width)
  line83(line_width)
  line84(line_width)
  line85(line_width)
  line86(line_width)
  line87(line_width)
  line88(line_width)
  line89(line_width)
  line90(line_width)
  line91(line_width)
  line92(line_width)
  line93(line_width)
  line94(line_width)
  line95(line_width)
  line96(line_width)
  line97(line_width)
  line98(line_width)
  line99(line_width)

  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour collidingParticle }
  finish { matteFinish }
}

// 3d separatrix surface.
mesh2
{
  separatrix
  // Only show part of separatrix contained within our axes.
  bounded_by { box { v0, vx+vy+vz } }
  clipped_by { bounded_by }

  pigment { colour collidingParticle }
  finish { shinyFinish }
}
