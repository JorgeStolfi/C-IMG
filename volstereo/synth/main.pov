// Last edited on 2003-04-17 12:44:28 by stolfi
// Imitation of the m21 fragment image

// global_settings { assumed_gamma 2.2 }

#declare bgcolor = <0.306, 0.261, 0.250>;
background { color rgb bgcolor }

// All dimensions are in millimeters unless said otherwise.
// All angles are in degrees unless said otherwise.

// Fudge factor to avoid coincident surfaces:
#declare eps = 0.01;

#include "frag-stand.inc"
#include "frag-camera-whole.inc"
#include "frag-lights.inc"
#include "frag-monolith.inc"
#include "frag-frag.inc"

// Z coordinate of mirror surface
#declare mirroring_Z = mirror_min.z;

// Center of monolith bottom
#declare mono_1_tilt = 0.0;
#declare mono_1_turn = 0.0;
#declare mono_1_pos = < +050.1, +225.0, 0000.0 > + (frame_max.z + eps) * z;

#declare mono_2_tilt = 0.0;
#declare mono_2_turn = -3.0;
#declare mono_2_pos = < +130.1, +225.0, 0000.0 > + (frame_max.z + eps) * z;

// Center of fragment bottom
#declare fragment_pos = < +047.0, +190.0, mirror_max.z >;

#declare scene = 
  union {
    object { stand }
    object { monolith 
      rotate mono_1_tilt * y 
      rotate mono_1_turn * z 
      translate mono_1_pos
    }
    object { monolith 
      rotate mono_2_tilt * y 
      rotate mono_2_turn * z 
      translate mono_2_pos
    }
    object { lights }
    // Must add the virtual lights explicitly, since POV-Ray 
    //   does not reflect light off the mirror.   Note  that
    //   the mirror must be declared "no_shadow" for 
    //   this to work.
    object { lights 
      translate -mirroring_Z * z 
      scale <1,1,-1> 
      translate mirroring_Z * z 
    }
    object { fragment translate fragment_pos }
  }

object { scene }
