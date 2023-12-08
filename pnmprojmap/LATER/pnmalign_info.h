/* Documentation of {pnmalign.c} */
/* Last edited on 2023-10-14 10:57:51 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -radius {XR} {YR} \\\n" \
  "  -maxDisp {XMAX} {YMAX} \\\n" \
  "  [ -maxTrials {MAX_TRIALS} ] \\\n" \
  "  " imgc_parse_haxis_HELP " \\\n" \
  "  " imgc_parse_vaxis_HELP " \\\n" \
  "    -image {PNMFILE[0]} [ -center ] {X[0]} {Y[0]} \\\n" \
  "  { -image {PNMFILE[1]} [ -center ] {X[1]} {Y[1]} }.. "

#define DEFAULT_MAX_TRIALS 10000

#define PROG_INFO_DESC \
  "  The program reads {N} images {IM[0..N-1]}, and takes" \
  " an initial reference point {P[i]} for" \
  " each image {IM[i]}.  It then computes a new reference" \
  " point {Q[i]} near each {P[i]}, so that the" \
  " neighborhoods of points {Q[0..N-1]} in images" \
  " {IM[0..N-1]} are maximally similar.  Then it" \
  " writes the new reference points {Q[0..N-1]}" \
  " to standard output, in the same format as the" \
  " \"-image\" command line arguments above.\n" \
  "\n" \
  "   Each neighborhood is a Hann-weighted window with" \
  " half-width {XR} and half-height {YR} centered at the" \
  " *original* reference point {P[i]}.  The neighborhoods are compared" \
  " with sample interpolation but without rotation.  The" \
  " similarity metric is {fimm_mismatch_var} applied separately" \
  " to each channel and summed over all channels.\n" \
  "\n" \
  "  " imgc_user_axes_intro_INFO "" \
  "  " imgc_pixel_axes_intro_INFO "" \
  "  " imgc_pixel_centers_intro_INFO "\n" \
  "\n" \
  "  The images are read from PNM image files" \
  " called \"{PNMFILE[0]}\", \"{PNMFILE[1]}\"."
  
#define PROG_INFO_OPTS \
  imgc_parse_haxis_HELP_INFO "" \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments, including" \
  " {XC} and {X[i,j]}." \
  "  " imgc_parse_haxis_pbm_default_INFO "\n" \
  "\n" \
  imgc_parse_vaxis_HELP_INFO "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments, including" \
  " {YC} and {Y[i,j]}." \
  "  " imgc_parse_vaxis_pbm_default_INFO "\n" \
  "\n" \
  "  -radius {XR} {YR}\n" \
  "    This required argument specifies the size of the neighborhood," \
  " around each alignment point, which is considered when comparing" \
  " the aligned images.  The neighborhood is a rectangle {2*XR}" \
  " pixels wide and {2*YR} pixels hight, centered at the alignment" \
  " point.\n" \
  "\n" \
  "  -maxDisp {XMAX} {YMAX}\n" \
  "    This required argument specifies the maximum adjustment" \
  " allowed on any reference point, in either direction, along" \
  " each axis.\n" \
  "\n" \
  "  -maxTrials {MAX_TRIALS}\n" \
  "    This optional argument specifies the maximum number of" \
  " possible alignments that will be tried. The default of" \
  " is " stringify(DEFAULT_MAX_TRIALS) ".\n" \
  "\n" \
  "  -image {PNMFILE[i]} [ -center ] {X[i]} {Y[i]}\n"  \
  "    Each instance of this argument specifies the name of the PBM/PGM/PPM file" \
  " that contains an additional input image {IM[i]}, and the coordinates" \
  " {(X[i],Y[i])} of the initial (approximate) reference point {P[i]}," \
  " in the image's domain.  The coordinate system is defined" \
  " by the \"-haxis\", \"-vaxis\", and \"-center\" arguments.  If this" \
  " argument is not given, the program does" \
  " nothing (since there are no images). \n" \
  "\n" \
  "  -center\n" \
  "    This keyword may appear immediately after the filename {PNMFILE[i]}" \
  " of an \"-image\" argument.  It sepcifies that the following" \
  " coordinates {X[i]} {Y[i]} are to be measured from the image's" \
  " center.  If omitted, the coordinates are measured from the" \
  " default origin (at the low-side edge). \n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics." \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmgtran(1), pnmscale(1), pnmrotate(1), pnmshear(1), pnm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created aug/2002 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  nov/2006 Rewritten to use sane PNM, argparser, etc. by J. Stolfi, IC-UNICAMP.\n" \
  "  aug/2007 Removed the guts to the library {float_image_align}.\n"               \
  "  aug/2007 Simplified so as to only adjust the alignment points.\n"               \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2002 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

