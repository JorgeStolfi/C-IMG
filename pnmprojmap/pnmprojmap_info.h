/* Manpage for {pnmprojmap.c}. */
#ifndef pnmprojmap_info_H
#define pnmprojmap_info_H

/* Last edited on 2024-10-31 10:20:50 by stolfi */

#include <image_coords.h>
#include <sample_conv.h>
#include <float_image_transform.h>
#include <argparser.h>

#define img_ops_HELP \
  "    " imgc_parse_center_org_HELP("center","org","In[k]") " \\\n" \
  "        " imgc_parse_unit_HELP("unit","In[k]") " \\\n" \
  "        [ matrix {mapName[k]} ] \\\n" \
  "        [ points {ptsName[k]} ]"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " imgc_parse_x_axis_HELP " \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
  "    -inPrefix {inPrefix} \\\n" \
  "    { -image {imgName[k]} {ext[k]} \\\n" \
  "    " img_ops_HELP " \\\n" \
  "    }.. \\\n" \
  "    [ -fromMatrix | -fromPoints ] \\\n" \
  "    " imgc_parse_unit_HELP("-oUnit","Out") " \\\n" \
  "    " imgc_parse_size_HELP("-oSize","Out") " \\\n" \
  "    " imgc_parse_center_org_HELP("-oCenter","-oOrg","Out") " \\\n" \
  "    [ -interpolate {intOrder} ] \\\n" \
  "    [ -defaultVal {defaultVal} | -extend ] \\\n" \
  "    [ -maxval {maxvalOut} ] \\\n" \
  "    [ -isMask {isMask} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -noImages ] \\\n" \
  "    [ -debug {debugCol} {debugRow} ] \\\n" \
  "    -outPrefix {outPrefix}"

#define PROG_INFO_COORDS_INTRO \
  "  " imgc_user_axes_intro_INFO "\n" \
  "\n" \
  "  The origin of each user coordinate system is" \
  " specified individually for each input image through" \
  " the \"center\" and \"org\" image attributes, and globally" \
  " for all output images through the \"-oCenter\" and \"-oOrg\" options.  Similarly," \
  " the unit of length in each user coordinate system, input or output," \
  " is a /user unit/ that is multiple or fraction of a pixel, specified by" \
  " the \"unit\" image attribute or the \"-oUnit\" option. See details in" \
  " the OPTIONS section below.  The coordinate system parameters affect" \
  " the interpretation of the projective map specification, whether given" \
  " as pairs of points or as a projective matrix.\n" \
  "\n" \
  "  " imgc_pixel_axes_intro_INFO "" \
  "  " imgc_pixel_centers_intro_INFO

#define PROG_INFO_DESC \
  "  The main function of the program is to read one or more image" \
  " files \"{inPrefix}{imgName[k]}.{ext[k]}\", for {k=0,1,...} and apply to" \
  " each of them a specified projective transformation {M[k]}; so that the" \
  " pixel value at a point {p} of the input image is copied to" \
  " point {M[k](p)} the output image. Then it writes that" \
  " output image to file \"{outPrefix}{imgName[k]}\".\n" \
  "\n" \
  "  The maps {M[k]} can may be given explicitly or determined" \
  " indirectly by matching lists of points.  In any case, the maps can be " \
  " writen out in a format that can be used by subsequent runs of this program."

#define xaxis_def imgc_parse_x_axis_INFO_OPTS_default_pbm
#define yaxis_def imgc_parse_y_axis_INFO_OPTS_default_pbm
#define iorg_def imgc_parse_center_org_INFO_OPTS_default_zero("org")
#define osize_def \
  "If not specified, the size of each output image, in output user units, will be" \
  " the same as that of the corresponding input image, in its input user units."
#define oorg_def \
  "If not specified, the origin of the output coordinate system" \
  " will be the center of the image."

#define PROG_INFO_MAP_FILES \
  "  If the \"-fromMatrix\" option is given (see below), the projective" \
  " transformation {M[k]} that will be applied to input image" \
  " number {k} is read from a" \
  " file specified with the \"matrix\" image" \
  " attribute (see the \"-image\" option below).\n" \
  "\n" \
  "  The file must contain three lines, each with three numbers, being the elements of" \
  " the matrix {M[k].dir} as described in {hr2.h}, one row per line.  A '#' character" \
  " indicates that the rest of the lines is a comment, which is" \
  " ignored.  Lines that are blank or contain only '#'-comments are" \
  " allowed only before the first data line.  Each line, including the last one, must be" \
  " properly terminated by an end-of-line (ASCII NL) character."

#define PROG_INFO_POINTS_FILES \
  "  If the \"-fromPoints\" option is given (see below), the projective" \
  " transformation {M[k]} to apply to each input image is computed by comparing" \
  " a list of /key feature points/ associated to image {k} with those" \
  " associated to the first image (index {k=0}). The key feature points are read" \
  " from a \".pts\" file specified by the \"points\" image" \
  " attribute (see the \"-image\" option below).\n" \
  "\n" \
  "  Blank lines and #-comments in this file are ignored.  Each data line {j} of" \
  " this file has one of the three forms below:\n" \
  "\n" \
  "      {tag[j]} {XF[k,j]} {YF[k,j]} P\n" \
  "      {tag[j]} {XF[k,j]} {YF[k,j]} C {rad[k,j]} {ang[k,j]}\n" \
  "      {tag[j]} {XF[k,j]} {YF[k,j]} E {XU[k,j]} {YU[k,j]}  {XV[k,j]} {YV[k,j]}\n" \
  "\n" \
  "    where {tag[j]} is an identifier (that starts with a letter" \
  " and contains no blanks), the third column" \
  " is a capital letter as shown, and the other fields are" \
  " numbers, possibly with a decimal period.\n" \
  "\n" \
  "  In this case, the command line must specify at least" \
  " two images.  The map {M[0]} will" \
  " always be the identity, apart from scaling factor determined by the" \
  " ratio of input and output user units.  For every other {k}, the program" \
  " will compute the map {M[k]} so that it takes each" \
  " point {P[j]=(XF[k,j],YF[k,j])} to" \
  " a point as close as possible to {Q[j]=(XF[0,j],YF[0,j])}, in the least mean" \
  " square error sense.  For this purpose, the program will consider only lines" \
  " of the two \".pts\" files with the same {tag}, and will implicitly sort" \
  " and renumber them from 0; so that, effectively, {tag[k,j] = tag[0,j]} for" \
  " all {j}.\n" \
  "\n" \
  "  The map will take each point {P[j]} precisely to {Q[j]} if the number" \
  " of matched pairs does not exceed a limit that depends on the specified" \
  " map type.  The limit is 1 pair for a translation, 2 pairs for a" \
  " similarity, 3 for an affine map, and 4 for a general" \
  " transformation.  If the number of matching pairs exceed this limit, or" \
  " the specified mapping type is a congruence, each" \
  " point {P[j]} will usually /not/ match the corresponding {Q[j]}; it will" \
  " only be as close as possible to it.\n" \
  "\n" \
  "  The other fields in the \".pts\" data files are currently ignored." \

#define PROG_INFO_OPTS \
  imgc_parse_x_axis_INFO_OPTS(xaxis_def) \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments.\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(yaxis_def) "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments.\n" \
  "\n" \
  "  -inPrefix {inPrefix}\n" \
  "    This optional argument defines the common prefix for all input" \
  " file names.  If omitted, assumes \"\" (that is, read files from" \
  " the current directory).\n" \
  "\n" \
  "  -image {imgName[k]} {ext[k]} [ {ATTRIBUTE}.. ]\n" \
  "    Each occurrence of this keyword specifies the file" \
  " name and extension for an input image, and optionally zero" \
  " or more attributes that" \
  " affect the handling of that image, including the interpretation" \
  " of coordinates on its domain.  The image" \
  " file name will be \"{inPrefix}{imgName[k]}.{ext}\". Currenty the extension" \
  " must be \"ppm\", \"pgm\", or \"pbm\".  See also" \
  " the \"-fromPoints\" and \"-fromMatrix\" options below.\n" \
  "\n" \
  "    Each image {ATTRIBUTE} may be any or all of the following, in any order:\n" \
  "\n" \
  "      unit {unitIn[k]}\n" \
  "        This optional attribute specifies the size in pixels of the unit for" \
  " user coordinate on image {imgName[k]}.  If this argument is not specified, the" \
  " program assumes \"unit 1\" (that is, user units are pixels).\n" \
  "\n" \
  "      center\n" \
  "      org {CXIn[k]} {CYIn[k]}\n" \
  "        These mutually exclusive optional attributes specify the origin of the" \
  " user coordinate system for the image {imgName[k]}.  The \"center\" option sets the" \
  " origin at the center of the image domain.  The \"org\" option sets the origin" \
  " displaced {CXIn[k]} and {CYIn[k]} /user/ units in the directions of" \
  " the X and Y axes, from the corner of the image  opposite to those" \
  " directions." iorg_def "  " \
  ""  imgc_unit_affects_org_INFO_OPTS("unit","org","the input image") "\n" \
  "\n" \
  "      matrix {mapName[k]}\n" \
  "        This optional attribute specifies that, if the \"-fromMatrix\" option" \
  " is given, the projective map {M[k]} to apply to this image is to be" \
  " read from file \"{inPrefix}{mapName[k]}.map\".  See the" \
  " section PROJECTIVE MAP MATRIX FILES for the format of this" \
  " file.  If this attribute is not specified, {mapName[k]} defaults" \
  " to be the same as {imgName[k]}.\n" \
  "\n" \
  "      points {ptsName[k]}\n" \
  "        This optional attribute specifies that, if the \"-fromPoints\" option" \
  " is given, the list of key feature points for this image is to be read" \
  " from the file \"{inPrefix}{ptsName[k]}.pts\".  See the" \
  " section KEY FEATURE POINTS FILES for the format of this" \
  " file.  If this attribute is not specified, {ptsName[k]} defaults" \
  " to be the same as {imgName[k]}.\n" \
  "\n" \
  "      type {mapType[k]}\n" \
  "        This optional attribute specifies the type of the projective map {M[k]} used to transform this image.  It is either \".\n" \
  "\n" \
  "  -fromMatrix\n" \
  "  -fromPoints\n" \
  "    Exactly one of these two options must be present. They specify" \
  " how the projective map {M[k]} for each input" \
  " image is determined.\n" \
  "\n" \
  "    The option \"-fromMatrix\" specifies that {M[k]} is" \
  " read from the file \"{inPrefix}{mapName[k]}.map\" (see the \"matrix\" attribute" \
  " of \"-image\").\n" \
  "\n" \
  "    The option \"-fromPoints\" specifies instead that the map {M[k]} is to" \
  " be computed from a list of key feature points that is to be" \
  " read from file \"{inPrefix}{ptsName[k]}.pts\" (see the \"points\" attribute" \
  " of \"-image\"), which shall be matched with the list of key features the" \
  " first image.\n" \
  "\n" \
  imgc_parse_unit_INFO_OPTS("-oUnit","Out","every output image") "   If not" \
  " specified, the output user unit for each image will be the same as" \
  " its input user unit.\n" \
  "\n" \
  "    " imgc_unit_affects_org_INFO_OPTS("-oUnit","org","the corresponding image") "\n" \
  "\n" \
  imgc_parse_size_INFO_OPTS("-oSize","Out","all output images",osize_def) "\n" \
  "\n" \
  imgc_parse_center_org_INFO_OPTS("-oCenter","-oOrg","Out",oorg_def) "\n" \
  "\n" \
  "    " imgc_unit_affects_org_INFO_OPTS("-oUnit","-oOrg","the output image") "\n" \
  "\n" \
  "  -interpolate {intOrder} \n" \
  "    This optional argument specifies the way inpt pixels are" \
  " interpolated.  " float_image_transform_interpolation_HELP_INFO "  The default" \
  " is \"-interpolate 0\" (C0 bilinear interpolation).\n" \
  "\n" \
  "  -defaultVal {defaultVal}\n" \
  "  -extend \n" \
  "    These mutually exclusive optional flags specify the handling" \
  " of source pixels that fall outside the input image's" \
  " domain.  If \"-defaultVal\" is used, any such pixel is assumed to" \
  " have value {defaultVal}, in a scale from 0 to 1.  If \"-extend\" is" \
  " used, the input image" \
  " will be implicitly extended to an infinite image, before" \
  " being transformed, by replicating the pixels" \
  " along its borders.  If neither option is specified," \
  " the program assumes \"-defaultVal 0.5\".\n" \
  "\n" \
  "  -maxval {maxvalOut}\n" \
  "    Specifies {maxvalOut} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to to the" \
  " input image's {maxval}, but at least 255 for color images" \
  " and 65535 for grayscale ones.\n" \
  "\n" \
  "  -isMask {isMask}\n" \
  "    This optional Boolean argument modifies the interpretation" \
  " of integer sample values, especially 0 and {MAXVAL}, in the" \
  " input and output files (see {sample_conv.h}).  If {isMask} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If {isMask} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMask F\".\n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics.\n" \
  "\n" \
  "  -noImages\n" \
  "    If this option is present, the program does not actually compute" \
  " nor write the output images.  It still reads or computes the perspective" \
  " map, and writes the merged feature point list.  ??? It still reads the input" \
  " images, in order to determine the input user-to-pixel mappings, and stil requires" \
  " the \"-oSize\" parameters.\n" \
  "\n" \
  "  -outPrefix {outPrefix}\n" \
  "    This mandatory argument defines the common prefix for all output file names.\n" \
  "\n" \
  "  -debug {debugCol} {debugRow}\n" \
  "    If this option is present, the program prints out" \
  " debugging information about the computation of the" \
  " output pixel on column {debugCol} and row {debugRow}.  This option ignores" \
  " all parameters defining the user coordinate systems and the" \
  " input-to-output map."

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
  "USER COORDINATE SYSTEMS\n" \
  PROG_INFO_COORDS_INTRO "\n" \
  "\n" \
  "PROJECTIVE MAP MATRIX FILES\n" \
  PROG_INFO_MAP_FILES "\n" \
  "\n" \
  "KEY FEATURE POINTS FILES\n" \
  PROG_INFO_POINTS_FILES "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "BUGS\n" \
  "  Sampling is not very scientific; it may blur more" \
  " than necessary, and may not work properly if the map is too warped.\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmscale(1), pnmradist(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created aug/2002 by Jorge Stolfi, IC-UNICAMP as \"pnmgtran\".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  All changes by J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "  nov/2006 Rewritten to use sane PNM, argparser, etc.\n" \
  "  jun/2007 Added \"-points\" option.\n" \
  "  jul/2007 Added \"-xAxis\", \"-yAxis\" options.\n" \
  "  aug/2007 Rearranged the points in \"-points\".\n" \
  "  jan/2008 Moved radial distortion to \"pnmradist\".\n" \
  "  jan/2008 Renamed from \"pnmgtran\" to \"pnmprojmap\".\n" \
  "  jan/2008 Added the \"-extend\" option.\n" \
  "  aug/2010 Added the \"-interpolate\" option.\n" \
  "  aug/2010 Added the \"-isMask\" option.\n" \
  "  mar/2017 Moved matrix help and parsing to {argparser_geo.h}\n" \
  "  mar/2023 Converted {int} to {int32_t}.\n" \
  "  aug/2023 Added the \"-scale\" option.\n" \
  "  aug/2023 Added the \"-mapFile\", \"-map\", and \"-unmap\" options.\n" \
  "  aug/2023 Replaced the \"-scale\" by \"-iUnit\" and \"-oUnit\".\n" \
  "  aug/2023 Updated for changes in {image_input_output_coords.h}.\n" \
  "  sep/2023 Added the \"-noImage\" option.\n" \
  "  oct/2023 Radical reform of options to prepare for multi-image processing.\n" \
  "  aug/2024 Bug fixes, got it to work again.\n" \
  "  aug/2024 Map type option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2006 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#endif
