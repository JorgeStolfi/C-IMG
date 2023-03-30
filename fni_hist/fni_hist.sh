#! /bin/bash 
# Last edited on 2012-12-08 21:39:12 by stolfilocal
 
PROG_NAME=${0##*/}
PROG_DESC="make a Postscript histogram of a selected channel of a FNI file"
PROG_HELP=(
  "${PROG_NAME} \\"
  "\n    [ -channel {CHNUM} ] \\"
  "\n    [ -range {VMIN} {VMAX} ] \\"
  "\n    [ -title {TITLESTRING} ] \\"
  "\n    -step {STEP} \\"
  "\n    < {INFILE}.fni > {OUTFILE}.eps"
)
PROG_INFO=(
  "\nNAME"
  "\n  ${PROG_NAME} - ${PROG_DESC}."
  "\n"
  "\nSYNOPSIS"
  "\n  ${PROG_HELP[@]}"
  "\n"
  "\nDESCRIPTION"
  "\n  Reads a multichannel float-valued image (\".fni\") file,"
  "\n  generates an Encapsulated Postscript histogram of a"
  "\n  selected channel using {gnuplot}."
  "\n"
  "\n  The histogram has bins of width {STEP}."
  "\n"
  "\n  If the \"-range\" option is given, the program ignores"
  "\n  any data values outside the range [{VMIN} _ {VMAX}]."
  "\n"
  "\nAUTHOR"
  "\n  Created 2006-04-02 by Jorge Stolfi, Unicamp"
)

# ----------------------------------------------------------------------
# INTERNAL OPTIONS

# ----------------------------------------------------------------------
# COMMAND LINE PARSING

# Parse command line switches: 
channel=0
vmin="+1.0"
vmax="-1.0"
step="0.0"
title=""
while [[ ( $# -ge 1 ) && ( "/$1" =~ /-.* ) ]]; do
  if [[ ( $# -ge 2 ) && ( $1 == "-channel" ) ]]; then 
    channel="$2"; shift; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-step" ) ]]; then 
    step="$2"; shift; shift;
  elif [[ ( $# -ge 3 ) && ( $1 == "-range" ) ]]; then 
    vmin="$2"; vmax="$3"; shift; shift; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-title" ) ]]; then 
    title="$2"; shift; shift;
  elif [[ ( $# -ge 1 ) && ( ( "/$1" == "/-help" ) || ( "/$1" == "/--help" ) ) ]]; then 
    echo -e "usage:\n  ${PROG_HELP[@]}"; exit 0;
  elif [[ ( $# -ge 1 ) && ( ( "/$1" == "/-info" ) || ( "/$1" == "/--info" ) ) ]]; then 
    echo -e "${PROG_INFO[@]}"; exit 0;
  else
    echo "unknown option $1" 1>&2 ;
    echo -e "usage:\n  ${PROG_HELP[@]}" 1>&2 ; exit 1 
  fi
done 

if [[ "/${title}" == "/" ]]; then
  title="Channel ${channel}"
fi

# Get positional parameters

# Check for leftover arguments:
if [[ $# -gt 0 ]]; then
  echo 'excess arguments "'"$1"'" ...' 1>&2 ;
  echo -e "usage:\n  ${PROG_HELP[@]}" 1>&2 ; exit 1 
fi

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

# Prefix for temporary file names
tmp="/tmp/$$"

# Save input to a temporary file:
fnifile=${tmp}.fni
cat > ${fnifile}

# Extract image size:
size=( `fni_size ${tmp}.fni` )
nx="${size[1]}"; shift;
ny="${size[2]}"; shift;

# Generate histogram:
hisfile=${tmp}.his

gawk \
  -v chan="${channel}" \
  -v step="${step}" \
  -v vmin="${vmin}" \
  -v vmax="${vmax}" \
  ' BEGIN {
      abort = -1;
      chan += 0; step += 0; vmin += 0; vmax += 0;
      if (step == 0) { arg_error("must specify a numeric \"-step\""); }
      split("", ct);
      kmin = +1.0e+100; kmax = -1.0e+100; nout = 0;
      koff = int(vmin/step);
      voff = koff*step;
      maxbins = 10001;
      NC = 0;
    }
    (abort >= 0) { exit abort; }
    /[=]/ { gsub(/[=]/, " = ", $0); }
    /^NC[ ]*[=]/ { NC = $3; }
    /^[ ]*[-+]?[0-9]/ { 
      z = ((chan < 0) || (chan >= NC) ? 0.0 : $(3+chan));
      if ((vmin < vmax) && ((z < vmin) || (z > vmax))) { nout++; next; }
      zz = (z - voff)/step; 
      k = int(zz);
      if (k < kmin) { kmin = k; }
      if (k > kmax) { kmax = k; }
      if (kmax - kmin + 1 > maxbins) { data_error("too many bins"); }
      if (! (k in ct)) { ct[k] = 0; }
      ct[k]++;
      next;
    }
    END { 
      if (abort >= 0) { exit abort; }
      if (nout > 0) 
        { printf "%d values outside the given range\n", nout > "/dev/stderr"; }
      if (kmin > kmax) { kmin = 0; kmax = 0; }
      for (k = kmin-1;  k <= kmax+1; k++)
        { md = voff + k*step; lo = md - step/2; hi = md + step/2;
          printf " %+14.7f %+14.7f %+14.7f %10d\n", lo, md, hi, ct[k]; 
        }
    }
    function arg_error(msg)
      { printf "** %s\n", msg > "/dev/stderr"; abort = 1; exit abort; }
    function data_error(msg)
      { printf "** line %d: %s\n", FNR, msg > "/dev/stderr"; abort = 1; exit abort; }
  ' \
  ${fnifile} \
  > ${hisfile}

# Get maxmum bin count:
maxfile=${tmp}.max
gawk \
  ' BEGIN {m=0;} 
    //{ if ($4 > m) { m = $4; } }
    END { print m; }
  ' \
  < ${hisfile} \
  > ${maxfile}
  
ymax="`cat ${maxfile}`"

gnuplot <<EOF
set terminal postscript eps color 
set output
set size 1.50,0.75
set nokey
ymax=${ymax}
set yrange[(-0.02*ymax):(1.02*ymax)]
set title "${title}"
plot "${hisfile}" using 2:4 with histeps
quit
EOF

# /bin/rm -f ${fnifile} ${hisfile}
