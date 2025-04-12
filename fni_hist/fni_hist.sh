#! /bin/bash 
# Last edited on 2025-04-09 15:11:09 by stolfi
 
PROG_NAME=${0##*/}
PROG_DESC="make a Postscript histogram of a selected channel of a FNI file"
PROG_HELP=(
  "${PROG_NAME} \\"
  "\n    [ -channel {CHNUM} ] \\"
  "\n    [ -range {VMIN} {VMAX} ] \\"
  "\n    [ -excludeRange {EMIN} {EMAX} ] \\"
  "\n    [ -logScale ] \\"
  "\n    [ -binRound { -1 | 0 | +1 } ] \\"
  "\n    [ -title {TITLESTRING} ] \\"
  "\n    [ -step {STEP} | -bins {BINS} ] \\"
  "\n    [ -verbose ] \\"
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
  "\n  If the \"-excludeRange\" option is given, the program ignores"
  "\n  any data values inside the range [{EMIN} _ {EMAX}]."
  "\n" 
  "\n The \"-binRound\" option determines how to handle samples that fall"
  "\n exactly on a bin boundary.  With \"-binRound -1\","
  "\n those values will be counted in the next lower bin."
  "\n With \"-binRround +1\", they will be counted in the next higher"
  "\n  bin. With \"-binRound 0\" (the default),"
  "\n each value will be counted as {0.5} in both bins."
  "\n" 
  "\n  If the \"-logScale\" option is given, the"
  "\n  plot uses log scale in the {Y} axis."
  "\n"
  "\n  If \"-step {STEP}\" is specified, each bin will have width {STEP}."
  "\n  If \"-bins {BINS}\" is specified instead, the width of each bin"
  "\n  will be a roundish number such that the total number"
  "\n  of bins as small as possible but no less"
  "\n  than {BINS}.  The default is \"-bins 200\"."
  "\n"
  "\nAUTHOR"
  "\n  Created 2006-04-02 by Jorge Stolfi, Unicamp"
  "\nMODIFICATION HISTORY"
  "\n  By J.Stolfi if not said otherwise"
  "\n  2025-02-12 Added \"-excludeRange\" and \"-logScale\""
  "\n  2025-04-08 Added \"-bins\""
  "\n  2025-04-09 Made \"-bins 200\" the default."
  "\n  2025-04-09 Added \"-verbose\"."
)

# ----------------------------------------------------------------------
# INTERNAL OPTIONS

# ----------------------------------------------------------------------
# COMMAND LINE PARSING

# Parse command line switches: 
channel=0
vmin="+99999999.0"
vmax="-99999999.0"
emin="+1.0"
emax="-1.0"
step="0"
bins="0"
binRound=0
logScale=0
verbose=0
title=""
while [[ ( $# -ge 1 ) && ( "/$1" =~ /-.* ) ]]; do
  if [[ ( $# -ge 2 ) && ( $1 == "-channel" ) ]]; then 
    channel="$2"; shift; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-step" ) ]]; then 
    step="$2"; shift; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-bins" ) ]]; then 
    bins="$2"; shift; shift;
  elif [[ ( $# -ge 3 ) && ( $1 == "-range" ) ]]; then 
    vmin="$2"; vmax="$3"; shift; shift; shift;
  elif [[ ( $# -ge 3 ) && ( $1 == "-excludeRange" ) ]]; then 
    emin="$2"; emax="$3"; shift; shift; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-binRound" ) ]]; then 
    binRound="$2"; shift; shift;
  elif [[ ( $# -ge 1 ) && ( $1 == "-logScale" ) ]]; then 
    logScale=1; shift;
  elif [[ ( $# -ge 2 ) && ( $1 == "-title" ) ]]; then 
    title="$2"; shift; shift;
  elif [[ ( $# -ge 1 ) && ( $1 == "-verbose" ) ]]; then 
    verbose=1; shift;
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

if [[ ( "/${bins}" == "/0" ) && ( "/${step}" == "/0" ) ]]; then
  bins=200
elif [[ ( "/${bins}" != "/0" ) && ( "/${step}" != "/0" ) ]]; then
  echo "** must specify at most one of \"-step\" or \"-bins\"" 1>&2; exit 1
fi

binRound=$(( ${binRound} + 0 ))

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

col=$(( 10#${channel} + 3 ))

make_histogram.gawk \
    -v col="${col}" \
    -v step="${step}" \
    -v bins="${bins}" \
    -v vmin="${vmin}" \
    -v vmax="${vmax}" \
    -v emin="${emin}" \
    -v emax="${emax}" \
    -v bround="${binRound}" \
    -v verbose="${verbose}" \
    ${fnifile} \
  > ${hisfile}
  
if [[ ! ( -s ${hisfile} ) ]]; then 
  echo "** {make_histogram.gawk} failed, ${hisfile} not generated" 1>&2; exit 1;
fi

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
echo "ymax = ${ymax}" 1>&2

gnuplot <<EOF
set terminal postscript eps color 
set output
set size 1.50,0.75
set nokey 
ymax=${ymax}
if (${logScale} != 0) {
  set yrange [ 0.40:(1.02*ymax)]; set logscale y 
} else { 
  set yrange[(-0.02*ymax):(1.02*ymax)]; unset logscale y
}
set title "${title}"
plot "${hisfile}" using 2:4 with histeps
quit
EOF

/bin/rm -f ${fnifile} ${hisfile} ${maxfile}
