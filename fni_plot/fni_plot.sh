#! /bin/bash 
# Last edited on 2023-01-14 10:47:28 by stolfi
 
PROG_NAME=${0##*/}
PROG_DESC="make a 3D Postscript plot of a selected channel of a FNI file"
PROG_HELP=(
  "${PROG_NAME} \\"
  "\n    [ -channel {CHNUM} ] \\"
  "\n    [ -title {TITLESRING} ] \\"
  "\n    [ -rows ] [ -kodak ] [ -kodakLog ] \\"
  "\n    [ -range {VMIN} {VMAX} ] \\"
  "\n    [ -ztics {VTICS} ] \\"
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
  "\n generates an Encapsulated Postscript 3D plot of a"
  "\n selected channel using gnuplot(1)."
  "\n"
  "\n  The \"-range\" option specifies the nominal range of values"
  "\n for vertical scale setting."
  "\n"
  "\n  The \"-ztics\" option specifies the value interval"
  "\n for vertical axis tics."
  "\n"
  "\n  If \"-rows\" is used, outputs a 2D plot where each"
  "\n row is treated as a curve."
  "\n"
  "\n  If \"-kodak\" is used with a 2D plot, adds a plot of"
  "\n the KODAK Q-13 gray calibration scale."
  "\n"
  "\n  If \"-kodakLog\" is used with a 2D plot, adds a plot of"
  "\n the logarithm of the KODAK Q-13 gray calibration scale."
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
title=""
rows=0
vrange=0
vmin=
vmax=
vtics=0
vstep=
kodak=0
kodakLog=0
while [[ ( $# -ge 1 ) && ( "/$1" =~ ^[/][-]. ) ]]; do
  if [[ ( $# -ge 2 ) && ( "$1" == "-channel" ) ]]; then 
    channel="$2"; shift; shift;
  elif [[ ( $# -ge 1 ) && ( "$1" == "-rows" ) ]]; then 
    rows=1; shift;
  elif [[ ( $# -ge 1 ) && ( "$1" == "-kodakLog" ) ]]; then 
    kodakLog=1; shift;
  elif [[ ( $# -ge 1 ) && ( "$1" == "-kodak" ) ]]; then 
    kodak=1; shift;
  elif [[ ( $# -ge 2 ) && ( "$1" == "-title" ) ]]; then 
    title="$2"; shift; shift;
  elif [[ ( $# -ge 3 ) && ( "$1" == "-range" ) ]]; then 
    vrange=1; vmin="$2"; vmax="$3"; shift; shift; shift;
  elif [[ ( $# -ge 2 ) && ( "$1" == "-ztics" ) ]]; then 
    vtics=1; vstep="$2"; shift; shift;
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

# Save input data in FNI temp file:
tmpfni=${tmp}.fni
cat > ${tmpfni}

# Extract image size:
size=( `fni_size.sh ${tmp}.fni` )
nx="${size[1]}";
ny="${size[2]}";

# Extract plottable data, save in another temp file:
tmpdat=${tmp}.dat
cat ${tmpfni} \
  | gawk \
      ' /^[ ]*[-+]?[0-9]/ { 
          y=$2+0; 
          if(oy!=y) { print ""; } 
          print; 
          oy=y;
        } 
      ' \
  > ${tmpdat}

# Preparations for vertical range:
if [[ ${vrange} -ne 0 ]]; then
  RangeVars="vmin=(${vmin}); vmax=(${vmax}); dv=0.02*(vmax-vmin)"
  if [[ ${vtics} -ne 0 ]]; then
    ZTicsVars="vstep=(${vstep})"
  else
    ZTicsVars="vstep=0.25*(vmax-vmin)"
  fi
else
  RangeVars=""
  if [[ ${vtics} -ne 0 ]]; then
    ZTicsVars="vstep=(${vstep})"
  else
    ZTicsVars=""
  fi
fi

if [[ $rows -ne 0 ]]; then
  
  # Two-dimensional plot
  echo "making two-dimensional plot..." 1>&2
  if [[ ${vrange} -ne 0 ]]; then
    SetRange="set yrange [(vmin-dv):(vmax+dv)]"
  else
    SetRange=""
  fi
  
  SetYTics="set ytics 0.1"
  
  if [[ ${kodak} -ne 0 ]]; then
    KodakPlot=', "'"${tmpdat}"'" using 1:(10.0**-(0.05 + 0.10*(int(20*column(1)/nx)))) with steps lt 3'
  else
    KodakPlot=""
  fi
  
  if [[ ${kodakLog} -ne 0 ]]; then
    KodakLogPlot=', "'"${tmpdat}"'" using 1:(log(10.0**-(0.05 + 0.10*(int(20*column(1)/nx))))) with steps lt 3'
    SetYTics="set ytics 1.0"
  else
    KodakLogPlot=""
  fi
  
gnuplot <<EOF
set terminal postscript eps color 
set output
set size 2,1.5
set xrange [-2:(${nx}+2)]
${RangeVars}
${SetRange}
set nokey
nx=${nx}; ny=${ny}
${SetYTics}
set mytics 5
set grid ytics lt 3, lt 0
#set grid mytics lt 3, lt 0
set title "${title}"
plot "${tmpdat}" using 1:(column(3+${channel})) with lines ${KodakPlot} ${KodakLogPlot}
quit
EOF

else

  # Three-dimensional plot
  echo "making three-dimensional plot..." 1>&2
  if [[ ${vrange} -ne 0 ]]; then
    SetRange="set zrange [(vmin-dv):(vmax+dv)]"
  else
    SetRange=""
  fi
  if [[ ${vtics} -ne 0 ]]; then
    SetZTics="set ztics vstep"
  else
    SetZtics=""
  fi
  
gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 28
set output
set view 30,330,1,1.25
set size 2,2
nx=${nx}; ny=${ny}
sz = (nx>ny?nx:ny)
set xrange [-2:(sz+2)]
set yrange [-2:(sz+2)]
${RangeVars}
${SetRange}
${ZTicsVars}
${SetZTics}
set nokey
set title "${title}"
set hidden3d
splot "${tmpdat}" using 1:2:(column(3+${channel})) with lines linetype 0 linecolor rgb '#0000ff'
quit
EOF

fi

/bin/rm -f ${tmpdat} ${tmpfni}