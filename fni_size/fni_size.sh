#! /bin/bash 
# Last edited on 2023-01-14 10:46:55 by stolfi
 
PROG_NAME=${0##*/}
PROG_DESC="print the number of channels, columns, and rows of a FNI file"
PROG_HELP=(
  "${PROG_NAME} \\"
  "\n    [ -format {FMTSTRING} ] \\"
  "\n    [<] {INFILE}.fni"
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
  "\n  prints out its dimensions as three integers \"{NC} {NX} {NY}\""
  "\n  (number of channels, number of columns, number of rows)."
  "\n"
  "\n  If the \"-format\" option is given, the integers are"
  "\n  printed with the given {FMTSTRING}."
  "\n"
  "\n  Prints \"0 0 0\" if it fails to parse the input file header."
  "\n"
  "\nAUTHOR"
  "\n  Created 2006-04-02 by Jorge Stolfi, Unicamp"
)

# Parse command line options: 

format="%d %d %d"
filename=()
while [[ ( $# -ge 1 ) && ( "/$1" =~ /-.* ) ]]; do
  if [[ ( $# -ge 2 ) && ( "/$1" == "/-format" ) ]]; then 
    format=$2; shift; shift;
  elif [[ ( $# -ge 1 ) && ( ( "/$1" == "/-help" ) || ( "/$1" == "/--help" ) ) ]]; then 
    echo -e "usage:\n  ${PROG_HELP[@]}"; exit 0;
  elif [[ ( $# -ge 1 ) && ( ( "/$1" == "/-info" ) || ( "/$1" == "/--info" ) ) ]]; then 
    echo -e "${PROG_INFO[@]}"; exit 0;
  else
    echo "unknown option $1" 1>&2 ;
    echo -e "usage:\n  ${PROG_HELP[@]}" 1>&2 ; exit 1 
  fi
done 

if [[ $# -ge 1 ]]; then
  filename=( $1 ); shift;
fi

if [[ $# -gt 0 ]]; then
  echo 'excess arguments "'"$1"'" ...' 1>&2 ;
  echo -e "usage:\n  ${PROG_HELP[@]}" 1>&2 ; exit 1 
fi

gawk \
  -v fmt="${format}" \
  ' \
    BEGIN { NC = -1; NX = -1; NY = -1; nlin = 0 } \
    /^[>|#]/ { next; } \
    // { nlin++; } \
    /^begin/{ next; } \
    /[=]/ { gsub(/[=]/, " = ", $0); } \
    /^NC *[=]/{ NC = $3; next; } \
    /^NX *[=]/{ NX = $3; next; } \
    /^NY *[=]/{ NY = $3; next; } \
    /^([ ]*[-+]?[0-9]|end)/ { exit 0; } \
    (nlin > 10) { exit 0; } \
    END { \
      if ((NC < 0) || (NX < 0) || (NY < 0)) { NC = 0; NX = 0; NY = 0; } \
      printf fmt, NC, NX, NY; exit 0; \
    } \
  ' \
  "${filename[@]}"

exit 0
