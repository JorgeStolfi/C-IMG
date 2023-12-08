#! /bin/bash
# Last edited on 2023-10-15 12:22:54 by stolfi

OT_FILE="$1"; shift
EX_FILE="$1"; shift
OT_EX_DIF_FILE="$1"; shift
SHOW="$1"; shift

# Compares a text file {OT_FILE} containing key feature points or a map matrix file with a reference file {EX_FILE},
# writes the differences to {OT_EX_DIF_FILE}, prints them to {stderr} if {SHOW} is "YES".

tmp="/tmp/$$"

function prt(){
  title="$1"; shift
  file="$1"; shift
  if [[ "/${SHOW}" == "/YES" ]]; then
    echo "${title}p:" 1>&2
    echo ".................................................." 1>&2
    cat "${file}" 1>&2
    echo ".................................................." 1>&2
  fi
}

rm -f -v "${OT_EX_DIF_FILE}"
if [[ ! ( -s "${OT_FILE}" ) ]]; then \
  echo "** ${OT_FILE} not generated" 1>&2
else
  prt "Output projective map" "${OT_FILE}"
fi
if [[ ! ( -s "${EX_FILE}" ) ]]; then \
  echo "** ${EX_FILE} not found" 1>&2
else
  prt "Expected output projective map" "${EX_FILE}"
fi
if [[ ( -s "${OT_FILE}" ) && ( -s "${EX_FILE}" ) ]]; then
  prdiff -Bb "${OT_FILE}" "${EX_FILE}" > "${OT_EX_DIF_FILE}"
  if [[ ! ( -s "${OT_EX_DIF_FILE}" ) ]]; then
    prt "Differences" "${OT_EX_DIF_FILE}"
  else
    echo "No differences." 1>&2
  fi
fi
