#! /bin/bash
# Last edited on 2024-10-31 22:34:32 by stolfi

OT_FILE="$1"; shift
EX_FILE="$1"; shift
OT_EX_DIF_FILE="$1"; shift
SHOW="$1"; shift

# Compares a mapped image {OT_FILE} with a reference image {EX_FILE},
# computes the difference image {OT_EX_DIF_FILE},
# displays them if {SHOW} is "YES".

tmp="/tmp/$$"

rm -f -v "${OT_EX_DIF_FILE}"

imgs=()
if [[ ! ( -s "${OT_FILE}" ) ]]; then \
  echo "** ${OT_FILE} not generated" 1>&2
else
  imgs+=( "${OT_FILE}" )
fi
if [[ ! ( -s "${EX_FILE}" ) ]]; then \
  echo "** ${EX_FILE} not found" 1>&2
else
  imgs+=( "${EX_FILE}" )
fi
if [[ ( -s "${OT_FILE}" ) && ( -s "${EX_FILE}" ) ]]; then
  pnmxarith -subtract -offset 0.5 "${OT_FILE}" "${EX_FILE}" > "${OT_EX_DIF_FILE}"
  imgs+=( "${OT_EX_DIF_FILE}" )
fi
if [[ "/${SHOW}" == "/YES" ]]; then
  tfile="${tmp}.png"
  convert "${imgs[@]}" -background black +append ${tfile}
  display -title "${OT_FILE}" -resize '>x800' ${tfile}
  rm ${tfile}
fi
