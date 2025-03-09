#! /bin/bash
# Last edited on 2025-02-17 02:20:20 by stolfi

inp_pgm="$1"; shift
tmp_prefix="/tmp/$$"
tmp_fni="${tmp_prefix}.fni"
tmp_png="${tmp_prefix}.png"

pnm_to_fni -min 0 -max 1 < ${inp_pgm} > ${tmp_fni}

fni_plot.sh -range 0.0 1.0 -ztics 0.2 < ${tmp_fni} > ${tmp_png}

display ${tmp_png}

rm -f ${tmp_fni} ${tmp_png} 
