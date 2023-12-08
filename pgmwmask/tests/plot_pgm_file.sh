#! /bin/bash
# Last edited on 2012-12-15 10:11:35 by stolfilocal

inp_pgm="$1"; shift
tmp_prefix="/tmp/$$"
tmp_fni="${tmp_prefix}.fni"
tmp_eps="${tmp_prefix}.eps"

pnm_to_fni -min 0 -max 1 < ${inp_pgm} > ${tmp_fni}

fni_plot.sh -range 0.0 1.0 -ztics 0.2 < ${tmp_fni} > ${tmp_eps}

gv ${tmp_eps}

rm -f ${tmp_fni} ${tmp_eps} 
