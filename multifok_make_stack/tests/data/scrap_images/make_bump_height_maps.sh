#! /bin/bash 
# Last edited on 2018-07-02 09:23:42 by stolfilocal

htdir="programs/c/PST/make_test_slope_maps/tests/out/0512x0512"

# Not copied:
#   11 low freq wave sin
#   12 low freq wave cos
#   13 mid freq wave sin
#   14 mid freq wave cos
#   16 hi freq wave cos
#   19 sphere again

nnhts=( \
  00:flat05 \
  01:gradix \
  02:gradiy \
  03:gradxy \
  04:parabo \
  05:sphere \
  06:pyram5 \
  07:conemt \
  08:ripple \
  09:babels \
  10:scisso \
  15:hiwave \
  17:babelv \
  18:trampv \
  20:arenas \
  21:circ3s \
  22:trians \
  23:pentas \
  24:circ1v \
  25:circ1s \
  26:qtbump \
  27:frcone \
  )
  
for nnht in ${nnhts[@]}; do 
    nn="${nnht%:*}"
    ht="${nnht#*:}"
    ifile="`ls ~/${htdir}/test-${nn}-1-05-2-11-A-*-0000-0000-Z.pgm | tail`"
    ofile="ht_${ht}_01.pgm"
    cat ${ifile} | pgmnorm -bsingle -wsingle >${ofile}
  done
  
# Create some special files:

pnmxarith -multiply ht_gradix_01.pgm ht_circ1s_01.pgm > ht_gradix_02.pgm
pnmxarith -multiply ht_gradiy_01.pgm ht_circ1s_01.pgm > ht_gradiy_02.pgm
pnmxarith -multiply ht_gradxy_01.pgm ht_circ1s_01.pgm > ht_gradxy_02.pgm
pnmxarith -multiply ht_scisso_01.pgm ht_circ1s_01.pgm > ht_scisso_02.pgm
pnmxarith -multiply ht_sphere_01.pgm ht_hiwave_01.pgm > ht_hiwave_02.pgm

convert ht_arenas_01.pgm -negate ht_arenas_02.pgm
display ht_arenas_02.pgm

# Create a composite of seven bumps:

convert +append ht_babelv_01.pgm ht_pyram5_01.pgm PGM:- \
  | pnmpad -black -left=256 -right=257 -top=69 -bottom=957 \
  > .seven_top.pgm
  
convert +append ht_arenas_02.pgm ht_qtbump_01.pgm ht_hiwave_02.pgm PGM:- \
  | pnmpad -black -left=0 -right=0 -top=513 -bottom=513 \
  > .seven_mid.pgm
  
convert +append ht_sphere_01.pgm ht_gradxy_02.pgm PGM:- \
  | pnmpad -black -left=256 -right=257 -top=957 -bottom=69 \
  > .seven_bot.pgm
  
pnmxarith -maximum .seven_top.pgm .seven_bot.pgm > .seven_tb.pgm
pnmxarith -maximum .seven_tb.pgm .seven_mid.pgm > .seven_all.pgm

convert .seven_all.pgm -resize '513x' ht_sevens_01.pgm
identify ht_sevens_01.pgm
display ht_sevens_01.pgm 

