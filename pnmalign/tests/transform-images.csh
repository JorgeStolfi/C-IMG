#! /bin/csh -f
# Last edited on 2002-10-04 17:30:58 by stolfi

set cmd = "$0"; set cmd = "${cmd:t}"
set usage = "${cmd} DSPFILE IMGFILE..."

if ( $#argv < 2 ) then
  echo "usage: ${usage}"; exit 1
endif

set dspfile = "$1"; shift
set infiles = ( $* )
set shifts = ( `cat ${dspfile}` )

set path = ( ${STOLFIHOME}/bin/${PLATFORM} ${path} )

set otfiles = ( )
@ i = 2
@ j = 3
foreach inf ( ${infiles} )
  set tx = "${shifts[$i]}"
  set ty = "${shifts[$j]}"
  set otf = "${inf:r}-shf.${inf:e}"
  echo "translating ${inf} by ${tx} ${ty} --> ${otf}"
  pnmgtran \
      -matrix 1 ${tx} ${ty} 0 1 0  0 0 1 ${inf} \
    > ${otf}
  set otfiles = ( ${otfiles} ${otf} )
  @ i = $i + 3
  @ j = $j + 3
end

set nimg = $#otfiles
if ( ${nimg} > 1 ) then
  set wt = `gawk -v n=${nimg} 'BEGIN{printf "%7.5f", 1.0/n;}'`
  set previmg = "${otfiles[1]}"
  set prevwt = ${wt}
  set ext = ${previmg:e}
  set tmpa = "/tmp/$$-a.${ext}"
  set tmpb = "/tmp/$$-b.${ext}"
  foreach otf ( ${otfiles[2-]} )
    if ( "/${otf:e}" != "/${ext} ) then
      echo "warning - inhomogeneous images ${previmg} ${otf}"
    endif
    pnmxarith -mix ${prevwt} ${wt} ${previmg} ${otf} > ${tmpb}
    mv ${tmpb} ${tmpa}
    set previmg = "${tmpa}"
    set prevwt = 1.0
  end
  mv ${previmg} ${avgfile}
  /bin/rm -f ${tmpa} ${tmpb}
endif

