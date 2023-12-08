#! /bin/csh -f
# Last edited on 2002-10-04 14:10:21 by stolfi

set cmd = "$0"; set cmd = "${cmd:t}"
set usage = "${cmd} DSPFILE IMGFILE..."

if ( $#argv < 2 ) then
  echo "usage: ${usage}"; exit 1
endif

set dspfile = "$1"; shift
set imgfiles = ( $* )
set shifts = ( `cat ${dspfile}` )

@ i = 2
@ j = 3
foreach f ( ${imgfiles} )
  echo " $f ${shifts[$i]} ${shifts[$j]}"
  @ i = $i + 3
  @ j = $j + 3
end


