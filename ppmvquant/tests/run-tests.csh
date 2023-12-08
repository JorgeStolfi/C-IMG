#! /bin/csh -f
# Last edited on 2011-06-05 23:59:05 by stolfi

# Tests of ppmvquant

set lib = "${STOLFIHOME}/PUB/colormaps"

rm -v -f out/*.{ppm,pgm,eps,hist}

foreach dist ( rgb yuv )
  if ( "/$dist" == "/rgb" ) then
    set distop = "-rgb"
  else
    set distop = "-yuv"
  endif

  foreach map ( 002 008 016 064 666 222 by2 )
    if ( "/$map" == "/666" ) then
      set ncop = ( -map ${lib}/standard-6x6x6-colormap-gm10.ppm )
    else if ( "/$map" == "/222" ) then
      set ncop = ( -map ${lib}/standard-2x2x2-colormap-gm10.ppm )
    else if ( "/$map" == "/by2" ) then
      set ncop = ( -map ${lib}/standard-blu-yel-colormap-gm10.ppm )
    else 
      set ncop = "$map"
    endif
    
    foreach mode ( n f )
      if ( "/$mode" == "/f" ) then
        set modeop = "-fs"
      else
        set modeop = "-nofs"
      endif

      foreach prog ( pbm old new )
        if ( "/$prog" == "/pbm" ) then
          set run = ( "ppmquant" )
        else if ( "/$prog" == "/old" ) then
          set run = ( "oldppmvquant" "-${dist}" )
        else
          set run = ( "ppmvquant" "-verbose" "-${dist}" )
        endif
        set tag = "${mode}-${dist}-${map}-${prog}"
        set oname = "out/${tag}"
        echo '===' ${tag} '==='
        set echo
        cat p.ppm | ${run} ${modeop} ${ncop} > "${oname}.ppm"
        unset echo
        pnmarith -difference p.ppm ${oname}.ppm | ppmtopgm > ${oname}-dif.pgm
        pgmhist ${oname}-dif.pgm > ${oname}-dif.hist
        plot-pgm-hist ${tag} 255 ${oname}-dif.hist > ${oname}-dif-hist.eps
        convert ${oname}-dif-hist.eps ${oname}-dif-hist.ppm
      end
      
      set jtag = "${mode}-${dist}-${map}"
      set jout = "out/${jtag}"
      pnmcat -lr ${jout}-{pbm,old,new}.ppm > ${jout}.ppm
      pnmcat -lr ${jout}-{pbm,old,new}-dif.pgm > ${jout}-dif.pgm
      pnmcat -lr ${jout}-{pbm,old,new}-dif-hist.ppm > ${jout}-dif-hist.ppm
    end
  end
end

display `ls out/?-???-???{,-dif,dif-hist}.{ppm,pgm} | sort -t. -k1,1 -k2,2r` &

