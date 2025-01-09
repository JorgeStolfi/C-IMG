#! /bin/bash
# Last edited on 2024-12-31 07:36:54 by stolfi

# Tests of ppmvquant

(cd out && rm -f *.{ppm,pgm,eps,hist} )

inppm="in/laclogop.ppm"

for space in rgb yuv; do
  if [[ "/$space" == "/rgb" ]]; then
    space_op="-rgb"
  else
    space_op="-yuv"
  fi

  for palette in 002 008 016 064 666 222 by2 ; do
    if [[  "/$palette" == "/666"  ]]; then
      mapFile="in/standard-6x6x6-colormap-gm10.ppm"; maxColors=""
    elif [[  "/$palette" == "/222"  ]]; then
      mapFile="in/standard-2x2x2-colormap-gm10.ppm"; maxColors=""
    elif [[  "/$palette" == "/by2"  ]]; then
      mapFile="in/standard-blu-yel-colormap-gm10.ppm"; maxColors=""
    else 
      mapFile=""; maxColors="${palette}"
    fi
    
    for floyd in F T; do 
      for prog in pbm old new ; do
        if [[  "/$prog" == "/pbm"  ]]; then
          run=( "ppmquant" )
          if [[ "/${maxColors}" == "/" ]]; then
            palette_op=( "-map" ${mapFile} )
          else
            palette_op=( ${maxColors} )
          fi
          if [[ "/${floyd}" == "/T" ]]; then
            floyd_op=( "-floyd" )
          else
            floyd_op=(  )
          fi
        elif [[  "/$prog" == "/old"  ]]; then
          run=( "oldppmvquant" "${space_op}" )
          if [[ "/${maxColors}" == "/" ]]; then
            palette_op=( "-map" ${mapFile} )
          else
            palette_op=( ${maxColors} )
          fi
          if [[ "/${floyd}" == "/T" ]]; then
            floyd_op=( "-floyd" )
          else
            floyd_op=(  )
          fi
        else
          run=( "ppmvquant" "-verbose" "${space_op}" )
          if [[ "/${maxColors}" == "/" ]]; then
            palette_op=( "-map" ${mapFile} )
          else
            palette_op=( "-maxColors" ${maxColors} )
          fi
          floyd_op=( "-floyd" ${floyd} )
        fi
        tag="${space}-fs${floyd}-${palette}-${prog}"
        oname="out/${tag}"
        
        echo '===' ${tag} '===' 1>&2
        echo "cat ${inppm} | ${run} ${floyd_op[@]} ${palette_op[@]}" 1>&2
        cat ${inppm} | ${run[@]} ${floyd_op[@]} ${palette_op[@]} > "${oname}.ppm"
        if [[ -s ${oname}.ppm ]]; then
          pnmarith -difference ${inppm} ${oname}.ppm | ppmtopgm > ${oname}-dif.pgm
          pnmxhist ${oname}-dif.pgm > ${oname}-dif.hist
          ./plot_pnm_hist.sh ${oname}-dif.hist ${oname}-dif-hist.png 255
        fi
      done
      
      jtag="${space}-fs${floyd}-${palette}"
      jout="out/${jtag}"
      pnmcat -td ${jout}-{pbm,old,new}.ppm > ${jout}.ppm
      pnmcat -td ${jout}-{pbm,old,new}-dif.pgm > ${jout}-dif.pgm
      pnmcat -td ${jout}-{pbm,old,new}-dif-hist.png > ${jout}-dif-hist.png
    done
  done
done

display `ls out/?-???-???{,-dif,dif-hist}.{ppm,pgm} | sort -t. -k1,1 -k2,2r` &

