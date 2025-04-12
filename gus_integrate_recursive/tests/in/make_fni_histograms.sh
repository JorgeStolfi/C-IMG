#! /bin/bash
# Last edited on 2025-04-09 15:07:11 by stolfi

show=1

for dir in */[0-9][0-9][0-9][0-9]x[0-9][0-9][0-9][0-9] ; do
  pushd ${dir}
  for ffile in *.fni ; do 
    if [[ -L ${ffile} ]]; then
      echo "!! ${dir}/${ffile} is a symbolic link, skipped" 1>&2
    elif [[ -f ${ffile} ]]; then
      nc="`cat ${ffile} | grep -e "NC" | sed -e 's:^.*= *::g' -e 's: *$::g'`"
      nc=$(( 10#${nc} + 0 ))
      for ch in 0 1 2 3 4 5 6 7 9 ; do
        if [[ $ch -lt $nc ]]; then
          hfile="${ffile/.fni}-ch${ch}-hist.eps"
          echo "plotting histogram of ${ffile} channel ${ch}" 1>&2
          title="${dir}/${ffile} channel ${ch}"
          fni_hist.sh -channel ${ch} -logScale -title "${title}" -bins 200 -verbose < ${ffile} > ${hfile}
          if [[ ${show} -ne 0 ]]; then evince ${hfile} ; fi
        fi
      done
    else
      echo "** invalid file type ${dir}/${ffile}" 1>&2 ; exit 1
    fi
  done
  
  popd
done
