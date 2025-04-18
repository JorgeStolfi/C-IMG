#! /bin/bash
# Last edited on 2025-04-18 16:54:41 by stolfi

echo "comparing pickled library files with originals in JSLIBS ..." 1>&2
ndiff=0
for libfil in `cat 00-SOURCES.txt | sed -e 's:[ ]*[#].*$::g' -e 's:[ ][ ]*$::g' -e 's:[ ][ ]*:@:g'` ; do
  lib=${libfil/@*}
  file=${libfil/*@}
  ofile="../../../JSLIBS/${lib}/${file}"
  
  if ! cmp -s ${ofile} ${file} ; then
    printf "%s comparing ${ofile} ${file} %s\n" "---" "---"
    prdiff -Bb ${ofile} ${file}
    ndiff=$(( ${ndiff} + 1 ));
  fi
done

if [[ ${ndiff} -gt 0 ]]; then
  echo "** differences detected - fix and recompile" 1>&2
  exit 1
fi
