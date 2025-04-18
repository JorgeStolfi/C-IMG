#! /bin/bash
# Last edited on 2019-12-07 13:16:51 by jstolfi

cmd=${0##*/}
usage="${cmd} [ -I{DIR} ]... [ {NAME}.c | {NAME}.h ]..."

# Reads a list of C source files (".h" or ".c") and
# outputs a list of dependency lines suitable for "make",
# using the ".ho" mechanism to avoid quadratic cost.

# For each source "{NAME}.{SEXT}", the script outputs
# a list of dependency pairs, one per line:
# 
#    The first line is "{NAME}.{OEXT}: {NAME}.{SEXT}", where
#    {OEXT} is ".o" or ".ho" if {SEXT} is ".c" or ".h", respectively.
# 
#    There follows one line for each header "{INCL}.h"
#    imported by "{NAME}.{SEXT}".  Normally the line 
#    says "{NAME}.{OEXT}: {DIR}/{INC}.ho". The directory
#    name {DIR} is obtained by searching "{INCL}.h" in
#    the list of directories provided with the "-I"
#    directive; files that are not found are flagged.
#    However, if {DIR} starts with a slash "/", the script
#    assumes that it is a system's include, so the dependency
#    line has ".h" instead of ".ho".
# 
# Unlike "makedepend" and "gcc -MM", this script considers
# only the headers directly included by the named source files;
# i.e. the included files themselves are not scanned for 
# indirect includes.  It also considers all "#include" commands,
# even those that are supressed by being inside "#if... #endif"
# or "/* ... */".

ipath=( ) 
while [ $# -gt 0 ]; do
  case "$1" in
    -I* ) 
      ipath=( ${ipath[@]} `echo "$1" | sed -e 's:^-I::' -e 's:[\/]*$::'` )
      shift ;;
    -* )
      error "unrecognized option $1" ;
      echo "usage: ${usage}" 1>&2 ; exit 1 ;;
    * ) break;;
  esac;
done

sources=( $@ ); 

echo "ipath = ( ${ipath[@]} )" 1>&2
# echo "sources = ( ${sources[*]} )"  1>&2

tmp="/tmp/$$"

for src in ${sources[@]} ; do
  name="${src%.*}"
  ext=${src##*.}
  if [ ".${ext}" == ".c" ]; then
    object="${name}.o"
  elif [ ".${ext}" == ".h" ]; then
    object="${name}.ho"
  else 
    echo "invalid extension '${ext}'" 1>&2
    continue
  fi
  
  echo "${object}: ${src}"
  
  cat ${src} \
    | gawk \
        -v src=${src} \
        ' /^ *[#]include *[<"][-a-zA-Z0-9_/+.]+[.][h][">] *$/ { \
            gsub(/^ *[#]include *["<]/, ""); \
            gsub(/[">] *$/, ""); \
            printf "%s:%s\n", $0, FNR; next; \
          } \
          /^ *[#]include/ { \
            printf "%s:%d: bad include \"%s\"\n", src, FNR, $0 > "/dev/stderr" ; \
            next; \
          } \
        ' \
    | egrep -v -e '^stdarg[.]h:' \
    > ${tmp}
    
  # ( echo "includes of ${src}:" ; cat ${tmp} ) 1>&2
    
  for incline in `cat ${tmp}` ; do
    inc="${incline%:*}"; 
    line="${incline##*:}"
    # echo "incline = ${incline}  inc = ${inc}  line = ${line}" 1>&2
    # Ignore include files that are actually internal to GCC or exist in MSDOS only:
    if [[ "/${inc}" == "/float.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/process.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/sys/time.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/sys/times.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/sys/types.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/sys/resource.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/sys/time.h" ]]; then continue; fi;
    if [[ "/${inc}" == "/fpu_control.h" ]]; then continue; fi;
    if [[ "${inc}X" == /* ]]; then
      # The included file's name starts with "/".
      # The dependency is to that file verbatim:
      dep="${inc}"
    else
      # The included file's name does not start with "/".
      # Look in the include path for a directory {dir} that has tha ".ho" or ".h" file:
      dir="/dev/null"
      for adir in ${ipath[@]} /usr/include ; do
        dep="${adir}/${inc}"
        if [[ ( -r "${dep}" ) || ( -r "${dep}o" ) ]]; then
          dir="${adir}"
          break
        fi
      done
      if [[ "${dir}" == "/dev/null" ]]; then
        # Neither ".h"nor ".ho" were found in the include list.
        # Print error and stop:
        printf "%s:%s: ** no %so or %s found\n" \
          "${src}" "${line}" "${inc}" "${inc}"  1>&2
        exit 1
      fi
      dep="${dir}/${inc}"
      if [ -r "${dep}o" ]; then
        # Found the ".ho" file, okay:
        dep="${dir}/${inc}o"
      elif [ -r "${dep}" ]; then
        # The ".h" file (but not the ".ho") was found
        # in directory "${dir}".
        # Check nature of directory: 
        if [[ "${dir}/" == "./" ]]; then
          # The ".h" file was found in the current directory.
          # The dependency is to the "ho" file, even if it doesn't exist:
          dep="${dir}/${inc}o"
        elif { [[ "${dir}/" == /${STOLFIHOME}/* ]] || 
             [[ "${dir}/" == */${STOLFIHOME}/* ]]; }; then
          # The directory seems to be an user's header repository.
          # There must be a ".ho" file next to it:
          printf "%s:%s: ** file %s found but %so file is missing\n" \
           "${src}" "${line}" "${dir}/${inc}" "${dir}/${inc}"  1>&2
          exit 1
        else
          # The directory apparently does not belong to the user.
          # The dependency is to the ".h" file itself:
          dep="${dir}/${inc}"
        fi
      fi
    fi
    # Remove any leading "./":
    dep="${dep#./}"
    printf "%s: %s\n" "${object}" "${dep}"
  done
  printf '\n' 
done

/bin/rm -f ${tmp}
