# Last edited on 2025-04-18 17:38:38 by stolfi

SHELL := /bin/bash

LIBNAME := libgus

IGNORE :=

.PHONY: debug-make uninstall pickle-libs compare-libs build-lib install
 
#---------------------------------------------------------------------
# Default target for "make" with no target: 

all: debug-make uninstall pickle-libs compare-libs build-lib install

# ----------------------------------------------------------------------
# Install folders for derived files and precompiled foreign libs:
GUS_LIB_DIR := ../lib
GUS_INC_DIR := ../include
GUS_BIN_DIR := ../bin

# ----------------------------------------------------------------------
# Package components and directories

# Library file name:
LIBFILE = ${LIBNAME}.a

# Get list of sources:
HFILES := ${filter-out ${IGNORE}, ${wildcard *.h}}
CFILES := ${filter-out ${IGNORE}, ${wildcard *.c}}

# Derived object files:
HOFILES := ${subst .h,.ho,${HFILES}}
OFILES := ${subst .c,.o,${CFILES}}

# X11 headers:
X11_INC_DIR := /usr/include/X11
X11_I_FLAGS := -I${X11_INC_DIR}

# OpenGL headers:
GL_INC_DIR := /usr/include/GL
GL_I_FLAGS := \
  -I${GUS_INC_DIR}/GL \
  -I${GL_INC_DIR}

#---------------------------------------------------------------------
# Directories to search for #include files:

I_FLAGS := \
  -I. \
  -I${GUS_INC_DIR} \
  ${GL_I_FLAGS} \
  ${X11_I_FLAGS} \
  
  
# ----------------------------------------------------------------------
# "make pickle-libs" makes local copies of the JSLIB sources listed in 00-SOURCES.txt

pickle-libs: 00-SOURCES.txt
	cat 00-SOURCES.txt | gawk '/^lib/ { printf "%s/%s\n", $$1, $$2; }' > .libs
	for f in `cat .libs` ; do \
          cp -avu ../../../JSLIBS/$$f ./ ; \
	  chmod a-w $${f/*\/} ; \
        done

# ----------------------------------------------------------------------
# "make compare-libs" compares the pickled sources in this folder 
# with the official sources in JSLIBS:

compare-libs: 00-SOURCES.txt compare_libs.sh
	./compare_libs.sh 

#---------------------------------------------------------------------
# "make depend" recreates the source dependency file ${DEPFILE}.

DEPFILE := Deps.make

depend: ${DEPFILE}

${DEPFILE}: ${HFILES} ${CFILES}
	@echo "creating dependency file ${DEPFILE} ..." 1>&2
	@/bin/rm -f ${DEPFILE}
	@${GUS_BIN_DIR}/extract_ho_deps.sh ${I_FLAGS} -I/usr/include \
            ${HFILES} ${CFILES} \
          | egrep -v ': /usr/include' \
          > ${DEPFILE}

# ----------------------------------------------------------------------
# Common C compiler flags for library files:

CC := /usr/bin/gcc

C_FLAGS := \
  -Wall \
  -Wcast-align \
  -Wchar-subscripts \
  -Wclobbered \
  -Wconversion \
  -Wmissing-field-initializers \
  -Wmissing-parameter-type \
  -Wmissing-prototypes \
  -Wpointer-arith \
  -Wstrict-prototypes \
  -Wuninitialized \
  -Wundef \
  -Werror \
  -ggdb \
  -fpcc-struct-return \
  -ffloat-store \
  -frounding-math \
  -std=gnu99 \
  -D_USE_MATH_DEFINES

# -Wno-sign-conversion \
# -Wno-type-limits \

#---------------------------------------------------------------------
# "make build" and "make build-lib" make sure that the ${DEPFILE}
# is up to date, and then recreate the object files from the source files
# as necessary.

build: build-lib

build-lib: ${DEPFILE} ${LIBFILE}
  
${LIBFILE}: ${HOFILES} ${OFILES}
	@echo "creating ${LIBFILE} ..." 1>&2
	-rm -f $*.a
	ar crv $*.a ${OFILES}
	ranlib $*.a

%.o: %.c
	${CC} -c ${C_FLAGS} ${I_FLAGS} $*.c
  
%.ho: %.h
	${CC} -o $*.ho -c ${C_FLAGS} ${I_FLAGS} -x c $*.h \
          || /bin/rm -f $*.ho

ifneq "/${wildcard ${DEPFILE}}" "/"
  # Include specific dependencies extracted by "make depend"
  include ${DEPFILE}
endif

#---------------------------------------------------------------------
# "make install" copies program and makefile to the public dir.

install: install-lib

install-lib: ${HFILES} ${HOFILES} ${LIBFILE}
	@echo "installing headers and library ..." 1>&2
	@for hho in ${HFILES} ${HOFILES} ; do \
          compdir="$${PWD}" ; \
          ( cd ${GUS_INC_DIR}/ && rm -f $${hho} ; ln -s $${compdir}/$${hho} ; ls -ld $${hho} ) ; \
        done
	@compdir="$${PWD}" ; ( cd ${GUS_LIB_DIR}/ && rm -f ${LIBFILE} ; ln -s $${compdir}/${LIBFILE} ; ls -ld ${LIBFILE} )

#---------------------------------------------------------------------
# "make uninstall" deletes the exported program and makefile:

uninstall: uninstall-lib

uninstall-lib:
	@echo "uninstalling headers and library ..." 1>&2
	( cd ${GUS_INC_DIR}/. && rm -f ${HFILES} ${HOFILES} )
	( cd ${GUS_LIB_DIR}/. && rm -f ${LIBFILE} )

#---------------------------------------------------------------------
# "make debug"-make prints some debugging info:
        
debug-make:
	@echo "--- make variables ---" 1>&2
	@echo "LIBFILE =       ${LIBFILE}" 1>&2
	@echo "LIBNAME =       ${LIBNAME}" 1>&2
	@echo "DEPFILE =       ${DEPFILE}" 1>&2
	@echo "GUS_BIN_DIR =   ${GUS_BIN_DIR}" 1>&2
	@echo "GUS_INC_DIR =   ${GUS_INC_DIR}" 1>&2
	@echo "GUS_LIB_DIR =   ${GUS_LIB_DIR}" 1>&2
	@echo "IGNORE =        ${IGNORE}" 1>&2
	@echo "CFILES =        ${CFILES}" 1>&2
	@echo "HFILES =        ${HFILES}" 1>&2
	@echo "HOFILES =       ${HOFILES}" 1>&2
	@echo "OFILES =        ${OFILES}" 1>&2
	@echo "GL_INC_DIR =    ${GL_INC_DIR}" 1>&2
	@echo "GL_I_FLAGS =    ${GL_I_FLAGS}" 1>&2
	@echo "X11_INC_DIR =   ${X11_INC_DIR}" 1>&2
	@echo "X11_I_FLAGS =   ${X11_I_FLAGS}" 1>&2
	@echo "I_FLAGS =       ${I_FLAGS}" 1>&2
	@echo "C_FLAGS =       ${C_FLAGS}" 1>&2
	@echo "------" 1>&2

#---------------------------------------------------------------------
# "make clean" deletes all derived files.
# It is a double-colon rule so that clients may add more actions.
#---------------------------------------------------------------------
# "make clean" deletes all derived files.

clean: clean-lib

clean-lib:
	@echo "removing derived files ..." 1>&2
	-/bin/rm -f *.o *.ho *.a core
