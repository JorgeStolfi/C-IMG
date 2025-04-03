# Last edited on 2025-04-03 16:59:34 by stolfi

SHELL := /bin/bash

PROG = gus_integrate_recursive

TESTS_DIR := tests

GUS_LIB_SRC_DIR := lib-src
GUS_LIB := libgus.a

# ----------------------------------------------------------------------
# Install folders for derived files and precompiled foreign libs:
GUS_LIB_DIR := lib
GUS_INC_DIR := include
GUS_BIN_DIR := bin

# Where to install scripts and programs:

# ----------------------------------------------------------------------
# Package components and folders

# Get list of all local ".h" files, and their derived objects:
HFILES := 
HOFILES := ${subst .h,.ho,${HFILES}}

# Get list of all local ".c" files (including the main sources of ${PROG}):
CFILES := 
OFILES := ${subst .c,.o,${CFILES}}

# Get list of ".c" files of good main programs and their objects:
PROG_CFILES := ${addsuffix .c,${PROG}}
PROG_OFILES := ${subst .c,.o,${PROG_CFILES}}

# ----------------------------------------------------------------------
# Locations of imported packages

# X11 libraries and headers:
X11_INC_DIR := /usr/include/X11 
X11_LIB_DIR := /usr/lib/x86_64-linux-gnu
X11_LIBRARIES := \
  ${X11_LIB_DIR}/libX11.so \
  ${X11_LIB_DIR}/libXext.so \
  ${X11_LIB_DIR}/libXmu.so.6 \
  ${X11_LIB_DIR}/libXt.so.6 \
  ${X11_LIB_DIR}/libXi.so.6 \
  ${X11_LIB_DIR}/libSM.so.6 \
  ${X11_LIB_DIR}/libICE.so.6
X11_I_FLAGS :=  -I${X11_INC_DIR}

# OpenGL libraries and headers:
GL_INC_DIR := /usr/include/GL
GL_LIB_DIR := /usr/lib/x86_64-linux-gnu
GL_LIBRARIES := \
  ${GUS_LIB_DIR}/libglut.so \
  ${GL_LIB_DIR}/libGLU.so \
  ${GL_LIB_DIR}/libGL.so
GL_I_FLAGS := \
  -I ${GUS_INC_DIR}/GL \
  -I ${GL_INC_DIR}

#---------------------------------------------------------------------
# Phony "make" targets:

.PHONY:: \
  all debug-make \
  uninstall-prog \
  build build-prog \
  run-tests \
  install-prog \
  clean

# ----------------------------------------------------------------------
# Common C compiler flags for all libraries and programs:

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

LD_FLAGS := \
  -ggdb

#---------------------------------------------------------------------
# Default target for "make" with no target: 
all: build-prog run-tests

#---------------------------------------------------------------------
# Folders to search for #include files:

I_FLAGS := \
  -I${GUS_INC_DIR} \
  ${GL_I_FLAGS} \
  ${X11_I_FLAGS}

build-prog: ${PROG}

LIBRARIES := \
  ${GUS_LIB_DIR}/${GUS_LIB} \
  ${addprefix ${GUS_LIB_DIR}/, ${GUS_LIB}} \
  ${GL_LIBRARIES} \
  ${X11_LIBRARIES}

%.o: %.c
	${CC} -c ${C_FLAGS} ${I_FLAGS} $*.c
  
%.ho: %.h
	${CC} -o $*.ho -c ${C_FLAGS} ${I_FLAGS} -x c $*.h \
          || /bin/rm -f $*.ho
  
${PROG}: ${PROG_OFILES} ${LIB_OFILES} ${LIBRARIES}
	@echo '# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 1>&2
	@echo "PROG_OFILES = ${PROG_OFILES}" 1>&2
	@echo "LIB_OFILES = ${LIB_OFILES}" 1>&2
	@echo "LIBRARIES = ${LIBRARIES}" 1>&2
	@echo building $@ ... 1>&2
	-rm -f $@
	${CC} -o $@ ${LD_FLAGS} $@.o \
          ${LIB_OFILES} \
          ${LIBRARIES} -lm -lrt
	@echo '# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        
${GUS_LIB_DIR}/${GUS_LIB}: \
	  ${wildcard ${GUS_LIB_SRC_DIR}/*.h ${GUS_LIB_SRC_DIR}/*.c} \
	  ${wildcard ${GUS_INC_DIR}/*.ho}
	cd ${GUS_LIB_SRC_DIR} && ${MAKE} ${GUS_LIB}

#---------------------------------------------------------------------
# "make install" is a no-op:

install:
	( cd ${GUS_BIN_DIR} && ln -s ../${PROG} ) 

#---------------------------------------------------------------------
# "make uninstall" is a no-op:

uninstall: 
	( cd ${GUS_BIN_DIR} && rm ${PROG} ) 

#---------------------------------------------------------------------
# "make clean" deletes all derived files.
# It is a double-colon rule so that clients may add more actions.

clean::
	-/bin/rm -f *.o core ${PROG}
	( cd ${GUS_LIB_SRC_DIR} && make clean )
	( cd ${TESTS_DIR}  & make clean )

#---------------------------------------------------------------------
# "make debug-make" prints some debugging info:
debug-make:
	@echo "PROG = ${PROG}" 1>&2
	@echo "CC = ${CC}" 1>&2
	@echo "GUS_BIN_DIR = ${GUS_BIN_DIR}" 1>&2
	@echo "GUS_INC_DIR = ${GUS_INC_DIR}" 1>&2
	@echo "GUS_LIB_DIR = ${GUS_LIB_DIR}" 1>&2
	@echo "GUS_LIB_SRC_DIR = ${GUS_LIB_SRC_DIR}" 1>&2
	@echo "GUS_LIB = ${GUS_LIB}" 1>&2
	@echo "HFILES = ${HFILES}" 1>&2
	@echo "HOFILES = ${HFILES}" 1>&2
	@echo "CFILES = ${CFILES}" 1>&2
	@echo "OFILES = ${CFILES}" 1>&2
	@echo "PROG_CFILES = ${PROG_CFILES}" 1>&2
	@echo "PROG_OFILES = ${PROG_OFILES}" 1>&2
	@echo "LIBRARIES = ${LIBRARIES}" 1>&2
	@echo "LIB_OFILES = ${LIB_OFILES}" 1>&2
	@echo "GL_LIB_DIR = ${GL_LIB_DIR}" 1>&2
	@echo "GL_LIBRARIES = ${GL_LIBRARIES}" 1>&2
	@echo "X11_LIB_DIR = ${X11_LIB_DIR}" 1>&2
	@echo "X11_LIBRARIES = ${X11_LIBRARIES}" 1>&2
	@echo "GL_I_FLAGS = ${GL_I_FLAGS}" 1>&2
	@echo "X11_INC_DIR = ${X11_INC_DIR}" 1>&2
	@echo "X11_I_FLAGS = ${X11_I_FLAGS}" 1>&2
	@echo "I_FLAGS = ${I_FLAGS}" 1>&2
	@echo "LD_FLAGS = ${LD_FLAGS}" 1>&2
	@echo "C_FLAGS = ${C_FLAGS}" 1>&2
	@echo "TESTS_DIR = ${TESTS_DIR}" 1>&2
