# Last edited on 2024-08-30 23:24:54 by stolfi

PROG := ppvip

TEST_LIB := libppv.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  do-test

do-test: ${PROG}
	${PROG}

