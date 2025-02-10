# Last edited on 2025-02-04 19:38:12 by stolfi

PROG = make_sphere_textures

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	${PROG}
	display -title '%f' `ls out/test-*.png | sort`

clean::
	/bin/rm -fv out/test-*.ppm
