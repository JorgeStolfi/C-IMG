# Last edited on 2003-09-25 16:19:02 by stolfi


INC = /home/staff/stolfi/PUB/include
LIB = /home/staff/stolfi/PUB/${PLATFORM}/lib

HFILES =

HOFILES =

OFILES =
  
PROG = \
  make-test-image
  
MAINOBJ = \
  ${PROG}.o
  
IMAGE = \
  test-image
  
LIBS = \
  ${LIB}/libjs.a
  
GCCFLAGS = \
  -I${INC} \
  -g \
  -ansi \
  -Wall -Wtraditional -Wpointer-arith -Wmissing-prototypes

all: cleanup ${PROG}

cleanup: ;\
  rm -f ${PROG}
  
%.o: %.c ;\
  gcc -c ${GCCFLAGS} $*.c
  
%.ho: %.h ;\
  gcc -o $*.ho -c ${GCCFLAGS} -x c $*.h \
  || /bin/rm -f $*.ho
  
${PROG}: ${OFILES} ${MAINOBJ} ${LIBS} ;\
  gcc -o ${PROG} ${MAINOBJ} ${OFILES} ${LIBS} -lm &&\
  ${PROG} &&\
  display ${IMAGE}.pgm

install:

# Dependencies of .h files: 
  
# Dependencies for .c files:

make-test-image.o::  make-test-image.c ${INC}/js.ho ${INC}/ps.ho

