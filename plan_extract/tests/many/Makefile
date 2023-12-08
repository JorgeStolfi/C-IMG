# Last edited on 2017-07-16 13:41:31 by stolfilocal

PROG := plan_extract
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TEST := \
  2007-03-01-google

TESTS_NOT := \
  2007-04-07-stolfi

SAMPLE := smp

DURFILE := ${TEST}/in/${SAMPLE}.dur
ARGFILE := ${TEST}/out/${SAMPLE}.arg
OUTFILE := ${TEST}/out/${SAMPLE}.out

all: ${OUTFILE}
# all: ${ARGFILE}

clean::
	rm -f ${ARGFILE} ${OUTFILE}

SHOW := display
HRULE := echo ' --------------------------------------------- '

${ARGFILE}: ${DURFILE} Makefile
	cat ${DURFILE} \
	  | gawk -v tt="${TEST}" \
             ' /./{ \
                 fn=$$2; nx=$$6; ny = $$7; \
                 pf=fn; gsub(/[\/]/, "", pf); gsub(/[.].*$$/, "", pf); \
                 printf "  -input out/%s/%s %5s %5s", tt, pf, nx, ny; \
                 printf "  -take %5s %5s 0 0", nx, ny; \
                 printf "  -pad tlbr"; \
                 printf "\n"; \
               } \
             ' \
          > ${ARGFILE}
	@${HRULE}

${OUTFILE}: ${ARGFILE} ${PROGDIR}/${PROG} Makefile
	head -10 ${ARGFILE} \
          | xargs \
            --max-lines=1 \
            --no-run-if-empty \
            ${PROGDIR}/${PROG} \
              -minSize 320 240 \
              -maxSize 320 240 \
              -maxScale 0.75 \
              -maxPad 0.1 \
              -debug 1
	touch ${OUTFILE}
	cat $@
	@${HRULE}
