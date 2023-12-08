# Last edited on 2017-07-28 02:39:24 by jstolfi
 
PROG := mavica-calib
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

PLOT := do-plot

IMAGES := \
  dsc00001 \
  dsc00002 \
  dsc00003 \
  dsc00004 \
  dsc00005 \
  dsc00006 \
  dsc00007 \
  dsc00008
  
# Input files:
MED_FILES := ${addprefix{in/, ${addsuffix -med.txt,${IMAGES}}}

# Output files:
RES_FILES := ${addprefix{out/, ${addsuffix -res.txt,${IMAGES}}}
PLO_FILES := ${addprefix{out/, ${addsuffix -plo.txt,${IMAGES}}}
PLO_PS_FILES := ${addprefix{out/, ${addsuffix -plo.ps,${IMAGES}}}
  
all: ${PROGDIR}/${PROG} ${PLO_PS_FILES} show-blacks show-gammas
	@echo ${PLO_PS_FILES}

clean::

out/%-plo.ps: in/%-med.txt ${PROGDIR}/${PROG} ${PROGDIR}/${PLOT}
	${PROGDIR}/${PROG} -r $* >& $*.log
	${PROGDIR}/${PLOT} $*

show-blacks: ${PLO_FILES}
	egrep '[#] *Preto' ${RES_FILES}

show-gammas: ${PLO_FILES}
	egrep '[#] *Gama' ${RES_FILES}
