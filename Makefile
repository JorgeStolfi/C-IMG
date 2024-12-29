# Last edited on 2024-08-30 23:40:53 by stolfi

PROGS_BUG := \
  pnmpairtopng \

PROGS_TO_FINISH := \
  asus_capture \
  fni_std_extract \
  geostereo \
  image_match_features \
  image_stitch \
  image_stitch_n \
  langev \
  ms_img_match \
  multifok \
  multifok_analyze \
  multifok_make_stack \
  optdither \
  pnm_asus_eee_pc_701 \
  pnmalign \
  pnmbend \
  pnmclassif \
  pnmhog \
  pnmhough \
  pnmmatch \
  pnmpralign \
  pnmxhist \
  ppvip \
  stereo_adjust \
  volstereo \
  
IGNOREDIRS := \
  ${PROGS_BUG} \
  ${PROGS_TO_FINISH}
   
all: build install

include ${STOLFIHOME}/programs/c/GENERIC-ROOT-DIR.make

