# Last edited on 2024-08-30 23:27:25 by stolfi

PROG := asus_capture

OTHER_I_FLAGS := \
  -I/usr/include/ffmpeg-compat \
  -I/usr/include/ffmpeg-compat/libavformat \
  -I/usr/include/ffmpeg-compat/libavdevice \
  -I/usr/include/ffmpeg-compat/libavcodec \
  -I/usr/include/ffmpeg-compat/libswscale \
  -I/usr/include/gstreamer-1.0/gst/rtsp \
  -I/usr/include/gstreamer-1.0 \
  -I/usr/include/glib-2.0 \
  -I/usr/lib64/glib-2.0/include \

BOTHER_I_FLAGS := \
  -I/usr/include \
  -I/usr/include/ffmpeg-compat/libavutil \

BOTHER_I_FLAGS :=

JS_LIBS := \
  libimg.a  \
  libgeo.a\
  libjs.a
  
OTHER_LIBS := \
  /usr/lib64/libavutil.so \
  

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

