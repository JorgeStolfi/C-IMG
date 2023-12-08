#define PROG_NAME "asus-capture"
#define PROG_DESC "captures a sequence of frames fom the ASUS eeePC camera"
#define PROG_VERS "1.0"

#define asus_capture_C_COPYRIGHT "Copyright © 2008 by the State University of Campinas (UNICAMP)"

/* Last edited on 2017-06-24 01:46:26 by stolfilocal */

#define PROG_HELP \
  " " PROG_NAME " \\\n" \
  "    [ -format [ rgb | yuv ] ] \\\n" \
  "    [ -prefix {PREFIX} ] \\\n" \
  "    {NUM_FRAMES} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program captures a specified number of frames from the" \
  " {ASUS eeePC} built-in camera, as fast as it can, and writes" \
  " them to disk in the PPM or PGM format --- hopefully avoiding" \
  " the artifacts of JPEG compression.\n" \
  "\n" \
  "  In order to maximize the frame rate, the program saves all" \
  " frames in memory during the capture, in the camera's native" \
  " format, and writes them out only at the end.  Therefore, the" \
  " number of frames is limited by the amount of real memory" \
  " available.  The use of virtual memory may have an impact" \
  " on the frame rate.\n" \
  "\n" \
  "  Each frame can be saved to disk as a single PPM file in" \
  " the RGB color system, called \"{PREFIX}{NNNNNN}-rgb.ppm\"; or" \
  " as three separate PGM files in the YUV color system," \
  " called \"{PREFIX}{NNNNNNN}-{CHAN}.pgm\"; where {CHAN} is" \
  " either \"Y\", \"U\", or \"V\", {NNNNNN} is a six-digit frame" \
  " index (starting from \"000000\"), and {PREFIX} is a pathname" \
  " prefix specified by the \"-prefix\" option." \
  "\n" \
  "  The RGB and Y images have the nominal SVGA resolution," \
  " 640 by 480 pixels. The U and V images have only half as" \
  " many columns, namely 320 by 480, and need to be stretched" \
  " horizontally by a factor of 2.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -format [ RGB | YUV ]\n" \
  "    This optional argument specifies the output format.  The" \
  " possible values are \"RGB\" (one PPM file in RGB encoding) or" \
  " \"YUV\" (three separate PGM files with the Y, U, and V" \
  " channels).  The default is \"-format RGB\".\n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    This optional argument specifies a common prefix for all" \
  " output file names.  If the {PREFIX} contains slashes, all" \
  " relevant directories must exist.  The" \
  " default is \"-prefix \'\'\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmtopng(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2008-11 by Rafael Saracchini and Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " asus_capture_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>
/* #include <unistd.h> */
/* #include <sys/times.h> */

typedef time_t time_foo_t;

#include <avformat.h>

typedef AVFormatParameters foo_t;

#include <avdevice.h>
#include <rtsp.h>
#include <avcodec.h>
/* #include <opt.h> */
#include <swscale.h>

#undef FALSE
#undef TRUE
#include <bool.h>
#include <jstime.h>
#include <argparser.h>

typedef struct options_t 
  { int numFrames;    /* Number of frames to capture. */
    bool_t rgbFormat; /* TRUE write RGB files, FALSE write Y/U/V files. */
    char *prefix;     /* Prefix for output file names. */
  } options_t;
  
options_t *acp_parse_options(int argc, char *argv[]);
  /* Parses the command-line arguments and returns them as an {options_t} record. */

void acp_save_frame(AVFrame *pFrame, int width, int height, int iFrame, char *dst_dir);

void acp_init_ffmpeg(void);
   /* Register all {ffmpeg} formats and codecs. */
    
void acp_init_device(char *deviceName, AVFormatContext **pFormatCtx, AVCodecContext  **pCodecCtx, AVCodec **pCodec, AVInputFormat **iFormat, AVFormatParameters *formatParams, int *streamId);
  /* Fails if cannoy open the camera and decoder. */

void acp_alloc_frame_and_buffer(AVCodecContext  *pCodecCtx, int pixel_fmt, AVFrame **frame, uint8_t **buffer);

void acp_free_frame(AVFrame *frame, uint8_t *buffer);

void acp_convert_frame(AVFrame* srcFrame, int srcWidth, int srcHeight, int srcFmt, AVFrame* dstFrame, int dstWidth, int dstHeight, int dstFmt);
  /* Converts a frame format to another, rescaling it if needed. */
    
int main(int argc, char *argv[]);

void acp_init_ffmpeg(void)
  { avcodec_register_all();
    avdevice_register_all();
    av_register_all();
  }

void acp_save_frame(AVFrame *pFrame, int width, int height, int iFrame, char *dst_dir)
  { char *fname = NULL;
    asprintf(&fname, "%s/frame%03d.ppm", dst_dir, iFrame);
    FILE *pFile = (FILE*)notnull(fopen(fname, "wb"), "could not open the output file");
    // Write header
    fprintf(pFile, "P6\n%d %d\n255\n", width, height);
    // Write pixel data
    int  y;
    for (y = 0; y < height; y++)
      { fwrite(pFrame->data[0]+y*pFrame->linesize[0], 1, width*3, pFile); }
    // Close file
    fclose(pFile);
  }

void acp_init_device(char *devName, AVFormatContext **pFormatCtx, AVCodecContext  **pCodecCtx, AVCodec **pCodec, AVInputFormat **iFormat, AVFormatParameters *formatParams, int *streamId)
  {
    //Init params with default options
    (*formatParams).standard = NULL;
    (*formatParams).width = 640;
    (*formatParams).height = 480;
    (*formatParams).time_base = (AVRational){1, 30};
    fprintf(stderr, "checking whether the input format \"video4linux22\" is recognized... ");
    (*iFormat) = (AVInputFormat *)notnull(av_find_input_format("video4linux2"), "not recognized");
    fprintf(stderr, "OK\n");

    fprintf(stderr,"opening camera device... ");
    int device_err = av_open_input_file(pFormatCtx, devName, *iFormat, 0, formatParams);
    demand(device_err == 0, "cannot open device!\n");
    fprintf(stderr,"OK\n");

    fprintf(stderr,"retrieving the camera stream information... ");
    int stream_info_err = av_find_stream_info(*pFormatCtx);
    demand(stream_info_err >= 0, "not found!\n");
    fprintf(stderr,"OK\n");

    // Dump information about file onto standard error
    fprintf(stderr,"*********************Device Information**********************\n");
    dump_format(*pFormatCtx, 0, devName, 0);
    fprintf(stderr,"*************************************************************\n");

    fprintf(stderr, "finding the first available video stream... ");
    *streamId = -1;
    int i;
    for (i = 0; (i < (*pFormatCtx)->nb_streams) && ((*streamId) < 0); i++)
      { if ((*pFormatCtx)->streams[i]->codec->codec_type == CODEC_TYPE_VIDEO)
          { (*streamId) = i; }
      }
    demand((*streamId) >= 0, "not found !");
    fprintf(stderr,"OK\n");

    // Get a pointer to the codec context for the video stream
    (*pCodecCtx) = (*pFormatCtx)->streams[(*streamId)]->codec;
    
    fprintf(stderr,"locating the device decoder... ");
    *pCodec = avcodec_find_decoder((*pCodecCtx)->codec_id);
    demand((*pCodec) != NULL,  "unsupported codec!\n");
    fprintf(stderr,"OK\n");

    fprintf(stderr,"opening the device decoder stream... ");
    int codec_open_err = avcodec_open(*pCodecCtx, *pCodec);
    demand(codec_open_err >= 0,  "codec could not be opened!\n");
    fprintf(stderr,"OK\n");
  }

void acp_alloc_frame_and_buffer(AVCodecContext  *pCodecCtx, int pixel_fmt, AVFrame **frame, uint8_t **buffer)
  { //Alocate basic data structures, buffer is allocated in a different section (why ?)
    (*frame) = (AVFrame *)notnull(avcodec_alloc_frame(), "no space for frame");

    // Determine required buffer size and allocate buffer
    int numBytes = avpicture_get_size(pixel_fmt, pCodecCtx->width, pCodecCtx->height);
    (*buffer) = (uint8_t *)notnull(av_malloc(numBytes*sizeof(uint8_t)), "no space for buffer");
    
    // Assign appropriate parts of buffer to image planes in pFrameRGB
    // Note that pFrameRGB is an AVFrame, but AVFrame is a superset of AVPicture
    avpicture_fill((AVPicture *)(*frame), (*buffer), pixel_fmt, pCodecCtx->width, pCodecCtx->height);
  }

void acp_free_frame(AVFrame *frame, uint8_t *buffer)
  { //Just free all needed structures
    if (buffer != NULL) { av_free(buffer); }
    av_free(frame);
  }

void acp_convert_frame(AVFrame *srcFrame, int srcWidth, int srcHeight, int srcFmt, AVFrame *dstFrame, int dstWidth, int dstHeight, int dstFmt)
  { struct SwsContext* encoderSwsContext;
    encoderSwsContext = sws_getContext(srcWidth, srcHeight, srcFmt, dstWidth, dstHeight, dstFmt, SWS_BICUBIC, NULL, NULL, NULL);
    sws_scale(encoderSwsContext, srcFrame->data, srcFrame->linesize, 0, srcHeight, dstFrame->data, dstFrame->linesize);
  }


int main(int argc, char *argv[])
  { 
    // Parse command line arguments:
    options_t *o = acp_parse_options(argc, argv);

    // Initialize the ffmpeg codecs, compatible devices and so on: 	
    acp_init_ffmpeg();

    // Start the video capture: 
    char *devName = "/dev/video0"; // The ASUS eeePC camera device.
    AVFormatContext *pFormatCtx;
    AVCodecContext *pCodecCtx;
    AVCodec *pCodec;
    AVInputFormat *iFormat;
    AVFormatParameters formatParams;
    int videoStreamId;
    acp_init_device(devName, &pFormatCtx, &pCodecCtx, &pCodec, &iFormat, &formatParams, &videoStreamId);

    // Allocate frame header and data for packets got from video stream:
    AVFrame *pFrameYUV = avcodec_alloc_frame();
    
    fprintf(stderr,"alocating storage for grabbed frames... ");
    AVFrame **vecFrameYUV = (AVFrame **)notnull(malloc(sizeof(AVFrame *)*o->numFrames), "no mem for frames");
    uint8_t **vecBufferYUV = (uint8_t **)notnull(malloc(sizeof(uint8_t *)*o->numFrames), "no mem for buffers");
    int i;
    for (i = 0; i < o->numFrames; i++)
      { acp_alloc_frame_and_buffer(pCodecCtx, pCodecCtx->pix_fmt, &(vecFrameYUV[i]), &(vecBufferYUV[i])); }
    fprintf(stderr,"OK\n");

    fprintf(stderr,"capturing %d frames... \n", o->numFrames);
    double startCaptureTime = real_time_usec();
    int count = 0;
    int npk_good = 0;
    int npk_bad = 0;
    int npk_other = 0;
    while (count < o->numFrames)
      { // Read the next packet from the AV stream:
	AVPacket packet;
	int frame_read_err = av_read_frame(pFormatCtx, &packet);
        if (frame_read_err < 0)
          { if (verbose) { fprintf(stderr, "(D)\n"); }
            npk_bad++;
          }
        else if (packet.stream_index != videoStreamId)
          { if (verbose) { fprintf(stderr, "(O)\n"); }
            npk_other++;
          }
        else 
          { if (verbose) { fprintf(stderr, "."); }
            // Decode the video frame:
            int frameFinished;
            avcodec_decode_video(pCodecCtx, pFrameYUV, &frameFinished, packet.data, packet.size);
            // Did we get a video frame?
            if (frameFinished)
              { if (verbose) { fprintf(stderr, "[%d]", count); }
                //Copy frame from buffer to our vector
                av_picture_copy((AVPicture *)vecFrameYUV[count],(AVPicture *)pFrameYUV, pCodecCtx->pix_fmt, pCodecCtx->width, pCodecCtx->height);
                count++;
              }
            npk_good++;
          }
        // Free the packet that was allocated by av_read_frame
        av_free_packet(&packet);
      }
    double stopCaptureTime = real_time_usec();	
    fprintf(stderr,"captured %d frames\n", count);
    double fps = 1.0e6 * count/(stopCaptureTime - startCaptureTime);
    double mspf = (stopCaptureTime - startCaptureTime)/count / 1.0e3;
    fprintf(stderr,"mean rate = %.6lf frames/sec (%.2lf ms/frame)\n", fps, mspf);
    fprintf(stderr,"%d good packets\n", npk_good);
    fprintf(stderr,"%d corrupted packets (discarded)\n", npk_bad);
    fprintf(stderr,"%d packets with wrong stream id (ignored)\n", npk_other);

    // Save frames to disk:
    AVFrame *pFrameRGB;
    uint8_t *bufferRGB;
    acp_alloc_frame_and_buffer(pCodecCtx, PIX_FMT_RGB24, &pFrameRGB, &bufferRGB);
    for (i = 0; i < numFrames; i++)
      { acp_convert_frame(vecFrameYUV[i], pCodecCtx->width, pCodecCtx->height, pCodecCtx->pix_fmt, pFrameRGB, pCodecCtx->width, pCodecCtx->height, PIX_FMT_RGB24);
        acp_save_frame(pFrameRGB, pCodecCtx->width, pCodecCtx->height, i, dst_dir);
        acp_free_frame(vecFrameYUV[i], vecBufferYUV[i]);
      }
      
    //free the remaining structures
    free(vecFrameYUV);
    free(vecBufferYUV);
    // Free the YUV frame
    av_free(pFrameYUV);
    // Close the codec
    avcodec_close(pCodecCtx);
    //Close the video file
    av_close_input_file(pFormatCtx);
    //Byeeee
    return 0;
  }
options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    if (argparser_keyword_present(pp, "-format"))
      { char *fmt = argparser_get_next_non_keyword(pp); 
        if (strcmp(fmt, "RGB") == 0)
          { o->rgbFormat = TRUE; }
        else if (strcmp(fmt, "YUV") == 0)
          { o->rgbFormat = FALSE; }
        else 
          { argparser_error(pp, "unrecognized \"-format\""); }
      }
    else
      { o->rgbFormat = TRUE; }

    if (argparser_keyword_present(pp, "-prefix"))
      { o->prefix = argparser_get_next(pp); }
    else
      { o->prefix = ""; }

    argparser_skip_parsed(pp);

    o->numFrames = argparser_get_next_int(pp, 0, MAX_NUM_FRAMES);
    
    argparser_finish(pp);
    return o;
  }
