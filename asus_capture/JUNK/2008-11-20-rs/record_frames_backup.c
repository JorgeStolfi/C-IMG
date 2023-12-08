
#include <math.h>
#include <limits.h>
#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libavformat/rtsp.h>
#include <libavdevice/avdevice.h>
#include <libswscale/swscale.h>
#include <libavcodec/opt.h>
#include <libswscale/swscale.h>
#include <stdio.h>


#include <sys/times.h>
#include <unistd.h>

double get_current_time(){
	 struct timespec buf;
    	 clock_gettime(CLOCK_REALTIME, &buf);
    	 return (((double)buf.tv_sec)*1.0) + (((double)buf.tv_nsec)/1.0e9);
	
}

void SaveFrame(AVFrame *pFrame, int width, int height, int iFrame,char* dest_dir) {
	FILE *pFile;
	char szFilename[32];
	int  y;
	// Open file
	sprintf(szFilename, "%s/frame%03d.ppm",dest_dir, iFrame);
	pFile=fopen(szFilename, "wb");
	if(pFile==NULL)
	return;

	// Write header
	fprintf(pFile, "P6\n%d %d\n255\n", width, height);
	// Write pixel data
	for(y=0; y<height; y++)
		fwrite(pFrame->data[0]+y*pFrame->linesize[0], 1, width*3, pFile);
	// Close file
	fclose(pFile);
}

int main(int argc, char *argv[]) {
	AVFormatContext *pFormatCtx;
	int             i, videoStream;
	AVCodecContext  *pCodecCtx;
	AVCodec         *pCodec;
	AVFrame         *pFrame; 
	AVFrame         *pFrameRGB;
	AVPacket        packet;
	int             frameFinished;
	int             numBytes;
	uint8_t         *buffer;
	int numFrames;
	
	if(argc < 3) {
		printf("Program Usage: record_frames <number of frames> <output dir>\n");
		return -1;
	}
	if(sscanf(argv[1],"%d",&numFrames) != 1){
		printf("Please provide a VALID number of  frames!\n");
		return -1;
	}
	char* dest_dir = argv[2];
	// Register all formats and codecs
	avcodec_register_all();
	avdevice_register_all();
	av_register_all();
	AVFormatParameters formatParams ;
	AVInputFormat *iformat = NULL;
	formatParams.standard = NULL; 
	formatParams.width = 640;
	formatParams.height = 480;
	formatParams.time_base = (AVRational){1, 30};
	iformat = av_find_input_format("video4linux2");
	fprintf(stderr,"Deteting Device Format...");
	if (iformat == NULL) {
        	printf("No format found\n");
        	return -2;
	}
	fprintf(stderr,"OK\n");

	fprintf(stderr,"Opening Device...");
	// Open video file
	if(av_open_input_file(&pFormatCtx,"/dev/video0", iformat, 0, &formatParams)!=0){
		fprintf(stderr,"couldn't open device !\n");
		return -1; // Couldn't open file
	}
	fprintf(stderr,"OK\n");
	// Retrieve stream information
	fprintf(stderr,"Retrieving stream information...");
	if(av_find_stream_info(pFormatCtx)<0){
		fprintf(stderr,"not found !\n");
		return -1; // Couldn't find stream information
	}
	fprintf(stderr,"OK\n");
	fprintf(stderr,"*********************Device Information**********************\n");
	// Dump information about file onto standard error
	dump_format(pFormatCtx, 0, "/dev/video0", 0);
	fprintf(stderr,"*************************************************************\n");
	// Find the first video stream
	fprintf(stderr,"Detecting video stream...");
	videoStream=-1;
	for(i=0; i<pFormatCtx->nb_streams; i++){
		if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) {
			videoStream=i;
			break;
		}
	}
	if(videoStream==-1){
		fprintf(stderr,"not found !");
		return -1; // Didn't find a video stream
	}
	fprintf(stderr,"OK\n");
	// Get a pointer to the codec context for the video stream
	pCodecCtx=pFormatCtx->streams[videoStream]->codec;
	// Find the decoder for the video stream
	fprintf(stderr,"Locating device decoder...");
	pCodec=avcodec_find_decoder(pCodecCtx->codec_id);
	if(pCodec==NULL) {
		fprintf(stderr, "unsupported codec !\n");
		return -1; // Codec not found
	}
	// Open codec
	if(avcodec_open(pCodecCtx, pCodec)<0){
		fprintf(stderr, "codec could not be opened !\n");
		return -1; // Could not open codec
	}
	fprintf(stderr,"OK\n");
	// Allocate video frame
	pFrame=avcodec_alloc_frame();
	// Allocate an AVFrame structure
	fprintf(stderr,"Alocating frame structures...");
	pFrameRGB=avcodec_alloc_frame();
	if(pFrameRGB==NULL){
		fprintf(stderr,"failed !");
		return -1;
	}
	fprintf(stderr,"OK\n");
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(PIX_FMT_RGB24, pCodecCtx->width,pCodecCtx->height);
	buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));
	// Assign appropriate parts of buffer to image planes in pFrameRGB
	// Note that pFrameRGB is an AVFrame, but AVFrame is a superset
	// of AVPicture
	avpicture_fill((AVPicture *)pFrameRGB, buffer, PIX_FMT_RGB24,pCodecCtx->width, pCodecCtx->height);
	
	// Read frames and save first five frames to disk
	i=0;
	int count = 0;
	int failed = 0;
	fprintf(stderr,"Starting Capture\n");
	double startCaptureTime = get_current_time();
	//while(av_read_frame(pFormatCtx, &packet)>=0) {
	while(1){
		if(av_read_frame(pFormatCtx, &packet) < 0){
			fprintf(stderr, "WARNING... corrupted packet discarded\n");
			failed++;
			continue;
		}
		// Is this a packet from the video stream?
		if(packet.stream_index==videoStream) {
			// Decode video frame
			avcodec_decode_video(pCodecCtx, pFrame, &frameFinished,packet.data, packet.size);
			// Did we get a video frame?
			if(frameFinished) {
				struct SwsContext* encoderSwsContext;
				encoderSwsContext = sws_getContext(640, 480,pCodecCtx->pix_fmt , 640, 480,PIX_FMT_RGB24 , SWS_BICUBIC, NULL, NULL, NULL); 
			        sws_scale(encoderSwsContext, pFrame->data, pFrame->linesize, 0, 640, pFrameRGB->data, pFrameRGB->linesize); 
				// Convert the image from its native format to RGB
				//img_convert((AVPicture *)pFrameRGB, PIX_FMT_RGB24, (AVPicture*)pFrame, pCodecCtx->pix_fmt, pCodecCtx->width, pCodecCtx->height);
				// Save the frame to disk
				//if(++i<=5){
				SaveFrame(pFrameRGB, pCodecCtx->width, pCodecCtx->height,count,dest_dir);
				//count++;
				//}*/
				count++;
				//fprintf(stderr,"%d ",count);
			}
		}
		//fprintf(stderr,"%d frames produced against %d\n",count,numFrames);
		// Free the packet that was allocated by av_read_frame
		av_free_packet(&packet);
		if(count > numFrames) break;
	}
	double stopCaptureTime = get_current_time();	
	double fps = count/(stopCaptureTime - startCaptureTime);
	fprintf(stderr,"\nFinished Capture - %lf FPS\n",fps);
	fprintf(stderr,"%d corrupted packets discarded\n",failed);
	// Free the RGB image
	av_free(buffer);
	av_free(pFrameRGB);
	// Free the YUV frame
	av_free(pFrame);
	// Close the codec
	avcodec_close(pCodecCtx);
	//Close the video file
	av_close_input_file(pFormatCtx);
	return 0;
}
