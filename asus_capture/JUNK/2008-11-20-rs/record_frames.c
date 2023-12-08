
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
#include <stdlib.h>


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

void InitFFMPEG(void){
	// Register all formats and codecs
	avcodec_register_all();
	avdevice_register_all();
	av_register_all();
}

int InitCaptureDevice(char* deviceName,AVFormatContext **pFormatCtx,AVCodecContext  **pCodecCtx,AVCodec **pCodec,AVInputFormat **iformat,AVFormatParameters *formatParams,int* streamId){
	

	int  videoStream,i;
	//Init params with default options
	(*formatParams).standard = NULL; 
	(*formatParams).width = 640;
	(*formatParams).height = 480;
	(*formatParams).time_base = (AVRational){1, 30};
	(*iformat) = av_find_input_format("video4linux2");

	//Detects if Video4Linux2 format is reconized
	fprintf(stderr,"Detecting Device Format...");
	if (*iformat == NULL) {
        	printf("No format found\n");
        	return -2;
	}
	fprintf(stderr,"OK\n");

	// Open video input device
	fprintf(stderr,"Opening Device...");
	if(av_open_input_file(pFormatCtx,deviceName, *iformat, 0, formatParams)!=0){
		fprintf(stderr,"couldn't open device !\n");
		return -1; // Couldn't open file
	}
	fprintf(stderr,"OK\n");

	// Retrieve stream information
	fprintf(stderr,"Retrieving stream information...");
	if(av_find_stream_info(*pFormatCtx)<0){
		fprintf(stderr,"not found !\n");
		return -1; // Couldn't find stream information
	}
	fprintf(stderr,"OK\n");

	// Dump information about file onto standard error
	fprintf(stderr,"*********************Device Information**********************\n");
	dump_format(*pFormatCtx, 0, deviceName, 0);
	fprintf(stderr,"*************************************************************\n");

	// Find the first video stream avaiable
	fprintf(stderr,"Detecting video stream...");
	videoStream=-1;
	for(i=0; i<(*pFormatCtx)->nb_streams; i++){
		if((*pFormatCtx)->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) {
			videoStream=i;
			break;
		}
	}
	if(videoStream==-1){
		fprintf(stderr,"not found !");
		return -1; // Didn't find a video stream
	}
	*streamId = videoStream;
	fprintf(stderr,"OK\n");

	// Get a pointer to the codec context for the video stream
	(*pCodecCtx) = (*pFormatCtx)->streams[videoStream]->codec;
	// Find the decoder for the video stream
	fprintf(stderr,"Locating device decoder...");
	*pCodec=avcodec_find_decoder((*pCodecCtx)->codec_id);
	if(*pCodec==NULL) {
		fprintf(stderr, "unsupported codec !\n");
		return -1; // Codec not found
	}

	// Open codec
	if(avcodec_open(*pCodecCtx, *pCodec)<0){
		fprintf(stderr, "codec could not be opened !\n");
		return -1; // Could not open codec
	}
	fprintf(stderr,"OK\n");

	//OK, all damn things finished
	return 1;
}

AVFrame*  AllocateFrame(AVCodecContext  *pCodecCtx,int pixel_fmt,uint8_t **buffer){
	AVFrame* pFrameRGB;
	int numBytes;
	//Alocate basic data structures, buffer is allocated in a different section (why ?)
	pFrameRGB=avcodec_alloc_frame();
	if(pFrameRGB==NULL){
		return NULL;
	}
	
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(pixel_fmt, pCodecCtx->width,pCodecCtx->height);
	*buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));
	if(*buffer == NULL){
		return NULL;
	}
	// Assign appropriate parts of buffer to image planes in pFrameRGB
	// Note that pFrameRGB is an AVFrame, but AVFrame is a superset of AVPicture
	avpicture_fill((AVPicture *)pFrameRGB, *buffer, pixel_fmt,pCodecCtx->width, pCodecCtx->height);
	
	
	return pFrameRGB;
}

void FreeFrame(AVFrame* frame,uint8_t* buffer){
	//Just free all needed structures
	if(buffer != NULL){
		av_free(buffer);
	}
	av_free(frame);
}

void ConvertFrameFormat(AVFrame* sourceFrame,int sourceWidth, int sourceHeight,int sourceFmt, AVFrame* destFrame,int destWidth, int destHeight,int destFmt){
	//this function converts to a frame format to another, reescaling it if needed
	struct SwsContext* encoderSwsContext;
	encoderSwsContext = sws_getContext(sourceWidth, sourceHeight,sourceFmt , destWidth, destHeight,destFmt , SWS_BICUBIC, NULL, NULL, NULL); 
	sws_scale(encoderSwsContext, sourceFrame->data, sourceFrame->linesize, 0, sourceHeight, destFrame->data, destFrame->linesize);
}


int main(int argc, char *argv[]) { 
	

	AVFormatContext *pFormatCtx;
	AVInputFormat *iformat = NULL;
	AVFormatParameters formatParams;
	int             i;
	AVCodecContext  *pCodecCtx;
	AVCodec         *pCodec;
	AVFrame         *pFrame; 
	//AVFrame         *pFrameRGB;
	AVPacket        packet;
	int             frameFinished;
	int             numBytes;
	uint8_t         *buffer;
	int numFrames,videoStream;
	
	if(argc < 3) {
		printf("Program Usage: record_frames <number of frames> <output dir>\n");
		return -1;
	}
	if(sscanf(argv[1],"%d",&numFrames) != 1){
		printf("Please provide a VALID number of  frames!\n");
		return -1;
	}
	char* dest_dir = argv[2];

	//Init FFMPEG codecs, compatible devices and so on	
	InitFFMPEG();
	//Do all those damn things needed to start a video capture
	fprintf(stderr,"Initializing Capture Device\n");
	if(InitCaptureDevice("/dev/video0",&pFormatCtx,&pCodecCtx,&pCodec,&iformat,&formatParams,&videoStream) < 0){
		fprintf(stderr,"Capture Device Setup failure !");
		return -1;	
	}
	fprintf(stderr,"Capture Device Setup sucessfull !\n");
	
	// Allocate video frame buffer packets got from video stream
	pFrame=avcodec_alloc_frame();
	//Alocate a vector that will contain frames to be grabbed
	fprintf(stderr,"Alocating frame structures...");
	AVFrame** vecFrame = (AVFrame**)malloc(sizeof(AVFrame*)*numFrames);
	//Data buffer to each frame must be allocated separated
	uint8_t** vecBuffer = (uint8_t**)malloc(sizeof(uint8_t*)*numFrames);
	for(i = 0; i < numFrames; i++){
		vecFrame[i] = AllocateFrame(pCodecCtx,pCodecCtx->pix_fmt,&(vecBuffer[i]));
		if(vecFrame[i] == NULL){
			fprintf(stderr,"failed !");
			return -1;
		}
	}
	fprintf(stderr,"OK\n");


	i=0;
	int count = 0;
	int failed = 0;
	fprintf(stderr,"Starting Capture\n");
	double startCaptureTime = get_current_time();
	
	while(1){
		//Read frame... if it fails it means a corrupted frame, so we just discard it
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
				//Copy frame from buffer to our vector
				av_picture_copy((AVPicture *)vecFrame[count],(AVPicture *)pFrame,pCodecCtx->pix_fmt,pCodecCtx->width, pCodecCtx->height);
				count++;
			}
		}
		
		// Free the packet that was allocated by av_read_frame
		av_free_packet(&packet);
		if(count >= numFrames) break;
	}
	//Show usefull(?) data
	double stopCaptureTime = get_current_time();	
	double fps = count/(stopCaptureTime - startCaptureTime);
	fprintf(stderr,"\nFinished Capture - %lf FPS\n",fps);
	fprintf(stderr,"%d corrupted packets discarded\n",failed);
	fprintf(stderr,"Saving frames (please wait)... ");
	
	//Save all captured frames
	AVFrame* pFrameRGB = AllocateFrame(pCodecCtx,PIX_FMT_RGB24,&buffer);
	for(i = 0; i < numFrames; i++){
		ConvertFrameFormat(vecFrame[i],pCodecCtx->width, pCodecCtx->height,pCodecCtx->pix_fmt, pFrameRGB,pCodecCtx->width, pCodecCtx->height,PIX_FMT_RGB24);
		SaveFrame(pFrameRGB,pCodecCtx->width, pCodecCtx->height,i,dest_dir);
		FreeFrame(vecFrame[i],vecBuffer[i]);
		
	}
	//free the remaining structures
	free(vecFrame);
	free(vecBuffer);
	fprintf(stderr,"OK");
	// Free the YUV frame
	av_free(pFrame);
	// Close the codec
	avcodec_close(pCodecCtx);
	//Close the video file
	av_close_input_file(pFormatCtx);
	//Byeeee
	return 0;
}
