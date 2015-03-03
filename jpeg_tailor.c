#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <syslog.h>
#include <sys/stat.h>
#include <jpeglib.h>
#include <string.h>
#include <time.h>
#include <math.h>
/* Jpeg-Tailor, by Oren Watson */
#define max(a,b) (a>b?a:b)


/*global variables to allow gradient function*/
unsigned char *bmp_buffer;
int row_stride, width, height, pixsiz, seams;

/* function computes gradient at i,j. Gradient is
magnitude of X derivative + Y derivative*/
float gradient(int i,int j){
	int k,offset=width*j+i;
	float xdiff=0,ydiff=0;
	for(k=0;k<pixsiz;k++)xdiff += fabs(bmp_buffer[offset*pixsiz+k]
		- bmp_buffer[(offset+1)*pixsiz+k]);
	if(j>=height-1)return xdiff;
	for(k=0;k<pixsiz;k++)ydiff = fabs(bmp_buffer[offset*pixsiz+k]
		- bmp_buffer[(offset+width)*pixsiz+k]);
	return xdiff+ydiff;
}

/*function to output gradient as grayscale image.*/
void output_energy(float * energy, int curwid){
		int i,j;
		/*first convert to bitmap in memory*/
		unsigned char *nrg_bmp = malloc(height*(curwid-1));
		for(j=0;j<height;j++)
		for(i=0;i<curwid-1;i++){
			nrg_bmp[j*(curwid-1)+i]=(energy[j*width+i]/(6*255)*255);
		}
		FILE *nrgf = fopen("nrg.pgm","w");
		/* NetPGM format. See
		http://netpbm.sourceforge.net/doc/pgm.html */
		fprintf(nrgf,"P5 %d %d 255\n",curwid-1,height);
		fwrite(nrg_bmp,1,height*(curwid-1),nrgf);
		fclose(nrgf);
		free(nrg_bmp);
}


int main (int argc, char *argv[]) {
	int rc, i, j, k;
	
	if (argc < 3 || argc > 4) {
		fprintf(stderr, "USAGE: %s filename.jpg N\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	
	srand(time(0));
	
	/* Some code copied/adapted from libjpeg documentation, see
	http://libjpeg.cvs.sourceforge.net/viewvc/libjpeg/libjpeg/example.c?view=markup */
	struct jpeg_decompress_struct dinfo;
	struct jpeg_error_mgr jerr;
	unsigned long bmp_size;
	seams = strtol(argv[2],0,0);
	FILE* file = fopen(argv[1], "r");
	dinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&dinfo);
	jpeg_stdio_src(&dinfo, file);
	(void) jpeg_read_header(&dinfo, TRUE);
	printf("Reading in jpeg...\n");
	jpeg_start_decompress(&dinfo);
	width = dinfo.output_width;
	height = dinfo.output_height;
	pixsiz = dinfo.output_components;
	bmp_size = width * height * pixsiz;
	bmp_buffer = (unsigned char*) malloc(bmp_size);
	row_stride = width * pixsiz;
	while (dinfo.output_scanline < dinfo.output_height) {
		unsigned char *buffer_array[1];
		buffer_array[0] = bmp_buffer + \
			(dinfo.output_scanline) * row_stride;
		jpeg_read_scanlines(&dinfo, buffer_array, 1);
	}
	jpeg_finish_decompress(&dinfo);
	jpeg_destroy_decompress(&dinfo);
	
	float *energy = malloc(width*height*sizeof(float));
	/*** Calculate Energy function (gradient)*/
	int q=0;
	for(j=0;j<height;j++)
	for(i=0;i<width-1;i++){
		energy[j*width+i] = gradient(i,j);
	}
	
	/*output gradient as grayscale, uncomment if needed*/
	/*output_energy(energy,width);*/
	
	/*allocate some buffer for summing*/
	float *minsum = malloc(width*height*sizeof(float));
	void *scrap = malloc(width*max(sizeof(float),pixsiz));
	for(k=0;k<seams;k++){
		int curwid = width-k-1;
		/*we are going to sum down the minimum path to each point*/
		for(i=0;i<curwid;i++){
			minsum[i]=energy[i];
		}
		for(j=0;j<height;j++){
			minsum[j*width+curwid]=INFINITY;
			/*this is necessary to prevent problems as the
			image becomes less wide.*/
		}
		
		/** This loop is highly optimized (as it should be, because it 
		is a part of the program that is O(N^3).)
		I optimized by trying various ways, as it turns out this
		way is the best because of SSE instructions on modern processors.*/
		for(j=1;j<height;j++){
			float a,b,c;
			a = INFINITY;/*we don't want the first pixel to go left*/
			b = minsum[width*j-width];
			int oft = width*j;/*offset begins at the start of the scanline*/
			int end = oft + curwid;
			for(;oft<end;oft++){
				c = minsum[oft+1-width];/*get the above-right pixel*/
				if(c < a) a = c;
				if(b < a) a = b;
				minsum[oft]=energy[oft]+a;
				if(++oft>=end)break;
				a = minsum[oft+1-width];/*get the above-right pixel*/
				if(a < b) b = a;
				if(c < b) b = c;
				minsum[oft]=energy[oft]+b;
				if(++oft>=end)break;
				b = minsum[oft+1-width];/*get the above-right pixel*/
				if(b < c) c = b;
				if(a < c) c = a;
				minsum[oft]=energy[oft]+c;
			}
		}
		
		/* find minimal seam */
		float *last = minsum + width*(height-1);
		int sp=0;
		float enr = INFINITY;
		for(i=0;i<curwid;i++){
			if(enr > last[i]){
				sp = i;
				enr = last[i];
			}
		}
		
		//printf("\nselected %d/%d at %d enr %f",k,seams,sp,enr);
		for(j=height-1;;){
			/* modify bitmap, removing one pixel at seam */
			int start = (width*j + sp);
			int len = (curwid-sp);
			void *dest = bmp_buffer + start*pixsiz;
			void *src = bmp_buffer + (start+1)*pixsiz;
			memcpy(scrap,src,len*pixsiz);
			memcpy(dest,scrap,len*pixsiz);
			/* modify energy in parallel, recomputing values where removal
			of seam will cause gradient to change */
			dest = energy + start;
			src = energy + start + 1;
			memcpy(scrap,src,len*sizeof(float));
			memcpy(dest,scrap,len*sizeof(float));
			
			if(sp>0)energy[width*j+sp-1]=gradient(sp-1,j);
			if(sp<curwid)energy[width*j+sp]=gradient(sp,j);
			if(sp<curwid-1)energy[width*j+sp+1]=gradient(sp+1,j);
			/* go up the seam */
			j--;
			if(j==0)break;
			int up = -1;
			float min = -INFINITY;
			if(sp!=0)min = minsum[width*j+sp-1];
			float cng = minsum[width*j+sp];
			if(min>cng)min = cng, up = 0;
			float rng = minsum[width*j+sp+1];
			if(min>rng)up = 1;
			sp += up;
		}
		
	}
	free(minsum);
	free(energy);
	
	printf("\noutputting carved...\n");
	
	struct jpeg_compress_struct cinfo;
	FILE *outfile=fopen("out.jpeg","w");
	row_stride = (width-seams)*pixsiz;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo,outfile);
	cinfo.image_width = width-seams;
	cinfo.image_height = height;
	cinfo.input_components = pixsiz;
	cinfo.in_color_space = dinfo.out_color_space;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo,100,TRUE);
	jpeg_start_compress(&cinfo,TRUE);
	while((j=cinfo.next_scanline) < height){
		unsigned char *row_pointer[1];
		row_pointer[0] = bmp_buffer + j*width*pixsiz;
		(void)jpeg_write_scanlines(&cinfo,row_pointer,1);
	}
	jpeg_finish_compress(&cinfo);
	if(1){
	}
	return EXIT_SUCCESS;
}
