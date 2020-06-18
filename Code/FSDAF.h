#pragma once

#include "cuda_runtime.h"

struct parameter
{
	char InputF1[200];
	char InputC1[200];
	char InputC2[200];
	char InputClass[200];
	size_t w = 20;								//window's size for searching similar pixels
	size_t num_similar_pixel = 20;				//number of choosen pixels in final calculation
	size_t num_pure = 100;						//number of most purest coarse pixels in each class selected fro change calculation
	float DN_min = 0.0;							//set the range of DN(Digital Number) value of the image, If byte, 0 and 255
	float DN_max = 10000.0;
	size_t scale_factor = 16;					//one coarse resolution has scale_factor * scale_factor fine resolutions
	int background = 0;							//value of background pixels. 0 means that pixels will be considered as background if one of its bands = 0
	size_t background_band = 1;					//which band with value = background indicating background pixels.

	size_t IDWSearchRadius = 2 * scale_factor;	//search radius for IDW
	int IDWPower = 2;							//power of distance in IDW

	dim3 dimGrid;
	dim3 dimBlock;

};

int parseParameters(char *fname, parameter* par);

int readImgSize(char argv[], size_t &ns, size_t &nl, size_t &nb);

int writeImg(char InputC2[], short* Img);

int adaptiveDomainDecomposition(int ns, int nl, int nb, parameter* p1, int * &location_block, size_t &nBlock, size_t &nNeighborWidth);

int locationCompute(int *location_block, int iblock, int scale_factor, int ns, int nl, int nNeighborWidth, int location_neighbor[4]);

int resultSummary(int ns, int nl, int nb, int iblock, int* location_block, int location_block_neighbor[4], float* fineImg2_block, short* fineImg2);

int blockComputing(parameter* p1, int ns, int nl, int nb, int iblock, int* location_block, int location_block_neighbor[4], float* fineImg2_block);

/* CUDA version of some steps */
int getBlockLine(parameter* p1, int &nL, int ns, int nb, int num_similar_pixel);

int cuGetHetIndex(parameter* p1, int scale_d, int ns_block, int nl_block, short* L1_class, float* het_index);

int cuInterpolate_IDW(parameter* p1, int ns_block, int nl_block, int nb, int ns_c, int nl_c, int *col_c, int *row_c, float *coarse_c2, float* L2_IDW);

int cuFinalCalculation(parameter* p1, float* fine2, float *FineImg1, float* CoarseImg1, float* CoarseImg2, float *change_21, float *D_D_all, float* similar_th,
	size_t ns_block, size_t nl_block, size_t nb_block);
