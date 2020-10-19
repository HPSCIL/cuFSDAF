#include <iostream>
#include <stdio.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_functions.h"

#include "FSDAF.h"

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}


__host__ __device__
void sortSimilarPixel(float arr[], size_t id[], size_t len, size_t num_cand)
{
	for (size_t i = 0; i < num_cand; i++)
	{
		float fMaxArr = arr[i];
		size_t nID = id[i];
		size_t _iCount = i;
		for (size_t j = i + 1; j < len; j++)
		{
			if (arr[j] < fMaxArr)
			{
				fMaxArr = arr[j];
				nID = id[j];
				_iCount = j;
			}
		}
		arr[_iCount] = arr[i];
		id[_iCount] = id[i];
		arr[i] = fMaxArr;
		id[i] = nID;
	}
}


int getBlockLine(parameter* p1, int &nL, int ns, int nb, int num_similar_pixel)
{
	size_t a = 0;

	cudaDeviceProp deviceProp;
	int deviceCount;
	CHECK(cudaGetDeviceCount(&deviceCount));
	for (size_t i = 0; i < deviceCount; i++)
	{
		CHECK(cudaGetDeviceProperties(&deviceProp, i));
		a = deviceProp.totalGlobalMem;
	}

	int nThread = p1->dimGrid.x * p1->dimBlock.x;			//size of 1 cycle
	int nWindowSize = (p1->w + 1 + p1->w) * (p1->w + 1 + p1->w);	//size of window, and w is the length of the searching window

	//Constant memory of final computing
	int nConstantMemory = sizeof(float) * nWindowSize + sizeof(int) * (p1->w + 1) * (p1->w + 1) * nWindowSize * 2
		+ sizeof(int) * nWindowSize * nThread * 2 + sizeof(float) * nWindowSize * nThread + sizeof(double) * num_similar_pixel * nThread;

	//Memory of each line
	int nSizeOfEachLine = (sizeof(float) * 3 + sizeof(short) * 3) * ns * nb;

	nL = (a * 0.8 - nConstantMemory) / nSizeOfEachLine;

	return 0;
}

__global__
void cuHet(int iCycle, int nThread, int scale_d, int ns_block, int nl_block, short* L1_class_d, float* het_index_d)
{
	size_t iThread = blockIdx.x * blockDim.x + threadIdx.x;

	if (iThread < nThread)
	{
		if (iThread + iCycle * nThread < ns_block * nl_block)
		{
			size_t i = (iThread + iCycle * nThread) % ns_block;
			size_t j = (iThread + iCycle * nThread) / ns_block;

			size_t ai = 0;
			size_t bi = 0;
			size_t aj = 0;
			size_t bj = 0;

			if (i < scale_d) ai = 0; else ai = i - scale_d;
			if (i + scale_d > ns_block - 1) bi = ns_block - 1; else bi = i + scale_d;
			if (j < scale_d) aj = 0; else aj = j - scale_d;
			if (j + scale_d > nl_block - 1) bj = nl_block - 1; else bj = j + scale_d;

			size_t num_sameclass = 0;
			for (size_t l = aj; l < bj; l++)
				for (size_t k = ai; k < bi; k++)
					if (L1_class_d[k + l * ns_block] == L1_class_d[i + j * ns_block])
						num_sameclass++;
			het_index_d[i + j * ns_block] = float(num_sameclass) / ((bi - ai)*(bj - aj));

		}
	}
}

int cuGetHetIndex(parameter* p1, int scale_d, int ns_block, int nl_block, short* L1_class, float* het_index)
{
	dim3 dimGrid = p1->dimGrid;
	dim3 dimBlock = p1->dimBlock;

	size_t nThread = dimGrid.x * dimBlock.x;			//size of 1 cycle

	// Malloc the memory in device(GPU)
	short *L1_class_d;
	float *het_index_d;

	cudaMalloc((void**)&L1_class_d, sizeof(short) * ns_block * nl_block);
	cudaMalloc((void**)&het_index_d, sizeof(float) * ns_block * nl_block);

	// Copy data from host to device (CPU to GPU)
	cudaMemcpy(L1_class_d, L1_class, sizeof(short) * ns_block * nl_block, cudaMemcpyHostToDevice);
	cudaMemcpy(het_index_d, het_index, sizeof(float) * ns_block * nl_block, cudaMemcpyHostToDevice);

	for (int iCycle = 0; iCycle * dimGrid.x * dimBlock.x < ns_block * nl_block; iCycle++)
	{
		cuHet << <dimGrid, dimBlock >> > (iCycle, nThread, scale_d, ns_block, nl_block, L1_class_d, het_index_d);
		CHECK(cudaThreadSynchronize());
	}

	cudaMemcpy(het_index, het_index_d, sizeof(float) *  ns_block * nl_block, cudaMemcpyDeviceToHost);

	cudaFree(L1_class_d);
	cudaFree(het_index_d);

	CHECK(cudaDeviceReset());

	return 0;
}

__global__
void cuIDW(int icycle, int nThread, int scale_factor, int background_band, int ns_block, int nl_block, int nb, int ns_c, int nl_c, int IDWSearchRadius, int IDWPower, int *col_c_D, int *row_c_D, float *coarse_c2_D, float* L2_IDW_D)
{
	size_t iThread = blockIdx.x * blockDim.x + threadIdx.x;

	if (iThread < nThread)
	{
		size_t i = iThread + icycle * nThread;
		if (i < ns_block * nl_block)
		{
			for (size_t ib = 0; ib < nb; ib++)
			{
				double wSum = 0;	// weight sum

				//coordinate indexes of the minimum enclosing rectangle for the searching circle 
				int index0 = (i % ns_block + scale_factor / 2);
				int index1 = (i % ns_block + IDWSearchRadius - (scale_factor - 1) / 2) / scale_factor;
				int index2 = (i / ns_block + scale_factor / 2);
				int index3 = (i / ns_block + IDWSearchRadius - (scale_factor - 1) / 2) / scale_factor;

				if (index0 < IDWSearchRadius) index0 = 0; else index0 = (index0 - IDWSearchRadius) / scale_factor;
				if (index1 > ns_c - 1) index1 = ns_c - 1;
				if (index2 < IDWSearchRadius) index2 = 0; else index2 = (index2 - IDWSearchRadius) / scale_factor;
				if (index3 > nl_c - 1) index3 = nl_c - 1;

				for (size_t _j = index2; _j <= index3; _j++)
				{
					if (abs(wSum - 1) <= 1.0e-8)
						break;
					for (size_t _i = index0; _i <= index1; _i++)
					{
						size_t iKnownPoint = _i + _j * ns_c;
						double d = sqrt(pow((double)(col_c_D[iKnownPoint] - (int(i) % ns_block)), 2) + pow((double)(row_c_D[iKnownPoint] - (int(i) / ns_block)), 2));
						d = rsqrt(d);

						if (abs(d) >= 1.0e-8)		//if sqrt((knownPoint[nKnownPointIndex].x - (i % ns_block))^2 + (y - yi)^2 ) <= IDWSearchRadius, ...
						{
							if (d <= IDWSearchRadius)
							{
								wSum += pow(d, -IDWPower);
								L2_IDW_D[i + ib * ns_block * nl_block] += pow(d, -IDWPower) * coarse_c2_D[iKnownPoint + ib * ns_c * nl_c];
							}
						}
						else
						{
							L2_IDW_D[i + ib * ns_block * nl_block] = coarse_c2_D[iKnownPoint + ib * ns_c * nl_c];
							wSum = 1;
							break;
						}
					}
				}

				if (abs(wSum) > 1.0e-8)
					L2_IDW_D[i + ib * ns_block * nl_block] /= wSum;
				else
				{
					L2_IDW_D[i + ib * ns_block * nl_block] = 0;
					printf("IDW=0");
				}
			}
		}
	}
}

int cuInterpolate_IDW(parameter* p1, int ns_block, int nl_block, int nb, int ns_c, int nl_c, int *col_c, int *row_c, float *coarse_c2, float* L2_IDW)
{
	dim3 dimGrid = p1->dimGrid;
	dim3 dimBlock = p1->dimBlock;
	size_t nThread = dimGrid.x * dimBlock.x;			//size of 1 cycle

	// Malloc the memory in device(GPU)
	int *col_c_d, *row_c_d;
	float *coarse_c2_d, *L2_IDW_d;
	cudaMalloc((void**)&col_c_d, sizeof(int) * ns_c * nl_c);
	cudaMalloc((void**)&row_c_d, sizeof(int) * ns_c * nl_c);
	cudaMalloc((void**)&coarse_c2_d, sizeof(float) * ns_c * nl_c * nb);
	cudaMalloc((void**)&L2_IDW_d, sizeof(float) * ns_block * nl_block * nb);

	// Copy data from host to device (CPU to GPU)
	cudaMemcpy(col_c_d, col_c, sizeof(int) * ns_c * nl_c, cudaMemcpyHostToDevice);
	cudaMemcpy(row_c_d, row_c, sizeof(int) * ns_c * nl_c, cudaMemcpyHostToDevice);
	cudaMemcpy(coarse_c2_d, coarse_c2, sizeof(float) * ns_c * nl_c * nb, cudaMemcpyHostToDevice);
	cudaMemcpy(L2_IDW_d, L2_IDW, sizeof(float) * ns_block * nl_block * nb, cudaMemcpyHostToDevice);
	
	for (size_t iCycle = 0; iCycle * dimGrid.x * dimBlock.x < ns_block * nl_block; iCycle++)
	{
		cuIDW << <dimGrid, dimBlock >> > (iCycle, nThread, p1->scale_factor, p1->background_band, ns_block, nl_block, nb, ns_c, nl_c, p1->IDWSearchRadius, p1->IDWPower, col_c_d, row_c_d, coarse_c2_d, L2_IDW_d);
		CHECK(cudaThreadSynchronize());
	}

	cudaMemcpy(L2_IDW, L2_IDW_d, sizeof(float) *  ns_block * nl_block * nb, cudaMemcpyDeviceToHost);

	cudaFree(col_c_d);
	cudaFree(row_c_d);
	cudaFree(coarse_c2_d);
	cudaFree(L2_IDW_d);

	CHECK(cudaDeviceReset());

	return 0;
}

__global__
void cuFine2_1(float* FineImg1_d, float *CoarseImg1_d, float *CoarseImg2_d, float *change_21_d, float *fine2_d,
	float *D_D_all_d, float *similar_th_d, size_t *col_wind_d, size_t *row_wind_d, size_t *positionCand_orderDis_d, float *mmap_order_dis_d, double *D_D_cand_d, size_t ns, size_t nl, int nb, int background_band,
	int background, size_t w, float DN_max, float DN_min, size_t num_similar_pixel, size_t cycle_size, size_t cycle_time)
{
	size_t i_thread = blockIdx.x * blockDim.x + threadIdx.x;

	if (i_thread < cycle_size)
	{
		size_t i = (i_thread + cycle_size * cycle_time) % ns;
		size_t j = (i_thread + cycle_size * cycle_time) / ns;
		if (i < ns && j < nl &&
			(CoarseImg1_d[i + j * ns + (background_band - 1) * ns * nl] - background) > 1e-6)
		{
			// searching range
			size_t ai = max(i, w) - w;
			size_t bi = min(ns - 1, i + w);
			size_t aj = max(j, w) - w;
			size_t bj = min(nl - 1, j + w);
			size_t ci = i - ai;
			size_t cj = j - aj;
			size_t nI = w, nJ = w;
			if (i - w < 0) nI = i;
			if (j - w < 0) nJ = j;
			if (i + w > ns - 1) nI = ns - 1 - i;
			if (j + w > nl - 1) nJ = nl - 1 - j;
		
			//search similar pixels within the window
			size_t nWindowSize = (w + 1 + w) * (w + 1 + w);
			size_t number_cand0 = 0;
			for (size_t k = 0; k < (bi - ai + 1) * (bj - aj + 1); k++)
			{
				positionCand_orderDis_d[k + i_thread * nWindowSize] = 1;
				number_cand0++;
				for (size_t iband = 0; iband < nb; iband++)
				{
					size_t nGlbIdx = (k % (bi - ai + 1)) + ai + (k / (bi - ai + 1)) * ns + aj * ns + iband * ns * nl;
					size_t nGlbIdx_c = ci + ai + cj * ns + aj * ns + iband * ns * nl;
					if (abs(FineImg1_d[nGlbIdx] - FineImg1_d[nGlbIdx_c]) >= similar_th_d[nGlbIdx])
					{
						positionCand_orderDis_d[k + i_thread * nWindowSize] = 0;
						number_cand0--;
						break;
					}
				}
			}

			size_t order_dis_count = 0;
			for (size_t k = 0; k < (bi - ai + 1) * (bj - aj + 1); k++)
			{
				if (positionCand_orderDis_d[k + i_thread * nWindowSize] != 0)
				{
					double similar_cand_k = 0.0;
					for (size_t ib = 0; ib < nb; ib++)
					{
						size_t nGlbIdx = (k % (bi - ai + 1)) + ai + (k / (bi - ai + 1)) * ns + aj * ns + ib * ns * nl;
						size_t nGlbIdx_c = ci + ai + cj * ns + aj * ns + ib * ns * nl;
						similar_cand_k += abs(FineImg1_d[nGlbIdx] - FineImg1_d[nGlbIdx_c]) / (double)FineImg1_d[nGlbIdx_c];
					}
					mmap_order_dis_d[order_dis_count + i_thread * nWindowSize] = similar_cand_k;
					positionCand_orderDis_d[order_dis_count + i_thread * nWindowSize] = k;
					++order_dis_count;
				}
			}
			
			size_t number_cand = min(number_cand0, num_similar_pixel);

			sortSimilarPixel(mmap_order_dis_d + i_thread * nWindowSize, positionCand_orderDis_d + i_thread * nWindowSize, number_cand0, number_cand);

			if ((bi - ai + 1)*(bj - aj + 1) < (w*2.0 + 1)*(w*2.0 + 1))
			{
				for (size_t k = 0; k < number_cand; k++)
				{
					size_t nGlbIdx = positionCand_orderDis_d[k + i_thread * nWindowSize] + (nI + nJ * (w + 1)) * (w + 1 + w) * (w + 1 + w);
					D_D_cand_d[k + num_similar_pixel * i_thread] = sqrtf(double((ci - col_wind_d[nGlbIdx]) * (ci - col_wind_d[nGlbIdx]) + (cj - row_wind_d[nGlbIdx]) * (cj - row_wind_d[nGlbIdx]))) + 0.0000001;
				}
			}
			else
			{
				for (size_t k = 0; k < number_cand; k++)
					D_D_cand_d[k + num_similar_pixel * i_thread] = D_D_all_d[positionCand_orderDis_d[k + i_thread * nWindowSize]];
			}
			for (size_t k = 0; k < number_cand; k++)
			{
				double similar_cand_ind_same_class_k = 0.0;
				for (size_t ib = 0; ib < nb; ib++)
				{
					size_t nGlbIdx = (positionCand_orderDis_d[k + i_thread * nWindowSize] % (bi - ai + 1)) + ai + (positionCand_orderDis_d[k + i_thread * nWindowSize] / (bi - ai + 1)) * ns + aj * ns + ib * ns * nl;
					size_t nGlbIdx_c = ci + ai + cj * ns + aj * ns + ib * ns * nl;
					similar_cand_ind_same_class_k += abs(FineImg1_d[nGlbIdx] - FineImg1_d[nGlbIdx_c]) / (double)FineImg1_d[nGlbIdx_c];
				}
				D_D_cand_d[k + num_similar_pixel * i_thread] = (1.0 + D_D_cand_d[k + num_similar_pixel * i_thread] / w)*(similar_cand_ind_same_class_k + 1.0);
			}
			double temp_Total_C_D = 0.0;
			for (size_t k = 0; k < number_cand; k++)
			{
				temp_Total_C_D += 1.0 / D_D_cand_d[k + num_similar_pixel * i_thread];
			}
			//predict the value (formula 24)
			for (size_t iband = 0; iband < nb; iband++)
			{
				double temp_total_weight = 0;
				for (size_t k = 0; k < number_cand; k++)
				{
					size_t glbIdx = (positionCand_orderDis_d[k + i_thread * nWindowSize] % (bi - ai + 1)) + ai + (positionCand_orderDis_d[k + i_thread * nWindowSize] / (bi - ai + 1)) * ns + aj * ns + iband * ns * nl;
					temp_total_weight += (1.0 / D_D_cand_d[k + num_similar_pixel * i_thread] / temp_Total_C_D) * change_21_d[glbIdx];
				}
				//int _nIdxFine2 = i + j * ns + iband * ns * nl;
				size_t _nIdx = i + j * ns + iband * ns * nl;
				fine2_d[_nIdx] = FineImg1_d[_nIdx] + temp_total_weight;
				//revise the abnormal prediction

				if (fine2_d[_nIdx] <= DN_min || fine2_d[_nIdx] >= DN_max)
				{
					fine2_d[_nIdx] = FineImg1_d[_nIdx] + (CoarseImg2_d[_nIdx] - CoarseImg1_d[_nIdx]);
					if (fine2_d[_nIdx] > DN_max) fine2_d[_nIdx] = DN_max;
					if (fine2_d[_nIdx] < DN_min) fine2_d[_nIdx] = DN_min;
				}
			}
		}
	}
}

int cuFinalCalculation(parameter* p1, float* fine2, float *FineImg1, float* CoarseImg1, float* CoarseImg2, float *change_21, float *D_D_all, float* similar_th,
	size_t ns_block, size_t nl_block, size_t nb)
{
	size_t w = p1->w;
	dim3 dimGrid = p1->dimGrid;
	dim3 dimBlock = p1->dimBlock;
	size_t imgSize = ns_block * nl_block * nb;				//size of image, ns, nl, nb indicate the width, length, and the bands of img
	size_t nThread = dimGrid.x * dimBlock.x;				//size of 1 cycle
	size_t nWindowSize = (w + 1 + w) * (w + 1 + w);			//size of window, and w is the length of the searching window

	size_t *col_wind = new size_t[(w + 1) * (w + 1) * (w + 1 + w) * (w + 1 + w)]();
	size_t *row_wind = new size_t[(w + 1) * (w + 1) * (w + 1 + w) * (w + 1 + w)]();
	for (size_t i = 0; i < w + 1; i++)
	{
		for (size_t j = 0; j < w + 1; j++)
		{
			size_t ai, bi, aj, bj;
			if (i > w) ai = i - w; else ai = 0;
			bi = i + w;
			if (j > w) aj = j - w; else aj = 0;
			bj = j + w;

			for (size_t k = 0; k < bi - ai + 1; k++)
			{
				for (size_t l = 0; l < bj - aj + 1; l++)
				{
					col_wind[k + l * (bi - ai + 1) + (i + j * (w + 1)) * (w + 1 + w) * (w + 1 + w)] = k;
					row_wind[k + l * (bi - ai + 1) + (i + j * (w + 1)) * (w + 1 + w) * (w + 1 + w)] = l;
				}
			}
		}
	}

	float *CoarseImg1_d, *CoarseImg2_d;
	float *fine2_d;
	float *FineImg1_d, *change_21_d;
	float *D_D_all_d, *similar_th_d;
	size_t *col_wind_d, *row_wind_d;
	size_t *positionCand_orderDis_d;
	float *mmap_order_dis_d;
	double *D_D_cand_d;

	// Malloc the memory in device(GPU)
	cudaMalloc((void**)&fine2_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&FineImg1_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&CoarseImg1_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&CoarseImg2_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&change_21_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&similar_th_d, sizeof(float) * imgSize);
	cudaMalloc((void**)&D_D_all_d, sizeof(float) * nWindowSize);
	cudaMalloc((void**)&col_wind_d, sizeof(size_t) * (w + 1) * (w + 1) * nWindowSize);
	cudaMalloc((void**)&row_wind_d, sizeof(size_t) * (w + 1) * (w + 1) * nWindowSize);
	cudaMalloc((void**)&positionCand_orderDis_d, sizeof(size_t) * nWindowSize * nThread);
	cudaMalloc((void**)&mmap_order_dis_d, sizeof(float) * nWindowSize * nThread);
	cudaMalloc((void**)&D_D_cand_d, sizeof(double) * p1->num_similar_pixel * nThread);


	// Copy data from host to device (CPU to GPU)
	cudaMemcpy(fine2_d, fine2, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(FineImg1_d, FineImg1, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(CoarseImg1_d, CoarseImg1, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(CoarseImg2_d, CoarseImg2, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(change_21_d, change_21, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(similar_th_d, similar_th, sizeof(float) * imgSize, cudaMemcpyHostToDevice);
	cudaMemcpy(D_D_all_d, D_D_all, sizeof(float) * nWindowSize, cudaMemcpyHostToDevice);
	cudaMemcpy(col_wind_d, col_wind, sizeof(size_t) * (w + 1) * (w + 1) * nWindowSize, cudaMemcpyHostToDevice);
	cudaMemcpy(row_wind_d, row_wind, sizeof(size_t) * (w + 1) * (w + 1) * nWindowSize, cudaMemcpyHostToDevice);

	for (size_t iCycle = 0; iCycle * dimGrid.x * dimBlock.x < ns_block * nl_block; iCycle++)
	{
		size_t nCyl = (ns_block * nl_block / dimGrid.x / dimBlock.x + 20 - 1) / 20;
		if (nCyl != 0 && (iCycle % nCyl) == 0)
		{
			printf("\r*");
			for (size_t i = 0; i <= (iCycle / nCyl); i++)
				printf("-");
			for (size_t i = 0; i <= (ns_block * nl_block / (dimGrid.x * dimBlock.x) / nCyl) - (iCycle / nCyl); i++)
				printf(" ");
			printf("*");
		}
		cuFine2_1 << <dimGrid, dimBlock >> > (FineImg1_d, CoarseImg1_d, CoarseImg2_d, change_21_d, fine2_d, D_D_all_d, similar_th_d,
			col_wind_d, row_wind_d, positionCand_orderDis_d, mmap_order_dis_d, D_D_cand_d, ns_block, nl_block, nb, p1->background_band,
			p1->background, w, p1->DN_max, p1->DN_min, p1->num_similar_pixel, nThread, iCycle);

		CHECK(cudaThreadSynchronize());
	}


	// Copy result from device to host (GPU to CPU)
	cudaMemcpy(fine2, fine2_d, sizeof(float) * imgSize, cudaMemcpyDeviceToHost);

	cudaFree(fine2_d);
	cudaFree(FineImg1_d);
	cudaFree(CoarseImg1_d);
	cudaFree(CoarseImg2_d);
	cudaFree(change_21_d);
	cudaFree(D_D_all_d);
	cudaFree(similar_th_d);
	cudaFree(col_wind_d);
	cudaFree(row_wind_d);
	cudaFree(positionCand_orderDis_d);
	cudaFree(mmap_order_dis_d);
	cudaFree(D_D_cand_d);

	CHECK(cudaDeviceReset());

	delete[]col_wind;
	delete[]row_wind;
	col_wind = NULL;
	row_wind = NULL;

	return 0;
}


