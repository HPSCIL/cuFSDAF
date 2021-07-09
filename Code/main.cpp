/**************************************************************************
* cuFSDAF v1.1
* Author:		Huan Gao
* E-mail:		ghcug14@cug.edu.cn
* update date:2021/07/09
* Update history:
  1. 20210709 Update codes handling background values
	(1) Add a new parameter to point the images for masking fusion results
	(2) update the function for removing background values
***************************************************************************
* NOTE: this algorithm can ONLY be used for EDUCATIONAL and SCIENTIFIC
* purposes, NO COMMERCIAL usages are allowed unless the author is
* contacted and a permission is granted
***************************************************************************/


#include <stdio.h>
#include <string>
#include <time.h> 
#include "FSDAF.h"

int main(int argc, char* argv[])
{
	/* Parameters setting */
	parameter *p1 = new parameter[1];
	parseParameters(argv[1], p1);
	/************************************************************	
	 Dimension of threads for CUDA, bigger number for better GPU.
	 Recommended number: 128, 256, 512.
	 If not sure, do not change it.
	************************************************************/
	dim3 iDimGrid(512);
	dim3 iDimBlock(128);
	p1->dimGrid = iDimGrid;
	p1->dimBlock = iDimBlock;

	/* Read image size */
	size_t ns, nl, nb;
	readImgSize(p1->InputF1, ns, nl, nb);
	printf("Input images have %u x %u pixels and %u bands.\n", ns, nl, nb);

	/* Adaptive Domain-decompisition */
	int* location_block;
	size_t nBlock = 0, nNeighborWidth = 0;
	adaptiveDomainDecomposition(ns, nl, nb, p1, location_block, nBlock, nNeighborWidth);
	printf("\nThere are %u blocks: ", nBlock);

	/* For each sub-domain */
	short* fineImg2 = new short[ns * nl * nb]();
	for (size_t iblock = 0; iblock < nBlock; iblock++)
	{
		/* Location including neighborhood */
		int location_block_neighbor[4];
		locationCompute(location_block, iblock, p1->scale_factor, ns, nl, nNeighborWidth, location_block_neighbor);
		
		/* Sub-domain size: ns_block * nl_block */
		int ns_block = location_block_neighbor[1] - location_block_neighbor[0] + 1;			//columns of each block
		int nl_block = location_block_neighbor[3] - location_block_neighbor[2] + 1;			//rows of each block
		printf("\n\nBlock %d has %d x %d pixels.\n", iblock + 1, ns_block, nl_block);
		printf("The index of matching fine pixels: [%4d, %4d], [%4d, %4d]\n", location_block_neighbor[0], location_block_neighbor[1], location_block_neighbor[2], location_block_neighbor[3]);

		float *fineImg2_block = new float[ns_block * nl_block * nb]();
		blockComputing(p1, ns, nl, nb, iblock, location_block, location_block_neighbor, fineImg2_block);
		
		/* Result summary */
		resultSummary(ns, nl, nb, iblock, location_block, location_block_neighbor, fineImg2_block, fineImg2);

		delete[]fineImg2_block;
		fineImg2_block = NULL;
	}

	/* Result output */
	writeImg(p1->InputC2, fineImg2);

	delete[] location_block;
	delete[] fineImg2;
	delete[] p1;
	location_block = NULL;
	fineImg2 = NULL;
	p1 = NULL;

	return 0;
}
