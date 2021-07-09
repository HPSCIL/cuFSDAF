
#include <algorithm>
#include <fstream>
#include "gdal_priv.h"
#include "ogr_core.h"
#include "ogr_spatialref.h"
#include "dataanalysis.h"

#include <iostream>
using namespace std;

#include "FSDAF.h"


int getPerPixSize(size_t &mnDatalength, GDALDataType datatype)
{
	switch (datatype)
	{
	case GDT_Byte:
		mnDatalength = sizeof(unsigned char);
		break;
	case GDT_UInt16:
		mnDatalength = sizeof(unsigned short);
		break;
	case GDT_Int16:
		mnDatalength = sizeof(short);
		break;
	case GDT_UInt32:
		mnDatalength = sizeof(unsigned int);
		break;
	case GDT_Int32:
		mnDatalength = sizeof(int);
		break;
	case GDT_Float32:
		mnDatalength = sizeof(float);
		break;
	case GDT_Float64:
		mnDatalength = sizeof(double);
		break;
	default:
		printf("WriteImg(): Unknown data type!\n");
		return -1;
	}
	return 0;
}

/* Read and parse parameters from file 'fname' and save into par */
int parseParameters(char *fname, parameter* par)
{
	int   i, k, total_pairs = 0;
	char  buffer[1000] = "\0";
	const char  *label = NULL;

	const char  *tokenptr = NULL;
	char readpath[1000];
	std::string argName;
	std::string t;
	char  *separator = "= ,";
	FILE  *in;
	if ((in = fopen(fname, "r")) == NULL)
	{
		printf("Can't open input %s\n", fname);
		return -1;
	}
	fscanf(in, "%s", buffer);
	if (strcmp(buffer, "cuFSDAF_PARAMETER_START") != 0)
	{
		printf("This is not a valid input file\n");
		return -1;
	}
	int nn = 0;
	while (1)
	{
		nn++;
		if (nn > 1000)
		{
			std::cerr << "Parameter file should end using \"cuFSDAF_PARAMETER_END\".\n";
			break;
		}
		memset(buffer, 0, 1000);
		if (fgets(buffer, 1000, in) == NULL)
			continue;
		tokenptr = strtok(buffer, separator);
		label = tokenptr;
		if (strcmp(label, "cuFSDAF_PARAMETER_END") == 0)
			break;
		if (strcmp(label, "#") == 0) continue;
		while (tokenptr != NULL)
		{
			tokenptr = strtok(NULL, separator);
			if (strcmp(label, "IN_F1_NAME") == 0)
				sscanf(tokenptr, "%s", par->InputF1);
			else if (strcmp(label, "IN_C1_NAME") == 0)
				sscanf(tokenptr, "%s", par->InputC1);
			else if (strcmp(label, "IN_C2_NAME") == 0)
				sscanf(tokenptr, "%s", par->InputC2);
			else if (strcmp(label, "IN_F1_CLASS_NAME") == 0)
				sscanf(tokenptr, "%s", par->InputClass);
			else if (strcmp(label, "W") == 0)
				par->w = atoi(tokenptr);
			else if (strcmp(label, "NUM_SIMILAR_PIXEL") == 0)
				par->num_similar_pixel = atoi(tokenptr);
			else if (strcmp(label, "NUM_PURE") == 0)
				par->num_pure = atoi(tokenptr);
			else if (strcmp(label, "DN_MIN") == 0)
				par->DN_min = atof(tokenptr);
			else if (strcmp(label, "DN_MAX") == 0)
				par->DN_max = atof(tokenptr);
			else if (strcmp(label, "SCALE_FACTOR") == 0)
				par->scale_factor = atoi(tokenptr);
			else if (strcmp(label, "BACKGROUND") == 0)
				par->background = atof(tokenptr);
			else if (strcmp(label, "BACKGROUND_IMG") == 0)
				par->background_img = atoi(tokenptr);
			else if (strcmp(label, "BACKGROUND_BAND") == 0)
				par->background_band = atoi(tokenptr);
			else if (strcmp(label, "IDW_SEARCH_RADIUS") == 0)
				par->IDWSearchRadius = atoi(tokenptr);
			else if (strcmp(label, "IDW_POWER") == 0)
				par->IDWPower = atoi(tokenptr);
			label = tokenptr;
		}
	}
	fclose(in);
	return 0;
}

/* Calculate Standard Deviation */
double stddev(double *dif, int a_size)
{
	double totaldif = 0;
	double result = 0;
	double totaldev = 0;
	for (int i = 0; i < a_size; i++)
		totaldif += dif[i];
	for (int i = 0; i < a_size; i++)
		totaldev += pow((dif[i] - totaldif / double(a_size - 1)), 2);
	result = sqrt(totaldev / double(a_size));
	return result;
}

/* Read size and projection information of image InputF1 */
int readImgSize(char InputF1[], size_t &ns, size_t &nl, size_t &nb)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");	//support Chinese path
	GDALDataset *ImgBef = (GDALDataset*)GDALOpen(InputF1, GA_ReadOnly);
	if (ImgBef == NULL)
	{
		printf("load error! The file name \" %s \" may be wrong.\n", InputF1);
		exit(0);
	}
	ns = ImgBef->GetRasterXSize();			//samples of each line
	nl = ImgBef->GetRasterYSize();			//number of lines
	nb = ImgBef->GetRasterCount();			//band number
	GDALClose(ImgBef);
	ImgBef = NULL;
	return 0;
}

/* Write fusion image */
template<typename  T> int writeOut(T _val, GDALDataset* pwrite_out, short* Img, unsigned char* mpData)
{
	int ns = pwrite_out->GetRasterXSize();
	int nl = pwrite_out->GetRasterYSize();
	int nb = pwrite_out->GetRasterCount();
	for (int i = 0; i < nl; i++)
	{
		for (int j = 0; j < ns; j++)
		{
			for (int k = 0; k <nb; k++)
			{
				_val = Img[k * ns * nl + i * ns + j];
				size_t nloc = (k * ns * nl + i * ns + j)* sizeof(T);
				memcpy(mpData + nloc, &_val, sizeof(T));
			}
		}
	}
	return 0;
}

/* Write fusion image */
int writeImg(char InputC2[], short* Img)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");	//support Chinese path
	GDALDataset *ImgBef = (GDALDataset*)GDALOpen(InputC2, GA_ReadOnly);
	if (ImgBef == NULL)
	{
		printf("load error! The file name may be wrong.\n");
		exit(0);
	}
	GDALDataType datatype = ImgBef->GetRasterBand(1)->GetRasterDataType();
	size_t ns = ImgBef->GetRasterXSize();
	size_t nl = ImgBef->GetRasterYSize();
	size_t nb = ImgBef->GetRasterCount();
	double* fNoDataV = new double[nb]();
	for (size_t i = 0; i < nb; i++)
		fNoDataV[i] = ImgBef->GetRasterBand(i + 1)->GetNoDataValue();
	double geoTransform[6];
	char projref[2048];

	ImgBef->GetGeoTransform(geoTransform);
	const char *sProRef = ImgBef->GetProjectionRef();
	strcpy(projref, sProRef);
	GDALClose(ImgBef);
	ImgBef = NULL;

	char str3[100];
	char str2[] = "_cuFSDAF";
	char *cPtr = NULL;
	char *separator = ".";
	cPtr = strtok(InputC2, separator);
	cPtr = strtok(NULL, separator);
	sprintf(str3, "%s%s", InputC2, str2);
	if (cPtr != NULL)
		sprintf(str3, "%s.%s", str3, cPtr);

	GDALDriver* mpoDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* ImgOut = mpoDriver->Create(str3, ns, nl, nb, datatype, NULL);
	if (ImgOut == NULL)
	{
		printf("CGDALWrite::init : Create poDataset Failed.\n");
		GDALClose(ImgOut);
		ImgOut = NULL;
		return -1;
	}
	ImgOut->SetGeoTransform(geoTransform);
	ImgOut->SetProjection(projref);
	for (size_t i = 0; i < nb; i++)
		ImgOut->GetRasterBand(i + 1)->SetNoDataValue(fNoDataV[i]);
	size_t nPerPixSize = 0;
	getPerPixSize(nPerPixSize, datatype);
	size_t mnDatalength = ns * nl * nb * nPerPixSize;
	unsigned char* mpData = new unsigned char[mnDatalength];
	memset(mpData, 0, mnDatalength);
	int _iNs = ImgOut->GetRasterXSize();
	unsigned char _ucV = 0;
	unsigned short _usV = 0;
	short _sV = 0;
	unsigned int _uiV = 0;
	int _iV = 0;
	float _fV = 0.0;
	double _dV = 0.0;
	switch (datatype)
	{
	case 1:
		writeOut(_ucV, ImgOut, Img, mpData);
		break;
	case 2:
		writeOut(_usV, ImgOut, Img, mpData);
		break;
	case 3:
		writeOut(_sV, ImgOut, Img, mpData);
		break;
	case 4:
		writeOut(_uiV, ImgOut, Img, mpData);
		break;
	case 5:
		writeOut(_iV, ImgOut, Img, mpData);
		break;
	case 6:
		writeOut(_fV, ImgOut, Img, mpData);
		break;
	case 7:
		writeOut(_dV, ImgOut, Img, mpData);
		break;
	default:
		break;
	}
	if (ImgOut != NULL && mpData != NULL)
	{
		ImgOut->RasterIO(GF_Write, 0, 0, ns, nl, \
			mpData, ns, nl, datatype, nb, 0, 0, 0, 0);
		ImgOut->FlushCache();
	}
	printf("\nwrite \"%s\" successfully!\n", str3);

	delete[]fNoDataV;
	fNoDataV = NULL;
	GDALClose(ImgBef);
	GDALClose(ImgOut);
	ImgBef = NULL;
	ImgOut = NULL;

	return 0;
}

/* remove background values based on background band */
int imgNoback(int ns, int nl, int nb, int background_band, float background, float* Img)
{
	std::vector <int> _background_whole(ns * nl, 0);
	int _num_back = 0;
	for (int i = 0; i < ns * nl; i++)
	{
		if (abs(Img[i + ns * nl * (background_band - 1)] - background) <= 1.0e-8)
		{
			_background_whole[i] = 1;
			_num_back++;
		}
	}
	if (_num_back > 0)
	{
		for (int iband = 0; iband < nb; iband++)
		{
			double count = 0;
			for (int i = 0; i < ns * nl; i++)
				if (_background_whole[i] == 0)
					count += Img[i + ns * nl * iband];
			for (int i = 0; i < ns * nl; i++)
				if (_background_whole[i] == 1)
					Img[i + ns * nl * iband] = count / (ns * nl - _num_back);
		}
	}
	return 0;
}

/* remove background values in full image */
int imgNoback1(int ns, int nl, int nb, float background, float* Img)
{
	for (int iband = 0; iband < nb; iband++)
	{
		std::vector <int> _background_whole(ns * nl, 0);
		int _num_back = 0;
		for (int i = 0; i < ns * nl; i++)
		{
			if (abs(Img[i + ns * nl * iband] - background) <= 1.0e-8)
			{
				_background_whole[i] = 1;
				_num_back++;
			}
		}
		if (_num_back > 0)
		{
			double count = 0;
			for (int i = 0; i < ns * nl; i++)
				if (_background_whole[i] == 0)
					count += Img[i + ns * nl * iband];
			for (int i = 0; i < ns * nl; i++)
				if (_background_whole[i] == 1)
					Img[i + ns * nl * iband] = count / (ns * nl - _num_back);
		}


	}
	return 0;
}

/* extreme values correction */
int extremeCorrect(int ns, int nl, int nb, int DN_min, int DN_max, float* Img)
{
	std::vector<float> _FineImg1_sort;
	for (int i = 0; i < nl * ns * nb; i++)
		_FineImg1_sort.push_back(Img[i]);
	float data_1_4[2];		//max and min in [0.0001, 0.9999]
	for (int ib = 0; ib < nb; ib++)
	{
		std::sort(_FineImg1_sort.begin() + ns * nl * ib, _FineImg1_sort.begin() + ns * nl * (ib + 1));
		int i1, i2;
		for (int i = 0; i < ns * nl + 1; i++)
		{
			if (float(i) / (ns * nl) < 0.0001)
			{
				i1 = i;
				data_1_4[0] = _FineImg1_sort[i + ns * nl * ib];
			}
			if (float(i) / (ns * nl) < 0.9999)
			{
				data_1_4[1] = _FineImg1_sort[i + ns * nl * ib];
				i2 = i;
			}
		}
		for (size_t i = i1; i < ns * nl; i++)
		{
			if (_FineImg1_sort[i + ns * nl * ib] == data_1_4[0])
				continue;
			if (_FineImg1_sort[i + ns * nl * ib] >= DN_min)
			{
				data_1_4[0] = _FineImg1_sort[i + ns * nl * ib];
				break;
			}
		}
		for (size_t i = i2; i > 0; i--)
		{
			if (_FineImg1_sort[i + ns * nl * ib] == data_1_4[1])
				continue;
			if (_FineImg1_sort[i + ns * nl * ib] <= DN_max)
			{
				data_1_4[1] = _FineImg1_sort[i + ns * nl * ib];
				break;
			}
		}
		for (int i = 0; i < ns * nl; i++)
		{
			Img[i + ns * nl * ib] = std::max(Img[i + ns * nl * ib], data_1_4[0]);
			Img[i + ns * nl * ib] = std::min(Img[i + ns * nl * ib], data_1_4[1]);
		}
	}
	return 0;
}

/* decompose images into sub-domains adaptively */
int adaptiveDomainDecomposition(int ns, int nl, int nb, parameter* p1, int * &location_block, size_t &nBlock, size_t &nNeighborWidth)
{
	int nSimilarPixels = p1->num_similar_pixel;
	nNeighborWidth = std::max(p1->scale_factor, p1->IDWSearchRadius);
	nNeighborWidth = (nNeighborWidth + p1->w + p1->scale_factor - 1) / p1->scale_factor;
	/* Each sub-domain has nL lines at most */
	int nL = 0;
	getBlockLine(p1, nL, ns, nb, nSimilarPixels);
	if (nL > 2 * nNeighborWidth * p1->scale_factor)
		nL = nL - 2 * nNeighborWidth * p1->scale_factor;
	else{
		nNeighborWidth = 0;
		printf("Caution in adaptiveDomainDecomposition().\n");
		printf("Very few lines in one block, results in edge of block may be noncontinuous.\n");
		printf("Consider replacing GPU with larger video memory.\n");
	}
	int nl_block = std::min(nL, nl);
	int ns_block = ns;
	/* number of sub-domains */
	int n_ns = (ns - 1 + ns_block) / ns_block;
	int n_nl = (nl - 1 + nl_block) / nl_block;
	nBlock = n_ns * n_nl;
	/* Location of sub-domains */
	location_block = new int[4 * nBlock]();	//locations of start point and end point of each block
	for (int i_ns = 0; i_ns < n_ns; i_ns++)
	{
		for (int i_nl = 0; i_nl < n_nl; i_nl++)
		{
			location_block[0 + (i_ns + n_ns * i_nl) * 4] = i_ns * ns_block;
			location_block[1 + (i_ns + n_ns * i_nl) * 4] = std::min(ns - 1, (i_ns + 1) * ns_block - 1);
			location_block[2 + (i_ns + n_ns * i_nl) * 4] = i_nl * nl_block;
			location_block[3 + (i_ns + n_ns * i_nl) * 4] = std::min(nl - 1, (i_nl + 1) * nl_block - 1);
		}
	}
	return 0;
}

/* Compute locations of each sub-domain including neighborhood */
int locationCompute(int *location_block, int iblock, int scale_factor, int ns, int nl, int nNeighborWidth, int location_neighbor[4])
{
	location_neighbor[0] = std::max(0, location_block[iblock * 4] - nNeighborWidth * scale_factor);
	location_neighbor[1] = std::min(ns - 1, location_block[1 + iblock * 4] + nNeighborWidth * scale_factor);
	location_neighbor[2] = std::max(0, location_block[2 + iblock * 4] - nNeighborWidth * scale_factor);
	location_neighbor[3] = std::min(nl - 1, location_block[3 + iblock * 4] + nNeighborWidth * scale_factor);

	return 0;
}

/* Summarize results of each sub-domain */
int resultSummary(int ns, int nl, int nb, int iblock, int* location_block, int location_block_neighbor[4], float* fineImg2_block, short* fineImg2)
{
	int ns_block = location_block_neighbor[1] - location_block_neighbor[0] + 1;
	int nl_block = location_block_neighbor[3] - location_block_neighbor[2] + 1;
	for (size_t i_ns = 0; i_ns < (location_block[1 + iblock * 4] - location_block[0 + iblock * 4]) + 1; i_ns++)
	{
		for (size_t i_nl = 0; i_nl < (location_block[3 + iblock * 4] - location_block[2 + iblock * 4]) + 1; i_nl++)
		{
			for (size_t i_nb = 0; i_nb < nb; i_nb++)
			{
				size_t _blockIdx = i_ns + (location_block[0 + iblock * 4] - location_block_neighbor[0]) + (i_nl + location_block[2 + iblock * 4] - location_block_neighbor[2]) * ns_block + i_nb * ns_block * nl_block;
				size_t _globalIdx = location_block[0 + iblock * 4] + i_ns + (location_block[2 + iblock * 4] + i_nl) * ns + i_nb * ns * nl;
				fineImg2[_globalIdx] = (short)round(fineImg2_block[_blockIdx]);
			}
		}
	}
	return 0;
}

/* Read input images from pointer *p according to data type */
template<typename  T> void pointerConversion(T& nV, unsigned char *p, int datatype)
{
	switch (datatype)
	{
	case 0:
		printf("Unknown data type!\n");
		break;
	case 1:

		nV = *(unsigned char*)p;
		break;
	case 2:
		nV = *(unsigned short*)p;
		break;
	case 3:
		nV = *(short*)p;
		break;
	case 4:
		nV = *(unsigned int*)p;
		break;
	case 5:
		nV = *(int*)p;
		break;
	case 6:
		nV = *(float*)p;
		break;
	case 7:
		nV = *(double*)p;
		break;
	default:
		nV = *(unsigned char*)p;
		break;
	}
}


template<typename  T>
int LoadImg(char* inputFile, int location_block_neighbor[4], T* Img_block)
{
	int ns_block = location_block_neighbor[1] - location_block_neighbor[0] + 1;			//columns of each block(including neighborhood)
	int nl_block = location_block_neighbor[3] - location_block_neighbor[2] + 1;			//rows of each block(including neighborhood)
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");	//support Chinese path
	GDALDataset *pread = (GDALDataset*)GDALOpen(inputFile, GA_ReadOnly);
	int ns = pread->GetRasterXSize();						//samples of each line
	int nl = pread->GetRasterYSize();						//lines of fineImg1
	int nb = pread->GetRasterCount();					//band number
	GDALDataType datatype = pread->GetRasterBand(1)->GetRasterDataType();
	size_t mnPerPixSize = 0;
	getPerPixSize(mnPerPixSize, datatype);
	size_t mnDatalength = ns * nl * nb * mnPerPixSize;
	unsigned char* pData = new unsigned char[mnDatalength];
	memset(pData, 0, mnDatalength);
	pread->RasterIO(GF_Read, 0, 0, ns, nl, pData, ns, nl, datatype, nb, 0, 0, 0, 0);
	for (int j = location_block_neighbor[2]; j < location_block_neighbor[3] + 1; j++)
	{
		for (int i = location_block_neighbor[0]; i < location_block_neighbor[1] + 1; i++)
		{
			for (int k = 0; k < nb; k++)
			{
				T _sV = 0;
				pointerConversion(_sV, &(pData[(k * ns * nl + j * ns + i)*mnPerPixSize]), datatype);
				Img_block[i - location_block_neighbor[0] + (j - location_block_neighbor[2]) * ns_block + k * ns_block * nl_block] = _sV;
			}
		}
	}
	GDALClose(pread);
	pread = NULL;
	delete[] pData;
	pData = NULL;
	return 0;
}

/* Read images of each block */
int readBlockImage(parameter* p1, int iblock, int &num_class, int location_block_neighbor[4],
		float* fineImg1_block, float* coarseImg1_block, float* coarseImg2_block, short* L1_class_block)
{
	LoadImg(p1->InputF1, location_block_neighbor, fineImg1_block);
	LoadImg(p1->InputC1, location_block_neighbor, coarseImg1_block);
	LoadImg(p1->InputC2, location_block_neighbor, coarseImg2_block);
	LoadImg(p1->InputClass, location_block_neighbor, L1_class_block);

	int ns_block = location_block_neighbor[1] - location_block_neighbor[0] + 1;			//columns of each block(including neighborhood)
	int nl_block = location_block_neighbor[3] - location_block_neighbor[2] + 1;			//rows of each block(including neighborhood)

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");	//support Chinese path

	GDALDataset *pread = (GDALDataset*)GDALOpen(p1->InputF1, GA_ReadOnly);
	int nb = pread->GetRasterCount();
	int nMaxClassNum = 0;
	for (int i = 0; i < ns_block * nl_block; i++)
		if (L1_class_block[i] >= nMaxClassNum)
			nMaxClassNum = L1_class_block[i];
	if (nMaxClassNum <= 0)
	{
		printf("Number of class is wrong!\n");
		exit(0);
	}
	for (size_t iclass = 0; iclass < nMaxClassNum; iclass++)
	{
		size_t num_ic = 0;
		std::vector<size_t> ind_ic;
		ind_ic.reserve(ns_block * nl_block);
		for (size_t _i = 0; _i < ns_block * nl_block; _i++)
		{
			size_t nGlbIdx = _i % ns_block + (_i / ns_block) * ns_block + (p1->background_band - 1) * ns_block * nl_block;
			if (nGlbIdx >= ns_block * nl_block * nb)
			{
				printf("\nERROR: The parameter \"background band\" may be bigger than the number of bands, Please check your input.\n");
				exit(0);
			}
			if ( (L1_class_block[_i] == iclass + 1) && (abs(fineImg1_block[nGlbIdx] - p1->background) >= 1.0e-8) )
			{
				ind_ic.push_back(_i);
				num_ic++;
			}
		}
		if (num_ic > 0)
		{
			for (size_t _ic = 0; _ic < num_ic; _ic++)
			{
				L1_class_block[ind_ic[_ic]] = num_class + 1;
			}
			num_class++;
		}
	}
	if (num_class >= 0)
	{
		printf("load \" %s \" successfully!\n", p1->InputF1);
		printf("load \" %s \" successfully!\n", p1->InputC1);
		printf("load \" %s \" successfully!\n", p1->InputC2);
		printf("load \" %s \" successfully!\n", p1->InputClass);
	}
	GDALClose(pread);
	pread = NULL;
	return 0;
}

/* Calculate fraction of each class in each coarse pixel */
int getFraction(int ns_c, int nl_c, int num_class, int ns_block, int nl_block, int scale_factor, short* L1_class, std::vector<float> &Fraction1)
{
	for (size_t ic = 0; ic < ns_c; ic++)
	{
		for (size_t jc = 0; jc < nl_c; jc++)
		{
			int *num_ic = new int[num_class]();

			size_t I1 = ic * scale_factor;
			size_t I2 = std::min((ic + 1) * scale_factor, (size_t)ns_block);
			size_t J1 = jc * scale_factor;
			size_t J2 = std::min((jc + 1) * scale_factor, (size_t)nl_block);
			size_t num_c = (I2 - I1) * (J2 - J1);
			for (size_t i = I1; i < I2; i++)
				for (size_t j = J1; j < J2; j++)
					for (int _ic = 0; _ic < num_class; _ic++)
						if (_ic + 1 == L1_class[i])
							num_ic[_ic]++;
			float Fraction_total = 0;
			for (int iclass = 0; iclass < num_class; iclass++)
			{
				Fraction1[ic + jc * ns_c + iclass * ns_c * nl_c] = float(num_ic[iclass]) / float(num_c);
				Fraction_total += Fraction1[ic + jc * ns_c + iclass * ns_c * nl_c];
			}
			if (Fraction_total <= 0.999)
				for (int iclass = 0; iclass < num_class; iclass++)
					Fraction1[ic + jc * ns_c + iclass * ns_c * nl_c] = 0;

			delete[]num_ic;
			num_ic = NULL;
		}
	}
	return 0;
}

/* Get change rate of each class from t1 to t2*/
int getChangeRate(int ns_c, int nl_c, int nb, int num_class, int num_pure, std::vector<float> coarse_c1, float *coarse_c2,
	std::vector<float> Fraction1, std::vector<double> &c_rate)
{
	for (int ib = 0; ib < nb; ib++)
	{
		float *X_matrix0 = new float[num_class * num_pure * num_class]();
		float *Y_matrix0 = new float[num_class * num_pure]();

		int ii = 0;
		for (int ic = 0; ic < num_class; ic++)
		{
			std::vector<int> order(ns_c * nl_c, 0);
			std::multimap<float, int>mmap1;
			for (int i = 0; i < ns_c * nl_c; i++)
			{
				mmap1.insert(std::make_pair(Fraction1[i + ic * ns_c * nl_c], i));
			}
			std::multimap<float, int>::iterator mmap1_it = mmap1.begin();
			for (int i = 0; i < mmap1.size(); i++)
			{
				order[i] = mmap1_it->second;
				mmap1_it++;
			}
			std::reverse(order.begin(), order.end());
			int num_f = 0;
			for (int i = 0; i < ns_c * nl_c; i++)
				if (Fraction1[i + ic * ns_c * nl_c] > 0.01)
					num_f++;
			//make sure the ordered samples(coarse pixel)' class weight > 0.1 and the sum < num_pure
			int num_pure1 = std::min(num_f, num_pure);
			std::multimap<float, int>mmap_change_c;
			for (int i = 0; i < num_pure1; i++)
			{
				//pick first num_pure1 samples and sort them by change_c
				mmap_change_c.insert(std::make_pair(coarse_c2[order[i] + ib * ns_c * nl_c] - coarse_c1[order[i] + ib * ns_c * nl_c], order[i]));
			}
			std::multimap<float, int>::iterator mmap_change_c_it = mmap_change_c.begin();
			std::vector<int> sortIndex(num_pure1, 0);				//index of ordered samples
			std::vector<float> sortIndexValue(num_pure1, 0);		//the change value of each ordered samples(change_c[])
			for (int i = 0; i < num_pure1; i++)
			{
				sortIndex[i] = mmap_change_c_it->second;
				sortIndexValue[i] = mmap_change_c_it->first;
				mmap_change_c_it++;
			}
			float data_1_4[2];		//index of point in change_c[] which locates in 0.25 and 0.75
			for (int i = 0; i < num_pure1 + 1; i++)
			{
				if (float(i) / num_pure1 <= 0.25)
					data_1_4[0] = sortIndexValue[i];
				if (float(i) / num_pure1 <= 0.75)
					data_1_4[1] = sortIndexValue[i];
			}
			//regard samples which locate in [0.25,0.75] as class type nochange
			int num_nonc = 0;
			for (int i = 0; i < num_pure1; i++)
				if (sortIndexValue[i] >= data_1_4[0] && sortIndexValue[i] <= data_1_4[1])
					num_nonc++;
			std::vector<int> ind_nochange(num_nonc, 0);
			size_t count_temp = 0;
			for (int i = 0; i < num_pure1; i++)
				if (sortIndexValue[i] >= data_1_4[0] && sortIndexValue[i] <= data_1_4[1])
					ind_nochange[count_temp++] = i;
			for (int i = 0; i < num_nonc; i++)
				Y_matrix0[0 + ii + i] = sortIndexValue[ind_nochange[i]];
			for (int icc = 0; icc < num_class; icc++)
				for (int i = 0; i < num_nonc; i++)
					X_matrix0[icc + (i + ii) * num_class] = Fraction1[sortIndex[ind_nochange[i]] + icc * ns_c * nl_c];
			ii = ii + num_nonc;
		}
		std::vector<float> X_matrix1(num_class * ii, 0);
		std::vector<float> Y_matrix1(ii, 0);
		for (int i_ii = 0; i_ii < ii; i_ii++)
		{
			Y_matrix1[i_ii] = Y_matrix0[i_ii];
			for (int ic = 0; ic < num_class; ic++)
				X_matrix1[ic + i_ii * num_class] = X_matrix0[ic + i_ii * num_class];
		}

		float *XY_matrix = new float[(num_class + 1) * ii]();
		for (int i_ii = 0; i_ii < ii; i_ii++)
		{
			for (int ic = 0; ic < num_class; ic++)
			{
				XY_matrix[ic + i_ii * (num_class + 1)] = X_matrix1[ic + i_ii * num_class];
			}
			XY_matrix[num_class + i_ii * (num_class + 1)] = Y_matrix1[i_ii];
		}
		//linear regression
		std::string strOut = "[";
		for (int i_ii = 0; i_ii < ii - 1; i_ii++)
		{
			strOut += "[";
			for (int ic = 0; ic < num_class; ic++)
			{
				strOut += std::to_string(XY_matrix[ic + i_ii * (num_class + 1)]) + ",";
			}
			strOut += std::to_string(XY_matrix[num_class + i_ii * (num_class + 1)]) + "],";
		}
		strOut += "[";
		for (int ic = 0; ic < num_class; ic++)
		{
			strOut += std::to_string(XY_matrix[ic + (ii - 1) * (num_class + 1)]) + ",";
		}
		strOut += std::to_string(XY_matrix[num_class + (ii - 1) * (num_class + 1)]) + "]]";
		const char *xy_chr = strOut.c_str();
		alglib::real_2d_array xy = xy_chr;
		alglib::ae_int_t info;
		alglib::ae_int_t nvars;
		alglib::linearmodel model;
		alglib::lrreport rep;
		alglib::real_1d_array c;
		lrbuildz(xy, xy.rows(), num_class, info, model, rep);
		lrunpack(model, c, nvars);
		for (int i = 0; i < num_class; i++)
		{
			c_rate[i + ib * num_class] = c[i];
		}

		delete[] X_matrix0;
		delete[] Y_matrix0;
		delete[] XY_matrix;
		X_matrix0 = NULL;
		Y_matrix0 = NULL;
		XY_matrix = NULL;
	}
	return 0;
}

/* Distributes residuals of two predictions */
int residualRedistribution(int ns_c, int nl_c, int ns_block, int nl_block, int nb, int scale_factor, std::vector<float> coarse_c2_p, std::vector<float> fine_c1, std::vector<float> coarse_c1, float *coarse_c2,
	float* L2_Interpolation, std::vector<double> L2_1, float *het_index, float* change_21)
{
	std::vector<double> predict_change_c(ns_c * nl_c * nb, 0);
	std::vector<double> real_change_c(ns_c * nl_c * nb, 0);
	std::vector<double> change_R(ns_c * nl_c * nb, 0);
	for (size_t i = 0; i < ns_c * nl_c * nb; i++)
	{
		predict_change_c[i] = coarse_c2_p[i] - fine_c1[i];
		real_change_c[i] = coarse_c2[i] - coarse_c1[i];
		change_R[i] = real_change_c[i] - predict_change_c[i];
	}
	for (size_t ic = 0; ic < ns_c; ic++)
	{
		for (size_t jc = 0; jc < nl_c; jc++)
		{
			size_t I1 = ic * scale_factor;
			size_t I2 = std::min((ic + 1) * scale_factor, (size_t)ns_block);
			size_t J1 = jc * scale_factor;
			size_t J2 = std::min((jc + 1) * scale_factor, (size_t)nl_block);
			size_t num_ii = (I2 - I1) * (J2 - J1);
			std::vector<int> ind_c(num_ii, 0);
			size_t count_temp = 0;
			for (size_t i = I1; i < I2; i++)
				for (size_t j = J1; j < J2; j++)
					ind_c[count_temp++] = i + j * ns_block;

			for (size_t ib = 0; ib < nb; ib++)
			{
				double diff_change = change_R[ic + jc * ns_c + ib * ns_c * nl_c];
				std::vector<double> w_change_tps(num_ii, 0);
				for (size_t i = 0; i < num_ii; i++)
				{
					w_change_tps[i] = L2_Interpolation[ind_c[i] + ib * ns_block * nl_block] - L2_1[ind_c[i] + ib * ns_block * nl_block];
				}
				if (diff_change <= 0)
				{
					int num_noc = 0;
					for (size_t i = 0; i < num_ii; i++)
						if (w_change_tps[i] > 0)
							num_noc++;
					std::vector<int> ind_noc(num_noc, 0);	//if w_change_tps > 0, save the index
					size_t count_temp = 0;
					for (size_t i = 0; i < num_ii; i++)
						if (w_change_tps[i] > 0)
						{
							ind_noc[count_temp] = i;
							count_temp++;
						}
					if (num_noc > 0)
					{
						for (size_t i = 0; i < num_noc; i++)
						{
							w_change_tps[ind_noc[i]] = 0;
						}
					}
				}
				else
				{
					int num_noc = 0;
					for (size_t i = 0; i < num_ii; i++)
						if (w_change_tps[i] < 0)
							num_noc++;
					std::vector<int> ind_noc(num_noc, 0);		//if w_change_tps < 0, save the index
					size_t count_temp = 0;
					for (size_t i = 0; i < num_ii; i++)
						if (w_change_tps[i] < 0)
						{
							ind_noc[count_temp] = i;
							count_temp++;
						}
					if (num_noc > 0)
					{
						for (size_t i = 0; i < num_noc; i++)
						{
							w_change_tps[ind_noc[i]] = 0;
						}
					}
				}
				for (size_t i = 0; i < num_ii; i++)
				{
					w_change_tps[i] = abs(w_change_tps[i]);
				}
				std::vector<double> w_unform(num_ii, 0);

				for (size_t i = 0; i < num_ii; i++)
					w_unform[i] = abs(diff_change);
				std::vector<double> w_change(num_ii, 0);
				double total_w_change = 0;
				for (size_t i = 0; i < num_ii; i++)
				{
					w_change[i] = w_change_tps[i] * het_index[ind_c[i]] + w_unform[i] * (1 - het_index[ind_c[i]]) + 0.000000001;
					total_w_change += w_change[i];
				}
				for (size_t i = 0; i < num_ii; i++)
					w_change[i] = w_change[i] / (total_w_change / num_ii);		//weight normalization
				//change extreme weight to the mean of every weight
				int num_extrem = 0;
				for (size_t i = 0; i < num_ii; i++)
				{
					if (w_change[i] > 10)
					{
						num_extrem++;
					}
				}
				std::vector<int> ind_extrem(num_extrem, 0);
				size_t temp_count = 0;
				for (size_t i = 0; i < num_ii; i++)
				{
					if (w_change[i] > 10)
					{
						ind_extrem[temp_count] = i;
						temp_count++;
					}
				}
				total_w_change = 0;
				for (size_t i = 0; i < num_ii; i++)
					total_w_change += w_change[i];
				if (num_extrem > 0)
				{
					for (size_t i = 0; i < num_extrem; i++)
					{
						w_change[ind_extrem[i]] = total_w_change / num_ii;
					}
				}
				total_w_change = 0;
				for (size_t i = 0; i < num_ii; i++)
				{
					total_w_change += w_change[i];
				}
				for (size_t i = 0; i < num_ii; i++)
				{
					w_change[i] = w_change[i] / (total_w_change / num_ii);
					change_21[ind_c[i] + ib * ns_block * nl_block] = w_change[i] * diff_change;
				}
			}
		}
	}
	return 0;
}

/* Get change rate of each pixel from t1 to t2*/
int getChange_21(float * change_21_block, int ns_block, int nl_block, int nb, parameter* p1, int num_class,
	float *fineImg1_block, float *coarseImg1_block, float *coarseImg2_block, short *L1_class)
{
	int w = p1->w;
	int num_pure = p1->num_pure;
	float background = p1->background;
	int background_band = p1->background_band;
	int num_similar_pixel = p1->num_similar_pixel;
	float DN_max = p1->DN_max;
	float DN_min = p1->DN_min;
	//the index of Image
	int ns_c = int((ns_block + p1->scale_factor - 1) / p1->scale_factor);				//coarse pixels of each row
	int nl_c = int((nl_block + p1->scale_factor - 1) / p1->scale_factor);				//coarse pixels of each column
	//the index of interpolation samples
	std::vector<int> row_ind(ns_block * nl_block, 0);
	std::vector<int> col_ind(ns_block * nl_block, 0);
	for (int i = 0; i < ns_block; i++)
		for (int j = 0; j < nl_block; j++)
		{
			col_ind[i + j * ns_block] = i;
			row_ind[i + j * ns_block] = j;
		}
	//resample images to coarse image
	std::vector<float> fine_c1(ns_c * nl_c * nb, 0);
	std::vector<float> coarse_c1(ns_c * nl_c * nb, 0);
	float* coarse_c2 = new float[ns_c * nl_c * nb]();
	int *row_c = new int[ns_c * nl_c]();					//coarse pixel's row index in fine images
	int *col_c = new int[ns_c * nl_c]();					//coarse pixel's col index in fine images
	for (size_t ic = 0; ic < ns_c; ic++)
	{
		for (size_t jc = 0; jc < nl_c; jc++)
		{
			size_t total_row = 0;
			size_t total_col = 0;
			float * total_fine1 = new float[nb]();
			float * total_coarse1 = new float[nb]();
			float * total_coarse2 = new float[nb]();


			size_t I1 = ic * p1->scale_factor;
			size_t I2 = std::min((ic + 1) * p1->scale_factor, (size_t)ns_block);
			size_t J1 = jc * p1->scale_factor;
			size_t J2 = std::min((jc + 1) * p1->scale_factor, (size_t)nl_block);
			size_t Count_total_rc = (I2 - I1) * (J2 - J1);

			for (size_t i = I1; i < I2; i++)
			{
				for (size_t j = J1; j < J2; j++)
				{
					total_col += i;
					total_row += j;
					for (size_t _iband = 0; _iband < nb; _iband++)
					{
						total_fine1[_iband] += fineImg1_block[i + j * ns_block + _iband * ns_block * nl_block];
						total_coarse1[_iband] += coarseImg1_block[i + j * ns_block + _iband * ns_block * nl_block];
						total_coarse2[_iband] += coarseImg2_block[i + j * ns_block + _iband * ns_block * nl_block];
					}
				}
			}

			row_c[ic + jc * ns_c] = total_row / Count_total_rc;
			col_c[ic + jc * ns_c] = total_col / Count_total_rc;
			for (int _iband = 0; _iband < nb; _iband++)
			{
				fine_c1[ic + jc * ns_c + _iband * ns_c * nl_c] = total_fine1[_iband] / Count_total_rc;
				coarse_c1[ic + jc * ns_c + _iband * ns_c * nl_c] = total_coarse1[_iband] / Count_total_rc;
				coarse_c2[ic + jc * ns_c + _iband * ns_c * nl_c] = total_coarse2[_iband] / Count_total_rc;
			}
			delete[]total_fine1;
			delete[]total_coarse1;
			delete[]total_coarse2;
			total_fine1 = NULL;
			total_coarse1 = NULL;
			total_coarse2 = NULL;
		}
	}

	//step2:get the fraction of each class in coarse pixel at t1 
	std::vector<float> Fraction1(ns_c * nl_c * num_class, 0);
	getFraction(ns_c, nl_c, num_class, ns_block, nl_block, p1->scale_factor, L1_class, Fraction1);

	//get the homogeneity of each fine pixel   
	float *het_index = new float[ns_block * nl_block]();
	int scale_d = w;
	cuGetHetIndex(p1, scale_d, ns_block, nl_block, L1_class, het_index);

	//change rate of each class from t1 to t2
	std::vector<double> c_rate(num_class * nb, 0);
	getChangeRate(ns_c, nl_c, nb, num_class, num_pure, coarse_c1, coarse_c2, Fraction1, c_rate);

	//step4: prediction without considering spatial change
	std::vector<double> L2_1(ns_block * nl_block * nb, 0);
	for (int i = 0; i < ns_block * nl_block * nb; i++)
		L2_1[i] = fineImg1_block[i];
	for (int i = 0; i < ns_block * nl_block; i++)
		for (int ic = 0; ic < num_class; ic++)
			if (ic + 1 == L1_class[i])
				for (int ib = 0; ib < nb; ib++)
					L2_1[i + ib * ns_block * nl_block] = L2_1[i + ib * ns_block * nl_block] + c_rate[ic + ib * num_class];

	//resample L2_1[] to coarse image coarse_c2_p[]
	std::vector<float> coarse_c2_p(ns_c * nl_c * nb, 0);
	for (size_t i_ns = 0; i_ns < ns_c; i_ns++)
	{
		for (size_t i_nl = 0; i_nl < nl_c; i_nl++)
		{
			for (size_t ib = 0; ib < nb; ib++)
			{
				float total_temp = 0;
				size_t I1 = i_ns * p1->scale_factor;
				size_t I2 = std::min((i_ns + 1) * p1->scale_factor, (size_t)ns_block);
				size_t J1 = i_nl * p1->scale_factor;
				size_t J2 = std::min((i_nl + 1) * p1->scale_factor, (size_t)nl_block);
				size_t count_temp = (I2 - I1) * (J2 - J1);

				for (int k = I1; k < I2; k++)
					for (int l = J1; l < J2; l++)
						total_temp += L2_1[k + l * ns_block + ib * ns_block * nl_block];
				coarse_c2_p[i_ns + i_nl * ns_c + ib * ns_c * nl_c] = total_temp / count_temp;
			}
		}
	}
	std::vector<float> min_allow_c2(nb, 0);
	std::vector<float> max_allow_c2(nb, 0);
	for (int ib = 0; ib < nb; ib++)
	{
		float min_allow0 = DN_max;
		float max_allow0 = -DN_max;
		for (int i = 0; i < ns_c * nl_c; i++)
		{
			if (coarse_c2[i + ib * ns_c * nl_c] < min_allow0)
				min_allow0 = coarse_c2[i + ib * ns_c * nl_c];
			if (coarse_c2[i + ib * ns_c * nl_c] > max_allow0)
				max_allow0 = coarse_c2[i + ib * ns_c * nl_c];
		}
		for (int i = 0; i < ns_block * nl_block; i++)
		{
			if (L2_1[i + ib * ns_block * nl_block] < min_allow0)
				min_allow0 = L2_1[i + ib * ns_block * nl_block];
			if (L2_1[i + ib * ns_block * nl_block] > max_allow0)
				max_allow0 = L2_1[i + ib * ns_block * nl_block];
		}
		if (min_allow0 < DN_min)
			min_allow_c2[ib] = DN_min;
		else
			min_allow_c2[ib] = min_allow0;
		if (max_allow0 > DN_max)
			max_allow_c2[ib] = DN_max;
		else
			max_allow_c2[ib] = max_allow0;
	}

	/* predict L2 using IDW interpolation */
	float* L2_IDW = new float[ns_block * nl_block * nb]();
	cuInterpolate_IDW(p1, ns_block, nl_block, nb, ns_c, nl_c, col_c, row_c, coarse_c2, L2_IDW);
	for (size_t i = 0; i < ns_block * nl_block * nb; i++)
	{
		if (abs(L2_IDW[i] - 0) <= 1.0e-8)
		{
			printf("Some pixels cannot search neighboring pixels in IDW.\nUse larger IDWSearchRadius or check input images.\n");
			printf("i: %d\t Band: %d\t Index: [%d, %d].\n", i, i / (ns_block * nl_block), i % ns_block, (i % (ns_block * nl_block) / ns_block));
			return -1;
		}

	}

	//residual redistribution
	residualRedistribution(ns_c, nl_c, ns_block, nl_block, nb, p1->scale_factor, coarse_c2_p, fine_c1, coarse_c1, coarse_c2, L2_IDW, L2_1, het_index, change_21_block);

	//second prediction from TPS 
	std::vector<float> fine2_2(ns_block * nl_block * nb, 0);
	for (int i = 0; i < ns_block * nl_block * nb; i++)
		fine2_2[i] = L2_1[i] + change_21_block[i];
	for (int ib = 0; ib < nb; ib++)
	{
		for (int i = 0; i < ns_block * nl_block; i++)
		{
			if (fine2_2[i + ib * ns_block * nl_block] < min_allow_c2[ib])
				fine2_2[i + ib * ns_block * nl_block] = min_allow_c2[ib];
			if (fine2_2[i + ib * ns_block * nl_block] > max_allow_c2[ib])
				fine2_2[i + ib * ns_block * nl_block] = max_allow_c2[ib];
		}
	}
	for (int i = 0; i < ns_block * nl_block * nb; i++)
		change_21_block[i] = fine2_2[i] - fineImg1_block[i];

	delete[] row_c;
	delete[] col_c;
	delete[] coarse_c2;
	delete[]het_index;
	delete[]L2_IDW;

	row_c = NULL;
	col_c = NULL;
	coarse_c2 = NULL;
	het_index = NULL;
	L2_IDW = NULL;

	return 0;
}

/* Get thresholds of similarity for searching similar pixels */
int getSimilarTh(int iblock, int ns_block, int nl_block, int scale_factor, int nb, int num_class,
	int nSubLines, int nSubSamples, int nSubNgbWidth, int location_block_neighbor[4], float* fineImg1_block, float* similar_th)
{
	size_t n_ns_sub1 = (ns_block / scale_factor - 1 + nSubSamples) / nSubSamples;				//each row has n_ns subs
	size_t n_nl_sub1 = (nl_block / scale_factor - 1 + nSubLines) / nSubLines;					//each column has n_nl subs
	int *location_sub1 = new int[4 * n_ns_sub1 * n_nl_sub1]();									//location of start point and end point

	for (int i_ns = 0; i_ns < n_ns_sub1; i_ns++)
	{
		for (int i_nl = 0; i_nl < n_nl_sub1; i_nl++)
		{
			location_sub1[0 + (i_ns + n_ns_sub1 * i_nl) * 4] = std::max(0, i_ns * nSubSamples * scale_factor);
			location_sub1[1 + (i_ns + n_ns_sub1 * i_nl) * 4] = std::min(ns_block - 1, (i_ns + 1) * nSubSamples * scale_factor - 1);
			location_sub1[2 + (i_ns + n_ns_sub1 * i_nl) * 4] = std::max(0, i_nl * nSubLines * scale_factor);
			location_sub1[3 + (i_ns + n_ns_sub1 * i_nl) * 4] = std::min(nl_block - 1, (i_nl + 1) * nSubLines * scale_factor - 1);
		}
	}

	for (size_t isub = 0; isub < n_ns_sub1 * n_nl_sub1; isub++)
	{
		//location including neighborhood
		int location_sub_neighbor1[4] = { location_sub1[0 + isub * 4] , location_sub1[1 + isub * 4] , location_sub1[2 + isub * 4] , location_sub1[3 + isub * 4] };
		location_sub_neighbor1[0] = std::max(0, location_sub_neighbor1[0] - nSubNgbWidth * scale_factor);
		location_sub_neighbor1[1] = std::min(ns_block - 1, location_sub_neighbor1[1] + nSubNgbWidth * scale_factor);
		location_sub_neighbor1[2] = std::max(0, location_sub_neighbor1[2] - nSubNgbWidth * scale_factor);
		location_sub_neighbor1[3] = std::min(nl_block - 1, location_sub_neighbor1[3] + nSubNgbWidth * scale_factor);

		int ns_sub = location_sub_neighbor1[1] - location_sub_neighbor1[0] + 1;	//columns of each sub(including neighborhood)
		int nl_sub = location_sub_neighbor1[3] - location_sub_neighbor1[2] + 1;	//rows of each sub(including neighborhood)
		float* fineImg1_sub = new float[ns_sub * nl_sub * nb]();
		for (size_t i = location_sub_neighbor1[0]; i < location_sub_neighbor1[1] + 1; i++)
		{
			for (size_t j = location_sub_neighbor1[2]; j < location_sub_neighbor1[3] + 1; j++)
			{
				int nBlkIdx = i + j * ns_block;
				int nSubIdx = i - location_sub_neighbor1[0] + (j - location_sub_neighbor1[2]) * ns_sub;
				for (size_t k = 0; k < nb; k++)
					fineImg1_sub[nSubIdx + k * ns_sub * nl_sub] = fineImg1_block[nBlkIdx + k * ns_block * nl_block];
			}
		}
		for (size_t i_nb = 0; i_nb < nb; i_nb++)
		{
			std::vector<float> similar_th_sub(nb, 0);
			double *fine1_iband = new double[ns_sub * nl_sub]();
			for (size_t i = 0; i < ns_sub * nl_sub; i++)
				fine1_iband[i] = fineImg1_sub[i + i_nb * ns_sub * nl_sub];

			similar_th_sub[i_nb] = stddev(fine1_iband, ns_sub * nl_sub) * 2.0 / num_class;
			delete[]fine1_iband;
			fine1_iband = NULL;
			for (size_t i_nl = 0; i_nl < (location_sub1[3 + isub * 4] - location_sub1[2 + isub * 4] + 1); i_nl++)
			{
				for (size_t i_ns = 0; i_ns < (location_sub1[1 + isub * 4] - location_sub1[0 + isub * 4] + 1); i_ns++)
				{
					size_t nGlbIdx = location_sub1[0 + isub * 4] + i_ns + (location_sub1[2 + isub * 4] + i_nl) * ns_block + i_nb * ns_block * nl_block;
					similar_th[nGlbIdx] = similar_th_sub[i_nb];
				}
			}
		}
		delete[]fineImg1_sub;
		fineImg1_sub = NULL;
	}
	return 0;
}

/* Mask final images using the background band in one image (0, 1, 2: F1, C1, C2) */
int getBackgArea(int ns, int nl, int nb, int iBackgImg, int background_band, float background, float* F1, float* C1, float* C2, unsigned char* BackgArea)
{
	float * Img = NULL;
	switch (iBackgImg)
	{
	case 0:
		Img = F1;
		break;
	case 1:
		Img = C1;
		break;
	case 2:
		Img = C2;
		break;
	default:
		break;
	}
	std::vector <int> _background_whole(ns * nl, 0);
	int _num_back = 0;
	for (int i = 0; i < ns * nl; i++)
	{
		if (abs(Img[i + ns * nl * (background_band - 1)] - background) <= 1.0e-8)
		{
			BackgArea[i] = 1;
		}
	}
	return 0;
}


int finalImgMask(int ns, int nl, int nb, float fBackgVal, unsigned char* BackgArea, float* Img)
{
	for (int iband = 0; iband < nb; iband++)
		for (int i = 0; i < ns * nl; i++)
			if (BackgArea[i] == 1)
				Img[i + ns * nl * iband] = fBackgVal;
	return 0;
}

/* Compute in each sub-domain */
int blockComputing(parameter* p1, int ns, int nl, int nb, int iblock,
	int* location_block, int location_block_neighbor[4], float* fineImg2_block)
{
	int ns_block = location_block_neighbor[1] - location_block_neighbor[0] + 1;			//columns of each block
	int nl_block = location_block_neighbor[3] - location_block_neighbor[2] + 1;			//rows of each block
	float* fineImg1_block = new float[ns_block * nl_block * nb]();
	float* coarseImg1_block = new float[ns_block * nl_block * nb]();
	float* coarseImg2_block = new float[ns_block * nl_block * nb]();
	short* L1_class_block = new short[ns_block * nl_block]();

	/* Read block Images */
	int num_class = 0;
	readBlockImage(p1, iblock, num_class, location_block_neighbor, fineImg1_block, coarseImg1_block, coarseImg2_block, L1_class_block);
	printf("\nPlease confirm the data is read correctly:\n");
	printf("fineImg1[0] = %.3f  coarseImg1[0] = %.3f  coarseImg2[0] = %.3f  L1_class[0] = %d num_class = %d\n", fineImg1_block[0], coarseImg1_block[0], coarseImg2_block[0], L1_class_block[0], num_class);

	///* Replace background values by average values */
	//imgNoback(ns_block, nl_block, nb, p1->background_band, p1->background, fineImg1_block);
	//imgNoback(ns_block, nl_block, nb, p1->background_band, p1->background, coarseImg1_block);
	//imgNoback(ns_block, nl_block, nb, p1->background_band, p1->background, coarseImg2_block);
	
	// get background locations
	unsigned char* ucharBackgArea = new unsigned char[ns_block * nl_block]();
	getBackgArea(ns_block, nl_block, nb, p1->background_img, p1->background_band, p1->background, fineImg1_block, coarseImg1_block, coarseImg2_block, ucharBackgArea);


	// Replace all background vals in every band by average values
	imgNoback1(ns_block, nl_block, nb, p1->background, fineImg1_block);
	imgNoback1(ns_block, nl_block, nb, p1->background, coarseImg1_block);
	imgNoback1(ns_block, nl_block, nb, p1->background, coarseImg2_block);
	
	/* Correct extreme values */
	extremeCorrect(ns_block, nl_block, nb, p1->DN_min, p1->DN_max, fineImg1_block);

	/* Compute changes from t1 to t2 */
	float* change_21_block = new float[ns_block * nl_block * nb]();
	getChange_21(change_21_block, ns_block, nl_block, nb, p1, num_class, fineImg1_block, coarseImg1_block, coarseImg2_block, L1_class_block);

	/* Mitigate errors using neighborhood information */
	float* D_D_all = new float[(p1->w * 2 + 1) * (p1->w * 2 + 1)]();
	for (int i = 0; i < p1->w * 2 + 1; i++)
		for (int j = 0; j < p1->w * 2 + 1; j++)
			D_D_all[i + j * (p1->w * 2 + 1)] = sqrtf(pow((i - 20.0), 2) + pow((j - 20.0),2));
	float *similar_th = new float[ns_block * nl_block * nb]();
	getSimilarTh(iblock, ns_block, nl_block, p1->scale_factor, nb, num_class, ns_block, nl_block, 0, location_block_neighbor, fineImg1_block, similar_th);
	cuFinalCalculation(p1, fineImg2_block, fineImg1_block, coarseImg1_block, coarseImg2_block, change_21_block, D_D_all, similar_th,
		ns_block, nl_block, nb);

	/* Masking F2 using backg area */
	finalImgMask(ns_block, nl_block, nb, p1->background, ucharBackgArea, fineImg2_block);


	delete[] fineImg1_block;
	delete[] coarseImg1_block;
	delete[] coarseImg2_block;
	delete[] L1_class_block;
	delete[] change_21_block;
	delete[] D_D_all;
	delete[] similar_th;
	delete[] ucharBackgArea;
	
	fineImg1_block = NULL;
	coarseImg1_block = NULL;
	coarseImg2_block = NULL;
	L1_class_block = NULL;
	change_21_block = NULL;
	D_D_all = NULL;
	similar_th = NULL;
	ucharBackgArea = NULL;

	return 0;
}