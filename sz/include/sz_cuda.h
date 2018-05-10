#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "sz.h"

unsigned char * SZ_compress_float_1D_MDQ_cuda(float *oriData, size_t r1, double realPrecision, size_t * comp_size);

unsigned char * SZ_compress_float_2D_MDQ_cuda(float *oriData, size_t r1, size_t r2, double realPrecision, size_t * comp_size);

unsigned char * SZ_compress_float_3D_MDQ_cuda(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size);

void decompressDataSeries_float_1D_cuda(float** data, size_t r1, unsigned char* comp_data);

void decompressDataSeries_float_3D_cuda(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data);

void decompressDataSeries_float_2D_cuda(float** data, size_t r1, size_t r2, unsigned char* comp_data);

void Huffman_init_cuda(int *s, size_t length, int thread_num);
