/**
 *  @file TightDataPointStorageF.h
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief Header file for the tight data point storage (TDPS).
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _TightDataPointStorageF_H
#define _TightDataPointStorageF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h> 

typedef struct TightDataPointStorageF
{
	size_t dataSeriesLength;
	int allSameData;
	double realPrecision; //it's used as the pwrErrBoundRatio when errBoundMode==PW_REL
	float medianValue;
	char reqLength;
	char radExpo; //used to compute reqLength based on segmented precisions in "pw_rel_compression"
	
	size_t exactDataNum;
	float reservedValue;
	
	unsigned char* rtypeArray;
	size_t rtypeArray_size;
	
	unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
	size_t typeArray_size;
	
	unsigned char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
	size_t leadNumArray_size;
	
	unsigned char* exactMidBytes;
	size_t exactMidBytes_size;
	
	unsigned char* residualMidBits;
	size_t residualMidBits_size;
	
	unsigned int intervals; //quantization_intervals
	
	unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)
	
	size_t segment_size;
	
	unsigned char* pwrErrBoundBytes;
	int pwrErrBoundBytes_size;
} TightDataPointStorageF;

void new_TightDataPointStorageF_Empty(TightDataPointStorageF **self);
int new_TightDataPointStorageF_fromFlatBytes(TightDataPointStorageF **self, unsigned char* flatBytes, size_t flatBytesLength);
void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps);
void decompressDataSeries_float_1D_pwr(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps);

float* extractRealPrecision_2D_float(size_t R1, size_t R2, int blockSize, TightDataPointStorageF* tdps);
void decompressDataSeries_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps);
void decompressDataSeries_float_2D_pwr(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps);

float* extractRealPrecision_3D_float(size_t R1, size_t R2, size_t R3, int blockSize, TightDataPointStorageF* tdps);
void decompressDataSeries_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps);
void decompressDataSeries_float_3D_pwr(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps);

void decompressDataSeries_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps);

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps, int errBoundMode);

void new_TightDataPointStorageF(TightDataPointStorageF **self,
		size_t dataSeriesLength, size_t exactDataNum,
		int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
		unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		unsigned char* resiMidBits, size_t resiMidBits_size,
		unsigned char* resiBitLength, size_t resiBitLengthSize, 
		double realPrecision, float medianValue, char reqLength, unsigned int intervals, 
		unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo);

void convertTDPStoBytes_float(TightDataPointStorageF* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte);
void convertTDPStoBytes_float_reserve(TightDataPointStorageF* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte);
void convertTDPStoFlatBytes_float(TightDataPointStorageF *tdps, unsigned char** bytes, size_t *size);
void convertTDPStoFlatBytes_float_args(TightDataPointStorageF *tdps, unsigned char* bytes, size_t *size);

void free_TightDataPointStorageF(TightDataPointStorageF *tdps);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _TightDataPointStorageF_H  ----- */
