/**
 *  @file TightDataPointStorageD.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for the tight data point storage (TDPS).
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _TightDataPointStorageD_H
#define _TightDataPointStorageD_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TightDataPointStorageD
{
	int dataSeriesLength;
	int allSameData;
	int bestChoice;
	double realPrecision;
	
	int exactDataNum;
	double reservedValue;
	
	char* rtypeArray;
	int rtypeArray_size;
	
	char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
	int typeArray_size;
	
	char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
	int leadNumArray_size;
	
	char* exactMidBytes;
	int exactMidBytes_size;
	
	char* escBytes;
	int escBytes_size;
	
	char* residualMidBits;
	int residualMidBits_size;

} TightDataPointStorageD;

void new_TightDataPointStorageD_Empty(TightDataPointStorageD **this);
void new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **this, char* flatBytes, int flatBytesLength);
void decompressDataSeries_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps);
void decompressDataSeries_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps);
void decompressDataSeries_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps);
void getSnapshotData_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps);
void getSnapshotData_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps);
void getSnapshotData_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps);
void new_TightDataPointStorageD(TightDataPointStorageD **this, 
		int dataSeriesLength, int exactDataNum, 
		char* type, char* exactMidBytes, int exactMidBytes_size,
		char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		char* resiMidBits, int resiMidBits_size,
		char* escBytes, int escBytes_size,
		char* resiBitLength, int resiBitLengthSize, double realPrecision);
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, char** bytes, int *size);
void free_TightDataPointStorageD(TightDataPointStorageD *tdps);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _TightDataPointStorageD_H  ----- */
