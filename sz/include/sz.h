/**
 *  @file sz.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole detector.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_H
#define _SZ_H

#include <stdio.h>
#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */
#include "iniparser.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "VarSet.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#define HZ 102
#define SZ 101

#define ABS 0
#define REL 1
#define ABS_AND_REL 2
#define ABS_OR_REL 3

#define SZ_FLOAT 0
#define SZ_DOUBLE 1

#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1

#define DynArrayInitLen 1024

int sysEndianType; //endian type of the system
int dataEndianType; //endian type of the data
int maxSegmentNum;

char maxHeap[10];
 
long status;

int sol_ID;
int errorBoundMode; //ABS, REL, ABS_AND_REL, or ABS_OR_REL

int gzipMode; //four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION

char *sz_cfgFile;

int offset;

double absErrBound;
double relBoundRatio;

int versionNumber[3];

SZ_VarSet* sz_varset;

typedef union lshort
{
	unsigned short svalue;
	char byte[2];
} lshort;

typedef union ldouble
{
    double value;
    unsigned long lvalue;
    char byte[8];
} ldouble;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    char byte[4];
} lfloat;

//dataCompression.c
void computeRangeSize_double(double* oriData, int size, double* valueRangeSize, double* medianValue);
void computeRangeSize_float(float* oriData, int size, float* valueRangeSize, float* medianValue);
double min_d(double a, double b);
double max_d(double a, double b);
float min_f(float a, float b);
float max_f(float a, float b);
double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio);
float getRealPrecision_float(float valueRangeSize, int errBoundMode, float absErrBound, float relBoundRatio);
void symTransform_8bytes(char data[8]);
void flush_to_bigEndian_8bytes(char data[8], int dataEndianType);
void symTransform_4bytes(char data[4]);
void flush_to_bigEndian_4bytes(char data[4]);
void bigEndian_to_OSEndian_double(char data[8]);
void bigEndian_to_OSEndian_float(char data[4]);
void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
int compIdenticalLeadingBytesCount_double(char* preBytes, char* curBytes);
int compIdenticalLeadingBytesCount_float(char* preBytes, char* curBytes);
void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce);

//ByteToolkit.c
inline int bytesToInt_bigEndian(char* bytes);
inline void intToBytes_bigEndian(char *b, unsigned int num);
inline long bytesToLong_bigEndian(char* b);
inline void longToBytes_bigEndian(char *b, unsigned long num);
long doubleToOSEndianLong(double value);
int floatToOSEndianInt(float value);
inline short getExponent_float(float value);
short getPrecisionReqLength_float(float precision);
inline short getExponent_double(double value);
short getPrecisionReqLength_double(double precision);
inline int numberOfLeadingZeros_Int(int i);
inline int numberOfLeadingZeros_Long(long i);
inline char getLeadingNumbers_Int(int v1, int v2);
inline char getLeadingNumbers_Long(long v1, long v2);
short bytesToShort(char* bytes);
int bytesToInt(char* bytes);
long bytesToLong(char* bytes);
float bytesToFloat(char* bytes);
void floatToBytes(char *b, float num);
double bytesToDouble(char* bytes);
void doubleToBytes(char *b, double num);
int extractBytes(char* byteArray, int k, int validLength);
inline int getMaskRightCode(int m);
inline int getLeftMovingCode(int kMod8);
inline int getRightMovingSteps(int kMod8, int resiBitLength);
inline int getRightMovingCode(int kMod8, int resiBitLength);

//TypeManager.c
int convertIntArray2ByteArray_fast_2b(char* timeStepType, int timeStepTypeLength, char **result);
void convertByteArray2IntArray_fast_2b(int stepLength, char* byteArray, int byteArrayLength, int **intArray);
int getLeftMovingSteps(int k, char resiBitLength);
int convertIntArray2ByteArray_fast_dynamic(char* timeStepType, char* resiBitLength, int resiBitLengthLength, char **bytes);
int computeBitNumRequired(int dataLength);
void decompressBitArraybySimpleLZ77(int** result, char* bytes, int bytesLength, int totalLength, int validLength);

//test_zlib.c
ulong zlib_compress(char* data, ulong dataLength, char** compressBytes, int level);
ulong zlib_compress2(char* data, ulong dataLength, char** compressBytes, int level);
ulong zlib_uncompress(char* compressBytes, ulong cmpSize, char** oriData, ulong targetOriSize);
ulong zlib_uncompress2(char* compressBytes, ulong cmpSize, char** oriData, ulong targetOriSize);
ulong zlib_uncompress3(char* compressBytes, ulong cmpSize, char** oriData, ulong targetOriSize);

//szf.c
void sz_init_c_(char *configFile,int *len,int *ierr);
void sz_finalize_c_();
void SZ_writeData_inBinary_d1_Float_(float* data, char *fileName, int *len);
void sz_compress_d1_float_(float* data, char *bytes, int *outSize, int *r1);
void sz_compress_d1_float_rev_(float* data, float *reservedValue, char *bytes, int *outSize, int *r1);
void sz_compress_d2_float_(float* data, char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_float_rev_(float* data, float *reservedValue, char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_float_(float* data, char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_float_rev_(float* data, float *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_float_(float* data, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_float_rev_(float* data, float *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_(float* data, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_float_rev_(float* data, float *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_double_(double* data, char *bytes, int *outSize, int *r1);
void sz_compress_d1_double_rev_(double* data, double *reservedValue, char *bytes, int *outSize, int *r1);
void sz_compress_d2_double_(double* data, char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_double_rev_(double* data, double *reservedValue, char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_double_(double* data, char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_double_rev_(double* data, double *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_double_(double* data, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_double_rev_(double* data, double *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_(double* data, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_double_rev_(double* data, double *reservedValue, char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_args_(float* data, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_args_(float* data, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_args_(float* data, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_args_(float* data, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_args_(float* data, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_args_(double* data, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_args_(double* data, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_args_(double* data, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_args_(double* data, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_args_(double* data, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_rev_args_(float* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_rev_args_(float* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_rev_args_(float* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_rev_args_(float* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_rev_args_(float* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_rev_args_(double* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_rev_args_(double* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_rev_args_(double* data, float *reservedValue, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_rev_args_(double* data, double *reservedValue, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_rev_args_(double* data, double *reservedValue, char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_decompress_d1_float_(char *bytes, int *byteLength, float *data, int *r1);
void sz_decompress_d2_float_(char *bytes, int *byteLength, float *data, int *r1, int *r2);
void sz_decompress_d3_float_(char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_float_(char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_float_(char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_decompress_d1_double_(char *bytes, int *byteLength, double *data, int *r1);
void sz_decompress_d2_double_(char *bytes, int *byteLength, double *data, int *r1, int *r2);
void sz_decompress_d3_double_(char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_double_(char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_double_(char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_batchaddVar_d1_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_batchaddvar_d2_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchaddvar_d1_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_batchaddvar_d2_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchdelvar_c_(char* varName, int *len, int *errState);
void sz_batch_compress_c_(char* bytes, int *outSize);
void sz_batch_decompress_c_(char* bytes, int *byteLength);
void sz_getvardata_float_(char* varName, int *len, float* data, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_getvardata_double_(char* varName, int *len, double* data, int *r1, int *r2, int *r3, int *r4, int *r5);


//sz.h
int SZ_Init(char *configFilePath);
int computeDataLength(int r5, int r4, int r3, int r2, int r1);
int computeDimension(int r5, int r4, int r3, int r2, int r1);

void SZ_compress_args_float_NoCkRngeNoGzip(char** newByteData, float *oriData, int dataLength, float realPrecision, int *outSize);
void SZ_compress_args_double_NoCkRngeNoGzip(char** newByteData, double *oriData, int dataLength, double realPrecision, int *outSize);
void SZ_compress_args_float_withinRange(char** newByteData, float *oriData, int dataLength, int *outSize);
void SZ_compress_args_double_withinRange(char** newByteData, double *oriData, int dataLength, int *outSize);

void SZ_compress_args_float(char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, float absErr_Bound, float rel_BoundRatio);
void SZ_compress_args_double(char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

void SZ_compress_args_float_wRngeNoGzip(char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, float absErr_Bound, float rel_BoundRatio);
void SZ_compress_args_double_wRngeNoGzip(char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

char *SZ_compress(int dataType, void *data, int *outSize, int r5, int r4, int r3, int r2, int r1);
char *SZ_compress_args(int dataType, void *data, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_args2(int dataType, void *data, char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);

char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
char *SZ_compress_rev(int dataType, void *data, void *reservedValue, int *outSize, int r5, int r4, int r3, int r2, int r1);

void SZ_decompress_args_float(float** newData, int r5, int r4, int r3, int r2, int r1, char* cmpBytes, int cmpSize);
void SZ_decompress_args_double(double** newData, int r5, int r4, int r3, int r2, int r1, char* cmpBytes, int cmpSize);
void *SZ_decompress(int dataType, char *bytes, int byteLength, int r5, int r4, int r3, int r2, int r1);
int SZ_decompress_args(int dataType, char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1);

void filloutDimArray(int* dim, int r5, int r4, int r3, int r2, int r1);
char* SZ_batch_compress(int *outSize);
SZ_VarSet* SZ_batch_decompress(char* compressedStream, int length);
void SZ_Finalize();

#endif /* ----- #ifndef _SZ_H  ----- */
