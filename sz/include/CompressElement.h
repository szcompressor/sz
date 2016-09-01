/**
 *  @file CompressElement.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Compress Elements such as DoubleCompressELement.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _CompressElement_H
#define _CompressElement_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DoubleValueCompressElement
{
	double data;
	long curValue;
	char curBytes[8]; //big_endian
	int reqBytesLength;
	int resiBitsLength;
} DoubleValueCompressElement;

typedef struct FloatValueCompressElement
{
	float data;
	int curValue;
	char curBytes[4]; //big_endian
	int reqBytesLength;
	int resiBitsLength;
} FloatValueCompressElement;

typedef struct LossyCompressionElement
{
	int leadingZeroBytes; //0,1,2,or 3
	char integerMidBytes[8];
	int integerMidBytes_Length; //they are mid_bits actually
	//char curBytes[8];
	//int curBytes_Length; //4 for single_precision or 8 for double_precision	
	int resMidBitsLength;
	int residualMidBits;
} LossyCompressionElement;

inline void listAdd_double(double last3CmprsData[3], double value);
inline void listAdd_float(float last3CmprsData[3], float value);
inline int validPrediction_double(double minErr, double precision);
inline int validPrediction_float(float minErr, float precision);
void new_LossyCompressionElement(LossyCompressionElement *lce, int leadingNum, char* intMidBytes, 
		int intMidBytes_Length, int resiMidBitsLength, int resiBits);
inline void updateLossyCompElement_Double(char* curBytes, char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce);
inline void updateLossyCompElement_Float(char* curBytes, char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _CompressElement_H  ----- */
