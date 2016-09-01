/**
 *  @file ExpSegmentConstructor.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for the exponential segment constructor.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ExpSegmentConstructor_H
#define _ExpSegmentConstructor_H

#include "ExpSegment.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ExpSegmentConstructor
{	
	ExpSegment *header;
	short reqExpo;
	char precisionBitNum; //xx=32 for single-precision and xx=64 for double-precision
	float* floatData;
	double* doubleData;
	float floatRealPrecision;
	double doubleRealPrecision;
	ExpSegment *curExp;

} ExpSegmentConstructor;

void new_ExpSegmentConstructor_float(ExpSegmentConstructor **this, float* data, short reqExpo, float realPrecision);
void new_ExpSegmentConstructor_double(ExpSegmentConstructor **this, double* data, short reqExpo, double realPrecision);
void new_ExpSegmentConstructor_escbytes(ExpSegmentConstructor **this, int precisionBitNum, char* flatBytes, 
			int flatBytesLength, int totalNumOfData);
void free_ExpSegmentConstructor(ExpSegmentConstructor *esc);
ExpSegment* createNewExpSegment(ExpSegment* prevES, short expo, int i);
ExpSegment* backTrackCompute(ExpSegment* curES, short nextLevel, int lengthBound);
void computeMinMaxMedian(ExpSegment* header, void* data, int dataLength, int DATA_TYPE);
void computeOffset_3orders(ExpSegment* header, void* data_, int dataLength, float errBound, int DATA_TYPE);
void computeOffsetMedianValue(ExpSegment* header, void* data_, int DATA_TYPE);
void computeReqLength_double(ExpSegment* es, short reqExpo);
void computeReqLength_float(ExpSegment* es, short reqExpo);
void computeReqLength(ExpSegmentConstructor* this, short reqExpo);
void getExpSegment_fast(ExpSegmentConstructor *esc, int index);
void cleanESwithNoUnpredictable(ExpSegmentConstructor *esc);
int convertESCToBytes(ExpSegmentConstructor *esc, char **bytes);
void construct(ExpSegmentConstructor *esc, void* data_, int dataLength, int lengthBound, short reqExpo, float errBound, int SZ_DATA_Type);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ExpSegmentConstructor_H  ----- */
