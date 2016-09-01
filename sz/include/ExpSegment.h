/**
 *  @file ExpSegment.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for the exponential segment constructor.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ExpSegment_H
#define _ExpSegment_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ExpSegment
{	
	int fixed;
	short level;
	int startStep;
	int length;
	int endStep;
	int unpredNum;
	
	double min_value_d; 
	double max_value_d;
	float min_value_f; 
	float max_value_f;
	
	double medianValue_d0;
	float medianValue_f0;
	
	double medianValue_d;
	float medianValue_f;
	
	int reqLength;
	int reqBytesLength;
	int resiBitsLength;
	char offset; //offset is equal to [8 - beta]+2 (or 1), where beta  alpha - [alpha/8]*8, where alpha average length of leadingzero number
	
	struct ExpSegment* prev; 
	struct ExpSegment* next;
	
} ExpSegment;

void new_ExpSegment(ExpSegment **this, short level, int startStep);
void new_ExpSegment_float(ExpSegment **this, int startStep, float medianValue0, int reqLength);
void new_ExpSegment_double(ExpSegment **this, short startStep, double medianValue0, int reqLength);
void toLatestES(ExpSegment **es);
void removeSelf(ExpSegment* this);
void mergeSelf_bestfit(ExpSegment* this, ExpSegment** merged, short prevLevel, short nextLevel, int lengthBound);
short mergeSelf_toprev(ExpSegment* this);
short mergeSelf_tonext(ExpSegment* this);
int getEndStep(ExpSegment* this);
void cleanUp(ExpSegment* this, void* data, long dataLength, float errBound, int DATA_TYPE);
void updateEndStep(ExpSegment* this);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ExpSegment_H  ----- */
