/**
 *  @file ExpSegmentConstructor.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Constructor of Exponential Segment
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "sz.h"
#include "ExpSegment.h"
#include "ExpSegmentConstructor.h"
#include "DynamicByteArray.h"

void new_ExpSegmentConstructor_float(ExpSegmentConstructor **this, float* data, short reqExpo, float realPrecision)
{
	*this = (ExpSegmentConstructor*) malloc(sizeof(ExpSegmentConstructor));
	(*this)->precisionBitNum = 32;
	(*this)->floatData = data;
	(*this)->reqExpo = reqExpo;
	(*this)->floatRealPrecision = realPrecision;
} 

void new_ExpSegmentConstructor_double(ExpSegmentConstructor **this, double* data, short reqExpo, double realPrecision)
{
	*this = (ExpSegmentConstructor*) malloc(sizeof(ExpSegmentConstructor));
	(*this)->precisionBitNum = 64;
	(*this)->doubleData = data;
	(*this)->reqExpo = reqExpo;
	(*this)->doubleRealPrecision = realPrecision;
}

void new_ExpSegmentConstructor_escbytes(ExpSegmentConstructor **this, int precisionBitNum, char* flatBytes, 
			int flatBytesLength, int totalNumOfData)
{
	*this = (ExpSegmentConstructor*)malloc(sizeof(ExpSegmentConstructor));
	(*this)->floatData = NULL;
	(*this)->doubleData = NULL;
	(*this)->precisionBitNum = (unsigned char)precisionBitNum;
	if(precisionBitNum == 32)
	{
		new_ExpSegment(&((*this)->header), (short)-128, -1);
		ExpSegment* prevES = (*this)->header;
		
		(*this)->header->fixed = 1;
		int unitSize = 4+4+1;
		
		if(flatBytesLength%unitSize!=0)
		{
			printf("flatBytes.length(%d) mod unitSize(%d)!=0!!\n", flatBytesLength, unitSize);
			exit(0);
		}
		
		int startStep, num = flatBytesLength/unitSize;
		char startStepBytes[4];
		char reqLength, medValueBytes[4];
		float medValue;
		int i, j, k = 0;
		for(i = 0;i<num;i++)
		{
			memcpy(startStepBytes, &(flatBytes[k]),4);
			k+=4;
			//for(j = 0;j<4;j++)
			//	startStepBytes[j] = flatBytes[k++];
			memcpy(medValueBytes, &(flatBytes[k]),4);
			k+=4;		
			//for(j = 0;j<4;j++)
			//	medValueBytes[j] = flatBytes[k++];
			reqLength = flatBytes[k++];
			
			startStep = bytesToInt_bigEndian(startStepBytes);
			medValue = bytesToFloat(medValueBytes);
			
			ExpSegment *es;
			new_ExpSegment_float(&es, startStep, medValue, reqLength);
			es->medianValue_f = medValue;
			prevES->next = es;
			es->prev = prevES;
			if(i!=0)
			{
				prevES->endStep = startStep-1;
			}
			prevES = es;
		}
		
		prevES->endStep = totalNumOfData-1;
	}
	else //64
	{
		new_ExpSegment(&((*this)->header), (short)-1024, -1);
		ExpSegment* prevES = (*this)->header;

		(*this)->header->fixed = 1;
		int unitSize = 4+8+1;
		
		if(flatBytesLength%unitSize!=0)
		{
			printf("flatBytes.length(%d) mod unitSize(%d)!=0!!\n", flatBytesLength, unitSize);
			exit(0);
		}
		
		int startStep, num = flatBytesLength/unitSize;
		char startStepBytes[4];
		char reqLength, medValueBytes[8];
		double medValue;
					
		int i, j, k = 0;
		for(i = 0;i<num;i++)
		{
			memcpy(startStepBytes, &(flatBytes[k]),4);
			k+=4;			
			//for(j = 0;j<4;j++)
			//	startStepBytes[j] = flatBytes[k++];
			memcpy(medValueBytes, &(flatBytes[k]),8);
			k+=8;				
			//for(j = 0;j<8;j++)
			//	medValueBytes[j] = flatBytes[k++];
			reqLength = flatBytes[k++];
			
			startStep = bytesToInt_bigEndian(startStepBytes);
			medValue = bytesToDouble(medValueBytes);
			
			ExpSegment *es;
			new_ExpSegment_double(&es, startStep, medValue, reqLength);
			es->medianValue_d = medValue;
			prevES->next = es;
			es->prev = prevES;
			if(i!=0)
			{
				prevES->endStep = startStep-1;
			}
			prevES = es;
		}
		
		prevES->endStep = totalNumOfData-1;
	}
	(*this)->curExp = (*this)->header->next;
}

void free_ExpSegmentConstructor(ExpSegmentConstructor *esc)
{
	//free expSegments....
	ExpSegment *p = esc->header;
	ExpSegment *q;
	while(p!=NULL)
	{
		q = p->next;
		free(p);
		p = q;
	}
	//note that esc->floatData or esc->doubleData should not be freed here because they point to the original data.
	//if(esc->precisionBitNum == 32 && esc->floatData!=NULL)//single-precision
	//	free(esc->floatData);
	//else if(esc->precisionBitNum == 64 && esc->doubleData!=NULL)
	//	free(esc->doubleData);

	free(esc);
}

ExpSegment* createNewExpSegment(ExpSegment* prevES, short expo, int i)
{
	ExpSegment* es;
	new_ExpSegment(&es, expo, i);
	prevES->next = es;
	es->prev = prevES;
	return es;
}

ExpSegment* backTrackCompute(ExpSegment* curES, short nextLevel, int lengthBound)
{
	short next_level = nextLevel;
	if(curES->fixed==1 || curES->startStep==-1) //isHeader
		return NULL;
	ExpSegment* preES = curES->prev;
	ExpSegment* es; 
	mergeSelf_bestfit(curES, &es, preES->level, next_level, lengthBound);
	if(es==NULL) // a bump here
	{
		curES->fixed = 1;
		next_level = curES->level;
	}
	else
	{
		next_level = es->level;
		preES = es->prev;
	}
	ExpSegment* mergedES = backTrackCompute(preES, next_level, lengthBound);
	if(mergedES == NULL)
	{
		if(es==NULL) //this means es is a bump
			return curES;
		else
			return es;
	}
	else
		return mergedES;
}

void computeMinMaxMedian(ExpSegment* header, void* data, int dataLength, int DATA_TYPE)
{
	int i;
	ExpSegment* curES = header;
	if(DATA_TYPE==SZ_DOUBLE)
	{
		int startStep, endStep;
		double value = 0, minValue, maxValue;
		double *doubleData = (double *)data;
		while(curES->next!=NULL)
		{
			curES = curES->next;
			startStep = curES->startStep;
			endStep = getEndStep(curES);
			minValue = doubleData[startStep];
			maxValue = doubleData[startStep];
			for(i = startStep;i<=endStep;i++)
			{
				value = doubleData[i];
				if(value<minValue)
					minValue = value;
				else if(value>maxValue)
					maxValue = value;
			}
			curES->min_value_d = minValue;
			curES->max_value_d = maxValue;
			curES->medianValue_d0 = minValue + (maxValue - minValue)/2;
		}		
	}
	else if(DATA_TYPE==SZ_FLOAT)
	{
		int startStep, endStep;
		float value = 0, minValue, maxValue;
		float *floatData = (float *)data;
		while(curES->next!=NULL)
		{
			curES = curES->next;
			startStep = curES->startStep;
			endStep = getEndStep(curES);
			minValue = floatData[startStep];
			maxValue = floatData[startStep];
			for(i = startStep;i<=endStep;i++)
			{
				value = floatData[i];
				if(value<minValue)
					minValue = value;
				else if(value>maxValue)
					maxValue = value;
			}
			curES->min_value_f = minValue;
			curES->max_value_f = maxValue;
			curES->medianValue_f0 = minValue + (maxValue - minValue)/2;
		}		
	}
	else
	{
		printf("Error: wrong data type\n");
		exit(0);
	}
}

void computeOffset_3orders(ExpSegment* header, void* data_, int dataLength, float errBound, int DATA_TYPE)
{
	int i;
	ExpSegment* curExp = header->next;
	if(DATA_TYPE==SZ_FLOAT)
	{
		float *floatData = (float *)data_;
		while(curExp!=NULL)
		{
			int start = curExp->startStep;
			int end = curExp->endStep;
			
			float preData = floatData[start];
			float prepreData = preData;
			float preprepreData = prepreData;

			int preUnpreData = floatToOSEndianInt(preData);
			int leadingNumSum_mod8 = 0, count = 0;
			for(i = start;i<=end;i++)
			{
				float value = floatData[i];
				if(i==1)
				{
					if(fabs(value - preData) > errBound) //actually, it is equal to (value-medianValue)-(preData-medianValue)
					{
						float medianValue = curExp->medianValue_f0; //here, it's also known as lambda0
						//float medianValue = 0;
						int newUnpreData = floatToOSEndianInt(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Int(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				else if(i==2)//i>=2
				{
					if(fabs(value - preData)>errBound && fabs(value - (2*preData - prepreData))>errBound)
					{
						float medianValue = curExp->medianValue_f0;
						//float medianValue = 0;
						int newUnpreData = floatToOSEndianInt(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Int(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				else //i>=3
				{
					if(fabs(value - preData)>errBound && fabs(value - (2*preData - prepreData))>errBound && fabs(value - (3*preData - 3*prepreData + preprepreData))>errBound)
					{
						float medianValue = curExp->medianValue_f0;
						//float medianValue = 0;
						int newUnpreData = floatToOSEndianInt(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Int(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				preprepreData = prepreData;
				prepreData = preData;
				preData = value;
			}
			int meanLeadNum_mod8 = (int)(((float)leadingNumSum_mod8)/count);
			//Change places / Changing places, compared to SZ-0.5.11
			//float meanLeadNum_mod = meanLeadNum - ((int)(meanLeadNum/8))*8; 
			if(meanLeadNum_mod8 > 2) //4->2 is better for single-precision data
			{
				int offset_ = (8 - meanLeadNum_mod8)+offset;
				curExp->offset = (char)offset_;
			}

			curExp = curExp->next;
		}
	}
	else if(DATA_TYPE==SZ_DOUBLE)
	{
		double *doubleData = (double *)data_;
		while(curExp!=NULL)
		{
			int start = curExp->startStep;
			int end = curExp->endStep;
			
			double preData = doubleData[start];
			double prepreData = preData;
			double preprepreData = prepreData;
			
			long preUnpreData = doubleToOSEndianLong(preData);
			int leadingNumSum_mod8 = 0, count = 0;
			
			for(i = start;i<=end;i++)
			{
				double value = doubleData[i];
				if(i==1)
				{
					if(fabs(value - preData) > errBound) //actually, it is equal to (value-medianValue)-(preData-medianValue)
					{
						double medianValue = 0;//curExp->medianValue_d0; //here, it's also known as lambda0
						long newUnpreData = doubleToOSEndianLong(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Long(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				else if(i==2)//i>=2
				{
					if(fabs(value - preData)>errBound && fabs(value - (2*preData - prepreData))>errBound)
					{
						double medianValue = 0;//curExp->medianValue_d0;
						long newUnpreData = doubleToOSEndianLong(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Long(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				else //i>=3
				{
					if(fabs(value - preData)>errBound && fabs(value - (2*preData - prepreData))>errBound && fabs(value - (3*preData - 3*prepreData + preprepreData))>errBound)
					{
						double medianValue = 0;//curExp->medianValue_d0;
						long newUnpreData = doubleToOSEndianLong(value - medianValue);
						leadingNumSum_mod8 += getLeadingNumbers_Long(preUnpreData, newUnpreData)%8;
						preUnpreData = newUnpreData;
						count++;
					}
				}
				preprepreData = prepreData;
				prepreData = preData;
				preData = value;
			}
			int meanLeadNum_mod8 = (int)(((float)leadingNumSum_mod8)/count);
			//Change places / Changing places, compared to SZ-0.5.11
			//float meanLeadNum_mod = meanLeadNum - ((int)(meanLeadNum/8))*8; 
			if(meanLeadNum_mod8 > 4)
			{
				int offset_ = (8 - meanLeadNum_mod8)+offset;
				curExp->offset = (char)offset_;
			}

			curExp = curExp->next;
		}
	}
	else
	{
		printf("Error: wrong data type.\n");
		exit(0);
	}
}

void computeOffsetMedianValue(ExpSegment* header, void* data_, int DATA_TYPE)
{
	int i, j = 0;
	ExpSegment* curExp = header->next;
	if(DATA_TYPE == SZ_DOUBLE)
	{
		double *data = (double *) data_;
		while(curExp!=NULL)
		{
			char offset = curExp->offset;
			if(offset==0)
			{
				//curExp->medianValue_d = curExp->medianValue_d0;
				curExp->medianValue_d = 0;
			}
			else
			{
				//double lambda0 = curExp->medianValue_d0;
				double lambda0 = 0;
				int start  = curExp->startStep;
				int end = curExp->endStep;
				
				double sum = 0;
				for(i = start; i<=end;i++)
					sum += data[i];
				double avg = sum/(end-start+1);

				int exp = getExponent_double(avg - lambda0) + offset;
				if(exp >= 0)
				{
					double lambda = lambda0 - (1 << exp);
					curExp->medianValue_d = lambda;
				}
				else
				{
					double lambda = lambda0 - 1.0/(1 << -exp);
					curExp->medianValue_d = lambda;
				}
			}
			curExp = curExp->next;
		}
	}
	else if(DATA_TYPE == SZ_FLOAT)
	{
		float *data = (float *) data_;
		while(curExp!=NULL)
		{
			char offset = curExp->offset;
			if(offset==0)
			{
				curExp->medianValue_f = curExp->medianValue_f0;
				//curExp->medianValue_f = 0;
			}
			else
			{
				float lambda0 = curExp->medianValue_f0;
				//float lambda0 = 0;
				int start  = curExp->startStep;
				int end = curExp->endStep;
				
				float sum = 0;
				for(i = start; i<=end;i++)
					sum += data[i];
				float avg = sum/(end-start+1);

				int exp = getExponent_float(avg - lambda0) + offset;
				if(exp >= 0)
				{
					float lambda = lambda0 - (1 << exp);
					curExp->medianValue_f = lambda;
				}
				else
				{
					float lambda = lambda0 - 1.0/(1 << -exp);
					curExp->medianValue_f = lambda;
				}
			}
			curExp = curExp->next;
		}		
	}
	else
	{
		printf("Error: Wrong data type.\n");
		exit(0);
	}
}


void computeReqLength_double(ExpSegment* es, short reqExpo)
{
	double radius;
	double medianValue_d = es->medianValue_d;
	double medianValue0 = es->medianValue_d0;
	double minValue = es->min_value_d;
	double maxValue = es->max_value_d;	
	if(medianValue_d > medianValue0)
		radius = medianValue_d - minValue;
	else //<=
		radius = maxValue - medianValue_d;
	
	//retrieve the exponent of required precision.
	short radExpo = getExponent_double(radius);
	
	int reqMantLength = radExpo - reqExpo;
	int reqLength = 12+reqMantLength;
	if(reqLength>64)
	{
		medianValue_d = medianValue0;
		radius = medianValue_d -  minValue;
		radExpo = getExponent_double(radius);
		reqMantLength = radExpo - reqExpo;
		reqLength = 12+reqMantLength;
	}

	if(reqLength<12)
		reqLength = 12;
	if(reqLength>64)
	{
		reqLength = 64;
		medianValue_d = 0;
	}
	es->medianValue_d = medianValue_d;
	es->reqLength = reqLength;
	es->reqBytesLength = reqLength/8;
	es->resiBitsLength = reqLength%8;
}



void computeReqLength_float(ExpSegment* es, short reqExpo)
{
	float radius;
	float medianValue_f = es->medianValue_f;
	float medianValue0 = es->medianValue_f0;	
	float minValue = es->min_value_f;
	float maxValue = es->max_value_f;		
	if(medianValue_f > medianValue0)
		radius = medianValue_f - minValue;
	else //<=
		radius = maxValue - medianValue_f;
	
	//retrieve the exponent of required precision.
	short radExpo = getExponent_float(radius);
	
	int reqMantLength = radExpo - reqExpo;
	int reqLength = 9+reqMantLength;
	if(reqLength>32)
	{
		medianValue_f = medianValue0;
		radius = medianValue0 - minValue;
		//TODO: SHENGDI
		radExpo = getExponent_float(radius);
		reqMantLength = radExpo - reqExpo;
		reqLength = 9+reqMantLength;
	}
	if(reqLength<9)
		reqLength = 9;
	if(reqLength>32)
	{
		reqLength = 32;
		medianValue_f = 0;
	}
	
	es->medianValue_f = medianValue_f;
	es->reqLength = reqLength;
	es->reqBytesLength = reqLength/8;
	es->resiBitsLength = reqLength%8;
}

void computeReqLength(ExpSegmentConstructor* this, short reqExpo)
{
	if(this->precisionBitNum == 32)
	{
		ExpSegment* curES = this->header;
		while(curES->next!=NULL)
		{
			curES = curES->next;			
			computeReqLength_float(curES, reqExpo);
		}
	}
	else //64 double data
	{
		ExpSegment* curES = this->header;
		while(curES->next!=NULL)
		{
			curES = curES->next;
			computeReqLength_double(curES, reqExpo);
		}
	}
}

inline void getExpSegment_fast(ExpSegmentConstructor *esc, int index)
{
	while(index > esc->curExp->endStep)
		esc->curExp = (esc->curExp)->next;
}

void cleanESwithNoUnpredictable(ExpSegmentConstructor *esc)
{
	ExpSegment *nextES = esc->header->next;
	while(nextES!=NULL)
	{
		if(nextES->unpredNum==0)
		{
			mergeSelf_toprev(nextES);
			nextES = nextES->prev->next;
		}
		else
			nextES = nextES->next;
	}
	esc->curExp = esc->header->next;
}

int convertESCToBytes(ExpSegmentConstructor *esc, char **bytes)
{
	int i, size_;
	if(esc->precisionBitNum == 32)
	{
		char reqLength;
		char startBytes[4], mediaBytes[4];
		DynamicByteArray* dba;
		new_DBA(&dba, 1024);

		ExpSegment *curES = esc->header;
		while(curES->next!=NULL)
		{
			curES = curES->next;
			reqLength = (char)curES->reqLength;
			//TODO
			intToBytes_bigEndian(startBytes, curES->startStep); //int
			floatToBytes(mediaBytes, curES->medianValue_f); //int
			
			for(i = 0;i<4;i++)
				addDBA_Data(dba, startBytes[i]);
			for(i = 0;i<4;i++)
				addDBA_Data(dba, mediaBytes[i]);

			addDBA_Data(dba, reqLength);
		}
		size_ = dba->size;
		convertDBAtoBytes(dba, bytes);
		free_DBA(dba);
	}
	else //64 double data
	{
		char reqLength;
		char startBytes[8], mediaBytes[8];	
		DynamicByteArray *dba;
		new_DBA(&dba, 1024);
		//int k = 0,j=0;
		ExpSegment *curES = esc->header;
		while(curES->next!=NULL)
		{
			curES = curES->next;
			reqLength = (char)curES->reqLength;
			//printf("%d (%d,%d) level=%d median=%f median0=%f\n", k++, startStep, curES->startStep+curES->length-1, curES->level, curES->medianValue_d, curES->medianValue_d0);

			intToBytes_bigEndian(startBytes, curES->startStep); //int
			doubleToBytes(mediaBytes, curES->medianValue_d);			
			
			for(i = 0;i<4;i++)
				addDBA_Data(dba, startBytes[i]);
			for(i = 0;i<8;i++)
				addDBA_Data(dba, mediaBytes[i]);

			addDBA_Data(dba, reqLength);
		}
		size_ = dba->size;
		convertDBAtoBytes(dba, bytes);
		free_DBA(dba);
	}
	return size_;
}

void construct(ExpSegmentConstructor *esc, void* data_, int dataLength, int lengthBound, short reqExpo, float errBound, int SZ_DATA_Type)
{
	int debug_counter = 0;
	if(SZ_DATA_Type == SZ_DOUBLE)
	{
		double* data = (double *)data_;
		ExpSegment* esHeader;
		new_ExpSegment(&esHeader, (short)-1024,-1); //-1024 must be smaller than any possible exponent (such as -1023)
		esHeader->fixed = 1;
		short preExpo = getExponent_double(data[0]);
		if(preExpo < reqExpo)
			preExpo = reqExpo; //because mantiReqLength = expo - reqExpo, which should be >=0
		ExpSegment* initES;
		new_ExpSegment(&initES, preExpo, 0);
		ExpSegment* prevES = initES;
		esHeader->next = initES;
		initES->prev = esHeader;
		int i;
		double value;
		short expo;
		for(i = 1;i<dataLength;i++)
		{
			value = data[i];
			expo = getExponent_double(value);
			if(expo < reqExpo)
				expo = reqExpo; //because mantiReqLength = expo - reqExpo, which should be >=0
				
			if(expo < preExpo) //record
			{
				prevES = createNewExpSegment(prevES, expo, i);
			}
			else if(expo > preExpo) //compute
			{
				ExpSegment* mergedES = backTrackCompute(prevES, expo, lengthBound);
				if(mergedES==NULL) //seems only related to header
					mergedES = prevES;
				else
				{
					toLatestES(&mergedES); //convert it to the latest/last ES
				}
				if(prevES->fixed==1 && prevES->level < expo)
				{
					prevES = createNewExpSegment(prevES, expo, i);
				}
				else
				{
					prevES = mergedES;
					prevES->length += 1;
				}
			}
			else //==
				prevES->length += 1;
			
			preExpo = expo;
		}
		
		backTrackCompute(prevES, prevES->level, lengthBound);
		
		cleanUp(esHeader, data, dataLength, errBound, SZ_DOUBLE);
		
		updateEndStep(esHeader);
		
		computeMinMaxMedian(esHeader, data, dataLength, SZ_DOUBLE);
		
		computeOffset_3orders(esHeader, data, dataLength, errBound, SZ_DOUBLE);
		
		computeOffsetMedianValue(esHeader, data, SZ_DOUBLE);
		
		esc->header = esHeader;		
	}
	else if(SZ_DATA_Type == SZ_FLOAT)
	{
		float* data = (float *)data_;
		ExpSegment* esHeader; 
		new_ExpSegment(&esHeader, (short)-128,-1); //-1024 must be smaller than any possible exponent (such as -1023)
		esHeader->fixed = 1;
		
		short preExpo = getExponent_float(data[0]);
		if(preExpo < reqExpo)
			preExpo = reqExpo; //because mantiReqLength = expo - reqExpo, which should be >=0
		ExpSegment* initES;
		new_ExpSegment(&initES, preExpo, 0);
		ExpSegment* prevES = initES;
		esHeader->next = initES;
		initES->prev = esHeader;
		int i;
		float value;
		short expo;
		for(i = 1;i<dataLength;i++)
		{
			value = data[i];
			expo = getExponent_float(value);
			if(expo < reqExpo)
				expo = reqExpo; //because mantiReqLength = expo - reqExpo, which should be >=0
			if(expo < preExpo) //record
			{
				prevES = createNewExpSegment(prevES, expo, i);
			}
			else if(expo > preExpo) //compute
			{
				ExpSegment* mergedES = backTrackCompute(prevES, expo, lengthBound);
				if(mergedES==NULL) //seems only related to header
					mergedES = prevES;
				else
				{
					toLatestES(&mergedES); //convert it to the latest/last ES
				}
				if(prevES->fixed == 1 && prevES->level < expo)
				{
					prevES = createNewExpSegment(prevES, expo, i);
				}
				else
				{
					prevES = mergedES;
					prevES->length += 1;
				}
			}
			else //==
				prevES->length += 1;
			
			preExpo = expo;
		}
		
		backTrackCompute(prevES, prevES->level, lengthBound);
		
		cleanUp(esHeader, data, dataLength, errBound, SZ_FLOAT);
		
		updateEndStep(esHeader);
		
		computeMinMaxMedian(esHeader, data, dataLength, SZ_FLOAT);
		
		computeOffset_3orders(esHeader, data, dataLength, errBound, SZ_FLOAT);
	
		computeOffsetMedianValue(esHeader, data, SZ_FLOAT);
		
		esc->header = esHeader;
	}
	else
	{
		printf("error type\n");
	}
}

