/**
 *  @file sz_float.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
#include "ExpSegmentConstructor.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "TightDataPointStorageF.h"
#include "zlib.h"
#include "rw.h"

void SZ_compress_args_float_NoCkRngeNoGzip_1D(char** newByteData, float *oriData, int dataLength, float realPrecision, int *outSize)
{
	int i;
	short reqExpo = getPrecisionReqLength_float(realPrecision);
	
	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;
		
	float* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
		
	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_float(&esc, spaceFillingValue, reqExpo, realPrecision);
			
	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_FLOAT);
	esc->curExp = esc->header->next;
	
	computeReqLength(esc, reqExpo);
			
	//DynamicByteArray *reqByteLengthArray;
	//new_DBA(&reqByteLengthArray, DynArrayInitLen);
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);
	
	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	type[0] = 0;
	
	char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;
			
	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	float medianValue = curExp->medianValue_f;
	float last3CmprsData[3] = {0};

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
				
	//add the first data	
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
		
	//add the second data
	type[1] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);			
	compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
	
	unsigned char state;
	float lcf, qcf;		
	float interval, checkRadius, curData;
	float pred;
	float predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;		
	checkRadius = 255*realPrecision;
	interval = 2*realPrecision;
	
	for(i=2;i<dataLength;i++)
	{				
		curData = spaceFillingValue[i];
		pred = 2*last3CmprsData[0] - last3CmprsData[1];
		predAbsErr = fabs(curData - pred);	
		if(predAbsErr<=checkRadius)
		{
			state = (unsigned char)((predAbsErr/realPrecision+1)/2);
			if(curData>=pred)
			{
				type[i] = 128+state;
				pred = pred + state*interval;
			}
			else //curData<pred
			{
				type[i] = 128-state;
				pred = pred - state*interval;
			}
			listAdd_float(last3CmprsData, pred);					
			continue;
		}
		
		//unpredictable data processing
		getExpSegment_fast(esc, i);
		curExp = esc->curExp;			
		type[i] = 0;
		curExp->unpredNum++;
		medianValue = curExp->medianValue_f;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;
		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
						
		listAdd_float(last3CmprsData, vce->data);	
	}//end of for
			
	cleanESwithNoUnpredictable(esc);
		
	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			expSegmentsInBytes, expSegmentsInBytes_size, 
			resiBitLengthArray->array, resiBitLengthArray->size, realPrecision);

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	//free(type);
	
	free_ExpSegmentConstructor(esc);
		
	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
//	TightDataPointStorageF* tdps2;
//	new_TightDataPointStorageF_fromFlatBytes(&tdps2, *newByteData, outSize);

	free_DBA(exactMidByteArray);	
	
	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);	
}

/**
 * 
 * Note: @r1 is high dimension
 * 		 @r2 is low dimension 
 * */
void SZ_compress_args_float_NoCkRngeNoGzip_2D(char** newByteData, float *oriData, int r1, int r2, float realPrecision, int *outSize)
{
	int i,j;
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;
		
	int dataLength = r1*r2;	
	
	P0 = (float*)malloc(r2*sizeof(float));
	P1 = (float*)malloc(r2*sizeof(float));
		
	short reqExpo = getPrecisionReqLength_float(realPrecision);
	
	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;
		
	float* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
				
	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_float(&esc, spaceFillingValue, reqExpo, realPrecision);
			
	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_FLOAT);
	esc->curExp = esc->header->next;
	
	computeReqLength(esc, reqExpo);
			
	//DynamicByteArray *reqByteLengthArray;
	//new_DBA(&reqByteLengthArray, DynArrayInitLen);
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);
	
	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	type[0] = 0;
	
	char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;
			
	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	float medianValue = curExp->medianValue_f;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
			
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < 256)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + 128;
		P1[1] = pred1D + 2 * (type[1] - 128) * realPrecision;
	}
	else
	{
		type[1] = 0;

		getExpSegment_fast(esc, 1);
		curExp = esc->curExp;
		curExp->unpredNum++;
		medianValue = curExp->medianValue_f;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;

		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < 256)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + 128;
			P1[j] = pred1D + 2 * (type[j] - 128) * realPrecision;
		}
		else
		{
			type[j] = 0;

			getExpSegment_fast(esc, j);
			curExp = esc->curExp;
			curExp->unpredNum++;
			medianValue = curExp->medianValue_f;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r1-1 */
	int index;
	for (i = 1; i < r1; i++)
	{	
		/* Process row-i data 0 */
		index = i*r2;
		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < 256)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + 128;
			P0[0] = pred1D + 2 * (type[index] - 128) * realPrecision;
		}
		else
		{
			type[index] = 0;

			getExpSegment_fast(esc, index);
			curExp = esc->curExp;
			curExp->unpredNum++;
			medianValue = curExp->medianValue_f;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}
									
		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < 256)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + 128;
				P0[j] = pred2D + 2 * (type[index] - 128) * realPrecision;
			}
			else
			{
				type[index] = 0;

				getExpSegment_fast(esc, index);
				curExp = esc->curExp;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);			
	cleanESwithNoUnpredictable(esc);
		
	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			expSegmentsInBytes, expSegmentsInBytes_size, 
			resiBitLengthArray->array, resiBitLengthArray->size, realPrecision);

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	//free(type);
	
	free_ExpSegmentConstructor(esc);
		
	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
//	TightDataPointStorageF* tdps2;
//	new_TightDataPointStorageF_fromFlatBytes(&tdps2, *newByteData, outSize);

	free_DBA(exactMidByteArray);	
	
	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);	
}

void SZ_compress_args_float_NoCkRngeNoGzip_3D(char** newByteData, float *oriData, int r1, int r2, int r3, float realPrecision, int *outSize)
{
	int i,j,k;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int dataLength = r1*r2*r3;

	P0 = (float*)malloc(r2*r3*sizeof(float));
	P1 = (float*)malloc(r2*r3*sizeof(float));

	short reqExpo = getPrecisionReqLength_float(realPrecision);

	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;

	float* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;

	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_float(&esc, spaceFillingValue, reqExpo, realPrecision);

	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_FLOAT);
	esc->curExp = esc->header->next;

	computeReqLength(esc, reqExpo);

	//DynamicByteArray *reqByteLengthArray;
	//new_DBA(&reqByteLengthArray, DynArrayInitLen);

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;

	char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;

	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	float medianValue = curExp->medianValue_f;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < 256)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + 128;
		P1[1] = pred1D + 2 * (type[1] - 128) * realPrecision;
	}
	else
	{
		type[1] = 0;

		getExpSegment_fast(esc, 1);
		curExp = esc->curExp;
		curExp->unpredNum++;
		medianValue = curExp->medianValue_f;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;

		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < 256)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + 128;
			P1[j] = pred1D + 2 * (type[j] - 128) * realPrecision;
		}
		else
		{
			type[j] = 0;

			getExpSegment_fast(esc, j);
			curExp = esc->curExp;
			curExp->unpredNum++;
			medianValue = curExp->medianValue_f;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r2-1 */
	int index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;
		pred1D = P1[index-r3];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < 256)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + 128;
			P1[index] = pred1D + 2 * (type[index] - 128) * realPrecision;
		}
		else
		{
			type[index] = 0;

			getExpSegment_fast(esc, index);
			curExp = esc->curExp;
			curExp->unpredNum++;
			medianValue = curExp->medianValue_f;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < 256)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + 128;
				P1[index] = pred2D + 2 * (type[index] - 128) * realPrecision;
			}
			else
			{
				type[index] = 0;

				getExpSegment_fast(esc, index);
				curExp = esc->curExp;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r2*r3;
		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < 256)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + 128;
			P0[0] = pred1D + 2 * (type[index] - 128) * realPrecision;
		}
		else
		{
			type[index] = 0;

			getExpSegment_fast(esc, index);
			curExp = esc->curExp;
			curExp->unpredNum++;
			medianValue = curExp->medianValue_f;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}


	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			index = k*r2*r3+j;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < 256)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + 128;
				P0[j] = pred2D + 2 * (type[index] - 128) * realPrecision;
			}
			else
			{
				type[index] = 0;

				getExpSegment_fast(esc, index);
				curExp = esc->curExp;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

	    /* Process Row-1 --> Row-r2-1 */
		int index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r2*r3 + i*r3;
			index2D = i*r3;
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < 256)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + 128;
				P0[index2D] = pred2D + 2 * (type[index] - 128) * realPrecision;
			}
			else
			{
				type[index] = 0;

				getExpSegment_fast(esc, index);
				curExp = esc->curExp;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				index = k*r2*r3 + i*r3 + j;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				diff = spaceFillingValue[index] - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < 256)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + 128;
					P0[index2D] = pred3D + 2 * (type[index] - 128) * realPrecision;
				}
				else
				{
					type[index] = 0;

					getExpSegment_fast(esc, index);
					curExp = esc->curExp;
					curExp->unpredNum++;
					medianValue = curExp->medianValue_f;
					reqLength = curExp->reqLength;
					reqBytesLength = curExp->reqBytesLength;
					resiBitsLength = curExp->resiBitsLength;

					addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);
	cleanESwithNoUnpredictable(esc);

	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			expSegmentsInBytes, expSegmentsInBytes_size,
			resiBitLengthArray->array, resiBitLengthArray->size, realPrecision);

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n",
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);

//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	//free(type);

	free_ExpSegmentConstructor(esc);

	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

//	TightDataPointStorageF* tdps2;
//	new_TightDataPointStorageF_fromFlatBytes(&tdps2, *newByteData, outSize);

	free_DBA(exactMidByteArray);

	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_withinRange(char** newByteData, float *oriData, int dataLength, int *outSize)
{
	TightDataPointStorageF* tdps = (TightDataPointStorageF*) malloc(sizeof(TightDataPointStorageF));
	tdps->rtypeArray = NULL;
	tdps->typeArray = NULL;	
	tdps->leadNumArray = NULL;
	tdps->escBytes = NULL;
	tdps->residualMidBits = NULL;
	
	tdps->allSameData = 1;
	tdps->dataSeriesLength = dataLength;
	tdps->exactMidBytes = (char*)malloc(sizeof(char)*4);
	float value = oriData[0];
	floatToBytes(tdps->exactMidBytes, value);
	tdps->exactMidBytes_size = 4;
	
	int tmpOutSize;
	char *tmpByteData;
	convertTDPStoFlatBytes_float(tdps, &tmpByteData, &tmpOutSize);

	*newByteData = (char*)malloc(sizeof(char)*12); //for floating-point data (1+3+4+4)
	memcpy(*newByteData, tmpByteData, 12);
	*outSize = 12;
	free_TightDataPointStorageF(tdps);	
}

void SZ_compress_args_float_wRngeNoGzip(char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, float absErr_Bound, float rel_BoundRatio)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, rel_BoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
//		SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize);
		if(r5==0&&r4==0&&r3==0&&r2==0)
			SZ_compress_args_float_NoCkRngeNoGzip_1D(newByteData, oriData, r1, realPrecision, outSize);
		else if(r5==0&&r4==0&&r3==0)
			SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize);
		else if(r5==0&&r4==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r3, r2, r1, realPrecision, outSize);
		else if(r5==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r4*r3, r2, r1, realPrecision, outSize);
	}
}

void SZ_compress_args_float(char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, float absErr_Bound, float relBoundRatio)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		int tmpOutSize = 0;
		char* tmpByteData;
		if (r2==0)
			SZ_compress_args_float_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize);
		else
		if (r3==0)
			SZ_compress_args_float_NoCkRngeNoGzip_2D(&tmpByteData, oriData, r2, r1, realPrecision, &tmpOutSize);
		else
		if (r4==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r3, r2, r1, realPrecision, &tmpOutSize);
		else
		if (r5==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r4*r3, r2, r1, realPrecision, &tmpOutSize);			
		else
		{
			printf("Error: doesn't support 5 dimensions for now.\n");
			exit(0);
		}
		//TODO: call Gzip to do the further compression.
		*outSize = (int)zlib_compress2(tmpByteData, tmpOutSize, newByteData, gzipMode);
		free(tmpByteData);
	}
}

void SZ_decompress_args_float(float** newData, int r5, int r4, int r3, int r2, int r1, char* cmpBytes, int cmpSize)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	char* tmpBytes;
	int targetUncompressSize = dataLength <<2; //i.e., *4
	int tmpSize = 12;
	if(cmpSize!=12)
		tmpSize = zlib_uncompress3(cmpBytes, (ulong)cmpSize, &tmpBytes, (ulong)targetUncompressSize);
	else
		tmpBytes = cmpBytes;
	//tmpSize must be "much" smaller than dataLength
	
	//tmpBytes = cmpBytes;
	//ulong tmpSize = cmpSize;
	char* szTmpBytes = (char*)malloc(sizeof(char)*tmpSize);
	memcpy(szTmpBytes, tmpBytes, tmpSize);
	free(tmpBytes); //release useless memory
	
	//TODO: convert szTmpBytes to data array.
	TightDataPointStorageF* tdps;
	new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	if (dataLength == r1)
		getSnapshotData_float_1D(newData,r1,tdps);
	else
	if (dataLength == r1*r2)
		getSnapshotData_float_2D(newData,r2,r1,tdps);
	else
	if (dataLength == r1*r2*r3)
		getSnapshotData_float_3D(newData,r3,r2,r1,tdps);
	else
	if (dataLength == r1*r2*r3*r4)
		getSnapshotData_float_3D(newData,r4*r3,r2,r1,tdps);
	free_TightDataPointStorageF(tdps);
	free(szTmpBytes);
}
