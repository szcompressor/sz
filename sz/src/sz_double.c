/**
 *  @file sz_double.c
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
#include "zlib.h"
#include "rw.h"

void SZ_compress_args_double_NoCkRngeNoGzip_1D(char** newByteData, double *oriData, int dataLength, double realPrecision, int *outSize)
{
	int i;
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	
	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;
		
	double* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
		
	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_double(&esc, spaceFillingValue, reqExpo, realPrecision);
			
	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_DOUBLE);
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
	
	char preDataBytes[8];
	longToBytes_bigEndian(preDataBytes, 0);
	
	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;
			
	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	double medianValue = curExp->medianValue_d;
	double last3CmprsData[3] = {0};

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
				
	//add the first data	
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_double(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);		
		
	//add the second data
	type[1] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);			
	compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_double(last3CmprsData, vce->data);
	
	unsigned char state;
	double lcf, qcf;		
	double interval, checkRadius, curData;
	double pred;
	double predAbsErr;
	double min_pred, minErr, minIndex;
	int a = 0;		
	checkRadius = 255*realPrecision;
	interval = 2*realPrecision;

	for(i=2;i<dataLength;i++)
	{				
		//printf("%.30G\n",last3CmprsData[0]);
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
			listAdd_double(last3CmprsData, pred);
			continue;
		}
		
		//unpredictable data processing
		getExpSegment_fast(esc, i);
		curExp = esc->curExp;			
		type[i] = 0;
		curExp->unpredNum++;
		medianValue = curExp->medianValue_d;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;
		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		
		compressSingleDoubleValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,8);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
							
		listAdd_double(last3CmprsData, vce->data);
	}//end of for		
					
	cleanESwithNoUnpredictable(esc);		
		
	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageD* tdps;
			
	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			expSegmentsInBytes, expSegmentsInBytes_size, 
			resiBitLengthArray->array, resiBitLengthArray->size, realPrecision);
	
//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	//free(type);
	
	free_ExpSegmentConstructor(esc);
		
	//TODO: return bytes....
	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	free_DBA(exactMidByteArray);	
	
	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);	
}


/**
 * 
 * Note: @r1 is high dimension
 * 		 @r2 is low dimension 
 * */
void SZ_compress_args_double_NoCkRngeNoGzip_2D(char** newByteData, double *oriData, int r1, int r2, double realPrecision, int *outSize)
{
	int i,j;
	double pred1D, pred2D;
	double diff = 0.0;
	long itvNum = 0;
	double *P0, *P1;
		
	int dataLength = r1*r2;	
	
	P0 = (double*)malloc(r2*sizeof(double));
	P1 = (double*)malloc(r2*sizeof(double));
		
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	
	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;
		
	double* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
				
	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_double(&esc, spaceFillingValue, reqExpo, realPrecision);
			
	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_DOUBLE);
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
	
	char preDataBytes[8];
	longToBytes_bigEndian(preDataBytes, 0);
	
	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;
			
	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	double medianValue = curExp->medianValue_d;

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
			
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
		medianValue = curExp->medianValue_d;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;

		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,8);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
			medianValue = curExp->medianValue_d;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
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

		itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
			medianValue = curExp->medianValue_d;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}
									
		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		double *Pt;
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
	
	TightDataPointStorageD* tdps;
			
	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum, 
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
	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	free_DBA(exactMidByteArray);	
	
	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);	
}

void SZ_compress_args_double_NoCkRngeNoGzip_3D(char** newByteData, double *oriData, int r1, int r2, int r3, double realPrecision, int *outSize)
{
	int i,j,k;
	double pred1D, pred2D, pred3D;
	double diff = 0.0;
	long itvNum = 0;
	double *P0, *P1;

	int dataLength = r1*r2*r3;

	P0 = (double*)malloc(r2*r3*sizeof(double));
	P1 = (double*)malloc(r2*r3*sizeof(double));

	short reqExpo = getPrecisionReqLength_double(realPrecision);

	unsigned char* type = (unsigned char*) malloc(dataLength*sizeof(char));
	//type[dataLength]=0;

	double* spaceFillingValue = oriData; //
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;

	ExpSegmentConstructor* esc;
	new_ExpSegmentConstructor_double(&esc, spaceFillingValue, reqExpo, realPrecision);

	construct(esc, spaceFillingValue, dataLength, lengthBound, reqExpo, realPrecision, SZ_DOUBLE);
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

	char preDataBytes[8];
	longToBytes_bigEndian(preDataBytes, 0);

	getExpSegment_fast(esc, 0);
	ExpSegment* curExp = esc->curExp;
	curExp->unpredNum ++;

	int reqLength = curExp->reqLength;
	int reqBytesLength = curExp->reqBytesLength;
	int resiBitsLength = curExp->resiBitsLength;
	double medianValue = curExp->medianValue_d;

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
		medianValue = curExp->medianValue_d;
		reqLength = curExp->reqLength;
		reqBytesLength = curExp->reqBytesLength;
		resiBitsLength = curExp->resiBitsLength;

		addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
		compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,8);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
			medianValue = curExp->medianValue_d;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
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

		itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
			medianValue = curExp->medianValue_d;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
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

		itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
			medianValue = curExp->medianValue_d;
			reqLength = curExp->reqLength;
			reqBytesLength = curExp->reqBytesLength;
			resiBitsLength = curExp->resiBitsLength;

			addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}


	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			index = k*r2*r3+j;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
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

			itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;

				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
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

				itvNum = (long) (fabs(diff)/realPrecision) + 1;

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
					medianValue = curExp->medianValue_d;
					reqLength = curExp->reqLength;
					reqBytesLength = curExp->reqBytesLength;
					resiBitsLength = curExp->resiBitsLength;

					addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
					compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,8);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		double *Pt;
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

	TightDataPointStorageD* tdps;

	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum,
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
	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	free_DBA(exactMidByteArray);

	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);
}

void SZ_compress_args_double_withinRange(char** newByteData, double *oriData, int dataLength, int *outSize)
{
	TightDataPointStorageD* tdps = (TightDataPointStorageD*) malloc(sizeof(TightDataPointStorageD));
	tdps->rtypeArray = NULL;
	tdps->typeArray = NULL;
	tdps->leadNumArray = NULL;
	tdps->escBytes = NULL;
	tdps->residualMidBits = NULL;
	
	tdps->allSameData = 1;
	tdps->dataSeriesLength = dataLength;
	tdps->exactMidBytes = (char*)malloc(sizeof(char)*8);
	double value = oriData[0];
	doubleToBytes(tdps->exactMidBytes, value);
	tdps->exactMidBytes_size = 8;
	
	int tmpOutSize;
	char *tmpByteData;
	convertTDPStoFlatBytes_double(tdps, &tmpByteData, &tmpOutSize);

	*newByteData = (char*)malloc(sizeof(char)*16); //for floating-point data (1+3+4+4)
	memcpy(*newByteData, tmpByteData, 16);
	*outSize = 16;
	free_TightDataPointStorageD(tdps);	
}

void SZ_compress_args_double_wRngeNoGzip(char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double rel_BoundRatio)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	double valueRangeSize = 0, medianValue = 0;
	
	computeRangeSize_double(oriData, dataLength, &valueRangeSize, &medianValue);
	double realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, rel_BoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_double_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		if(r5==0&&r4==0&&r3==0&&r2==0)
			SZ_compress_args_double_NoCkRngeNoGzip_1D(newByteData, oriData, r1, realPrecision, outSize);
		else if(r5==0&&r4==0&&r3==0)
			SZ_compress_args_double_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize);
		else if(r5==0&&r4==0)
			SZ_compress_args_double_NoCkRngeNoGzip_3D(newByteData, oriData, r3, r2, r1, realPrecision, outSize);
	}
}

void SZ_compress_args_double(char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	double valueRangeSize = 0, medianValue = 0;
	
	computeRangeSize_double(oriData, dataLength, &valueRangeSize, &medianValue);
	double realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_double_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		int tmpOutSize = 0;
		char* tmpByteData;
		if (r2==0)
			SZ_compress_args_double_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize);
		else
		if (r3==0)
			SZ_compress_args_double_NoCkRngeNoGzip_2D(&tmpByteData, oriData, r2, r1, realPrecision, &tmpOutSize);
		else
		if (r4==0)
			SZ_compress_args_double_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r3, r2, r1, realPrecision, &tmpOutSize);
		else
		if (r5==0)
			SZ_compress_args_double_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r4*r3, r2, r1, realPrecision, &tmpOutSize);
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

void SZ_decompress_args_double(double** newData, int r5, int r4, int r3, int r2, int r1, char* cmpBytes, int cmpSize)
{
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	char* tmpBytes;
	int targetUncompressSize = dataLength <<3; //i.e., *8
	int tmpSize = 16;
	if(cmpSize!=16)
		tmpSize = zlib_uncompress3(cmpBytes, (ulong)cmpSize, &tmpBytes, (ulong)targetUncompressSize);
	else
		tmpBytes = cmpBytes;
	//tmpSize must be "much" smaller than dataLength
	
	//tmpBytes = cmpBytes;
	//ulong tmpSize = cmpSize;
	char* szTmpBytes = (char*)malloc(sizeof(char)*tmpSize);
	memcpy(szTmpBytes, tmpBytes, tmpSize);
	free(tmpBytes); //release useless memory
	//TODO: convert szTmpBytes to double array.
	TightDataPointStorageD* tdps;
	new_TightDataPointStorageD_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	if (dataLength == r1)
		getSnapshotData_double_1D(newData,r1,tdps);
	else
	if (dataLength == r1*r2)
		getSnapshotData_double_2D(newData,r2,r1,tdps);
	else
	if (dataLength == r1*r2*r3)
		getSnapshotData_double_3D(newData,r3,r2,r1,tdps);
	else
	if (dataLength == r1*r2*r3*r4)
		getSnapshotData_double_3D(newData,r4*r3,r2,r1,tdps);		
	
	free_TightDataPointStorageD(tdps);
	free(szTmpBytes);
}
