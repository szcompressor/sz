/**
 *  @file sz.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
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
//#include "CurveFillingCompressStorage.h"

int SZ_Init(char *configFilePath)
{
	char str[512]="", str2[512]="", str3[512]="";
	sz_cfgFile = configFilePath;
	int loadFileResult = SZ_LoadConf();
	if(loadFileResult==0)
		exit(0);
	sz_varset = (SZ_VarSet*)malloc(sizeof(SZ_VarSet));
	sz_varset->header = (SZ_Variable*)malloc(sizeof(SZ_Variable));
	sz_varset->lastVar = sz_varset->header;
	sz_varset->count = 0;
}

int computeDimension(int r5, int r4, int r3, int r2, int r1)
{
	int dimension;
	if(r1==0) 
	{
		dimension = 0;
	}
	else if(r2==0) 
	{
		dimension = 1;
	}
	else if(r3==0) 
	{
		dimension = 2;
	}
	else if(r4==0) 
	{
		dimension = 3;
	}
	else if(r5==0) 
	{
		dimension = 4;
	}
	else 
	{
		dimension = 5;
	}
	return dimension;	
}

int computeDataLength(int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;
	if(r1==0) 
	{
		dataLength = 0;
	}
	else if(r2==0) 
	{
		dataLength = r1;
	}
	else if(r3==0) 
	{
		dataLength = r1*r2;
	}
	else if(r4==0) 
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0) 
	{
		dataLength = r1*r2*r3*r4;
	}
	else 
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

void SZ_compress_args_float_NoCkRngeNoGzip(char** newByteData, float *oriData, int dataLength, float realPrecision, int *outSize)
{
	int i;
	short reqExpo = getPrecisionReqLength_float(realPrecision);
	
	char* type = (char*) malloc(dataLength*sizeof(char));
		
	float* spaceFillingValue = oriData; //
	float* values = spaceFillingValue;
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
	//TODO:
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
	//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	
	float medianValue = curExp->medianValue_f;
	
	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	
	float last3CmprsData[3];
	float* recentData;		
	float pred_1, pred_2, pred_3, curData, pred_1AbsErr, pred_2AbsErr, pred_3AbsErr;
	float min_pred, minErr;
	int a = 0;		
			
	listAdd_float(last3CmprsData, vce->data);
	for(i=1;i<dataLength;i++)
	{	
		getExpSegment_fast(esc, i);
		curExp = esc->curExp;
		if(i==1)
		{
			pred_1 = last3CmprsData[0];//last3CmprsData.get(0);
			pred_1AbsErr = fabs(pred_1 - values[1]);
			if(validPrediction_float(pred_1AbsErr, realPrecision))
			{
				type[i] = 1;
				listAdd_float(last3CmprsData, pred_1);
			}
			else
			{
				type[i] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleFloatValue(vce, values[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				
				listAdd_float(last3CmprsData, vce->data);
			}
		}
		else if(i==2)
		{
			pred_1 = last3CmprsData[0];//last3CmprsData.get(0);
			pred_2 = 2*last3CmprsData[0] - last3CmprsData[1];//2*last3CmprsData.get(0) - last3CmprsData.get(1);
			
			pred_1AbsErr = fabs(pred_1 - values[2]);
			pred_2AbsErr = fabs(pred_2 - values[2]);
			minErr = pred_1AbsErr<pred_2AbsErr ? pred_1AbsErr : pred_2AbsErr;
			if(validPrediction_float(minErr, realPrecision))
			{
				if(pred_1AbsErr <= pred_2AbsErr)
				{
					type[2] = 1;
					listAdd_float(last3CmprsData, pred_1);
				}
				else
				{
					type[2] = 2;
					listAdd_float(last3CmprsData, pred_2);
				}
			}
			else
			{
				type[2] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleFloatValue(vce, values[2], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				
				listAdd_float(last3CmprsData, vce->data);
			}					
		}
		else //i>=3
		{
			recentData = last3CmprsData;
			
			pred_1 = recentData[0];
			pred_2 = 2*recentData[0] - recentData[1];
			pred_3 = 3*recentData[0] - 3*recentData[1] + recentData[2];
			
			curData = values[i];
			pred_1AbsErr = fabs(pred_1 - curData);
			pred_2AbsErr = fabs(pred_2 - curData);	
			pred_3AbsErr = fabs(pred_3 - curData);
			min_pred = pred_1AbsErr<pred_2AbsErr?pred_1AbsErr:pred_2AbsErr;	
			minErr = min_pred<pred_3AbsErr?min_pred:pred_3AbsErr;
			
			if(validPrediction_float(minErr, realPrecision))
			{
				if(pred_1AbsErr<=pred_2AbsErr&&pred_1AbsErr<=pred_3AbsErr)
				{
					type[i] = 1;
					listAdd_float(last3CmprsData, pred_1);
				}
				else if(pred_2AbsErr<pred_1AbsErr&&pred_2AbsErr<=pred_3AbsErr)
				{
					type[i] = 2;
					listAdd_float(last3CmprsData, pred_2);
				}
				else //pred_3 is the smallest (pred_3<pred_2 and pred_3<pred_1)
				{
					type[i] = 3;
					listAdd_float(last3CmprsData, pred_3);
				}					
			}
			else
			{
				type[i] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_f;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
								
				listAdd_float(last3CmprsData, vce->data);
			}
		}//end of if() else....
	}//end of for			
											
	cleanESwithNoUnpredictable(esc);		
	
	//char* resiBitLengthBytes;
	//convertDBAtoBytes(resiBitLengthArray, &resiBitLengthBytes);
	
	//char* exactLeadNumIntArray;
	//convertDIAtoInts(exactLeadNumArray, &exactLeadNumIntArray);

	
	//char* exactMidBytes;
	//convertDBAtoBytes(exactMidByteArray, &exactMidBytes);
	
	//char* resiBitIntArray;
	//convertDIAtoInts(resiBitArray, &resiBitIntArray);
		
	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			expSegmentsInBytes, expSegmentsInBytes_size, 
			resiBitLengthArray->array, resiBitLengthArray->size);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	
	free_ExpSegmentConstructor(esc);
	
	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	free_DBA(exactMidByteArray);	
	//*outSize= tmpOutSize;
	//*newByteData = tmpByteData;
	
	//free(preDataBytes);  //no need to release preDataBytes, because it's in dvce or lce?
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
		SZ_compress_args_float_NoCkRngeNoGzip(newByteData, oriData, dataLength, realPrecision, outSize);
	}
}

void SZ_compress_args_float(char** newByteData, float *oriData, 
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
		int tmpOutSize;
		char* tmpByteData;
		SZ_compress_args_float_NoCkRngeNoGzip(&tmpByteData, oriData, dataLength, realPrecision, &tmpOutSize);
		
		//TODO: call Gzip to do the further compression.
		*outSize = (int)zlib_compress2(tmpByteData, tmpOutSize, newByteData, gzipMode);
		free(tmpByteData);
	}
}

void SZ_compress_args_double_NoCkRngeNoGzip(char** newByteData, double *oriData, int dataLength, double realPrecision, int *outSize)
{
	int i;
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	char* type = (char*) malloc(dataLength*sizeof(char));
		
	double* spaceFillingValue = oriData; //
	double* values = spaceFillingValue;
	int lengthBound = dataLength/maxSegmentNum;
	if(lengthBound < 128)
		lengthBound = 128;
	//TODO:
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
	//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
	addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
	
	double medianValue = curExp->medianValue_d;
	
	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	
	double last3CmprsData[3];
	double* recentData;		
	double pred_1, pred_2, pred_3, curData, pred_1AbsErr, pred_2AbsErr, pred_3AbsErr;
	double min_pred, minErr;
			
	listAdd_double(last3CmprsData, vce->data);
	for(i=1;i<dataLength;i++)
	{
		getExpSegment_fast(esc, i);
		curExp = esc->curExp;
		if(i==1)
		{
			pred_1 = last3CmprsData[0];
			pred_1AbsErr = fabs(pred_1 - values[1]);
			if(validPrediction_double(pred_1AbsErr, realPrecision))
			{
				type[i] = 1;
				listAdd_double(last3CmprsData, pred_1);
			}
			else
			{
				type[i] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleDoubleValue(vce, values[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				
				listAdd_double(last3CmprsData, vce->data);
			}
		}
		else if(i==2)
		{
			pred_1 = last3CmprsData[0];//last3CmprsData.get(0);
			pred_2 = 2*last3CmprsData[0] - last3CmprsData[1];//2*last3CmprsData.get(0) - last3CmprsData.get(1);
			
			pred_1AbsErr = fabs(pred_1 - values[2]);
			pred_2AbsErr = fabs(pred_2 - values[2]);
			minErr = pred_1AbsErr<pred_2AbsErr ? pred_1AbsErr : pred_2AbsErr;
			if(validPrediction_double(minErr, realPrecision))
			{
				if(pred_1AbsErr <= pred_2AbsErr)
				{
					type[2] = 1;
					listAdd_double(last3CmprsData, pred_1);
				}
				else
				{
					type[2] = 2;
					listAdd_double(last3CmprsData, pred_2);
				}
			}
			else
			{
				type[2] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleDoubleValue(vce, values[2], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				
				listAdd_double(last3CmprsData, vce->data);
			}					
		}
		else //i>=3
		{
			recentData = last3CmprsData;
			
			pred_1 = recentData[0];
			pred_2 = 2*recentData[0] - recentData[1];
			pred_3 = 3*recentData[0] - 3*recentData[1] + recentData[2];
			
			curData = values[i];
			pred_1AbsErr = fabs(pred_1 - curData);
			pred_2AbsErr = fabs(pred_2 - curData);	
			pred_3AbsErr = fabs(pred_3 - curData);
			min_pred = pred_1AbsErr<pred_2AbsErr?pred_1AbsErr:pred_2AbsErr;	
			minErr = min_pred<pred_3AbsErr?min_pred:pred_3AbsErr;
			
			if(validPrediction_double(minErr, realPrecision))
			{
				if(pred_1AbsErr<=pred_2AbsErr&&pred_1AbsErr<=pred_3AbsErr)
				{
					type[i] = 1;
					listAdd_double(last3CmprsData, pred_1);
				}
				else if(pred_2AbsErr<pred_1AbsErr&&pred_2AbsErr<=pred_3AbsErr)
				{
					type[i] = 2;
					listAdd_double(last3CmprsData, pred_2);
				}
				else //pred_3 is the smallest (pred_3<pred_2 and pred_3<pred_1)
				{
					type[i] = 3;
					listAdd_double(last3CmprsData, pred_3);
				}					
			}
			else
			{
				type[i] = 0;
				curExp->unpredNum++;
				medianValue = curExp->medianValue_d;
				reqLength = curExp->reqLength;
				reqBytesLength = curExp->reqBytesLength;
				resiBitsLength = curExp->resiBitsLength;
				//addDBA_Data(reqByteLengthArray, (char)reqBytesLength);
				addDBA_Data(resiBitLengthArray, (char)resiBitsLength);
				
				compressSingleDoubleValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
								
				listAdd_double(last3CmprsData, vce->data);
			}
		}//end of if() else....
	}//end of for			
											
	cleanESwithNoUnpredictable(esc);		
	
	char* resiBitLengthBytes;
	convertDBAtoBytes(resiBitLengthArray, &resiBitLengthBytes);
	
	char* exactLeadNumIntArray;
	convertDIAtoInts(exactLeadNumArray, &exactLeadNumIntArray);

	
	char* exactMidBytes;
	convertDBAtoBytes(exactMidByteArray, &exactMidBytes);
	
	char* resiBitIntArray;
	convertDIAtoInts(resiBitArray, &resiBitIntArray);
		
	char* expSegmentsInBytes;
	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageD* tdps;
			
	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum, 
			type, exactMidBytes, exactMidByteArray->size,  
			exactLeadNumIntArray,  
			resiBitIntArray, resiBitArray->size, 
			expSegmentsInBytes, expSegmentsInBytes_size, 
			resiBitLengthBytes, resiBitLengthArray->size);
	
	//free memory
	free_DBA(exactMidByteArray);	
	free_DBA(resiBitLengthArray);
	free_DIA(resiBitArray);
	free_DIA(exactLeadNumArray);
	
	free_ExpSegmentConstructor(esc);
	
	//TODO: return bytes....
	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	//*outSize= tmpOutSize;
	//*newByteData = tmpByteData;
	
	//free(preDataBytes);  //no need to release preDataBytes, because it's in dvce or lce?
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
int errBoundMode, double absErr_Bound, double relBoundRatio)
{

	int dataLength = computeDataLength(r5,r4,r3,r2,r1);

	double valueRangeSize, medianValue;

	computeRangeSize_double(oriData, dataLength, &valueRangeSize, &medianValue);
	double realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_double_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		SZ_compress_args_double_NoCkRngeNoGzip(newByteData, oriData, dataLength, realPrecision, outSize);
	}	
}

void SZ_compress_args_double(char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{

	int dataLength = computeDataLength(r5,r4,r3,r2,r1);

	double valueRangeSize, medianValue;

	computeRangeSize_double(oriData, dataLength, &valueRangeSize, &medianValue);
	double realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_double_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		int tmpOutSize;
		char* tmpByteData;
		SZ_compress_args_double_NoCkRngeNoGzip(&tmpByteData, oriData, dataLength, realPrecision, &tmpOutSize);
		
		//TODO: call Gzip to do the further compression.
		*outSize = (int)zlib_compress2(tmpByteData, tmpOutSize, newByteData, gzipMode);
		free(tmpByteData);
	}	
}

/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
char* SZ_compress_args(int dataType, void *data, int *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	//TODO
	if(dataType==SZ_FLOAT)
	{
		char *newByteData;
		
		SZ_compress_args_float(&newByteData, (float *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, (float)absErrBound, (float)relBoundRatio);
		
		return newByteData;
	}
	else if(dataType==SZ_DOUBLE)
	{
		char *newByteData;
		SZ_compress_args_double(&newByteData, (double *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return newByteData;
	}
	else
	{
		printf("Error: dataType can only be SZ_FLOAT or SZ_DOUBLE.\n");
		exit(1);
		return NULL;
	}
}

int SZ_compress_args2(int dataType, void *data, char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	char* bytes = SZ_compress_args(dataType, data, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
    memcpy(compressed_bytes, bytes, *outSize);
    free(bytes); 
	return 0;
}

char *SZ_compress(int dataType, void *data, int *outSize, int r5, int r4, int r3, int r2, int r1)
{	
	char *newByteData = SZ_compress_args(dataType, data, outSize, errorBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	
	return newByteData;
}

//////////////////
/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param		reservedValue  the reserved value
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;
	char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;	
}

int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	char* bytes = SZ_compress_rev_args(dataType, data, reservedValue, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	memcpy(compressed_bytes, bytes, *outSize);
	free(bytes); //free(bytes) is removed , because of dump error at MIRA system (PPC architecture), fixed?
	return 0;
}

char *SZ_compress_rev(int dataType, void *data, void *reservedValue, int *outSize, int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;

	char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);	
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;
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
	char* szTmpBytes = (char*)malloc(sizeof(char)*tmpSize);
	memcpy(szTmpBytes, tmpBytes, tmpSize);
	free(tmpBytes); //release useless memory
	
	//TODO: convert szTmpBytes to data array.
	TightDataPointStorageF* tdps;
	new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	getSnapshotData_float(newData,dataLength,tdps);
	
	free_TightDataPointStorageF(tdps);
	free(szTmpBytes);
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
	getSnapshotData_double(newData,dataLength,tdps);
	
	free_TightDataPointStorageD(tdps);
}

void *SZ_decompress(int dataType, char *bytes, int byteLength, int r5, int r4, int r3, int r2, int r1)
{
	if(dataType == SZ_FLOAT)
	{
		float *newFloatData;
		SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newFloatData;
	}
	else if(dataType == SZ_DOUBLE)
	{
		double *newDoubleData;
		SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newDoubleData;
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);		
	}
}

int SZ_decompress_args(int dataType, char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1)
{
	int i;
	int nbEle = computeDataLength(r5,r4,r3,r2,r1);
	if(dataType == SZ_FLOAT)
	{
		float* data = (float *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		float* data_array = (float *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(float));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];	
		free(data); //this free operation seems to not work with BlueG/Q system.
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* data = (double *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		double* data_array = (double *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(double));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];
		free(data); //this free operation seems to not work with BlueG/Q system.
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);				
	}

	return nbEle;
}

/*-----------------------------------batch data compression--------------------------------------*/

void filloutDimArray(int* dim, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
		dim[0] = r1;
	else if(r3==0)
	{
		dim[0] = r2;
		dim[1] = r1;
	}
	else if(r4==0)
	{
		dim[0] = r3;
		dim[1] = r2;
		dim[2] = r1;
	}
	else if(r5==0)
	{
		dim[0] = r4;
		dim[1] = r3;
		dim[2] = r2;
		dim[3] = r1;
	}
	else
	{
		dim[0] = r5;
		dim[1] = r4;
		dim[2] = r3;
		dim[3] = r2;
		dim[4] = r1;		
	}
}

char* SZ_batch_compress(int *outSize)
{
	int dataLength;
	DynamicByteArray* dba; 
	new_DBA(&dba, 32768);
	
	//number of variables
	int varCount = sz_varset->count;
	char  countBufBytes[4];
	intToBytes_bigEndian(countBufBytes, varCount);
	//add only the lats two bytes, i.e., the maximum # variables is supposed to be less than 32768
	addDBA_Data(dba, countBufBytes[2]);
	addDBA_Data(dba, countBufBytes[3]);
	
	int i, j, k = 0;
	SZ_Variable* p = sz_varset->header->next;
	while(p!=NULL)
	{
		if(p->dataType==SZ_FLOAT)
		{
			char *newByteData;
			int outSize;
			SZ_compress_args_float_wRngeNoGzip(&newByteData, (float *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, (float)(p->absErrBound), (float)(p->relBoundRatio));
			
			p->compressedBytes = newByteData;
		}
		else if(p->dataType==SZ_DOUBLE)
		{
			char *newByteData;
			int outSize;
			SZ_compress_args_double_wRngeNoGzip(&newByteData, (double *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, p->absErrBound, p->relBoundRatio);
			
			p->compressedBytes = newByteData;
		}
		
		//TODO: copy metadata information (totally 1+4*x+4+(1+y) bytes)
		//byte format of variable storage: meta 1 bytes (000000xy) x=0 SZ/1 HZ, y=0 float/1 double ; 
		//4*x: x integers indicating the dimension sizes; 
		//compressed length (int) 4 bytes;
		//1+y: 1 means the length of varname, y bytes represents the varName string. 
		int meta = 0;
		if(p->dataType == SZ_FLOAT)
			meta = 2; //10
		else if(p->dataType == SZ_DOUBLE)
			meta = 3; //11
		
		//keep dimension information
		int dimNum = computeDimension(p->r5, p->r4, p->r3, p->r2, p->r1);
		int dimSize[dimNum];
		memset(dimSize, 0, dimNum*sizeof(int));
		meta = meta | dimNum << 2; //---aaabc: aaa indicates dim, b indicates HZ, c indicates dataType
		
		addDBA_Data(dba, (char)meta);
		
		filloutDimArray(dimSize, p->r5, p->r4, p->r3, p->r2, p->r1);
		
		for(j=0;j<dimNum;j++)
		{
			intToBytes_bigEndian(countBufBytes, dimSize[j]);
			for(i = 0;i<4;i++)
				addDBA_Data(dba, countBufBytes[i]);
		}
			 
		//Keep compressed size information	 
		intToBytes_bigEndian(countBufBytes, p->compressedSize);
		
		for(i = 0;i<4;i++)
			addDBA_Data(dba, countBufBytes[i]);			 
			 
		//Keep varName information	 
		char varNameLength = (char)strlen(p->varName);
		addDBA_Data(dba, varNameLength);
		memcpyDBA_Data(dba, p->varName, varNameLength);
			 
		//copy the compressed stream
		memcpyDBA_Data(dba, p->compressedBytes, p->compressedSize);
		free(p->compressedBytes);
		p->compressedBytes = NULL;
			 
		p = p->next;
	}
	
	char* tmpFinalCompressedBytes;
	convertDBAtoBytes(dba, &tmpFinalCompressedBytes);
	
	char* tmpCompressedBytes2;
	int tmpGzipSize = (int)zlib_compress2(tmpFinalCompressedBytes, dba->size, &tmpCompressedBytes2, gzipMode);
	free(tmpFinalCompressedBytes);
	
	char* finalCompressedBytes = (char*) malloc(sizeof(char)*(4+tmpGzipSize));
	
	intToBytes_bigEndian(countBufBytes, dba->size);
	
	memcpy(finalCompressedBytes, countBufBytes, 4);
	
	memcpy(&(finalCompressedBytes[4]), tmpCompressedBytes2, tmpGzipSize);
	free(tmpCompressedBytes2);
	
	*outSize = 4+tmpGzipSize;
	return finalCompressedBytes;
}

SZ_VarSet* SZ_batch_decompress(char* compressedStream, int compressedLength)
{
	int i, j, k = 0;
	char intByteBuf[4];
	
	//get target decompression size for Gzip (zlib)
	intByteBuf[0] = compressedStream[0];
	intByteBuf[1] = compressedStream[1];
	intByteBuf[2] = compressedStream[2];
	intByteBuf[3] = compressedStream[3];
	
	int targetUncompressSize = bytesToInt_bigEndian(intByteBuf);
	
	//Gzip decompression
	char* gzipDecpressBytes;
	int gzipDecpressSize = zlib_uncompress3(&(compressedStream[4]), (ulong)compressedLength, &gzipDecpressBytes, (ulong)targetUncompressSize);
	
	if(gzipDecpressSize!=targetUncompressSize)
	{
		printf("Error: incorrect decompression in zlib_uncompress3: gzipDecpressSize!=targetUncompressSize\n");
		exit(0);
	}
	
	//Start analyzing the byte stream for further decompression	
	intByteBuf[0] = 0;
	intByteBuf[1] = 0; 
	intByteBuf[2] = gzipDecpressBytes[k++];
	intByteBuf[3] = gzipDecpressBytes[k++];
	
	int varCount = bytesToInt_bigEndian(intByteBuf);	
		
	int varNum = sz_varset->count;
	int dataLength, cpressedLength;
	
	if(varNum==0)
	{
		SZ_Variable* lastVar = sz_varset->header; 
		if(lastVar->next!=NULL)
		{
			printf("Error: sz_varset.count is inconsistent with the number of variables in sz_varset->header.\n");
			exit(0);
		}
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				exit(0);
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			int dimSize[5];
			memset(dimSize, 0, 5*sizeof(int));
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
				k+=4;
				dimSize[start_dim+j] = bytesToInt_bigEndian(intByteBuf);
			}	
			
			//get compressed length
			memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
			k+=4;
			cpressedLength = bytesToInt_bigEndian(intByteBuf);	
					
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
		
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);				
				getSnapshotData_float(&newData,dataLength,tdps);
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);
				free_TightDataPointStorageF(tdps);			
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);				
				getSnapshotData_double(&newData,dataLength,tdps);
				free_TightDataPointStorageD(tdps);
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);				
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				exit(0);
			}
	
			k+=cpressedLength;
		}
	}
	else if(varNum!=varCount)
	{
		printf("Error: Inconsistency of the # variables between sz_varset and decompressed stream.\n");
		printf("Note sz_varset.count should be either 0 or the correct number of variables stored in the decompressed stream.\n");
		printf("Currently, sz_varset.count=%d, expected number of variables = %d\n", varNum, varCount);
		exit(0);
	}
	else //if(varNum>0 && varNum==varCount)
	{
		SZ_Variable* p = sz_varset->header; 
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				exit(0);
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			int dimSize[5];
			memset(dimSize, 0, 5*sizeof(int));
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
				k+=4;
				dimSize[start_dim+j] = bytesToInt_bigEndian(intByteBuf);
			}
			
			//get compressed length
			memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
			k+=4;
			cpressedLength = bytesToInt_bigEndian(intByteBuf);	
			
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
			
			SZ_Variable* curVar = p->next;
			
			int checkVarName = strcmp(varNameString, curVar->varName);
			if(checkVarName!=0)
			{
				printf("Error: the varNames in the compressed stream are inconsistent with those in the sz_varset\n");
				exit(100);
			}
						
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);					
				getSnapshotData_float(&newData,dataLength,tdps);
				free_TightDataPointStorageF(tdps);	
				curVar->data = newData;
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);					
				getSnapshotData_double(&newData,dataLength,tdps);
				free_TightDataPointStorageD(tdps);	
				curVar->data = newData;
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				exit(0);
			}		
			
			k+=cpressedLength;
			p = p->next;
		}
	}
	
	return sz_varset;
}

void SZ_Finalize()
{
	//TODO
}
