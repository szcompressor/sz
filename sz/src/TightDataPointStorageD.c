/**
 *  @file TightPointDataStorageD.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief The functions used to construct the tightPointDataStorage element for storing compressed bytes.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "TightDataPointStorageD.h"
#include "sz.h"
#include "ExpSegmentConstructor.h"

void new_TightDataPointStorageD_Empty(TightDataPointStorageD **this)
{
	*this = (TightDataPointStorageD*)malloc(sizeof(TightDataPointStorageD));
	(*this)->dataSeriesLength = 0;
	(*this)->allSameData = 0;
	(*this)->exactDataNum = 0;
	(*this)->reservedValue = 0;
	
	(*this)->rtypeArray = NULL;
	(*this)->rtypeArray_size = 0;
	
	(*this)->typeArray = NULL; //its size is dataSeriesLength/4 (or xxx/4+1) 
	(*this)->typeArray_size = 0;
	
	(*this)->leadNumArray = NULL; //its size is exactDataNum/4 (or exactDataNum/4+1)
	(*this)->leadNumArray_size = 0;
	
	(*this)->exactMidBytes = NULL;
	(*this)->exactMidBytes_size = 0;
	
	(*this)->escBytes = NULL;
	(*this)->escBytes_size = 0;
	
	(*this)->residualMidBits = NULL;
	(*this)->residualMidBits_size = 0;
}

void new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **this, char* flatBytes, int flatBytesLength)
{
	new_TightDataPointStorageD_Empty(this);
	int i, index = 0;
	char version[3];
	for (i = 0; i < 3; i++)
		version[i] = flatBytes[index++];
	if(checkVersion(version)!=1)
	{
		//wrong version
		printf("Wrong version: \nCompressed-data version (%d.%d.%d)\n",version[0], version[1], version[2]);
		printf("Current sz version: (%d.%d.%d)\n", versionNumber[0], versionNumber[1], versionNumber[2]);
		exit(0);
	}
	
	char dsLengthBytes[4];
	for (i = 0; i < 4; i++)
		dsLengthBytes[i] = flatBytes[index++];
		
	//TODO
	(*this)->dataSeriesLength = bytesToInt_bigEndian(dsLengthBytes);// 4
	char sameRByte = flatBytes[index++];
	int same = sameRByte & 0x01;
	
	if(same==1)
	{
		(*this)->allSameData = 1;
		int exactMidBytesLength = flatBytesLength - 3 - 4 -1;
		if(exactMidBytesLength>0)
			(*this)->exactMidBytes = (char*)malloc(sizeof(char)*exactMidBytesLength);
		else
			(*this)->exactMidBytes = NULL;
		for(i = 0;i<exactMidBytesLength;i++)
			(*this)->exactMidBytes[i] = flatBytes[index++];
		return;
	}
	else
		(*this)->allSameData = 0;
	
	int rtype_ = sameRByte & 0x02;
	if(rtype_!=0)
	{
		char rTypeLengthBytes[4];
		for(i = 0;i<4;i++) //4
			rTypeLengthBytes[i] = flatBytes[index++];
		(*this)->rtypeArray_size = bytesToInt_bigEndian(rTypeLengthBytes);//4		
	}
	else
		(*this)->rtypeArray_size = 0;	
	
	char exactLengthBytes[4];
	for (i = 0; i < 4; i++)
		exactLengthBytes[i] = flatBytes[index++];
	(*this)->exactDataNum = bytesToInt_bigEndian(exactLengthBytes);// 4

	char escByteLengthBytes[4];
	for (i = 0; i < 4; i++)
		escByteLengthBytes[i] = flatBytes[index++];
	(*this)->escBytes_size = bytesToInt_bigEndian(escByteLengthBytes);//4
	if((*this)->escBytes_size>0)
		(*this)->escBytes = (char*)malloc(sizeof(char)*(*this)->escBytes_size);
	else
		(*this)->escBytes = NULL;
	char exactMidBytesLengthBytes[4];
	for (i = 0; i < 4; i++)
		exactMidBytesLengthBytes[i] = flatBytes[index++];// 4
	(*this)->exactMidBytes_size = bytesToInt_bigEndian(exactMidBytesLengthBytes);// 4

	int typeArrayLength = 0;
	if (rtype_ != 0) {
		if((*this)->rtypeArray_size>0)
			(*this)->rtypeArray = (char*)malloc(sizeof(char)*(*this)->rtypeArray_size);
		else
			(*this)->rtypeArray = NULL;
		char reservedValueBytes[8];
		for (i = 0; i < 8; i++)
			reservedValueBytes[i] = flatBytes[index++];
		(*this)->reservedValue = bytesToDouble(reservedValueBytes);
	}

	int logicTypeLengthBitsNum = (*this)->dataSeriesLength * 2;
	if (logicTypeLengthBitsNum % 8 == 0)
	{
		(*this)->typeArray_size = logicTypeLengthBitsNum >> 3;//divided by 8
	}	
	else
	{
		(*this)->typeArray_size = (logicTypeLengthBitsNum >> 3) + 1;
	}
	if((*this)->typeArray_size>0)
		(*this)->typeArray = (char*)malloc(sizeof(char)*(*this)->typeArray_size);
	else
		(*this)->typeArray = NULL;

	int logicLeadNumBitsNum = (*this)->exactDataNum * 2;
	if (logicLeadNumBitsNum % 8 == 0)
	{
		(*this)->leadNumArray_size = logicLeadNumBitsNum >> 3;
	}
	else
	{
		(*this)->leadNumArray_size = (logicLeadNumBitsNum >> 3) + 1;
	}
	if((*this)->leadNumArray_size>0)
		(*this)->leadNumArray = (char*)malloc(sizeof(char)*(*this)->leadNumArray_size);
	else
		(*this)->leadNumArray = NULL;
	
	if((*this)->exactMidBytes_size>0)
		(*this)->exactMidBytes = (char*)malloc(sizeof(char)*(*this)->exactMidBytes_size);
	else
		(*this)->exactMidBytes = NULL;

	if ((*this)->rtypeArray != NULL) 
	{
		(*this)->residualMidBits_size = flatBytesLength - 3 - 4 - 1 - 4 - 4 - 4 - 4
				- 8 - (*this)->rtypeArray_size - (*this)->escBytes_size
				- (*this)->typeArray_size - (*this)->leadNumArray_size
				- (*this)->exactMidBytes_size;
		for (i = 0; i < (*this)->rtypeArray_size; i++)
			(*this)->rtypeArray[i] = flatBytes[index++];
	}
	else
	{
		(*this)->residualMidBits_size = flatBytesLength - 3 - 4 - 1 - 4 - 4
				- 4 - (*this)->escBytes_size - (*this)->typeArray_size
				- (*this)->leadNumArray_size - (*this)->exactMidBytes_size;
	}	
	if((*this)->residualMidBits_size>0)
		(*this)->residualMidBits = (char*)malloc(sizeof(char)*(*this)->residualMidBits_size);
	else
		(*this)->residualMidBits = NULL;
	for (i = 0; i < (*this)->escBytes_size; i++)
		(*this)->escBytes[i] = flatBytes[index++];
	for (i = 0; i < (*this)->typeArray_size; i++)
		(*this)->typeArray[i] = flatBytes[index++];
	for (i = 0; i < (*this)->leadNumArray_size; i++)
		(*this)->leadNumArray[i] = flatBytes[index++];
	for (i = 0; i < (*this)->exactMidBytes_size; i++)
	{
		(*this)->exactMidBytes[i] = flatBytes[index++];
	}	
	for (i = 0; i < (*this)->residualMidBits_size; i++)
		(*this)->residualMidBits[i] = flatBytes[index++];	
}


//TODO
void decompressDataSeries_double(double** data, int dataSeriesLength, TightDataPointStorageD* tdps) {
	int* leadNum;
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type;
	convertByteArray2IntArray_fast_2b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	ExpSegmentConstructor *esc;
	new_ExpSegmentConstructor_escbytes(&esc, 64, tdps->escBytes, tdps->escBytes_size, dataSeriesLength);

	char preBytes[8];
	char curBytes[8];
	
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum
	int curByteIndex = 0;	
	ExpSegment* curES;	
	int reqBytesLength, resiBitsLength, resiBits, leadingNum;	
	double medianValue, exactData;
	for (i = 0; i < dataSeriesLength; i++) {
		switch (type[i]) {
		case 0:
			getExpSegment_fast(esc,i);
			curES = esc->curExp;
			reqBytesLength = curES->reqBytesLength;
			resiBitsLength = curES->resiBitsLength;
			medianValue = curES->medianValue_d;

			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			int leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				char resiByte = (char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			double exactData = bytesToDouble(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
			break;
		case 1:
			(*data)[i] = (*data)[i - 1];
			break;
		case 2:
			(*data)[i] = 2 * (*data)[i - 1] - (*data)[i - 2];
			break;
		case 3:
			(*data)[i] = 3 * (*data)[i - 1] - 3
					* (*data)[i - 2] + (*data)[i - 3];
			break;
		default:
			printf("Error: type[%d] cannot be %d, but only 0, 1, 2, 3\n", i, type[i]);
			exit(0);
		}
	}
	free_ExpSegmentConstructor(esc);
	return;
}

//TODO
void getSnapshotData_double(double** data, int dataSeriesLength, TightDataPointStorageD* tdps) {
	
	int i;
	if (tdps->allSameData) {
		double value = bytesToDouble(tdps->exactMidBytes);
		*data = (double*)malloc(sizeof(double)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			
			decompressDataSeries_double(data, dataSeriesLength, tdps);
			return;
		} else {
			*data = (double*)malloc(sizeof(double)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			int count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			double* decmpData;
			decompressDataSeries_double(&decmpData, count, tdps);
			// insert the decompressed data
			int k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
		}
	}
}

/**
 * 
 * type's length == dataSeriesLength
 * exactMidBytes's length == exactMidBytes_size
 * leadNumIntArray's length == exactDataNum
 * escBytes's length == escBytes_size
 * resiBitLength's length == resiBitLengthSize
 * */
void new_TightDataPointStorageD(TightDataPointStorageD **this, int dataSeriesLength, int exactDataNum, 
		char* type, char* exactMidBytes, int exactMidBytes_size,
		char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		char* resiMidBits, int resiMidBits_size,
		char* escBytes, int escBytes_size,
		char* resiBitLength, int resiBitLengthSize) {
	*this = (TightDataPointStorageD *)malloc(sizeof(TightDataPointStorageD));
	(*this)->allSameData = 0;
	(*this)->dataSeriesLength = dataSeriesLength;
	(*this)->exactDataNum = exactDataNum;
	
	(*this)->rtypeArray = NULL;
	(*this)->rtypeArray_size = 0;
	
	(*this)->typeArray_size = convertIntArray2ByteArray_fast_2b(type, dataSeriesLength, &((*this)->typeArray));
	
	(*this)->exactMidBytes = exactMidBytes;
	(*this)->exactMidBytes_size = exactMidBytes_size;
	
	(*this)->leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &((*this)->leadNumArray));

	//(*this)->residualMidBits = resiMidBits;
	//(*this)->residualMidBits_size = resiMidBits_size;

	(*this)->escBytes = escBytes;
	(*this)->escBytes_size = escBytes_size;
	
	(*this)->residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiMidBits, resiBitLength, resiBitLengthSize, &((*this)->residualMidBits));
}

//TODO: convert TightDataPointStorageD to bytes...
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, char** bytes, int *size) 
{
	int i, k = 0;
	char rTypeLengthBytes[4];
	char dsLengthBytes[4];
	char exactLengthBytes[4];
	char escBytesLength[4];
	char exactMidBytesLength[4];
	char reservedValueBytes[8];
	
	intToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//4
	char sameByte = tdps->allSameData==1?(char)1:(char)0;
	
	if(tdps->allSameData==1)
	{
		int totalByteLength = 3 + 4 + 1 + tdps->exactMidBytes_size;
		*bytes = (char *)malloc(sizeof(char)*totalByteLength);
	
		for (i = 0; i < 3; i++)//3
			(*bytes)[k++] = versionNumber[i];
		for (i = 0; i < 4; i++)
			(*bytes)[k++] = dsLengthBytes[i];
		(*bytes)[k++] = sameByte;
		for (i = 0; i < tdps->exactMidBytes_size; i++)
			(*bytes)[k++] = tdps->exactMidBytes[i];
		
		*size = totalByteLength;
	}
	else if (tdps->rtypeArray == NULL) 
	{
		int residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
		int totalByteLength = 3 + 4 + 1 + 4 + 4 + 4 + tdps->escBytes_size
				+ tdps->typeArray_size + tdps->leadNumArray_size + tdps->exactMidBytes_size + residualMidBitsLength;

		*bytes = (char *)malloc(sizeof(char)*totalByteLength);

		for(i = 0;i<3;i++)//3
			(*bytes)[k++] = versionNumber[i];
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = dsLengthBytes[i];
		(*bytes)[k++] = sameByte;	//1			
		
		intToBytes_bigEndian(exactLengthBytes, tdps->exactDataNum);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactLengthBytes[i];
		
		intToBytes_bigEndian(escBytesLength, tdps->escBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = escBytesLength[i];
		
		intToBytes_bigEndian(exactMidBytesLength, tdps->exactMidBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactMidBytesLength[i];

		memcpy(&((*bytes)[k]), tdps->escBytes, tdps->escBytes_size);
		k += tdps->escBytes_size;
		memcpy(&((*bytes)[k]), tdps->typeArray, tdps->typeArray_size);
		k += tdps->typeArray_size;
		memcpy(&((*bytes)[k]), tdps->leadNumArray, tdps->leadNumArray_size);
		k += tdps->leadNumArray_size;
		memcpy(&((*bytes)[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
		k += tdps->exactMidBytes_size;		
		if(tdps->residualMidBits!=NULL)
		{
			memcpy(&((*bytes)[k]), tdps->residualMidBits, tdps->residualMidBits_size);
			k += tdps->residualMidBits_size;	
		}

		*size = totalByteLength;
	}
	else //the case with reserved value
	{
		int residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
		int totalByteLength = 3 + 4 + 1 + 4 + 4 + 4 + 4 + 8 + tdps->rtypeArray_size
		+ tdps->escBytes_size + tdps->typeArray_size + tdps->leadNumArray_size 
		+ tdps->exactMidBytes_size + residualMidBitsLength;

		sameByte = (char) (sameByte | 0x02); // 00000010, the second bit
												// denotes whether it is
												// with "reserved value"

		*bytes = (char*)malloc(sizeof(char)*totalByteLength);
		
		for(i = 0;i<3;i++)//3
			(*bytes)[k++] = versionNumber[i];
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = dsLengthBytes[i];
		(*bytes)[k++] = sameByte;						//1
		
		intToBytes_bigEndian(rTypeLengthBytes, tdps->rtypeArray_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = rTypeLengthBytes[i];	
		
		intToBytes_bigEndian(exactLengthBytes, tdps->exactDataNum);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactLengthBytes[i];
		
		intToBytes_bigEndian(escBytesLength, tdps->escBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = escBytesLength[i];
		
		intToBytes_bigEndian(exactMidBytesLength, tdps->exactMidBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactMidBytesLength[i];

		doubleToBytes(reservedValueBytes, tdps->reservedValue);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = reservedValueBytes[i];
		
		memcpy(&((*bytes)[k]), tdps->rtypeArray, tdps->rtypeArray_size);
		k += tdps->rtypeArray_size;		
		memcpy(&((*bytes)[k]), tdps->escBytes, tdps->escBytes_size);
		k += tdps->escBytes_size;
		memcpy(&((*bytes)[k]), tdps->typeArray, tdps->typeArray_size);
		k += tdps->typeArray_size;
		memcpy(&((*bytes)[k]), tdps->leadNumArray, tdps->leadNumArray_size);
		k += tdps->leadNumArray_size;
		memcpy(&((*bytes)[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
		k += tdps->exactMidBytes_size;		
		if(tdps->residualMidBits!=NULL)
		{
			memcpy(&((*bytes)[k]), tdps->residualMidBits, tdps->residualMidBits_size);
			k += tdps->residualMidBits_size;	
		}
		
		*size = totalByteLength;
	}
}

void free_TightDataPointStorageD(TightDataPointStorageD *tdps)
{
	if(tdps->rtypeArray!=NULL)
		free(tdps->rtypeArray);
	if(tdps->typeArray!=NULL)
		free(tdps->typeArray);
	if(tdps->leadNumArray!=NULL)
		free(tdps->leadNumArray);
	if(tdps->exactMidBytes!=NULL)
		free(tdps->exactMidBytes);
	if(tdps->escBytes!=NULL)
		free(tdps->escBytes);
	if(tdps->residualMidBits!=NULL)
		free(tdps->residualMidBits);
	free(tdps);
}
