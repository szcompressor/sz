/**
 *  @file TypeManager.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief TypeManager is used to manage the type array: parsing of the bytes and other types in between.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include "DynamicByteArray.h"

/**
 * little endian
 * [01|10|11|00|....]-->[01|10|11|00][....]
 * @param timeStepType
 * @return
 */
int convertIntArray2ByteArray_fast_2b(char* timeStepType, int timeStepTypeLength, char **result)
{
	int i, j, byteLength = 0;
	if(timeStepTypeLength%4==0)
		byteLength = timeStepTypeLength*2/8;
	else
		byteLength = timeStepTypeLength*2/8+1;
	if(byteLength>0)
		*result = (char*)malloc(byteLength*sizeof(char));
	else
		*result = NULL;
	int n = 0;
	for(i = 0;i<byteLength;i++)
	{
		int tmp = 0;
		for(j = 0;j<4&&n<timeStepTypeLength;j++)
		{
			int type = timeStepType[n];
			switch(type)
			{
			case 0: 
				
				break;
			case 1:
				tmp = (tmp | (1 << (6-j*2)));
				break;
			case 2:
				tmp = (tmp | (2 << (6-j*2)));
				break;
			case 3:
				tmp = (tmp | (3 << (6-j*2)));
				break;
			default:
				printf("Error: wrong timestep type...: type[%d]=%d\n", n, type);
				exit(0);
			}
			n++;
		}
		(*result)[i] = (char)tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_2b(int stepLength, char* byteArray, int byteArrayLength, int **intArray)
{
	if(stepLength > byteArrayLength*4)
	{
		printf("Error: stepLength > byteArray.length*4\n");
		printf("stepLength=%d, byteArray.length=%d\n", stepLength, byteArrayLength);
		exit(0);
	}
	if(stepLength>0)
		*intArray = (int*)malloc(stepLength*sizeof(int));
	else
		*intArray = NULL;
	int i, n = 0;

	for (i = 0; i < byteArrayLength; i++) {
		int tmp = byteArray[i];
		(*intArray)[n++] = (tmp & 0xC0) >> 6;
		if(n==stepLength)
			break;
		(*intArray)[n++] = (tmp & 0x30) >> 4;
		if(n==stepLength)
			break;
		(*intArray)[n++] = (tmp & 0x0C) >> 2;
		if(n==stepLength)
			break;
		(*intArray)[n++] = tmp & 0x03;
		if(n==stepLength)
			break;
	}
}

int getLeftMovingSteps(int k, char resiBitLength)
{
	return 8 - k%8 - resiBitLength;
}

/**
 * 
 * @param timeStepType is the resiMidBits
 * @param resiBitLength is the length of resiMidBits for each element, (the number of resiBitLength == the # of unpredictable elements
 * @return
 */
int convertIntArray2ByteArray_fast_dynamic(char* timeStepType, char* resiBitLength, int resiBitLengthLength, char **bytes)
{
	int i = 0, j = 0, k = 0, value;
	DynamicByteArray* dba;
	new_DBA(&dba, 1024);
	int tmp = 0, leftMovSteps = 0;
	for(j = 0;j<resiBitLengthLength;j++)
	{
		char rbl = resiBitLength[j];
		if(rbl==0)
			continue;
		value = timeStepType[i];
		leftMovSteps = getLeftMovingSteps(k, rbl);
		if(leftMovSteps < 0)
		{
			tmp = tmp | (value >> (-leftMovSteps));
			addDBA_Data(dba, (char)tmp);
			tmp = 0 | (value << (8+leftMovSteps));
		}
		else if(leftMovSteps > 0)
		{
			tmp = tmp | (value << leftMovSteps);
		}
		else //==0
		{
			tmp = tmp | value;
			addDBA_Data(dba, (char)tmp);
			tmp = 0;
		}
		i++;
		k += rbl;
	}
	if(leftMovSteps != 0)
		addDBA_Data(dba, (char)tmp);
	convertDBAtoBytes(dba, bytes);
	int size = dba->size;
	free_DBA(dba);
	return size;
}

int computeBitNumRequired(int dataLength)
{
	return 32 - numberOfLeadingZeros_Int(dataLength);
}

void decompressBitArraybySimpleLZ77(int** result, char* bytes, int bytesLength, int totalLength, int validLength)
{
	int pairLength = (bytesLength*8)/(validLength+1);
	int tmpLength = pairLength*2;
	int tmpResult[tmpLength];
	int i, j, k = 0;
	for(i = 0;i<tmpLength;i+=2)
	{
		int outIndex = k/8;
		int innerIndex = k%8;

		char curByte = bytes[outIndex];
		tmpResult[i] = (curByte >> (8-1-innerIndex)) & 0x01;
		k++;
		
		int numResult = extractBytes(bytes, k, validLength);
		
		tmpResult[i+1] = numResult;
		k = k + validLength;
	}
	
	*result = (int*)malloc(sizeof(int)*totalLength);
	k = 0;
	for(i = 0;i<tmpLength;i=i+2)
	{
		int state = tmpResult[i];
		int num = tmpResult[i+1];
		for(j = 0;j<num;j++)
			(*result)[k++] = state;
	}
}
