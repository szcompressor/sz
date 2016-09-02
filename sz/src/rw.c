/**
 *  @file rw.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief io interface for fortrance
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rw.h"
#include "sz.h"

int checkFileSize(char *srcFilePath)
{
	int filesize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		exit(1);
	}
	fseek(pFile, 0, SEEK_END);
    filesize = (int)ftell(pFile);
    return filesize;
}

char * readByteData(char *srcFilePath, ulong *byteLength)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = (int)ftell(pFile);
    fclose(pFile);
    
    char *byteBuf = (char *)malloc((*byteLength)*sizeof(char)); //sizeof(char)==1
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    return byteBuf;
}

double *readDoubleData(char *srcFilePath, ulong *nbEle)
{
	if(dataEndianType==sysEndianType)
	{
		double *daBuf = readDoubleData_systemEndian(srcFilePath, nbEle);
		return daBuf;
	}
	else
	{
		ulong i,j;
		
		ulong byteLength;
		char* bytes = readByteData(srcFilePath, &byteLength);
		double *daBuf = (double *)malloc(byteLength);
		*nbEle = byteLength/8;
		
		ldouble buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*8;
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

float *readFloatData(char *srcFilePath, ulong *nbEle)
{
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = readFloatData_systemEndian(srcFilePath, nbEle);
		return daBuf;
	}
	else
	{
		ulong i,j;
		
		ulong byteLength;
		char* bytes = readByteData(srcFilePath, &byteLength);
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;
		
		lfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

double * readDoubleData_systemEndian(char *srcFilePath, ulong *nbEle)
{
	ulong inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/8; //only support double in this version
    fclose(pFile);
    
    double *daBuf = (double *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(daBuf, 8, *nbEle, pFile);
    fclose(pFile);
    return daBuf;
}

float * readFloatData_systemEndian(char *srcFilePath, ulong *nbEle)
{
	ulong inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/4; 
    fclose(pFile);
    
    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		exit(0);
	}
    
    float *daBuf = (float *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    return daBuf;
}

void writeByteData(char *bytes, ulong byteLength, char *tgtFilePath)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
    
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
}

void writeDoubleData(double *data, ulong nbEle, char *tgtFilePath)
{
	ulong i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
    
    for(i = 0;i<nbEle;i++)
	{
		sprintf(s,"%.20G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
}

void writeFloatData(float *data, ulong nbEle, char *tgtFilePath)
{
	ulong i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
   
    for(i = 0;i<nbEle;i++)
	{
		//printf("i=%d\n",i);
		//printf("data[i]=%f\n",data[i]);
		sprintf(s,"%.30G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
}

void writeData(void *data, int dataType, int nbEle, char *tgtFilePath)
{
	if(dataType == SZ_FLOAT)
	{
		float* dataArray = (float *)data;
		writeFloatData(dataArray, nbEle, tgtFilePath);
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* dataArray = (double *)data;
		writeDoubleData(dataArray, nbEle, tgtFilePath);
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);	
	}
}
