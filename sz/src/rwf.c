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

void checkfilesizec_(char *srcFilePath, int *len, int *filesize)
{
	int i; 
	int status;
	char s[*len+1];
	for(i=0;i<*len;i++)
		s[i]=srcFilePath[i];
	s[*len]='\0';
	*filesize = checkFileSize(s, &status);
}

void readbytefile_(char *srcFilePath, size_t *len, unsigned char *bytes, int *byteLength)
{
	size_t i, ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';
    unsigned char *tmp_bytes = readByteData(s, byteLength, &ierr);
    memcpy(bytes, tmp_bytes, *byteLength);
    free(tmp_bytes);
}

void readdoublefile_(char *srcFilePath, size_t *len, double *data, int *nbEle)
{
	size_t i; 
	int ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';	
	double *tmp_data = readDoubleData(s, nbEle, &ierr);
	memcpy(data, tmp_data, *nbEle);
	free(tmp_data);
}

void readfloatfile_(char *srcFilePath, size_t *len, float *data, int *nbEle)
{
	size_t i; 
	int ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';
	float *tmp_data = readFloatData(s, nbEle, &ierr);
	memcpy(data, tmp_data, *nbEle);
	free(tmp_data);
}

void writebytefile_(unsigned char *bytes, int *byteLength, char *tgtFilePath, size_t *len)
{
	size_t i; 
	int ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';
	writeByteData(bytes, *byteLength, s, &ierr);
}

void writedoublefile_(double *data, size_t *nbEle, char *tgtFilePath, size_t *len)
{
	size_t i
	int ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';	
	writeDoubleData(data, *nbEle, s, &ierr);
}

void writefloatfile_(float *data, size_t *nbEle, char *tgtFilePath, size_t *len)
{
	size_t i; 
	int ierr;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';
	writeFloatData(data, *nbEle, s, &ierr);
}
