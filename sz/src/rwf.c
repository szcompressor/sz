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

void checkfilesizec_(char *srcFilePath, ulong *len, ulong *filesize)
{
	ulong i;
	char s[*len+1];
	for(i=0;i<*len;i++)
		s[i]=srcFilePath[i];
	s[*len]='\0';	
	*filesize = checkFileSize(s);
}

void readbytefile_(char *srcFilePath, ulong *len, char *bytes, ulong *byteLength)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';
    char *tmp_bytes = readByteData(s, byteLength);
    memcpy(bytes, tmp_bytes, *byteLength);
}

void readdoublefile_(char *srcFilePath, ulong *len, double *data, ulong *nbEle)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';	
	double *tmp_data = readDoubleData(s, nbEle);
	memcpy(data, tmp_data, *nbEle);
}

void readfloatfile_(char *srcFilePath, ulong *len, float *data, ulong *nbEle)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=srcFilePath[i];
    s[*len]='\0';
	float *tmp_data = readFloatData(s, nbEle);
	memcpy(data, tmp_data, *nbEle);
}

void writebytefile_(char *bytes, ulong *byteLength, char *tgtFilePath, ulong *len)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';
	writeByteData(bytes, *byteLength, s);
}

void writedoublefile_(double *data, ulong *nbEle, char *tgtFilePath, ulong *len)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';	
	writeDoubleData(data, *nbEle, s);
}

void writefloatfile_(float *data, ulong *nbEle, char *tgtFilePath, ulong *len)
{
	ulong i;
    char s[*len+1];
    for(i=0;i<*len;i++)
        s[i]=tgtFilePath[i];
    s[*len]='\0';
	writeFloatData(data, *nbEle, s);
}
