/**
 *  @file io.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole io interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _IO_H
#define _IO_H

#include <stdio.h>
#include <stdio.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

float** create2DArray_float(int m, int n);
float*** create3DArray_float(int p, int m, int n);
void free2DArray_float(float** data, int m);
void free3DArray_float(float*** data, int p, int m);

double** create2DArray_double(int m, int n);
double*** create3DArray_double(int p, int m, int n);
void free2DArray_double(double** data, int m);
void free3DArray_double(double*** data, int p, int m);

int checkFileSize(char *srcFilePath);
unsigned char *readByteData(char *srcFilePath, int *byteLength);
double *readDoubleData_systemEndian(char *srcFilePath, int *nbEle);
float *readFloatData_systemEndian(char *srcFilePath, int *nbEle);
double *readDoubleData(char *srcFilePath, int *nbEle);
float *readFloatData(char *srcFilePath, int *nbEle);
void writeByteData(unsigned char *bytes, int outSize, char *tgtFilePath);
void writeDoubleData(double *data, int nbEle, char *tgtFilePath);
void writeFloatData(float *data, int nbEle, char *tgtFilePath);
void writeData(void *data, int dataType, int nbEle, char *tgtFilePath);
void writeFloatData_inBytes(float *data, int nbEle, char* tgtFilePath);
void writeDoubleData_inBytes(double *data, int nbEle, char* tgtFilePath);
void writeShortData(unsigned short *states, int stateLength, char *tgtFilePath);
unsigned short* readShortData(char *srcFilePath, int *dataLength);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _IO_H  ----- */
