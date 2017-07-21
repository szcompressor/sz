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

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

float** create2DArray_float(size_t m, size_t n);
float*** create3DArray_float(size_t p, size_t m, size_t n);
void free2DArray_float(float** data, size_t m);
void free3DArray_float(float*** data, size_t p, size_t m);

double** create2DArray_double(size_t m, size_t n);
double*** create3DArray_double(size_t p, size_t m, size_t n);
void free2DArray_double(double** data, size_t m);
void free3DArray_double(double*** data, size_t p, size_t m);

size_t checkFileSize(char *srcFilePath, int *status);
unsigned char *readByteData(char *srcFilePath, size_t *byteLength, int *status);
double *readDoubleData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
float *readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
double *readDoubleData(char *srcFilePath, size_t *nbEle, int *status);
float *readFloatData(char *srcFilePath, size_t *nbEle, int *status);
void writeByteData(unsigned char *bytes, size_t outSize, char *tgtFilePath, int *status);
void writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status);
void writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status);
void writeData(void *data, int dataType, size_t nbEle, char *tgtFilePath, int *status);
void writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status);
void writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status);
void writeShortData(unsigned short *states, size_t stateLength, char *tgtFilePath, int *status);
unsigned short* readShortData(char *srcFilePath, size_t *dataLength, int *status);


void checkfilesizec_(char *srcFilePath, size_t *len, int *filesize);
void readbytefile_(char *srcFilePath, size_t *len, unsigned char *bytes, int *byteLength);
void readdoublefile_(char *srcFilePath, size_t *len, double *data, int *nbEle);
void readfloatfile_(char *srcFilePath, size_t *len, float *data, int *nbEle);
void writebytefile_(unsigned char *bytes, int *byteLength, char *tgtFilePath, size_t *len);
void writedoublefile_(double *data, size_t *nbEle, char *tgtFilePath, size_t *len);
void writefloatfile_(float *data, size_t *nbEle, char *tgtFilePath, size_t *len);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _IO_H  ----- */
