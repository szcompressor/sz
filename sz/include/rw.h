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

int checkFileSize(char *srcFilePath);
char *readByteData(char *srcFilePath, ulong *byteLength);
double *readDoubleData_systemEndian(char *srcFilePath, ulong *nbEle);
float *readFloatData_systemEndian(char *srcFilePath, ulong *nbEle);
double *readDoubleData(char *srcFilePath, ulong *nbEle);
float *readFloatData(char *srcFilePath, ulong *nbEle);
void writeByteData(char *bytes, ulong outSize, char *tgtFilePath);
void writeDoubleData(double *data, ulong nbEle, char *tgtFilePath);
void writeFloatData(float *data, ulong nbEle, char *tgtFilePath);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _IO_H  ----- */
