/**
 *  @file DynamicIntArray.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Dynamic Int Array.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _DynamicIntArray_H
#define _DynamicIntArray_H

typedef struct DynamicIntArray
{	
	char* array; //char* (one byte) is enough, don't have to be int*
	int size;
	int capacity;
} DynamicIntArray;

void new_DIA(DynamicIntArray **dia, int cap);
void convertDIAtoInts(DynamicIntArray *dia, char **data);
void free_DIA(DynamicIntArray *dia);
int getDIA_Data(DynamicIntArray *dia, int pos);
inline void addDIA_Data(DynamicIntArray *dia, int value);

#endif /* ----- #ifndef _DynamicIntArray_H  ----- */
