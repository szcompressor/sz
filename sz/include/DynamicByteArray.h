/**
 *  @file DynamicByteArray.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Dynamic Byte Array.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _DynamicByteArray_H
#define _DynamicByteArray_H

typedef struct DynamicByteArray
{	
	char* array;
	int size;
	int capacity;
} DynamicByteArray;

void new_DBA(DynamicByteArray **dba, int cap);
void convertDBAtoBytes(DynamicByteArray *dba, char** bytes);
void free_DBA(DynamicByteArray *dba);
int getDBA_Data(DynamicByteArray *dba, int pos);
inline void addDBA_Data(DynamicByteArray *dba, char value);

#endif /* ----- #ifndef _DynamicByteArray_H  ----- */
