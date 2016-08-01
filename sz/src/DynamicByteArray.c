/**
 *  @file DynamicByteArray.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief Dynamic Byte Array
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "DynamicByteArray.h"

void new_DBA(DynamicByteArray **dba, int cap) {
		*dba = (DynamicByteArray *)malloc(sizeof(DynamicByteArray));
        (*dba)->size = 0;
        (*dba)->capacity = cap;
        (*dba)->array = (char*)malloc(sizeof(char)*cap);
    }

void convertDBAtoBytes(DynamicByteArray *dba, char** bytes)
{
	int size = dba->size;
	if(size>0)
		*bytes = (char*)malloc(size * sizeof(char));
	else
		*bytes = NULL;
	memcpy(*bytes, dba->array, size*sizeof(char));	
}

void free_DBA(DynamicByteArray *dba)
{
	free(dba->array);
	free(dba);
}

int getDBA_Data(DynamicByteArray *dba, int pos)
{
	if(pos>=dba->size)
	{
		printf("Error: wrong position of DBA.\n");
		exit(0);
	}
	return dba->array[pos];
}

inline void addDBA_Data(DynamicByteArray *dba, char value)
{
	if(dba->size==dba->capacity)
	{
		dba->capacity = dba->capacity << 1;
		dba->array = (char *)realloc(dba->array, dba->capacity*sizeof(char));
	}
	dba->array[dba->size] = value;
	dba->size ++;
}

void memcpyDBA_Data(DynamicByteArray *dba, char* data, int length)
{
	if(dba->size + length > dba->capacity)
	{
		dba->capacity = dba->size + length;
		dba->array = (char *)realloc(dba->array, dba->capacity*sizeof(char));
	}
	memcpy(&(dba->array[dba->size]), data, length);
	dba->size += length;
}
