/**
 *  @file DynamicIntArray.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief Dynamic Int Array
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "DynamicIntArray.h"

void new_DIA(DynamicIntArray **dia, int cap) {
		*dia = (DynamicIntArray *)malloc(sizeof(DynamicIntArray));
        (*dia)->size = 0;
        (*dia)->capacity = cap;
        (*dia)->array = (char*)malloc(sizeof(char)*cap);
    }

void convertDIAtoInts(DynamicIntArray *dia, char **data)
{
	int size = dia->size;
	if(size>0)
		*data = (char*)malloc(size * sizeof(char));
	else
		*data = NULL;
	memcpy(*data, dia->array, size*sizeof(char));	
}

void free_DIA(DynamicIntArray *dia)
{
	free(dia->array);
	free(dia);
}

int getDIA_Data(DynamicIntArray *dia, int pos)
{
	if(pos>=dia->size)
	{
		printf("Error: wrong position of DIA.\n");
		exit(0);
	}
	return dia->array[pos];
}

inline void addDIA_Data(DynamicIntArray *dia, int value)
{
	if(dia->size==dia->capacity)
	{
		dia->capacity = dia->capacity << 1;
		dia->array = (char *)realloc(dia->array, dia->capacity*sizeof(char));
	}
	dia->array[dia->size] = (char)value;
	dia->size ++;
}
