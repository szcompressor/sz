/**
 *  @file ExpSegment.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief 
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sz.h"
#include "ExpSegment.h"

void new_ExpSegment(ExpSegment **this, short level, int startStep)
{
	*this = (ExpSegment*)malloc(sizeof(ExpSegment));
	(*this)->fixed = 0;
	(*this)->unpredNum = 0;
	(*this)->level = level;
	(*this)->startStep = startStep;
	(*this)->length = 1;
	(*this)->medianValue_d = 0;
	(*this)->medianValue_f = 0;
	(*this)->medianValue_d0 = 0;
	(*this)->medianValue_f0 = 0;
	(*this)->offset = 0;
	(*this)->prev = NULL;
	(*this)->next = NULL;
}

void new_ExpSegment_float(ExpSegment **this, int startStep, float medianValue0, int reqLength)
{
	*this = (ExpSegment*)malloc(sizeof(ExpSegment));
	(*this)->fixed = 0;
	(*this)->level = 0;	
	(*this)->unpredNum = 0;
	(*this)->startStep = startStep;
	(*this)->medianValue_f0 = medianValue0;
	(*this)->medianValue_f = 0;
	(*this)->medianValue_d0 = 0;
	(*this)->medianValue_d = 0;	
	(*this)->offset = 0;
	(*this)->reqLength = reqLength;
	(*this)->reqBytesLength = reqLength/8;
	(*this)->resiBitsLength = reqLength%8;
	(*this)->prev = NULL;
	(*this)->next = NULL;
}

void new_ExpSegment_double(ExpSegment **this, short startStep, double medianValue0, int reqLength)
{
	*this = (ExpSegment*)malloc(sizeof(ExpSegment));
	(*this)->fixed = 0;
	(*this)->unpredNum = 0;
	(*this)->level = 0;
	(*this)->startStep = startStep;
	(*this)->medianValue_d0 = medianValue0;
	(*this)->medianValue_d = 0;
	(*this)->medianValue_f0 = 0;
	(*this)->medianValue_f = 0;	
	(*this)->offset = 0;
	(*this)->reqLength = reqLength;
	(*this)->reqBytesLength = reqLength/8;
	(*this)->resiBitsLength = reqLength%8;
	(*this)->prev = NULL;
	(*this)->next = NULL;	
}

void toLatestES(ExpSegment **es)
{
	while((*es)->next!=NULL)
		(*es) = (*es)->next;
}

void removeSelf(ExpSegment* this)
{
	if(this->startStep!=-1)
		this->prev->next = this->next;
	if(this->next!=NULL)
		this->next->prev = this->prev;
	free(this);
}

void mergeSelf_bestfit(ExpSegment* this, ExpSegment** merged, short prevLevel, short nextLevel, int lengthBound)
{
	short level = this->level;
	if(level>prevLevel && level>nextLevel) // this is a bump
	{	
		*merged = NULL; //don't fit merging requirement
		return;
	}
	if(this->fixed==1 || level==nextLevel)
	{
		*merged = this;
		return;
	}
	
	short bestLevel; 
	if(this->fixed==1 || this->prev->startStep==-1 || prevLevel <= level)
		bestLevel = nextLevel;
	else if(nextLevel <= level)
		bestLevel = prevLevel;
	else
		bestLevel = prevLevel < nextLevel ? prevLevel:nextLevel;
	
	if(this->length*(bestLevel-level) < lengthBound) //should merge_bestfit
	//if(this->length < lengthBound) //should merge_bestfit
	{
		if((prevLevel<level && level<=nextLevel) || (prevLevel > nextLevel && nextLevel >= level)) //right merge
		{
			this->level = nextLevel;
			//this.length++;
			*merged = this;				
		}
		else//left merge
		{
			this->prev->length += this->length;
			*merged = this->prev;
			removeSelf(this);
			if((*merged)->fixed==0 && (*merged)->prev->startStep!=-1)//at the moment, merged here is actually the "current" segment
			{
				mergeSelf_bestfit(*merged, merged, (*merged)->prev->level, nextLevel, lengthBound);
				if(*merged==NULL) //this is a bump
					*merged = this->prev;
			}
			if((*merged)->fixed==1)
				(*merged)->fixed = 0;
		}
	}
	else
	{
		this->fixed = 1;
		*merged = this; //no merge
	}
}

short mergeSelf_toprev(ExpSegment* this)
{
	this->prev->length += this->length;
	removeSelf(this);
	return this->prev->level;
}

short mergeSelf_tonext(ExpSegment* this)
{
	this->next->length += this->length;
	this->next->startStep = this->startStep;
	removeSelf(this);
	return this->next->level;
}

int getEndStep(ExpSegment* this)
{
	return this->startStep+this->length-1;
}

void cleanUp(ExpSegment* this, void* data, long dataLength, float errBound, int DATA_TYPE)
{
	int endStep = 0;
	if(DATA_TYPE==SZ_DOUBLE)
	{
		double* values = (double*) data;
		ExpSegment* curES = this->next;//this is header
		if(curES==NULL)
			return;
		ExpSegment* nextES = curES->next;
		while(nextES!=NULL)
		{
			short nextLevel = nextES->level;
			if(curES->level==nextLevel)
			{
				mergeSelf_toprev(nextES);
				//curES is the same.
			}
			else 
			{
				int i = nextES->startStep;
				double max = values[i], min = max;

				i++;
				for(;i<getEndStep(nextES);i++)
				{
					double v = values[i];
					if(max < v)
						max = v;
					if(min > v)
						min = v;
				}
				if(max - min < errBound)
				{
					mergeSelf_toprev(nextES);
				}
				else
				{
					curES = nextES;
				}
			}
			
			nextES = curES->next;
		}
	}
	else if(DATA_TYPE==SZ_FLOAT)
	{
		float* values = (float*) data;
		ExpSegment* curES = this->next;//this is header
		if(curES==NULL)
			return;
		ExpSegment* nextES = curES->next;
		while(nextES!=NULL)
		{
			short nextLevel = nextES->level;
			if(curES->level==nextLevel)
			{
				mergeSelf_toprev(nextES);
				//curES is the same.
			}
			else 
			{
				int i = nextES->startStep;
				float max = values[i], min = max;

				i++;
				for(;i<getEndStep(nextES);i++)
				{
					float v = values[i];
					if(max < v)
						max = v;
					if(min > v)
						min = v;
				}
				if(max - min < errBound)
				{
					mergeSelf_toprev(nextES);
				}
				else
				{
					curES = nextES;
				}
			}
			
			nextES = curES->next;
		}
	}
	else
	{
		printf("Error: wrong type of data\n");
		exit(0);
	}
}

void updateEndStep(ExpSegment* this)
{
	ExpSegment* curES = this->next;//this is the first one
	while(curES!=NULL)
	{
		curES->endStep = getEndStep(curES);
		curES = curES->next;
	}
}
