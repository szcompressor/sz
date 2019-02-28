//
// Created by borelset on 2018/12/24.
//

#ifndef _MULTILEVELCACHETABLEWIDEINTERVAL_H
#define _MULTILEVELCACHETABLEWIDEINTERVAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <memory.h>
#include <stdlib.h>
#include "stdio.h"

struct SubLevelTableWideInterval{
    uint64_t baseIndex;
    uint64_t topIndex;
    uint16_t* table;
    uint16_t expoIndex;
} SubLevelTableWideInterval;

struct TopLevelTableWideInterval{
    uint16_t bits;
    uint16_t baseIndex;
    uint16_t topIndex;
    struct SubLevelTableWideInterval* subTables;
    double bottomBoundary;
    double topBoundary;
} TopLevelTableWideInterval;

uint16_t MLCTWI_GetExpoIndex(double value);
uint16_t MLCTWI_GetRequiredBits(double precision);
uint64_t MLCTWI_GetMantiIndex(double value, int bits);

double MLTCWI_RebuildDouble(uint16_t expo, uint64_t manti, int bits);
void MultiLevelCacheTableWideIntervalBuild(struct TopLevelTableWideInterval* topTable, double* precisionTable, int count, double precision, int plus_bits);
uint32_t MultiLevelCacheTableWideIntervalGetIndex(double value, struct TopLevelTableWideInterval* topLevelTable);
void MultiLevelCacheTableWideIntervalFree(struct TopLevelTableWideInterval* table);

#ifdef __cplusplus
}
#endif

#endif //_MULTILEVELCACHETABLEWIDEINTERVAL_H
