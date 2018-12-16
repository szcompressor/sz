/**
 *  @file ArithmeticCoding.h
 *  @author Sheng Di
 *  @date Dec, 2018
 *  @brief Header file for the ArithmeticCoding.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ArithmeticCoding_H
#define _ArithmeticCoding_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#define ONE_FOURTH (0x4000000000000000)
#define ONE_HALF (0x8000000000000000)
#define THREE_FOURTHS (0xC000000000000000)
#define MAX_CODE (0xFFFFFFFFFFFFFFFF)

typedef struct Prob {
    size_t low;
    size_t high;
    int state;
} Prob;

typedef struct AriCoder
{
	int numOfRealStates; //the # real states menas the number of states after the optimization of # intervals
	int numOfValidStates; //the # valid states means the number of non-zero frequency cells (some states/codes actually didn't appear)
	size_t total_frequency;	
	Prob* cumulative_frequency; //used to encode data more efficiencly
} AriCoder;

void output_bit_1(int* buf);
void output_bit_0(int* buf);
int output_bit_1_plus_pending(int pending_bits);
int output_bit_0_plus_pending(int pending_bits);

AriCoder *createAriCoder(int numOfStates, int *s, size_t length);
void ari_init(AriCoder *ariCoder, int *s, size_t length);
unsigned int pad_ariCoder(AriCoder* ariCoder, unsigned char** out);
void unpad_ariCoder(AriCoder** ariCoder, unsigned char* bytes);

inline unsigned char get_bit(unsigned char* p, int offset);

void ari_encode(AriCoder *ariCoder, int *s, size_t length, unsigned char *out, size_t *outSize);
void ari_decode(AriCoder *ariCoder, unsigned char *s, size_t s_len, size_t targetLength, int *out);

Prob* getCode(AriCoder *ariCoder, size_t scaled_value);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ArithmeticCoding_H  ----- */

