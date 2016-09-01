/**
 *  @file Huffman.h
 *  @author Sheng Di
 *  @date Aug., 2016
 *  @brief Header file for the exponential segment constructor.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _Huffman_H
#define _Huffman_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct node_t {
	struct node_t *left, *right;
	int freq;
	char c;
} *node;
 
//for multi-thread version (in the future), these global variables must be carefully handled
//e.g., share by passing addresses instead.
struct node_t pool[128];
node qqq[127], *qq;
int n_nodes, qend; //n_nodes is for compression
char bufByte;
unsigned char code[128];
unsigned char cout[128];
char n_inode; //n_inode is for decompression

node huff_new_node(int freq, char c, node a, node b);
node huff_new_node2(char c);
void huff_qinsert(node n);
node huff_qremove();
void Huff_build_code(node n, int len, unsigned int out);
void huff_init(const char *s);
int huff_encode(const char *s, char *out);
void huff_decode(char *s, int targetLength, node t, char *out);
void pad_tree(char*L, char* R, char* C, char i, node root);
int convert_HuffTree_to_bytes_8states(char* out);
int convert_HuffTree_to_bytes_16states(char* out);
int convert_HuffTree_to_bytes_32states(char* out);
int convert_HuffTree_and_Encode_8states(char* type, int dataLength, char** out);
int convert_HuffTree_and_Encode_16states(char* type, int dataLength, char** out);
int convert_HuffTree_and_Encode_32states(char* type, int dataLength, char** out);
void unpad_tree(char* L, char* R, char* C, char i, node root);
node reconstruct_HuffTree_from_bytes_8states(char* bytes);
node reconstruct_HuffTree_from_bytes_16states(char* bytes);
node reconstruct_HuffTree_from_bytes_32states(char* bytes);
void reconstruct_HuffTree_and_Decode_8states(char *bytes, int targetDataLength, char** out);
void reconstruct_HuffTree_and_Decode_16states(char *bytes, int targetDataLength, char** out);
void reconstruct_HuffTree_and_Decode_32states(char *bytes, int targetDataLength, char** out);

#ifdef __cplusplus
}
#endif

#endif
