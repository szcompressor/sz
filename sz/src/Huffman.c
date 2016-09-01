/**
 *  @file Huffman.c
 *  @author Sheng Di
 *  @date Aug., 2015
 *  @brief Huffman Encoding, Compression and Decompression functions
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Huffman.h"

node huff_new_node(int freq, char c, node a, node b)
{
	node n = pool + n_nodes++;
	if (freq) n->c = c, n->freq = freq;
	else {
		n->left = a, n->right = b;
		n->freq = a->freq + b->freq;
	}
	return n;
}
 
node huff_new_node2(char c)
{
	pool[n_nodes].c = c;
	return pool+n_nodes++;
} 
 
/* priority queue */
void huff_qinsert(node n)
{
	int j, i = qend++;
	while ((j = i / 2)) {
		if (qq[j]->freq <= n->freq) break;
		qq[i] = qq[j], i = j;
	}
	qq[i] = n;
}
 
node huff_qremove()
{
	int i, l;
	node n = qq[i = 1];
 
	if (qend < 2) return 0;
	qend--;
	while ((l = i * 2) < qend) {
		if (l + 1 < qend && qq[l + 1]->freq < qq[l]->freq) l++;
		qq[i] = qq[l], i = l;
	}
	qq[i] = qq[qend];
	return n;
}
 
/* walk the tree and put 0s and 1s */
void huff_build_code(node n, int len, unsigned int out)
{
	if (n->c) {
		code[n->c] = out >> 1;
		cout[n->c] = len;
		return;
	}

	out = (out | 0) << 1;
	huff_build_code(n->left, len + 1, out);
	out = (out>>1 | 1) << 1;
	huff_build_code(n->right, len + 1, out);
}
 
void huff_init(const char *s)
{
	int i, freq[128] = {0};
	unsigned int c = 0;
 
	while(*s) freq[(int)*s++]++;
	for (i = 0; i < 128; i++)
		if (freq[i]) huff_qinsert(huff_new_node(freq[i], i, 0, 0));
	while (qend > 2) 
		huff_qinsert(huff_new_node(0, 0, huff_qremove(), huff_qremove()));
	huff_build_code(qq[1], 0, c);
}
 
/**
 * Huffman encoding
 * @par s: input char array
 * @par out: output byte stream
 * 
 * @return the size of output byte stream (in byte) 
 * 
 * */
int huff_encode(const char *s, char *out)
{
	int i=0, r=8; //i:index of out, r: residual index
	char c; //length of code
	unsigned char codeByte;
	while (*s) 
	{
		codeByte = code[*s];
		c = cout[*s];
		if(r>=c)
		{
			out[i] |= (codeByte << (r-c));
			r -= c;
		}
		else
		{
			out[i] |= (codeByte >> (c-r));
			i++;
			out[i] |= (codeByte << (r-c+8));
			r = r - c + 8;
		}
		s++;
	}
	return i+1; //length
}
 
/**
 * Huffman decoding
 * @par s: input encoded byte stream
 * @par targetLength: the target size of decoded char array
 * @par out: the decoded char array
 * 
 * */ 
void huff_decode(char *s, int targetLength, node t, char *out)
{
	int i = 0, r, byteIndex= 0, count=0;
	node n = t;
	char byte;
	for(i=0;count<targetLength;i++)
	{
		byteIndex = i>>3; //i/8
		r = i%8;
		if(((s[byteIndex] >> (7-r)) & 0x01) == 0)
			n = n->left;
		else
			n = n->right;
		if (n->c) {
			//putchar(n->c); 
			out[count] = n->c - 48;
			n = t; 
			count++;
		}
	}
//	putchar('\n');
	if (t != n) printf("garbage input\n");
	return;
}
	 
void pad_tree(char*L, char* R, char* C, char i, node root)
{
	C[i] = root->c;
	node lroot = root->left;
	if(lroot!=0)
	{
		n_inode++;
		L[i] = n_inode;
		pad_tree(L,R,C,n_inode, lroot);
	}
	node rroot = root->right;
	if(rroot!=0)
	{
		n_inode++;
		R[i] = n_inode;
		pad_tree(L,R,C,n_inode, rroot);
	}
} 
 
int convert_HuffTree_to_bytes_8states(char* out) 
{
	char L[15] = {0}, R[15] = {0}, C[15] = {0};
	pad_tree(L,R,C,0,qq[1]);
	memcpy(out, L, 15);
	memcpy(out+15,R,15);
	memcpy(out+30,C,15);
	return 45;
}

int convert_HuffTree_to_bytes_16states(char* out) 
{
	char L[31] = {0}, R[31] = {0}, C[31] = {0};
	pad_tree(L,R,C,0,qq[1]);
	memcpy(out, L, 31);
	memcpy(out+31,R,31);
	memcpy(out+62,C,31);
	return 93;
}

int convert_HuffTree_to_bytes_32states(char* out) 
{
	char L[63] = {0}, R[63] = {0}, C[63] = {0};
	pad_tree(L,R,C,0,qq[1]);
	memcpy(out, L, 63);
	memcpy(out+63,R,63);
	memcpy(out+126,C,63);
	return 189;
}

int convert_HuffTree_and_Encode_8states(char* type, int dataLength, char** out)
{
	int size = 45+dataLength*3/8+1;
	*out = (char*)malloc(sizeof(char)*size);
	memset(*out, 0, size);
	convert_HuffTree_to_bytes_8states(*out);
	int encodeSize = huff_encode(type, (*out)+45); //treeSize==45
	return 45+encodeSize;
}

int convert_HuffTree_and_Encode_16states(char* type, int dataLength, char** out)
{
	int size = 93+dataLength*4/8+1;
	*out = (char*)malloc(sizeof(char)*size);
	memset(*out, 0, size);
	convert_HuffTree_to_bytes_16states(*out);
	int encodeSize = huff_encode(type, (*out)+93); //treeSize==45
	return 93+encodeSize;
}

int convert_HuffTree_and_Encode_32states(char* type, int dataLength, char** out)
{
	int size = 189+dataLength*5/8+1;
	*out = (char*)malloc(sizeof(char)*size);
	memset(*out, 0, size);
	convert_HuffTree_to_bytes_32states(*out);
	int encodeSize = huff_encode(type, (*out)+189); //treeSize==45
	return 189+encodeSize;
}
 
void unpad_tree(char* L, char* R, char* C, char i, node root)
{
	//root->c = C[i];
	if(root->c==0)
	{
		char l, r;
		l = L[i];
		if(l!=0)
		{
			node lroot = huff_new_node2(C[l]);
			root->left = lroot;
			unpad_tree(L,R,C,l,lroot);
		}
		r = R[i];
		if(r!=0)
		{
			node rroot = huff_new_node2(C[r]);
			root->right = rroot;
			unpad_tree(L,R,C,r,rroot);
		}
	}
}
 
/**
 * The length of the bytes must be 45.
 * @par bytes: [L(15bytes)][R(15bytes)][C(15bytes)]
 * */
node reconstruct_HuffTree_from_bytes_8states(char* bytes)
{
	char L[15], R[15], C[15];
	memcpy(L, bytes, 15);
	memcpy(R, bytes+15, 15);
	memcpy(C, bytes+30, 15);
	
	node root = huff_new_node2(0);
	unpad_tree(L,R,C,0,root);
	return root;
}

node reconstruct_HuffTree_from_bytes_16states(char* bytes)
{
	char L[31], R[31], C[31];
	memcpy(L, bytes, 31);
	memcpy(R, bytes+31, 31);
	memcpy(C, bytes+62, 31);
	
	node root = huff_new_node2(0);
	unpad_tree(L,R,C,0,root);
	return root;
}

node reconstruct_HuffTree_from_bytes_32states(char* bytes)
{
	char L[63], R[63], C[63];
	memcpy(L, bytes, 63);
	memcpy(R, bytes+63, 63);
	memcpy(C, bytes+126, 63);
	
	node root = huff_new_node2(0);
	unpad_tree(L,R,C,0,root);
	return root;
}

void reconstruct_HuffTree_and_Decode_8states(char *bytes, int targetDataLength, char** out)
{
	n_nodes = 0;
	node root = reconstruct_HuffTree_from_bytes_8states(bytes);
	*out = (char*)malloc(sizeof(char)*targetDataLength);
	huff_decode(bytes+45, targetDataLength, root, *out);
}

void reconstruct_HuffTree_and_Decode_16states(char *bytes, int targetDataLength, char** out)
{
	n_nodes = 0;
	node root = reconstruct_HuffTree_from_bytes_16states(bytes);
	*out = (char*)malloc(sizeof(char)*targetDataLength);
	huff_decode(bytes+93, targetDataLength, root, *out);
}

void reconstruct_HuffTree_and_Decode_32states(char *bytes, int targetDataLength, char** out)
{
	n_nodes = 0;
	node root = reconstruct_HuffTree_from_bytes_32states(bytes);
	*out = (char*)malloc(sizeof(char)*targetDataLength);
	huff_decode(bytes+189, targetDataLength, root, *out);
}
