/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*
*   Author: Amir Jalali                     ajalali2016@fau.edu
*                       
*                       All rights reserved   
*****************************************************************************/
#ifndef ARITH_H
#define ARITH_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
//////////////////  Parameters ///////////////////////////

#define SMALL_PRIMES_COUNT  74
#define NWORDS_64   8
#define MAX_EXPONENT    5

uint64_t one_Mont[NWORDS_64];
const uint64_t smallprimes[SMALL_PRIMES_COUNT];
const uint64_t four_sqrt_p[8];
uint64_t zero[NWORDS_64];
//////////////////  Datatypes  ///////////////////////////

// Datatype for representing 512-bit integer
typedef uint64_t UINT512_t[NWORDS_64];

// Datatype for a pointer representing 512-bit finite field element
typedef uint64_t felm_t[NWORDS_64];

// Datatype for presenting projective points
typedef struct {
    felm_t X;
    felm_t Z;
} proj_point;

// Datatype for presenting projective curve coefficients
typedef struct {
    felm_t A;
    felm_t C;
} proj_coeff;

// Static defenition of a pointer to a projective point 
typedef proj_point proj_point_t[1];

// Static definition of a pointer to a projective curve coefficient
typedef proj_coeff proj_coeff_t[1];

///////////////////  Integer Arithmetic ////////////////////

bool mp_add_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

unsigned int mp_sub_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void mp_mul_u64(const uint64_t *a, const uint64_t b, uint64_t *c);

void mp_U512_set_zero(uint64_t *a);

void mp_U512_set_one(uint64_t *a);

int mp_U512_bit(const uint64_t *a, uint64_t k);

///////////////////  Field Arithmetic  /////////////////////

void fp_random_512(uint64_t *a);

void fp_add_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_sub_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_mul_mont_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_sqr_mont_512(const uint64_t *a, uint64_t *c);

void fp_inv(uint64_t *a);

bool fp_issquare(const uint64_t *a);

void to_mont(const uint64_t *in, uint64_t *out);

void from_mont(const uint64_t *in, uint64_t *out);

void fp_cpy(const uint64_t *in, uint64_t *out);

void fp_init_zero(uint64_t *a);

void fp_init_one(uint64_t *a);

void fp_print(uint64_t *a);


///////////////////  Group Arithmetic  //////////////////////
void swap_points(proj_point_t P, proj_point_t Q, const uint64_t mask);

void xDBL(proj_point_t Q, const proj_point_t A, const proj_point_t P);

void xADD(proj_point_t S, const proj_point_t P, const proj_point_t Q, const proj_point_t PQ);

void xDBLADD(proj_point_t R, proj_point_t S, const proj_point_t P, const proj_point_t Q, const proj_point_t PQ, const proj_point_t A);

void xMUL(proj_point_t Q, const proj_point_t A, proj_point_t P, const UINT512_t k);

void xISOG(proj_point_t A, proj_point_t P, const proj_point_t K, uint64_t k);


#endif