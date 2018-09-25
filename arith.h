/*****************************************************************************
 *      Constant-time Arithmetic Library for CSIDH 
 * 
 *      Author: Amir Jalali             ajalali2016@fau.edu
 * ***************************************************************************/
#ifndef ARITH_H
#define ARITH_H

#include <stdbool.h>
#include <stdint.h>

//////////////////  Datatypes  ///////////////////////////

// Datatype for representing 512-bit integer
typedef uint64_t UINT512_t[8];

// Datatype for a pointer representing 512-bit finite field element
typedef uint64_t felm_t[8];

// Datatype for presenting projective points
typedef struct {
    felm_t X;
    felm_t Z;
} proj_point;

// Static defenition of a pointer to a projective point 
typedef proj_point proj_point_t[1];
///////////////////  Integer Arithmetic ////////////////////

bool mp_add_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

bool mp_sub_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void mp_mul_u64(const uint64_t *a, const uint64_t b, uint64_t *c);

///////////////////  Field Arithmetic  /////////////////////

void fp_add_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_sub_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_mul_mont_512(const uint64_t *a, const uint64_t *b, uint64_t *c);

void fp_sqr_mont_512(const uint64_t *a, uint64_t *c);

void fp_inv(uint64_t *a);

void to_mont(const uint64_t *in, uint64_t *out);

void from_mont(const uint64_t *in, uint64_t *out);

void fp_cpy(const uint64_t *in, uint64_t *out);
/*
///////////////////  Group Arithmetic  //////////////////////

void xDBL(const proj_point_t P, const felm_t A24, const felm_t C24, proj_point_t Q);

void xADD(const proj_point_t P, proj_point_t Q, const proj_point_t xPQ);

void xDBLADD(proj_point_t P, proj_point_t Q, const felm_t xPQ, const felm_t A24);

void mont_ladder();

void xISOG(proj_point_t Q, const proj_point_t R, uint64_t k, felm_t A, felm_t C);
*/

#endif