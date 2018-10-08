/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*
*   Author: Modified Amir Jalali                     ajalali2016@fau.edu
*                       
*                       All rights reserved   
*****************************************************************************/
#ifndef CSIDH_API_H
#define CSIDH_API_H

#include "arith.h"

typedef struct private_key {
    int8_t exponents[(SMALL_PRIMES_COUNT + 1) / 2];
} private_key;

typedef struct public_key {
    felm_t A; 
} public_key;

typedef struct shared_secret {
    felm_t A;
} shared_secret;

typedef private_key private_key_t[1];

typedef public_key public_key_t[1];

typedef shared_secret shared_secret_t[1];

bool csidh_validate(const public_key_t in);

void action(const public_key_t in, const private_key_t priv, public_key_t out);

void csidh_keypair(private_key_t priv, public_key_t pub);

void csidh_sharedsecret(const public_key_t in, const private_key_t priv, shared_secret_t out);


#endif