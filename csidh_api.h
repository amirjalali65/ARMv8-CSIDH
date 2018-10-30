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

////////////////////////// Main API //////////////////////////////////////////
/*
The validate fucntion gets a pointer to the CSIDH public key and return a boolean value
corresponding to the validation of the public key. If the input public key is valid, the
return value is true, otherwise it is false. Since the validation algorithm does not contain
any secret value, the implementation is exactly the one which is proposed by Castryck et al.
and it is totally non-constant time.
*/
bool csidh_validate(const public_key_t in);

/*
The keypair function generates CSIDH public key and private key for each party. The constant-time
implementation of this function is included in this library. The generated public key is computed 
by performing an action on the starting curve. The action operation is implemented both in constant-time. 
*/
void csidh_keypair(private_key_t priv, public_key_t pub);

/*
The shared secrete generation is basically a wrapper around the action operation which gets 
Alice's public-key and Bob's private key, or Alice's private key and Bob's public key to generate
the shared secret using CSIDH for each party. This function is implemented in constant time.
*/
void csidh_sharedsecret(const public_key_t in, const private_key_t priv, shared_secret_t out);

#endif