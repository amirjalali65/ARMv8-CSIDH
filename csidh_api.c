/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*   The changes made to the original POC code of CSIDH by Castryck et al.
*   Author: Created by Amir Jalali                     ajalali2016@fau.edu
*                       
*                       All rights reserved   
*****************************************************************************/
#include <string.h>
#include <assert.h>
#include "csidh_api.h"
#include "rng.h"


/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj_point_t *P, const proj_point_t A, size_t lower, size_t upper)
{
    // Since this function is only called by csidh_validate, it does not need to 
    // be constant-time from the security point of view 
    assert(lower < upper);  

    if (upper - lower == 1)
        return; 

    size_t mid = lower + (upper - lower + 1) / 2;

    UINT512_t cl, cu;
    mp_U512_set_one(cl);
    mp_U512_set_one(cu);
    for (size_t i = lower; i < mid; ++i)
        mp_mul_u64(cu, smallprimes[i], cu);
    for (size_t i = mid; i < upper; ++i)
        mp_mul_u64(cl, smallprimes[i], cl);

    xMUL_non_const(P[mid], A, P[lower], cu);
    xMUL_non_const(P[lower], A, P[lower], cl);

    cofactor_multiples(P, A, lower, mid);
    cofactor_multiples(P, A, mid, upper);
}

bool csidh_validate(const public_key_t in)
{
    // Since validation does not any secret information, the non-constant time
    // implementation does not seem to expose any vulnerability to the scheme

    proj_point_t A, P[SMALL_PRIMES_COUNT];
    fp_cpy(in->A, A->X);
    fp_cpy(one_Mont, A->Z);

    UINT512_t order, t;
    do {
      
        fp_random_512(P[0]->X);
        fp_cpy(one_Mont, P[0]->Z);
        
        /* maximal 2-power in p+1 */
        xDBL(P[0], A, P[0]);
        xDBL(P[0], A, P[0]);

        cofactor_multiples(P, A, 0, SMALL_PRIMES_COUNT);

    
        mp_U512_set_one(order);

        for (size_t i = SMALL_PRIMES_COUNT - 1; i < SMALL_PRIMES_COUNT; --i) {

            /* we only gain information if [(p+1)/l] P is non-zero */
            if (memcmp(P[i]->Z, zero, sizeof(felm_t))) {

                mp_U512_set_zero(t);
                t[0] = smallprimes[i];
                xMUL_non_const(P[i], A, P[i], t);

                if (memcmp(P[i]->Z, zero, sizeof(felm_t)))
                    /* P does not have order dividing p+1. */
                    return false;

                mp_mul_u64(order, smallprimes[i], order);

                if (mp_sub_512(four_sqrt_p, order, t))
                    /* order > 4 sqrt(p), hence definitely supersingular */
                    return true;
            }
        }

    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}

static void get_mont_rhs(const felm_t A, const felm_t x, felm_t rhs)
{
    felm_t t;
    fp_cpy(x, rhs);
    fp_sqr_mont_512(rhs, rhs);
    fp_mul_mont_512(A, x, t);
    fp_add_512(t, rhs, rhs);
    fp_add_512(one_Mont, rhs, rhs);
    fp_mul_mont_512(rhs, x, rhs);
}

// non-constant and constant-time implementation of action
static void action(const public_key_t in, const private_key_t priv, public_key_t out)
{
    UINT512_t k[2] = {{0}};
    k[0][0] = 4; 
    k[1][0] = 4;

    uint8_t e[2][SMALL_PRIMES_COUNT];
    int8_t t = 0;
    
#ifdef _CONSTANT_ 
    uint8_t t_sign;
    bool is_nonzero;

    for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) 
    {
        t = (int8_t) (priv->exponents[i / 2] << i % 2 * 4) >> 4;
        t_sign = ((t & 0x80) >> 7 | !t);
        is_nonzero = (bool)t;

        e[t_sign][i] = t - (2 * t_sign) * t;
        e[!t_sign][i] = 0;
        mp_mul_u64(k[!t_sign], smallprimes[i], k[!t_sign]);
        mp_mul_u64(k[!is_nonzero], (smallprimes[i] - ((is_nonzero)*(smallprimes[i]-1))), k[!is_nonzero]);
    }
#else
    for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) 
    {
        t = (int8_t) (priv->exponents[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) 
        {
            e[0][i] = t;
            e[1][i] = 0;
            mp_mul_u64(k[1], smallprimes[i], k[1]);
        }
        else if (t < 0)
        {
            e[1][i] = -t;
            e[0][i] = 0;
            mp_mul_u64(k[0], smallprimes[i], k[0]);
        }
        else 
        {
            e[0][i] = 0;
            e[1][i] = 0;
            mp_mul_u64(k[0], smallprimes[i], k[0]);
            mp_mul_u64(k[1], smallprimes[i], k[1]);
        }
    }
#endif
    proj_point_t A, P; UINT512_t one, cof; felm_t rhs;
    fp_cpy(in->A, A->X);
    fp_cpy(one_Mont, A->Z);

    mp_U512_set_one(one);

    bool  done[2] = {false, false};
#ifdef _CONSTANT_
    int count;
    bool donemask, sign, mask, esign_mask;
    proj_point_t bigA, AA, PP, K;
    unsigned int z_is_zero;
    uint64_t correction;

    for(count = 0; count <= UPPER_BOUND; count++) 
    {
        fp_cpy(A->X, bigA->X);
        fp_random_512(P->X);
        fp_cpy(one_Mont, P->Z);
        
        get_mont_rhs(A->X, P->X, rhs);
        sign = !fp_issquare(rhs);

        xMUL(P, A, P, k[sign]);

        done[sign] = true;
    
    
        for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) 
        {
            fp_cpy(A->X, AA->X);
            fp_cpy(A->Z, AA->Z);
            fp_cpy(P->X, PP->X);
            fp_cpy(P->Z, PP->Z);

            esign_mask = e[sign][i];
            mp_U512_set_one(cof);
            for (size_t j = i + 1; j < SMALL_PRIMES_COUNT; ++j)
            {
                mask = !e[sign][j];
                correction = mask * (smallprimes[j] - 1);
                mp_mul_u64(cof, (smallprimes[j] - correction), cof);
            }
            xMUL(K, A, P, cof);

            z_is_zero = !memcmp(K->Z, zero, sizeof(felm_t));

            xISOG(A, P, K, smallprimes[i]);
            cswap(A, AA, (0 - (uint64_t)(z_is_zero | !esign_mask)));
            cswap(P, PP, (0 - (uint64_t)(z_is_zero | !esign_mask)));

            mask = (--e[sign][i] | (bool)z_is_zero);
            mask = (mask | !esign_mask);
            e[sign][i] += z_is_zero;
            correction = mask * (smallprimes[i] - 1);
                
            mp_mul_u64(k[sign], (smallprimes[i] - correction), k[sign]);
            done[sign] &= !e[sign][i];
        }

        fp_inv(A->Z);
        fp_mul_mont_512(A->X, A->Z, A->X);
        fp_cpy(one_Mont, A->Z);     
        donemask ^= donemask;   
        cswap(A, bigA, (0 - (uint64_t)donemask));
        donemask = (done[0] & done[1]);
    } 
#else
        do
        {   
            fp_random_512(P->X);
            fp_cpy(one_Mont, P->Z);
            
            get_mont_rhs(A->X, P->X, rhs);
            bool sign = !fp_issquare(rhs);
            
            if (done[sign])
                continue;
            
            xMUL(P, A, P, k[sign]);

            done[sign] = true;
            for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) 
            {
                if (e[sign][i]) 
                {

                    
                    mp_U512_set_one(cof);
                    for (size_t j = i + 1; j < SMALL_PRIMES_COUNT; ++j)
                        if (e[sign][j])
                            mp_mul_u64(cof, smallprimes[j], cof);
                    proj_point_t K;
                    xMUL(K, A, P, cof);

                    if (memcmp(K->Z, zero, sizeof(felm_t))) {

                        xISOG(A, P, K, smallprimes[i]);

                        if (!--e[sign][i])
                            mp_mul_u64(k[sign], smallprimes[i], k[sign]);

                    }

                }
                done[sign] &= !e[sign][i];
            }
            fp_inv(A->Z);
            fp_mul_mont_512(A->X, A->Z, A->X);
            fp_cpy(one_Mont, A->Z);     
        }
        while(!(done[0] && done[1]));
#endif
    fp_cpy(A->X, out->A);
}

void csidh_keypair(private_key_t priv, public_key_t pub)
{
    int i, j;
    public_key_t base_curve;

    fp_init_zero(base_curve->A);
    memset(&priv->exponents, 0, sizeof(priv->exponents)); 

#ifdef _CONSTANT_
    for (i = 0; i < SMALL_PRIMES_COUNT; i++) 
    {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (j = 0; j < sizeof(buf); ++j) 
        {
            uint8_t compare_mask, full_mask;
            int8_t compare_lower, compare_upper, new_val, tmp;

            compare_upper = (MAX_EXPONENT + 1) - buf[j];
            compare_mask = !((compare_upper & 0x80) >> 7 | !compare_upper);
            compare_lower = buf[j] - (-MAX_EXPONENT - 1);
            compare_mask &= !((compare_lower & 0x80) >> 7 | !compare_lower);
            full_mask = 0 - compare_mask;

            new_val = priv->exponents[i/2] | (buf[j] & 0xf) << i % 2 * 4;     
            tmp = full_mask & (new_val ^ priv->exponents[i/2]);
            priv->exponents[i/2] = tmp ^ priv->exponents[i/2];
        }
    }
#else
    for (i = 0; i < SMALL_PRIMES_COUNT;) 
    {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (j = 0; j < sizeof(buf); ++j) 
        {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->exponents[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= SMALL_PRIMES_COUNT)
                    break;
            }
        }
    }
#endif
    // Generate Public-key
    action(base_curve, priv, pub);
}

void csidh_sharedsecret(const public_key_t in, const private_key_t priv, shared_secret_t out)
{
    public_key_t tmp;

    action(in, priv, tmp);
    fp_cpy(tmp->A, out->A);
}

