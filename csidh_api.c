#include <string.h>
#include <assert.h>
#include "csidh_api.h"
#include "rng.h"

void csidh_keypair(private_key_t priv, public_key_t pub)
{
    int i;
    public_key_t base_curve;
    fp_init_zero(base_curve->A);
    memset(&priv->exponents, 0, sizeof(priv->exponents)); 
    // TODO: Private-key generation should be constant time?! 
    for (i = 0; i < SMALL_PRIMES_COUNT; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->exponents[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= SMALL_PRIMES_COUNT)
                    break;
            }
        }
    }
    action(base_curve, priv, pub);
}

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj_point_t *P, const proj_point_t A, size_t lower, size_t upper)
{
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

    xMUL(P[mid], A, P[lower], cu);
    xMUL(P[lower], A, P[lower], cl);

    cofactor_multiples(P, A, lower, mid);
    cofactor_multiples(P, A, mid, upper);
}

/* never accepts invalid keys. */
bool validate(const public_key_t in)
{
    proj_point_t A;
    fp_cpy(in->A, A->X);
    fp_cpy(one_Mont, A->Z);

    do {
        proj_point_t P[SMALL_PRIMES_COUNT];
        fp_random_512(P[0]->X);
        fp_cpy(one_Mont, P[0]->Z);
        
        /* maximal 2-power in p+1 */
        xDBL(P[0], A, P[0]);
        xDBL(P[0], A, P[0]);

        cofactor_multiples(P, A, 0, SMALL_PRIMES_COUNT);

        UINT512_t order;
        mp_U512_set_one(order);

        for (size_t i = SMALL_PRIMES_COUNT - 1; i < SMALL_PRIMES_COUNT; --i) {

            /* we only gain information if [(p+1)/l] P is non-zero */
            if (memcmp(P[i]->Z, zero, sizeof(felm_t))) {

                UINT512_t t;
                mp_U512_set_zero(t);
                t[0] = smallprimes[i];
                xMUL(P[i], A, P[i], t);

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

/* compute x^3 + Ax^2 + x */
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

/* totally not constant-time. */
void action(const public_key_t in, const private_key_t priv, public_key_t out)
{
    UINT512_t k[2];
    mp_U512_set_zero(k[0]);
    mp_U512_set_zero(k[1]);
    k[0][0] = 4; 
    k[1][0] = 4;

    uint8_t e[2][SMALL_PRIMES_COUNT];
    int8_t t = 0;

    for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) {

        t = (int8_t) (priv->exponents[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            mp_mul_u64(k[1], smallprimes[i], k[1]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            mp_mul_u64(k[0], smallprimes[i], k[0]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            mp_mul_u64(k[0], smallprimes[i], k[0]);
            mp_mul_u64(k[1], smallprimes[i], k[1]);
        }
    }
    proj_point_t A;
    fp_cpy(in->A, A->X);
    fp_cpy(one_Mont, A->Z);
    
    UINT512_t one;
    mp_U512_set_one(one);
    
    // Definition of variable should be outside of the loop
    proj_point_t P;
    felm_t rhs;
    UINT512_t cof;

    bool done[2] = {false, false};

    do {

        assert(!memcmp(A->Z, one_Mont, sizeof(felm_t)));
        
        fp_random_512(P->X);
        fp_cpy(one_Mont, P->Z);
        
        get_mont_rhs(A->X, P->X, rhs);
        bool sign = !fp_issquare(rhs);

        if (done[sign])
            continue;

        xMUL(P, A, P, k[sign]);

        done[sign] = true;

        for (size_t i = 0; i < SMALL_PRIMES_COUNT; ++i) {

            if (e[sign][i]) {

                
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

    } while (!(done[0] && done[1]));

    fp_cpy(A->X, out->A);
}

/* includes public-key validation. */
bool csidh_sharedsecret(const public_key_t in, const private_key_t priv, shared_secret_t out)
{
    public_key_t tmp;
    bool valid = true;
    fp_init_zero(tmp->A);
    if (!validate(in)) {
        fp_random_512(out->A);
        valid = false;
    }
    action(in, priv, tmp);
    fp_cpy(tmp->A, out->A);
    return valid;
}

