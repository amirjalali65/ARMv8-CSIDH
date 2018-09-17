
#include <string.h>
#include <assert.h>

#include "csidh.h"
#include "rng.h"

/* specific to p, should perhaps be somewhere else */
const unsigned primes[num_primes] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
    317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
};

const u512 four_sqrt_p = {{
    0x85e2579c786882cf, 0x4e3433657e18da95, 0x850ae5507965a0b3, 0xa15bc4e676475964,
}};


const public_key base = {0}; /* A = 0 */

void csidh_private(private_key *priv)
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < num_primes; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= max_exponent && buf[j] >= -max_exponent) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= num_primes)
                    break;
            }
        }
    }
}

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj *P, const proj *A, size_t lower, size_t upper)
{
    assert(lower < upper);

    if (upper - lower == 1)
        return;

    size_t mid = lower + (upper - lower + 1) / 2;

    u512 cl = u512_1, cu = u512_1;
    for (size_t i = lower; i < mid; ++i)
        u512_mul3_64(&cu, &cu, primes[i]);
    for (size_t i = mid; i < upper; ++i)
        u512_mul3_64(&cl, &cl, primes[i]);

    xMUL(&P[mid], A, &P[lower], &cu);
    xMUL(&P[lower], A, &P[lower], &cl);

    cofactor_multiples(P, A, lower, mid);
    cofactor_multiples(P, A, mid, upper);
}

/* never accepts invalid keys. */
bool validate(public_key const *in)
{
    const proj A = {in->A, fp_1};

    do {

        proj P[num_primes];
        fp_random(&P->x);
        P->z = fp_1;

        /* maximal 2-power in p+1 */
        xDBL(P, &A, P);
        xDBL(P, &A, P);

        cofactor_multiples(P, &A, 0, num_primes);

        u512 order = u512_1;

        for (size_t i = num_primes - 1; i < num_primes; --i) {

            /* we only gain information if [(p+1)/l] P is non-zero */
            if (memcmp(&P[i].z, &fp_0, sizeof(fp))) {

                u512 tmp;
                u512_set(&tmp, primes[i]);
                xMUL(&P[i], &A, &P[i], &tmp);

                if (memcmp(&P[i].z, &fp_0, sizeof(fp)))
                    /* P does not have order dividing p+1. */
                    return false;

                u512_mul3_64(&order, &order, primes[i]);

                if (u512_sub3(&tmp, &four_sqrt_p, &order)) /* returns borrow */
                    /* order > 4 sqrt(p), hence definitely supersingular */
                    return true;
            }
        }

    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}

/* totally not constant-time. */
void action(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */

    uint8_t e[2][num_primes];

    for (size_t i = 0; i < num_primes; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = 0; i < num_primes; ++i) {

            if (e[sign][i]) {

                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j)
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
}

/* includes public-key validation. */
bool csidh(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    action(out, in, priv);
    return true;
}

