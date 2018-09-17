#ifndef FP_H
#define FP_H

#include "u512.h"

/* fp is in the Montgomery domain, so interpreting that
   as an integer should never make sense.
   enable compiler warnings when mixing up u512 and fp. */
typedef struct fp {
    u512 x;
} fp;

extern const fp fp_0;
extern const fp fp_1;

void fp_set(fp *x, uint64_t y);
void fp_cswap(fp *x, fp *y, bool c);

void fp_enc(fp *x, u512 const *y); /* encode to Montgomery representation */
void fp_dec(u512 *x, fp const *y); /* decode from Montgomery representation */

void fp_add2(fp *x, fp const *y);
void fp_sub2(fp *x, fp const *y);
void fp_mul2(fp *x, fp const *y);

void fp_add3(fp *x, fp const *y, fp const *z);
void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_mul3(fp *x, fp const *y, fp const *z);

void fp_sq1(fp *x);
void fp_sq2(fp *x, fp const *y);
void fp_inv(fp *x);
bool fp_issquare(fp const *x);

void fp_random(fp *x);

#endif
