#ifndef UINT_H
#define UINT_H

#include <stdbool.h>
#include <stdint.h>

typedef struct u512 {
    uint64_t c[8];
} u512;

extern const u512 u512_1;

void u512_set(u512 *x, uint64_t y);

bool u512_bit(u512 const *x, uint64_t k);

bool u512_add3(u512 *x, u512 const *y, u512 const *z); /* returns carry */
bool u512_sub3(u512 *x, u512 const *y, u512 const *z); /* returns borrow */

void u512_mul3_64(u512 *x, u512 const *y, uint64_t z);

#endif
