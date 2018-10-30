/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*
*   Author: Modified by Amir Jalali                     ajalali2016@fau.edu
*                       
*                       All rights reserved   
*****************************************************************************/

#include "arith.h"
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>



///////////////////  Field Arithmetic  /////////////////////

uint64_t prime511[NWORDS_64] = { 0x1b81b90533c6c87b, 0xc2721bf457aca835,
                         0x516730cc1f0b4f25, 0xa7aac6c567f35507,
                         0x5afbfcc69322c9cd, 0xb42d083aedc88c42,
                         0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf };

uint64_t r2_Mont[NWORDS_64] = { 0x36905b572ffc1724, 0x67086f4525f1f27d,
                        0x4faf3fbfd22370ca, 0x192ea214bcc584b1,
                        0x5dae03ee2f5de3d0, 0x1e9248731776b371,
                        0xad5f166e20e4f52d, 0x4ed759aea6f3917e };

uint64_t one_Mont[NWORDS_64] = { 0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 
                         0x5d319e67c1e961b4, 0xb0aa7275301955f1,
                         0x4a080672d9ba6c64, 0x97a5ef8a246ee77b,
                         0x06ea9e5d4383676a, 0x3496e2e117e0ec80 };

uint64_t zero[NWORDS_64] = {0};

const uint64_t smallprimes[SMALL_PRIMES_COUNT] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
    317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
};

const uint64_t four_sqrt_p[8] = { 0x85e2579c786882cf, 0x4e3433657e18da95,
                                  0x850ae5507965a0b3, 0xa15bc4e676475964,
                                  0x0000000000000000, 0x0000000000000000,
                                  0x0000000000000000, 0x0000000000000000};

void mp_U512_set_zero(uint64_t *a)
{
    int i;
    for(i = 0; i < NWORDS_64; i++)
    {
        a[i] = 0;
    }
}

void mp_U512_set_one(uint64_t *a)
{
    int i;
    for(i = 0; i < NWORDS_64; i++)
    {
        a[i] = 0;
    }
    a[0] = 1;
}

void fp_random_512(uint64_t *a)
{
    static int fd = -1;
    int n, i;
    if (fd < 0 && 0 > (fd = open("/dev/urandom", O_RDONLY)))
        exit(1);
    for (i = 0; i < NWORDS_64 * 8; i += n)
        if (0 >= (n = read(fd, (char *) a + i, (NWORDS_64 * 8) - i)))
            exit(2);
    a[7] &= 0x3FFFFFFFFFFFFFFF;
}


void fp_sqr_mont_512(const uint64_t *a, uint64_t *c){
    fp_mul_mont_512(a, a, c);
}

void fp_cpy(const uint64_t *a, uint64_t *c)
{
    int i;
    for(i = 0; i < 8; i++)
        c[i] = a[i];
}

void fp_init_zero(uint64_t *a)
{
    int i;
    for(i = 0; i < NWORDS_64; i++)
        a[i] = 0;
}

void fp_init_one(uint64_t *a)
{
    int i;
    for(i = 0 ; i < NWORDS_64; i++)
        a[i] = 0;
    a[0] = 1;
}

void fp_inv(uint64_t *a)
{
    // Field inversion using addition chain
    felm_t tmp[28], t;
    int i;

    // Pre-computed Table
    // tmp0->a3     tmp1->a5    tmp2->a11   tmp3->a13   tmp4->a15
    // tmp5->a17    tmp6->a19   tmp7->a21   tmp8->a25   tmp9->a27
    // tmp10->a29   tmp11->a31  tmp12->a33  tmp13->a35  tmp14->a37
    // tmp15->a39   tmp16->a41  tmp17->a43  tmp18->a45  tmp19->a47
    // tmp20->a49   tmp21->a51  tmp22->a53  tmp23->a55  tmp24->a57
    // tmp25->a59   tmp26->a61  tmp27->a63    
    fp_sqr_mont_512(a, t);
    fp_mul_mont_512(a, t, tmp[0]);
    fp_mul_mont_512(t, tmp[0], tmp[1]);
    fp_sqr_mont_512(tmp[1], tmp[2]);
    fp_mul_mont_512(tmp[2], a, tmp[2]);
    for(i = 2; i <= 7; i++ )
        fp_mul_mont_512(t, tmp[i], tmp[i+1]);
    fp_mul_mont_512(tmp[8], t, tmp[8]);
    for(i = 8; i <= 26; i++)
        fp_mul_mont_512(t, tmp[i], tmp[i+1]);

    fp_cpy(a, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[18], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 3; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[1], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[11], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[27], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[27], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[20], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[2], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[10], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[25], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[24], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[18], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[11], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[27], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[21], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[8], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[24], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[18], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[9], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[27], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[7], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[10], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[21], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[21], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[11], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[18], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[0], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[15], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[19], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[22], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[22], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[23], t, t);
    for(i = 0; i < 12; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[23], t, t);
    for(i = 0; i < 3; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[15], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[8], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 3; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    fp_cpy(t, a);
}

bool fp_issquare(const uint64_t *a)
{
    // Square-root check using addition chain
    felm_t tmp[27], t;
    int i;
    
    // Pre-computed Table
    // tmp0->a3     tmp1->a5    tmp2->a13   tmp3->a15   tmp4->a17
    // tmp5->a19    tmp6->a21   tmp7->a25   tmp8->a27   tmp9->a29
    // tmp10->a31   tmp11->a33  tmp12->a35  tmp13->a37  tmp14->a39
    // tmp15->a41   tmp16->a43  tmp17->a45  tmp18->a47  tmp19->a49
    // tmp20->a51   tmp21->a53  tmp22->a55  tmp23->a57  tmp24->a59
    // tmp25->a61   tmp26->a63
    fp_sqr_mont_512(a, t);
    fp_mul_mont_512(a, t, tmp[0]);
    fp_mul_mont_512(t, tmp[0], tmp[1]);
    fp_sqr_mont_512(tmp[1], tmp[2]);
    fp_mul_mont_512(tmp[2], tmp[0], tmp[2]);
    for(i = 2; i <= 5; i++ )
        fp_mul_mont_512(t, tmp[i], tmp[i+1]);
    fp_mul_mont_512(tmp[6], t, tmp[7]);
    fp_mul_mont_512(tmp[7], t, tmp[7]);
    for(i = 7; i <= 25; i++)
        fp_mul_mont_512(t, tmp[i], tmp[i+1]);
    
    fp_cpy(a, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 3; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[1], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[25], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[2], t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[10], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[2], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[19], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[9], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[24], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[23], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[11], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[10], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[20], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[2], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[7], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[23], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[25], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[8], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[26], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[6], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[9], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[4], t, t);
    for(i = 0; i < 7; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[20], t, t);
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[5], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[20], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[10], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[17], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[3], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[0], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 8; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[11], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[18], t, t);
    for(i = 0; i < 4; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[2], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[16], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[21], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[13], t, t);
    for(i = 0; i < 2; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[21], t, t);
    for(i = 0; i < 9; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[22], t, t);
    for(i = 0; i < 12; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[22], t, t);
    for(i = 0; i < 3; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(a, t, t);
    for(i = 0; i < 11; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[15], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[14], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[12], t, t);
    for(i = 0; i < 6; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[7], t, t);
    for(i = 0; i < 10; i++) fp_sqr_mont_512(t, t);
    fp_mul_mont_512(tmp[25], t, t);

    return (memcmp(t, one_Mont, sizeof(felm_t)) == 0) ? true : false;
}   

void to_mont(const uint64_t *in, uint64_t *out)
{
    fp_mul_mont_512(in, r2_Mont, out);
}

void from_mont(const uint64_t *in, uint64_t *out)
{
    uint64_t one[8] = { 0 };
    one[0] = 1;
    fp_mul_mont_512(in, one, out);
}

void fp_print(uint64_t *a)
{
    int i;
    printf("\n");
    for(i = 7; i >= 0; i--)
    {
        printf("%lx", a[i]);
    }
    printf("\n");    
}

//////////////// Group Arithmetic ////////////////////////
void xDBLADD(proj_point_t R, proj_point_t S, const proj_point_t P, const proj_point_t Q, const proj_point_t PQ, const proj_point_t A24)
{
    felm_t t0, t1, t2;

    fp_add_512(P->X, P->Z, t0);
    fp_sub_512(P->X, P->Z, t1);
    fp_sqr_mont_512(t0, R->X);
    fp_sub_512(Q->X, Q->Z, t2);
    fp_add_512(Q->X, Q->Z, S->X);
    fp_mul_mont_512(t0, t2, t0);
    fp_sqr_mont_512(t1, R->Z);
    fp_mul_mont_512(t1, S->X, t1);
    fp_sub_512(R->X, R->Z, t2);
    fp_mul_mont_512(A24->Z, R->Z, R->Z);
    fp_mul_mont_512(R->X, R->Z, R->X);
    fp_mul_mont_512(t2, A24->X, S->X);
    fp_sub_512(t0, t1, S->Z);
    fp_add_512(R->Z, S->X, R->Z);
    fp_add_512(t0, t1, S->X);
    fp_mul_mont_512(t2, R->Z, R->Z);
    fp_sqr_mont_512(S->Z, S->Z);
    fp_sqr_mont_512(S->X, S->X);
    fp_mul_mont_512(PQ->X, S->Z, S->Z);
    fp_mul_mont_512(PQ->Z, S->X, S->X);
}

void xDBL(proj_point_t Q, const proj_point_t A, const proj_point_t P)
{
    felm_t t0, t1, t2;

    fp_add_512(P->X, P->Z, t0);
    fp_sqr_mont_512(t0, t0);
    fp_sub_512(P->X, P->Z, t1);
    fp_sqr_mont_512(t1, t1);
    fp_sub_512(t0, t1, t2);
    fp_add_512(t1, t1, t1); 
    fp_add_512(t1, t1, t1);
    fp_mul_mont_512(t1, A->Z, t1);
    fp_mul_mont_512(t0, t1, Q->X);
    fp_add_512(A->Z, A->Z, t0); 
    fp_add_512(A->X, t0, t0);
    fp_mul_mont_512(t0, t2, t0);
    fp_add_512(t0, t1, t0);
    fp_mul_mont_512(t0, t2, Q->Z);
}

void xADD(proj_point_t S, const proj_point_t P, const proj_point_t Q, const proj_point_t PQ)
{
    felm_t t0, t1, t2, t3;

    fp_add_512(P->X, P->Z, t0);
    fp_sub_512(P->X, P->Z, t1);
    fp_add_512(Q->X, Q->Z, t2);
    fp_sub_512(Q->X, Q->Z, t3);
    fp_mul_mont_512(t0, t3, t0);
    fp_mul_mont_512(t1, t2, t1);
    fp_add_512(t0, t1, t2);
    fp_sub_512(t0, t1, t3);
    fp_sqr_mont_512(t2, t2);
    fp_sqr_mont_512(t3, t3);
    fp_mul_mont_512(PQ->Z, t2, S->X);
    fp_mul_mont_512(PQ->X, t3, S->Z);
}


void cswap(proj_point_t P, proj_point_t Q, const uint64_t mask)
{
  // Constant-time point swap   
  // mask = 0 then P <- P and Q <- Q
  // mask = 0xFF...FF then P <- Q and Q <- P
    uint64_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_64; i++) 
    {
        temp = mask & (P->X[i] ^ Q->X[i]);
        P->X[i] = temp ^ P->X[i]; 
        Q->X[i] = temp ^ Q->X[i]; 
        temp = mask & (P->Z[i] ^ Q->Z[i]);
        P->Z[i] = temp ^ P->Z[i]; 
        Q->Z[i] = temp ^ Q->Z[i]; 
    }
}

int mp_U512_bit(const uint64_t *a, uint64_t k)
{
    return (a[k >> 6] >> (k & (63))) & 1;
}

// Montgomery ladder implementation 
// Constant-time and non constant-time
void xMUL(proj_point_t Q, const proj_point_t A,  proj_point_t P, const UINT512_t k)
{
    proj_point_t R, tmp, A24, Pcopy;
    fp_cpy(P->X, R->X);
    fp_cpy(P->Z, R->Z);

    fp_cpy(P->X, Pcopy->X);
    fp_cpy(P->Z, Pcopy->Z);

    fp_cpy(one_Mont, Q->X);
    fp_cpy(zero, Q->Z);

    fp_add_512(A->Z, A->Z, A24->X);
    fp_add_512(A24->X, A24->X, A24->Z); // 4C
    fp_add_512(A24->X, A->X, A24->X);   // A + 2C

    int bit, nbits = 511;

#ifdef _CONSTANT_

    #ifdef _FASTLADDER_
        while(--nbits && !mp_U512_bit(k, nbits));
        xDBL(Q, A, P);
    #endif

    int i, swap, bprev = 0;
    uint64_t mask;

    for(i = nbits-1; i >= 0; i--)
    {
        bit = mp_U512_bit(k, i);
        swap = bit ^ bprev;
        bprev = bit;
        mask = 0 - (uint64_t)swap;
        cswap(Q, R, mask);
        xDBLADD(Q, R, Q, R, Pcopy, A24);
    }

    cswap(Q, R, (0 - (uint64_t)bit));

#else

    while(--nbits && !mp_U512_bit(k, nbits));

    do
    {
        bit = mp_U512_bit(k, nbits);
        if(bit)
        {
            fp_cpy(Q->X, tmp->X);fp_cpy(Q->Z, tmp->Z);
            fp_cpy(R->X, Q->X);fp_cpy(R->Z, Q->Z);
            fp_cpy(tmp->X, R->X);fp_cpy(tmp->Z, R->Z);
        }
        xDBLADD(Q, R, Q, R, Pcopy, A24);
        if(bit)
        {
            fp_cpy(Q->X, tmp->X);fp_cpy(Q->Z, tmp->Z);
            fp_cpy(R->X, Q->X);fp_cpy(R->Z, Q->Z);
            fp_cpy(tmp->X, R->X);fp_cpy(tmp->Z, R->Z);
        }
    } while (nbits--);
#endif
}

void xMUL_non_const(proj_point_t Q, const proj_point_t A,  proj_point_t P, const UINT512_t k)
{
    proj_point_t R, tmp, A24, Pcopy;
    fp_cpy(P->X, R->X);
    fp_cpy(P->Z, R->Z);

    fp_cpy(P->X, Pcopy->X);
    fp_cpy(P->Z, Pcopy->Z);

    fp_cpy(one_Mont, Q->X);
    fp_cpy(zero, Q->Z);

    fp_add_512(A->Z, A->Z, A24->X);
    fp_add_512(A24->X, A24->X, A24->Z);
    fp_add_512(A24->X, A->X, A24->X);

    int bit, nbits = 512;
    
    while(--nbits && !mp_U512_bit(k, nbits));

    do
    {
        bit = mp_U512_bit(k, nbits);
        if(bit)
        {
            fp_cpy(Q->X, tmp->X);fp_cpy(Q->Z, tmp->Z);
            fp_cpy(R->X, Q->X);fp_cpy(R->Z, Q->Z);
            fp_cpy(tmp->X, R->X);fp_cpy(tmp->Z, R->Z);
        }
        xDBLADD(Q, R, Q, R, Pcopy, A24);
        if(bit)
        {
            fp_cpy(Q->X, tmp->X);fp_cpy(Q->Z, tmp->Z);
            fp_cpy(R->X, Q->X);fp_cpy(R->Z, Q->Z);
            fp_cpy(tmp->X, R->X);fp_cpy(tmp->Z, R->Z);
        }
    } while (nbits--);

}

void xISOG(proj_point_t A, proj_point_t P, const proj_point_t K, const uint64_t k)
{
    felm_t tmp0, tmp1;
    felm_t T[4];
    fp_cpy(K->Z, T[0]);
    fp_cpy(K->X, T[1]);
    fp_cpy(K->X, T[2]);
    fp_cpy(K->Z, T[3]);
    proj_point_t Q;

    fp_mul_mont_512(P->X, K->X, Q->X);
    fp_mul_mont_512(P->Z, K->Z, tmp0);
    fp_sub_512(Q->X, tmp0, Q->X);

    fp_mul_mont_512(P->X, K->Z, Q->Z);
    fp_mul_mont_512(P->Z, K->X, tmp0);
    fp_sub_512(Q->Z, tmp0, Q->Z);

    proj_point_t M[3];
    int i;
    for (i = 0; i < 3; i++)
    {
        fp_cpy(K->X, M[i]->X);
        fp_cpy(K->Z, M[i]->Z);
    }
    xDBL(M[1], A, K);

    for (i = 1; i < k / 2; ++i) {
        if (i >= 2)
            xADD(M[i % 3], M[(i - 1) % 3], K, M[(i - 2) % 3]);

        fp_mul_mont_512(M[i % 3]->X, T[0], tmp0);
        fp_mul_mont_512(M[i % 3]->Z, T[1], tmp1);
        fp_add_512(tmp0, tmp1, T[0]);

        fp_mul_mont_512(M[i % 3]->X, T[1], T[1]);

        fp_mul_mont_512(M[i % 3]->Z, T[2], tmp0);
        fp_mul_mont_512(M[i % 3]->X, T[3], tmp1);
        fp_add_512(tmp0, tmp1, T[2]);

        fp_mul_mont_512(M[i % 3]->Z, T[3], T[3]);


        fp_mul_mont_512(P->X, M[i % 3]->X, tmp0);
        fp_mul_mont_512(P->Z, M[i % 3]->Z, tmp1);
        fp_sub_512(tmp0, tmp1, tmp0);
        fp_mul_mont_512(Q->X, tmp0, Q->X);

        fp_mul_mont_512(P->X, M[i % 3]->Z, tmp0);
        fp_mul_mont_512(P->Z, M[i % 3]->X, tmp1);
        fp_sub_512(tmp0, tmp1, tmp0);
        fp_mul_mont_512(Q->Z, tmp0, Q->Z);
    }

    fp_mul_mont_512(T[0], T[1], T[0]);
    fp_add_512(T[0], T[0], T[0]); 
    fp_sqr_mont_512(T[1], T[1]);
    fp_mul_mont_512(T[2], T[3], T[2]);
    fp_add_512(T[2], T[2], T[2]); 
    fp_sqr_mont_512(T[3], T[3]);
    fp_mul_mont_512(T[1], T[2], tmp0);
    fp_mul_mont_512(T[0], T[3], tmp1);
    fp_sub_512(tmp0, tmp1, tmp0);
    fp_mul_mont_512(tmp0, A->Z, tmp0);
    fp_add_512(tmp0, tmp0, tmp1); 
    fp_add_512(tmp0, tmp1, tmp0);
    fp_mul_mont_512(T[1], T[3], tmp1);
    fp_mul_mont_512(tmp1, A->X, tmp1);
    fp_sub_512(tmp1, tmp0, A->X);
    fp_sqr_mont_512(T[3], T[3]);
    fp_mul_mont_512(A->Z, T[3], A->Z);
    fp_sqr_mont_512(Q->X, Q->X);
    fp_sqr_mont_512(Q->Z, Q->Z);
    fp_mul_mont_512(P->X, Q->X, P->X);
    fp_mul_mont_512(P->Z, Q->Z, P->Z);
}
