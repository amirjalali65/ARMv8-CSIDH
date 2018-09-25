#include "arith.h"
#include <assert.h>

///////////////////  Field Arithmetic  /////////////////////

uint64_t prime511[8] = { 0x1b81b90533c6c87b, 0xc2721bf457aca835,
                         0x516730cc1f0b4f25, 0xa7aac6c567f35507,
                         0x5afbfcc69322c9cd, 0xb42d083aedc88c42,
                         0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf };

uint64_t r2_Mont[8] = { 0x36905b572ffc1724, 0x67086f4525f1f27d,
                        0x4faf3fbfd22370ca, 0x192ea214bcc584b1,
                        0x5dae03ee2f5de3d0, 0x1e9248731776b371,
                        0xad5f166e20e4f52d, 0x4ed759aea6f3917e };

uint64_t one_Mont[8] = { 0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 
                         0x5d319e67c1e961b4, 0xb0aa7275301955f1,
                         0x4a080672d9ba6c64, 0x97a5ef8a246ee77b,
                         0x06ea9e5d4383676a, 0x3496e2e117e0ec80 };

uint64_t zero[8] = {0};

void fp_sqr_mont_512(const uint64_t *a, uint64_t *c){
    fp_mul_mont_512(a, a, c);
}

void fp_cpy(const uint64_t *a, uint64_t *c)
{
    int i;
    for(i = 0; i < 8; i++)
        c[i] = a[i];
}

void fp_inv(uint64_t *a)
{
    // Field inversion using FLT chain
    felm_t tmp[28], t;
    int i;

    // Table
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
    for(i = 0; i < 5; i++) fp_sqr_mont_512(t, t);
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

//////////////// Group Arithmetic ////////////////////////

void xDBLADD(proj_point_t R, proj_point_t S, const proj_point_t P, const proj_point_t Q, const proj_point_t PQ, const proj_point_t A)
{
    felm_t t0, t1, t2, t3;

    fp_add_512(Q->X, Q->Z, t0);
    fp_sub_512(Q->X, Q->Z, t1);
    fp_add_512(P->X, P->Z, t2);
    fp_sub_512(P->X, P->Z, t3);
    fp_sqr_mont_512(R->X, t2);
    fp_sqr_mont_512(S->X, t3);
    fp_mul_mont_512(t1, t2, t2);
    fp_mul_mont_512(t0, t3, t3);
    fp_sub_512(R->X, S->X, t1);
    fp_add_512(A->Z, A->Z, t0); /* multiplication by 2 */
    fp_mul_mont_512(t0, S->X, R->Z);
    fp_add_512(A->X, t0, S->X);
    fp_add_512(R->Z, R->Z, R->Z); /* multiplication by 2 */
    fp_mul_mont_512(R->X, R->Z, R->X);
    fp_mul_mont_512(t1, S->X, S->X);
    fp_sub_512(t2, t3, S->Z);
    fp_add_512(R->Z, S->X, R->Z);
    fp_add_512(t2, t3, S->X);
    fp_mul_mont_512(t1, R->Z, R->Z);
    fp_sqr_mont_512(S->Z, t3);
    fp_sqr_mont_512(S->X, t1);
    fp_mul_mont_512(PQ->Z, t1, S->X);
    fp_mul_mont_512(PQ->X, t3, S->Z);
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
    fp_add_512(A->Z, A->Z, t0); /* multiplication by 2 */
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

void point_swap(proj_point_t R, proj_point_t Q, bool swap_bit)
{
    proj_point_t tmp;
    if(swap_bit)
    {
        fp_cpy(Q->X, tmp->X);
        fp_cpy(Q->Z, tmp->Z);
        fp_cpy(R->X, Q->X);
        fp_cpy(R->Z, Q->Z);
        fp_cpy(tmp->X, R->X);
        fp_cpy(tmp->Z, R->Z);
    }
}


/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
/*
void xMUL(proj_point_t Q, const proj_point_t A, const proj_point_t P, const UINT512_t k)
{
    proj_point_t R;
    proj_point_t Pcopy; // in case Q = P 

    fp_cpy(P->X, R->X);
    fp_cpy(P->Z, R->Z);
    fp_cpy(P->X, Pcopy->X);
    fp_cpy(P->Z, Pcopy->Z);


    fp_cpy(one_Mont, Q->X);
    fp_cpy(zero, Q->Z);

    unsigned long i = 512;

    // TODO: make it using a foor loop
    // TODO: implement u512_bit: it returns the bit value of k at position i
    while (--i && !u512_bit(k, i));

    do {
        bool bit = u512_bit(k, i);
        point_swap(R, Q, bit);
        xDBLADD(Q, R, Q, R, Pcopy, A);
        point_swap(R, Q, bit);
    } while (i--);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k 
void xISOG(proj_point_t A, proj_point_t P, const proj_point_t K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

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
    fp_add_512(T[0], T[0], T[0]); // multiplication by 2 

    fp_sqr_mont_512(T[1], T[1]);

    fp_mul_mont_512(T[2], T[3], T[2]);
    fp_add_512(T[2], T[2], T[2]); // multiplication by 2 

    fp_sqr_mont_512(T[3], T[3]);

    // Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) 
    fp_mul_mont_512(T[1], T[2], tmp0);
    fp_mul_mont_512(T[0], T[3], tmp1);
    fp_sub_512(tmp0, tmp1, tmp0);
    fp_mul_mont_512(tmp0, A->Z, tmp0);
    fp_add_512(tmp0, tmp0, tmp1); 
    fp_add_512(tmp0, tmp1, tmp0); // multiplication by 3 

    fp_mul_mont_512(T[1], T[3], tmp1);
    fp_mul_mont_512(tmp1, A->X, tmp1);

    fp_sub_512(tmp1, tmp0, A->X);

    // Az := Az * T[3]^2 
    fp_sqr_mont_512(T[3], T[3]);
    fp_mul_mont_512(A->Z, T[3], A->Z);

    // X := X * Xim^2, Z := Z * Zim^2 
    fp_sqr_mont_512(Q->X, Q->X);
    fp_sqr_mont_512(Q->Z, Q->Z);
    fp_mul_mont_512(P->X, Q->X, P->X);
    fp_mul_mont_512(P->Z, Q->Z, P->Z);
}
*/