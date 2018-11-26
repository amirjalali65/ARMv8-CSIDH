/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*
*   Author: Amir Jalali                     ajalali2016@fau.edu
*                       
*                       All rights reserved   
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "arith.h"
#include "rng.h"
#include <string.h>

#define TEST_LOOP 1000

int test_fp_arithmetic()
{
    int i, passed = 1;
    felm_t a={0}, b={0}, acpy={0}, bcpy={0}, c={0}, d={0}, c1={0}, c2={0}, c3={0}, c4={0};

    for(i = 0; i < TEST_LOOP; i++)
    {
        fp_random_512(a);fp_random_512(b);
        fp_random_512(c);fp_random_512(d);
        fp_cpy(a, acpy);fp_cpy(b, bcpy);
        to_mont(a, a);to_mont(b, b);
        to_mont(c, c);to_mont(d, d);
        
        fp_add_512(a, b, c1);           
        fp_add_512(c, d, c2);           
        fp_mul_mont_512(c1, c2, c3);    
        fp_mul_mont_512(a, c, c1);
        fp_mul_mont_512(a, d, c2);
        fp_add_512(c1, c2, c1);
        fp_mul_mont_512(b, c, c2);
        fp_mul_mont_512(b, d, c4);
        fp_add_512(c2, c4, c2);
        fp_add_512(c2, c1, c2);
        if(memcmp(c2, c3, 64) != 0)
            passed = 0;

        fp_sub_512(a, b, c1);           
        fp_sub_512(c, d, c2);           
        fp_mul_mont_512(c1, c2, c3);    
        fp_mul_mont_512(a, c, c1);  
        fp_mul_mont_512(a, d, c2);  
        fp_sub_512(c1, c2, c1);     
        fp_mul_mont_512(b, c, c2);  
        fp_mul_mont_512(b, d, c4);  
        fp_sub_512(c2, c4, c2);     
        fp_sub_512(c1, c2, c2);     
        if(memcmp(c2, c3, 64) != 0)
            passed = 0;

        fp_init_one(c);
        fp_cpy(c, d);
        fp_inv(c);
        fp_mul_mont_512(c, d, c1);
        fp_mul_mont_512(c1, one_Mont, c1);
        if(memcmp(c1, one_Mont, 64) != 0)
        {
            fp_print(c1);
            fp_print(one_Mont);
            passed = 0;
        }
    }



    return passed;
}

int main()
{
    int passed;

    passed = test_fp_arithmetic();

    if(passed)
    {
        printf("\nfp arithmetic tests passed\n");
    }else{
        printf("\nfp arithmetic tests failed\n");
    }
    return 0;
}