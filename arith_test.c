#include <stdio.h>
#include <string.h>
#include "arith.h"
#include "rng.h"
#include <string.h>

#define TEST_LOOP 1000


void print_big(uint64_t *a)
{
    int i;
    printf("\n");
    for(i = 7; i >= 0; i--)
    {
        printf("%llx", a[i]);
    }
    printf("\n");
}

int test_fp_arithmetic()
{
    int i, passed = 1;
    felm_t a, b, c, d, c1, c2, c3, c4;

    for(i = 0; i < TEST_LOOP; i++)
    {
        randombytes(a, 64);randombytes(b, 64);
        randombytes(c, 64);randombytes(d, 64);
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

        fp_add_512(a, b, c1);   
        fp_inv(c);              
        fp_mul_mont_512(c, c1, c1); 

        fp_mul_mont_512(a, c, c2);  
        fp_mul_mont_512(b, c, c3);  
        fp_add_512(c2, c3, c2);     

        if(memcmp(c1, c2, 64) != 0)
            passed = 0;
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
}