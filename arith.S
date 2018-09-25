/****************************************************************************
*   Efficient implementation of finite field arithmetic over p511 on ARMv8
*                   Constant-time Implementation of CSIDH
*
*   Amir Jalali                     ajalali2016@fau.edu
*****************************************************************************/
.data
p511:
.quad 0x1b81b90533c6c87b
.quad 0xc2721bf457aca835
.quad 0x516730cc1f0b4f25
.quad 0xa7aac6c567f35507
.quad 0x5afbfcc69322c9cd
.quad 0xb42d083aedc88c42
.quad 0xfc8ab0d15e3e4c4a
.quad 0x65b48e8f740f89bf

minus_p511_inverse:
.quad 0x66c1301f632e294d

.text

.macro schoolbook_compact_registers_first
    mul     x10, x3, x7
    umulh   x12, x3, x7
    
    mul     x14, x4, x7
    umulh   x16, x4, x7
    
    mul     x18, x5, x7
    umulh   x20, x5, x7

    mul     x22, x6, x7
    umulh   x24, x6, x7

    adds    x12, x12, x14
    adcs    x16, x16, x18
    adcs    x20, x20, x22

    ldp     x3, x4, [x0, #32]
    ldp     x5, x6, [x0, #48]

    mul     x25, x3, x7
    umulh   x26, x3, x7

    mul     x27, x4, x7
    umulh   x14, x4, x7

    mul     x18, x5, x7
    umulh   x22, x5, x7

    mul     x3, x6, x7
    umulh   x4, x6, x7

    adcs    x24, x24, x25
    adcs    x26, x26, x27
    adcs    x14, x14, x18
    adcs    x22, x22, x3
    adc     x4, x4, xzr
.endm

.macro schoolbook_compact_registers_second
    ldr     x3, p511
    ldr     x4, p511 + 8
    ldr     x5, p511 + 16
    ldr     x6, p511 + 24

    mul     x10, x3, x28
    umulh   x12, x3, x28
    
    mul     x14, x4, x28
    umulh   x16, x4, x28
    
    mul     x18, x5, x28
    umulh   x20, x5, x28

    mul     x22, x6, x28
    umulh   x24, x6, x28

    adds    x12, x12, x14
    adcs    x16, x16, x18
    adcs    x20, x20, x22

    ldr     x3, p511 + 32
    ldr     x4, p511 + 40
    ldr     x5, p511 + 48
    ldr     x6, p511 + 56
    
    mul     x25, x3, x28
    umulh   x26, x3, x28

    mul     x27, x4, x28
    umulh   x14, x4, x28

    mul     x18, x5, x28
    umulh   x22, x5, x28 
    
    mul     x3, x6, x28
    umulh   x4, x6, x28 

    adcs    x24, x24, x25
    adcs    x26, x26, x27
    adcs    x14, x14, x18
    adcs    x22, x22, x3
    adc     x4, x4, xzr
.endm

.global fp_add_512
.global fp_sub_512
.global mp_add_512
.global mp_sub_512
.global mp_mul_u64
.global fp_mul_mont_512


fp_add_512:
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0,#16]
    ldp     x7, x8,   [x0,#32]
    ldp     x9, x10,  [x0,#48]
     
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1,#16]
    ldp     x15, x16, [x1,#32]
    ldp     x17, x18, [x1,#48]  

    adds    x3, x3, x11
    adcs    x4, x4, x12
    adcs    x5, x5, x13
    adcs    x6, x6, x14
    adcs    x7, x7, x15
    adcs    x8, x8, x16
    adcs    x9, x9, x17
    adcs    x10, x10, x18

    ldr     x11, p511
    ldr     x12, p511 + 8
    ldr     x13, p511 + 16
    ldr     x14, p511 + 24
    ldr     x15, p511 + 32
    ldr     x16, p511 + 40
    ldr     x17, p511 + 48
    ldr     x18, p511 + 56

    subs    x3, x3, x11
    sbcs   x4, x4, x12
    sbcs   x5, x5, x13
    sbcs   x6, x6, x14
    sbcs   x7, x7, x15
    sbcs   x8, x8, x16
    sbcs   x9, x9, x17
    sbcs   x10, x10, x18
    sbc     x19, xzr, xzr

    and     x11, x11, x19
    adds    x3, x3, x11
    and     x12, x12, x19
    adcs    x4, x4, x12
    and     x13, x13, x19
    adcs    x5, x5, x13
    and     x14, x14, x19
    adcs    x6, x6, x14
    and     x15, x15, x19
    adcs    x7, x7, x15
    and     x16, x16, x19
    adcs    x8, x8, x16
    and     x17, x17, x19
    adcs    x9, x9, x17 
    and     x18, x18, x19
    adcs    x10, x10, x18

    stp     x3, x4, [x2]
    stp     x5, x6, [x2,#16]
    stp     x7, x8, [x2,#32]
    stp     x9, x10, [x2,#48]
    
    ret

fp_sub_512:
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0,#16]
    ldp     x7, x8,   [x0,#32]
    ldp     x9, x10,  [x0,#48]
     
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1,#16]
    ldp     x15, x16, [x1,#32]
    ldp     x17, x18, [x1,#48]

    subs    x3, x3, x11
    sbcs    x4, x4, x12
    sbcs    x5, x5, x13
    sbcs    x6, x6, x14
    sbcs    x7, x7, x15
    sbcs    x8, x8, x16
    sbcs    x9, x9, x17
    sbcs    x10, x10, x18
    sbc     x19, xzr, xzr

    ldr     x11, p511
    ldr     x12, p511 + 8
    ldr     x13, p511 + 16
    ldr     x14, p511 + 24
    ldr     x15, p511 + 32
    ldr     x16, p511 + 40
    ldr     x17, p511 + 48
    ldr     x18, p511 + 56

    and     x11, x11, x19
    adds    x3, x3, x11
    and     x12, x12, x19
    adcs    x4, x4, x12
    and     x13, x13, x19
    adcs    x5, x5, x13
    and     x14, x14, x19
    adcs    x6, x6, x14
    and     x15, x15, x19
    adcs    x7, x7, x15
    and     x16, x16, x19
    adcs    x8, x8, x16
    and     x17, x17, x19
    adcs    x9, x9, x17 
    and     x18, x18, x19
    adcs    x10, x10, x18
    
    stp     x3, x4, [x2]
    stp     x5, x6, [x2,#16]
    stp     x7, x8, [x2,#32]
    stp     x9, x10, [x2,#48]
    
    ret

mp_add_512:
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0,#16]
    ldp     x7, x8,   [x0,#32]
    ldp     x9, x10,  [x0,#48]
     
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1,#16]
    ldp     x15, x16, [x1,#32]
    ldp     x17, x18, [x1,#48]

    adds    x3, x3, x11
    adcs    x4, x4, x12
    adcs    x5, x5, x13
    adcs    x6, x6, x14
    adcs    x7, x7, x15
    adcs    x8, x8, x16
    adcs    x9, x9, x17
    adc    x10, x10, x18

    stp     x3, x4, [x2]
    stp     x5, x6, [x2,#16]
    stp     x7, x8, [x2,#32]
    stp     x9, x10, [x2,#48]
    
    ret

// TODO: Check if negative results need to be implemented in constant-time
mp_sub_512:
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0,#16]
    ldp     x7, x8,   [x0,#32]
    ldp     x9, x10,  [x0,#48]
     
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1,#16]
    ldp     x15, x16, [x1,#32]
    ldp     x17, x18, [x1,#48]

    subs    x3, x3, x11
    sbcs    x4, x4, x12
    sbcs    x5, x5, x13
    sbcs    x6, x6, x14
    sbcs    x7, x7, x15
    sbcs    x8, x8, x16
    sbcs    x9, x9, x17
    sbc    x10, x10, x18

    stp     x3, x4, [x2]
    stp     x5, x6, [x2,#16]
    stp     x7, x8, [x2,#32]
    stp     x9, x10, [x2,#48]
    
    ret

mp_mul_u64:
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0,#16]
    ldp     x7, x8,   [x0,#32]
    ldp     x9, x10,  [x0,#48]
    
    ldr     x11, [x1]

    mul     x12, x3, x11
    umulh   x13, x3, x11

    mul     x14, x4, x11
    umulh   x15, x4, x11         

    mul     x16, x5, x11
    umulh   x17, x5, x11         

    mul     x18, x6, x11
    umulh   x19, x6, x11         

    mul     x20, x7, x11
    umulh   x21, x7, x11         

    mul     x22, x8, x11
    umulh   x23, x8, x11         

    mul     x24, x9, x11
    umulh   x25, x9, x11         

    mul     x26, x10, x11
    umulh   x27, x10, x11

    adds    x13, x13, x14
    adcs    x15, x15, x16
    adcs    x17, x17, x18
    adcs    x19, x19, x20
    adcs    x21, x21, x22
    adcs    x23, x23, x24
    adcs    x25, x25, x26
    adc     x27, xzr, xzr 

    stp     x13, x15, [x2]
    stp     x17, x19, [x2, #16]
    stp     x21, x23, [x2, #32]
    stp     x25, x27, [x2, #48]

    ret        

fp_mul_mont_512:
    sub     sp, sp, #48
    stp     x15, x16, [sp]
    stp     x17, x18, [sp, #16]
    stp     x19, x20, [sp, #32]

    // 0
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1]

    mul     x8, x3, x7
    umulh   x9, x3, x7
    
    mul     x10, x4, x7
    umulh   x11, x4, x7
    
    mul     x12, x5, x7
    umulh   x13, x5, x7

    mul     x14, x6, x7
    umulh   x15, x6, x7

    adds    x9, x9, x10
    adcs    x11, x11, x12
    adcs    x13, x13, x14

    ldp     x3, x4, [x0, #32]
    ldp     x5, x6, [x0, #48]

    mul     x16, x3, x7
    umulh   x17, x3, x7

    mul     x18, x4, x7
    umulh   x19, x4, x7

    mul     x20, x5, x7
    umulh   x21, x5, x7

    mul     x22, x6, x7
    umulh   x23, x6, x7

    adcs    x15, x15, x16
    adcs    x17, x17, x18
    adcs    x19, x19, x20
    adcs    x21, x21, x22
    adc     x23, x23, xzr

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x8

    ldr     x3, p511
    ldr     x4, p511 + 8
    ldr     x5, p511 + 16
    ldr     x6, p511 + 24

    mul     x10, x3, x28
    umulh   x12, x3, x28

    mul     x14, x4, x28
    umulh   x16, x4, x28

    mul     x18, x5, x28
    umulh   x20, x5, x28

    mul     x22, x6, x28
    umulh   x24, x6, x28

    adds    x12, x12, x14
    adcs    x16, x16, x18
    adcs    x20, x20, x22

    ldr     x3, p511 + 32
    ldr     x4, p511 + 40
    ldr     x5, p511 + 48
    ldr     x6, p511 + 56
    
    mul     x25, x3, x28
    umulh   x26, x3, x28

    mul     x27, x4, x28
    umulh   x14, x4, x28

    mul     x18, x5, x28
    umulh   x22, x5, x28 
    
    mul     x3, x6, x28
    umulh   x4, x6, x28 

    adcs    x24, x24, x25
    adcs    x26, x26, x27
    adcs    x14, x14, x18
    adcs    x22, x22, x3
    adc     x4, x4, xzr

    adds    x8, x8, x10
    adcs    x9, x9, x12
    adcs    x11, x11, x16
    adcs    x13, x13, x20
    adcs    x15, x15, x24
    adcs    x17, x17, x26
    adcs    x19, x19, x14
    adcs    x21, x21, x22
    adcs    x23, x23, x4
    adc     x8, xzr, xzr    
    
    // 1
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #8]

    schoolbook_compact_registers_first

    adds    x9, x9, x10
    adcs    x11, x11, x12
    adcs    x13, x13, x16
    adcs    x15, x15, x20
    adcs    x17, x17, x24
    adcs    x19, x19, x26
    adcs    x21, x21, x14
    adcs    x23, x23, x22
    adc     x8, x8, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x9
    
    schoolbook_compact_registers_second

    adds    x9, x9, x10
    adcs    x11, x11, x12
    adcs    x13, x13, x16
    adcs    x15, x15, x20
    adcs    x17, x17, x24
    adcs    x19, x19, x26
    adcs    x21, x21, x14
    adcs    x23, x23, x22
    adcs    x8, x8, x4
    adc     x9, xzr, xzr   
    
    // 2
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #16]

    schoolbook_compact_registers_first

    adds    x11, x11, x10
    adcs    x13, x13, x12
    adcs    x15, x15, x16
    adcs    x17, x17, x20
    adcs    x19, x19, x24
    adcs    x21, x21, x26
    adcs    x23, x23, x14
    adcs    x8, x8, x22
    adc     x9, x9, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x11

    schoolbook_compact_registers_second

    adds    x11, x11, x10
    adcs    x13, x13, x12
    adcs    x15, x15, x16
    adcs    x17, x17, x20
    adcs    x19, x19, x24
    adcs    x21, x21, x26
    adcs    x23, x23, x14
    adcs    x8, x8, x22
    adcs    x9, x9, x4
    adc     x11, xzr, xzr  

    // 3
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #24]

    schoolbook_compact_registers_first

    adds    x13, x13, x10
    adcs    x15, x15, x12
    adcs    x17, x17, x16
    adcs    x19, x19, x20
    adcs    x21, x21, x24
    adcs    x23, x23, x26
    adcs    x8, x8, x14
    adcs    x9, x9, x22
    adc     x11, x11, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x13

    schoolbook_compact_registers_second

    adds    x13, x13, x10
    adcs    x15, x15, x12
    adcs    x17, x17, x16
    adcs    x19, x19, x20
    adcs    x21, x21, x24
    adcs    x23, x23, x26
    adcs    x8, x8, x14
    adcs    x9, x9, x22
    adcs    x11, x11, x4
    adc     x13, xzr, xzr   

    // 4
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #32]

    schoolbook_compact_registers_first

    adds    x15, x15, x10
    adcs    x17, x17, x12
    adcs    x19, x19, x16
    adcs    x21, x21, x20
    adcs    x23, x23, x24
    adcs    x8, x8, x26
    adcs    x9, x9, x14
    adcs    x11, x11, x22
    adc     x13, x13, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x15

    schoolbook_compact_registers_second
  
    adds    x15, x15, x10
    adcs    x17, x17, x12
    adcs    x19, x19, x16
    adcs    x21, x21, x20
    adcs    x23, x23, x24
    adcs    x8, x8, x26
    adcs    x9, x9, x14
    adcs    x11, x11, x22
    adcs    x13, x13, x4
    adc     x15, xzr, xzr   

    // 5
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #40]

    schoolbook_compact_registers_first

    adds    x17, x17, x10
    adcs    x19, x19, x12
    adcs    x21, x21, x16
    adcs    x23, x23, x20
    adcs    x8, x8, x24
    adcs    x9, x9, x26
    adcs    x11, x11, x14
    adcs    x13, x13, x22
    adc     x15, x15, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x17

    schoolbook_compact_registers_second 

    adds    x17, x17, x10
    adcs    x19, x19, x12
    adcs    x21, x21, x16
    adcs    x23, x23, x20
    adcs    x8, x8, x24
    adcs    x9, x9, x26
    adcs    x11, x11, x14
    adcs    x13, x13, x22
    adcs    x15, x15, x4
    adc     x17, xzr, xzr   

    // 6
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #48]

    schoolbook_compact_registers_first

    adds    x19, x19, x10
    adcs    x21, x21, x12
    adcs    x23, x23, x16
    adcs    x8, x8, x20
    adcs    x9, x9, x24
    adcs    x11, x11, x26
    adcs    x13, x13, x14
    adcs    x15, x15, x22
    adc     x17, x17, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x19

    schoolbook_compact_registers_second
  
    adds    x19, x19, x10
    adcs    x21, x21, x12
    adcs    x23, x23, x16
    adcs    x8, x8, x20
    adcs    x9, x9, x24
    adcs    x11, x11, x26
    adcs    x13, x13, x14
    adcs    x15, x15, x22
    adcs    x17, x17, x4
    adc     x19, xzr, xzr   

    // 7
    ldp     x3, x4, [x0]
    ldp     x5, x6, [x0, #16]
    ldr     x7, [x1, #56]

    schoolbook_compact_registers_first

    adds    x21, x21, x10
    adcs    x23, x23, x12
    adcs    x8, x8, x16
    adcs    x9, x9, x20
    adcs    x11, x11, x24
    adcs    x13, x13, x26
    adcs    x15, x15, x14
    adcs    x17, x17, x22
    adc     x19, x19, x4

    ldr     x28, minus_p511_inverse
    mul     x28, x28, x21

    schoolbook_compact_registers_second
    
    adds    x21, x21, x10
    adcs    x23, x23, x12   // mc0
    adcs    x8, x8, x16     // mc1
    adcs    x9, x9, x20     // mc2
    adcs    x11, x11, x24   // mc3
    adcs    x13, x13, x26   // mc4
    adcs    x15, x15, x14   // mc5
    adcs    x17, x17, x22   // mc6
    adc    x19, x19, x4     // mc7

    ldr     x10, p511
    ldr     x12, p511 + 8
    ldr     x14, p511 + 16
    ldr     x16, p511 + 24

    ldr     x3, p511 + 32
    ldr     x4, p511 + 40
    ldr     x5, p511 + 48
    ldr     x6, p511 + 56
    
    subs    x23, x23, x10
    sbcs    x8, x8, x12
    sbcs    x9, x9, x14
    sbcs    x11, x11, x16
    sbcs    x13, x13, x3
    sbcs    x15, x15, x4
    sbcs    x17, x17, x5
    sbcs    x19, x19, x6
    sbc     x28, xzr, xzr

    and     x10, x10, x28
    and     x12, x12, x28
    and     x14, x14, x28
    and     x16, x16, x28
    and     x3, x3, x28
    and     x4, x4, x28
    and     x5, x5, x28
    and     x6, x6, x28

    adds    x23, x23, x10
    adcs    x8, x8, x12
    adcs    x9, x9, x14
    adcs    x11, x11, x16
    adcs    x13, x13, x3
    adcs    x15, x15, x4
    adcs    x17, x17, x5
    adcs    x19, x19, x6
    
    stp     x23, x8, [x2]
    stp     x9, x11, [x2, #16]
    stp     x13, x15, [x2, #32]
    stp     x17, x19, [x2, #48]

    ldp     x15, x16, [sp]
    ldp     x17, x18, [sp, #16]
    ldp     x19, x20, [sp, #32]
    add     sp, sp, #48

    ret


