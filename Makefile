# Makefile for cross-compilation on LINUX 
# Target processor: aarch64 ARMv8  

OPTIMIZATION = -O3

# Default cross-sompiler 
CC=aarch64-linux-gnu-gcc-7
CROSS_FLAGS= -static

# Default value set to non-constant time implementation
ifeq "$(CONSTANT)" "TRUE"
	CONST=-D _CONSTANT_
endif

ifeq "$(DEBUG)" "TRUE"
	DEB=-g
endif

CFLAGS= -c $(DEB) $(OPTIMIZATION) $(CROSS_FLAGS) $(CONST) 

OBJECTS=arith.o arith_asm.o csidh_api.o rng.o csidh_test.o

CSIDH_TEST: $(OBJECTS)
	$(CC) $(CROSS_FLAGS) $(OPTIMIZATION) $(ADDITIONAL_FLAGS) -o CSIDH_TEST $(OBJECTS) $(TEST_OBJECTS)

arith.o: arith.c arith.h
	$(CC) $(CFLAGS) arith.c

arith_asm.o: arith_asm.S
	$(CC) $(CFLAGS) arith_asm.S

csidh_api.o: csidh_api.c csidh_api.h
	$(CC) $(CFLAGS) csidh_api.c

rng.o: rng.c rng.h
	$(CC) $(CFLAGS) rng.c

csidh_test.o: csidh_test.c
	$(CC) $(CFLAGS) csidh_test.c

.PHONY: clean

clean:
	rm $(OBJECTS) $(TEST_OBJECT) CSIDH_TEST