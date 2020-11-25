# Makefile for cross-compilation on LINUX 
# Target processor: aarch64 ARMv8  

OPTIMIZATION = -O3

BITS?= 512
CONSTANT?=TRUE
FASTLADDER?=FALSE

# Default cross-sompiler 
CC=aarch64-linux-gnu-gcc
CROSS_FLAGS= -static -fcommon

# Default value set to non-constant time implementation
ifeq "$(CONSTANT)" "TRUE"
	ifeq "$(FASTLADDER)" "TRUE"
	CONST=-D _CONSTANT_ -D _FASTLADDER_
	else
	CONST=-D _CONSTANT_
	endif
endif

ifeq "$(DEBUG)" "TRUE"
	DEB=-g
endif

CFLAGS= -c $(DEB) $(OPTIMIZATION) $(CROSS_FLAGS) $(CONST)

COMMON_OBJECTS=arith.o arith_asm.o csidh_api.o rng.o
OBJECTS=$(COMMON_OBJECTS) csidh_test.o
OBJECTS_UTIL=$(COMMON_OBJECTS) csidh_util.o

test: CSIDH_TEST
CSIDH_TEST: $(OBJECTS)
	$(CC) $(CROSS_FLAGS) $(OPTIMIZATION) $(ADDITIONAL_FLAGS) -o CSIDH_TEST $(OBJECTS) $(TEST_OBJECTS)

util: CSIDH_UTIL
CSIDH_UTIL: $(OBJECTS_UTIL)
	$(CC) $(CROSS_FLAGS) $(OPTIMIZATION) $(ADDITIONAL_FLAGS) -o csidh-p$(BITS)-util $(OBJECTS_UTIL)

regenerate_test_vectors: util
	./csidh-p512-util -g -p sample-keys/1.montgomery.le.pk -s sample-keys/1.montgomery.le.sk
	./csidh-p512-util -g -p sample-keys/2.montgomery.le.pk -s sample-keys/2.montgomery.le.sk
	./csidh-p512-util -g -p sample-keys/3.montgomery.le.pk -s sample-keys/3.montgomery.le.sk
	./csidh-p512-util -g -p sample-keys/4.montgomery.le.pk -s sample-keys/4.montgomery.le.sk
	./csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/1.montgomery.le.sk > sample-keys/1-2.ss
	./csidh-p512-util -d -p sample-keys/1.montgomery.le.pk -s sample-keys/2.montgomery.le.sk > sample-keys/2-1.ss
	./csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/3.montgomery.le.sk > sample-keys/3-2.ss
	./csidh-p512-util -d -p sample-keys/3.montgomery.le.pk -s sample-keys/4.montgomery.le.sk > sample-keys/4-3.ss

util_test: util
	echo "BEGIN util-test"
	./csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/1.montgomery.le.sk > sample-keys/1-2.ss.test_result
	diff sample-keys/1-2.ss.test_result sample-keys/1-2.ss
	./csidh-p512-util -d -p sample-keys/1.montgomery.le.pk -s sample-keys/2.montgomery.le.sk > sample-keys/2-1.ss.test_result
	diff sample-keys/2-1.ss.test_result sample-keys/2-1.ss
	./csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/3.montgomery.le.sk > sample-keys/3-2.ss.test_result
	diff sample-keys/3-2.ss.test_result sample-keys/3-2.ss
	./csidh-p512-util -d -p sample-keys/3.montgomery.le.pk -s sample-keys/4.montgomery.le.sk > sample-keys/4-3.ss.test_result
	diff sample-keys/4-3.ss.test_result sample-keys/4-3.ss
	rm sample-keys/*.test_result
	echo "END util-test"

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

csidh_util.o: csidh_util.c csidh_util.h
	$(CC) $(CFLAGS) -D BITS=$(BITS) csidh_util.c

.PHONY: clean

clean:
	-rm $(OBJECTS) $(TEST_OBJECT) $(OBJECTS_UTIL) CSIDH_TEST csidh-p$(BITS)-util sample-keys/*.test_result
