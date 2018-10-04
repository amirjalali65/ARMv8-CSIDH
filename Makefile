all:
	@cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-O2 -funroll-loops \
		rng.c \
		u512.s fp.s \
		mont.c \
		csidh.c \
		main.c \
		-o main

debug:
	cc \
		-std=c99 -pedantic \
		-Wall -Wextra \
		-g \
		rng.c \
		u512.s fp.s \
		mont.c \
		csidh.c \
		main.c \
		-o main
arm_64:
	aarch64-linux-gnu-gcc \
		-std=c99 -static \
		-Wall -Wextra	\
		-O3	-funroll-loops	\
		rng.c \
		arith.c \
		arith.S	\
		csidh_api.c \
		csidh_test.c \
		-o CSIDH_test

clean:
	rm -f main 

