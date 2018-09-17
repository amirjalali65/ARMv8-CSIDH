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

clean:
	rm -f main

