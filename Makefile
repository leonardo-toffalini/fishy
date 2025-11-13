CC = clang
CFLAGS = -Wextra -Wall

.PHONY = run, clean

all: main run

main:
	CC src/main.c -o bin/main

run:
	./bin/main

clean:
	rm bin/*
