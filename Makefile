CC=gcc
CPP=g++
CFLAGS=-fsanitize=address -Wall -Wextra -Wno-reorder
CPPFLAGS=-std=c++17

NAME=useful

all: ${NAME}

%: %.c
	${CC} ${CFLAGS} $< -o $@

%: %.cpp
	${CPP} ${CPPFLAGS} ${CFLAGS} $< -o $@

clean:
	rm -rf ${NAME}