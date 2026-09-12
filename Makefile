CC = gcc
CFLAGS = -Wall -Wextra -O3 -Iinclude
LDFLAGS = -lm

SRC = $(wildcard src/*.c)
OBJ = $(SRC:src/%.c=build/%.o)
TARGET = bin/nullo_extender

$(TARGET): $(OBJ)
	@mkdir -p bin
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build bin/*

.PHONY: clean
