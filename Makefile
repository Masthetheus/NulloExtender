CC = gcc
BASE_CFLAGS = -Wall -Wextra -Iinclude
BASE_LDFLAGS = -lm

SAN_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer -g

SRC = $(wildcard src/*.c)
OBJ = $(SRC:src/%.c=build/%.o)
TARGET = bin/nullo_extender

ifeq ($(DEBUG),1)
  CFLAGS = $(BASE_CFLAGS) $(SAN_FLAGS) -O1
  LDFLAGS = $(BASE_LDFLAGS) $(SAN_FLAGS)
else
  CFLAGS = $(BASE_CFLAGS) -O3
  LDFLAGS = $(BASE_LDFLAGS)
endif

$(TARGET): $(OBJ)
	@mkdir -p bin
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build bin/*

.PHONY: clean
