# Example

## 1. Description

This folder contains the compact nullomer set for the B. subtilis organism at a k of 12, consisting of both MAWs and trivial extensions for testing the script baseline usage.

## 2. Usage

Firstly, make sure the script is correctly compiled inside bin.

```
gcc -O3 -std=c99 -o nullo_extender src/main.c -lm
```

Then, call the binary pointing to the example file and with parametrs of your choice, such as:

```
./nullo_extender example/bsubtilis_k12_example 60 3 12 20
```

This command will output directly to the terminal 12 sequences of size 20 that are nullomers, via extension rule, of the B. subtilis with at maximum 60% GC composition and no more than 3 subsequential equal bases.
