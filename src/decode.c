#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

void decode_kmer(uint64_t idx, int k, char *seq) {
    for (int i = 0; i < k; i++) {
        uint64_t base = (idx >> ((k-i-1)*2)) & 3;
        seq[i] = BITS_TO_BASE[base];
    }
    seq[k] = '\0';
}
