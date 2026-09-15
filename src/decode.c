#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "seed_operations.h"

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

void decode_kmer(seed *s, int k, char *seq) {
        for (int i = 0; i < s->length; i++) {
                seq[i] = BITS_TO_BASE[seed_get_base(s, i)];
        }
        seq[s->length] = '\0';
}
