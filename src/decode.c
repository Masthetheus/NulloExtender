#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "seed_operations.h"

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

void decode_kmer(seed *s, int k, char *seq) {
        fprintf(stderr,"%d\n",s->length);
        for (int i = 0; i < s->length; i++) {
                fprintf(stderr,"BITS TO BASE %d\n",BITS_TO_BASE[seed_get_base(s, i)]);
                seq[i] = BITS_TO_BASE[seed_get_base(s, i)];
                fprintf(stderr,"%s\n", seq[i]);
        }
        seq[s->length] = '\0';
}
