#ifndef FILTERS_H
#define FILTERS_H

#include "nn_table.h"

bool gc_check(seed *s, uint8_t base, int max_gc_count);

bool hp_check(seed *s, uint8_t base, int hp_max);

int hamming_check(const seed *a, const seed *b);

void tm_add_base(TmAccumulator *acc, uint8_t prev_base, uint8_t curr_base);

double tm_finalize(TmAccumulator *acc, uint8_t first_base, uint8_t last_base, double na_conc, int seq_len);

#endif
