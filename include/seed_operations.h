#ifndef SEED_OPERATIONS_H
#define SEED_OPERATIONS_H

#include "seed.h"

void initialize_seeds(seed *seeds, int ext_k, int n);

void seed_push_base(seed *s, uint64_t base);

uint64_t seed_get_base(const seed *s, int pos);

void seed_destroy(seed *seeds, int n);

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int k);

bool seed_check(seed *s, int k, int gc_max_count, int hp_max);

void extend_seed(seed *s, int diff_k,int k, float gc_max, int hp_max);

#endif
