#ifndef SEED_OPERATIONS_H
#define SEED_OPERATIONS_H

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int ext_k);

bool seed_check(seed *s, int k, int gc_max_count, int gc_min_count, int hp_max);

#endif
