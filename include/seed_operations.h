#ifndef SEED_OPERATIONS_H
#define SEED_OPERATIONS_H

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int ext_k);

bool seed_check(seed *s, int k, int gc_max_count, int hp_max);

uint64_t extend_seed(seed *s, int diff_k,int k, float gc_max, int hp_max);

#endif
