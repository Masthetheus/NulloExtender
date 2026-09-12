#ifndef FILTERS_H
#define FILTERS_H

bool gc_check(seed *s, uint64_t base, int max_gc_count);

bool hp_check(seed *s, uint64_t base, int hp_max);

int hamming_check(uint64_t a, uint64_t b, int ext_k);

#endif
