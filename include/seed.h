#ifndef SEED_H
#define SEED_H

#include "nn_table.h"

typedef struct seed{
        uint64_t *blocks;
        int length;
        int n_blocks;
	int gc;
	int hp_count;
	uint64_t last_base;
	TmAccumulator acc;
} seed;

#endif
