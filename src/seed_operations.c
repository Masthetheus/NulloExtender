#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "seed.h"
#include "seed_operations.h"
#include "filters.h"

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int ext_k){
	uint64_t random = rand() % *count;
	if (nullomers[random]){
		seeds[i].idx = nullomers[random];
	} else {
		printf("No more suitable seeds are available, try lowering seed number.");
		exit(EXIT_FAILURE);	
	}
	uint64_t holder = nullomers[random];
	nullomers[random] = nullomers[*count - 1];
	nullomers[*count - 1] = holder;
	*count--;
}

bool seed_check(seed *s, int k, int gc_max_count, int gc_min_count, int hp_max){
	int gc = 0;
	int hp_count = 1;
	uint64_t idx = s->idx;
	
	uint64_t first_base = (idx >> ((k-1)*2)) & 3;
	
	if (first_base == 1 || first_base == 3){
		s->gc++;
	}	
	s->last_base = first_base;

	for (int i = 1; i < k; i++){
		uint64_t base = (idx >> ((k-i-1)*2)) & 3;
		bool gc_checker = gc_check(s, base, gc_max_count, gc_min_count);
		bool hp_checker = hp_check(s, base, hp_max);
		if (gc_checker == false || hp_checker == false){
			return false;
		}
	}

	return true;
}

