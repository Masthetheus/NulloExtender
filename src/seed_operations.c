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
#include "nn_table.h"

#define ALPHABET_SIZE 4 // A, C, T AND G

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int ext_k){
	uint64_t random = rand() % *count;
	if (nullomers[random]){
		seeds[i].idx = nullomers[random];
		seeds[i].gc = 0;
		seeds[i].hp_count = 1;
		seeds[i].last_base = (seeds[i].idx >>((ext_k-1)*2)) & 3;
		seeds[i].acc = (TmAccumulator){0};
	} else {
		printf("No more suitable seeds are available, try lowering seed number.");
		exit(EXIT_FAILURE);	
	}
	uint64_t holder = nullomers[random];
	nullomers[random] = nullomers[*count - 1];
	nullomers[*count - 1] = holder;
	(*count)--;

}

bool seed_check(seed *s, int k, int gc_max_count, int hp_max){
	uint64_t idx = s->idx;
	
	uint64_t first_base = (idx >> ((k-1)*2)) & 3;
	
	if (first_base == 1 || first_base == 3){
		s->gc++;
	}	

	for (int i = 1; i < k; i++){
		uint64_t base = (idx >> ((k-i-1)*2)) & 3;
		uint64_t prev_base = s->last_base;
		bool gc_checker = gc_check(s, base, gc_max_count);
		bool hp_checker = hp_check(s, base, hp_max);
		if (gc_checker == false || hp_checker == false){
			return false;
		}
		tm_add_base(&s->acc, prev_base, base);
	}

	return true;
}

uint64_t extend_seed(seed *s, int diff_k,int k, float gc_max, int hp_max){
	uint64_t curr_idx = s->idx;
	int gc_max_count = (int)((gc_max/100)*k);

	for(int i = 0; i < diff_k; i++){

		uint64_t base = rand() % ALPHABET_SIZE;
		uint64_t prev_base = s->last_base;
		bool quality_checker = false;

		while (quality_checker == false){
			bool hp_checker = hp_check(s, base, hp_max);
			bool gc_checker = gc_check(s, base, gc_max_count);

			if(gc_checker == false || hp_checker == false){
				base = (base + 1 + (rand() % (ALPHABET_SIZE - 1))) % ALPHABET_SIZE;
			} else {
				quality_checker = true;
			}
		}

		tm_add_base(&s->acc, prev_base, base);
		int current_k = k + i + 1;
		gc_max_count = (int)((gc_max/100)*current_k);
		curr_idx = (curr_idx << 2) | base;
	}

	return curr_idx;
}


