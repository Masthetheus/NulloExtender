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

void initialize_seeds(seed *seeds, int ext_k, int n) {
        int n_blocks = (ext_k + 31)/32;
        for (int i = 0; i < n; i++){
                seeds[i].blocks = calloc(n_blocks, sizeof(uint64_t));
                seeds[i].length = 0;
                seeds[i].n_blocks = n_blocks;
        }
}

void seed_push_base(seed *s, uint64_t base){
        int block_idx = s->length/32;
        int offset = (s->length % 32) * 2;
        s->blocks[block_idx] |= (base << offset);
        s->length++;
}

uint8_t seed_get_base(const seed *s, int pos){
        int block_idx = pos/32;
        int offset = (pos%32) * 2;
        return (s->blocks[block_idx] >> offset) & 3;
}

void seed_destroy(seed *seeds, int n){
        for (int i = 0; i < n; i++){
                free(seeds[i].blocks);
        }
}

void select_seed(seed *seeds, uint64_t *nullomers, int i, int *count, int k){
        if (*count == 0){
                fprintf(stderr, "No more suitable seeds available.\n");
                exit(EXIT_FAILURE);
        }
	uint64_t random = rand() % *count;
	if (nullomers[random]){
                seeds[i].length = 0; 
                memset(seeds[i].blocks, 0, seeds[i].n_blocks * sizeof(uint64_t));
                for (int j = 0; j < k; j++){
                        uint64_t base = (nullomers[random] >> ((k-1-j)*2)) & 3;
                        seed_push_base(&seeds[i], base);
                }
		seeds[i].gc = 0;
		seeds[i].hp_count = 1;
		seeds[i].last_base = seed_get_base(&seeds[i],0);
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

void seed_debug_print(const seed *s, int k) {
    fprintf(stderr, "length=%d gc=%d hp_count=%d last_base=%d | ", s->length, s->gc, s->hp_count, (int)s->last_base);
    for (int i = 0; i < s->length; i++) {
        fprintf(stderr, "%d", (int)seed_get_base(s, i));
    }
    fprintf(stderr, "\n");
}

bool seed_check(seed *s, int k, int gc_max_count, int hp_max){
        uint8_t first_base = seed_get_base(s,0);
	
	if (first_base == 1 || first_base == 3){
		s->gc++;
	}	

	for (int i = 1; i < k; i++){
		uint64_t base = seed_get_base(s, i);
		uint64_t prev_base = s->last_base;
		bool gc_checker = gc_check(s, base, gc_max_count);
		bool hp_checker = hp_check(s, base, hp_max);
		if (gc_checker == false || hp_checker == false){
			return false;
		}
		tm_add_base(&s->acc, prev_base, base);
	}

        //seed_debug_print(s, k);

	return true;
}

void extend_seed(seed *s, int diff_k,int k, float gc_max, int hp_max){
	uint64_t curr_idx = s->blocks[0];
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
                seed_push_base(s, base);
	}
}


