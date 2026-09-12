#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include <seed.h>

bool gc_check(seed *s, uint64_t base, int max_gc_count){
	if (base == 1 || base == 3){
		s->gc++;
		if (s->gc >= max_gc_count){
			s->gc--;
			return false;
		} else {
			return true;
		}
	}

	return true;
}

bool hp_check(seed *s, uint64_t base, int hp_max){
	if (s->last_base == base){
		s->hp_count++;
		if (s->hp_count >= hp_max){
			s->hp_count--;
			return false;
		}
	} else {
		s->hp_count = 1;
	}
	s->last_base = base;
	return true;
}

int hamming_check(uint64_t a, uint64_t b, int ext_k){
	uint64_t diff = a ^ b;
	uint64_t merged = (diff | (diff >> 1)) & 0x5555555555555555ULL;
	return __builtin_popcountll(merged);
}


