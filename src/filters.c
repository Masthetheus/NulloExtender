#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "seed.h"
#include "nn_table.h"

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

void tm_add_base(TmAccumulator *acc, uint8_t prev_base, uint8_t curr_base) {
    acc->delta_h_acc += nn_table[prev_base][curr_base].delta_h;
    acc->delta_s_acc += nn_table[prev_base][curr_base].delta_s;
}

double tm_finalize(TmAccumulator *acc, uint8_t first_base, uint8_t last_base,
                    double na_conc, int seq_len) {
    double h = acc->delta_h_acc;
    double s = acc->delta_s_acc;

    // salt correction
    s += 0.368 * (seq_len - 1) * log(na_conc);

    const double R = 1.987; // cal/(mol*K)
    double Ct = 250e-9;     // 250nM, default
    double tm_kelvin = (h * 1000.0) / (s + R * log(Ct / 4.0));
    return tm_kelvin - 273.15;
}


