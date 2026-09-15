#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "seed.h"
#include "nn_table.h"

bool gc_check(seed *s, uint8_t base, int max_gc_count){
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

bool hp_check(seed *s, uint8_t base, int hp_max){
	int candidate_hp_count;

	if (s->last_base == base){
		candidate_hp_count = s->hp_count + 1;
		if (candidate_hp_count >= hp_max){
			return false;
		}
	} else {
		candidate_hp_count = 1;
	}

	s->hp_count = candidate_hp_count;
	s->last_base = base;
	return true;
}

int hamming_check(const seed *a, const seed *b){
    int dist = 0;
    for (int i = 0; i < a->n_blocks; i++){
        uint64_t diff = a->blocks[i] ^ b->blocks[i];

        if (i == a->n_blocks - 1){
            int used_bits = (a->length % 32 == 0) ? 64 : (a->length % 32) * 2;
            uint64_t mask = (used_bits == 64) ? ~0ULL : ((1ULL << used_bits) - 1);
            diff &= mask;
        }

        uint64_t merged = (diff | (diff >> 1)) & 0x5555555555555555ULL;
        dist += __builtin_popcountll(merged);
    }
    return dist;
}

void tm_add_base(TmAccumulator *acc, uint8_t prev_base, uint8_t curr_base) {
    acc->delta_h_acc += nn_table[prev_base][curr_base].delta_h;
    acc->delta_s_acc += nn_table[prev_base][curr_base].delta_s;
}

double tm_finalize(TmAccumulator *acc, uint8_t first_base, uint8_t last_base,
                    double na_conc, int seq_len) {
    double h = acc->delta_h_acc;
    double s = acc->delta_s_acc;

    // terminal initiation correction
    if (first_base == 1 || first_base == 3) { // C or G
        h += init_params.gc_init_h;
        s += init_params.gc_init_s;
    } else { // A or T
        h += init_params.at_init_h;
        s += init_params.at_init_s;
    }

    if (last_base == 1 || last_base == 3) {
        h += init_params.gc_init_h;
        s += init_params.gc_init_s;
    } else {
        h += init_params.at_init_h;
        s += init_params.at_init_s;
    }

    // salt correction
    s += 0.368 * (seq_len - 1) * log(na_conc);

    const double R = 1.987;
    double Ct = 250e-9;
    double tm_kelvin = (h * 1000.0) / (s + R * log(Ct / 4.0));
    return tm_kelvin - 273.15;
}

