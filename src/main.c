#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "file_parser.h"
#include "seed.h"
#include "seed_operations.h"
#include "filters.h"
#include "decode.h"

#define ALPHABET_SIZE 4 // A, C, T AND G
#define INITIAL_CAPACITY 100000


typedef struct {
	double delta_h;
	double delta_s;
} NNParams;

static const NNParams nn_table[4][4] = {
    /*            A              C            T             G        */
    /* A */ {{-7.9,-22.2}, {-8.4,-22.4}, {-7.2,-20.4}, {-8.2,-22.2}},
    /* C */ {{-8.5,-22.7}, {-8.0,-19.9}, {-7.8,-21.0}, {-9.8,-24.4}},
    /* T */ {{-7.2,-21.3}, {-8.2,-22.2}, {-7.9,-22.2}, {-8.4,-22.4}},
    /* G */ {{-8.2,-22.2}, {-9.8,-24.4}, {-8.5,-22.7}, {-8.0,-19.9}}
};

typedef struct {
	double gc_init_h, gc_init_s;
	double at_init_h, at_init_s;
} InitParams;

static const InitParams init_params = {
    .gc_init_h = 0.1,  .gc_init_s = -2.8,
    .at_init_h = 2.3,  .at_init_s = 4.1
};

typedef struct {
    double delta_h_acc;
    double delta_s_acc;
    int initialized;
} TmAccumulator;

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

int main(int argc, char *argv[]){
        srand(time(NULL));
        if (argc != 7){
                fprintf(stderr, "Usage: %s <nullomer_file> <gc_max> <gc_min> <homopolymer_max> <number_of_seeds> <target_length>\n", argv[0]);
                return 1;
        }

        FILE *f = fopen(argv[1], "rb");
        if (!f){
                perror("Error opening file!");
                return 1;
        }

	char *endptr;
	float gc_max = strtof(argv[2], &endptr);
	if (endptr == argv[2]) {
		fprintf(stderr, "Error: No valid numeric digits found.\n");
		return 1;
	}
	float gc_min = strtof(argv[3], &endptr);
	if (endptr == argv[3]) {
		fprintf(stderr, "Error: No valid numeric digits found.\n");
		return 1;
	}

	if (gc_min >= gc_max){
		fprintf(stderr, "Error: Min GC must be smaller than Max GC!");
		return 1;
	}	

	int hp_max = atoi(argv[4]);
	int n = atoi(argv[5]);
        int ext_k = atoi(argv[6]);
		
	int k = 0;
	size_t count = 0;
	uint64_t *nullomers = parse_nullomer_file(f, &k, &count);

        seed seeds[n];

	int livecount = count;

        for (int i = 0; i < n; i++){
		select_seed(seeds, nullomers, i, &livecount, k);
        }

	for (int i = 0; i < n; i++){
		bool passed = true;
		int max_gc_count = (gc_max/100)*k;
		int min_gc_count = (gc_min/100)*k;

		bool seed_checker = seed_check(&seeds[i], k, max_gc_count, hp_max);

		if (seed_checker == false){
			select_seed(seeds, nullomers, i, &livecount, k);
			i--;	
		}
	}

	printf("All seeds checked\n");
	uint64_t final_seq[n];
	int diff_k = ext_k - k;

	int max_gc_count = (gc_max/100)*k;
        for (int i = 0; i < n; i++){
                uint64_t curr_idx = seeds[i].idx;
		int min_gc_count = (gc_min/100)*ext_k;
		curr_idx = extend_seed(&seeds[i], diff_k, k, gc_max, hp_max);	
		if (seeds[i].gc < min_gc_count){
			select_seed(seeds, nullomers, i, &livecount, k);
			seed_check(&seeds[i], k, max_gc_count, hp_max);
			i--;
			continue;
		}	
		final_seq[i] = curr_idx;
                char *seq = malloc((ext_k+1) * sizeof(char));
                decode_kmer(curr_idx, ext_k, seq);
                printf("%s\n", seq);
                
                free(seq);
        }

	for (int i = 0; i < n; i++){
		uint64_t a = final_seq[i];
		for (int j = i+1; j < n; j++){
			uint64_t b = final_seq[j];
			int hamming = hamming_check(a,b,ext_k);
			if (hamming < 3){
				select_seed(seeds, nullomers, i, &livecount, k);
				seed_check(&seeds[i], k, max_gc_count, hp_max);
				uint64_t curr_idx = extend_seed(&seeds[i], diff_k, k, gc_max, hp_max);
				final_seq[i] = curr_idx;
				i--;
				break;
			}
		}
	}

	for (int i = 0; i < n; i++){
		uint64_t curr_idx = final_seq[i];
                char *seq = malloc((ext_k+1) * sizeof(char));
                decode_kmer(curr_idx, ext_k, seq);
                printf("%s\n", seq);
                
                free(seq);

	}
               

        fclose(f);
        free(nullomers);

        return 0;
}
