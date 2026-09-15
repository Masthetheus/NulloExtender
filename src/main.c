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
#include "nn_table.h"

#define ALPHABET_SIZE 4 // A, C, T AND G
#define INITIAL_CAPACITY 100000
#define NA_CONC 0.05

int main(int argc, char *argv[]){
        srand(time(NULL));
        if (argc != 9){
                fprintf(stderr, "Usage: %s <nullomer_file> <gc_max> <gc_min> <homopolymer_max> <number_of_seeds> <target_length> <min_tm> <max_tm>\n", argv[0]);
                return 1;
        }

        FILE *f = fopen(argv[1], "rb");
        if (!f){
                perror("Error opening file!");
                return 1;
        }
	
	// loading of filtering parameters
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
	double tm_min = strtof(argv[7], &endptr);
	double tm_max = strtof(argv[8], &endptr);
	if (tm_min >= tm_max){
	    fprintf(stderr, "Error: Min Tm must be smaller than Max Tm!");
	    return 1;
	}
		
	// gather null set from file
	int k = 0;
	size_t count = 0;
	uint64_t *nullomers = parse_nullomer_file(f, &k, &count);


	// obtain first draft of seeds
        seed seeds[n];
        initialize_seeds(seeds, ext_k, n);
	int livecount = count;

        for (int i = 0; i < n; i++){
		select_seed(seeds, nullomers, i, &livecount, k);
        }

	// check if initial seeds passes filtering, selecting new ones if needed
	for (int i = 0; i < n; i++){
		int max_gc_count = (gc_max/100)*k;

		bool seed_checker = seed_check(&seeds[i], k, max_gc_count, hp_max);

		if (seed_checker == false){
			select_seed(seeds, nullomers, i, &livecount, k);
			i--;	
		}
	}

	// extend all approved seeds
	// making sure they are still on par with the filters
	seed final_seq[n];
	int diff_k = ext_k - k;
	int max_gc_count = (gc_max/100)*k;
        for (int i = 0; i < n; i++){
		int min_gc_count = (gc_min/100)*ext_k;

		extend_seed(&seeds[i], diff_k, k, gc_max, hp_max);	
                uint8_t first_base = seed_get_base(&seeds[i], 0);
		double final_tm = tm_finalize(&seeds[i].acc, first_base, seeds[i].last_base, NA_CONC, ext_k);
		if (seeds[i].gc < min_gc_count || final_tm < tm_min || final_tm > tm_max){
			select_seed(seeds, nullomers, i, &livecount, k);
			seed_check(&seeds[i], k, max_gc_count, hp_max);
			i--;
			continue;
		}	
        }

	// check orthogonality, randomizing new seeds
	// when it fails the minimum hamming distance
        for (int i = 0; i < n; i++){
                for (int j = i+1; j < n; j++){
                        int hamming = hamming_check(&seeds[i], &seeds[j]);
                        if (hamming < 3){
                                select_seed(seeds, nullomers, i, &livecount, k);
                                seed_check(&seeds[i], k, max_gc_count, hp_max);
                                extend_seed(&seeds[i], diff_k, k, gc_max, hp_max);
                                i--;
                                break;
                        }
                }
        }

	// output final primers
	for (int i = 0; i < n; i++){
                char *seq = malloc((ext_k+1) * sizeof(char));
                decode_kmer(&seeds[i], ext_k, seq);
                printf("%s\n", seq);
                free(seq);
	}
               

        seed_destroy(seeds, n);
        fclose(f);
        free(nullomers);

        return 0;
}
