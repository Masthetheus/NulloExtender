#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define ALPHABET_SIZE 4 // A, C, T AND G
#define INITIAL_CAPACITY 100000

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

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

typedef struct seed{
	uint64_t idx;
	int gc;
	int hp_count;
	uint64_t last_base;
	TmAccumulator *acc;
} seed;

void decode_kmer(uint64_t idx, int k, char *seq) {
    for (int i = 0; i < k; i++) {
        uint64_t base = (idx >> ((k-i-1)*2)) & 3;
        seq[i] = BITS_TO_BASE[base];
    }
    seq[k] = '\0';
}

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

        char *seq = malloc((ext_k+1) * sizeof(char));
	decode_kmer(nullomers[random], ext_k, seq);
	printf("New seed:\n%s\n", seq);
	
	free(seq);
}

bool gc_check(seed *s, uint64_t base, int max_gc_count, int min_gc_count){
	if (base == 1 || base == 3){
		s->gc++;
		if (s->gc >= max_gc_count){
			s->gc--;
			return false;
		} else {
			return true;
		}
	}

	if (s->gc < min_gc_count){
		return false;
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
	return true;
}

bool seed_check(seed *s, int k, int gc_max_count, int gc_min_count, int hp_max){
	int gc = 0;
	int hp_count = 1;
	uint64_t idx = s->idx;
	
	uint64_t first_base = (idx >> ((k-1)*2)) & 3;
	
	if (first_base == 1 || first_base == 3){
		gc++;
	}

	for (int j = 1; j < k; j++){
		uint64_t base = (idx >> ((k-j-1)*2)) & 3;
		
		if (j == 1){
			s->last_base = base;
		}

		if (base == 1 || base == 3){
			gc++;
			if (gc > gc_max_count){
				return false;
			}
		}

		if (first_base == base){
			hp_count++;
			if (hp_count > hp_max){
				return false;
			}
		} else {
			hp_count = 1;
			first_base = base;
		}
		s->last_base = base;
	}
	
	if (gc < gc_min_count){
		return false;
	}
	
	uint64_t idx_holder = s->idx;
	s->gc = gc;
	s->hp_count = hp_count;
	s->last_base = idx_holder & 3;	
	return true;
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
		printf("Error: No valid numeric digits found.\n");
		return 1;
	}
	float gc_min = strtof(argv[3], &endptr);
	if (endptr == argv[3]) {
		printf("Error: No valid numeric digits found.\n");
		return 1;
	}
	int hp_max = atoi(argv[4]);
	int n = atoi(argv[5]);
        int ext_k = atoi(argv[6]);

        uint64_t *nullomers = NULL;
        size_t count = 0;
        size_t capacity = 0;


        fseek(f, 6, SEEK_SET);

        uint16_t parameters[4];
        fread(parameters,2,2, f);

        int k = parameters[0];
        int half_k = parameters[1];

        uint8_t codes[2];

        if (fread(codes, 1, 2, f) != 2){
                fprintf(stderr, "Error in file parsing.");
                fclose(f);
                exit(EXIT_FAILURE);
        }
        unsigned char byte_size = codes[0];
        unsigned char counter_size = codes[1];

        uint64_t v1 = 0;
        while (1) {
                v1 = 0;
                if (fread(&v1, byte_size, 1, f) != 1) break;

                uint64_t counter = 0;
                fread(&counter, counter_size, 1, f);
                uint64_t null_count = counter;

                uint8_t *null_collection = malloc(null_count * byte_size);
                if (fread(null_collection, byte_size, null_count, f) != null_count){
                        fprintf(stderr, "Error in file composition.\n");
                        fclose(f);
                        exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < null_count; i++) {
                        uint64_t v2 = 0;
                        memcpy(&v2, null_collection + i * byte_size, byte_size);
                        uint64_t idx = (v1 << ((uint64_t) (k - half_k) * 2)) | v2;
                        if (count == capacity){
                                size_t new_capacity = (capacity==0) ? INITIAL_CAPACITY : capacity * 2;
                                uint64_t *tmp = realloc(nullomers, new_capacity * sizeof(uint64_t));
                                if (!tmp) {
                                        perror("Realloc failed");
                                        free(nullomers);
                                        fclose(f);
                                        return 1;
                                }
                                nullomers = tmp;
                                capacity = new_capacity;
                        }
                        nullomers[count++] = idx;
                }
                free(null_collection);
        }

        seed seeds[n];

	int livecount = count;

        for (int i = 0; i < n; i++){
		select_seed(seeds, nullomers, i, &livecount, k);
        }

	for (int i = 0; i < n; i++){
		bool passed = true;
		int max_gc_count = (gc_max/100)*k;
		int min_gc_count = (gc_min/100)*k;

		bool seed_checker = seed_check(&seeds[i], k, max_gc_count, min_gc_count, hp_max);

		if (seed_checker == false){
			select_seed(seeds, nullomers, i, &livecount, k);
			i--;	
		}
	}

	uint64_t final_seq[n];

        for (int i = 0; i < n; i++){
                uint64_t curr_idx = seeds[i].idx;
		int max_gc_count = (gc_max/100)*k;
		int min_gc_count = (gc_min/100)*k;

                for(int j = 0; j < ext_k - k; j++){
                        uint64_t base = rand() % ALPHABET_SIZE;
			bool quality_checker = false;
	
			while (quality_checker == false){
				bool hp_checker = hp_check(&seeds[i], base, hp_max);
				bool gc_checker = gc_check(&seeds[i], base, max_gc_count, min_gc_count);

				if(gc_checker == false || hp_checker == false){
					base = (base + 1 + (rand() % (ALPHABET_SIZE - 1))) % ALPHABET_SIZE;
				} else {
					quality_checker = true;
				}
			}

			//tm_add_base(seeds[i]->acc, seeds[i].last_base, base);
			seeds[i].last_base = base;
			
			int current_k = k + j + 1;
			max_gc_count = (int)((gc_max/100)*current_k);
			min_gc_count = (int)((gc_min/100)*current_k);
                        curr_idx = (curr_idx << 2) | base;
                }
                char *seq = malloc((ext_k+1) * sizeof(char));
                decode_kmer(curr_idx, ext_k, seq);
                printf("%s\n", seq);
                
                free(seq);
        }
               

        fclose(f);
        free(nullomers);

        return 0;
}
