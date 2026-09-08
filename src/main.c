#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#define ALPHABET_SIZE 4 // A, C, T AND G
#define INITIAL_CAPACITY 100000

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

typedef struct seed{
	uint64_t idx;
	int gc;
	int hp_count;
	uint64_t last_base;
} seed;

void decode_kmer(uint64_t idx, int k, char *seq) {
    for (int i = 0; i < k; i++) {
        uint64_t base = (idx >> ((k-i-1)*2)) & 3;
        seq[i] = BITS_TO_BASE[base];
    }
    seq[k] = '\0';
}

void select_new_seed(seed *seeds, uint64_t *nullomers, int n, int count, int ext_k){
	uint64_t random = rand() % count;
	seeds[n].idx = nullomers[random];

        char *seq = malloc((ext_k+1) * sizeof(char));
	decode_kmer(nullomers[random], ext_k, seq);
	printf("New seed:\n%s\n", seq);
	
	free(seq);
}

int main(int argc, char *argv[]){
        srand(time(NULL));
        if (argc != 5){
                fprintf(stderr, "Usage: %s <nullomer_file> <gc_max> <homopolymer_max> <number_of_seeds>\n", argv[0]);
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
	int homopolymer_max = atoi(argv[3]);
	int n = atoi(argv[4]);

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

        int ext_k = k + 5;
	int livecount = count;

        for (int i = 0; i < n; i++){
                uint64_t random = rand() % count;
                seeds[i].idx = nullomers[random];

		uint64_t holder = nullomers[random];
		nullomers[random] = nullomers[livecount - 1];
		nullomers[livecount - 1] = holder;
		livecount--;

                char *seq = malloc((k+1) * sizeof(char));

                decode_kmer(seeds[i].idx, k, seq);
                printf("%s\n", seq);

                free(seq);
        }

        
	//select_new_seed(seeds, nullomers, 1, count, ext_k);

	uint64_t hp_max = 3;

	for (int i = 0; i < n; i++){
		bool passed = true;
		int max_gc_count = (gc_max/100)*k;
		uint64_t idx = seeds[i].idx;
		uint64_t first_base = (idx >> ((k-1)*2)) & 3;
		int hp_count = 1;
		int gc = 0;

		if (first_base == 1 || first_base == 3){
			 gc++;
		}
		for (int j = 1; j < k; j++){
			uint64_t base = (idx >> ((k-j-1)*2)) & 3;
			if (base == 1 || base == 3){
				gc++;
				if (gc > max_gc_count){
					passed = false;
				}
			}
			if (first_base == base){
				hp_count++;
				if (hp_count > hp_max){
					passed = false;
				}
			} else{
				hp_count = 1;
				first_base = base;
			}
			if (passed == false){
				select_new_seed(seeds, nullomers, i, livecount, k);
				i--;	
				break;
			}
		}
		if (passed == true){
			uint64_t idx_holder = seeds[i].idx;
			seeds[i].gc = gc;
			seeds[i].hp_count = hp_count;
			seeds[i].last_base = idx_holder & 3;
		}
	}

	uint64_t final_seq[n];

        for (int i = 0; i < n; i++){
                uint64_t curr_idx = seeds[i].idx;
		int max_gc_count = (gc_max/100)*k;

                for(int j = 0; j < ext_k - k; j++){
                        uint64_t base = rand() % ALPHABET_SIZE;
			if (base == 1 || base == 3){
				seeds[i].gc++;
				if (seeds[i].gc >= max_gc_count){
					seeds[i].gc--;
					base--;
				}
			}
			if (seeds[i].last_base == base){
				seeds[i].hp_count++;
				if (seeds[i].hp_count >= hp_max){
					seeds[i].hp_count--;
					bool flag = true;
					while (flag){
						uint64_t new_base = rand() % ALPHABET_SIZE;
						if (new_base != base){
							base = new_base;
							flag = false;
						}
					}
				}
			}
			seeds[i].last_base = base;

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
