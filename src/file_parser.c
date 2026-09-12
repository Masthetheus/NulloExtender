#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "file_parser.h"

#define INITIAL_CAPACITY 100000

uint64_t *parse_nullomer_file(FILE *f, int *k, size_t *count){
	uint64_t *nullomers = NULL;
	size_t capacity = 0;

	// skips basic header
	fseek(f, 6, SEEK_SET);

	// obtain k and half_k values
	uint16_t parameters[2];
	fread(parameters, 2, 2, f);
	*k = parameters[0];
	int half_k = parameters[1];

	// reads bit encoding of file
	uint8_t codes[2];	
	
	if (fread(codes, 1, 2, f) != 2){
		fprintf(stderr, "Error in file parsing.");
		fclose(f);
		exit(EXIT_FAILURE);
	}

	unsigned char byte_size = codes[0];
	unsigned char counter_size = codes[1];

	// starts reading the v1 and v2 groups
	uint64_t v1;
	while(fread(&v1, byte_size, 1, f) == 1) {

		// obtain how many v2s for current v1
                uint64_t counter = 0;
                if (fread(&counter, counter_size, 1, f) != 1) break;
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

                        uint64_t idx = (v1 << ((uint64_t) (*k - half_k) * 2)) | v2;

                        if (*count == capacity){
                                size_t new_capacity = (capacity==0) ? INITIAL_CAPACITY : capacity * 2;
                                uint64_t *tmp = realloc(nullomers, new_capacity * sizeof(uint64_t));
                                if (!tmp) {
                                        perror("Realloc failed");
                                        free(nullomers);
                                        fclose(f);
					exit(EXIT_FAILURE);
                                }
                                nullomers = tmp;
                                capacity = new_capacity;
                        }

                        nullomers[*count] = idx;
			(*count)++;
                }
                free(null_collection);
        }
	return nullomers;
}
