#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "seed.h"
#include "seed_operations.h"

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
}


