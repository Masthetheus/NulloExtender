#ifndef NN_TABLE_H
#define NN_TABLE_H

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

#endif
