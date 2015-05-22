#ifndef SCALER_H_
#define SCALER_H_


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <inttypes.h>
#include <math.h>
/* utility functions */

void error(void);
void release(void);

/* input and output */

void define_settings();
void define_vectors();
void compute_engine();


// #include "vector.h"
typedef struct  {
    int64_t sum;
    int sumFlag;
    int64_t mode;
    int modeFlag;
    int64_t median;
    int medianFlag;
    int64_t minimum;
    int minFlag;
    int64_t maximum;
    int maxFlag;
} props;
/* utility functions */

int64_t fast_rand(void);

void set_seed(int64_t value);
void set_length(int64_t value);

// void * prelimSieveofEratosthnes (int64_t * primes);
// void * printlulz(int64_t * vector);

/* vector operations */

int64_t* new_vector(void);
int64_t* prime_vector(int64_t start);
int64_t* random_vector(int64_t seed);
int64_t* uniform_vector(int64_t value);
int64_t* sequence_vector(int64_t start, int64_t step);

int64_t* cloned(int64_t* vector);
int64_t* reversed(int64_t* vector);
int64_t* ascending(int64_t* vector);
int64_t* descending(int64_t* vector);

int64_t* scalar_add(int64_t* vector, int64_t scalar);
int64_t* scalar_mul(int64_t* vector, int64_t scalar);

int64_t* vector_add(int64_t* vector_a, int64_t* vector_b);
int64_t* vector_mul(int64_t* vector_a, int64_t* vector_b);

/* compute operations */

int64_t get_sum(int64_t* vector);
int64_t get_mode(int64_t* vector);
int64_t get_median(int64_t* vector);
int64_t get_minimum(int64_t* vector);
int64_t get_maximum(int64_t* vector);
int64_t get_element(int64_t* vector, int64_t index);
int64_t get_frequency(int64_t* vector, int64_t value);

void display(int64_t* vector, int64_t label);

#endif
