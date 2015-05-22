#ifndef VECTOR_H_
#define VECTOR_H_

#include <stdint.h>

/* utility functions */

int64_t fast_rand(void);

void set_seed(int64_t value);
void set_length(int64_t value);
void vecInfo(int64_t field, int64_t data);

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
