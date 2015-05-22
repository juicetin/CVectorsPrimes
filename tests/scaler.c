////////////////////////////////
///      CORRECT CONTROL     ///
////////////////////////////////

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <inttypes.h>

#include "scaler.h"
#include "vector.h"

#define UNSET -1
#define MAX_BUFFER 64

#define OUTPUT_INT32_T "%" PRId32 "\n"
#define OUTPUT_INT64_T "%" PRId64 "\n"

static int64_t g_ncores        = 0; /* 1 <= n <= 128 */
static int64_t g_nthreads      = 0; /* 1 <= n <= 256 */

static int64_t g_length        = 0; /* 1 <= n <= 1,000,000,000 */
static int64_t g_nstored       = 0; /* 1 <= n <= 1,000,000 */
static int64_t g_nvectors      = 0; /* 1 <= n <= 1,000,000 */
static int64_t g_ncomputations = 0; /* 1 <= n <= 1,000,000 */

static int64_t** g_vectors     = NULL;

/**
 * Release all dynamically allocated memory
 */
void release(void) {

    if (g_vectors == NULL) {
        return;
    }

    /* free any vectors that have been allocated */

    for (int64_t i = 0; i < g_nvectors; i++) {
        if (g_vectors[i] != NULL) {
            free(g_vectors[i]);
        }
    }

    free(g_vectors);
}

/**
 * Outputs error and terminates program
 */
void error(void) {

    puts("error");
    release();
    exit(0);
}

/**
 * Defines setting based on given input
 */
void define_settings(void) {

    int parsed_args = 0;
    char input[MAX_BUFFER];

    /* <# cores> <# threads> */

    fgets(input, MAX_BUFFER, stdin);
    parsed_args = sscanf(input, "%" PRId64 " %" PRId64, &g_ncores, &g_nthreads);

    if (parsed_args != 2) {
        error();
    }

    /* <vector length> <# vectors> <# computations> */

    fgets(input, MAX_BUFFER, stdin);
    parsed_args = sscanf(input, "%" PRId64 " %" PRId64 " %" PRId64, &g_length, &g_nvectors, &g_ncomputations);

    if (parsed_args != 3) {
        error();
    }

    set_length(g_length);
}

/**
 * Defines vectors based on given input
 */
void define_vectors(void) {

    g_vectors = (int64_t **) calloc(g_nvectors, sizeof(int64_t *));

    for (int64_t i = 0; i < g_nvectors; i++) {

        char input[MAX_BUFFER];
        char function[MAX_BUFFER];

        int64_t arg1 = UNSET;
        int64_t arg2 = UNSET;

        /* vector= <function> [args] */

        fgets(input, MAX_BUFFER, stdin);
        sscanf(input, "vector= %s %" PRId64 " %" PRId64, function, &arg1, &arg2);

        if (strcmp(function, "primes") == 0) {
            g_vectors[g_nstored++] = prime_vector(arg1);
        } else if (strcmp(function, "random") == 0) {
            g_vectors[g_nstored++] = random_vector(arg1);
        } else if (strcmp(function, "uniform") == 0) {
            g_vectors[g_nstored++] = uniform_vector(arg1);
        } else if (strcmp(function, "sequence") == 0) {
            g_vectors[g_nstored++] = sequence_vector(arg1, arg2);
        } else if (strcmp(function, "cloned") == 0) {
            g_vectors[g_nstored++] = cloned(g_vectors[arg1]);
        } else if (strcmp(function, "reversed") == 0) {
            g_vectors[g_nstored++] = reversed(g_vectors[arg1]);
        } else if (strcmp(function, "ascending") == 0) {
            g_vectors[g_nstored++] = ascending(g_vectors[arg1]);
        } else if (strcmp(function, "descending") == 0) {
            g_vectors[g_nstored++] = descending(g_vectors[arg1]);
        } else if (strcmp(function, "scalar#add") == 0) {
            g_vectors[g_nstored++] = scalar_add(g_vectors[arg1], arg2);
        } else if (strcmp(function, "scalar#mul") == 0) {
            g_vectors[g_nstored++] = scalar_mul(g_vectors[arg1], arg2);
        } else if (strcmp(function, "vector#add") == 0) {
            g_vectors[g_nstored++] = vector_add(g_vectors[arg1], g_vectors[arg2]);
        } else if (strcmp(function, "vector#mul") == 0) {
            g_vectors[g_nstored++] = vector_mul(g_vectors[arg1], g_vectors[arg2]);
        } else {
            error();
        }
    }
}

/**
 * Runs computations based on given input
 */
void compute_engine(void) {

    for (int64_t i = 0; i < g_ncomputations; i++) {

        char input[MAX_BUFFER];
        char function[MAX_BUFFER];

        int64_t arg1 = UNSET;
        int64_t arg2 = UNSET;

        /* compute <function> [args] */

        fgets(input, MAX_BUFFER, stdin);
        sscanf(input, "compute %s %" PRId64 " %" PRId64, function, &arg1, &arg2);

        if (strcmp(function, "display") == 0) {
            display(g_vectors[arg1], arg1);
        } else {

            if (arg2 == UNSET) {
                printf("%s of vector~%" PRId64 " = ", function, arg1);
            }

            if (strcmp(function, "sum") == 0) {
                printf(OUTPUT_INT64_T, get_sum(g_vectors[arg1]));
            } else if (strcmp(function, "mode") == 0) {
                printf(OUTPUT_INT64_T, get_mode(g_vectors[arg1]));
            } else if (strcmp(function, "median") == 0) {
                printf(OUTPUT_INT64_T, get_median(g_vectors[arg1]));
            } else if (strcmp(function, "minimum") == 0) {
                printf(OUTPUT_INT64_T, get_minimum(g_vectors[arg1]));
            } else if (strcmp(function, "maximum") == 0) {
                printf(OUTPUT_INT64_T, get_maximum(g_vectors[arg1]));
            } else if (strcmp(function, "element") == 0) {
                printf("%s %" PRId64 " in vector~%" PRId64 " = ", function, arg2, arg1);
                printf(OUTPUT_INT64_T, get_element(g_vectors[arg1], arg2));
            } else if (strcmp(function, "frequency") == 0) {
                printf("%s of %" PRId64 " in vector~%" PRId64 " = ", function, arg2, arg1);
                printf(OUTPUT_INT64_T, get_frequency(g_vectors[arg1], arg2));
            } else {
                error();
            }
        }
    }
}

/**
 * Main function
 */
int main(void)
{
    define_settings();
    define_vectors();
    compute_engine();

    release();

    return 0;
}
