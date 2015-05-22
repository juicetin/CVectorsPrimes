#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>

#include "vector.h"

static int64_t g_seed = 0;
static int64_t g_length = 0;

////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////

int int64Ascend(const void *x, const void *y)
{
    const int64_t * a = (int64_t *) x;
    const int64_t * b = (int64_t *) y;

    if (*a > *b) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
    // return (*(int64_t *)x > *(int64_t *)y) - (*(int64_t *)x < *(int64_t *)y);
}

//Compare function for two int64s to sort in descending order
int int64Descend(const void *x, const void *y)
{
    const int64_t * a = (int64_t *) x;
    const int64_t * b = (int64_t *) y;

    if (*a < *b) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
    // return (*(int64_t *)x < *(int64_t *)y) - (*(int64_t *)x > *(int64_t *)y);
}

void swap(int64_t * vector, int64_t a, int64_t b)
{
    int64_t tmp = vector[a];
    vector[a] = vector[b];
    vector[b] = tmp;
}

int64_t quickSelect (int64_t * arr, int64_t n, int64_t k)
{
    // clock_t tick, tock;
    // tick = clock();
    int64_t i, ir, j, l, mid, a;
    l = 0;
    ir = n - 1;
    for(;;) {
        if (ir <= l+1) { 
            if (ir == l+1 && arr[ir] < arr[l]) {
                swap(arr, l, ir);
            }
            // tock = clock();
            // printf("Time elapsed: %f seconds\n", (double)(tock-tick)/CLOCKS_PER_SEC);
            return arr[k];
        }
        else 
        {
            // mid=(l+ir) >> 1;
            if ((l+ir) % 2 == 0)
                mid=((l+ir)-1) >> 1;
            else 
                mid=(l+ir) >> 1; 
            swap(arr, mid, l+1);
            if (arr[l] > arr[ir])
                swap(arr, l, ir);
            if (arr[l+1] > arr[ir]) 
                swap(arr, l+1, ir);
            if (arr[l] > arr[l+1]) 
                swap(arr, l, l+1);
            i=l+1; 
            j=ir;
            a=arr[l+1]; 
            for (;;) 
            { 
                do i++; while (arr[i] < a); 
                do j--; while (arr[j] > a); 
                if (j < i) break; 
                swap(arr, i, j);
            } 
            arr[l+1]=arr[j]; 
            arr[j]=a;
            if (j >= k) ir=j-1; 
            if (j <= k) l=i;
        }
    }
}

/**
 * Sets the number of elements that each vector will contain
 */
 void set_length(int64_t length) {

    g_length = length;
}

/**
 * Sets the seed used when generating pseudorandom numbers
 */
 void set_seed(int64_t seed) {

    g_seed = seed;
}

/**
 * Returns pseudorandom number determined by the seed
 */
 int64_t fast_rand(void) {

    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

////////////////////////////////
///   VECTOR INITALISATIONS  ///
////////////////////////////////

/**
 * Returns new vector, with all elements set to zero
 */
 int64_t* new_vector(void) {

    return (int64_t *) calloc(g_length, sizeof(int64_t));
}

/**
 * Returns new vector, with all elements set to given value
 */
 int64_t* uniform_vector(int64_t value) {

    int64_t* vector = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        vector[i] = value;
    }

    return vector;
}

/**
 * Returns new vector, with elements generated at random using given seed
 */
 int64_t* random_vector(int64_t seed) {

    int64_t* vector = new_vector();

    set_seed(seed);

    for (int64_t i = 0; i < g_length; i++) {
        vector[i] = fast_rand();
    }

    return vector;
}

/**
 * Returns whether given number is prime
 */
 bool is_prime(int64_t number) {

    if (number == 0) {
        return false;
    } else if (number == 1) {
        return false;
    }

    for (int64_t i = 2; i < number; i++) {
        if (number % i == 0) {
            return false;
        }
    }

    return true;
}

/**
 * Returns new vector, containing primes numbers in sequence from given start
 */
 int64_t* prime_vector(int64_t start) {

    int64_t* vector = new_vector();

    int64_t number = start;
    int64_t i = 0;

    while (i < g_length) {

        if (is_prime(number)) {
            vector[i] = number;
            i += 1;
        }

        number += 1;
    }

    return vector;
}

/**
 * Returns new vector, with elements in sequence from given start and step
 */
 int64_t* sequence_vector(int64_t start, int64_t step) {

    int64_t* vector = new_vector();

    int64_t current = start;

    for (int64_t i = 0; i < g_length; i++) {
        vector[i] = current;
        current += step;
    }

    return vector;
}

////////////////////////////////
///     VECTOR OPERATIONS    ///
////////////////////////////////

int64_t* cloned(int64_t* vector) {

    int64_t* clone = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        clone[i] = vector[i];
    }

    return clone;
}

int64_t* reversed(int64_t* vector) {

    int64_t* result = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[g_length - 1 - i];
    }

    return result;
}


int64_t* ascending(int64_t* vector) {
    int64_t* result = cloned(vector);
    qsort(result, g_length, sizeof(int64_t), int64Ascend);
    return result;
}

int64_t* descending(int64_t* vector) {
    int64_t* result = cloned(vector);
    qsort(result, g_length, sizeof(int64_t), int64Descend);
    return result;
}

int64_t* scalar_add(int64_t* vector, int64_t scalar) {

    int64_t* result = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[i] + scalar;
    }

    return result;
}

int64_t* scalar_mul(int64_t* vector, int64_t scalar) {

    int64_t* result = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[i] * scalar;
    }

    return result;
}

int64_t* vector_add(int64_t* vector_a, int64_t* vector_b) {

    int64_t* result = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector_a[i] + vector_b[i];
    }

    return result;
}

int64_t* vector_mul(int64_t* vector_a, int64_t* vector_b) {

    int64_t* result = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector_a[i] * vector_b[i];
    }

    return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

/**
 * Returns the sum of all elements
 */
 int64_t get_sum(int64_t* vector) {

    int64_t sum = 0;

    for (int64_t i = 0; i < g_length; i++) {
        sum += vector[i];
    }

    return sum;
}


// * Returns the most frequently occuring element
// * or -1 if there is no such unique element
int64_t get_mode(int64_t* vector) {
    int64_t* result = cloned(vector);
    qsort(result, g_length, sizeof(int64_t), int64Ascend);

    int64_t curSpanCount = 0, hiSpanCount = 0, 
    currentCheck = 0, repeatCount = 0, mode = 0;

    for (int64_t i = 0; i < g_length; i++)
    {
        if (currentCheck != result[i])
        {
            currentCheck = result[i];
            curSpanCount = 0;
        }
        curSpanCount++;

        if (curSpanCount > hiSpanCount)
        {
            repeatCount = 0;
            hiSpanCount = curSpanCount;
            mode = currentCheck;   
        }
        else if (currentCheck != result[i+1] && curSpanCount == hiSpanCount)
        {
            repeatCount = 1;
        }
    }

    if (repeatCount > 0)
        return -1;
    return mode;
}

 // Returns the lower median

int64_t get_median(int64_t* vector) {

    int64_t n = g_length/2;
    if (g_length % 2 == 0)
        n = g_length/2 - 1;

    int64_t * result = vector;
    return quickSelect (result, g_length, n);
}

/**
 * Returns the smallest value in the vector
 */
 int64_t get_minimum(int64_t* vector) {

    int64_t min = INT64_MAX;
    for (int64_t = 0; i < g_length; i++) {
        if (vector[i] < min)
            min = vector[i];
    }
    return min;
}

/**
 * Returns the largest value in the vector
 */
 int64_t get_maximum(int64_t* vector) {

    int64_t max = 0;
    for (int64_t i = 0; i < g_length; i++) {
        if (max < vector[i])
            max = vector[i];
    }
    return max;
}

/**
 * Returns the frequency of the value in the vector
 */
 int64_t get_frequency(int64_t* vector, int64_t value) {

    int64_t count = 0;

    for (int64_t i = 0; i < g_length; i++) {
        if (vector[i] == value) {
            count += 1;
        }
    }

    return count;
}

/**
 * Returns the value stored at the given element index
 */
 int64_t get_element(int64_t* vector, int64_t index) {

    return vector[index];
}

/**
 * Output given vector to standard output
 */
 void display(int64_t* vector, int64_t label) {

    printf("vector~%" PRId64 " ::", label);

    for (int64_t i = 0; i < g_length; i++) {
        printf(" %" PRId64, vector[i]);
    }

    puts("");
}

