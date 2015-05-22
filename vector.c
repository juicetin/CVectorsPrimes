/*Sublime text is recommended to read this file. The task was to implement a number of mathematical operations in C
to run as quickly as possible using algorithms and optimisation. My ranking was 8th in a cohort of 300 for this 
competitive assignment.*/

#include "vector.h"
#include "immintrin.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

/*
 * To change pre-sieve values, simply decide on a number
 * of pre-generated primes, e.g. 150,000 - the 150000th
 * prime is 2,015,177, and it's root is 1419.57.
 * 150000 primes, cap 2015177, root 1420
 * 200000 primes, cap 2750160, root 1659
*/
#define PRIMECAP 2015177
#define ROOTPCAP 1420
#define NUMPRIMES 150000

#define PARAFLOOR 1000 //g_length required before parallelising the various functions. Prevent overhead

#define UP 1 //Denotes prime, uniform, ascending, 0, and + step seqeuences
#define DOWN -1 //Denotes descending, and -step sequences
#define SEQUP 2
#define SEQDOWN -2
#define UNIFORM 3
#define PRIME 4


 int64_t ** g_vec_properties = NULL;
 static int64_t sieve_size = 0; 
 static int64_t g_seed = 0;

//Compare function for two int64s to sort in ascending order
 static inline int int64Ascend(const void *x, const void *y)
 {
    return (*(int64_t *)x > *(int64_t *)y) - (*(int64_t *)x < *(int64_t *)y);
}

//Compare function for two int64s to sort in descending order
static inline int int64Descend(const void *x, const void *y)
{
    return (*(int64_t *)x < *(int64_t *)y) - (*(int64_t *)x > *(int64_t *)y);
}

void vecInfo(int64_t field, int64_t data)
{
    g_vec_properties[g_nstored][field] = data;
}

////////////////////////////////
///    PARALLEL FUNCTIONS    ///
////////////////////////////////

typedef struct {
    int64_t * vectorA, * vectorB, *vectorC;
    int64_t start, end, result, info1;
    // char padding [8];
} wargs;

static inline void thread_ops (void *(*worker)(void*), wargs wargs[g_nthreads], int num_threads)
{
    pthread_t threads[num_threads];
    for (int64_t i = 0; i < num_threads; i++)
        pthread_create(&threads[i], NULL, worker, &wargs[i]);

    for (int64_t i = 0; i < num_threads; i++)
        pthread_join(threads[i], NULL);
}

//Get count of a function
void * freq_worker (void *arg) 
{
    wargs * counter = (wargs *) arg;
    register int64_t count = 0;
    register int64_t check = counter->info1;

    for (int i = counter->start; i < counter->end; i++)
    {
        if (counter->vectorA[i] == check)
            count++;
    }   
    counter->result = count;
    return NULL;
}

int64_t para_freq (int64_t * vector, int64_t toCheck) {
    wargs freq[g_nthreads];
    int64_t splitter = g_length/g_nthreads;

    for (int i = 0; i < g_nthreads; i++)
    {
        freq[i].vectorA = vector;
        freq[i].info1 = toCheck;
        freq[i].start = i * splitter;
        freq[i].end = (i+1) * splitter;
    }
    freq[g_nthreads-1].end = g_length;

    thread_ops (freq_worker, freq, g_nthreads);

    int frequency = 0;
    for (int i = 0; i < g_nthreads; i++)
    {
        frequency += freq[i].result;
    }
    return frequency;
}

//Get count of a function
void * max_worker (void *arg) 
{
    wargs * maximum = (wargs *) arg;
    register int64_t max = 0;
    for (int64_t i = maximum->start; i < maximum->end; i++)
    {
        if (maximum->vectorA[i] > max)
        {
            max = maximum->vectorA[i];
        }
    }
    maximum->result = max;
    return NULL;
}

int64_t para_max (int64_t * vector) {
    wargs maximum[g_nthreads];
    int64_t splitter = g_length/g_nthreads;

    for (int i = 0; i < g_nthreads; i++)
    {
        maximum[i].vectorA = vector;
        maximum[i].start = i * splitter;
        maximum[i].end = (i+1) * splitter;
    }
    maximum[g_nthreads-1].end = g_length;

    thread_ops(max_worker, maximum, g_nthreads);

    int64_t max = 0;
    for (int i = 0; i < g_nthreads; i++)
    {   
        if (maximum[i].result > max)
            max = maximum[i].result;
    }
    g_existCalcs[arg1].maximum = max;
    return max;
}

typedef struct {
    int64_t * vector;
    int64_t localMin, start, end;
} paraMin;

void * min_worker (void *arg)
{
    wargs * minimum = (wargs *) arg;
    register int64_t min = minimum->vectorA[minimum->start];

    for (int64_t i = minimum->start; i < minimum->end; i++)
    {
        if (minimum->vectorA[i] < min)
        {
            min = minimum->vectorA[i];
        }
    }
    minimum->result = min;
    return NULL;
}

int64_t para_min (int64_t * vector) {
    wargs minimum[g_nthreads];
    int64_t splitter = g_length/g_nthreads;
    for (int i = 0; i < g_nthreads; i++)
    {
        minimum[i].vectorA = vector;
        minimum[i].start = i * splitter;
        minimum[i].end = (i+1) * splitter;
    }
    minimum[g_nthreads-1].end = g_length;

    thread_ops(min_worker, minimum, g_nthreads);

    int64_t min = INT64_MAX;
    for (int i = 0; i < g_nthreads; i++)
    {   
        if (minimum[i].result < min)
            min = minimum[i].result;
    }
    g_existCalcs[arg1].minimum = min;
    return min;
}

static inline void * scalmul_worker (void *arg)
{
    wargs * scalmul = (wargs *) arg;
    int64_t scalar = scalmul->info1, end = scalmul->end;

    for (int64_t i = scalmul->start; i < end/6*6; i+=6)
    {
        scalmul->vectorB[i] = scalmul->vectorA[i] * scalar;
        scalmul->vectorB[i+1] = scalmul->vectorA[i+1] * scalar;
        scalmul->vectorB[i+2] = scalmul->vectorA[i+2] * scalar;
        scalmul->vectorB[i+3] = scalmul->vectorA[i+3] * scalar;
        scalmul->vectorB[i+4] = scalmul->vectorA[i+4] * scalar;
        scalmul->vectorB[i+5] = scalmul->vectorA[i+5] * scalar;
    }

    for (int64_t i = end/6*6; i < end; i++)
        scalmul->vectorB[i] = scalmul->vectorA[i] * scalar;

    return NULL;
}

void para_scalar_mul (int64_t * vector, int64_t * result, int64_t scalar) {
    wargs scalmul[g_nthreads];
    int64_t splitter = g_length/g_nthreads;

    for (int i = 0; i < g_nthreads; i++)
    {
        scalmul[i].vectorA = vector;
        scalmul[i].info1 = scalar;
        scalmul[i].vectorB = result;
        scalmul[i].start = i * splitter;
        scalmul[i].end = (i+1) * splitter;
    }
    scalmul[g_nthreads-1].end = g_length;

    thread_ops(scalmul_worker, scalmul, g_nthreads);
}

static inline void * scaladd_worker (void *arg)
{
    wargs * scaladd = (wargs *) arg;
    int64_t scalar = scaladd->info1, end = scaladd->end;

    for (int64_t i = scaladd->start; i < end/6*6; i+=6)
    {
        scaladd->vectorB[i] = scaladd->vectorA[i] + scalar;
        scaladd->vectorB[i+1] = scaladd->vectorA[i+1] + scalar;
        scaladd->vectorB[i+2] = scaladd->vectorA[i+2] + scalar;
        scaladd->vectorB[i+3] = scaladd->vectorA[i+3] + scalar;
        scaladd->vectorB[i+4] = scaladd->vectorA[i+4] + scalar;
        scaladd->vectorB[i+5] = scaladd->vectorA[i+5] + scalar;
    }

    for (int64_t i = end/6*6; i < end; i++)
        scaladd->vectorB[i] = scaladd->vectorA[i] + scalar;

    return NULL;
}

void para_scalar_add (int64_t * vector, int64_t * result, int64_t scalar) {
    wargs scaladd[g_nthreads];
    int64_t splitter = g_length/g_nthreads;

    for (int i = 0; i < g_nthreads; i++)
    {
        scaladd[i].vectorA = vector;
        scaladd[i].info1 = scalar;
        scaladd[i].vectorB = result;
        scaladd[i].start = i * splitter;
        scaladd[i].end = (i+1) * splitter;
    }
    scaladd[g_nthreads-1].end = g_length;

    thread_ops(scaladd_worker, scaladd, g_nthreads);
}

//Adds two non-pattern vectors
static inline void * vecAdd_worker (void *arg)
{
    wargs * vecAdd = (wargs *) arg;
    int64_t end = vecAdd->end;

    for (int64_t i = vecAdd->start; i < end/6*6; i+=6)
    {
        vecAdd->vectorC[i] = vecAdd->vectorA[i] + vecAdd->vectorB[i];
        vecAdd->vectorC[i+1] = vecAdd->vectorA[i+1] + vecAdd->vectorB[i+1];
        vecAdd->vectorC[i+2] = vecAdd->vectorA[i+2] + vecAdd->vectorB[i+2];
        vecAdd->vectorC[i+3] = vecAdd->vectorA[i+3] + vecAdd->vectorB[i+3];
        vecAdd->vectorC[i+4] = vecAdd->vectorA[i+4] + vecAdd->vectorB[i+4];
        vecAdd->vectorC[i+5] = vecAdd->vectorA[i+5] + vecAdd->vectorB[i+5];
    }

    for (int64_t i = end/6*6; i < end; i++)
        vecAdd->vectorC[i] = vecAdd->vectorA[i] + vecAdd->vectorB[i];
    
    return NULL;
}

//Adds two sequences - O(1) access to the original vectors
static inline void * vecAdd_worker_2SEQ (void *arg)
{
    wargs * vecAdd = (wargs *) arg;
    int64_t start = vecAdd->vectorA[0] + vecAdd->vectorB[1];
    int64_t step = g_vec_properties[arg1][2] + g_vec_properties[arg2][2];
    int64_t end = vecAdd->end;

    for (int64_t i = vecAdd->start; i < end/6*6; i+=6)
    {
        vecAdd->vectorC[i] = start + i * step;
        vecAdd->vectorC[i+1] = start + (i+1) * step;
        vecAdd->vectorC[i+2] = start + (i+2) * step;
        vecAdd->vectorC[i+3] = start + (i+3) * step;
        vecAdd->vectorC[i+4] = start + (i+4) * step;
        vecAdd->vectorC[i+5] = start + (i+5) * step;
    }

    for (int64_t i = end/6*6; i < end; i++)
        vecAdd->vectorC[i] = start + i * step;
    
    return NULL;
}

//Adds a sequence with something else - access to only one vector
static inline void * vecAdd_worker_ASEQ (void *arg)
{
    wargs * vecAdd = (wargs *) arg;
    int64_t start = vecAdd->vectorA[0];
    int64_t step = vecAdd->vectorA[1] - vecAdd->vectorA[0];
    int64_t end = vecAdd->end;

    for (int64_t i = vecAdd->start; i < end/6*6; i+=6)
    {
        vecAdd->vectorC[i] = vecAdd->vectorB[i] + start + i * step;
        vecAdd->vectorC[i+1] = vecAdd->vectorB[i+1] + start + (i+1) * step;
        vecAdd->vectorC[i+2] = vecAdd->vectorB[i+2] + start + (i+2) * step;
        vecAdd->vectorC[i+3] = vecAdd->vectorB[i+3] + start + (i+3) * step;
        vecAdd->vectorC[i+4] = vecAdd->vectorB[i+4] + start + (i+4) * step;
        vecAdd->vectorC[i+5] = vecAdd->vectorB[i+5] + start + (i+5) * step;
    }
    
    for (int64_t i = end/6*6; i < end; i++)
        vecAdd->vectorC[i] = vecAdd->vectorB[i] + start + i * step;
    
    return NULL;
}

//Adds a sequence with something else - access to only one vector
static inline void * vecAdd_worker_BSEQ(void *arg)
{
    wargs * vecAdd = (wargs *) arg;
    int64_t start = vecAdd->vectorB[0];
    int64_t step = vecAdd->vectorB[1] - vecAdd->vectorB[0];
    int64_t end = vecAdd->end;
    
    for (int64_t i = vecAdd->start; i < end/6*6; i+=6)
    {
        vecAdd->vectorC[i] = vecAdd->vectorA[i] + start + i * step;
        vecAdd->vectorC[i+1] = vecAdd->vectorA[i+2] + start + (i+1) * step;
        vecAdd->vectorC[i+2] = vecAdd->vectorA[i+3] + start + (i+2) * step;
        vecAdd->vectorC[i+3] = vecAdd->vectorA[i+4] + start + (i+3) * step;
        vecAdd->vectorC[i+4] = vecAdd->vectorA[i+5] + start + (i+4) * step;
        vecAdd->vectorC[i+5] = vecAdd->vectorA[i+5] + start + (i+5) * step;
    }

    for (int64_t i = end/6*6; i < end; i++)
        vecAdd->vectorC[i] = vecAdd->vectorA[i] + start + i * step;

    return NULL; 
}

void para_vector_add (void *(*worker)(void*), int64_t * vector1, int64_t * vector2, int64_t * result) {
    wargs vecAdd[g_nthreads];
    int64_t splitter = g_length/g_nthreads;
    
    for (int i = 0; i < g_nthreads; i++)
    {
        vecAdd[i].vectorA = vector1;
        vecAdd[i].vectorB = vector2;
        vecAdd[i].vectorC = result;
        vecAdd[i].start = i * splitter;
        vecAdd[i].end = (i+1) * splitter;
    }
    vecAdd[g_nthreads-1].end = g_length;

    thread_ops(worker, vecAdd, g_nthreads);
}

static inline void * vecMul_worker (void *arg)
{
    wargs * vecMul = (wargs *) arg;
    int64_t end = vecMul->end;

    for (int64_t i = vecMul->start; i < end/6*6; i+=6)
    {
        vecMul->vectorC[i] = vecMul->vectorA[i] * vecMul->vectorB[i];
        vecMul->vectorC[i+1] = vecMul->vectorA[i+1] * vecMul->vectorB[i+1];
        vecMul->vectorC[i+2] = vecMul->vectorA[i+2] * vecMul->vectorB[i+2];
        vecMul->vectorC[i+3] = vecMul->vectorA[i+3] * vecMul->vectorB[i+3];
        vecMul->vectorC[i+4] = vecMul->vectorA[i+4] * vecMul->vectorB[i+4];
        vecMul->vectorC[i+5] = vecMul->vectorA[i+5] * vecMul->vectorB[i+5];
    }

    for (int64_t i = end/6*6; i < end; i++)
        vecMul->vectorC[i] = vecMul->vectorA[i] * vecMul->vectorB[i];

    return NULL;
}

void para_vector_mul (int64_t * vector1, int64_t * vector2, int64_t * result) {
    wargs vecMul[g_nthreads];
    int64_t splitter = g_length/g_nthreads;
    for (int i = 0; i < g_nthreads; i++)
    {
        vecMul[i].vectorA = vector1;
        vecMul[i].vectorB = vector2;
        vecMul[i].vectorC = result;
        vecMul[i].start = i * splitter;
        vecMul[i].end = (i+1) * splitter;
    }
    vecMul[g_nthreads-1].end = g_length;

    thread_ops(vecMul_worker, vecMul, g_nthreads);
}

static inline void * sums_worker (void *arg)
{
    wargs * sums = (wargs *) arg;
    int64_t sum = 0;
    for (int64_t i = sums->start; i < sums->end; i++)
    {
        sum += sums->vectorA[i];
    }

    sums->result = sum;
    return NULL;
}

int64_t para_sum (int64_t * vector, int n_split)
{
    wargs sums[n_split];
    int64_t splitter = g_length/n_split;
    for (int i = 0; i < n_split; i++)
    {
        sums[i].vectorA = vector;
        sums[i].start = i * splitter;
        sums[i].end = (i+1) * splitter;
    }
    sums[n_split-1].end = g_length;

    thread_ops(sums_worker, sums, n_split);

    int64_t sum = 0;
    for (int i = 0; i < n_split; i++)
        sum += sums[i].result;

    g_existCalcs[arg1].sum = sum;
    return sum;
}

static inline void * clone_worker (void *arg)
{
    wargs * clone = (wargs *) arg;
    for (int64_t i = clone->start; i < clone->end; i++)
    {
        clone->vectorB[i] = clone->vectorA[i];
    }
    return NULL;
}

void para_clone (int64_t * clone, int64_t * original)
{
    wargs cloned[g_nthreads];

    int64_t splitter = g_length/g_nthreads;
    for(int i = 0; i < g_nthreads; i++)
    {
        cloned[i].vectorA = original;
        cloned[i].vectorB = clone;
        cloned[i].start = i * splitter;
        cloned[i].end = (i+1) * splitter;
    }
    cloned[g_nthreads-1].end = g_length;

    thread_ops(clone_worker, cloned, g_nthreads);
}

static inline void * reverse_worker (void *arg)
{
    wargs * reverse = (wargs *) arg;
    for (int64_t i = reverse->start; i < reverse->end; i++)
    {
        reverse->vectorB[i] = reverse->vectorA[g_length-1-i];
    }
    return NULL;
}

void para_reverse (int64_t * reversed, int64_t * original)
{
    wargs rev[g_nthreads];

    int64_t splitter = g_length/g_nthreads;
    for (int i = 0; i < g_nthreads; i++)
    {
        rev[i].vectorA = original;
        rev[i].vectorB = reversed;
        rev[i].start = i * splitter;
        rev[i].end = (i+1) * splitter;
    }

    thread_ops(reverse_worker, rev, g_nthreads);
}

int64_t * clonePart (int64_t * vector, int64_t start, int64_t end)
{
    int64_t* clone = (int64_t *)calloc(end - start, sizeof(int64_t));
    int64_t index = 0;
    for (int64_t i = start; i < end; i++) {
        clone[index++] = vector[i];
    }
    return clone;
}

bool powOf2(int x) 
{
    return x && !(x & (x - 1));
}

static inline void swap(int64_t * vector, int64_t a, int64_t b)
{
    int64_t tmp = vector[a];
    vector[a] = vector[b];
    vector[b] = tmp;
}

int64_t quickSelect (int64_t * vector, int64_t start, int64_t len, int goal)
{
    int64_t begin, end, pivL, pivR, mid, midChk;
    begin = start;
    end = len - 1;
    while (1) 
    {
        if (end <= begin + 1) 
        { 
            if (end == begin + 1 && vector[end] < vector[begin]) 
                swap(vector, begin, end);
            return vector[goal];
        }

        else 
        {
            if ((begin+end) % 2 == 0)
                mid=((begin + end) - 1) / 2;

            else
                mid=(begin + end) / 2; 

            swap(vector, mid, begin+1);

            if (vector[begin] > vector[end])
                swap(vector, begin, end);

            if (vector[begin + 1] > vector[end])
                swap(vector, begin + 1, end);

            if (vector[begin] > vector[begin + 1])
                swap(vector, begin, begin + 1);

            pivL = begin + 1; 
            pivR = end;
            midChk = vector[begin + 1]; 

            while (1)
            { 
                do pivL++ ;
                while (vector[pivL] < midChk);
                do pivR--;
                while (vector[pivR] > midChk);
                
                if (pivR < pivL)
                    break; 
                
                swap(vector, pivL, pivR);
            } 

            vector[begin + 1] = vector[pivR]; 
            vector[pivR] = midChk;

            if (pivR >= goal)
                end = pivR - 1; 

            if (pivR <= goal) 
                begin = pivL;
        }
    }
}


static inline int getdigit (int64_t num, int n)
{
    static int64_t tenpows[] = {1, 10, 100, 1000, 10000, 100000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000};
    return ((num / tenpows[n]) % 10);
}

void radix_sort (int64_t * vec, int64_t len)
{
    //////////////////////////////////////////////////
    //////////////IMPORANT: TO DO/////////////////////
    //////////////////////////////////////////////////
    ////REMOVE ASSUMPTION THAT VEC LEN = 5////////////
    //////////////////////////////////////////////////

    //Bucket of values based on digit check
    int64_t ** bucket = calloc(10, sizeof(int64_t *));
    //Count of values in each bucket
    int64_t * buckCnt = calloc(10, sizeof(int64_t *));

    int64_t max = get_maximum(vec);
   int64_t digits = 1;
   while (max)
   {
        max /= 10;
        digits++;
   }

    for (int i = 0; i < digits; i++)
    {

        buckCnt = calloc(10, sizeof(int64_t *));
        for (int m = 0; m <10; m++)
            bucket[m] = calloc(g_length, sizeof(int64_t));

        for (int64_t j = 0; j < len; j++)
        {
            int curDigit = getdigit(vec[j], i);
            bucket[curDigit][buckCnt[curDigit]] = vec[j];
            buckCnt[curDigit]++;
        }

        int64_t put = 0;
        for (int i = 0; i < 10; i++)
        {
            for (int64_t j = 0; j < buckCnt[i]; j++)
            {
                vec[put++] = bucket[i][j];
            }
        }
    }

    for (int i = 0; i < digits; i++)
        free(bucket[i]);
    free(bucket);
    free(buckCnt);
}

////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////

//Sets the number of elements that each vector will contain
void set_length(int64_t length) 
{
    g_length = length;
}

//Sets the seed used when generating pseudorandom numbers
void set_seed(int64_t seed) 
{
    g_seed = seed;
}

//Returns pseudorandom number determined by the seed
int64_t fast_rand(void) 
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

////////////////////////////////
///   VECTOR INITALISATIONS  ///
////////////////////////////////

//Returns new vector, with all elements set to zero
int64_t* new_vector(void) {
    return (int64_t *) calloc(g_length, sizeof(int64_t));
}

//Returns new vector, with all elements set to given value
int64_t* uniform_vector(int64_t value)
{
    int64_t* vector = new_vector();

    for (int64_t i = 0; i < g_length; i++) {
        vector[i] = value;
    }

    return vector;
}

//Returns new vector, with elements generated at random using given seed
int64_t* random_vector(int64_t seed) {
    int64_t* vector = new_vector();
    set_seed(seed);

    for (int64_t i = 0; i < g_length; i++) {
        vector[i] = fast_rand();
    }

    return vector;
}

bool is_prime_serial(int64_t number) {
    int64_t limit = sqrt(number) + 1;
    int64_t i = 0;
    while (primeList[i] < limit)
    {
        if (number % primeList[i++] == 0)
            return false;
    }

    if (PRIMECAP >= limit)
        return true;

    for (i = PRIMECAP + 2; i < limit; i+=2)
    {
        if(number % i == 0)
            return false;
    }

    return true;
}

//Returns new vector for bit primes, all elements set to zero
int * new_sieve(int64_t size)
{
    return (int *) calloc(size, sizeof(int));
}

//Returns new array of booleans
bool * new_bool_sieve(int64_t size)
{
    return (bool *) calloc(size, sizeof(bool));
}

//Performs primitive sieve of Eratosthenes on a given boolean array
void generateSieve (bool * sieve, int64_t sievesize, int64_t root)
{
    int64_t i = 0;
    for (i = 2; i < sievesize; i++)
        sieve[i] = true;

    int64_t j = 0;
    for (i = 2; i < root; i++)
    {
        if (sieve[i])
        {    
            for (j = i; i * j < sievesize; j++)
                sieve [i * j] = false;
        }
    }
}

//Primitive sieve of Eratosthenes
void sieveofEratosthenes (int64_t * primes, int64_t number, int64_t root)
{
    bool * sieve = new_bool_sieve(sieve_size);
    int64_t i = 0;
    generateSieve(sieve, sieve_size, root);

    int64_t j = 0;

    if (number < 0) number = 2;
    if (number % 2 == 0 && number != 2) 
        number++;
    else if (number == 2)
    { 
        primes[j] = 2; 
        number++; 
        j++;
    }
    for (i = number; i < sieve_size; i+=2)
    {
        if (sieve[i])
            primes[j++] = i;
        if (j == g_length)
            break;
    }
    free(sieve);    
}

    //Segmented sieve of Eratosthenes
void segSieveOfEratosthenes (int64_t * primes, int64_t root, int64_t begin, int64_t factor)
{
    //Generate all primes up to the root of the max prime needed (taking into account starting number)
    bool * initialSieve = (bool *) calloc(root, sizeof(bool));              //Bool sieve, all set to false
    int64_t iSsize = sqrt(root) + 1;                                        //Root (of root) for primitive SOE generation
    generateSieve (initialSieve, root, iSsize);                             //Complete the bool sieve and set true/false for primes
    int64_t * primesToSeg = (int64_t *) calloc(root, sizeof(int64_t));      
    int64_t segPrimesCount = 0;
    int64_t length_count = 0;

    //Fill primes list with all prime/true values from boolean sieve
    for (int64_t i = 3; i < root; i += 2)
    {
        if(initialSieve[i])
        {
            primesToSeg[segPrimesCount++] = i;
        }
    }
    primesToSeg = (int64_t *)realloc(primesToSeg, segPrimesCount * sizeof(int64_t));    //Reallocte prime sieve based on exact size

    /*Find a good block size to sieve in
    Even numbers for block size and start number because
    -0.5 * (start + 1 + prime) (see below for what this is)
    nees to always be an even number for this seg sieve to work*/
    int64_t blockSize = segPrimesCount * factor;
    
    while (!(blockSize % 2 == 0))   
        blockSize++;

    int64_t start = begin;

    if (start <= 2)
    {
        primes[length_count++] = 2;
        start = 3;                      //Account for start = 0 and 1
    }

    if (start % 2 == 0)
        start++;                        //Allow odd skips in next step

    //Check if starting number is within primes list!
    if (start <= primesToSeg[segPrimesCount-1])
    {
        for (int64_t i = start; i < root; i += 2)
        {
            if (initialSieve[i])
                primes[length_count++] = i;
            if (length_count == g_length)
            {
                free(initialSieve);
                free(primesToSeg);
                return; 
            }
        }
    }

    //Initial sieve is no longer needed
    free(initialSieve);

    if (start % 2 != 0)
        start--;

    //Begin sieving by blocks. THIS IS THE MAIN PART OF THE SEGMENTED SIEVE
    for (; start < sieve_size; start += blockSize)
    {
    //Create boolean array of size blocksize, initialise all to true
        bool * sieveBlock = (bool *)calloc(blockSize/2, sizeof(bool));
        for (int64_t j = 0; j < blockSize/2; j++)
            sieveBlock[j] = true;


        /*Create offset values. The sieve boolean array is half the block
        size instead of the full size because rather than using the generic
        -Lmodp offset which requires more memory space, the alternate formula
        -1/2(L+1+p)modp, where L is the start value of the block, and p is the
        prime number, should mean it can operate fractionally faster*/
        int64_t * primeOffsets = (int64_t *)calloc(segPrimesCount, sizeof(int64_t));
        for (int64_t j = 0; j < segPrimesCount; j++)
        {
            int64_t qval = (start+1+primesToSeg[j])/(-2) % primesToSeg[j];
            if (qval < 0)
                qval += primesToSeg[j];     //Q = -0.5 * (start + 1 + prime) % prime
            primeOffsets[j] = qval;
        }

        /*Go through the current block and set all the appropriate elements of the
        bool sieve to false according to the prime number and the q-value linked to it*/
        for (int64_t j = 0; j < segPrimesCount; j++)
        {
            for (int64_t k = primeOffsets[j]; k < blockSize/2; k += primesToSeg[j])
            {
                sieveBlock[k] = false;
            }
        }

        /*Go through the sieved block and take out all the primes until length is reached.*/
        for (int64_t i = 0; i < blockSize/2; i++)
        {
            if (length_count >= g_length)
            {
                free(sieveBlock); 
                free(primeOffsets);
                free(primesToSeg);
                return;
            }
            if (sieveBlock[i])
                primes[length_count++] = start + i * 2 + 1;  //The true values' indexes are extrated to give the true prime number in this way
            
        }
        free(sieveBlock); 
        free(primeOffsets);            
    }
}

// Returns new vector, containing primes numbers in sequence from given start
int64_t * prime_vector(int64_t start)
{
    int64_t * primes = new_vector();

    int64_t number = start;

    //Modified nth prime estimation - always overestimate
    int lenLog = 0;
    if (g_length > 1)
        lenLog = log(g_length) + log(log(g_length)) + 1;
    else lenLog = 1;

    //Modified X primes under n multiplier - always overestimate
    int64_t numPrimes = 4;
    if (number >= 10)
        numPrimes = number/(log(number) - 1);
    int numLog = log(numPrimes) + log(log(numPrimes)) + 1;

    //Determine an overestimate of a number higher than the higher prime needed
    sieve_size = (g_length * lenLog) + (numPrimes * numLog);
    if (sieve_size < 6) sieve_size = 6;

    int64_t root = ceil(sqrt(sieve_size)); //Pass root value into functions. Ceil + Sqrt operations cost


    if (g_length > 2000)
           segSieveOfEratosthenes(primes, root, start, 200);
    else
       segSieveOfEratosthenes(primes, root, start, 2);
   
    return primes;
}

//Returns new vector, with elements in sequence from given start and step
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

//Returns new vector, cloning elements from given vector
int64_t* cloned(int64_t* vector) 
{
    int64_t* clone = new_vector();
    if (g_length > PARAFLOOR)
    {
        para_clone(clone, vector);
        return clone;
    }

    for (int64_t i = 0; i < g_length; i++) {
        clone[i] = vector[i];
    }

    return clone;
}

//Returns new vector, with elements ordered in reverse
int64_t* reversed(int64_t* vector) 
{

    int64_t* result = new_vector();
    if (g_vec_properties[arg1][0] == UNIFORM)
        return cloned(vector);
    else if (g_length > PARAFLOOR)
    {
        para_reverse(result, vector);
        return result;
    }

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[g_length - 1 - i];
    }

    return result;
}

//Returns new vector, with elements ordered from smallest to largest
int64_t* ascending(int64_t* vector) 
{   
    int64_t* result = cloned(vector);

    if (g_vec_properties[arg1][0] == UP || g_vec_properties[arg1][0] == UNIFORM || 
        g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == PRIME)
    {
        return result;
    }
    else if (g_vec_properties[arg1][0] == DOWN || g_vec_properties[arg1][0] == SEQDOWN)
    {
        free(result);
        return reversed(vector);
    }
    
    
    qsort(result, g_length, sizeof(int64_t), int64Ascend);
    return result;
}

//Returns new vector, with elements ordered from largest to smallest
int64_t* descending(int64_t* vector) 
{
    int64_t* result = cloned(vector);
    
    if (g_vec_properties[arg1][0] == DOWN || g_vec_properties[arg1][0] == UNIFORM || g_vec_properties[arg1][0] == SEQDOWN)
    {
        return result;
    }
    else if (g_vec_properties[arg1][0] == UP || g_vec_properties[arg1][0] == SEQUP ||
        g_vec_properties[arg1][0] == PRIME)
    {
        free(result);
        return reversed(vector);
    }

    qsort(result, g_length, sizeof(int64_t), int64Descend);
    return result;
}

// Returns new vector, adding scalar to each element
int64_t* scalar_add(int64_t* vector, int64_t scalar) {
    int64_t* result = new_vector();
    if (g_length > PARAFLOOR)
    {
        para_scalar_add (vector, result, scalar);
        return result;
    }

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[i] + scalar;
    }
    
    return result;
}

//Returns new vector, multiplying scalar to each element
int64_t* scalar_mul(int64_t* vector, int64_t scalar)
{
    if (arg2 == 0)
        return zerovector;
    if (arg2 == 1)
        return vector;

    int64_t* result = new_vector();
    if (g_length > PARAFLOOR)
    {
        para_scalar_mul (vector, result, scalar);
        return result;
    }

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector[i] * scalar;
    }

    return result;
}

//Returns new vector, adding elements with the same index
int64_t* vector_add(int64_t* vector_a, int64_t* vector_b)
{
    //Returns vector_b if a is all 0s
    if (g_vec_properties[arg1][0] == UNIFORM && vector_a[0] == 0)
        return cloned(vector_b);
    
    //Returns vector_a if b is all 0s
    else if (g_vec_properties[arg2][0] == UNIFORM && vector_b[0] == 0)
        return cloned(vector_a);
    
    //Returns an array of 0s if the two are uniform inverses
    else if (g_vec_properties[arg1][0] == UNIFORM && g_vec_properties[arg2][0] == UNIFORM &&
        vector_a[0] == -vector_b[0])
        return cloned(zerovector);

    int64_t* result = new_vector();
    
    //Parallelise if above a certain threshold for g_length
    if (g_length > PARAFLOOR)
    {
         //If both vectors are sequences - no memory reading of either is required for the result
        if ((g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == SEQDOWN) &&
            (g_vec_properties[arg2][0] == SEQUP || g_vec_properties[arg2][0] == SEQDOWN))
        {
            para_vector_add (vecAdd_worker_2SEQ, vector_a, vector_b, result);
            return result;
        }

         //If only one vector is a sequence - only the non-seq needs to be read
        if (g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == SEQDOWN)
        {
            para_vector_add(vecAdd_worker_ASEQ, vector_a, vector_b, result);
            return result;
        }

        if (g_vec_properties[arg2][0] == SEQDOWN || g_vec_properties[arg1][0] == SEQUP)
        {
            para_vector_add(vecAdd_worker_BSEQ, vector_a, vector_b, result);
            return result;
        }

        para_vector_add (vecAdd_worker, vector_a, vector_b, result);
        return result;
    }

    //Simply sequential
    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector_a[i] + vector_b[i];
    }

    return result;
}

//Returns new vector, multiplying elements with the same index
int64_t* vector_mul(int64_t* vector_a, int64_t* vector_b)
{
    if (g_vec_properties[arg1][0] == UNIFORM)
    {
        if (vector_a[0] == 0)
            return cloned(zerovector);
        if (vector_a[0] == 1)
            return cloned(vector_b);
    }
    if (g_vec_properties[arg2][0] == UNIFORM)
    {
        if (vector_b[0] == 0)
            return cloned(zerovector);
        if (vector_a[0] == 1)
            return cloned(vector_a);
    }

    int64_t* result = new_vector();
    if (g_length > PARAFLOOR)
    {
        para_vector_mul (vector_a, vector_b, result);
        return result;
    }

    for (int64_t i = 0; i < g_length; i++) {
        result[i] = vector_a[i] * vector_b[i];
    }

    return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

int64_t get_sum_vectorised (int64_t * vector)
{

    __m128i sum = _mm_setzero_si128();
    int64_t actualSum = 0;

    
    for (int64_t i = 0; i < g_length/4*4; i += 4)
    {
        __m128i temp = _mm_loadu_si128((__m128i *)(vector + i));
        sum = _mm_add_epi64(sum, temp);
    }

    int64_t A[4] = {0,0,0,0};
    _mm_storeu_si128((__m128i *)A, sum);
    actualSum += A[0] + A[1] + A[2] + A[3];
    
    for (int64_t i = g_length/4*4; i < g_length; i++)
    {
        actualSum += vector[i];
    }

    return actualSum;   
}

//Returns the sum of all elements
int64_t get_sum(int64_t* vector) 
{

    // return get_sum_vectorised(vector);
    if (g_existCalcs[arg1].sumFlag == 1)
    {
        return g_existCalcs[arg1].sum;
    }
    g_existCalcs[arg1].sumFlag = 1;

    if (g_vec_properties[arg1][0] == UNIFORM)
    {
        int64_t sum = vector[0] * g_length;
        g_existCalcs[arg1].sum = sum;
        return sum;
    }

    if (g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == SEQDOWN)
    {
        int64_t sum = (g_length * (vector[0] + vector[g_length-1]))/2;
        g_existCalcs[arg1].sum = sum;
        return sum;
    }

    if (g_length > PARAFLOOR)
    {
        return para_sum (vector, g_nthreads);
    }

    int64_t sum = 0;

    for (int64_t i = 0; i < g_length; i++) {
        sum += vector[i];
    }
    
    g_existCalcs[arg1].sum = sum;
    return sum;
}

//Creates what is essentially a bucket holding counts of all 
//values in the array. Requires O(k) space and runs in O(n)
//time, where k = the maximum - minimum (+1) values in the array 
int64_t fast_mode(int64_t* vector, int64_t max, int64_t min)
{
    int64_t * bucket = calloc(max - min + 1, sizeof(int64_t));

    int64_t currMode = vector[0], count = 0;
    bool repeatMode = true;

    for (int64_t i = 0; i < g_length; i++)
    {
        int64_t bucketIndex = vector[i] - min;  //Offset by the minimum to get index of element in bucket
        bucket[bucketIndex]++;

        if (bucket[bucketIndex] > count)        //New mode. Resets mode contention condition
        {
            count = bucket[bucketIndex];
            currMode = vector[i];
            repeatMode = false;
        }
        else if (bucket[bucketIndex] == count)  //If there is a mode contention - temporary mode = -1
            repeatMode = true;
    }
    
    free(bucket);
    if (repeatMode)
    {
        g_existCalcs[arg1].mode = -1;
        return -1;
    }

    else
    {
        g_existCalcs[arg1].mode = currMode;
        return currMode;
    }
}

//Returns the most frequently occuring element
//or -1 if there is no such unique element
int64_t get_mode(int64_t* vector) 
{

    if (g_existCalcs[arg1].modeFlag == 1)
        return g_existCalcs[arg1].mode;
    g_existCalcs[arg1].modeFlag = 1;
    
    //Mode of a uniform vector is any element within it
    if (g_length == 1 || g_vec_properties[arg1][0] == UNIFORM)
    {
        g_existCalcs[arg1].mode = vector[0];
        return vector[0];
    }

    //All vectors with unique values have no mode
    else if (g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == SEQDOWN ||
        g_vec_properties[arg1][0] == PRIME)
    {
        g_existCalcs[arg1].mode = -1;
        return -1;
    }

    int64_t max = get_maximum(vector);
    int64_t min = get_minimum(vector);
    if ((max - min) <   2147483648) //16gb can hold 2147483648 int64s
                                    //12gb can hold 1610612736 int64s
                                    //8gb  can hold 1073741824 int64s
                                    //4gb  can hold 536870912  int64s
        return fast_mode(vector, max, min);

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
    {
        g_existCalcs[arg1].mode = -1;
        return -1;
    }

    g_existCalcs[arg1].mode = mode;
    return mode;
}

//Returns the lower median
int64_t get_median(int64_t* vector)
{

    if (g_existCalcs[arg1].medianFlag == 1)
        return g_existCalcs[arg1].median;
    g_existCalcs[arg1].medianFlag = 1;

    //Median of a uniform vector is any value in it - as it equals the lower median
    if (g_length == 1 || g_vec_properties[arg1][0] == UNIFORM)
    {
        g_existCalcs[arg1].median = vector[0];
        return vector[0];
    }

    //Calculate middle index
    int64_t n = g_length/2;
    if (g_length % 2 == 0)
        n = g_length/2 - 1;

    //Median of sorted vectors is simply at n
    if (g_vec_properties[arg1][0] == UP || g_vec_properties[arg1][0] == SEQUP || 
        g_vec_properties[arg1][0] == PRIME)
    {
        g_existCalcs[arg1].median = vector[n];
        return vector[n];
    }
    else if (g_vec_properties[arg1][0] == DOWN || g_vec_properties[arg1][0] == SEQDOWN)
    {
        g_existCalcs[arg1].median = vector[n-1];
        return vector[n-1];
    }

    int64_t * result = vector;

    int64_t median = quickSelect(result, 0, g_length, n);
    g_existCalcs[arg1].median = median;
    return median;
}

//Returns the smallest value in the vector
int64_t get_minimum(int64_t* vector) 
{

    if (g_existCalcs[arg1].minFlag == 1)
        return g_existCalcs[arg1].minimum;
    g_existCalcs[arg1].minFlag = 1;

    //Minimum value of all ascending vectors is the first element
    if (g_length == 1 || g_vec_properties[arg1][0] == UP || 
        g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == UNIFORM ||
        g_vec_properties[arg1][0] == PRIME)
    {
        g_existCalcs[arg1].minimum = vector[0];
        return vector[0];
    }

    //Mminimum value of all descnding vectors is the last element
    if (g_vec_properties[arg1][0] == DOWN || g_vec_properties[arg1][0] == SEQDOWN)
    {
        g_existCalcs[arg1].minimum = vector[g_length-1];
        return vector[g_length-1];
    }

    if (g_length > PARAFLOOR)
        return para_min (vector);

    int64_t min = INT64_MAX;
    for (int i = 0; i < g_length; i++) {
        if (vector[i] < min)
            min = vector[i];
    }
    g_existCalcs[arg1].minimum = min;
    return min;
}

//Returns the largest value in the vector
int64_t get_maximum(int64_t* vector)
{

    if (g_existCalcs[arg1].maxFlag == 1)
        return g_existCalcs[arg1].maximum;
    g_existCalcs[arg1].maxFlag = 1;

    //Max value for all ascending vectors is the last value
    if (g_length == 1 || g_vec_properties[arg1][0] == UP || 
        g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == UNIFORM ||
        g_vec_properties[arg1][0] == PRIME)
    {
        g_existCalcs[arg1].maximum = vector[g_length-1];
        return vector[g_length-1];
    }

    //Max value for all descending vectors is the first value
    if (g_vec_properties[arg1][0] == DOWN || g_vec_properties[arg1][0] == SEQDOWN)
    {
        g_existCalcs[arg1].maximum = vector[0];
        return vector[0];
    }

    if (g_length > PARAFLOOR)
        return para_max(vector);

    int64_t max = 0;
    for (int i = 0; i < g_length; i++) {
        if (max < vector[i])
            max = vector[i];
    }
    g_existCalcs[arg1].maximum = max;
    return max;
}

//Returns the frequency of the value in the vector
int64_t get_frequency(int64_t* vector, int64_t value) 
{
    //Checks uniform vectors. Frequency is either length or 0.
    if (g_vec_properties[arg1][0] == UNIFORM)
    {
        if (value == vector[0])
            return g_length;
        return 0;
    }

    //Checks sequences. Frequency is either 1 or 0. Performs divide operation and
    //checks for a whole number to determine if the number exists.
    if (g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == SEQDOWN)
    {
        long double wholeNumCheck = ((long double)value - (long double)g_vec_properties[arg1][1])/(long double)g_vec_properties[arg1][2];
        if (ceilf(wholeNumCheck) == wholeNumCheck && wholeNumCheck >= 0)
            return 1;
        else
            return 0;
    }

    if (g_length > PARAFLOOR)
        return para_freq(vector, value);

    int64_t count = 0;
    for (int64_t i = 0; i < g_length; i++) {
        if (vector[i] == value) {
            count ++;
        }
    }
    return count;
}

//Returns the value stored at the given element index
int64_t get_element(int64_t* vector, int64_t index) 
{
    return vector[index];
}

//Output given vector to standard output
void display(int64_t* vector, int64_t label) {

    printf("vector~%" PRId64 " ::", label);

    for (int64_t i = 0; i < g_length; i++) {
        printf(" %" PRId64, vector[i]);
    }

    puts("");
}