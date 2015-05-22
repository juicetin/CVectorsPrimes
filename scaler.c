/*Sublime text is recommended to read this file. The task was to implement a number of mathematical operations in C
to run as quickly as possible using algorithms and optimisation. My ranking was 8th in a cohort of 300 for this 
competitive assignment.*/


#include "scaler.h"
#include "vector.h"

#define UNSET -1
#define MAX_BUFFER 64

#define OUTPUT_INT32_T "%" PRId32 "\n"
#define OUTPUT_INT64_T "%" PRId64 "\n"

#define UP 1 //Ascending
#define DOWN -1 //Denotes descending 
#define SEQUP 2 //Positive step sequences - unique values
#define SEQDOWN -2 //Negative step sequences - unique values
#define UNIFORM 3 //Uniform vectors (including when multiplied by zero, etc.)
#define PRIME 4 //Non-sequential step unique values
#define RAND 5 //Random ONLY vectors. Specifically for radix sorting

int64_t g_nthreads               = 0; /* 1 <= n <= 256 */
int64_t g_length                 = 0; /* 1 <= n <= 1,000,000,000 */
int64_t arg1 			     = UNSET;
int64_t arg2                 = UNSET;
int64_t g_nstored                = 0; /* 1 <= n <= 1,000,000 */
int64_t * zerovector          = NULL;
int64_t g_ncores                 = 0; /* 1 <= n <= 128 */

int64_t * primeList           = NULL;

static int64_t** g_vectors     = NULL;
static int64_t g_nvectors               = 0; /* 1 <= n <= 1,000,000 */
static int64_t g_ncomputations          = 0; /* 1 <= n <= 1,000,000 */

props * g_existCalcs;

extern int64_t ** g_vec_properties;
// extern existCalcs * storeCalcs;

//Release all dynamically allocated memory
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

	for (int64_t i = 0; i < g_nvectors; i++)
		free(g_vec_properties[i]);
	free(g_vec_properties);

	free(g_vectors);
	free(zerovector);
	free(g_existCalcs); 
}

//Outputs error and terminates program
void error(void) {

	puts("error");
	release();
	exit(0);
}

//Defines setting based on given input
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

	g_existCalcs = calloc(g_nvectors, sizeof(props));

	set_length(g_length);

	g_vec_properties = (int64_t **) calloc (g_nvectors, sizeof(int64_t *));
	for (int64_t i = 0; i <g_nvectors; i++)
		g_vec_properties[i] = (int64_t *) calloc(3, sizeof(int64_t));

	zerovector = (int64_t *) calloc (g_length, sizeof(int64_t));
}

//Defines vectors based on given input
void define_vectors(void) 
{
	g_vectors = (int64_t **) calloc(g_nvectors, sizeof(int64_t *));
	arg1 = UNSET;
	arg2 = UNSET;

	for (int64_t i = 0; i < g_nvectors; i++) {

		char input[MAX_BUFFER];
		char function[MAX_BUFFER];

        /* vector= <function> [args] */

		fgets(input, MAX_BUFFER, stdin);
		sscanf(input, "vector= %s %" PRId64 " %" PRId64, function, &arg1, &arg2);

		if (strcmp(function, "primes") == 0)
		{
			vecInfo(0, PRIME);
			g_vectors[g_nstored++] = prime_vector(arg1);
		} 
		else if (strcmp(function, "random") == 0) 
		{
			vecInfo(0, RAND);
			g_vectors[g_nstored++] = random_vector(arg1);
		} 
		else if (strcmp(function, "uniform") == 0)
		{
			vecInfo(0, UNIFORM);
			g_vectors[g_nstored++] = uniform_vector(arg1);
		} 
		else if (strcmp(function, "sequence") == 0)
		{
			if (arg2 >= 0)
				vecInfo(0, SEQUP);
			else
				vecInfo(0, SEQDOWN);
			vecInfo(1, arg1);
			vecInfo(2, arg2);

			g_vectors[g_nstored++] = sequence_vector(arg1, arg2);
		} 
		else if (strcmp(function, "cloned") == 0) {
			vecInfo(0, g_vec_properties[arg1][0]);
			g_vectors[g_nstored++] = cloned(g_vectors[arg1]);
		} 
		else if (strcmp(function, "reversed") == 0) {
			vecInfo(0, g_vec_properties[arg1][0]);
			g_vectors[g_nstored++] = reversed(g_vectors[arg1]);
		} 
		else if (strcmp(function, "ascending") == 0) {
			vecInfo(0, UP);
			g_vectors[g_nstored++] = ascending(g_vectors[arg1]);
		} 
		else if (strcmp(function, "descending") == 0) {
			vecInfo(0, DOWN);
			g_vectors[g_nstored++] = descending(g_vectors[arg1]);
		} 
		else if (strcmp(function, "scalar#add") == 0) {
			//Retains the same value 
			vecInfo(0, g_vec_properties[arg1][0]);
			//TODO: Change start/step values
			g_vectors[g_nstored++] = scalar_add(g_vectors[arg1], arg2);
		} 
		else if (strcmp(function, "scalar#mul") == 0) {
			if (arg2 > 0)
				vecInfo(0, g_vec_properties[arg1][0]);
			else if (arg2 == 0)
				vecInfo(0, UNIFORM);
			else
				vecInfo(0, -g_vec_properties[arg1][0]);
			g_vectors[g_nstored++] = scalar_mul(g_vectors[arg1], arg2);
		} 
		else if (strcmp(function, "vector#add") == 0) {
			if (g_vec_properties[arg1][0] == SEQUP && g_vec_properties[arg2][0] == SEQUP)
			{
				vecInfo(0, SEQUP);
				vecInfo(1, g_vec_properties[arg1][1] + g_vec_properties[arg2][1]);
				vecInfo(2, g_vec_properties[arg1][2] + g_vec_properties[arg2][2]);
			}
			else if (g_vec_properties[arg1][0] == SEQDOWN && g_vec_properties[arg2][0] == SEQDOWN)
			{
				vecInfo(0, SEQDOWN);
				vecInfo(1, g_vec_properties[arg1][1] + g_vec_properties[arg2][1]);
				vecInfo(2, g_vec_properties[arg1][2] + g_vec_properties[arg2][2]);
			}
			else if ((g_vec_properties[arg1][0] == SEQUP || g_vec_properties[arg1][0] == UP) &&
					(g_vec_properties[arg2][0] == SEQUP || g_vec_properties[arg2][0] == UP))
			{
				vecInfo(0, UP);
			}
			else if ((g_vec_properties[arg1][0] == SEQDOWN || g_vec_properties[arg1][0] == DOWN) &&
					(g_vec_properties[arg2][0] == SEQDOWN || g_vec_properties[arg2][0] == DOWN))
			{
				vecInfo(0, DOWN);
			} 

			g_vectors[g_nstored++] = vector_add(g_vectors[arg1], g_vectors[arg2]);
		} 
		else if (strcmp(function, "vector#mul") == 0) {
			if ((g_vec_properties[arg1][0] == UNIFORM || g_vec_properties[arg1][0] == UNIFORM) && 
				(g_vectors[arg1][0] == 0 || g_vectors[arg1][0] == 0))
				vecInfo(0, UNIFORM);
			g_vectors[g_nstored++] = vector_mul(g_vectors[arg1], g_vectors[arg2]);
		} 
		else {
			error();
		}
	}
}

//Runs computations based on given input
void compute_engine(void) {

	for (int64_t i = 0; i < g_ncomputations; i++) {

		char input[MAX_BUFFER];
		char function[MAX_BUFFER];

		arg1 = UNSET;
		arg2 = UNSET;

        /* compute <function> [args] */

		fgets(input, MAX_BUFFER, stdin);
		sscanf(input, "compute %s %" PRId64 " %" PRId64, function, &arg1, &arg2);

		if (strcmp(function, "display") == 0) {
			display(g_vectors[arg1], arg1);
		} else {

			if (arg2 == UNSET) {
				printf("%s of vector~%" PRId64 " = ", function, arg1);
			}

			if (strcmp(function, "sum") == 0)
			{
				printf(OUTPUT_INT64_T, get_sum(g_vectors[arg1]));
			} 
			else if (strcmp(function, "mode") == 0) {
 				// printf("arg from scaler is %ld\n", arg1);
				printf(OUTPUT_INT64_T, get_mode(g_vectors[arg1]));
			} 
			else if (strcmp(function, "median") == 0) {
				printf(OUTPUT_INT64_T, get_median(g_vectors[arg1]));
			} 
			else if (strcmp(function, "minimum") == 0) {
				printf(OUTPUT_INT64_T, get_minimum(g_vectors[arg1]));
			} 
			else if (strcmp(function, "maximum") == 0) {
				printf(OUTPUT_INT64_T, get_maximum(g_vectors[arg1]));
			} 
			else if (strcmp(function, "element") == 0) {
				printf("%s %" PRId64 " in vector~%" PRId64 " = ", function, arg2, arg1);
				printf(OUTPUT_INT64_T, get_element(g_vectors[arg1], arg2));
			} 
			else if (strcmp(function, "frequency") == 0) {
				printf("%s of %" PRId64 " in vector~%" PRId64 " = ", function, arg2, arg1);
				printf(OUTPUT_INT64_T, get_frequency(g_vectors[arg1], arg2));
			} 
			else {
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

 	// typedef float floatSSE[4];
 	// floatSSE a = {0, 2, 4, 6};
 	// floatSSE b = {0, 2, 4, 6};
 	// __m128 ma = _mm_loadu_ps(&a[0]);
 	// __m128 mb = _mm_loadu_ps(&b[0]);

 // 	__m128 mprod = _mm_mul_ps(ma,mb);

 // 	SSEadd prod;
 // 	_mm_store_ps(prod,mprod);
 	
	// printf("Product is (%.2f,%.2f,%.2f,%2.f)\n",
	// 	prod[0],prod[1],prod[2],prod[3]);

 	// int64_t lol[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

 	// __m128i sum;
 	// for (int i = 0; i < 12; i+=4)
 	// {
 	// 	sum += _mm_loadu_si128(&lol[i]);
 	// }

 	// typedef int64_t sseADD[4];
 	// sseADD c = {100, 200, 300, 400};
 	// sseADD d = {500, 600, 700, 800};
 	// __m128i mc = _mm_load_si128(&c[0]);


 	// // printf("<%ld>\n", (int64_t)sum);

 	// int64_t j = 0;
 	// while (j < 12)
 	// {
 	// 	int temp_1 = _mm_load_si128()
 	// }

 	define_settings();
 	define_vectors();
 	compute_engine();

 	release();

 	return 0;
 }
