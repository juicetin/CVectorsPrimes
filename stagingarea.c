typedef struct {
	int64_t * vectorA, *vectorB, *vectorC;
	int64_t start, end, result, info1;
} wargs;

void * freq_worker (void *arg)
{
	wargs *counter = (wargs*) arg;
	int64_t count = 0;
	int64_t check = counter->info1;

	for (int i = counter->start; i < counter->end; ++i)
	{
		if (counter->vectorA[i] == check)
		{
			count++;
		}
	}
	counter->result = count;
	return NULL;
}

int64_t para_freq (int64_t *vector, int64_t to_check)
{
	wargs freq[g_nthreads];	// array to hold local thread info
	int64_t splitter = g_length/g_nthreads;	// array i separator

	// pass data to each thread to work on
	for (int i = 0; i < g_nthreads; ++i)
	{
		freq[i].vectorA = vector;
		freq[i].info1 = toCheck;
		freq[i].start = i * splitter;
		freq[i].end (i+1) * splitter;
	}
	freq[g_nthreads-1].end = g_length;

	pthread_t threads[g_nthreads];	// create threads

	// start all threads
	for (int i = 0; i < g_nthreads; ++i)
	{
		pthread_create(&threads[i], NULL, freq_worker, &freq[i]);
	}

	// wait until all threads are done
	for (int i = 0; i < g_nthreads; ++i)
	{
		pthread_join(threads[i], NULL);
	}

	// go through local frequency of each thread to get overall
	// frequency
	int frequency = 0;
	for (int i = 0; i < g_nthreads; ++i)
	{
		frequency += freq[i].result;
	}
	return frequency;
}


