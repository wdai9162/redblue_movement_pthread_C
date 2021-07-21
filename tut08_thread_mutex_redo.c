/******************************************************************************
* The exercise is to add a range of natural integer numbers from 0 to n
*   using a number of threads.  
******************************************************************************/
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

struct thrd_data{
	int id;
	int start;
	int end; /* the sub-range is from start to end-1 */
};

double  sum=0.0; /* global variable */
pthread_mutex_t sum_mutex; /* global mutex */


void *do_work(void *thrd_arg) 
{
	struct thrd_data *t_data;
	int i,start, end;
	int myid;
	double mysum=0.0;

	/* Initialize my part of the global array and keep local sum */
	t_data = (struct thrd_data *) thrd_arg; // cast the pointer type 
	start = t_data->start;
	end = t_data->end;
	myid = t_data->id;

	printf ("Thread %d doing sum from %d to %d\n", myid,start,end-1); 
	for (i=start; i < end ; i++) {
		mysum = mysum + i * 1.0;
    }

	/* Lock the mutex and update the global sum */
	pthread_mutex_lock(&sum_mutex);
	sum+=mysum;
	pthread_mutex_unlock(&sum_mutex);
	/* thread exit */
	pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
	int i, n, n_threads;
	int k, nq, nr;
	struct thrd_data *t_arg;
	pthread_t *thread_id;
	pthread_attr_t attr;

	/* Pthreads setup: initialize mutex and explicitly create */
	/* threads in a joinable state (for portability)  */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	/* initialize sum_mutex */
	pthread_mutex_init(&sum_mutex,NULL);
	
	/* ask to enter n and n_threads from the user */
	printf("Please enter n: \n");
	scanf("%d", &n);
	printf("Please enter n_threads:\n");
	scanf("%d", &n_threads);

	/* create arrays of thread_ids and thread t_args */
	thread_id = (pthread_t *)malloc(sizeof(pthread_t)*n_threads);
	t_arg = (struct thrd_data *)malloc(sizeof(struct thrd_data)*n_threads);

	/* distribute load and create threads for computation */
	
	nq = n/n_threads;
	nr = n%n_threads;
	k=1;
	for(i=0; i<n_threads; i++){
		t_arg[i].id=i;
		t_arg[i].start=k;
		if (i<nr) k=k+nq+1;
		else k = k+nq;
		t_arg[i].end=k;
		pthread_create(&thread_id[i], &attr, do_work, (void *)&t_arg[i]);
	}

	/* Wait for all threads to complete then print global sum */ 
	for (i=0; i<n_threads; i++) {
	  pthread_join(thread_id[i], NULL);
	}

	printf ("Done. Sum= %d \n", sum);

	sum=0.0;
	for (i=1;i<=n;i++)
		sum = sum + i * 1.0; 
	printf("Check Sum= %d\n",sum);

	/* Clean up and exit */
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&sum_mutex);
	pthread_exit (NULL);
}