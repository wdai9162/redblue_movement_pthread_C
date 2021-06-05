/* Jacobi iteration using pthreads
   two functions - mylib_barrier() and mylib_barrier_destroy() in this program are missing and 
   you need to implement them.

   to compile:
     gcc -o jacobi jacobithreadSol.c -lpthread 
   to run (assuming numWorkers divides gridSize):
     jacobi gridSize numWorkers numIters 
*/

#include <pthread.h>
#include <stdio.h>
#define MAXGRID 258   /* maximum grid size, including boundaries */
#define MAXWORKERS 16  /* maximum number of worker threads */
#include <stdlib.h>

void *Worker(void *); 
void InitializeGrids();

typedef struct {
  pthread_mutex_t		count_lock;		/* mutex semaphore for the barrier */
  pthread_cond_t		ok_to_proceed;	/* condition variable for leaving */
  int				count;		/* count of the number who have arrived */
} mylib_barrier_t;

void mylib_barrier_init(mylib_barrier_t *b){
  b->count = 0;
  pthread_mutex_init(&(b->count_lock), NULL);
  pthread_cond_init(&(b->ok_to_proceed), NULL);
}

void mylib_barrier(mylib_barrier_t *b, int num_threads)
{
  pthread_mutex_lock(&(b->count_lock));
  b->count++; 
  if (b->count == num_threads) {
    b->count = 0; 
    pthread_cond_broadcast(&(b->ok_to_proceed));
  } 
  else {
    while (pthread_cond_wait(&(b->ok_to_proceed), &(b->count_lock)) != 0);  
  }
  pthread_mutex_unlock(&(b->count_lock));
}

void mylib_barrier_destroy(mylib_barrier_t *b) {		/* please implement this function */  
 
  pthread_mutex_destroy(&(b->count_lock));
  pthread_cond_destroy(&(b->ok_to_proceed));

}

int gridSize, numWorkers, numIters, stripSize;
double maxDiff[MAXWORKERS];
double grid1[MAXGRID][MAXGRID], grid2[MAXGRID][MAXGRID];

/* barrier */
mylib_barrier_t	barrier;

/* main() -- read command line, initialize grids, and create threads
             when the threads are done, print the results */

int main(int argc, char *argv[]) {
  /* thread ids and attributes */
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;

  int i, j, tids[MAXWORKERS];
  double maxdiff = 0.0;
  FILE *results;

 
  /* make threads joinable */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* initialize barrier */
  mylib_barrier_init(&barrier);

  /* read command line and initialize grids */
  gridSize = atoi(argv[1]);
  numWorkers = atoi(argv[2]);
  numIters = atoi(argv[3]);
  stripSize = gridSize/numWorkers;
  InitializeGrids();

  /* create the workers, then wait for them to finish */
  for (i = 0; i < numWorkers; i++){
    tids[i] = i;
    pthread_create(&workerid[i], NULL, Worker, (void *) &tids[i]);
  }

  for (i = 0; i < numWorkers; i++)
    pthread_join(workerid[i], NULL);

  /* print the results */
  for (i = 0; i < numWorkers; i++)
    if (maxdiff < maxDiff[i])
      maxdiff = maxDiff[i];
  printf("number of iterations:  %d\nmaximum difference:  %e\n",
          numIters, maxdiff);

  results = fopen("results", "w");
  for (i = 1; i <= gridSize; i++) {
    for (j = 1; j <= gridSize; j++) {
      fprintf(results, "%f ", grid1[i][j]);
    }
   fprintf(results, "\n");
  }


  /* Clean up and exit */
  pthread_attr_destroy(&attr);
  mylib_barrier_destroy(&barrier); /* destroy barrier object */
  pthread_exit (NULL);

}


/* Each Worker computes values in one strip of the grids.
   The main worker loop does two computations to avoid copying from
   one grid to the other.  */

void *Worker(void *tid) {
  int *myid;
  double maxdiff, temp;
  int i, j, iters;
  int first, last;

  /* determine first and last rows of my strip of the grids */
  myid = (int *) tid;
  first = (*myid)*stripSize + 1;
  last = first + stripSize - 1;

  printf("myid = %d, first = %d, last = %d \n", *myid, first, last);

  for (iters = 1; iters <= numIters; iters++) {
    /* update my points */
    for (i = first; i <= last; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid2[i][j] = (grid1[i-1][j] + grid1[i+1][j] + 
                       grid1[i][j-1] + grid1[i][j+1]) * 0.25;
      }
    }

    mylib_barrier(&barrier, numWorkers);

    /* update my points again */
    for (i = first; i <= last; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid1[i][j] = (grid2[i-1][j] + grid2[i+1][j] +
               grid2[i][j-1] + grid2[i][j+1]) * 0.25;
      }
    }

    mylib_barrier(&barrier, numWorkers);

  }
  /* compute the maximum difference in my strip and set global variable */
  maxdiff = 0.0;
  for (i = first; i <= last; i++) {
    for (j = 1; j <= gridSize; j++) {
      temp = grid1[i][j]-grid2[i][j];
      if (temp < 0)
        temp = -temp;
      if (maxdiff < temp)
        maxdiff = temp;
    }
  }
  maxDiff[*myid] = maxdiff;

  /* thread exit */
  pthread_exit(NULL);
}

void InitializeGrids() {
  /* initialize the grids (grid1 and grid2)
     set boundaries to 1.0 and interior points to 0.0  */
  int i, j;
  for (i = 0; i <= gridSize+1; i++)
    for (j = 0; j <= gridSize+1; j++) {
      grid1[i][j] = 0.0;
      grid2[i][j] = 0.0;
    }
  for (i = 0; i <= gridSize+1; i++) {
    grid1[i][0] = 1.0;
    grid1[i][gridSize+1] = 1.0;
    grid2[i][0] = 1.0;
    grid2[i][gridSize+1] = 1.0;
  }
  for (j = 0; j <= gridSize+1; j++) {
    grid1[0][j] = 1.0;
    grid2[0][j] = 1.0;
    grid1[gridSize+1][j] = 1.0;
    grid2[gridSize+1][j] = 1.0;
  }
}
