#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#define NUM_THREADS 5 
 
struct thrd_args {
    int tid;
    char *message;
};


void *printHello(void *thrd_args){
    struct thrd_args *myarg;
    myarg = (struct thrd_args *) thrd_args;

    int my_tid = myarg->tid; 
    char *my_messsage = myarg->message;

    printf("%s", my_messsage);
    pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int rc, t, tids[NUM_THREADS];
  struct thrd_args *args; 
  args = (struct thrd_args *)malloc(NUM_THREADS*sizeof(struct thrd_args));

  for (t=0; t<NUM_THREADS; t++){
      printf("In main: creating thread %d\n", t);
      args[t].tid = t;
      args[t].message = "Hello World\n";
      rc = pthread_create(&threads[t], NULL, printHello, (void*)&args[t]);
      if(rc){
          printf("ERROR; return code from pthread_create() is %d\n", rc);
          exit(-1);
      }
  }

  pthread_exit(NULL);

}
