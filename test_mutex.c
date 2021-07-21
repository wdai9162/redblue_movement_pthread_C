#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#define NUM_THREADS 50

//mutex在global scope下进行声明 --> 在main function里进行初始化
pthread_mutex_t minimum_value_lock;
int best_cost; 


struct thrd_args {
    int tid;
    int my_cost;
};



void *findMin(void *thrd_args){
    struct thrd_args *myarg;
    int my_cost; 
    myarg = (struct thrd_args *) thrd_args;

    my_cost = myarg->my_cost;
    pthread_mutex_lock(&minimum_value_lock);
    if(my_cost>best_cost){
      best_cost = my_cost;
      } 
    pthread_mutex_unlock(&minimum_value_lock);
    printf("ending best_cost = %d\n", best_cost);
    pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
  //mutex在main function下进行初始化 --> 在global scope里声明
  pthread_mutex_init(&minimum_value_lock,NULL);

  pthread_t threads[NUM_THREADS];
  int rc, t, tids[NUM_THREADS];
  struct thrd_args *args; 
  args = (struct thrd_args *)malloc(NUM_THREADS*sizeof(struct thrd_args));

  for (t=0; t<NUM_THREADS; t++){
      printf("Main thread: creating thread %d\n", t);
      args[t].tid = t;
      args[t].my_cost = ((t+12)*(t+23)+283%(t+2))/3;
      rc = pthread_create(&threads[t], NULL, findMin, (void*)&args[t]);
      if(rc){
          printf("ERROR; return code from pthread_create() is %d\n", rc);
          exit(-1);
      }
  }

  /* Wait for all threads to complete then print global sum */
  for (int i=0; i < NUM_THREADS; i++) {
    pthread_join(threads[i],NULL);
  }

  printf("best cost is: %d \n", best_cost);
  pthread_exit(NULL);

}
