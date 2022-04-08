#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#define NUM_THREADS 5 
 
typedef struct {
  pthread_mutex_t thrd_count_lock;
  pthread_cond_t ok_to_proceed;
  int thrd_count; 
} mylib_barrier_t; 

//定义初始化函数
void mylib_barrier_init(mylib_barrier_t *b){
  b->thrd_count = 0;
  pthread_mutex_init(&(b->thrd_count_lock),NULL);
  pthread_cond_init(&(b->ok_to_proceed),NULL);
}

void mylib_barrier(mylib_barrier_t *b, int num_threads){
  //每一个call到这个barrier的thread，那么thrd_count加一表示该函数抵达barrier。该thread立即对比thrd_count和NUM_THREADS，如果小于NUM_THREADS进入等待状态
  //如果thrd_count等于NUM_THREADS，表示所有的thread都执行到了这个步奏，最后到达的thread那就signal ok_to_proceed
  
  pthread_mutex_lock(&(b->thrd_count_lock));
  b->thrd_count++;
  if(b->thrd_count<num_threads){
    pthread_cond_wait(&(b->ok_to_proceed), &(b->thrd_count_lock));
  }
  else{
    pthread_cond_signal(&(b->ok_to_proceed));
  }
  pthread_mutex_unlock(&(b->thrd_count_lock));
}

//定义销毁函数
void mylib_barrier_destory(mylib_barrier_t *b){

}



int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int rc, t, tids[NUM_THREADS];
  struct thrd_args *args; 


  pthread_exit(NULL);

}
