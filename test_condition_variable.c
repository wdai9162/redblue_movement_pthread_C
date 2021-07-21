#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#define NUM_THREADS 5 
 
pthread_cond_t cond_queue_is_empty, cond_queue_is_full;    //分别两个condition variable代表队列空和队列满两种条件
pthread_mutex_t task_queue_cond_lock;                      //一个lock对应queue的状态
int task_available; 

struct thrd_args {
    int tid;
    char *message;
};


void *producer(void *thrd_args){
    struct thrd_args *myarg;
    myarg = (struct thrd_args *) thrd_args;
    int my_tid = myarg->tid; 

    while(!done()){
      create_task();
      pthread_mutex_lock(&task_queue_cond_lock);
      if(task_available==1){
        pthread_cond_wait(&cond_queue_is_empty, &task_queue_cond_lock);   //这里面mutex task_queue_cond_lock就暂时释放掉了，这个thread等待condition的时候，这个mutex就可以被其他thread获取了
        //直到在这个cond被其他thread signal以后，这个mutex重新还给该thread，该thread重新lock mutex然后继续运算
      }
      insert_into_queue();
      task_available==1;
      pthread_cond_signal(&cond_queue_is_full);    //这里signal 在等待这个cond_queue_is_full的thread
      pthread_mutex_unlock(&task_queue_cond_lock); //释放mutex，thread运算结束
    pthread_exit(NULL);
}

void *consumer(void *consumer_thread_data){
  while(!done(){
    pthread_mutex_lock(&task_queue_cond_lock);
    if(task_available==0){
      pthread_cond_wait(&cond_queue_is_full,&task_queue_cond_lock);
    }
    my_task=extract_from_queue();
    task_available=0;
    pthread_cond_signal(cond_queue_is_empty);
    pthread_mutex_unlock(&task_queue_cond_lock);
    process_task(my_task);
  }
}

int main(int argc, char *argv[])
{ 
  /* declarations and initializations */
  task_available=0;
  pthread_cond_init(&cond_queue_is_empty, NULL);
  pthread_cond_init(&cond_queue_is_full, NULL);
  pthread_mutex_init(&task_queue_cond_lock, NULL);

  /* create and join producer and consumer threads */ 


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
