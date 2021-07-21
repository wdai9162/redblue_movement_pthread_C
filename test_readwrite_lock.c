#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#define NUM_THREADS 5 
 
typedef struct {
  int reader_count;                              //记录有多少个reader_lock
  int writer;                                    //0或者1来记录是否有writer存在
  int pending_writer;                            //记录有多少个writer在排队
  pthread_cond_t cond_reader_can_proceed;        //condition variable来针对reader可以访问的信号
  pthread_cond_t cond_writer_can_proceed;        //condition variable来针对writer可以访问的信号
  pthread_mutex_t read_write_lock;               //condition varaible和共享的struct的互斥访问锁
} mylib_rwlock_t; 

//定义一个函数来初始化这个struct
void mylib_rwlock_init(mylib_rwlock_t *rwlock){
  rwlock->reader_count=0;
  rwlock->writer=0;
  rwlock->pending_writer=0;
  pthread_mutex_init(&(rwlock->read_write_lock),NULL);
  pthread_cond_init(&(rwlock->cond_reader_can_proceed),NULL);
  pthread_cond_init(&(rwlock->cond_writer_can_proceed),NULL);
}

//定义reader lock的函数
void mylib_rwlock_rlock(mylib_rwlock_t *rwlock){
  //如果已经存在writer lock, 所有reader lock都得执行condition wait等待写入结束
  //如果不存在，那么授权reader lock并增加reader count
  pthread_mutex_lock(&(rwlock->read_write_lock));
  while(rwlock->writer>0||rwlock->pending_writer>0){
    pthread_cond_wait(&(rwlock->cond_reader_can_proceed),&(rwlock->read_write_lock));
  }
  //read访问授权成功，可以read数据
  rwlock->reader_count++;
  pthread_mutex_unlock(&(rwlock->read_write_lock));
}

//定义writer lock的函数
void mylib_rwlock_wlock(mylib_rwlock_t *rwlock){
  //如果已经存在reader或者writer lock, 增加pending_writer的数量and wait
  //被signal唤醒以后，减少pending_writer的数量并且增加writer count
  pthread_mutex_lock(&(rwlock->read_write_lock));
  rwlock->pending_writer++;
  while((rwlock->writer>0)||(rwlock->reader_count>0)){
    pthread_cond_wait(&(rwlock->cond_writer_can_proceed),&(rwlock->read_write_lock));
  }
  rwlock->pending_writer--; 
  rwlock->writer++;
  pthread_mutex_unlock(&(rwlock->read_write_lock));
}

//定义reader unlock函数
void mylib_rwlock_runlock(mylib_rwlock_t *rwlock){
  //减少reader数量，如果reader count为0且pending_writer大于0
  //那就signal让writer访问
  pthread_mutex_lock(&(rwlock->read_write_lock));
  rwlock->reader_count--;
  if(rwlock->reader_count==0){
    pthread_cond_signal(&(rwlock->cond_writer_can_proceed));
  }
  pthread_mutex_unlock(&(rwlock->read_write_lock));
}

//定义writer unlock函数
void mylib_rwlock_wunlock(mylib_rwlock_t *rwlock){
  //减少/重置writer数量为0
  //如果有pending_writer > 0，那就signal让队列中的writer访问
  //如果没有pending writer，那就signal所有等待中的reader访问
  rwlock->writer=0;
  if(rwlock->pending_writer>0){
    pthread_cond_signal(&(rwlock->cond_writer_can_proceed));
  }
  else {
    pthread_cond_broadcast(&(rwlock->cond_reader_can_proceed));
  }
  pthread_mutex_unlock(&(rwlock->read_write_lock));
}

//定义rwlock destroy


int main(int argc, char *argv[])
{
  pthread_t threads[NUM_THREADS];
  int rc, t, tids[NUM_THREADS];
  struct thrd_args *args; 


  pthread_exit(NULL);

}
