#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

typedef struct {
 int count; /* the integer variable */
 pthread_mutex_t s_m; /* mutex */
 pthread_cond_t s_cv; /* condition variable */
} mylib_sem_t;

mylib_sem_init(mylib_sem_t *s, int value){
    s->count=0;
    pthread_mutex_init(&(s->s_m),NULL);
    pthread_cond_init(&(s->s_cv),NULL);
}

mylib_sem_wait(mylib_sem_t *s){
    pthread_mutex_lock(&(s->s_m));
    s->count--;
    if(s->count<0){
        pthread_cond_wait(&(s->s_cv), &(s->s_m));
    }
    pthread_mutex_unlock(&(s->s_m));
}

mylib_sem_signal(mylib_sem_t *s){
    pthread_mutex_lock(&(s->s_m));
    s->count++;
    if(s->count<=0){
        pthread_cond_signal(&(s->s_cv));
    }
    pthread_mutex_unlock(&(s->s_m));
}


int main(void)
{

}