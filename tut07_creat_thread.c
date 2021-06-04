#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#define NUM_THREADS	8

char *messages[NUM_THREADS];

struct thread_data
{
	int	thread_id;
	char *message;
};

struct thread_data thread_data_array[NUM_THREADS];

void *PrintHello(void *threadarg)
{
	int taskid;
	char *hello_msg;
	struct thread_data *my_data;
    my_data = threadarg;
	sleep(1);

	/* retrieve message and print the id and the message. */
	printf("%s\n",my_data->message);
	pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
	pthread_t threads[NUM_THREADS];
	int rc, t;

	messages[0] = "English: Hello World!";
	messages[1] = "French: Bonjour, le monde!";
	messages[2] = "Spanish: Hola al mundo";
	messages[3] = "Klingon: Nuq neH!";
	messages[4] = "German: Guten Tag, Welt!"; 
	messages[5] = "Russian: Zdravstvytye, mir!";
	messages[6] = "Japan: Sekai e konnichiwa!";
	messages[7] = "Latin: Orbis, te saluto!";

	for(t=0;t<NUM_THREADS;t++) {
		thread_data_array[t].thread_id = t;
		thread_data_array[t].message = messages[t];

		/* create a thread */
		printf("Creating thread %d\n", t);
		rc = pthread_create(&threads[t],NULL,PrintHello,&thread_data_array[t]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	pthread_exit(NULL);
}