#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>

/* create struct to store thread arguments */
struct thread_data {
  int tid;
};

/* create struct to store barrier parameters */
typedef struct {
  pthread_mutex_t	count_lock;		       // mutex semaphore for the barrier
  pthread_cond_t all_thrd_compelte;	   // condition variable for leaving
  int	count;	              	         // count of the number who have arrived
} mylib_barrier_t;

struct timespec start, finish;

/* initialize global variables */ 
time_t now;
char time_local[100];
int **board;                           //the main red/blue board used for computation
int num_thrds;                         //num_thrds = total number of threads 
int n;                                 //n = board size n x n 
int t;                                 //t = tile grid size t x t
int c;                                 //c = terminating threshold c%
int max_iters;                         //max_iters = maximum number of iterations allowed

int tile_size;                         //tile size = n/t x n/t 
int cells_in_tile;                     //total number of cells in one tile
int tile_column;                       //number of tile columns

pthread_mutex_t		count_lock;		       //mutex semaphore for the barrier
mylib_barrier_t barrier;               //barrier for thread synchronization
bool global_finished; 


/* initialize all functions */ 
int** board_init(int row, int column);
void redMovement(int start_row, int finish_row); 
void blueMovement(int start_column, int finish_column); 
bool determineConvergence(int tile_row_start, int tile_row_finish, int tid);
void *threadComputation(void *thread_arg);
void mylib_barrier_init(mylib_barrier_t *b);
void mylib_barrier(mylib_barrier_t *b, int num_thrds);
void mylib_barrier_destroy(mylib_barrier_t *b);
void printBoard(int** board, int start_row, int finish_row);
void selfCheck(int** board);

int main(int argc, char **argv) {
  
  /* initialize global timestamp */
  now = time (0);
  strftime (time_local, 100, "%Y-%m-%d-%H:%M:%S", localtime (&now)); 

  /* user input preliminary check */
  if(argc < 6 || argc > 6 ) {
    printf("[Error]5 arguments expected!!\n");
    printf("Input Format: ./redblue_movement_v2_0_pthread <# of threads> <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: ./redblue_movement_v2_0_pthread 4 30 5 60 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
   }
  else if(atoi(argv[2])%atoi(argv[3])!=0){
    printf("[Error]!!Grid size n must be divisible by the tile size t!!\n");
    printf("SAMPLE: ./redblue_movement_v2_0_pthread 4 30 5 60 1000\n");
    exit(1);
  }
  else if(atoi(argv[1])<0||atoi(argv[2])<0||atoi(argv[3])<0||atoi(argv[4])<0||atoi(argv[5])<0||atoi(argv[4])>100){
    printf("[Error]!!Do NOT input negative numbers or threshold c larger than 100!!\n");
    printf("Input Format: ./redblue_movement_v2_0_pthread <# of threads> <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: ./redblue_movement_v2_0_pthread 4 30 5 60 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
  }

  num_thrds = atoi(argv[1]);           //number of threads nthrds
  n = atoi(argv[2]);                   //cell grid size n
  t = atoi(argv[3]);                   //tile grid size t
  c = atoi(argv[4]);                   //terminating threshold c
  max_iters = atoi(argv[5]);           //maximum number of iterations max_iters
  
  tile_size = n/t;                     //each tile is of size tile_size x tile_size cells
  cells_in_tile = n/t * n/t;
  tile_column = t;
  
  int i;
  pthread_t *thread_id;
  struct thread_data *thread_data_array; 
  mylib_barrier_init(&barrier);

  /* allocate memory for thread_ids and thread_data_array */
  thread_id = (pthread_t *)malloc(sizeof(pthread_t) * num_thrds);
  thread_data_array = (struct thread_data *) malloc (sizeof(struct thread_data) * num_thrds);  

  /* initialize board for computation and make a copy of original state for self checcking */
  board = board_init(n,n);
  int *board_scheck_cell = (int *)malloc(n * n * sizeof(int));
  int **board_scheck = (int **)malloc(n * sizeof(int*));
  /* assign the row pointer to the right place on the board  */
  for(i=0; i<n; i++) {
    board_scheck[i] = &board_scheck_cell[i*n];
  }
  /* copy values of the original board  */
  for(i=0; i<n; i++) {
    memcpy(board_scheck[i], board[i], sizeof(int) * n);
  }

  printf("\n=================================================\n");
  printf("The board size：n = %d (%dX%d)\n", n, n, n);
  printf("The tile grid size: t = %d (%dX%d)\n", t, t, t);
  printf("The terminating threshold: c = %d percent\n", c);
  printf("The maximum number of iterations: max_iters = %d\n", max_iters);
  printf("Current Local Time: %s\n", time_local);
  printf("=================================================\n");
  
  /* single thread messages */
  if (num_thrds==1) {

    printf("\n[Process #main] %s [single-thread]: INITIAL STATE of the board:\n\n", time_local);
    printBoard(board,0,n);
    printf("\n=================================================\n");

    printf("[Process #main] %s [single-thread]: sequential computation starting...\n", time_local);
    
    clock_gettime(CLOCK_REALTIME, &start);
    thread_data_array[0].tid = 0;
    pthread_create(&thread_id[0], NULL, threadComputation, &thread_data_array[0]);
    
    /* Wait for all threads to complete then print result */ 
    for (i=0; i < 1; i++) {
      pthread_join(thread_id[i],NULL);
    }
    
    clock_gettime(CLOCK_REALTIME, &finish);
    double time_spent_single = finish.tv_sec - start.tv_sec;
    time_spent_single += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("[Process #main] %s [single-thread]: FINAL STATE of the board:\n", time_local);
    printBoard(board,0,n);
    printf("[Process #main] %s [single-thread]: computation convergence time [%f]\n", time_local, time_spent_single);
  }
  
  else {

    printf("\n[Process #main] %s [multithreading]: INITIAL STATE of the board:\n\n", time_local);
    printBoard(board,0,n);
    printf("\n=================================================\n");

    printf("[Process #main] %s [multithreading]: multi threading computation starting...\n", time_local);
    
    clock_gettime(CLOCK_REALTIME, &start);
    
    /* create the threads, then wait for them to finish */
    for (i = 0; i < num_thrds; i++){
      thread_data_array[i].tid = i;
      pthread_create(&thread_id[i], NULL, threadComputation, &thread_data_array[i]);
    }

    /* Wait for all threads to complete then print global sum */ 
    for (i=0; i < 1; i++) {
      pthread_join(thread_id[i],NULL);
    }

    /* calculate time spent in total for multithread computation */ 
    clock_gettime(CLOCK_REALTIME, &finish);
    double time_spent_multi = finish.tv_sec - start.tv_sec;
    time_spent_multi += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    //sequentialComputation(board,n,t,c,max_iters);
    printf("[Process #main] %s [multithreading]: FINAL STATE of the board:\n", time_local);
    printBoard(board,0,n);
    printf("[Process #main] %s [multithreading]: computation convergence time [%f]\n", time_local, time_spent_multi);
    
    /* clean up the barrier */
    mylib_barrier_destroy(&barrier); /* destroy barrier object */
    
    /* Self checking serialization computation */
    printf("[Process #main] %s [multithreading][check]: self check starting in 9 seconds...\n", time_local);
    for (i=9;i>0;i--){
       printf("[Process #main] %s [multithreading][check]: %d...\r", time_local,i);
       sleep(1);
      fflush(stdout);
    }
    selfCheck(board_scheck);

  }
  pthread_exit (NULL);
}

/* this function is used to randomly generat a 2-d matrix for computation */ 
int** board_init(int row, int column) {

  /* allocate memory block for all the cells on the board */
  int* board_cell = (int *)malloc(row * column * sizeof(int));
  /* allocate memory block for the 1st pointers of each row */
  int** board =  (int **)malloc(row * sizeof(int*));
  if(board == NULL){
    fprintf(stderr, "**board out of memory\n");
    exit(1);
  }

  /* assign the row pointer to the right place on the board  */
  for(int i=0; i<row; i++) {
    board[i] = &board_cell[i*column];
  }

  /* set current time as seed for random generator */
  time_t t;
	srand((unsigned) time(&t));
  for (int i = 0; i < row; i++)
     for (int j = 0; j < column; j++)
        *(board[i] + j) = rand()%3;
  return board;
}
/* this function is used to compute horizontall red movement */ 
void redMovement(int start_row, int finish_row){
  
  int i,j;

  for (i = start_row; i < finish_row; i++){                   //row loop
    if (board[i][0] == 1 && board[i][1] == 0){                //check edge case where the first cell is red and can move LEFT, so the last red does NOT move.
      board[i][0] = 4;
      board[i][1] = 3;
    }
    for (j = 1; j < n; j++){                                  //column loop
      if (board[i][j] == 1 && board[i][(j+1)%n] == 0){        //when red can move right; (j+1)%n ensures if j = <last column>, the update wrap around to 1st column;
        board[i][j] = 0;
        board[i][(j+1)%n] = 3;
      }
      else if (board[i][j] == 3) board[i][j] = 1;
    }
    if (board[i][0] == 3) board[i][0] = 1;                    //2nd time to check first cell in this row iteration
    else if (board[i][0] == 4) board[i][0] = 0;
  }
}
/* this function is used to compute vertically blue movement */
void blueMovement(int start_column, int finish_column){

  int i,j;

  for (j = start_column; j < finish_column; j++){             //column loop
    if (board[0][j] == 2 && board[1][j] == 0){                //check edge case where the first cell is blue and can move DOWN, so the bottom blue does NOT move.
      board[0][j] = 4;
      board[1][j] = 3;
    }
    for (i = 1; i < n; i++){                                  //row loop
      if (board[i][j] == 2 && board[(i+1)%n][j] == 0){        //when blue can move down wraparound; (i+1)%n ensures if i = <last row>, the pointer moves back to first row;
        board[i][j] = 0;
        board[(i+1)%n][j] = 3;
      }
      else if (board[i][j] == 3) board[i][j] = 2;
    }
    if (board[0][j] == 3) board[0][j] = 2;                    //2nd time to check first cell in this column iteration
    else if (board[0][j] == 4 ) board[0][j] = 0;
  }
}
/* this function is used to analyze convergence by calculating red/blue percentage in each tile */ 
bool determineConvergence(int tile_row_start, int tile_row_finish, int tid){

  int red_1_count, blue_2_count;
  double red_percentage, blue_percentage;
  bool local_finished;
  int t_r, t_c, c_r, c_c;

  /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
  // Loop structure:
  //-tile_row(t_r) loop
  //   -tile_column(t_c) loop
  //         -cell_row(c_r) loop
  //              -cell_column(c_c) loop

  for (t_r = tile_row_start; t_r < tile_row_finish; t_r++) {
    for (t_c = 0; t_c < tile_column; t_c++) {
      /* initialize and reset red/blue count to Zero for each tile*/
      red_1_count = 0;
      blue_2_count = 0;
      for (c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {        //(t_r+1)*tile_size is the start row of the below tile
        for (c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){       //(t_c+1)*tile_size is the start column of the tile on the right
          if (board[c_r][c_c] == 1) {
            red_1_count+=1;                                                //if cell value is 1, red count plus 1;
          }
          else if (board[c_r][c_c] == 2){
            blue_2_count+=1;                                               //if cell value is 2, blue count plus 1;
          } 
        }
      }

      /* calculate color percentage and compare with threshold to decide whether to terminate */
      red_percentage = ((double)red_1_count/cells_in_tile)*100;
      blue_percentage = ((double)blue_2_count/cells_in_tile)*100;
      
      if (red_percentage > c) {
        local_finished = true;
        printf("[Thread #%d] %s [multithreading]: converged RED with %d red cells ==> %.2f percent ==> converged tile:\n", tid, time_local, red_1_count, red_percentage);
        /* print out the converged tile */
        for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
          for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
            printf("%d ", board[c_r][c_c]);
          }
          printf("\n");
        }
      }
      else if (blue_percentage > c) {
        local_finished = true;
        printf("[Thread #%d] %s [multithreading]: converged BLUE with %d blue cells ==> %.2f percent ==> converged tile:\n",tid, time_local, blue_2_count, blue_percentage);
        /* print out the converged tile */
        for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
          for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
            printf("%d ", board[c_r][c_c]);
          }
          printf("\n");
        }
      }
    }
  }
  return local_finished;
}
/* this is the main workload of each thread: 
   1. calculate the portion of its assigned sub-matrix; 
   2. call redMovement
   3. call blueMovement 
   4. call determineConvergence
   5. update global_finished if any thread converges and therefore terminate all threads 
*/
void *threadComputation(void *thread_arg){

  struct thread_data *t_data; 
  bool local_finished;
  int n_itrs = 0; 

  pthread_mutex_t	count_lock;
  t_data = (struct thread_data *)thread_arg;
  int tid = t_data->tid;                   // get my thread id from thread data

  int q = t / num_thrds;                   //minimum tile numbers each process handles
  int r = t % num_thrds;                   //tile remainders
  int ib, kn, K; 

  /* determine the portion of board for this thread */
  /* since tile unit is considered, multiply row pointers by tize_size n/t */
  if (tid < r){                            //evenly allocate remainder tiles when r != 0
    ib = tid * (q+1) * n/t;                //ib is the first row of each sub-board assigned to the thread, therefore q multiply by tile size
    K = (q+1) * n/t;
  }
  else{
    ib = tid * q * n/t + r * n/t;
    K = q * n/t;                           //K is the number of rows in this sub-board
  }
  kn = K * n;                              //K * N is the size of the assigned sub-board

  int tile_row_start = ib/tile_size;
  int tile_row_finish = (ib+K)/tile_size;
  
  while (!(global_finished) && n_itrs < max_iters){
    n_itrs++;
    redMovement(ib, ib+K);
    mylib_barrier(&barrier,num_thrds);
    blueMovement(ib, ib+K);
    mylib_barrier(&barrier,num_thrds);
    local_finished = determineConvergence(tile_row_start, tile_row_finish, tid);
    if (local_finished) global_finished = true;
    mylib_barrier(&barrier,num_thrds);
    if (global_finished) {
      printf("[Thread #%d] %s [multithreading]: terminated after [%d] iterations\n", tid, time_local, n_itrs);
      break;
    }
  }
  printf("[Thread #%d] %s [multithreading]: terminated after [%d] iterations\n", tid, time_local, n_itrs);
}
/* this function is used to define the barrier */
void mylib_barrier_init(mylib_barrier_t *b){
  b->count = 0;
  pthread_mutex_init(&(b->count_lock), NULL);
  pthread_cond_init(&(b->all_thrd_compelte), NULL);
}
/* this fucntion is used to call the barrier */
void mylib_barrier(mylib_barrier_t *b, int num_thrds)
{
  pthread_mutex_lock(&(b->count_lock));
  b->count++; 
  if (b->count == num_thrds) {
    b->count = 0; 
    pthread_cond_broadcast(&(b->all_thrd_compelte));
  } 
  else {
    while (pthread_cond_wait(&(b->all_thrd_compelte), &(b->count_lock)) != 0);  
  }
  pthread_mutex_unlock(&(b->count_lock));
}
/* this function is used to clean the barrier */ 
void mylib_barrier_destroy(mylib_barrier_t *b) {
 
  pthread_mutex_destroy(&(b->count_lock));
  pthread_cond_destroy(&(b->all_thrd_compelte));

}
/* this function is used to print any given matrix */ 
void printBoard(int** board, int start_row, int finish_row){
  
  int i,j; 
  for (int i = start_row; i < finish_row; i++) {
    for (int j = 0; j < n; j++) {
      printf("%d ", *(board[i] + j));
    }
    printf("\n");
  }
}
/* this function is used to do self check by sequential computation */
void selfCheck(int** board_scheck){

  bool finished_check = false;
  int n_itrs_check = 0;
  int red_1_count_check, blue_2_count_check;
  int tile_row = n/(n/t);
  int i,j;

  printf("[Process #main] %s [multithreading][check]: self checking computation in progress...\n", time_local);

  clock_gettime(CLOCK_REALTIME, &start);
  while (!finished_check && n_itrs_check < max_iters){
    n_itrs_check++;
    /***** Stage 1: Red Movement ******/
    for (int i = 0; i < n; i++){
      if (board_scheck[i][0] == 1 && board_scheck[i][1] == 0){
        board_scheck[i][0] = 4;
        board_scheck[i][1] = 3;
      }
      for (int j = 1; j < n; j++){
        if (board_scheck[i][j] == 1 && board_scheck[i][(j+1)%n] == 0){
        board_scheck[i][j] = 0;
        board_scheck[i][(j+1)%n] = 3;
        }
        else if (board_scheck[i][j] == 3) board_scheck[i][j] = 1;
      }
      if (board_scheck[i][0] == 3) board_scheck[i][0] = 1;
      else if (board_scheck[i][0] == 4) board_scheck[i][0] = 0;
    }
    /***** Stage 2: Blue Movement ******/
    for (int j = 0; j < n; j++){
      if (board_scheck[0][j] == 2 && board_scheck[1][j] == 0){
      board_scheck[0][j] = 4;
      board_scheck[1][j] = 3;
    }
    for (int i = 1; i < n; i++){
      if (board_scheck[i][j] == 2 && board_scheck[(i+1)%n][j] == 0){
        board_scheck[i][j] = 0;
        board_scheck[(i+1)%n][j] = 3;
      }
      else if (board_scheck[i][j] == 3) board_scheck[i][j] = 2;
    }
    if (board_scheck[0][j] == 3) board_scheck[0][j] = 2;
    else if (board_scheck[0][j] == 4 ) board_scheck[0][j] = 0;
  }
    /***** Stage 3: Determine if the computation has converged ******/
    for (int t_r = 0; t_r < tile_row; t_r++) {
      for (int t_c = 0; t_c < tile_column; t_c++) {

        /* initialize and reset red/blue count to Zero for each tile*/
        red_1_count_check = 0;
        blue_2_count_check = 0;

        for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
          for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){

            if (board_scheck[c_r][c_c] == 1) {
              red_1_count_check+=1;
            }
            else if (board_scheck[c_r][c_c] == 2){
              blue_2_count_check+=1;
            }
          }
        }

        /* calculate color percentage and compare with threshold to decide whether to terminate */
        double red_percentage_check = ((double)red_1_count_check/cells_in_tile)*100;
        double blue_percentage_check = ((double)blue_2_count_check/cells_in_tile)*100;

        if (red_percentage_check > c) {
          finished_check = true;
          printf("[Process #main] %s [multithreading][check]: converged RED with %d red cells ==> [%.2f] ==> converged tile:\n", time_local, red_1_count_check, red_percentage_check);
          /* print out the converged tile */
          for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
            for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
              printf("%d ", board_scheck[c_r][c_c]);
            }
            printf("\n");
          }
        }
        else if (blue_percentage_check > c) {
          finished_check = true;
          printf("[Process #main] %s [multithreading][check]: converged BLUE with %d blue cells ==> [%.2f] ==> converged tile:\n", time_local, blue_2_count_check, blue_percentage_check);
          /* print out the converged tile */
          for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
            for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
              printf("%d ", board_scheck[c_r][c_c]);
            }
            printf("\n");
          }
        }
      }
    }
  }
  clock_gettime(CLOCK_REALTIME, &finish);
  double time_spent_check = finish.tv_sec - start.tv_sec;
  time_spent_check += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  printf("[Process #main] %s [multithreading][check]: self checking completed after [%d] iterations. FINAL SELF CHECK STATE of the board_scheck:\n", time_local,n_itrs_check);
  printBoard(board_scheck,0,n);
  
  //validate each cell and count the difference
  int cell_dif = 0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      if(board_scheck[i][j]!=board[i][j]) {
        cell_dif++;
      }
    }
  }
  if (cell_dif == 0){
    printf("[Process #main] %s [multithreading][check]: multithreading computation matches self checking sequential computation result!!\n", time_local);
  }
  else {
    printf("[Process #main] %s [multithreading][check]: Oops...there are [%d] cells different between multithreading computation and sequential self checking computation result...\n", time_local, cell_dif);
  }
  printf("[Process #main] %s [multithreading][check]: sequential self checking computation convergence time [%f]\n", time_local, time_spent_check);
}