#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

char time_local[100];
time_t now;

/* create thread argument struct for redMove() and blueMove() */
struct thread_data {
  int tid;
  int cell_start;
  int cell_end;
  int grid_start;
  int grid_end; 
  int **board;
  int n;
  int t;
  int max_iters;
  int c;
};

/*define a function to initialize the board as required */
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

void redMovement(int** board, int n){
  
  int i,j;

  for (i = 0; i < n; i++){                                            //row loop
    if (board[i][0] == 1 && board[i][1] == 0){                //check edge case where the first cell is red and can move LEFT, so the last red does NOT move.
      board[i][0] = 4;
      board[i][1] = 3;
    }
    for (j = 1; j < n; j++){                                          //column loop
      if (board[i][j] == 1 && board[i][(j+1)%n] == 0){        //when red can move right; (j+1)%n ensures if j = <last column>, the pointer moves back to first column;
      board[i][j] = 0;
      board[i][(j+1)%n] = 3;
      }
      else if (board[i][j] == 3) board[i][j] = 1;
    }
    if (board[i][0] == 3) board[i][0] = 1;                    //2nd time to check first cell in this row iteration
    else if (board[i][0] == 4) board[i][0] = 0;
  }
}

void blueMovement(int** board, int n){

  int i,j;

  for (j = 0; j < n; j++){                                            //column loop
    if (board[0][j] == 2 && board[1][j] == 0){                //check edge case where the first cell is blue and can move DOWN, so the bottom blue does NOT move.
      board[0][j] = 4;
      board[1][j] = 3;
    }
    for (i = 1; i < n; i++){                                            //row loop
      if (board[i][j] == 2 && board[(i+1)%n][j] == 0){          //when blue can move down wraparound; (i+1)%n ensures if i = <last row>, the pointer moves back to first row;
        board[i][j] = 0;
        board[(i+1)%n][j] = 3;
      }
      else if (board[i][j] == 3) board[i][j] = 2;
    }
    if (board[0][j] == 3) board[0][j] = 2;                  //2nd time to check first cell in this column iteration
    else if (board[0][j] == 4 ) board[0][j] = 0;
  }
}

bool determineConvergence(int** board, int n, int tile_size, int tile_row, int tile_column, int cells_in_tile, int c , int tid){

  /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
  // Loop structure:
  //-tile_row(t_r) loop
  //   -tile_column(t_c) loop
  //         -cell_row(c_r) loop
  //              -cell_column(c_c) loop
  int red_1_count, blue_2_count;
  double red_percentage, blue_percentage;
  int t_r, t_c, c_r, c_c;
  bool local_finished;
  
  for (t_r = 0; t_r < tile_row; t_r++) {
    for (t_c = 0; t_c < tile_column; t_c++) {
      /* initialize and reset red/blue count to Zero for each tile*/
      red_1_count = 0;
      blue_2_count = 0;
      for (c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {    //(t_r+1)*t is the start row of the below tile
        for (c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){   //(t_c+1)*t is the start column of the tile on the right
          //printf("%d ", *(board[c_r]+c_c));          //Print Tiles
          if (board[c_r][c_c] == 1) {
            red_1_count+=1;
          }         //if cell value is 1, red count plus 1;
          else if (board[c_r][c_c] == 2){
            blue_2_count+=1;
          }        //if cell value is 2, blue count plus 1;
        }
      }

      /* calculate color percentage and compare with threshold to decide whether to terminate */
      red_percentage = ((double)red_1_count/cells_in_tile)*100;
      blue_percentage = ((double)blue_2_count/cells_in_tile)*100;
      
      if (red_percentage > c) {
        local_finished = true;
        printf("[Thread #%d] %s [serial]: converged RED with %d red cells ==> %.2f percent ==> converged tile:\n", tid, time_local, red_1_count, red_percentage);
        /* print out the converged tile */
        for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
          for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
            printf("%d ", board[c_r][c_c]);
          }
          printf("\n");
        }
        break;
      }
      else if (blue_percentage > c) {
        local_finished = true;
        printf("[Thread #%d] %s [serial]: converged BLUE with %d blue cells ==> %.2f percent ==> converged tile:\n",tid, time_local, blue_2_count, blue_percentage);
        /* print out the converged tile */
        for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
          for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
            printf("%d ", board[c_r][c_c]);
          }
          printf("\n");
        }
        break;
      }
    }
    if (local_finished == true); break;
  }
  return local_finished;
}
/*
void sequentialComputation(int** board, int n, int t, int c, int max_iters){
  

  int red_1_count, blue_2_count;
  int tile_size = n/t;
  int cells_in_tile = n/t * n/t;
  int tile_row = n/(n/t);
  int tile_column = n/(n/t);

  bool *finished = false;
  int n_itrs = 0; 

  while (!finished && n_itrs < max_iters){
    n_itrs++;
    redMovement(board,n);
    blueMovement(board,n);
    int tid=0;
    determineConvergence(board,n,tile_size,tile_row,tile_column,cells_in_tile,c,finished,tid);
    if (finished) break;
  }
    printf("[Process #%d] %s [serial]: serial computation interation = [%d]\n", myid, time_local, n_itrs);
    for(i=0; i<M; i++){
      for(j=0; j<N; j++){
        printf("%d ", board[i][j]);
        if(j == N - 1)
        printf("\n");
      }
    }
  } 
}
*/

void *threadComputation(void *thread_arg){

  struct thread_data *t_data; 
  int i;
  int cell_start;
  int cell_end;
  int grid_start;
  int grid_end;

  bool local_finished;
  int n_itrs = 0; 

  /* Initialize thread portion of the sub-board */
  t_data = (struct thread_data *)thread_arg;
  int n = t_data->n; 
  int t = t_data->t;
  int c = t_data->c; 
  int max_iters = t_data->max_iters;
  int **board = t_data->board;
  int tid = t_data->tid;

  int tile_size = n/t;
  int cells_in_tile = n/t * n/t;
  int tile_row = n/(n/t);
  int tile_column = n/(n/t);

  while (!(local_finished) && n_itrs < max_iters){
    n_itrs++;
    redMovement(board,n);
    blueMovement(board,n);
    local_finished = determineConvergence(board,n,tile_size,tile_row,tile_column,cells_in_tile,c,tid);
    if (local_finished) break;
  }
  printf("\nTotal iterations = %d\n", n_itrs);
}

int main(int argc, char **argv) {

  /* user input preliminary check */
  if(argc < 6 || argc > 6 ) {
    printf("[Error]!!4 arguments expected!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
   }
  else if(atoi(argv[2])%atoi(argv[3])!=0){
    printf("[Error]!!Grid size n must be divisible by the tile size t!!\n");
    printf("SAMPLE: mpirun –np 5 main 25 5 50 500\n");
    exit(1);
  }
  else if(atoi(argv[2])<0||atoi(argv[3])<0||atoi(argv[4])<0||atoi(argv[5])<0||atoi(argv[4])>100){
    printf("[Error]!!Do NOT input negative numbers or threshold c larger than 100!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
  }

  int num_thrds = atoi(argv[1]);   //number of threads nthrds
  int n = atoi(argv[2]);           //cell grid size n
  int t = atoi(argv[3]);           //tile grid size t
  int c = atoi(argv[4]);           //terminating threshold c
  int max_iters = atoi(argv[5]);   //maximum number of iterations max_iters

  int i,j;

  pthread_t *thread_id;
  pthread_attr_t *attr; 
  struct thread_data *thread_data_array; 

  /* initialize global time variable  */
  now = time (0);
  strftime (time_local, 100, "%Y-%m-%d-%H:%M:%S", localtime (&now)); 
  
  /* create arrays of thread_ids and thread thread_arg_array */
  thread_id = (pthread_t *)malloc(sizeof(pthread_t) * num_thrds);
  thread_data_array = (struct thread_data *) malloc (sizeof(struct thread_data) * num_thrds);  

  ///////////////
  int n_itrs = 0;
  int myid,numprocs;
  int *board_cell, **board, **board_final;
  int K, K_sub;
  int q, r;
  //int ib, kn, i, j;
  bool finished;
  bool local_finished = false;
  bool global_finished = false;

  if (num_thrds==1) {

    printf("\n=================================================\n");
    printf("The grid size：n = %d\n", n);
    printf("The tile grid size: t = %d\n", t);
    printf("The terminating threshold: c = %d percent\n", c);
    printf("The maximum number of iterations: max_iters = %d\n", max_iters);
    printf("Current Local Time: %s\n", time_local);
    printf("=================================================\n");

    board = board_init(n,n);
    printf("\n[Process #%d] %s [serial]: INITIAL STATE of the board:\n\n", myid, time_local);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        printf("%d ", *(board[i] + j));
      }
      printf("\n");
    }
    printf("\n=================================================\n");

    printf("[Process #%d] %s [serial]: serial computation starting...\n", myid, time_local);
    clock_t begin = clock();
    
    /* distribute load and create threads for computation */
    thread_data_array[0].tid = 0;
    thread_data_array[0].board = board; 
    thread_data_array[0].n = n;
    thread_data_array[0].t = t;
    thread_data_array[0].max_iters = max_iters;
    thread_data_array[0].c = c;
    pthread_create(&thread_id[0], NULL, threadComputation, &thread_data_array[0]);
    
    /* Wait for all threads to complete then print global sum */ 
    for (i=0; i < 1; i++) {
      pthread_join(thread_id[i],NULL);
    }
    
    //sequentialComputation(board,n,t,c,max_iters);
    clock_t end = clock();
    double time_spent_serial = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("[Process #%d] %s [serial]: FINAL STATE of the board:\n", myid, time_local);
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
        printf("%d ", board[i][j]);
        if(j == n - 1)
        printf("\n");
      }
    }
    printf("[Process #%d] %s [serial]: computation convergence time [%f]\n", myid, time_local, time_spent_serial);
  }
  
  
  else {
    
  }
  
  
  pthread_exit (NULL);
}

int* calculateSubBoard(int board, int num_thrds, int n, int t, int tid ) {
      //compute a sub-board for every other processes
      //evenly allocate tiles so that each process has row difference no greater than one tile row
			
      int q, r; 
      int K, kn;
      int i, ib;

      q = t / num_thrds;       //minimum tile numbers each process handles
			r = t % num_thrds;       //tile remainders
      /* calculate assigned sub-board rows for process 0 */
      if (tid != r){
        K = (q+1) * n/t;
      }
      else{
        K = q * n/t;                                                          //K is the number of rows for this sub-board
      }
      kn = K * n;                                                             //kn is the number of cells in this sub-board

      printf("[Process #%d] %s [parallel]: partitioning - process #0 starting at row ib = 0, No. tile rows t_r = %d, No. board rows K = %d.\n", tid, time_local, K/(n/t),K);
      /* process 0 allocate memory for its own sub-board */
      /* duplicate the values for sub-board assigned to process 0 for computation, 2 ghost rows introduced */
      /*
      int** p_board_row_0 = board_init(K+2,N); //Row 0 and Row (k-1) are ghost rows
      for(i=1; i<K+1; i++){
        for(j=0; j<N; j++){
          p_board_row_0[i][j] = p_board_row[i-1][j];
        }
      }
      
      /* calculate and distribut sub-board rows for other processes */
			for (i=1; i<num_thrds; i++){                                             //since tile unit is considered, multiply row pointers by tize_size n/t
				//evenly allocate remainder tiles when r != 0
				if (i < r){
					ib = i * (q+1) * n/t;          //ib is the first row of each block assigned to other processes, therefore q multiply by tile size
					K = (q+1) * n/t;
				}
				else{ //evenly allocate remainder tiles when r = 0
					ib = i * q * n/t + r * n/t;
					K = q * n/t;                   //K is the number of rows in this sub-board
				}
        kn = K * n;                      //K * N is the size of the assigned sub-board

				printf("[Process #%d] %s [parallel]: partitioning - process #%d starting at row ib = %d, No. tile rows t_r = %d, No. board rows K = %d.\n", tid, time_local, i, ib, K/(n/t), K);

        printf("[Process #%d] %s [parallel]: partition distributed to process #%d...\n", tid, time_local, i);
			}
      printf("[Process #%d] %s [parallel]: computation in process...please wait...\n", tid, time_local);
      /*Re-calculate assigned sub-board rows for process 0 to fix a BUG*/
      if (tid != r){
        K = (q+1) * n/t;
      }
      else{
        K = q * n/t;
      }
      kn = K * n;
}