#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

/*define a function to initialize the board as required */
int** board_init(int row, int column) {

  /* allocate memory block for all the cells on the board */
  int* p_board_cell = malloc(row * column * sizeof(int));
  /* allocate memory block for the 1st pointers of each row */
  int** p_board_row = malloc(row * sizeof(int*));
  if(p_board_row == NULL){
    fprintf(stderr, "**p_board_row out of memory\n");
    exit(1);
  }

  /* assign the row pointer to the right place on the board  */
  for(int i=0; i<row; i++) {
    p_board_row[i] = &p_board_cell[i*column];
  }

  /* set current time as seed for random generator */
  time_t t;
	srand((unsigned) time(&t));
  for (int i = 0; i < row; i++)
     for (int j = 0; j < column; j++)
        *(p_board_row[i] + j) = rand()%3; //fill the cells with random values 0,1,2
/*
  printf("Board initial state is displayed as below:\n");
  for (int i = 0; i < row; i++) {
     for (int j = 0; j < column; j++) {
        printf("%d ", *(p_board_row[i] + j));
     }
     printf("\n");
  }
  printf("===========\n");
*/
  return p_board_row;
}

int main(int argc, char **argv) {

  /* user input preliminary check */
  if(argc < 5 || argc > 5 ) {
    printf("[Error]!!4 arguments expected!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
   }
  else if(atoi(argv[1])%atoi(argv[2])!=0){
    printf("[Error]!!Grid size n must be divisible by the tile size t!!\n");
    printf("SAMPLE: mpirun –np 5 main 25 5 50 500\n");
    exit(1);
  }
  else if(atoi(argv[1])<0||atoi(argv[2])<0||atoi(argv[3])<0||atoi(argv[4])<0||atoi(argv[3])>100){
    printf("[Error]!!Do NOT input negative numbers or threshold c larger than 100!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
  }

  int n = atoi(argv[1]); /* grid size */
  int t = atoi(argv[2]); /* tile size */
  int c = atoi(argv[3]); /* terminating threshold c  */
  int max_iters = atoi(argv[4]); /* maximum number of iterations  */
  int n_itrs = 0;
  int myid,numprocs;
  int *p_board_cell, **p_board_row, **p_board_row_final;
  int M = atoi(argv[1]);
  int N = atoi(argv[1]);
  int K, K_sub;
  int q, r;
  int ib, kn, i, j;
  MPI_Status status;
  bool finished;
  bool local_finished = false;
  bool global_finished = false;

  char time_local[100];
  time_t now = time (0);
  strftime (time_local, 100, "%Y-%m-%d-%H:%M:%S", localtime (&now));

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);

  /* user input preliminary check */
  if(numprocs > t){
    printf("[Error]!!You are wasting computational resource!!Be frugal!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(percent): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(0);
  }

  /* serial execution */
  if (numprocs == 1) {
    printf("\n=================================================\n");
    printf("The grid size：n = %d\n", n);
    printf("The tile grid size: t = %d\n", t);
    printf("The terminating threshold: c = %d\n", c);
    printf("The maximum number of iterations: max_iters = %d\n", max_iters);
    printf("Current Local Time: %s\n", time_local);
    printf("=================================================\n");

    /*  process 0 initialize the board!  */
    clock_t begin = clock();
    p_board_row = board_init(M,N);
    printf("\n[Process #%d] %s [serial]: INITIAL STATE of the board:\n\n", myid, time_local);
    for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            printf("%d ", *(p_board_row[i] + j));
         }
         printf("\n");
      }
    printf("\n=================================================\n");

    printf("[Process #%d] %s [serial]: serial computation starting...\n", myid, time_local);

    int red_1_count, blue_2_count;
    int tile_size = M/t;
    int cells_in_tile = M/t * M/t;
    int tile_row = M/(M/t);
    int tile_column = N/(N/t);

    while (!finished && n_itrs < max_iters){

      n_itrs++;
      /***** Stage 1: Red Movement ******/
      for (int i = 0; i < M; i++){                                            //row loop
        if (p_board_row[i][0] == 1 && p_board_row[i][1] == 0){                //check edge case where the first cell is red and can move LEFT, so the last red does NOT move.
          p_board_row[i][0] = 4;
          p_board_row[i][1] = 3;
        }
        for (int j = 1; j < N; j++){                                          //column loop
          if (p_board_row[i][j] == 1 && p_board_row[i][(j+1)%n] == 0){        //when red can move right; (j+1)%n ensures if j = <last column>, the pointer moves back to first column;
          p_board_row[i][j] = 0;
          p_board_row[i][(j+1)%n] = 3;
          }
          else if (p_board_row[i][j] == 3) p_board_row[i][j] = 1;
        }
        if (p_board_row[i][0] == 3) p_board_row[i][0] = 1;                    //2nd time to check first cell in this row iteration
        else if (p_board_row[i][0] == 4) p_board_row[i][0] = 0;
      }
      /***** Stage 2: Blue Movement ******/
      for (int j = 0; j < n; j++){                                            //column loop
        if (p_board_row[0][j] == 2 && p_board_row[1][j] == 0){                //check edge case where the first cell is blue and can move DOWN, so the bottom blue does NOT move.
        p_board_row[0][j] = 4;
        p_board_row[1][j] = 3;
      }
      for (int i = 1; i < n; i++){                                            //row loop
        if (p_board_row[i][j] == 2 && p_board_row[(i+1)%n][j] == 0){          //when blue can move down wraparound; (i+1)%n ensures if i = <last row>, the pointer moves back to first row;
          p_board_row[i][j] = 0;
          p_board_row[(i+1)%n][j] = 3;
        }
        else if (p_board_row[i][j] == 3) p_board_row[i][j] = 2;
      }
      if (p_board_row[0][j] == 3) p_board_row[0][j] = 2;                  //2nd time to check first cell in this column iteration
      else if (p_board_row[0][j] == 4 ) p_board_row[0][j] = 0;
    }
      /***** Stage 3: Determine if the computation has converged ******/
      /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
      // Loop structure:
      //-tile_row(t_r) loop
      //   -tile_column(t_c) loop
      //         -cell_row(c_r) loop
      //              -cell_column(c_c) loop
      for (int t_r = 0; t_r < tile_row; t_r++) {
         for (int t_c = 0; t_c < tile_column; t_c++) {

           /* initialize and reset red/blue count to Zero for each tile*/
           red_1_count = 0;
           blue_2_count = 0;

           for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {    //(t_r+1)*t is the start row of the below tile
             for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){   //(t_c+1)*t is the start column of the tile on the right

               //printf("%d ", *(p_board_row[c_r]+c_c));          //Print Tiles
               if (p_board_row[c_r][c_c] == 1) {
                 red_1_count+=1;
               }         //if cell value is 1, red count plus 1;
               else if (p_board_row[c_r][c_c] == 2){
                 blue_2_count+=1;
               }        //if cell value is 2, blue count plus 1;
             }
             //printf("\n");                                  //Print Tiles - next row
           }

           /* calculate color percentage and compare with threshold to decide whether to terminate */
           double red_percentage = ((double)red_1_count/cells_in_tile)*100;
           double blue_percentage = ((double)blue_2_count/cells_in_tile)*100;

           if (red_percentage > c) {
             finished = true;
             printf("[Process #%d] %s [serial]: converged RED with %d red cells ==> [%.2f] ==> converged tile:\n", myid, time_local, red_1_count, red_percentage);
             /* print out the converged tile */
             for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
                for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                  printf("%d ", p_board_row[c_r][c_c]);
                }
                printf("\n");
              }
             break;
           }
           else if (blue_percentage > c) {
             finished = true;
             printf("[Process #%d] %s [serial]: converged BLUE with %d blue cells ==> [%.2f] ==> converged tile:\n",myid, time_local, blue_2_count, blue_percentage);
             /* print out the converged tile */
             for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
               for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                 printf("%d ", p_board_row[c_r][c_c]);
               }
               printf("\n");
             }
             break;
           }
           //printf("red = %d ==> %.2f\n",red_1_count_check, red_percentage_check);  //Print Tiles - result
           //printf("blu = %d ==> %.2f\n",blue_2_count_check, blue_percentage_check);//Print Tiles - result
           //printf("==========\n");
         }
         if (finished) break;
      }
      printf("[Process #%d] %s [serial]: serial computation interation = [%d]\n", myid, time_local, n_itrs);
      for(i=0; i<M; i++){
        for(j=0; j<N; j++){
          printf("%d ", p_board_row[i][j]);
          if(j == N - 1)
          printf("\n");
        }
      }
    }
      clock_t end = clock();
      double time_spent_serial = (double)(end - begin) / CLOCKS_PER_SEC;

      printf("[Process #%d] %s [serial]: FINAL STATE of the board:\n", myid, time_local);
      for(i=0; i<M; i++){
        for(j=0; j<N; j++){
          printf("%d ", p_board_row[i][j]);
          if(j == N - 1)
          printf("\n");
        }
      }
      printf("[Process #%d] %s [serial]: computation convergence time [%f]\n", myid, time_local, time_spent_serial);

    MPI_Finalize();
  }

  /* parallel execution */
  else {
    /* Process 0 actions */
    if (myid == 0) {

      printf("\n=================================================\n");
      printf("The grid size：n = %d\n", n);
      printf("The tile grid size: t = %d\n", t);
      printf("The terminating threshold: c = %d\n", c);
      printf("The maximum number of iterations: max_iters = %d\n", max_iters);
      printf("Current Local Time: %s\n", time_local);
      printf("=================================================\n");

      clock_t begin = clock();
      /*  process 0 initialize the board!  */
  		p_board_row = board_init(M,N);
      /*  process 0 allocate memory to hold the board FINAL STATE  */
      p_board_row_final = board_init(M,N);
  		if(p_board_row == NULL || p_board_row_final == NULL){
  			fprintf(stderr, "out of memory\n");
  			exit(1);
  		}

      printf("\n[Process #%d] %s [parallel]: INITIAL STATE of the board:\n\n", myid, time_local);
      /*
      for (int i = 0; i < n; i++) {
           for (int j = 0; j < n; j++) {
              printf("%d ", *(p_board_row[i] + j));
           }
           printf("\n");
        }
      printf("\n=================================================\n");
      */

			//compute a sub-board for every other processes
      //evenly allocate tiles so that each process has row difference no greater than one tile row
			q = t / numprocs;       //minimum tile numbers each process handles
			r = t % numprocs;       //tile remainders
      /* calculate assigned sub-board rows for process 0 */
      if (myid != r){
        K = (q+1) * n/t;
      }
      else{
        K = q * n/t;                                                          //K is the number of rows for this sub-board
      }
      kn = K * N;
      printf("[Process #%d] %s [parallel]: partitioning - process #0 starting at row ib = 0, No. tile rows t_r = %d, No. board rows K = %d.\n", myid, time_local, K/(M/t),K);
      /* process 0 allocate memory for its own sub-board */
      /* duplicate the values for sub-board assigned to process 0 for computation, 2 ghost rows introduced */
      int** p_board_row_0 = board_init(K+2,N); //Row 0 and Row (k-1) are ghost rows
      for(i=1; i<K+1; i++){
        for(j=0; j<N; j++){
          p_board_row_0[i][j] = p_board_row[i-1][j];
        }
      }
      /* calculate and distribut sub-board rows for other processes */
			for (i=1; i<numprocs; i++){                                             //since tile unit is considered, multiply row pointers by tize_size n/t
				//evenly allocate remainder tiles when r != 0
				if (i < r){
					ib = i * (q+1) * n/t;          //ib is the first row of each block assigned to other processes, therefore q multiply by tile size
					K = (q+1) * n/t;
				}
				else{ //evenly allocate remainder tiles when r = 0
					ib = i * q * n/t + r * n/t;
					K = q * n/t;                   //K is the number of rows in this sub-board
				}
        kn = K * N;                      //K * N is the size of the assigned sub-board

				printf("[Process #%d] %s [parallel]: partitioning - process #%d starting at row ib = %d, No. tile rows t_r = %d, No. board rows K = %d.\n", myid, time_local, i, ib, K/(M/t), K);
				MPI_Send(&p_board_row[ib][0], kn, MPI_INT, i, 1, MPI_COMM_WORLD); //send row pointer to the other processes
        printf("[Process #%d] %s [parallel]: partition distributed to process #%d...\n", myid, time_local, i);
			}
      printf("[Process #%d] %s [parallel]: computation in process...please wait...\n", myid, time_local);
      /*Re-calculate assigned sub-board rows for process 0 to fix a BUG*/
      if (myid != r){
        K = (q+1) * n/t;
      }
      else{
        K = q * n/t;
      }
      kn = K * N;

      while (!local_finished && n_itrs < max_iters){
        n_itrs++;
      /* Process 0 computation starts from here */
      MPI_Sendrecv(&p_board_row_0[K][0], n, MPI_INT, (myid+1)%numprocs, 2, &p_board_row_0[0][0], n, MPI_INT, (myid-1+numprocs)%numprocs, 2, MPI_COMM_WORLD, &status);
      /* send first data row pointer &p_board_row[1][0] to process myid-1 and receive bottom ghost row pointer &p_board_row[K+1][0] from process myid+1 and update the local sub-board */
      MPI_Sendrecv(&p_board_row_0[1][0], n, MPI_INT, (myid-1+numprocs)%numprocs, 3, &p_board_row_0[K+1][0], n, MPI_INT, (myid+1)%numprocs, 3, MPI_COMM_WORLD, &status);

      /***** Stage 1: Red Movement ******/
        /*row loop*/
        for (int i = 0; i < K+2; i++){                                          //include two ghost rows
          if (p_board_row_0[i][0] == 1 && p_board_row_0[i][1] == 0){            //evaluate edge case where the 1st cell is red and can move LEFT ==> the last red does NOT move.
            p_board_row_0[i][0] = 4;
            p_board_row_0[i][1] = 3;
          }
          /*column loop*/
          for (int j = 1; j < N; j++){
            if ( p_board_row_0[i][j] == 1 && p_board_row_0[i][(j+1)%N] == 0){       //when red can move RIGHT; (j+1)%N ensures if j = N-1, the pointer moves back to the first column;
            p_board_row_0[i][j] = 0;
            p_board_row_0[i][(j+1)%N] = 3;
            }
            else if ( p_board_row_0[i][j] == 3) p_board_row_0[i][j] = 1;
          }
          if (p_board_row_0[i][0] == 3) p_board_row_0[i][0] = 1;                    //2nd time to evaluate the 1st cell in this row iteration
          else if (p_board_row_0[i][0] == 4) p_board_row_0[i][0] = 0;               //flip place holder 3/4 back to 1 or 0. conclude this row iteration.
        }
      /***** Stage 2: Blue Movement ******/
        /*column loop*/
        for (int j = 0; j < N; j++){                                                //include ghost row 0 and K+1 for blue movement.
          if ( p_board_row_0[0][j] == 2 && p_board_row_0[1][j] == 0){               //evaluate edge case where the 1st cell is blue and can move DOWN ==> the bottom blue does NOT move.
          p_board_row_0[0][j] = 4;
          p_board_row_0[1][j] = 3;
        }
        /*row loop*/
        for (int i = 1; i < K+2; i++){
          if ( p_board_row_0[i][j] == 2 && p_board_row_0[(i+1)%(K+2)][j] == 0){         //when blue can move DOWN; (i+1)%(K+2) ensures if i = K+1, the pointer moves back to the first row;
            p_board_row_0[i][j] = 0;
            p_board_row_0[(i+1)%(K+2)][j] = 3;
          }
          else if ( p_board_row_0[i][j] == 3) p_board_row_0[i][j] = 2;
        }
        if ( p_board_row_0[0][j] == 3) p_board_row_0[0][j] = 2;                     //2nd time to evaluate the 1st cell in this column iteration
        else if ( p_board_row_0[0][j] == 4 ) p_board_row_0[0][j] = 0;               //flip place holder 3/4 back to 2 or 0. conclude this column iteration.
      }
      /***** Stage 3: Check Results ******/
      /* count the number of red and blue in each tile and check if the computation can be terminated*/
      /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
      // Loop structure:
      //-tile_row(t_r) loop
      //   -tile_column(t_c) loop
      //         -cell_row(c_r) loop
      //              -cell_column(c_c) loop
      int tile_size = M/t;                                                      // this is a n/t * n/t tile, i.e board is 6x6, t = 2, therefore tile is 6/2 x 6/2 = 3x3;
      int cells_in_tile = M/t * M/t;
      int red_1_count, blue_2_count;
      int tile_row = K/(M/t);
      int tile_column = N/(N/t);
      for (int t_r = 0; t_r < tile_row; t_r++) {
       for (int t_c = 0; t_c < tile_column; t_c++) {
         /* initialize and reset red/blue count to Zero for each tile*/
         red_1_count = 0;
         blue_2_count = 0;

         for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {//(t_r+1)*t is the start row of the below tile, +1 for the ghost row
           for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){     //(t_c+1)*t is the start column of the tile on the right
             //printf("t_r=%d, t_c=%d, c_r=%d, c_c=%d\n", t_r, t_c, c_r, c_c);
             //printf("%d ", p_board_row_0[c_r][c_c]);                              //Print Tiles
             if (p_board_row_0[c_r][c_c] == 1) red_1_count+=1;                    //if cell value is 1, red count plus 1;
             else if (p_board_row_0[c_r][c_c] == 2) blue_2_count+=1;              //if cell value is 2, blue count plus 1;
           }
           //printf("\n");                                                        //Print Tiles - next row
         }

         /* calculate color percentage and compare with threshold to decide whether to terminate */
         double red_percentage = ((double)red_1_count/cells_in_tile)*100;
         double blue_percentage = ((double)blue_2_count/cells_in_tile)*100;
         if (red_percentage > c) {
           local_finished = true;
           printf("[Process #%d] %s [parallel]: converged RED with %d red cells ==> [%.2f] ==> converged tile:\n", myid, time_local, red_1_count, red_percentage);
           /* print out the converged tile */
           for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {
              for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                printf("%d ", p_board_row_0[c_r][c_c]);
              }
              printf("\n");
            }
          break;
         }
         else if (blue_percentage > c) {
           local_finished = true;
           printf("[Process #%d] %s [parallel]: converged BLUE with %d blue cells ==> [%.2f] ==> converged tile:\n", myid, time_local, blue_2_count, blue_percentage);
           /* print out the converged tile */
           for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {
             for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
               printf("%d ", p_board_row_0[c_r][c_c]);
             }
             printf("\n");
           }
           break;
         }
         //printf("red = %d ==> %.2f\n",red_1_count, red_percentage);  //Print Tiles - result
         //printf("blu = %d ==> %.2f\n",blue_2_count, blue_percentage);//Print Tiles - result
         //printf("==========\n");
       }
       if (local_finished) break;
      }
      MPI_Allreduce(&local_finished, &global_finished, 1, MPI_C_BOOL, MPI_LOR,MPI_COMM_WORLD);
      if (global_finished&&local_finished){
        printf("[Process #%d] %s [parallel]: detected local_finished is TRUE, terminating...\n", myid, time_local);
        break;
        }
      else if (global_finished){
        printf("[Process #%d] %s [parallel]: detected global_finished is TRUE, terminating...\n", myid, time_local);
        break;
        }
      }

      //send sub_board to process 0 itself after movement computation
      MPI_Sendrecv(&p_board_row_0[1][0], kn, MPI_INT, 0, 4, &p_board_row_final[0][0], kn, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);

      //receive the calculated sub-board from every other processes
			for (i=1; i<numprocs; i++)	{
				//calculate the first row for each sub-board
				if (i < r){
					ib = i * (q+1) * n/t;          //ib is the first row of each block assigned to other processes, therefore q multiply by tile size
					K = (q+1) * n/t;
				}
				else{
					ib = i * q * n/t + r * n/t;
					K = q * n/t;
				}
				kn = K * N;
				MPI_Recv(&p_board_row_final[ib][0], kn, MPI_INT, i, 5, MPI_COMM_WORLD, &status);
			}
      clock_t end = clock();
      double time_spent_parallel = (double)(end - begin) / CLOCKS_PER_SEC;

      printf("=================================================\n");
      printf("\n[Process #%d] %s [parallel]: results received from all processes!\n", myid, time_local);
  		printf("[Process #%d] %s [parallel]: total iterations elapsed = [%d].\n", myid, time_local, n_itrs);
      printf("[Process #%d] %s [parallel]: FINAL STATE of the board:\n", myid, time_local);
  		/*
      for(i=0; i<M; i++){
  			for(j=0; j<N; j++){
  				printf("%d ", p_board_row_final[i][j]);
  				if(j == N - 1)
  				printf("\n");
  			}
  		}
      */

      /* Self checking serialization computation */
      printf("[Process #%d] %s [parallel][check]: self check starting in 5 seconds...\n", myid, time_local);
      printf("[Process #%d] %s [parallel][check]: 5...\n", myid, time_local);
      sleep(1);
      printf("[Process #%d] %s [parallel][check]: 4...\n", myid, time_local);
      sleep(1);
      printf("[Process #%d] %s [parallel][check]: 3...\n", myid, time_local);
      sleep(1);
      printf("[Process #%d] %s [parallel][check]: 2...\n", myid, time_local);
      sleep(1);
      printf("[Process #%d] %s [parallel][check]: 1...\n", myid, time_local);
      sleep(1);


      bool finished_check = false;
      int n_itrs_check = 0;
      int red_1_count_check, blue_2_count_check;
      int tile_size = M/t;
      int cells_in_tile = M/t * M/t;
      int tile_row = M/(M/t);
      int tile_column = N/(N/t);
      printf("[Process #%d] %s [parallel][check]: self checking computation in progress...\n", myid, time_local);

      clock_t begin_c = clock();
      while (!finished_check && n_itrs_check < max_iters){

        n_itrs_check++;
        /***** Stage 1: Red Movement ******/
        for (int i = 0; i < M; i++){                                            //row loop
          if (p_board_row[i][0] == 1 && p_board_row[i][1] == 0){                //check edge case where the first cell is red and can move LEFT, so the last red does NOT move.
            p_board_row[i][0] = 4;
            p_board_row[i][1] = 3;
          }
          for (int j = 1; j < N; j++){                                          //column loop
            if (p_board_row[i][j] == 1 && p_board_row[i][(j+1)%n] == 0){        //when red can move right; (j+1)%n ensures if j = <last column>, the pointer moves back to first column;
            p_board_row[i][j] = 0;
            p_board_row[i][(j+1)%n] = 3;
            }
            else if (p_board_row[i][j] == 3) p_board_row[i][j] = 1;
          }
          if (p_board_row[i][0] == 3) p_board_row[i][0] = 1;                    //2nd time to check first cell in this row iteration
          else if (p_board_row[i][0] == 4) p_board_row[i][0] = 0;
        }
        /***** Stage 2: Blue Movement ******/
        for (int j = 0; j < n; j++){                                            //column loop
          if (p_board_row[0][j] == 2 && p_board_row[1][j] == 0){                //check edge case where the first cell is blue and can move DOWN, so the bottom blue does NOT move.
          p_board_row[0][j] = 4;
          p_board_row[1][j] = 3;
        }
        for (int i = 1; i < n; i++){                                            //row loop
          if (p_board_row[i][j] == 2 && p_board_row[(i+1)%n][j] == 0){          //when blue can move down wraparound; (i+1)%n ensures if i = <last row>, the pointer moves back to first row;
            p_board_row[i][j] = 0;
            p_board_row[(i+1)%n][j] = 3;
          }
          else if (p_board_row[i][j] == 3) p_board_row[i][j] = 2;
        }
        if (p_board_row[0][j] == 3) p_board_row[0][j] = 2;                  //2nd time to check first cell in this column iteration
        else if (p_board_row[0][j] == 4 ) p_board_row[0][j] = 0;
      }
        /***** Stage 3: Determine if the computation has converged ******/
        /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
        // Loop structure:
        //-tile_row(t_r) loop
        //   -tile_column(t_c) loop
        //         -cell_row(c_r) loop
        //              -cell_column(c_c) loop
        for (int t_r = 0; t_r < tile_row; t_r++) {
           for (int t_c = 0; t_c < tile_column; t_c++) {

             /* initialize and reset red/blue count to Zero for each tile*/
             red_1_count_check = 0;
             blue_2_count_check = 0;

             for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {    //(t_r+1)*t is the start row of the below tile
               for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){   //(t_c+1)*t is the start column of the tile on the right

                 //printf("%d ", *(p_board_row[c_r]+c_c));          //Print Tiles
                 if (p_board_row[c_r][c_c] == 1) {
                   red_1_count_check+=1;
                 }         //if cell value is 1, red count plus 1;
                 else if (p_board_row[c_r][c_c] == 2){
                   blue_2_count_check+=1;
                 }        //if cell value is 2, blue count plus 1;
               }
               //printf("\n");                                  //Print Tiles - next row
             }

             /* calculate color percentage and compare with threshold to decide whether to terminate */
             double red_percentage_check = ((double)red_1_count_check/cells_in_tile)*100;
             double blue_percentage_check = ((double)blue_2_count_check/cells_in_tile)*100;

             if (red_percentage_check > c) {
               finished_check = true;
               printf("[Process #%d] %s [parallel][check]: converged RED with %d red cells ==> [%.2f] ==> converged tile:\n", myid, time_local, red_1_count_check, red_percentage_check);
               /* print out the converged tile */
               for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
                  for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                    printf("%d ", p_board_row[c_r][c_c]);
                  }
                  printf("\n");
                }
               break;
             }
             else if (blue_percentage_check > c) {
               finished_check = true;
               printf("[Process #%d] %s [parallel][check]: converged BLUE with %d blue cells ==> [%.2f] ==> converged tile:\n",myid, time_local, blue_2_count_check, blue_percentage_check);
               /* print out the converged tile */
               for (int c_r = t_r * tile_size; c_r < (t_r+1)*tile_size; c_r++) {
                 for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                   printf("%d ", p_board_row[c_r][c_c]);
                 }
                 printf("\n");
               }
               break;
             }
             //printf("red = %d ==> %.2f\n",red_1_count_check, red_percentage_check);  //Print Tiles - result
             //printf("blu = %d ==> %.2f\n",blue_2_count_check, blue_percentage_check);//Print Tiles - result
             //printf("==========\n");
           }
           if (finished_check) break;
        }
      }
      clock_t end_c = clock();
      double time_spent_check = (double)(end_c - begin_c) / CLOCKS_PER_SEC;

      printf("[Process #%d] %s [parallel][check]: self checking iteration = [%d] ==> FINAL CHECK STATE of the board:\n", myid, time_local,n_itrs_check);
      /*
      for(i=0; i<M; i++){
        for(j=0; j<N; j++){
          printf("%d ", p_board_row[i][j]);
          if(j == N - 1)
          printf("\n");
        }
      }
      */
      //validate each cell and count the difference
      int cell_dif = 0;
      for(i=0; i<M; i++){
        for(j=0; j<N; j++){
          if(p_board_row[i][j]!=p_board_row_final[i][j]) {
            cell_dif++;
          }
        }
      }
      if (cell_dif == 0){
        printf("[Process #%d] %s [parallel][check]: parallel computation matches self checking serial computation result!!\n", myid, time_local);
      }
      else {
        
      }
      printf("[Process #%d] %s [parallel]: parallel computation convergence time [%f]\n", myid, time_local, time_spent_parallel);
      printf("[Process #%d] %s [parallel]:  serial  computation convergence time [%f]\n", myid, time_local, time_spent_check);
    }

    /* all other processes actions */
    else {
  	//create a sub-board of size (K+2) X N.
    q = t / numprocs;       //evenly allocate number of tiles
    r = t % numprocs;       //tile remainders
    if (myid < r){
      ib = i * (q+1) * n/t;
      K = (q+1) * n/t;
    }
    else{
      ib = i * q * n/t + r * n/t;     //ib is the first row of each block assigned to other processes, multiply by tile size
      K = q * n/t;                    //K is the number of rows in this block
    }
		kn = K * N;
    //printf("myid = %d, K = %d.\n", myid, K);
    p_board_row = board_init(K+2,N); //Row 0 and Row (k-1) are ghost rows
  	/* recv the allocated sub-board from Process 0.*/
  	MPI_Recv(&p_board_row[1][0], kn, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);   //leaving &p_board_row[0][0] as pointer to ghost row

    printf("[Process #%d] %s [parallel]: computation in process...please wait...\n", myid, time_local);

    while (!local_finished && n_itrs < max_iters){
      n_itrs++;
      //printf("n_itrs = %d\n", n_itrs);

      /* send last data row pointer &p_board_row[K][0] to process myid+1 and receive top ghost row pointer &p_board_row[0][0] from process myid-1 and update the local sub-board */
      MPI_Sendrecv(&p_board_row[K][0], n, MPI_INT, (myid+1)%numprocs, 2, &p_board_row[0][0], n, MPI_INT, (myid-1+numprocs)%numprocs, 2, MPI_COMM_WORLD, &status);
      /* send first data row pointer &p_board_row[1][0] to process myid-1 and receive bottom ghost row pointer &p_board_row[K+1][0] from process myid+1 and update the local sub-board */
      MPI_Sendrecv(&p_board_row[1][0], n, MPI_INT, (myid-1+numprocs)%numprocs, 3, &p_board_row[K+1][0], n, MPI_INT, (myid+1)%numprocs, 3, MPI_COMM_WORLD, &status);

      /***** Stage 1: Red Movement ******/
        /*row loop*/
        for (int i = 0; i < K+2; i++){                                          //include ghost rows for computation
          if (p_board_row[i][0] == 1 && p_board_row[i][1] == 0){                //evaluate edge case where the 1st cell is red and can move LEFT ==> the last red does NOT move.
            p_board_row[i][0] = 4;
            p_board_row[i][1] = 3;
          }
          /*column loop*/
          for (int j = 1; j < N; j++){
            if ( p_board_row[i][j] == 1 && p_board_row[i][(j+1)%N] == 0){       //when red can move RIGHT; (j+1)%N ensures if j = N-1, the pointer moves back to the first column;
            p_board_row[i][j] = 0;
            p_board_row[i][(j+1)%N] = 3;
            }
            else if ( p_board_row[i][j] == 3) p_board_row[i][j] = 1;
          }
          if (p_board_row[i][0] == 3) p_board_row[i][0] = 1;                    //2nd time to evaluate the 1st cell in this row iteration
          else if (p_board_row[i][0] == 4) p_board_row[i][0] = 0;               //flip place holder 3/4 back to 1 or 0. conclude this row iteration.
        }
      /***** Stage 2: Blue Movement ******/
        /*column loop*/
        for (int j = 0; j < N; j++){                                            //include ghost row 0 and K+1 for blue movement.
          if ( p_board_row[0][j] == 2 && p_board_row[1][j] == 0){               //evaluate edge case where the 1st cell is blue and can move DOWN ==> the bottom blue does NOT move.
          p_board_row[0][j] = 4;
          p_board_row[1][j] = 3;
        }
        /*row loop*/
        for (int i = 1; i < K+2; i++){
          if ( p_board_row[i][j] == 2 && p_board_row[(i+1)%(K+2)][j] == 0){         //when blue can move DOWN; (i+1)%(K+2) ensures if i = K+1, the pointer moves back to the first row;
            p_board_row[i][j] = 0;
            p_board_row[(i+1)%(K+2)][j] = 3;
          }
          else if ( p_board_row[i][j] == 3) p_board_row[i][j] = 2;
        }
        if ( p_board_row[0][j] == 3) p_board_row[0][j] = 2;                     //2nd time to evaluate the 1st cell in this column iteration
        else if ( p_board_row[0][j] == 4 ) p_board_row[0][j] = 0;               //flip place holder 3/4 back to 2 or 0. conclude this column iteration.
      }
      /***** Stage 3: Check Results ******/
      /* count the number of red and blue in each tile and check if the computation can be terminated*/
      /* check every tiles to see if any tile’s colored cells are more than c% in one color (blue or red) */
      // Loop structure:
      //-tile_row(t_r) loop
      //   -tile_column(t_c) loop
      //         -cell_row(c_r) loop
      //              -cell_column(c_c) loop
      int tile_size = M/t;                                                      // this is a n/t * n/t tile, i.e board is 6x6, t = 2, therefore tile is 6/2 x 6/2 = 3x3;
      int cells_in_tile = M/t * M/t;
      int red_1_count, blue_2_count;
      int tile_row = K/(M/t);
      int tile_column = N/(N/t);
      for (int t_r = 0; t_r < tile_row; t_r++) {
       for (int t_c = 0; t_c < tile_column; t_c++) {
         /* initialize and reset red/blue count to Zero for each tile*/
         red_1_count = 0;
         blue_2_count = 0;

         for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {//(t_r+1)*t is the start row of the below tile, +1 for the ghost row
           for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){     //(t_c+1)*t is the start column of the tile on the right
             //printf("t_r=%d, t_c=%d, c_r=%d, c_c=%d\n", t_r, t_c, c_r, c_c);
             //printf("%d ", p_board_row[c_r][c_c]);                              //Print Tiles
             if (p_board_row[c_r][c_c] == 1) red_1_count+=1;                    //if cell value is 1, red count plus 1;
             else if (p_board_row[c_r][c_c] == 2) blue_2_count+=1;              //if cell value is 2, blue count plus 1;
           }
           //printf("\n");                                                        //Print Tiles - next row
         }

         /* calculate color percentage and compare with threshold to decide whether to terminate */
         double red_percentage = ((double)red_1_count/cells_in_tile)*100;
         double blue_percentage = ((double)blue_2_count/cells_in_tile)*100;
         if (red_percentage > c) {
           local_finished = true;
           printf("[Process #%d] %s [parallel]: converged RED with %d red cells ==> [%.2f] ==> converged tile:\n", myid, time_local, red_1_count, red_percentage);
           /* print out the converged tile */
           for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {
              for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
                printf("%d ", p_board_row[c_r][c_c]);
              }
              printf("\n");
            }
          break;
         }
         else if (blue_percentage > c) {
           local_finished = true;
           printf("[Process #%d] %s [parallel]: converged BLUE with %d blue cells ==> [%.2f] ==> converged tile:\n", myid, time_local, blue_2_count, blue_percentage);
           /* print out the converged tile */
           for (int c_r = t_r * tile_size + 1; c_r < (t_r+1)*tile_size+1; c_r++) {
             for (int c_c = t_c * tile_size; c_c < (t_c+1)*tile_size; c_c++){
               printf("%d ", p_board_row[c_r][c_c]);
             }
             printf("\n");
           }
           break;
         }
          /* print out result for every tile */
         //printf("red = %d ==> %.2f\n",red_1_count, red_percentage);
         //printf("blu = %d ==> %.2f\n",blue_2_count, blue_percentage);
         //printf("===============\n");
       }
       if (local_finished) break;
      }

    MPI_Allreduce(&local_finished, &global_finished, 1, MPI_C_BOOL, MPI_LOR,MPI_COMM_WORLD);
    if (global_finished&&local_finished){
      printf("[Process #%d] %s [parallel]: detected local_finished is TRUE, terminating...\n", myid, time_local);
      break;
      }
    else if (global_finished){
      printf("[Process #%d] %s [parallel]: detected global_finished is TRUE, terminating...\n", myid, time_local);
      break;
      }
    }
    MPI_Send(&p_board_row[1][0], kn, MPI_INT, 0, 5, MPI_COMM_WORLD);
    }

    MPI_Finalize();
  }
  return 0;
}
