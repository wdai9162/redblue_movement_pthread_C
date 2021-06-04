#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>

/* define a function to initialize the matrix */
int* matrixInit(int N_vector, int M_length) {
  int* matrix; 
  int i,j; 
  /* set current time as seed for random generator */
  time_t t;
	srand((unsigned) time(&t));

  /* allocate memory for the matrix and assign a random integer value to each element */
  matrix = (int *)malloc(N_vector * M_length * sizeof(int)); 
  for (i=0; i<N_vector; i++){
    for (j=0; j<M_length; j++){
      *(matrix + i*M_length + j) = rand()%10;
    }
  }
  return matrix; 
}

void printMatrix(int* matrix,int N_vector, int M_length){
  int i,j;
  for (i=0; i<N_vector; i++){
    for (j=0; j<M_length; j++){
      printf("%d ", *(matrix + i*M_length + j));
    }
   printf("\n");
   }
};

void print_result(int *result,int *result_row_index,int N){
  int index; 
  for (int i=0; i<N; i++){
    for (int j=i; j<N; j++){
      index = result_row_index[i]+j-i; 
      printf("%d ", result[index]);
    }
   printf("\n");
   }
};

int main(int argc, char **argv) {
  clock_t begin = clock();
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);
  int myid,numprocs;
  int h,i,j,k;
  int *matrix, *recvbuf, *comp_buf, *coms_buf;
  int *result, *result_row_index, index, num_results, c;               //c is the temp result of each inner product calculation
  int *output, output_num, *output_recv_count, *displs;
  int *c_i, *c_j, *output_ci, *output_cj;

  char time_local[100];
  time_t now = time (0);
  strftime (time_local, 100, "%Y-%m-%d-%H:%M:%S", localtime (&now));

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  
  /* user input preliminary check */
  if(argc < 3 || argc > 3 || N%numprocs != 0 || N < numprocs || N%2!=0) {
    printf("[Error]!!2 arguments expected and N must be divisible by numprocs and N must be eqaul or larger than numprocs and N must be a even number!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_2_v1_0 <No. of vectors: int N> <length of vector: int M>\n");
    printf("SAMPLE: xxxxxxxxxxxx (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
   }
  /*
  else if(atoi(argv[1])%atoi(argv[2])!=0){
    printf("[Error]!!Grid size n must be divisible by the tile size t!!\n");
    printf("SAMPLE: mpirun –np 5 main 25 5 50 500\n");
    exit(1);
  }
  else if(atoi(argv[1])<0||atoi(argv[2])<0||atoi(argv[3])<0||atoi(argv[4])<0||atoi(argv[3])>100){
    printf("[Error]!!Do NOT input negative numbers or threshold c larger than 100!!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_1_red_blue_movement_v1_0 <grid size: int n> <tile size: int t> <threshold(%): int c> <max iteration: int max_iters>\n");
    printf("SAMPLE: mpirun –np 5 main 50 5 85 1000 (n MUST be divisible by t)(# of process MUST NOT be greater than t)\n");
    exit(1);
  }
  */

  /* ❗ sequential computation starts from here ❗ */
  if (numprocs == 1) {

    matrix = matrixInit(N,M);                               //initiate the matrix NXM
    result = (int*)malloc(N*(N+1)/2 * sizeof(int));         //allocate memory for results in Compressed Row Storage format (V[k])
    result_row_index = (int*)malloc(N * sizeof(int));       //allocate memory for index array for the above results (R[i])
    result_row_index[0] = 0;                                //assign the first index to be 0
    for (i=1; i<N; i++){
      result_row_index[i] = result_row_index[i-1]+N-i+1;    //calculate the index for the first elemnt of each row and store in index array
    }
    printf("\n[Process #%d] %s [sequential]: INITIAL STATE of the matrix:\n\n", myid, time_local);
    printMatrix(matrix,N,M);
    printf("[Process #%d] %s [sequential]: The row index are: ", myid, time_local);
    for (int i = 0; i < N; ++i) {
      printf("%d ", result_row_index[i]);
      }
    printf("\n");

    /* start pairwise inner product calculation */
    for (i=0; i<N; i++){                                    //iterate through each vector
      for (j=i; j<N; j++){                                  //iterate through each other vectors starting from itself
        c = 0;                                              //reset result 'c' for the next pair
        for (k=0; k<M; k++){                                //iterate through each element of both vectors
          c += matrix[i*M+k]*matrix[j*M+k];                 //calculate inner product "Cij"
        }
        index = result_row_index[i]+j-i;                    //calculate the index for Cij in Compressed Row Storage format
        result[index] = c;                                  //store calculated result for each pair in Compressed Row Storage format
      }
    }
    print_result(result,result_row_index,N);
    printf("\n");
  }

  /* ❗ parallel computation starts from here ❗ */
  else {

    //The problem of even number of vectors is that after load balancing the second/bottom half of the processes has 1 pair/step less load.
    //Same applies to partitioned blocks as well. First half of the processes will calculate one more block pair.  (numprocs/2+1) vs (numprocs/2) block pair
    //==================================================================================================================
    /*N=12                                                     Step0   1     2     3     4     5     6  (N/2+1) OR (N/2)
    00	01	02	03	04	05	06	07	08	09	010	011              00	  01	  02	  03	  04	  05	  06
    11	12	13	14	15	16	17	18	19	110	111	                 11	  12	  13	  14	  15	  16	  17
    22	23	24	25	26	27	28	29	210	211		                   22	  23	  24	  25	  26	  27  	28
    33	34	35	36	37	38	39	310	311			                     33	  34	  35	  36	  37	  38	  39
    44	45	46	47	48	49	410	411				                       44	  45	  46	  47	  48	  49	  410
    55	56	57	58	59	510	511					             =======>    55	  56	  57	  58	  59	  510	  511
    66	67	68	69	610	611						                           66	  67	  68	  69	  610	  611
    77	78	79	710	711							                             77	  78	  79	  710	  711	  07
    88	89	810	811								                               88	  89	  810	  811	  08	  18
    99	910	911									                                 99	  910	  911	  09	  19	  29
    1010	1011										                               1010	1011  010	  110	  210	  310
    1111											                                   1111	011	  111	  211	  311	  411
    */ 
   
    int src, dst;
    int block_N = N/numprocs;                                               //block_N is the number of vectors each process is assigned with, N of the sub matrix
    MPI_Request* 	reqs = (MPI_Request*)calloc(1, sizeof(MPI_Request));
    MPI_Status * 	stat=(MPI_Status*)malloc(sizeof(MPI_Status)*1);
    
    /* root process to initialize the matrix and ready to be distributed */
    if (myid == 0) {
      matrix = matrixInit(N,M);                                             //initiate the matrix for calculation
      printf("\n[Process #%d] %s [parallel]: INITIAL STATE of the matrix:\n\n", myid, time_local);
      printMatrix(matrix,N,M);
      printf("\n");
      
      output_num = N*(N+1)/2; 
      displs = (int *)malloc(numprocs * sizeof(int));
      output_recv_count = (int *)malloc(numprocs * sizeof(int));

      for (i=0; i<numprocs; i++){
        if (i==0) {
          displs[i] = 0; 
          output_recv_count[i] = (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2);
        }
        else if (i>0 && i<numprocs/2+1) {
          displs[i] = displs[i-1] + (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2);
          output_recv_count[i] = (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2);
        } else {
          displs[i] = displs[i-1] + (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2);
          output_recv_count[i] = (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2);
        }
      }
      printf("displs is: ");
      for (i=0; i<numprocs; i++){
        printf("%d ",displs[i]);
      }
      printf("\noutput_recv_count is:");
      for (i=0; i<numprocs; i++){
        printf("%d ",output_recv_count[i]);
      }
      printf("\n");
      output = (int *)malloc(output_num * sizeof(int));                     //allocate memory space for receiving the results
      output_ci = (int *)malloc(output_num * sizeof(int));                     //allocate memory space for receiving the results index Ci 
      output_cj = (int *)malloc(output_num * sizeof(int));                     //allocate memory space for receiving the results index Cj 
    }

    recvbuf = (int *)malloc(M * block_N * sizeof(int));                     //allocate memory space for the receive buff in each process
    MPI_Scatter (matrix, M * block_N, MPI_INT, recvbuf, M * block_N, MPI_INT, 0, MPI_COMM_WORLD);              
    //printf ("process: %d, received sub matrix:\n", myid);
    //printMatrix(recvbuf,N/numprocs,M);

    comp_buf = (int *)malloc(M * block_N * sizeof(int));                    //allocate memory space for the computation buff in each process
    memcpy( comp_buf, recvbuf, M * block_N *sizeof(int) );                  //copy the values over to computation buff 

    coms_buf = (int *)malloc(M * block_N * sizeof(int));                    //duplicate the computation buffer to generate the communication buffer for sending/receiving asynchronously
    memcpy( coms_buf, comp_buf, M * block_N *sizeof(int) );                 //copy the values over to communication buff

    if(myid<numprocs/2){
      //result = (int*)malloc((numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2) * sizeof(int));      //allocate memory for results in Compressed Row Storage format (V[k])
      //result_row_index = (int*)malloc(block_N * sizeof(int));                                       //allocate memory for index array for the above results (R[i])
      //result_row_index[0] = 0;                                                                      //assign the first index to be 0
      
      /* print row index TBD */
      //for (i=1; i<block_N; i++){
      //  result_row_index[i] = result_row_index[i-1]+block_N-i+1;                                    //calculate the index for the first elemnt of each row and store in index array
      //}

      /*
      printf("[Process #%d] %s [parallel]: The row index are: ", myid, time_local);
      for (int i = 0; i < block_N; i++ ) {
        printf("%d ", result_row_index[i]);
        }
      printf("\n");
      */
      //first half of the processes calculate (numprocs/2+1) block pairs and generate (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2) results
      num_results = (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2);
      result=(int*)malloc(num_results * sizeof(int));
      c_i=(int*)malloc(num_results * sizeof(int));
      c_j=(int*)malloc(num_results * sizeof(int));
      index=0;
      //calculation - block loop:
      for (h=0; h<(numprocs/2+1); h++){                                    //loop for (numprocs/2+1) steps for fisrt half of the processes
        if (h==0){                                                         //treat step 0 individually because it only calculate upper diagonal	
          //printf ("process: %d, loop value: %d\n", myid,h );
          //start pairwise inner product calculation 
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=i; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
                c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
              }
              //ci=i+myid*block_N
              //cj=j+myid*block_N+h*block_N
              c_i[index]= i+myid*block_N;
              c_j[index] = j+myid*block_N+h*block_N;
              //index = result_row_index[i]+j-i;                           //calculate the index for Cij in Compressed Row Storage format
              result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
              index++;
            }
          }
          if (myid==0){
            dst=numprocs-1;
          }else dst=myid-1;
          MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
        }
        else{
          if (myid==numprocs-1){
            src=0;
          }else src=myid+1;
          MPI_Irecv(coms_buf, M * block_N, MPI_INT, src, 0, MPI_COMM_WORLD, reqs);  //open a socket for receiving
          //printf ("process: %d, loop value: %d, receiving...\n", myid,h );
          MPI_Wait(reqs, stat);
          //printf("process: %d, MPI_Irecv received\n", myid);
          //printMatrix(coms_buf,N/numprocs,M);
          memcpy( comp_buf, coms_buf, M * block_N *sizeof(int) );          //copy newly received data from coms_buf to comp_buf

          if (h<numprocs/2 && myid !=0){
            MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data 
          } 
          else if (myid==0 && h<(numprocs/2-1)){
            MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data 
          }

          //start pairwise inner product calculation 
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=0; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
                c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
              }

              c_i[index]= i+myid*block_N;
              c_j[index] = j+myid*block_N+h*block_N;
              //index = result_row_index[i]+j-i;                           //calculate the index for Cij in Compressed Row Storage format
              result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
              index++;
            }
          }
        }
      }
      //printf ("process: %d, computation completed and my results are:", myid);
      //for (i=0;i<num_results;i++){
        //printf("%d ", result[i]);
      //} 
      
    }
    else if(numprocs/2 <= myid && myid <numprocs) {
    //second half of the processes calculate (numprocs/2) block pairs and generate (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2) results
    num_results = (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2);
    result = (int*)malloc( num_results * sizeof(int) );    //allocate memory for results in Compressed Row Storage format (V[k])
    c_i=(int*)malloc(num_results * sizeof(int));
    c_j=(int*)malloc(num_results * sizeof(int));
    index=0;
    //result_row_index = (int*)malloc(block_N * sizeof(int));                                             //allocate memory for index array for the above results (R[i])
    //result_row_index[0] = 0;                                                                            //assign the first index to be 0
    //for (i=1; i<block_N; i++){
     // result_row_index[i] = result_row_index[i-1]+block_N-i+1;                                          //calculate the index for the first elemnt of each row and store in index array
    //}
    /*
    printf("[Process #%d] %s [parallel]: The row index are: ", myid, time_local);
    for (int i = 0; i < block_N; i++ ) {
      printf("%d ", result_row_index[i]);
      }
    printf("\n");
    */
    //block loop:
    for (h=0; h<(numprocs/2); h++){                                      //loop for (numprocs/2) steps for second half of the processes
      if (h==0){
        //printf ("process: %d, loop value: %d\n", myid,h );
        //start pairwise inner product calculation 
        for (i=0; i < block_N; i ++){                                    //iterate through each vector
          for (j=i; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
            c = 0;                                                       //reset result 'c' for the next pair
            for (k=0; k<M; k++){                                         //iterate through each element of both vectors
              c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
            }
            c_i[index]= i+myid*block_N;
            c_j[index] = j+myid*block_N+h*block_N;
            result[index] = c;
            index++;                                           //store calculated result for each pair in Compressed Row Storage format
          }
        }
        if (myid==0){
          dst=numprocs-1;
        }else dst=myid-1;
        MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
      }
      else{
        if (myid==numprocs-1){
          src=0;
        }else src=myid+1;
        MPI_Irecv(coms_buf, M * block_N, MPI_INT, src, 0, MPI_COMM_WORLD, reqs);  //open a socket for receiving
        //printf ("process: %d, loop value: %d, receiving...\n", myid,h );
        MPI_Wait(reqs, stat);
        //printf("process: %d, MPI_Irecv received\n", myid);
        //printMatrix(coms_buf,N/numprocs,M);
        memcpy( comp_buf, coms_buf, M * block_N *sizeof(int) );
        if (h<numprocs/2 && myid == numprocs/2){
          MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
        } 
        else if (h<(numprocs/2-1) && myid!=numprocs/2){
          MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
        }
        
        if(myid>numprocs/2 && h > (numprocs-1-myid)){                //remap Cij number for balanced load
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=0; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
              c += comp_buf[i*M+k]*recvbuf[j*M+k];                       //calculate inner product "Cij"
              }
              c_i[index]= i+(h-(numprocs-myid))*block_N;
              c_j[index] = j+myid*block_N;
              result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
              index++;
            }
          }
        }
        else {
          //start pairwise inner product calculation 
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=0; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
                c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
              }
              c_i[index]= i+myid*block_N;
              c_j[index] = j+myid*block_N+h*block_N;
              //index = result_row_index[i]+j-i;                           //calculate the index for Cij in Compressed Row Storage format
              result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
              index++;
            }
          }
        }  
      }
    }
    }
  }
  

  printf("\n");


  MPI_Gatherv (result, num_results, MPI_INT, output, output_recv_count, displs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gatherv (c_i, num_results, MPI_INT, output_ci, output_recv_count, displs, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Gatherv (c_j, num_results, MPI_INT, output_cj, output_recv_count, displs, MPI_INT, 0, MPI_COMM_WORLD); 

    //root process to aggregate the results and run checks 
    if (myid==0) { 

    result_row_index = (int*)malloc(N * sizeof(int));       //allocate memory for index array for the above results (R[i])
    result_row_index[0] = 0;                                //assign the first index to be 0
    for (i=1; i<N; i++){
      result_row_index[i] = result_row_index[i-1]+N-i+1;    //calculate the index for the first elemnt of each row and store in index array
    }
    printf("[Process #%d] %s [sequential]: The row index are: ", myid, time_local);
    for (int i = 0; i < N; ++i) {
      printf("%d ", result_row_index[i]);
      }
    printf("\n");
    int* final_result=(int*)malloc(N*(N+1)/2 * sizeof(int));

    for (i=0; i<N*(N+1)/2; i++){
      index = result_row_index[output_ci[i]]+output_cj[i]-output_ci[i]; 
      final_result[index]=output[i];
    }
    print_result(final_result,result_row_index,N);

    printf("Sequential computation self check starting...");
    sleep(3);

    }

  MPI_Finalize();
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("time =: [%f]\n", time_spent);
  return 0;
}