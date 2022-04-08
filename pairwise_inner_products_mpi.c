#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>
struct timespec start, finish;
/* initialize all the functions */
int* matrixInit(int N_vector, int M_length);
void printMatrix(int* matrix,int N_vector, int M_length);
void print_result(int *result,int *result_row_index,int N);
void calculateInnerProduct(int *matrix, int *result, int *result_row_index, int N_vector, int M_length);
void calculateInnerProductParallel(int *recvbuf, int *comp_buf, int *result, int *c_i, int *c_j, int block_N, int M, int index, int myid, int h, int c);
void calculateDisplacementRecvCount(int *displs, int *output_recv_count, int numprocs, int block_N);
int checkResults(int *final_result, int *result_check, int output_num, int dif); 

int main(int argc, char **argv) {
  
  char time_local[100];
  time_t now = time (0);
  strftime (time_local, 100, "%Y-%m-%d-%H:%M:%S", localtime (&now));
  int N = atoi(argv[1]);                                  // N vectors = N rows in the matrix
  int M = atoi(argv[2]);                                  // M is the length of each vector = M columns in the matrix
  int myid,numprocs;
  int h,i,j,k;
  
  int *matrix;                                            // matrix is randomly generated for computation 
  int *recvbuf;                                           // recvbuf holds the static sub-matrix assigned to each process
  int *comp_buf;                                          // comp_buf holds the dynamic sub-matrix for inner product computation
  int *coms_buf;                                          // coms_buf is used to exchange the dynamic sub-matrix between processes

  int *result, *result_row_index, index, num_results, c;
  int *output, output_num, *output_recv_count, *displs;
  int *c_i, *c_j, *output_ci, *output_cj;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  
  /* user input preliminary check */
  if(argc < 3 || argc > 3 || N%numprocs != 0 || N < numprocs || numprocs%2!=0) {
    printf("[Error]!!2 arguments expected; number of processes must be an even number; N must be divisible by numprocs and N must be eqaul or larger than numprocs!\n");
    printf("Input Format: mpirun –np <# of process> comp5426_wdai9162_assignment_2_v1_0 <No. of vectors: int N> <length of vector: int M>\n");
    printf("SAMPLE: mpirun -np 2 assignment_2 30 10\n");
    exit(1);
   }

  /* ❗ sequential computation starts from here */
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

    /* start pairwise inner product calculation */
    calculateInnerProduct(matrix, result, result_row_index, N, M);
    print_result(result,result_row_index,N);
  }

  /* ❗ parallel computation starts from here */
  else {

    //The problem of even number of vectors is that after load balancing the second/bottom half of the processes has 1 pair/step less load.
    //Same applies to partitioned blocks as well. First half of the processes will calculate one more block pair => (numprocs/2+1) vs (numprocs/2) block pair
    //==================================================================================================================
    /*For example N=12                                         Step0   1     2     3     4     5     6  (N/2+1) =7  VS (N/2) = 6
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
   
    int src, dst;                                                           //source and destiantion processes for MPI communication
    int block_N = N/numprocs;                                               //block_N is the number of vectors each process is assigned with, N of the sub matrix
    MPI_Request* 	reqs = (MPI_Request*)calloc(1, sizeof(MPI_Request));      //initialize MPI_Request for async communication MPI_Isend
    MPI_Status * 	stat=(MPI_Status*)malloc(sizeof(MPI_Status)*1);           //initialize MPI_Status for async communication MPI_Wait
    
    /* root process to initialize the matrix and ready to be distributed */
    if (myid == 0) {
      matrix = matrixInit(N,M);                                             //initiate the matrix for calculation
      printf("\n[Process #%d] %s [parallel]: INITIAL STATE of the matrix:\n\n", myid, time_local);
      printMatrix(matrix,N,M);
      printf("\n");
      
      output_num = N*(N+1)/2;                                               //total number of inner products results equals  N*(N+1)/2
      displs = (int *)malloc(numprocs * sizeof(int));                       //allocate displacement memory for MPI_Gatherv call to receive all results back to Process 0
      output_recv_count = (int *)malloc(numprocs * sizeof(int));            //allocate memeory for numbers of results sent back by each process 

      calculateDisplacementRecvCount(displs, output_recv_count, numprocs, block_N);
      output = (int *)malloc(output_num * sizeof(int));                     //allocate memory space for receiving the results
      output_ci = (int *)malloc(output_num * sizeof(int));                  //allocate memory space for receiving the results index Ci 
      output_cj = (int *)malloc(output_num * sizeof(int));                  //allocate memory space for receiving the results index Cj 
    }

    /* ❗logical for all processes starts from here */ 

    clock_gettime(CLOCK_REALTIME, &start);

    recvbuf = (int *)malloc(M * block_N * sizeof(int));                     //allocate memory space for the receive buff in each process
    MPI_Scatter (matrix, M * block_N, MPI_INT, recvbuf, M * block_N, MPI_INT, 0, MPI_COMM_WORLD);    //distribute the workload   

    comp_buf = (int *)malloc(M * block_N * sizeof(int));                    //allocate memory space for the computation buff in each process
    memcpy( comp_buf, recvbuf, M * block_N *sizeof(int) );                  //copy the values over to computation buff 

    coms_buf = (int *)malloc(M * block_N * sizeof(int));                    //duplicate the computation buffer to generate the comms buffer for sending/receiving asynchronously
    memcpy( coms_buf, comp_buf, M * block_N *sizeof(int) );                 //copy the values over to comms buff

    if(myid<numprocs/2){
    
      //first half of the processes calculate (numprocs/2+1) block pairs and generate (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2) results
      num_results = (numprocs/2)*block_N*block_N+(block_N*(block_N+1)/2);
      result=(int*)malloc(num_results * sizeof(int));
      c_i=(int*)malloc(num_results * sizeof(int));
      c_j=(int*)malloc(num_results * sizeof(int));
      index=0;
      //pairwise inner product calculation - block loop:
      for (h=0; h<(numprocs/2+1); h++){                                    //loop for (numprocs/2+1) steps for fisrt half of the processes
        if (h==0){                                                         //treat step 0 individually because it only calculate upper diagonal	
          //start pairwise inner product calculation
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=i; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
                c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
              }
              c_i[index]= i+myid*block_N;
              c_j[index] = j+myid*block_N+h*block_N;
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
          MPI_Wait(reqs, stat);
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
              result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
              index++;
            }
          }
        }
      }
    }
    else if(numprocs/2 <= myid && myid <numprocs) {
      //second half of the processes calculate (numprocs/2) block pairs and generate (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2) results
      num_results = (numprocs/2-1)*block_N*block_N + (block_N*(block_N+1)/2);
      result = (int*)malloc( num_results * sizeof(int) );                  //allocate memory for results in Compressed Row Storage format (V[k])
      c_i=(int*)malloc(num_results * sizeof(int));
      c_j=(int*)malloc(num_results * sizeof(int));
      index=0;
      //pairwise inner product calculation - block loop:
      for (h=0; h<(numprocs/2); h++){                                      //loop for (numprocs/2) steps for second half of the processes
        if (h==0){
          for (i=0; i < block_N; i ++){                                    //iterate through each vector
            for (j=i; j < block_N; j ++){                                  //iterate through each other vectors starting from itself
              c = 0;                                                       //reset result 'c' for the next pair
              for (k=0; k<M; k++){                                         //iterate through each element of both vectors
                c += recvbuf[i*M+k]*comp_buf[j*M+k];                       //calculate inner product "Cij"
              }
              c_i[index]= i+myid*block_N;
              c_j[index] = j+myid*block_N+h*block_N;
              result[index] = c;
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
          MPI_Wait(reqs, stat);
          memcpy( comp_buf, coms_buf, M * block_N *sizeof(int) );
          if (h<numprocs/2 && myid == numprocs/2){
            MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
          } 
          else if (h<(numprocs/2-1) && myid!=numprocs/2){
            MPI_Isend(coms_buf, M * block_N, MPI_INT, dst, 0, MPI_COMM_WORLD, reqs);  //send the data and gets into the next loop 
          }
          if(myid>numprocs/2 && h > (numprocs-1-myid)){                      //remap Cij number for balanced load
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
                result[index] = c;                                           //store calculated result for each pair in Compressed Row Storage format
                index++;
              }
            }
          }  
        }
    }
  }
  
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

      printf("[Process #%d] %s [parallel]: results row index are: ", myid, time_local);
      for (int i = 0; i < N; ++i) {
        printf("%d ", result_row_index[i]);
        }
      printf("\n\n");
      int* final_result=(int*)malloc(N*(N+1)/2 * sizeof(int));

      for (i=0; i<N*(N+1)/2; i++){
        index = result_row_index[output_ci[i]]+output_cj[i]-output_ci[i]; 
        final_result[index]=output[i];
      }
      clock_gettime(CLOCK_REALTIME, &finish);
      printf("[Process #%d] %s [parallel]: final results are: \n\n", myid, time_local);
      print_result(final_result,result_row_index,N);

      double time_spent_parallel = finish.tv_sec - start.tv_sec;
      time_spent_parallel += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      printf("[Process #%d] %s [parallel]: parallel computation time elapsed: [%f] seconds\n", myid, time_local, time_spent_parallel);
      
      /* Self checking serialization computation */
      int dif = 0;
      printf("[Process #%d] %s [parallel]: self check starting in 9 seconds...\n", myid, time_local);
      for (i=9;i>0;i--){
        printf("[Process #%d] %s [parallel]: %d...\r", myid, time_local,i);
        sleep(1);
        fflush(stdout);
      }
      printf("[Process #%d] %s [parallel][check]: self checking computation in progress...\n", myid, time_local);
      int *result_check = (int*)malloc(N*(N+1)/2 * sizeof(int));         //allocate memory for results in Compressed Row Storage format (V[k])
      int *result_check_row_index = (int*)malloc(N * sizeof(int));       //allocate memory for index array for the above results (R[i])
      result_check_row_index[0] = 0;                                //assign the first index to be 0
      for (i=1; i<N; i++){
        result_check_row_index[i] = result_check_row_index[i-1]+N-i+1;    //calculate the index for the first elemnt of each row and store in index array
      }

      /* start pairwise inner product calculation */
      clock_gettime(CLOCK_REALTIME, &start);
      calculateInnerProduct(matrix, result_check, result_check_row_index, N, M);
      checkResults(final_result,result_check,output_num,dif);
      clock_gettime(CLOCK_REALTIME, &finish);
      double time_spent_check = finish.tv_sec - start.tv_sec;
      time_spent_check += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      if (dif == 0){
        printf("[Process #%d] %s [parallel][check]: parallel computation matches self checking sequential computation result!\n", myid, time_local);
      }
      else {
        printf("[Process #%d] %s [parallel][check]: Oops...there are [%d] results different between parallel computation and sequential self checking computation result...\n", myid, time_local, dif);
      }
      printf("[Process #%d] %s [parallel][check]: sequential self checking computation time elapsed: [%f] seconds\n", myid, time_local, time_spent_check);
    }
  }
  MPI_Finalize();
  return 0;
}

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
void printMatrix(int *matrix,int N_vector, int M_length){
  int i,j;
  for (i=0; i<N_vector; i++){
    for (j=0; j<M_length; j++){
      printf("%d ", *(matrix + i*M_length + j));
    }
    printf("\n");
  }
  printf("\n");
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
  printf("\n");
};
void calculateInnerProduct(int *matrix, int *result, int *result_row_index, int N_vector, int M_length){
  int i,j,k;
  int index;                                              //index for result storage
  int c;                                                  //temp inner product result
  for (i=0; i<N_vector; i++){                             //iterate through each vector
    for (j=i; j<N_vector; j++){                           //iterate through each other vectors starting from itself
      c = 0;                                              //reset result 'c' for the next pair
      for (k=0; k<M_length; k++){                         //iterate through each element of both vectors
        c += matrix[i*M_length+k]*matrix[j*M_length+k];   //calculate inner product "Cij"
      }
      index = result_row_index[i]+j-i;                    //calculate the index for Cij in Compressed Row Storage format
      result[index] = c;                                  //store calculated result for each pair in Compressed Row Storage format
    }
  }
}
void calculateDisplacementRecvCount(int *displs, int *output_recv_count, int numprocs, int block_N){

  int i; 

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
}
int checkResults(int *final_result, int *result_check, int output_num, int dif){
 
  for(int i=0; i<output_num; i++){
      if(final_result[i]!=result_check[i]) {
        dif++;
      }
  }
  return dif;
}