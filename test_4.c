#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]){

int myid, numprocs;
int recvbuf;
int* sendbuf = (int *)malloc(3 * sizeof(int)); 

MPI_Init (&argc, &argv);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank (MPI_COMM_WORLD, &myid);
if(myid == 0){

int row = 4;
int column = 5;

  int* matrix = (int *)malloc(row * column * sizeof(int)); 
  int i,j,k; 

  for (i=0; i<row; i++){
    for (j=0; j<column; j++){
      *(matrix + i*column + j) = rand()%10;
    }
  }

  //*(matrix + i*column + j) 指代这个地址/pointer上的数值
  /*
  idd@IDD-PC-NR200P:~/comp5426$ mpirun -n 1 test_4
  0
  0
  0
  0
  1015985232
  22081
  1015988368
  22081
  1
  0
  10
  0
  1015993248
  22081
  -426438156
  32570
  -320447504
  32570
  -320447504
  32570
  The matrix elements are:
  0 0 0 0 1015985232
  22081 1015988368 22081 1 0
  10 0 1015993248 22081 -426438156
  32570 -320447504 32570 -320447504 32570
  

  //(matrix + i*column + j) 指代这个地址/pointer本身，连贯的内存地址
  idd@IDD-PC-NR200P:~/comp5426$ mpirun -n 1 test_4
  0x558fd5f15440
  0x558fd5f15444
  0x558fd5f15448
  0x558fd5f1544c
  0x558fd5f15450
  0x558fd5f15454
  0x558fd5f15458
  0x558fd5f1545c
  0x558fd5f15460
  0x558fd5f15464
  0x558fd5f15468
  0x558fd5f1546c
  0x558fd5f15470
  0x558fd5f15474
  0x558fd5f15478
  0x558fd5f1547c
  0x558fd5f15480
  0x558fd5f15484
  0x558fd5f15488
  0x558fd5f1548c
  */


  printf("The matrix elements are:\n");
   for (i = 0; i < row; i++) {
      for (j = 0; j < column; j++) {
         printf("%d ", *(matrix + i*column + j)); 
      }
      printf("\n");
   }


}

//MPI_Scatter (sendbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
//printf ("rank: %i, value: %i\n", myid, recvbuf);


MPI_Finalize();
return 0;
}




