#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int* matrix_init(int row, int column) {
  int* matrix; 
  int h,i,j; 

  matrix = (int *)malloc(row * column * sizeof(int)); 

  for (i=0; i<row; i++){
    for (j=0; j<column; j++){
      *(matrix + i*column + j) = rand()%10;
    }
  }


  printf("The INITIAL STATE of the matrix is:\n");
   for (i = 0; i < row; i++) {
      for (j = 0; j < column; j++) {
         printf("%d ", *(matrix + i*column + j)); 
      }
      printf("\n");
   }
  return matrix; 
}

int main(int argc, char *argv[]){

int myid, numprocs;
MPI_Init (&argc, &argv);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank (MPI_COMM_WORLD, &myid);
int* matrix,* recvbuf;

if(myid == 0){

matrix = matrix_init(4,5);

}

recvbuf = (int *)malloc(5 * sizeof(int)); 
MPI_Scatter (matrix, 5, MPI_INT, recvbuf, 5, MPI_INT, 0, MPI_COMM_WORLD);
printf ("rank: %d, value: %d\n", myid, recvbuf[0]);


MPI_Finalize();
return 0;
}




