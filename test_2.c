#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]){

int myid, numprocs;

MPI_Init (&argc, &argv);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  if(myid == 0)
	{
// declare a pointer variable to point to allocated heap space
int    *p_array;
double *d_array;

// call malloc to allocate that appropriate number of bytes for the array

p_array = (int *)malloc(sizeof(int)*50);      // allocate 50 ints
d_array = (double *)malloc(sizeof(double)*100);  // allocate 100 doubles


// use [] notation to access array buckets 
// (THIS IS THE PREFERED WAY TO DO IT)
for(int i=0; i < 50; i++) {
  p_array[i] = 0;
}

// you can use pointer arithmetic (but in general don't)
double *dptr = d_array;    // the value of d_array is equivalent to &(d_array[0])
for(int i=0; i < 50; i++) {
  *dptr = 0;
  dptr++;
}

   
	}
 	MPI_Finalize();
	return 0;
}

