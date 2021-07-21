#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


int main (int argc , char **argv) {
    int myid, numprocs;
    int token = 0;
    MPI_Status status;

    MPI_Init(&argc,&argv); /* Used to send the command line argumenys to all procs */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* Initializes the number of procs in the group specified by mpirun */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);	 /* Initialize the rank of this process in the group */
  
    if(myid == 0)
    {
	token = 1;
        printf("Hello World! This is process %i from %i - I have the token now!\n",myid,numprocs);
        MPI_Send(&token, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD);
        MPI_Recv(&token, 1, MPI_INT, numprocs-1, 1, MPI_COMM_WORLD, &status);
    } else {
        MPI_Recv(&token, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, &status);
        printf("Hello World! This is process %i from %i - I have the token now!\n",myid,numprocs);
        MPI_Send(&token, 1, MPI_INT, (myid+1)%numprocs, 1, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
