#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]){

int myid, numprocs;
int recvbuf;
int sendbuf[3];

MPI_Init (&argc, &argv);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank (MPI_COMM_WORLD, &myid);
if(myid == 0){

sendbuf[0] = 3;
sendbuf[1] = 5;
sendbuf[2] = 7;


}

MPI_Scatter (sendbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
printf ("rank: %i, value: %i\n", myid, recvbuf);

MPI_Finalize();
return 0;
}




