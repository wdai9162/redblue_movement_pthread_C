#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


int main (int argc , char *argv[])
{
	int myid, numProcs;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


	int x, y, z;
	switch(myid) {
		case 0: x=1; y=4; z=2;
		MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Send(&y, 1, MPI_INT, 2, 2, MPI_COMM_WORLD);
		MPI_Recv(&z, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		MPI_Allreduce(&y, &x, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		break;
		case 1: x=3; y=8; z=6;
		MPI_Bcast(&y, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Send(&z, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
		MPI_Send(&x, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Allreduce(&z, &y, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		break;
		case 2: x=6; y=7; z=8;
		MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Recv(&y, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		MPI_Recv(&z, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		MPI_Allreduce(&x, &y, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		break;
	}
    
	printf("Process #%d: x=%d y=%d z=%d \n",myid,x,y,z);

	MPI_Finalize();

	return 0;
}

