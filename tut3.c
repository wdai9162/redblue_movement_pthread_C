#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*
1D torus communication: 

0 s ---> 1 r 
1 s ---> 2 r 
2 s ---> 3 r
3 s ---> 4 r
4 s ---> 0 r

*/
int main(int argc, char **argv) {
    int myid;
    int numprocs;
    int *buf;
	MPI_Status status;  
	
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
   
    buf = malloc(4*sizeof(int));
    buf[0] = buf[3] = -1; 
    buf[1] = myid; 
    buf[2] = myid + numprocs;
    
    /*
    if(myid == 0)
    {   
        MPI_Send(&buf[2], 1, MPI_INT, (myid+1)%numprocs, 1, MPI_COMM_WORLD);
        MPI_Recv(&buf[0], 1, MPI_INT, (myid-1+numprocs)%numprocs, 1, MPI_COMM_WORLD, &status);

        MPI_Send(&buf[1], 1, MPI_INT, (myid-1+numprocs)%numprocs, 1, MPI_COMM_WORLD);
        MPI_Recv(&buf[3], 1, MPI_INT, (myid+1)%numprocs, 1, MPI_COMM_WORLD, &status);
    } else {
        MPI_Recv(&buf[0], 1, MPI_INT,(myid-1+numprocs)%numprocs, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&buf[2], 1, MPI_INT, (myid+1)%numprocs, 1, MPI_COMM_WORLD);

        MPI_Recv(&buf[3], 1, MPI_INT, (myid+1)%numprocs, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&buf[1], 1, MPI_INT, (myid-1+numprocs)%numprocs, 1, MPI_COMM_WORLD);

    }
    */

    //用MPI_Sendrecv重写这个逻辑就不需要区分myid=0了，所有process都自动收发
    MPI_Sendrecv(&buf[2], 1, MPI_INT, (myid+1)%numprocs, 1, &buf[0], 1, MPI_INT, (myid-1+numprocs)%numprocs, 1, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&buf[1], 1, MPI_INT, (myid-1+numprocs)%numprocs, 1, &buf[3], 1, MPI_INT, (myid+1)%numprocs , 1, MPI_COMM_WORLD, &status);

    printf("process %i: %i, %i, %i, %i\n", myid, buf[0], buf[1], buf[2], buf[3]);
    MPI_Finalize();
    return 0;
}
