#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    int myid;
    int numprocs;
	int inputSize;
	int *sendArray;
	int *recvbuf;
	int commscount;
	int *result;
	
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

	
    if(myid == 0)
    {   
		printf("Please enter the number of elements to read in\n");
		//read from the standard input an int to numElements
		scanf ("%d", &inputSize);
		sendArray = (int *)malloc(inputSize*sizeof(int));
		result = (int *)malloc(inputSize*sizeof(int));
		for (int i=0; i < inputSize; i++) sendArray[i]=0;
		for (int i=0; i < inputSize; i++) printf("%d \n",sendArray[i]);
		commscount = inputSize/numprocs;
    }
	MPI_Bcast(&commscount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("process#%d: commscount is %d\n",myid, commscount);
	
	recvbuf = (int *)malloc(commscount * sizeof(int));
	printf("process#%d: recvbuf address is %p\n",myid, recvbuf);

	MPI_Scatter(sendArray, commscount, MPI_INT, recvbuf, commscount, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i=0; i < commscount; i++) recvbuf[i]=myid;

	MPI_Gather(recvbuf, commscount, MPI_INT, result, commscount, MPI_INT, 0, MPI_COMM_WORLD);
	if(myid==0) for (int i=0; i < inputSize; i++) printf("result %d \n",result[i]);


    MPI_Finalize();
    return 0;
}
