#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int *createMatrix (int nrows, int ncols) {
    int *matrix;
    int h, i, j;

    matrix = malloc(nrows*ncols*sizeof(int));

    for (h=0; h<nrows*ncols; h++) {
        matrix[h] = h+1;
    }

    printf("%p",matrix);
    return matrix;
}

void printArray (int *row, int nElements) {
    int i;
    for (i=0; i<nElements; i++) {
        printf("%d ", row[i]);
    }
    printf("\n");
}

int main (int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int p, id;
    MPI_Comm_size(MPI_COMM_WORLD, &p); // Get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &id); // Get own ID

    int *matrix;

    if (id == 0) {
        matrix = createMatrix(p, p); // Master process creates matrix
        printf("Initial matrix:\n");
        printArray(matrix, p*p);
    }

    int *procRow = malloc(sizeof(int) * p); // received row will contain p integers
    MPI_Scatter(matrix, p, MPI_INT, procRow, p, MPI_INT, 0, MPI_COMM_WORLD);

    printf("Process %d received elements: ", id);
    printf("%d",procRow[0]);
    printArray(procRow, p);

    MPI_Finalize();

    return 0;
}