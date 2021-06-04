#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
   int var = 20;   /* actual variable declaration */
   int c = 100;
   int  *ip;        /* pointer variable declaration */

   ip = (int*)malloc( 2 * sizeof(int) );;  /* store address of var in pointer variable*/

   /* address stored in pointer variable */
   printf("Address stored in ip variable: %p\n", &ip[1] );

   /* access the value using the pointer */
   printf("Value of *ip variable: %d\n", ip[1] );

   ip[1]=c; 
   /* address stored in pointer variable */
   printf("Address stored in ip variable: %p\n", &ip[0] );
   /* address stored in pointer variable */
   printf("Address stored in ip variable: %p\n", &ip[1] );

   /* access the value using the pointer */
   printf("Value of *ip variable: %d\n", ip[1] );

   return 0;

}

