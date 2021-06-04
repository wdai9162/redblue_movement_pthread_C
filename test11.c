#include <stdlib.h>
#include <stdio.h>

void addOne(int* n);

int main(int argc, char *argv[])
{
   int y=0;
   int *ip;
   ip = &y; 

   addOne(&y);
   printf("%d \n", y);
   addOne(&y);
   printf("%d \n", y);

   return 0;
}

void addOne(int* n)
{
    *n+=1;
}


