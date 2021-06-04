#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

bool changeBool(bool result);
int main(int argc, char *argv[])
{
   bool result=false;
   bool *ip;
   ip = &result; 


   result = changeBool(result);
   fputs(result ? "true\n" : "false\n", stdout);
   result = changeBool(result);
   fputs(result ? "true\n" : "false\n", stdout);

   return 0;
}

bool changeBool(bool result)
{
    if (result) {
        result = false;
    } else result = true;
    return result;
}


