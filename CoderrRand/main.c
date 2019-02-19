#include <stdio.h>
#include <stdlib.h>
#include "coderRand.h"
#include "coderRand_initialize.h"
#include "coderRand_terminate.h"
int main()
{
    coderRand_initialize();
    printf("coderRand=%g\n", coderRand());
    coderRand_terminate();
    
    puts("Press enter to quit:");
    getchar();

    return 0;
}