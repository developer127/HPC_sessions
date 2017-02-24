#include <stdio.h>
#include "vec.h"    //basic vector functions

double v[8];

int main()
{
     printf("Define the length of the vector: ");

    size_t length;
    if(scanf("%zu", &length) != 1){
        return 1;
    }
    printf("\n confirmed the length is: %zu\n",length);

    init_vector(v, length, 1);
    print_vector(v, length, 1);
    return 0;
}
