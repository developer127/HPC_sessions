#include <stdio.h>

void print_vector(double *vec, size_t len; ptrdiff_t inc)
{
    for(size_t i = 0; i < len; ++i){
        printf("%4f ",vec[i*inc]);
    }
    printf("\n");
}
