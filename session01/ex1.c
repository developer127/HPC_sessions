#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

void init_vector(double* vec, size_t len, ptrdiff_t incr)
{
    for(size_t i=0; i<len; ++i){
        vec[i*incr] = i;
    }
}

void print_vector(double* vec, size_t len, ptrdiff_t inc)
{
    for(size_t i = 0; i < len; ++i){
        printf("%4.1lf ",vec[i*inc]);
    }
    printf("\n");
}

/* Main program */

//double v[8];

int main()
{
    printf("Define the length of the vector: ");

    size_t length;
    if(scanf("%zu", &length) != 1){
        return 1;
    }
    printf("\n confirmed the length is: %zu\n",length);

    double* v = (double*) malloc(sizeof(double) * length);
    if (!v){
        printf("Error allocating  memory\n");
        return 1;
    }
    init_vector(v, length, 1);
    print_vector(v, length, 1);

    init_vector(v, length/2, 1);
    print_vector(v, length/2, 1);

    init_vector(&v[length/2], length/2, 1);
    print_vector(&v[length/2], length/2, 1);

    print_vector(v, length, 1);

    init_vector(v, length/2, 2);
    print_vector(v, length/2, 2);

    init_vector(&v[1], length/2, 2);
    print_vector(&v[1], length/2, 2);

    print_vector(v, length, 1);

    free(v);

    return 0;
}
