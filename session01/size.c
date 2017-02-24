#include <stdlib.h>
#include <stdio.h>

int main()
{
    double x[8];
    double *y;

    for(size_t i = 0; i<8; ++i){
        x[i] = i*i;
    }

    y = x;

    printf("Array x has size: %zu\n", sizeof(x));
    printf("The pointer pointing to the start of the array doesn't know\n"
           "anything about how many numbers there exist: %zu\n", sizeof(y));

    printf("y = ");
    for(size_t i = 0; i<8; ++i){
        printf("%4.1f ", y[i]);
    }
    printf("\n");
}
