#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "vec.h"
#include "walltime.h"


int main() {
    double t0, t1, t2;
    size_t length = 8192;
    printf("     len t1 (separate) t2 (interleaving)      t2/t1\n");

    while(length <= 67108864){
        printf("%8zd", length);
        double* vector1 = (double*) malloc(sizeof(double) * length * 2);
        if (!vector1){
            printf("Error allocating memory\n");
            return 1;
        }
        double* vector2 = (double*) malloc(sizeof(double) * length * 2);
        if (!vector2){
            printf("Error allocating memory\n");
            return 1;
        }
        // separate vectors
        t0 = walltime();
        init_vector(vector1, length, 1);
        init_vector(vector1+length, length, 1);
        t1 = walltime() - t0;
        printf(" %12.4lf", t1);
        // interleaved vectors
        t0 = walltime();
        init_vector(vector2, length, 2);
        init_vector(vector2+1, length, 2);
        t2 = walltime() - t0;
        free(vector1);
        free(vector2);
        printf(" %12.4lf %16.2lf\n", t2, t2/t1);
        length *=2;
    }
    return 0;
}
