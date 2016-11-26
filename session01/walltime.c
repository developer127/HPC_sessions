#include <stdio.h>
#include <stddef.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdlib.h>


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


/* return real time in seconds since start of the process */
double walltime() {
    static int ticks_per_second = 0;
    if (!ticks_per_second) {
        ticks_per_second = sysconf(_SC_CLK_TCK);
    }
    struct tms timebuf;
    /* times returns the number of real time ticks passed since start */
    return (double) times(&timebuf) / ticks_per_second;
}

int main() {
    double t0, t1;
    size_t length = 8192;

    while(length <= 67108864){
        t0 = walltime();
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
        init_vector(vector1, length, 1);
        init_vector(vector1+length, length, 1);
        init_vector(vector2, length, 2);
        init_vector(vector2+1, length, 2);
        t1 = walltime() - t0;
        free(vector1);
        free(vector2);
        printf("for the vector of size %zu it took %.4lf seconds\n",
                length*2,t1);
        length *=2;
    }
    return 0;
}
