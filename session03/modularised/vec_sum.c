#include <stddef.h>

double vec_sum(double * vec, size_t length, ptrdiff_t inc)
{
    double sum = 0;
    for(size_t i = 0; i<length; ++i){
        sum += vec[i*inc];
    }
    return sum;
}
