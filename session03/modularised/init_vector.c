#include <stddef.h>

void init_vector(double *vec, size_t len, ptrdiff_t incr)
{
    for(size_t i=0; i<len; ++i){
        vec[(ptrdiff_t)i*incr] = (double)i;
    }
}
