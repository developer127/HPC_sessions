#include <stdio.h>
#include "vec.h"    //basic vector functions

double v[8];

int main()
{
    init_vector(v,sizeof(v),1);
    print_vector(v, sizeof(v),1);
    return 0;
}
