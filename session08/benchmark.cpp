#include "benchmark.hpp"

//------------------------------------------------------------------------------
// Functions for benchmarking
//------------------------------------------------------------------------------
namespace bench{


double
walltime()
{
    struct tms    ts;
    static double ClockTick=0.0;

    if (ClockTick==0.0) {
        ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
    }
    return ((double) times(&ts)) * ClockTick;
}

} // namespace bench
