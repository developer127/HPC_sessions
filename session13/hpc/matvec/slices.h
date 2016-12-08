#ifndef INC_MATVEC_SLICES_H
#define INC_MATVEC_SLICES_H 1

#include <cassert>

namespace hpc { namespace matvec {

template <typename Index>
struct Slices
{
    const unsigned int nof_threads;
    const Index intervalSize;
    const Index offset;

    Slices(unsigned int nof_threads, Index isize):
        nof_threads(nof_threads), intervalSize(isize),
        offset ((intervalSize + nof_threads -1) / nof_threads)
    {
    }

    Index size(Index i)
    {
        assert(i<nof_threads);
        if(i < nof_threads -1) {
            return offset;
        } else {
            return (intervalSize % offset == 0)? offset : intervalSize % offset;
        }
    }
};

}}          // namespace hpc::matvec
#endif      // INC_MATVEC_SLICES_H
