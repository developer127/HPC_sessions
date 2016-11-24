#ifndef INC_ULMBLAS_HPP
#define INC_ULMBLAS_HPP

namespace ulmBLAS {

//------------------------------------------------------------------------------
// A <- B
//------------------------------------------------------------------------------

template <typename T, typename Size, typename Index>
void
gecopy(Size m, Size n,
       const T *A, Index incRowA, Index incColA,
       T *B, Index incRowB, Index incColB)
{
    for (Size j=0; j<n; ++j) {
        for (Size i=0; i<m; ++i) {
            B[i*incRowB+j*incColB] = A[i*incRowA+j*incColA];
        }
    }
}

//------------------------------------------------------------------------------
// Y <- Y + alpha*X
//------------------------------------------------------------------------------

template <typename T, typename Size, typename Index>
void
geaxpy(Size m, Size n,
       T alpha,
       const T *X, Index incRowX, Index incColX,
       T       *Y, Index incRowY, Index incColY)
{
    for (Size i=0; i<m; ++i) {
        for (Size j=0; j<n; ++j) {
            Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
        }
    }
}

//------------------------------------------------------------------------------
// A <- alpha * A
//------------------------------------------------------------------------------

template <typename T, typename Size, typename Index>
void
gescal(Size m, Size n,
       T alpha,
       T *X, Index incRowX, Index incColX)
{
    if (alpha!=T(1)) {
        for (Size i=0; i<m; ++i) {
            for (Size j=0; j<n; ++j) {
                X[i*incRowX+j*incColX] *= alpha;
            }
        }
    }
}

}   // ulmBLAS namespace

#endif
