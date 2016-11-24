#ifndef INC_GEMM_BLOCKED_HPP
#define INC_GEMM_BLOCKED_HPP 1

#include <complex>
#include <type_traits>
#include "ulmblas.hpp"

namespace blocked {

/* Defining the hard coded block sizes
 * for different datatypes
 */

template <typename T>
struct BlockSize
{
    static constexpr int MC = 64;
    static constexpr int KC = 64;
    static constexpr int NC = 256;
    static constexpr int MR = 2;
    static constexpr int NR = 2;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<float>
{
    static constexpr int MC = 256;
    static constexpr int KC = 512;
    static constexpr int NC = 4096;
    static constexpr int MR = 8;
    static constexpr int NR = 8;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<double>
{
    static constexpr int MC = 256;
    static constexpr int KC = 256;
    static constexpr int NC = 4096;
    static constexpr int MR = 4;
    static constexpr int NR = 8;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<std::complex<float>>
{
    static constexpr int MC = 256;
    static constexpr int KC = 256;
    static constexpr int NC = 4096;
    static constexpr int MR = 4;
    static constexpr int NR = 8;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<std::complex<double>>
{
    static constexpr int MC = 256;
    static constexpr int KC = 128;
    static constexpr int NC = 4096;
    static constexpr int MR = 4;
    static constexpr int NR = 4;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};


/* The actual block algorithms */

template <typename T, typename Size, typename Index>
void
pack_A(Size mc, Size kc,
       const T *A, Index incRowA, Index incColA,
       T *p)
{
    static constexpr int MR = BlockSize<T>::MR;
    Size mp = (mc+MR-1) / MR;

    for (Size j=0; j<kc; ++j) {
        for (Size l=0; l<mp; ++l) {
            for (Size i0=0; i0<MR; ++i0) {
                Size i  = l*MR + i0;
                Size nu = l*MR*kc + j*MR + i0;
                p[nu]          = (i<mc) ? A[i*incRowA+j*incColA]
                                        : T(0);
            }
        }
    }
}

template <typename T, typename Size, typename Index>
void
pack_B(Size kc, Size nc,
       const T *B, Index incRowB, Index incColB,
       T *p)
{
    static constexpr int NR = BlockSize<T>::NR;
    Size np = (nc+NR-1) / NR;

    for (Size l=0; l<np; ++l) {
        for (Size j0=0; j0<NR; ++j0) {
            for (Size i=0; i<kc; ++i) {
                Size j  = l*NR+j0;
                Size nu = l*NR*kc + i*NR + j0;
                p[nu]          = (j<nc) ? B[i*incRowB+j*incColB]
                                        : T(0);
            }
        }
    }
}

template <typename T, typename Size, typename Index>
void
ugemm(Size kc, T alpha,
      const T *A, const T *B,
      T beta,
      T *C, Index incRowC, Index incColC)
{
    const Size MR = BlockSize<T>::MR;
    const Size NR = BlockSize<T>::NR;
    T P[BlockSize<T>::MR*BlockSize<T>::NR];

    for (Size l=0; l<MR*NR; ++l) {
        P[l] = T(0);
    }
    for (Size l=0; l<kc; ++l) {
        for (Size j=0; j<NR; ++j) {
            for (Size i=0; i<MR; ++i) {
                P[i+j*MR] += A[i+l*MR]*B[l*NR+j];
            }
        }
    }
    if (beta!=T(1)) {
        if (beta!=T(0)) {
            for (Size j=0; j<NR; ++j) {
                for (Size i=0; i<MR; ++i) {
                    C[i*incRowC+j*incColC] *= beta;
                }
            }
        } else {
            for (Size j=0; j<NR; ++j) {
                for (Size i=0; i<MR; ++i) {
                    C[i*incRowC+j*incColC] = T(0);
                }
            }
        }
    }
    for (Size j=0; j<NR; ++j) {
        for (Size i=0; i<MR; ++i) {
            C[i*incRowC+j*incColC] += alpha*P[i+j*MR];
        }
    }
}

template <typename T, typename Size, typename Index>
void
mgemm(Size mc, Size nc, Size kc,
      T alpha,
      const T *A, const T *B,
      T beta,
      T *C, Index incRowC, Index incColC)
{
    const Size MR = BlockSize<T>::MR;
    const Size NR = BlockSize<T>::NR;
    T C_[BlockSize<T>::MR*BlockSize<T>::NR];

    const Size mp  = (mc+MR-1) / MR;
    const Size np  = (nc+NR-1) / NR;
    const Size mr_ = mc % MR;
    const Size nr_ = nc % NR;

    for (Size j=0; j<np; ++j) {
        const Size nr = (j!=np-1 || nr_==Size(0)) ? NR : nr_;

        for (Size i=0; i<mp; ++i) {
            const Size mr = (i!=mp-1 || mr_==Size(0)) ? MR : mr_;

            if (mr==MR && nr==NR) {
                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR],
                      beta,
                      &C[i*MR*incRowC+j*NR*incColC],
                      incRowC, incColC);
            } else {
                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR],
                      T(0),
                      C_, Index(1), MR);
                ulmBLAS::gescal(mr, nr, beta,
                                &C[i*MR*incRowC+j*NR*incColC],
                                incRowC, incColC);
                ulmBLAS::geaxpy(mr, nr, T(1), C_, Index(1), MR,
                                &C[i*MR*incRowC+j*NR*incColC],
                                incRowC, incColC);
            }
        }
    }
}

template <typename T, typename Size, typename Index>
void
gemm(Size m, Size n, Size k,
     T alpha,
     const T *A, Index incRowA, Index incColA,
     const T *B, Index incRowB, Index incColB,
     T beta,
     T *C, Index incRowC, Index incColC)
{
    const Size MC = BlockSize<T>::MC;
    const Size NC = BlockSize<T>::NC;
    const Size KC = BlockSize<T>::KC;

    const Size mb = (m+MC-1) / MC;
    const Size nb = (n+NC-1) / NC;
    const Size kb = (k+KC-1) / KC;

    const Size mc_ = m % MC;
    const Size nc_ = n % NC;
    const Size kc_ = k % KC;

    if (alpha==T(0) || k==T(0)) {
        ulmBLAS::gescal(m, n, beta, C, incRowC, incColC);
        return;
    }

    T *A_ = new T[MC*KC];
    T *B_ = new T[KC*NC];

    for (Size j=0; j<nb; ++j) {
        const Size nc = (j!=nb-1 || nc_==Size(0)) ? NC : nc_;

        for (Size l=0; l<kb; ++l) {
            const Size   kc    = (l!=kb-1 || kc_==Size(0)) ? KC : kc_;
            const T beta_ = (l==Size(0)) ? beta : T(1);

            pack_B(kc, nc,
                   &B[l*KC*incRowB+j*NC*incColB],
                   incRowB, incColB,
                   B_);

            for (Size i=0; i<mb; ++i) {
                Size mc = (i!=mb-1 || mc_==Size(0)) ? MC : mc_;

                pack_A(mc, kc,
                       &A[i*MC*incRowA+l*KC*incColA],
                       incRowA, incColA,
                       A_);

                mgemm(mc, nc, kc,
                      alpha, A_, B_, beta_,
                      &C[i*MC*incRowC+j*NC*incColC],
                      incRowC, incColC);
            }
        }
    }
    delete [] A_;
    delete [] B_;
}

} // namespace blocked

#endif // INK_GEMM_BLOCKED_HPP
