#include <mpi.h>
#include <iostream>
#include <cassert>
#include <printf.hpp>
#include "primes.hpp"
#include <hpc/gmp/integer.h>
#include <hpc/aux/slices.h>

using namespace hpc::aux;
using namespace hpc::gmp;

void send_integer(const Integer& value, int dest, int tag) {
   assert(dest >= 0);
   ExportedInteger exp_value(value);
   int len = (int) exp_value.length();
   int header[2] = {exp_value.sgn(), len};
   MPI_Send(header, 2, MPI_INT, dest, tag, MPI_COMM_WORLD);
   MPI_Send(exp_value.words, len, MPI_INT, dest, tag, MPI_COMM_WORLD);
}

void send_finish(int dest) {
   MPI_Send(nullptr, 0, MPI_INT, dest, 1, MPI_COMM_WORLD);
}

bool receive_integer(Integer& value, int& source, int& tag) {
   int lenval;
   MPI_Status status;
   int header[2];
   MPI_Recv(header, 2, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
   tag = status.MPI_TAG;
   source = status.MPI_SOURCE;
   if (tag) return false;
   int sgn = header[0]; unsigned int len = header[1];
   ExportedInteger exp_value(sgn, len);
   MPI_Recv(exp_value.words, len, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
   value = exp_value.get();
   return true;
}

char* progname;

void usage() {
   fmt::printf(std::cerr, "Usage: %s N1 N2 {n_i}\n", progname);
   exit(1);
}

void primes_master(int nofworkers, int argc, char** argv) {
   progname = *argv++; --argc;
   if (argc < 3) usage();
   Integer start(*argv++); --argc;
   Integer end(*argv++); --argc;

   int k = argc + 1;
   unsigned int* offsets = new unsigned int[k-1];
   for (int i = 0; i < k-1; ++i) {
      offsets[i] = atoi(*argv++); --argc;
   }

   MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
   MPI_Bcast(offsets, k-1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

   for (int worker = 1; worker <= nofworkers; ++worker) {
      send_integer(start, worker, 0);
      send_integer(end, worker, 0);
   }
   int running_workers = nofworkers;
   while (running_workers > 0) {
      int source = MPI_ANY_SOURCE; int tag = MPI_ANY_TAG;
      Integer prime;
      bool ok = receive_integer(prime, source, tag);
      if (!ok) {
	 fmt::printf("%d has finished\n", source);
	 --running_workers; continue;
      }
      fmt::printf("%d: %d\n", source, prime);
   }
}

void primes_worker(int nofworkers, int rank) {
   int k;
   MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
   unsigned int* offsets = new unsigned int[k-1];
   MPI_Bcast(offsets, k-1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

   Integer global_start, global_end;
   int source = 0, tag = 0;
   receive_integer(global_start, source, tag);
   receive_integer(global_end, source, tag);

   Slices<Integer> slices(nofworkers, global_end - global_start);
   Integer start = slices.offset(rank);
   Integer end = start + slices.size(rank);

   Integer prime;
   while (search_prime_constellation(start, end, k, offsets, prime)) {
      send_integer(prime, 0, 0);
      start = prime;
      start += 1;
   }
   send_finish(0);
}

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);

   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int nofworkers; MPI_Comm_size(MPI_COMM_WORLD, &nofworkers);
   --nofworkers; assert(nofworkers > 0);

   if (rank == 0) {
      primes_master(nofworkers, argc, argv);
   } else {
      primes_worker(nofworkers, rank-1);
   }

   MPI_Finalize();
}
