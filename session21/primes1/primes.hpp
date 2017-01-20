#ifndef PRIMES_HPP
#define PRIMES_HPP

#include <hpc/gmp/integer.h>

bool search_prime_constellation(
	 hpc::gmp::Integer& start_interval,
	 hpc::gmp::Integer& end_interval,
	 unsigned int k,           // size of prime constellation
	 unsigned int offsets[],   // k-1 offsets
	 hpc::gmp::Integer& first_prime);    // found result (if returning true)

#endif
