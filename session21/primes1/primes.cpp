#include "primes.hpp"

using namespace hpc::gmp;

bool search_prime_constellation(
        Integer& start_interval,
        Integer& end_interval,
        unsigned int k,         // size of prime constellation
        unsigned int offsets[], // k-1 offsets
        Integer& first_prime)   // found result (if returning true)
{
    Integer p = start_interval;
    Integer q;
    bool found = false;
    for(;;) {
        // lookout for the next prime in the interval
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        if (p > end_interval) break;
        // p is apparently prime, check for an constellation
        found = true;
        for (int i = 0; i < k-1; ++i) {
            unsigned int offset = offsets[i];
            q = p + offset;
            if (mpz_probab_prime_p(q.get_mpz_t(), 10) < 1) {
                found = false;
                break; // not a prime
            }
        }
        if (found) break;
    }
    if (found) {
        first_prime = p;
    }
    return found;
}
