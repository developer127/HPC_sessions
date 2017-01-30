#ifndef HPC_GMP_INTEGER_H
#define HPC_GMP_INTEGER_H

#include <gmp.h>
#include <gmpxx.h>

namespace hpc { namespace gmp {

using Integer = ::mpz_class;

struct ExportedInteger {
   using Word = unsigned int;
   static constexpr auto numb = 8 * sizeof(Word);
   const int sign;
   const unsigned int len;
   Word* words; /* exported representation of an GMP integer */

   ExportedInteger(const Integer& integer) :
	 sign(mpz_sgn(integer.get_mpz_t())),
	 len((mpz_sizeinbase(integer.get_mpz_t(), 2) + numb - 1) / numb),
	 words(new Word[len]) {
      mpz_export(words, nullptr, 1, sizeof(Word), 1, 0,
	 integer.get_mpz_t());
   }
   ExportedInteger(int sign, unsigned int len) :
      sign(sign), len(len), words(new Word[len]) {
   }
   unsigned int length() const {
      return len;
   }
   int sgn() const {
      return sign;
   }
   Integer get() const {
      Integer val;
      if (sign == 0) return val;
      mpz_import(val.get_mpz_t(), len, 1, sizeof(Word), 1, 0,
	 (const void*) words);
      if (sign > 0) {
	 return val;
      } else {
	 return -val;
      }
   }
   ~ExportedInteger() {
      delete words;
   }
};

} } // namespace gmp, hpc

#endif
