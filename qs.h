#ifndef QS_H
#define QS_H

#include <gmp.h>
#include "factor_list.h"
#include "settings.h"
#include "primes.h"

int quadratic_sieve(factor_list ** result, const mpz_t num);
factor_list * sieving(const mpz_t num);

#endif
