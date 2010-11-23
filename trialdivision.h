#ifndef TRIAL_DIVISION_H
#define TRIAL_DIVISION_H

#include <gmp.h>
#include "factor_list.h"

mpz_t * trial_division(factor_list ** f, unsigned int * primes, unsigned int prime_count, const mpz_t n);

#endif
