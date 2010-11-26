#ifndef QS_H
#define QS_H

#include <gmp.h>
#include "factor_list.h"

int quadratic_sieve(factor_list ** result, const mpz_t num);
void add_possible_factor_to_array(mpz_t factor, const mpz_t ofNumber, mpz_t array[], int index);

#endif
