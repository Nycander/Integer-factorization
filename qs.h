#ifndef QS_H
#define QS_H

#include <gmp.h>
#include "factor_list.h"

int quadratic_sieve(factor_list ** result, const mpz_t num);
int try_adding_factor_to_result(factor_list ** result, mpz_t factor, const mpz_t ofNumber, mpz_t visited[], int visited_length);

#endif
