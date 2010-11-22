#include <gmp.h>
#include "factor_list.h"

#ifndef POLLARD_H
#define POLLARD_H

int pollard(factor_list ** f, const mpz_t n);
int rho(mpz_t res, const mpz_t N);
void f(mpz_t result, mpz_t x, const mpz_t N);
int floyd(const mpz_t N, mpz_t divisor);
int brent(const mpz_t N, mpz_t divisor);
#endif
