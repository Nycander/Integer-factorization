#include <gmp.h>
#include "factor_list.h"

#ifndef POLLARD_H
#define POLLARD_H

int pollard(factor_list ** f, const mpz_t n);
int rho(mpz_t res, const mpz_t N);

#endif
