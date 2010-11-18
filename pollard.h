#include <gmp.h>
#include "factor.h"

#ifndef POLLARD_H
#define POLLARD_H

int pollard(factor ** f, const mpz_t n);
int rho(mpz_t res, const mpz_t N);

#endif
