#ifndef QS_H
#define QS_H

#include <gmp.h>


void quadraticSieve(const mpz_t num);
factor_list * sieving(const mpz_t num);


#endif
