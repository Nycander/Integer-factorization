#ifndef FACTOR_H
#define FACTOR_H

#include <gmp.h>

typedef struct factor
{
	mpz_t * value;
	struct factor * next;
} factor;

factor * factor_add(factor ** f, mpz_t * v);

#endif
