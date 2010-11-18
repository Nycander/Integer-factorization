#ifndef FACTOR_H
#define FACTOR_H

#include <gmp.h>

typedef struct factor_list
{
	mpz_t * value;
	struct factor_list * next;
} factor_list;

factor_list * factor_list_add(factor_list ** f, mpz_t * v);

#endif
