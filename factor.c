
#include <stdlib.h>
#include <gmp.h>
#include "factor.h"

factor * factor_add(factor ** f, mpz_t * v)
{
	factor * nF = malloc(sizeof(factor));
	nF->next = *f;
	nF->value = v;
	*f = nF;

	return *f;
};
