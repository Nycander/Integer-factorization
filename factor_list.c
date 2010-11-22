
#include <stdlib.h>
#include <gmp.h>
#include "factor_list.h"

factor_list * factor_list_add(factor_list ** f, mpz_t * v)
{
	factor_list * nF = malloc(sizeof(factor_list));
	nF->next = *f;
	nF->value = v;
	*f = nF;

	return *f;
};
