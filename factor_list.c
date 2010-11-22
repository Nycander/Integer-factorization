#include <stdlib.h>
#include <stdio.h>
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


void factor_list_print(factor_list * f)
{
	while(*(f->value) != NULL)
	{
		gmp_printf("%Zd\n", f->value);
		f = f->next;
	}
	printf("\n");
};