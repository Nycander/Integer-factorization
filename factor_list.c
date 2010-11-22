#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "factor_list.h"

void factor_list_print(factor_list * f)
{
	while(*(f->value) != NULL)
	{
		gmp_printf("%Zd\n", f->value);
		f = f->next;
	}
	printf("\n");
};

