#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#define VERBOSE 0

#include "factor.h"
#include "pollard.h"

int current_input_number = 0;

int main(int argc, char * argv[])
{
#if VERBOSE
	printf("Hyper mega global factoring program\n");
	printf("-----------------------------------\n");
#endif

	int limit = (argc == 2 ? atoi(argv[1]) : INT_MAX );

	mpz_t num;
	while(++current_input_number <= limit)
	{
#if VERBOSE
		printf("> ");
#endif

		mpz_init(num);
		if (gmp_scanf("%Zd", num) <= 0)
		{
			mpz_clear(num);
			break;
		}

		factor * factors = malloc(sizeof(factor));
		factors->value = NULL;
		factors->next = NULL;
		if (pollard(&factors, num))
		{
			while(*(factors->value) != NULL)
			{
				gmp_printf("%Zd\n", factors->value);
				factors = factors->next;
			}
			printf("\n");
		}
		else
		{
			printf("fail\n\n");
		}

		// TODO: free the linked list!

		mpz_clear(num);

	}
	return 0;
}
