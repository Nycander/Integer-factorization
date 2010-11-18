#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "factor.h"
#include "pollard.c"

#define VERBOSE 1

int current_input_number = 0;

int main(int argc, char * argv[])
{
#if VERBOSE
	printf("Hyper mega global factoring program\n");
	printf("-----------------------------------\n");
	printf("> ");
#endif

	int limit = (argc == 2 ? atoi(argv[1]) : INT_MAX );

	mpz_t threshold;
	mpz_init_set_str(threshold, "83003906250000000", 10); // Brute-forced value @ kattis

	mpz_t num;
	while(++current_input_number <= limit)
	{
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
		}
		else
		{
			printf("fail\n\n");
		}

		mpz_clear(num);

	#if VERBOSE
		printf("> ");
	#endif
	}
	return 0;
}