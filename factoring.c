#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "settings.h"
#include "factor_list.h"
#include "trialdivision.h"
#include "pollard.h"

int current_input_number = 0;

void factor(mpz_t n)
{
	factor_list * factors = malloc(sizeof(factor_list));
	factors->value = NULL;
	factors->next = NULL;

#if USE_TRIAL_DIVISION
	// Exhaust trivial primes with trial division
	int trivialPrimesCount = 0;
	while(1)
	{
		mpz_t * newPointer = trial_division(&factors, n);
		if (newPointer == 0)
			break;

		n = *newPointer;
		trivialPrimesCount++;
	}
#if VERBOSE
	gmp_printf("\tExhausted trivial primes after %d iterations; n = %Zd\n", trivialPrimesCount, n);
#endif
#endif

	if (pollard(&factors, n))
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
}


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

		factor(num);

		mpz_clear(num);
	}
	return 0;
}