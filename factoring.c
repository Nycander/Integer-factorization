#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "settings.h"
#include "factor_list.h"
#include "trialdivision.h"
#include "pollard.h"
#include "qs.h"
#include "primes.h"

int current_input_number = 0;

void factor(mpz_t n)
{
	factor_list * factors = malloc(sizeof(factor_list));
	factors->value = NULL;
	factors->next = NULL;

	mpz_t original_n;
	mpz_init_set(original_n, n);


	#if USE_TRIAL_DIVISION

	#if VERBOSE
	gmp_printf(" :: Using trial-division with the first %d primes on %Zd...\n\t%Zd", primes_count, n, n);
	#endif

	// Exhaust trivial primes with trial division
	n = *trial_division(&factors, primes, primes_count, n);

	#if VERBOSE
	factor_list * tmp = factors;
	while(tmp->value != NULL)
	{
		gmp_printf(" / %Zd", *(tmp->value));
		tmp = tmp->next;
	}
	gmp_printf(" = %Zd\n", n);

	gmp_printf(" :: Exhausted all trivial primes, the number is now %Zd\n", n);
	#endif
	#endif



	while (mpz_sizeinbase(n, 2) >= USE_QUADRATIC_SIEVE_BIT_THRESHOLD)
	{
		#if VERBOSE
		gmp_printf(" :: Using QS since %Zd is %d bits long. (Threshold for QS is %d)\n", n, mpz_sizeinbase(n, 2), USE_QUADRATIC_SIEVE_BIT_THRESHOLD);
		#endif

		int qs_result = quadratic_sieve(&factors, n);

		if (qs_result == 0)
		{
			printf("fail\n\n");
			return;
		}

		#if VERBOSE
		gmp_printf(" :: Dividing the number %Zd with all found factors from QS... \n \t%Zd", n, n);
		#endif

		mpz_set(n, original_n);
		factor_list * tmp = factors;
		while(tmp->value != NULL)
		{
			#if VERBOSE
			gmp_printf(" / %Zd", *(tmp->value));
			#endif
			mpz_divexact(n, n, *(tmp->value));
			tmp = tmp->next;
		}
		#if VERBOSE
		gmp_printf(" = %Zd\n\n", n);
		#endif
	}
	#if VERBOSE
	gmp_printf(" :: Letting Pollard's Rho solve %Zd...\n", n);
	#endif

	if (pollard(&factors, n))
	{
		factor_list_print(factors);
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
		printf("# ");
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
