#include <stdlib.h>
#include <gmp.h>

#include "trialdivision.h"
#include "factor_list.h"
#include "primes.h"

mpz_t * trial_division(factor_list ** f, int * primes, int primes_count, const mpz_t N)
{

	mpz_t * n = malloc(sizeof(mpz_t));
	mpz_init_set(*n, N);

	// Numbers below 2 should not be factored.
	if (mpz_cmp_ui(N, 1) <= 0)
	{
		return n;
	}

	for(int i = 0; i < primes_count; i++)
	{
		if (mpz_divisible_ui_p(*n, primes[i]))
		{
			mpz_t * prime = malloc(sizeof(mpz_t));
			mpz_init_set_ui(*prime, primes[i]);

			mpz_divexact(*n, *n, *prime);

			factor_list_add(f, prime);

			i--;
		}
	}

	return n;
}
