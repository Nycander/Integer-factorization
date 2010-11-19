#include <stdlib.h>
#include <gmp.h>

#include "factor_list.h"
#include "pollard.h"

/**
 * Factor numbers using the Pollard's rho algorithm.
 *
 * @return 0 if failed
 */
int pollard(factor_list ** f, const mpz_t n)
{
	// Numbers below 2 should not be factored.
	if (mpz_cmp_ui(n, 1) <= 0)
	{
		return 1;
	}

	// Base case: we have a prime number
	if (mpz_probab_prime_p(n, 10))
	{
		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, n);
		factor_list_add(f, v);
		return 1;
	}

#if VERBOSE
	gmp_printf("\tSearching for x * y = %Zd ...\n", n);
#endif

	mpz_t divisor;
	mpz_init(divisor);

	if (! rho(divisor, n))
		return 0;

	mpz_t divend;
	mpz_init(divend);
	mpz_divexact(divend, n, divisor);


#if VERBOSE
	gmp_printf("\tFound: %Zd * %Zd = %Zd\n", divisor, divend, n);
#endif

	int r = pollard(f, divisor);
	mpz_clear(divisor);

	if (! r)
		return 0;

	r = pollard(f, divend);
	mpz_clear(divend);

	if (! r)
		return 0;

	return 1;
}

int rho(mpz_t result, const mpz_t N)
{
	// Check if divided by 2
	if (mpz_even_p(N))
	{
		mpz_set_ui(result, 2);
		return 1;
	}

	mpz_t x;		mpz_init_set_ui(x, 2);
	mpz_t y;		mpz_init_set_ui(y, 2);
	mpz_t divisor;	mpz_init_set_ui(divisor, 1);

	unsigned int iterations = 0;

	while(mpz_cmp_ui(divisor, 1) == 0)
	{
		f(x, x, N);
		f(y, x, N);

		mpz_sub(x, x, y);
		mpz_abs(x, x);

		mpz_gcd(divisor, x, N);

		if (iterations++ == 100000)
		{
			mpz_clear(divisor);
			return 0;
		}
	}

	// Great success
	mpz_set(result, divisor);
	mpz_clear(divisor);
	mpz_clear(y);
	mpz_clear(x);
	return 1;
}

void f(mpz_t result, mpz_t x, const mpz_t N)
{
	mpz_clear(result);
	mpz_init_set(result, x);
	mpz_mul(result, result, result);
	mpz_add_ui(result, result, 1);
	mpz_mod(result, result, N);
}