#include <gmp.h>

#include "factor.h"
#include "pollard.h"

/**
 * Factor numbers using the Pollard's rho algorithm.
 *
 * @return 0 if failed
 */
int pollard(factor ** f, const mpz_t n)
{
	if (mpz_cmp_ui(n, 1) <= 0)
	{
		return 1;
	}

	// If we REALLY have a prime number
	if (mpz_probab_prime_p(n, 20) == 2)
	{
		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, n);

		factor * nF = malloc(sizeof(factor));
		nF->next = *f;
		nF->value = v;
		*f = nF;
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

	gmp_randstate_t rand_state;
	gmp_randinit_default(rand_state);

	mpz_t divisor;
	mpz_init_set_ui(divisor, 1);

	unsigned int iterations = 0;

	while(mpz_cmp_ui(divisor, 1) == 0 || mpz_cmp(divisor, N) == 0)
	{
		mpz_t x; 		mpz_init(x);
		mpz_urandomm(x, rand_state, N);

		mpz_t y; 		mpz_init(y);
		mpz_urandomm(y, rand_state, N);

		mpz_sub(x, x, y);
		mpz_abs(x, x);

		mpz_clear(y);

		mpz_gcd(divisor, x, N);

		mpz_clear(x);

		if (iterations++ == 1000000)
		{
			mpz_clear(divisor);
			return 0;
		}
	}

	// Great success
	mpz_set(result, divisor);
	mpz_clear(divisor);
	gmp_randclear(rand_state);
	return 1;
}
