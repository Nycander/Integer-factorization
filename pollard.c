#include <stdlib.h>
#include <gmp.h>

#include "settings.h"
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

	mpz_t divisor;	mpz_init_set_ui(divisor, 1);

	if(! floyd(N, divisor))
		return 0;

	// Great success
	mpz_set(result, divisor);
	mpz_clear(divisor);
	return 1;
}

int brent(const mpz_t N, mpz_t divisor)
{
	mpz_t x;		mpz_init_set_ui(x, 2);

	mpz_t tortoise, hare; mpz_init(tortoise); mpz_init(hare);
	unsigned int power = 1, lam = 1;
	mpz_set(tortoise, x);
	f(hare, tortoise, N);
	while(mpz_cmp(tortoise, hare) != 0)
	{
		if(power == lam)
		{
			mpz_set(tortoise, hare);
			power *= 2;
			lam = 0;
		}
		f(hare, hare, divisor);
		lam++;
	}
	unsigned int mu = 0;
	mpz_set(tortoise, x); mpz_set(hare, x);
	int i;
	for(i=0; i<lam;i++)
	{
		f(hare, hare, N);
	}
	while(mpz_cmp(tortoise, hare) != 0)
	{
		f(tortoise, tortoise, N);
		f(hare, hare, N);
		mu++;
	}
	mpz_sub(tortoise, tortoise, hare);
	mpz_abs(tortoise, tortoise);
	mpz_gcd(divisor, tortoise, N);
	mpz_clear(x);
	return 1;
}

int floyd(const mpz_t N, mpz_t divisor)
{
	mpz_t x;		mpz_init_set_ui(x, 2);
	mpz_t y;		mpz_init_set_ui(y, 2);

	unsigned int iterations = 0;

	while(mpz_cmp_ui(divisor, 1) == 0)
	{
		if (iterations++ == POLLARD_THRESHOLD)
		{
			#if VERBOSE
			gmp_printf("\tGave up after %d iterations on number: %Zd\n", iterations-1, N);
			#endif
			mpz_clear(divisor);
			return 0;
		}
		
		f(x, x, N);
		f(y, x, N);

	#if VERBOSE
		gmp_printf("x = %Zd, y = %Zd ...  ", x, y);
	#endif

		mpz_sub(x, x, y);
		mpz_abs(x, x);

		if (mpz_cmp_ui(x, 0) == 0)
		{
			continue;
		}

		mpz_gcd(divisor, x, N);
	}

	mpz_clear(x);
	mpz_clear(y);
	return 1;
}

void f(mpz_t result, mpz_t x, const mpz_t N)
{
	mpz_set(result, x);
	mpz_mul(result, result, result);
	mpz_add_ui(result, result, 1);
	mpz_mod(result, result, N);
}
