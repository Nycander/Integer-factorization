#include <stdlib.h>
#include <stdio.h>
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
	else if (mpz_probab_prime_p(n, 10))
	{
		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, n);
		factor_list_add(f, v);
		return 1;
	}


#if VERBOSE
	gmp_printf("\tSearching for x * y = %Zd ...\n", n);
#endif

	mpz_t divisor;		mpz_init(divisor);
	mpz_t divend;		mpz_init(divend);
	// Check for even square root
	if (mpz_perfect_square_p(n))
	{
		mpz_sqrt(divisor, n);
		mpz_set(divend, divisor);
	}
	else
	{
		if (! rho(divisor, n))
			return 0;

		mpz_divexact(divend, n, divisor);
	}

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

#if USE_BRENT_INSTEAD_OF_FLOYD
	if(! brent(N, divisor))
#else
	if(! floyd(N, divisor))
#endif
	{
		mpz_clear(divisor);
		return 0;
	}

	// Great success
	mpz_set(result, divisor);
	mpz_clear(divisor);
	return 1;
}

int brent(const mpz_t N, mpz_t divisor)
{
	mpz_t x, y, r, q, k, ys, g, m, tmp, i;
	mpz_inits(x, y, r, q, k, ys, g, m, tmp, i);
	mpz_set(y, x0); mpz_set_ui(r,1); mpz_set_uit(q,1);
	do
	{
		mpz_set(x,y);
		for(mpz_set_ui(i,0);mpz_cmp(r,i)>0;mpz_add_ui(i,i,1))
		{
			f(y,y,N);
		}
	}
	mpz_clears(x, y, r, q, k, ys, g, m, tmp, i);
}

int floyd(const mpz_t N, mpz_t divisor)
{
	mpz_t x;		mpz_init_set_ui(x, 2);
	mpz_t y;		mpz_init_set_ui(y, 2);

	unsigned long iterations = 0;
	while(mpz_cmp_ui(divisor, 1) == 0)
	{
		if (++iterations == POLLARD_THRESHOLD)
		{
		#if VERBOSE
			gmp_printf("\tGave up after %d iterations on number: %Zd\n", iterations-1, N);
		#endif
			mpz_clear(x);
			mpz_clear(y);
			return 0;
		}

		#if VERBOSE
		gmp_printf("\tx = f(%Zd)", x);
		#endif

		f(x, x, N); // x = f(x)

		#if VERBOSE
		gmp_printf(" = %Zd,\t", x);
		gmp_printf("y = f(f(%Zd))", y);
		#endif

		f(y, y, N); // y = f(x)
		f(y, y, N); // y = f(x)
		#if VERBOSE
		gmp_printf(" = %Zd", y);
		#endif

		mpz_t diff;	mpz_init(diff);
		mpz_sub(diff, x, y);
		mpz_abs(diff, diff);
		#if VERBOSE
		gmp_printf(",\t|x-y| = %Zd", diff);
		#endif

		if (mpz_cmp_ui(diff, 0) == 0)
		{
		#if VERBOSE
			gmp_printf("...\n\tVisited all numbers after %d iterations on number: %Zd\n", iterations-1, N);
		#endif
			mpz_clear(x);
			mpz_clear(y);
			return 0;
		}

		mpz_gcd(divisor, diff, N);
		mpz_clear(diff);
		#if VERBOSE
		gmp_printf(",\tdivisor = %Zd\n", divisor);
		#endif
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
