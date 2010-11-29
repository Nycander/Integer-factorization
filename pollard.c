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
	else if (mpz_probab_prime_p(n, 5))
	{
		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, n);
		factor_list_add(f, v);
		return 1;
	}
	// Check for perfect powers
	else if (mpz_perfect_power_p(n))
	{
		// Find the exponent
		mpz_t p, r;
		mpz_init(r);
		mpz_init_set_ui(p, 3);

		while(1)
		{
			if (mpz_root(r, n, mpz_get_ui(p)))
				break;

			mpz_nextprime(p, p);
		}
		// r^p = n
		mpz_t modified_n;
		mpz_init_set(modified_n, n);
		for(int i = 0; i < mpz_get_ui(p); i++)
		{
			mpz_t * v = malloc(sizeof(mpz_t));
			mpz_init_set(*v, r);
			factor_list_add(f, v);
		}
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
	// Check if divisable by 2
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
	/*1-20, 36 valid, none up to 48 is. */
	unsigned int iterations = 0;
	unsigned int x0 = 2;
	mpz_t power;	mpz_init_set_ui(power, 1);
	mpz_t lambda;	mpz_init_set_ui(lambda, 1);
	mpz_t tortoise;	mpz_init_set_ui(tortoise, x0);
	mpz_t hare;		mpz_init(hare); f(hare, tortoise, N);
	mpz_t diff;		mpz_init(diff);
	int ret = 1;

	while(mpz_cmp_ui(divisor,1)==0)
	{
		/*if(iterations++>POLLARD_THRESHOLD)
		{
		#if VERBOSE
			gmp_printf("Gave up on %Zd after %i iterations.\n",N,iterations);
		#endif
			ret = 0;
			break;
		}*/
		if(mpz_cmp(power, lambda)==0)
		{
			mpz_set(tortoise, hare);
			mpz_mul_ui(power, power, 2);
			mpz_set_ui(lambda, 0);
		}
		f(hare, hare, N);
		mpz_add_ui(lambda, lambda, 1);

		mpz_sub(diff, tortoise, hare);
		mpz_abs(diff, diff);
		if(mpz_cmp_ui(diff,0)==0)
		{
			x0++;
			iterations = 0;
			mpz_set_ui(power, 1);
			mpz_set_ui(lambda, 1);
			mpz_set_ui(tortoise, x0);
			mpz_set_ui(divisor,1);
			f(hare, tortoise, N);
			continue;
		}
		mpz_gcd(divisor,diff,N);

	}
	mpz_clear(diff);
	mpz_clear(power);
	mpz_clear(lambda);
	mpz_clear(tortoise);
	mpz_clear(hare);
	return mpz_cmp(divisor, N)==0?0:ret;
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
