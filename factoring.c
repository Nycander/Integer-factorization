#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#define VERBOSE 0

void factor(const mpz_t n);
void rho(mpz_t res, const mpz_t N);

gmp_randstate_t rand_state;
int current_input_number = 0;


int main()
{
#if VERBOSE
	printf("Hyper mega global factoring program\n");
	printf("-----------------------------------\n");
	printf("> ");
#endif

	gmp_randinit_default(rand_state);


	mpz_t num;
	mpz_init(num);
	while(gmp_scanf("%Zd", num) > 0)
	{
		current_input_number++;

		factor(num);
		printf("\n");
		mpz_init(num);

	#if VERBOSE
		printf("> ");
	#endif
	}
	mpz_clear(num);
	gmp_randclear(rand_state);
	return 0;
}



/**
 * Factor numbers using the Pollard's rho algorithm.
 */
void factor(const mpz_t n)
{
	if (mpz_cmp_ui(n, 1) <= 0)
	{
		return;
	}

	// If we REALLY have a prime number
	if (mpz_probab_prime_p(n, 20) == 2)
	{
		gmp_printf("%Zd\n", n);
		return;
	}

#if VERBOSE
	gmp_printf("\t%Zd is not a prime, finding factors for it...\n", n);
#endif

	mpz_t divisor;
	mpz_init(divisor);
	rho(divisor, n);

	mpz_t divend;
	mpz_init(divend);
	mpz_divexact(divend, n, divisor);


#if VERBOSE
	gmp_printf("\tFound: %Zd * %Zd = %Zd\n", divisor, divend, n);
#endif

	factor(divisor);
	mpz_clear(divisor);

	factor(divend);
	mpz_clear(divend);

	return;
}

void rho(mpz_t result, const mpz_t N)
{
	if (mpz_even_p(N))
	{
		mpz_set_ui(result, 2);
		return;
	}


	mpz_t c; 		mpz_init(c);
	mpz_urandomm(c, rand_state, N);
	mpz_mod(c, c, N);

	mpz_t x; 		mpz_init(x);
	mpz_urandomm(x, rand_state, N);
	mpz_mod(x, x, N);

	mpz_t y; 		mpz_init(y);
	mpz_t divisor;	mpz_init_set_ui(divisor, 1);

	unsigned int iterations = 0;

	while(mpz_cmp_ui(divisor, 1) == 0)
	{
		mpz_mul(x, x, x);
		mpz_add(x, x, c);
		mpz_mod(x, x, N);

		mpz_mul(y, x, x);
		mpz_add(y, y, c);
		mpz_mod(y, y, N);

		mpz_t xyDiff;
		mpz_init(xyDiff);
		mpz_sub(xyDiff, x, y);
		mpz_abs(xyDiff, xyDiff);

		mpz_gcd(divisor, xyDiff, N);

		if (iterations++ == 2000000)
		{
			exit(current_input_number);
		}

#if VERBOSE
		gmp_printf("\tx: %Zd, y: %Zd, diff: %Zd, gcd(diff, n) = %Zd\n", x, y, xyDiff, divisor);
#endif
	}

	// Fail?
	if (mpz_cmp(divisor, N) == 0)
	{
		rho(divisor, N);
		printf("Rho failed to find a good divisor. Exiting...");
		exit(1);
	}

	// Great success
	mpz_set(result, divisor);
	mpz_clear(divisor);
	mpz_clear(x);
	mpz_clear(y);
}
