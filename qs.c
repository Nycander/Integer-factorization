#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "qs.h"

#include "settings.h"
#include "primes.h"

int maxNumberOfSieving = 60;
int smoothnessBound = 500;

int quadratic_sieve(factor_list ** result, const mpz_t num)
{
	factor_list * list = sieving(num);

	// TODO: Factor the numbers

	// TODO: For each factor: % mod 1  and put in array

	// TODO: The array assembles a matrix
	return 0;
}

/*
* LOLZ fungerar det här e jag kung //HIRSCHEN
*/
factor_list * sieving(const mpz_t num){
	mpz_t sqrtN, tmp, mod;
	mpz_init(sqrtN), mpz_init(tmp), mpz_init(mod);
	mpz_sqrt(sqrtN, num);


	mpz_t numbers[maxNumberOfSieving];
	mpz_t copy[maxNumberOfSieving];


	// Generate numbers
	for(unsigned int i = 0; i < maxNumberOfSieving; i++)
	{
		mpz_add_ui(sqrtN, sqrtN, 1);
		mpz_mul(tmp, sqrtN, sqrtN);

		mpz_init(numbers[i]);
		mpz_sub(numbers[i],tmp,num);

		mpz_init(copy[i]);
		mpz_sub(copy[i],tmp,num);
	}


	// Time to find good prime numbers! :D

	// Find relevant primes to divide the numbers with
	int good_primes[smoothnessBound];
		good_primes[0] = 2;
	int good_primes_count = 1;

	for(unsigned int i = 1; i < smoothnessBound; i++)
	{
		mpz_set_ui (mod, primes[i]);
		mpz_powm_ui (tmp, num, (primes[i]-1)/2, mod);

		if(mpz_cmp_ui(tmp, 1) == 0)
		{
			good_primes[good_primes_count] = primes[i];
			good_primes_count++;
		}
	}

	// Initialize bit matrix
	char bit_matrix[maxNumberOfSieving][maxNumberOfSieving];
	for(int i = 0; i < maxNumberOfSieving; i++)
	{
		for(int j = 0; j < maxNumberOfSieving; j++)
		{
			bit_matrix[i][j] = 0;
		}
	}

	// Find the good prime numbers
	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		for(unsigned int p = 0; p < good_primes_count; p++) // wtf is i?
		{
			if(mpz_divisible_ui_p(numbers[i], good_primes[p]) != 0)
			{
				mpz_divexact_ui(numbers[i], numbers[i], good_primes[p]);
				bit_matrix[i][p] = (bit_matrix[i][p]+1) & 1;

				if (mpz_cmp_ui(numbers[i], 1) == 0)
				{
					break;
				}
				else
				{
					--p;
				}
			}
		}
	}

	// Select good number for return
	factor_list * ret = malloc(sizeof(factor_list));
	ret->value = NULL;
	ret->next = NULL;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++)
	{
		if(mpz_cmp_ui(numbers[i], 1) == 0)
		{
			mpz_t * goodPrimeNumber = malloc(sizeof(mpz_t));
			mpz_init_set(*goodPrimeNumber, copy[i]);

			factor_list_add(&ret, goodPrimeNumber);
		}
	}

	//Clear mpz
	mpz_clear(sqrtN), mpz_clear(tmp), mpz_clear(mod);
	for(unsigned int i = 0; i < maxNumberOfSieving; i++){
		mpz_clear(numbers[i]);
		mpz_clear(copy[i]);
	}

	return ret;
}
