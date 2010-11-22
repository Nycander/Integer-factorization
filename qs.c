#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "qs.h"

#include "settings.h"
#include "primes.h"

int maxNumberOfSieving = 60;
int smoothnessBound = 90;

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


	// Time to find good numbers! :D

	// Find relevant primes to divide the numbers with
	int good[smoothnessBound];
		good[0] = 2;
	int nrGoodNums = 1;

	for(unsigned int i = 1; i < smoothnessBound; i++)
	{
		mpz_set_ui (mod, primes[i]);
		mpz_powm_ui (tmp, num, (primes[i]-1)/2, mod);

		if(mpz_cmp_ui(tmp, 1) == 0)
		{
			good[nrGoodNums] = primes[i];
			nrGoodNums++;
		}
	}


	// Find the good numbers
	for(unsigned int i = 0; i < nrGoodNums; i++){
		for(unsigned int j = 0; j < maxNumberOfSieving; j++){
			if(mpz_divisible_ui_p(numbers[j], good[i])){
				mpz_divexact_ui(numbers[j], numbers[j], good[i]);
			}
		}
	}

	// Select good ones to return
	factor_list * ret = malloc(sizeof(factor_list));
	ret->value = NULL;
	ret->next = NULL;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++)
	{
		if(mpz_cmp_ui(numbers[i], 1) == 0)
		{
			mpz_t * goodNum = malloc(sizeof(mpz_t));
			mpz_init_set(*goodNum, copy[i]);

			factor_list_add(&ret, goodNum);
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
