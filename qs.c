#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "qs.h"
#include "factor_list.h"
#include "settings.h"
#include "primes.h"

int maxNumberOfSieving = 60;
int smoothnessBound = 90;

int quadratic_sieve(factor_list ** result, const mpz_t num)
{
	/*factor_list * list =*/ sieving(num);
	// TODO: Factor numbers % mod 1 in matrix
	return 0;
}

/*
* LOLZ fungerar det här e jag kung //HIRSCHEN
*/
factor_list * sieving(const mpz_t num){
	mpz_t sqrtN, tmp, mod;
	mpz_init(sqrtN), mpz_init(tmp), mpz_init(mod);
	mpz_sqrt(sqrtN, num);


	mpz_t numbers[maxNumberOfSieving], copy[maxNumberOfSieving];
	// Init
	for(unsigned int i = 0; i < maxNumberOfSieving; i++){
		mpz_init(numbers[i]);
		mpz_init(copy[i]);
	}

	// Generate numbers
	for(unsigned int i = 0; i < maxNumberOfSieving; i++, mpz_add_ui(sqrtN, sqrtN, 1))
	{
		mpz_mul_si(tmp, sqrtN, 2);
		mpz_sub(numbers[i],tmp,num);
		mpz_sub(copy[i],tmp,num);
	}


	// Compute to find good numbers
	// Find relevant primes to div
	int good[smoothnessBound];
	int nrGoodNums = 1;
	good[0] = 2;
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
			if(mpz_divisible_ui_p(numbers[j], good[i]) != 0){
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
