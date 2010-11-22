
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "primes.h"
#include "factor_list.h"

#define VERBOSE 1

void quadraticSieve(const mpz_t num);
factor_list * sieving(const mpz_t num);

gmp_randstate_t rand_state;
int current_input_number = 0;
int maxNumberOfSieving = 60;
int smoothnessBound = 90;



void quadraticSieve(const mpz_t num){
	//FUN STUFF
}

/*
* LOLZ fungerar det här e jag kung //HIRSCHEN
*/
factor_list * sieving(const mpz_t num){
	mpz_t sqrtN, tmp, mod;
	mpz_inits(sqrtN, tmp, mod);
	mpz_sqrt(sqrtN, num);


	mpz_t numbers[maxNumberOfSieving], copy[maxNumberOfSieving];
	//Init
	for(unsigned int i = 0; i < maxNumberOfSieving; i++){
		mpz_inits(numbers[i], copy[i]);
	}

	//Generate numbers
	for(unsigned int i = 0; i < maxNumberOfSieving; i++, mpz_add_ui(sqrtN, sqrtN, 1){
		mpz_mul_si(tmp, sqrtN, 2);
		mpz_sub(numbers[i],tmp,num);
		mpz_sub(copy[i],tmp,num);
	}


	//Compute to find good numbers
	//Find relevant primes to div
	int good[smoothnessBound], nrGoodNums = 1;
	int good[0] = 2;
	for(unsigned int i = 1; i < smoothnessBound; i++){
		mpz_set_ui (mod, primes[i]);
		mpz_powm_ui (tmp, num, (primes[i]-1)/2, mod);
		if(tmp == 1){
			good[nrGoodNums] = primes[i];
			nrGoodNums++;
		}
	}

	//Find the good numbers
	for(unsigned int i = 0; i < nrGoodNums; i++){
		for(unsigned int j = 0; j < maxNumberOfSieving; j++){
			if(mpz_divisible_ui_p(numbers[j], good[i]) != 0){
				mpz_divexact_ui(numbers[j], numbers[j], good[i]);
			}
		}
	}

	//Select good ones to return
	factor_list * ret = malloc(sizeof(factor_list));
	ret->value = NULL;
	ret->next = NULL;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++){
		if(numbers[i] == 1){
			mpz_t * goodNum = malloc(sizeof(mpz_t));
			mpz_init_set(*goodNum, copy[i]);

			factor_list_add(&ret, goodNum);
		}
	}

	//Clear mpz
	mpz_clears(sqrtN, tmp, mod);
	for(unsigned int i = 0; i < maxNumberOfSieving; i++){
		mpz_clears(numbers[i], copy[i]);
	}

	return ret;
}
