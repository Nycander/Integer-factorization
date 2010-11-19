#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "primes.h"
#include "factor_list.h"

#define VERBOSE 1

void factor(const mpz_t num);
factor_list * sieving(const mpz_t num);

gmp_randstate_t rand_state;
int current_input_number = 0;
int maxNumberOfSieving = 60;
int smoothnessBound = 90;


int main(int argc, char * argv[])
{
#if VERBOSE
	printf("Hyper mega global factoring program\n");
	printf("-----------------------------------\n");
	printf("> ");
#endif

	gmp_randinit_default(rand_state);

	int limit = (argc == 2 ? atoi(argv[1]) : INT_MAX );

	mpz_t threshold;
	mpz_init_set_str(threshold, "83003906250000000", 10); // Brute-forced value @ kattis

	mpz_t num;
	while(++current_input_number <= limit)
	{
		mpz_init(num);
		if (gmp_scanf("%Zd", num) <= 0)
		{
			mpz_clear(num);
			break;
		}

		if (mpz_cmp(num, threshold) <= 0)
		{
			factor(num);
			printf("\n");
		}
		else
		{
			printf("fail\n\n");
		}

		mpz_clear(num);

	#if VERBOSE
		printf("> ");
	#endif
	}
	gmp_randclear(rand_state);
	return 0;
}

void factor(const mpz_t num){
	//FUN STUFF
}

/*
* LOLZ fungerar det här e jag kung
*/
factor_list * sieving(const mpz_t num){
	mpz_t sqrtN, tmp, mod;
	mpz_inits(sqrtN, tmp, mod);
	mpz_sqrt(sqrtN, num);


	mpz_t numbers[maxNumberOfSieving], copy[maxNumberOfSieving];
	//Init
	for(int i = 0; i < maxNumberOfSieving; i++){
		mpz_inits(numbers[i], copy[i]);
	}

	//Generate numbers
	for(int i = 0; i < maxNumberOfSieving; i++, mpz_add_ui(sqrtN, sqrtN, 1){
		mpz_mul_si(tmp, sqrtN, 2);
		mpz_sub(numbers[i],tmp,num);
		mpz_sub(copy[i],tmp,num);
	}


	//Compute to find good numbers
	//Find relevant primes to div
	int good[smoothnessBound], nrGoodNums = 1;
	int good[0] = 2;
	for(int i = 1; i < smoothnessBound; i++){
		mpz_set_ui (mod, primes[i]);
		mpz_powm_ui (tmp, num, (primes[i]-1)/2, mod);
		if(tmp == 1){
			good[nrGoodNums] = primes[i];
			nrGoodNums++;
		}
	}

	//Find the good numbers
	for(int i = 0; i < nrGoodNums; i++){
		for(int j = 0; j < maxNumberOfSieving; j++){
			if(mpz_divisible_ui_p(numbers[j], good[i]) != 0){
				mpz_divexact_ui(numbers[j], numbers[j], good[i]);
			}
		}
	}

	//Select good ones to return
	factor_list * ret = malloc(sizeof(factor_list));
	ret->value = NULL;
	ret->next = NULL;
	for(int i = 0; i < maxNumberOfSieving; i++){
		if(numbers[i] == 1){
			mpz_t * goodNum = malloc(sizeof(mpz_t));
			mpz_init_set(*goodNum, copy[i]);

			factor_list_add(&ret, goodNum);
		}
	}

	//Clear mpz
	mpz_clears(sqrtN, tmp, mod);
	for(int i = 0; i < maxNumberOfSieving; i++){
		mpz_clears(numbers[i], copy[i]);
	}

	return ret;
}

