#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#define VERBOSE 1

void factor(const mpz_t n);

gmp_randstate_t rand_state;
int current_input_number = 0;


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

void factor(const mpz_t n){
	//FUN STUFF
}