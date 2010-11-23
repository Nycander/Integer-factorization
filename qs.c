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

		mpz_init_set(copy[i], numbers[i]);
	}


	// Time to find good prime numbers! :D

	// Find relevant primes to divide the numbers with
	int good_primes[smoothnessBound];
		good_primes[0] = 2;
	int good_primes_count = 1;

	#if VERBOSE
	printf("Finding good primes: \n\t2\n");
	#endif

	for(unsigned int i = 1; mpz_cmp_ui(num, primes[i]) > 0 && i < smoothnessBound; i++)
	{
		mpz_set_ui (mod, primes[i]);
		mpz_powm_ui (tmp, num, (primes[i]-1)/2, mod);

		if(mpz_cmp_ui(tmp, 1) == 0)
		{
			good_primes[good_primes_count] = primes[i];
			good_primes_count++;

			#if VERBOSE
			printf("\t%d\n", primes[i]);
			#endif
		}
	}

	// Initialize bit matrix
	char complete_bit_matrix[good_primes_count][maxNumberOfSieving];
	for(int i = 0; i < good_primes_count; i++)
	{
		for(int j = 0; j < maxNumberOfSieving; j++)
		{
			complete_bit_matrix[i][j] = 0;
		}
	}

	// Find the good prime numbers
	factor_list * ret = malloc(sizeof(factor_list));
	ret->value = NULL;
	ret->next = NULL;

	#if VERBOSE
	printf("Computing good numbers: \n");
	#endif
	int number_count = 0;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		#if VERBOSE
		gmp_printf("\t%Zd = ", numbers[i]);
		#endif

		// Let the trial division commence!
		for(unsigned int p = 0; p < good_primes_count; p++)
		{
			if(mpz_divisible_ui_p(numbers[i], good_primes[p]) != 0)
			{
				#if VERBOSE
				printf("%d * ", good_primes[p]);
				#endif

				mpz_divexact_ui(numbers[i], numbers[i], good_primes[p]);
				complete_bit_matrix[p][i] = (complete_bit_matrix[p][i]+1) & (char)1;

				if (mpz_cmp_ui(numbers[i], 1) == 0)
				{
					#if VERBOSE
					gmp_printf("1 = OK!");
					#endif
					mpz_t * goodPrimeNumber = malloc(sizeof(mpz_t));
					mpz_init_set(*goodPrimeNumber, copy[i]);

					factor_list_add(&ret, goodPrimeNumber);
					number_count++;

					break;
				}
				else
				{
					--p;
				}
			}
		}
		#if VERBOSE
		printf("\n");
		#endif
	}

	#if VERBOSE
	printf("\nComplete bit matrix:\n\n");
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < maxNumberOfSieving; row++)
		{
			printf("%d ", complete_bit_matrix[column][row]);
		}
		printf("\n");
	}
	#endif

	#if VERBOSE
	printf("\nMatrix will be %d x %d\n", number_count, good_primes_count);
	#endif

	// Initialize the real bit matrix
	char bit_matrix[good_primes_count][number_count];
	// Copy values from complete_matrix to bit_matrix
	int n = 0;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		if (mpz_cmp_ui(numbers[i], 1) == 0 && mpz_cmp_ui(copy[i], 1) != 0)
		{
			for(unsigned int p = 0; p < good_primes_count; p++)
			{
				bit_matrix[p][n] = complete_bit_matrix[p][i];
			}
			n++;
		}
	}

	#if VERBOSE
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < number_count; row++)
		{
			printf("%d ", bit_matrix[column][row]);
		}
		printf("\n");
	}
	#endif

	// Clear our variables!
	mpz_clear(sqrtN), mpz_clear(tmp), mpz_clear(mod);
	for(unsigned int i = 0; i < maxNumberOfSieving; i++)
	{
		mpz_clear(numbers[i]);
		mpz_clear(copy[i]);
	}

	// Guass elminiation
	for(int column = 0, row = 0; column < number_count; column++, row++)
	{
		#if VERBOSE
		printf("\n\tLooking at column %d amd row %d in matrix:\n", column, row);
		for(int row = 0; row < good_primes_count; row++)
		{
			printf("\t");
			for(int column = 0; column < number_count; column++)
			{
				printf("%d ", bit_matrix[row][column]);
			}
			printf("\n");
		}
		printf("\n");
		#endif

		// Find first 1 in column
		int maxRow = row;
		while(maxRow < good_primes_count-1)
		{
			if (bit_matrix[maxRow][column] == 1)
				break;

			maxRow++;
		}

		// If we couldn't find a 1
		if (bit_matrix[maxRow][column] == 0)
		{
			row--;
			continue;
		}

		#if VERBOSE
		printf("\tFound a 1 on row %d\n", maxRow);
		#endif

		// If we must replace the largest row to the top, swap them.
		if (maxRow != row)
		{
			#if VERBOSE
			printf("\tSwapping rows %d and %d...\n", row, maxRow);
			#endif

			// Swap row i and maxColumn
			char tmp[maxNumberOfSieving];
			for(int c = 0; c < number_count; c++)
			{
				tmp[c] = bit_matrix[row][c];
				bit_matrix[row][c] = bit_matrix[maxRow][c];
				bit_matrix[maxRow][c] = tmp[c];
			}
		}

		#if VERBOSE
		printf("\tXOR-ing row %d with rows... ", row);
		#endif

		// Make sure all rows below this row has an initial zero.
		for(int r = row+1; r < good_primes_count; r++)
		{
			if (bit_matrix[r][column] == 0)
				continue;

			#if VERBOSE
			printf("%d ", r);
			#endif

			// Subtract bit_matrix[k][row] * bit_matrix[column] from bit_matrix[k]
			for(int c = 0; c < number_count; c++)
			{
				bit_matrix[r][c] = bit_matrix[r][c] ^ bit_matrix[row][c];
			}
		}
		#if VERBOSE
		printf("\n");
		#endif
	}
	#if VERBOSE
	printf("\nAfter gauss:\n\n");
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < number_count; row++)
		{
			printf("%d ", bit_matrix[column][row]);
		}
		printf("\n");
	}
	#endif

	// TODO: Get solution vectors

	// TODO: Use epic math to calculate factors.

	return 0;
}
