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

	int number_count = 0;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		// Let the trial division commence!
		for(unsigned int p = 0; p < good_primes_count; p++)
		{
			if(mpz_divisible_ui_p(numbers[i], good_primes[p]) != 0)
			{
				mpz_divexact_ui(numbers[i], numbers[i], good_primes[p]);
				complete_bit_matrix[p][i] = (complete_bit_matrix[p][i]+1) & 1;

				if (mpz_cmp_ui(numbers[i], 1) == 0)
				{
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
	}
	char bit_matrix[good_primes_count][number_count];
	int n = 0;

	#if VERBOSE
	printf("Matrix will be %d x %d\n", number_count, good_primes_count);
	#endif

	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		if (mpz_cmp_ui(numbers[i], 1) == 0)
		{
			for(unsigned int p = 0; p < good_primes_count; p++)
			{
				bit_matrix[p][n++] = complete_bit_matrix[p][i];
			}
		}
	}

	// Clear our variables!
	mpz_clear(sqrtN), mpz_clear(tmp), mpz_clear(mod);
	for(unsigned int i = 0; i < maxNumberOfSieving; i++)
	{
		mpz_clear(numbers[i]);
		mpz_clear(copy[i]);
	}

	#if VERBOSE
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < number_count; row++)
		{
			printf("%d ", bit_matrix[row][column]);
		}
		printf("\n");
	}
	#endif
	/*
	// TODO: Guass elminiation
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < number_count; row++)
		{
			// Find first 1 in column
			int maxColumn = column;
			for(int r = column+1; bit_matrix[r][row] == 0 && r < good_primes_count; r++)
				maxColumn = r;

			// Swap row i and maxColumn
			char tmp[maxNumberOfSieving];
			for(int c = 0; c < number_count; c++)
			{
				tmp[c] = bit_matrix[column][c];
				bit_matrix[column][c] = bit_matrix[maxColumn][c];
				bit_matrix[maxColumn][c] = tmp[c];
			}

			// Make sure all rows below this row has an initial zero.
			for(int r = column+1; r < good_primes_count; r++)
			{
				if (bit_matrix[r][row] == 0)
					continue;

				// Subtract bit_matrix[k][row] * bit_matrix[column] from bit_matrix[k]
				for(int c = 0; c < number_count; c++)
				{
					bit_matrix[r][c] = bit_matrix[r][c] ^ bit_matrix[column][c];
				}
			}
		}
	}
	for(int column = 0; column < good_primes_count; column++)
	{
		for(int row = 0; row < number_count; row++)
		{
			printf("%d ", bit_matrix[column][row]);
		}
		printf("\n");
	}*/
	/*i := 1
	j := 1
	while (i ≤ m and j ≤ n) do
	  Find pivot in column j, starting in row i:
	  maxi := i
	  for k := i+1 to m do
	    if abs(A[k,j]) > abs(A[maxi,j]) then
	      maxi := k
	    end if
	  end for
	  if A[maxi,j] ≠ 0 then
	    swap rows i and maxi, but do not change the value of i
	    Now A[i,j] will contain the old value of A[maxi,j].
	    divide each entry in row i by A[i,j]
	    Now A[i,j] will have the value 1.
	    for u := i+1 to m do
	      subtract A[u,j] * row i from row u
	      Now A[u,j] will be 0, since A[u,j] - A[i,j] * A[u,j] = A[u,j] - 1 * A[u,j] = 0.
	    end for
	    i := i + 1
	  end if
	  j := j + 1
	end while
	*/
	// TODO: Get solution vectors

	// TODO: Use epic math to calculate factors.

	return 0;
}
