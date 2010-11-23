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


	// Time for some gauss!

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


	// Initialize the real bit matrix
	unsigned int bit_matrix_width = number_count;
	unsigned int bit_matrix_height = good_primes_count;
	#if VERBOSE
		printf("\nMatrix will be %d x %d\n", bit_matrix_width, bit_matrix_height);
	#endif
	unsigned char bit_matrix[bit_matrix_height][bit_matrix_width];
	// Copy values from complete_matrix to bit_matrix
	int n = 0;
	for(unsigned int i = 0; i < maxNumberOfSieving; i++) // numbers to factorize
	{
		if (mpz_cmp_ui(numbers[i], 1) == 0 && mpz_cmp_ui(copy[i], 1) != 0)
		{
			for(unsigned int p = 0; p < bit_matrix_height; p++)
			{
				bit_matrix[p][n] = complete_bit_matrix[p][i];
			}
			n++;
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
	for(int column = 0; column < bit_matrix_height; column++)
	{
		for(int row = 0; row < bit_matrix_width; row++)
		{
			printf("%d ", bit_matrix[column][row]);
		}
		printf("\n");
	}
	#endif

	// Gauss elimination
	for(int column = 0, row = 0; column < bit_matrix_width; column++, row++)
	{
		#if VERBOSE
		printf("\n\tLooking at column %d and row %d in matrix:\n", column, row);

		if (bit_matrix_height < 100 && bit_matrix_width < 100)
		{
			for(int column = 0; column < bit_matrix_height; column++)
			{
				printf("\t");
				for(int row = 0; row < bit_matrix_width; row++)
				{
					printf("%d ", bit_matrix[column][row]);
				}
				printf("\n");
			}
			printf("\n");
		}
		#endif

		// Find first 1 in column
		int maxRow = row;
		while(maxRow < bit_matrix_height-1)
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
			for(int c = 0; c < bit_matrix_width; c++)
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
		for(int r = row+1; r < bit_matrix_height; r++)
		{
			if (bit_matrix[r][column] == 0)
				continue;

			#if VERBOSE
			printf("%d ", r);
			#endif

			// Subtract bit_matrix[k][row] * bit_matrix[column] from bit_matrix[k]
			for(int c = 0; c < bit_matrix_width; c++)
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
	for(int column = 0; column < bit_matrix_height; column++)
	{
		for(int row = 0; row < bit_matrix_width; row++)
		{
			printf("%d ", bit_matrix[column][row]);
		}
		printf("\n");
	}
	#endif

	int known = 0;
	int knownIndexes[bit_matrix_width];
	for(int i = 0; i < bit_matrix_width; i++)
		knownIndexes[i] = -1;

	for(int y = 0; y < bit_matrix_height; y++)
	{
		int x = 0;
		while(x < bit_matrix_width)
		{
			if (bit_matrix[y][x] == 1)
			{
				break;
			}
			x++;
		}

		if (x < bit_matrix_width)
		{
			knownIndexes[x] = x;
			known++;
		}
	}

	#if VERBOSE
	printf("\nWe have %d known variables, thus there are %d unknowns\n", known, bit_matrix_width-known);
	#endif

	int unknowns = bit_matrix_width-known;

	// For all 2^unknowns permutations
	for(int i = 0; i < (1 << unknowns); i++)
	{
		// Idea: i can be used with masks to get the current value of the unknowns.
		char solution[bit_matrix_width];

		// Fill the unknowns first
		int unknown_cntr = 0;
		for(int s = 0; s < bit_matrix_width; s++)
		{
			if (knownIndexes[s] != s)
			{
				solution[s] = (i & (1 << unknown_cntr)) >> unknown_cntr;
				unknown_cntr++;
			}
			else
			{
				solution[s] = -1;
			}
		}

		// Fill the knowns bottom up
		for(int s = bit_matrix_width-1; s >= 0; s--)
		{
			if (solution[s] != -1)
				continue;

			solution[s] = 0;

			// Search the column for the "1" (should be only 1)
			int row = -1;
			for(int y = 0; y < bit_matrix_height; y++)
			{
				if (bit_matrix[y][s] == 1)
				{
					row = y;
					break;
				}
			}

			// Calculate value from equation vector
			for(int x = s; x < bit_matrix_width; x++)
			{
				solution[s] = solution[s] ^ bit_matrix[row][x];
			}
		}

		#if VERBOSE
		printf("Solution array: ");
		for(int s = 0; s < bit_matrix_width; s++)
		{
			printf("%d ", solution[s]);
		}
		printf("\n");
		#endif
	}


	return 0;
}
