#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "qs.h"

#include "pollard.h"
#include "settings.h"
#include "primes.h"

int smoothnessBound = 500;

// TEH EPIC QUADRATIC SIEVE!

int quadratic_sieve(factor_list ** result, const mpz_t num)
{

	mpz_t number_result;
	mpz_init_set(number_result, num);

	while(mpz_cmp_ui(number_result,1) != 0){

		#if VERBOSE
		gmp_printf("\nFactoring the number %Zd ...\n\n", number_result);
		#endif

		mpz_t sqrtN, tmp, mod,ret1, ret2;
		mpz_init(sqrtN), mpz_init(tmp), mpz_init(ret1), mpz_init(ret2), mpz_init(mod);
		mpz_sqrt(sqrtN, num);


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
		int num_sieved_count = good_primes_count+1;


		#if VERBOSE
		printf("Initializing bit matrix...\n");
		#endif

		// Initialize bit matrix
		char bit_matrix[good_primes_count][num_sieved_count];
		for(int i = 0; i < good_primes_count; i++)
		{
			for(int j = 0; j < num_sieved_count; j++)
			{
				bit_matrix[i][j] = 0;
			}
		}


		// Find the good prime numbers
		mpz_t nums[num_sieved_count];
		int number_count = 0;
		int OGindexes[num_sieved_count];

		#if VERBOSE
		gmp_printf("Sieving numbers starting at %Zd...\n", sqrtN);
		#endif

		// Generate numbers
		for(unsigned int i = 0; num_sieved_count > number_count; i++)
		{
			mpz_t smoothNumCand;
			mpz_init(smoothNumCand);

			mpz_add_ui(sqrtN, sqrtN, 1);
			mpz_mul(tmp, sqrtN, sqrtN);

			mpz_sub(smoothNumCand,tmp,num);
			mpz_set(tmp, smoothNumCand);

			// Let the trial division commence!
			for(unsigned int p = 0; p < good_primes_count; p++)
			{
				if(mpz_divisible_ui_p(tmp, good_primes[p]) != 0)
				{
					mpz_divexact_ui(tmp, tmp, good_primes[p]);

					if (mpz_cmp_ui(tmp, 1) == 0)
					{

						bit_matrix[p][number_count] = (bit_matrix[p][number_count]+1) & (char)1;
						mpz_init_set(nums[number_count], smoothNumCand);
						OGindexes[number_count++] = i;
						#if VERBOSE
						gmp_printf("Storing: %Zd \n", nums[number_count-1]);
						#endif
						break;
					}
					else
					{
						--p;
					}
				}
			}
			mpz_clear(smoothNumCand);
		}

		#if VERBOSE
		printf("\nFound the following numbers: ");
		for(int i = 0; i < number_count; i++)
		{
			gmp_printf("%Zd ", nums[i]);
		}
		printf("\n");
		#endif

		// Time for some gauss!

		// If we have an overdetermined matrix, we must fail.
		if (good_primes_count > number_count)
		{
			#if VERBOSE
			printf("\n!!!! We have an overdetermined matrix with %d equations and %d unknowns !!! \n", good_primes_count, number_count);
			#endif
			return 0;
		}

		// Initialize the real bit matrix
		unsigned int bit_matrix_width = number_count;
		unsigned int bit_matrix_height = good_primes_count;

		#if VERBOSE
		printf("\nMatrix is of size %d x %d\n", bit_matrix_width, bit_matrix_height);

		if (bit_matrix_height < 100 && bit_matrix_width < 100)
		{
			for(int column = 0; column < bit_matrix_height; column++)
			{
				for(int row = 0; row < bit_matrix_width; row++)
				{
					printf("%d ", bit_matrix[column][row]);
				}
				printf("\n");
			}
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
				char tmp[num_sieved_count];
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
		if (bit_matrix_height < 100 && bit_matrix_width < 100)
		{
			printf("\nAfter gauss:\n\n");
			for(int column = 0; column < bit_matrix_height; column++)
			{
				for(int row = 0; row < bit_matrix_width; row++)
				{
					printf("%d ", bit_matrix[column][row]);
				}
				printf("\n");
			}
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

			//Orginaltalens produkt
			mpz_set_ui(ret1,1);
			for(int s = 0; s < bit_matrix_width; s++)
			{
				if (solution[s] == 0)
					continue;

				mpz_mul(ret1, ret1, nums[s]);
			}
			mpz_sqrt(ret1, ret1);

			//Orginalfaktorernas produkt
			mpz_set_ui(tmp,1);
			mpz_set_ui(ret2,1);
			mpz_sqrt(sqrtN, num);
			for(int s = 0; s < bit_matrix_width; s++)
			{
				if (solution[s] == 0)
					continue;

				mpz_add_ui(tmp, sqrtN, (OGindexes[s]+1));
				mpz_mul(ret2, ret2, tmp);
			}

			//tmp save
			mpz_set(tmp, ret1);
			//num1
			mpz_add(ret1, ret2, ret1);
			//num2
			mpz_sub(ret2, ret2, tmp);

			//factor 1
			mpz_gcd(ret1, ret1, num);
			//factor 2
			mpz_gcd(ret2, ret2, num);

			// Try to store the factors
			if(try_adding_factor_to_result(ret1, number_result, result)){
				mpz_divexact(number_result, number_result, ret1);
			}
			if(try_adding_factor_to_result(ret2, number_result, result)){
				mpz_divexact(number_result, number_result, ret2);
			}
		}

		// Clear our variables!
		mpz_clear(sqrtN), mpz_clear(ret1), mpz_clear(ret2), mpz_clear(tmp), mpz_clear(mod);
	}

	return 1;
}

int try_adding_factor_to_result(mpz_t factor, const mpz_t ofNumber, factor_list ** result)
{
	if (mpz_cmp_ui(factor, 1) == 0)
	{
		return 0;
	}
	if (mpz_cmp(factor, ofNumber) == 0)
	{
		return 0;
	}

	#if VERBOSE
	gmp_printf("Storing factor: %Zd\n", factor);
	#endif

	pollard(result, factor);

	return 1;
	/*mpz_t * m = malloc(sizeof(mpz_t));
	mpz_init_set(*m, factor);
	factor_list_add(result, m);*/
}
