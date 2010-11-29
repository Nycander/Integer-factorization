#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>
#include <math.h>

#include "qs.h"

#include "pollard.h"
#include "settings.h"
#include "primes.h"
#include "shanks.h"

#define GOODPRIME_VERBOSE 0
#define SIEVE_VERBOSE 0
#define MATRIX_VERBOSE 0
#define SOLUTION_ARRAY_VERBOSE 0

int smoothness_bound = 500;

// TEH EPIC QUADRATIC SIEVE!

int quadratic_sieve(factor_list ** result, const mpz_t num)
{
	/**/
	int num_size = mpz_sizeinbase(num, 2);
	int ln_n = M_LN2 * (double)num_size;

	//smoothness_bound = (int) (0.63*pow(exp(sqrt(ln_n * log(ln_n))), 0.35355339059));
	long d = 1.5;
	smoothness_bound = (int) (0.63*pow(ln_n, d*sqrt(ln_n*log(ln_n))));


	#if VERBOSE
	printf(" :: Smoothness-bound = %d \n", smoothness_bound);
	#endif
	/**/
	//smoothness_bound = mpz_sizeinbase(num, 2);


	// Numbers below 2 should not be factored.
	if (mpz_cmp_ui(num, 1) <= 0)
	{
		return 1;
	}
	// Base case: we have a prime number
	else if (mpz_probab_prime_p(num, 10))
	{
		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, num);
		factor_list_add(result, v);
		return 1;
	}

	#if VERBOSE
	gmp_printf(" :: Factoring the number %Zd using QS: \n \n ", num);
	#endif

	mpz_t tmp;
	mpz_init(tmp);

	// Time to find good prime numbers! :D

	// Find the good prime numbers
	mpz_t nums[2*smoothness_bound];
	mpz_t nums_copy[2*smoothness_bound];
	mpz_t nums_p[2*smoothness_bound];

	#if VERBOSE
	gmp_printf("Finding %d numbers which satisfies the relation a^2 = %Zd (mod good_prime)\n", 2*smoothness_bound, num);
	#endif

	// Generate numbers
	int primes = 1;
	int num_index = 0;
	mpz_t prime;
	mpz_init_set_ui(prime, 2);
	while(num_index < 2*smoothness_bound)
	{
		mpz_nextprime(prime, prime); // side-effect: start at 3

		mpz_t * smooth_candidate = shanks_tonelli(num, prime);

		mpz_t r, R;
		mpz_init_set(R, *smooth_candidate);
		mpz_init(r);
		mpz_sub(r, prime, R);

		mpz_clear(*smooth_candidate);
		free(smooth_candidate);

		if (mpz_cmp(R, prime) == 0)
		{
			#if VERBOSE && SIEVE_VERBOSE
			gmp_printf(" %Zd ^2 = %Zd (mod %Zd): Useless, could not run Tonelli-Shanks?\n", R, num, prime);
			#endif
		}
		else
		{
			if (mpz_cmp_ui(R, 1) > 0)
			{
				#if VERBOSE && SIEVE_VERBOSE
				gmp_printf(" %Zd ^2 = %Zd (mod %Zd) : OK!\n", R, num, prime);
				#endif
				mpz_init_set(nums[num_index], R);
				mpz_init_set(nums_copy[num_index], R);
				mpz_init_set(nums_p[num_index], prime);
				num_index++;
			}

			if (mpz_cmp_ui(R, 1) > 0 && mpz_cmp(R, r) != 0)
			{
				#if VERBOSE && SIEVE_VERBOSE
				gmp_printf(" %Zd = %Zd - %Zd : OK!\n", r, R, prime);
				#endif
				mpz_init_set(nums[num_index], r);
				mpz_init_set(nums_copy[num_index], r);
				mpz_init_set(nums_p[num_index], prime);
				num_index++;
			}
		}
		mpz_clear(R);
		mpz_clear(r);
		primes++;
	}
	mpz_clear(prime);

	int number_count = num_index;

	#if VERBOSE
	printf("\n Visited %d prime numbers.\n", primes);
	printf("\n Found the following %d numbers: ", number_count);
	for(int i = 0; i < number_count; i++)
	{
		gmp_printf("%Zd ", nums[i]);
	}
	printf("\n ");

	#if MATRIX_VERBOSE
	printf("Initializing bit matrix...\n ");
	#endif
	#endif
	// Initialize bit matrix
	unsigned int bit_matrix_width = number_count;
	unsigned int bit_matrix_height = primes;


	char bit_matrix[bit_matrix_height][bit_matrix_width];

	// Populate matrix with trial division
	mpz_init_set_ui(prime, 2);
	mpz_t mod; mpz_init(mod);
	for(int i = 0, p = 0; i < primes; i++)
	{
		//int rowHasOne = 0;
		for(int j = 0; j < number_count; j++)
		{
			bit_matrix[p][j] = 0;

			// Factor the number!
			while (mpz_divisible_p(nums_copy[j], prime))
			{
				mpz_divexact(nums_copy[j], nums_copy[j], prime);
				bit_matrix[p][j] = bit_matrix[p][j]^1;
				//rowHasOne = rowHasOne ^ 1;
			}
		}
			p++;

		/*if (rowHasOne == 1)
		{
		}
		else
		{
			bit_matrix_height--;
		}*/

		mpz_nextprime(prime, prime);
	}

	// If we have an overdetermined matrix, we must fail.
/*
	if (bit_matrix_height > bit_matrix_width)
	{
		#if VERBOSE
		printf("\n !!!! We have an overdetermined matrix with %d equations and %d unknowns !!! \n ", bit_matrix_height, bit_matrix_width);
		#endif
		return 0;
	}*/

	#if VERBOSE && MATRIX_VERBOSE
	printf("All these numbers should be 1: ");
	for(int i = 0; i < number_count; i++)
	{
		gmp_printf("%Zd ", nums_copy[i]);
	}
	printf("\n");
	#endif

	#if VERBOSE
	printf("\n Matrix is of size %d x %d\n", bit_matrix_width, bit_matrix_height);
	#if MATRIX_VERBOSE
	for(int i = 0; i < bit_matrix_height; i++)
	{
		printf(" ");
		for(int j = 0; j < bit_matrix_width; j++)
		{
			printf("%d", bit_matrix[i][j]);
		}
		printf("\n");
	}
	#endif
	#endif

	mpz_clear(prime);

	#if VERBOSE
	printf("\n Will now solve the system of equations built from the factors...\n ");
	#endif

	// Gauss elimination
	for(int column = 0, row = 0; column < bit_matrix_width; column++, row++)
	{
		#if VERBOSE && MATRIX_VERBOSE
		printf("\n \tLooking at column %d and row %d in matrix:\n ", column, row);

		if (bit_matrix_height < 100 && bit_matrix_width < 100)
		{
			for(int column = 0; column < bit_matrix_height; column++)
			{
				printf("\t");
				for(int row = 0; row < bit_matrix_width; row++)
				{
					printf("%d", bit_matrix[column][row]);
				}
				printf("\n ");
			}
			printf("\n ");
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

		#if VERBOSE && MATRIX_VERBOSE
		printf("\tFound a 1 on row %d\n ", maxRow);
		#endif

		// If we must replace the largest row to the top, swap them.
		if (maxRow != row)
		{
			#if VERBOSE && MATRIX_VERBOSE
			printf("\tSwapping rows %d and %d...\n ", row, maxRow);
			#endif

			// Swap row i and maxColumn
			char tmp[smoothness_bound];
			for(int c = 0; c < bit_matrix_width; c++)
			{
				tmp[c] = bit_matrix[row][c];
				bit_matrix[row][c] = bit_matrix[maxRow][c];
				bit_matrix[maxRow][c] = tmp[c];
			}
		}

		#if VERBOSE && MATRIX_VERBOSE
		printf("\tXOR-ing row %d with rows... ", row);
		#endif

		// Make sure all rows below this row has an initial zero.
		for(int r = row+1; r < bit_matrix_height; r++)
		{
			if (bit_matrix[r][column] == 0)
				continue;

			#if VERBOSE && MATRIX_VERBOSE
			printf("%d ", r);
			#endif

			// Subtract bit_matrix[k][row] * bit_matrix[column] from bit_matrix[k]
			for(int c = 0; c < bit_matrix_width; c++)
			{
				bit_matrix[r][c] = bit_matrix[r][c] ^ bit_matrix[row][c];
			}
		}
		#if VERBOSE && MATRIX_VERBOSE
		printf("\n ");
		#endif
	}
	#if VERBOSE && MATRIX_VERBOSE
	if (bit_matrix_height < 100 && bit_matrix_width < 100)
	{
		printf("\n After gauss:\n \n ");
		for(int column = 0; column < bit_matrix_height; column++)
		{
			for(int row = 0; row < bit_matrix_width; row++)
			{
				printf("%d", bit_matrix[column][row]);
			}
			printf("\n ");
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
	printf("\n We have %d known variables, thus there are %d unknowns.\n", known, bit_matrix_width-known);
	#endif

	int unknowns = bit_matrix_width-known;

	int visited_threshold = 2;
	mpz_t visited[visited_threshold];
	int v_ptr = 0;

	mpz_t permutations_of_unknowns;
	mpz_init(permutations_of_unknowns);
	mpz_ui_pow_ui(permutations_of_unknowns, 2, (unknowns > MAX_NUMBER_OF_SOLUTION_VECTORS ? MAX_NUMBER_OF_SOLUTION_VECTORS : unknowns));

	mpz_t ret1, ret2;
	mpz_init(ret1), mpz_init(ret2);

	// For all 2^unknowns permutations
	mpz_t modified_n;
	mpz_init_set(modified_n, num);

	mpz_t i;
	for(mpz_init_set_ui(i, 1); mpz_cmp(i, permutations_of_unknowns) < 0; mpz_add_ui(i, i, 1))
	{
		// Idea: i can be used with masks to get the current value of the unknowns.
		char solution[bit_matrix_width];

		// Fill the unknowns first
		int unknown_cntr = 0;
		for(int s = 0; s < bit_matrix_width; s++)
		{
			if (knownIndexes[s] != s)
			{
				solution[s] = (char)mpz_tstbit(i, unknown_cntr);
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

		#if VERBOSE && SOLUTION_ARRAY_VERBOSE
		printf(" Solution array: ");
		for(int s = 0; s < bit_matrix_width; s++)
		{
			printf("%d", solution[s]);
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
		mpz_set_ui(ret2,1);
		for(int s = 0; s < bit_matrix_width; s++)
		{
			if (solution[s] == 0)
				continue;

			mpz_mul(ret2, ret2, nums_p[s]);
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
		try_adding_factor_to_result(result, ret1, &modified_n, visited, &v_ptr);

		try_adding_factor_to_result(result, ret2, &modified_n, visited, &v_ptr);

		if (v_ptr > 0)
		{
			break;
		}
	}

	for(int s = 0; s < bit_matrix_width; s++){
		mpz_clear(nums[s]);
		mpz_clear(nums_p[s]);
	}

	// Clear our variables!
	mpz_clear(ret1), mpz_clear(ret2), mpz_clear(tmp), mpz_clear(mod);

	#if VERBOSE
	printf(" :: QS over and out.\n\n");
	#endif

	if (v_ptr == 0)
		return 0;
	else
		return 1;
}

int try_adding_factor_to_result(factor_list ** result, mpz_t factor, mpz_t * ofNumber, mpz_t visited[], int * visited_length)
{
	if (mpz_cmp_ui(factor, 1) == 0)
	{
		return 0;
	}
	if (mpz_cmp(factor, *ofNumber) >= 0)
	{
		return 0;
	}
	for(int i = 0; i < *(visited_length); i++)
	{
		if (mpz_cmp(visited[i], factor) == 0)
		{
			return 0;
		}
	}

	if (mpz_probab_prime_p(factor, 5) && mpz_divisible_p(*ofNumber, factor) != 0)
	{
		#if VERBOSE
		gmp_printf(" Found factor %Zd, which is a prime number.\n", factor);
		#endif

		mpz_t * v = malloc(sizeof(mpz_t));
		mpz_init_set(*v, factor);
		factor_list_add(result, v);

		mpz_init_set(visited[*(visited_length)], factor);
		(*visited_length)++;

		mpz_divexact(*ofNumber, *ofNumber, factor);

		if (mpz_probab_prime_p(*ofNumber, 5))
		{
			mpz_t * v = malloc(sizeof(mpz_t));
			mpz_init_set(*v, *ofNumber);
			factor_list_add(result, v);

			mpz_init_set(visited[*(visited_length)], *ofNumber);
			(*visited_length)++;
		}
	}
	else
	{
		#if VERBOSE
		gmp_printf(" Found factor %Zd, using Pollard's Rho to find the prime factors.\n", factor, factor);
		#endif

		factor_list * pollards_factors = malloc(sizeof(factor_list));
		pollards_factors->value = NULL;
		pollards_factors->next = NULL;

		pollard(&pollards_factors, factor);

		#if VERBOSE
		gmp_printf("\n \t%Zd = {Pollard} = ", factor, factor);
		#endif
		// Go through all factors and try to add them
		while(pollards_factors->value != NULL)
		{
			int success = try_adding_factor_to_result(result, *(pollards_factors->value), ofNumber, visited, visited_length);

			#if VERBOSE
			if (success)
			{
				gmp_printf("%Zd %s", *(pollards_factors->value), pollards_factors->next->value == NULL ? "\n" : "* ");
			}
			else
			{
				gmp_printf("(%Zd) %s", *(pollards_factors->value), pollards_factors->next->value == NULL ? "\n" : "* ");
			}
			#endif

			pollards_factors = pollards_factors->next;
		}
	}

	return 1;
}
