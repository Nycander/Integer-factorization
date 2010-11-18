#include <gmp.h>

#include "trialdivision.h"
#include "factor.h"

#include "primes.h" // Should contain the array primes[]

int trial_division(factor ** f, const mpz_t n)
{
	for(int i = 0; i < primes.length; i++)
	{
		if (mpz_divisible_ui_p(n, primes[i]))
		{
			return 0;
		}
	}
	return 0;
}