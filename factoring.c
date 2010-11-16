#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void factor(mpz_t * num);

int main()
{
	mpz_t num;
	mpz_init(num);
	while(gmp_scanf("%Zd", num) > 0)
	{
		factor(&num);
		mpz_init(num);
	}
	return 0;
}


void factor(mpz_t * num)
{
	//gmp_printf("%Zd", num);
	printf("fail\n\n");
}