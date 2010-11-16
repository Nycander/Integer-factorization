#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void factor(mpz_t * num);

int main()
{
	mpz_t num;
	mpz_init(num);
	int i = 0;
	while(gmp_scanf("%Zd", num) > 0)
	{
		factor(&num);
		printf("\n");
		i++;
		mpz_init(num);
	}
	return i;
}


void factor(mpz_t * num)
{
	//gmp_printf("%Zd\n", num);
	printf("fail\n");
}
