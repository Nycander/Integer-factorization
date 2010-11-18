#ifndef __FACTOR_H_
#define __FACTOR_H_ 1

typedef struct factor
{
	mpz_t * value;
	struct factor * next;
} factor;

#endif