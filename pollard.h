#ifndef _POLLARD__H_
#define _POLLARD__H_

int pollard(factor ** f, const mpz_t n);
int rho(mpz_t res, const mpz_t N);

#endif
