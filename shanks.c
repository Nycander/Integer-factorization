#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <limits.h>

#include "shanks.h"


//Returns prime if fail otherwise return R satisfying R^2=n mod prime
mpz_t * shanks_tonelli(const mpz_t quad_resi, const mpz_t prime){
	if(mpz_legendre(quad_resi, prime) == 1){
		mpz_t Q, tmp;
		mpz_init(Q);
		mpz_init(tmp);
		mpz_t * R = malloc(sizeof(mpz_t));
		mpz_init(*R);
		mpz_sub_ui(Q, prime, 1);
		unsigned int S = 0;

		while(mpz_divisible_ui_p(Q, 2)){
			mpz_divexact_ui(Q, Q, 2);
			S++;
		}
		//If S == 1, R = +- n ^ ((p+1)/4)
		if(S == 1){
			mpz_add_ui(tmp, prime, 1);
			mpz_cdiv_q_ui(tmp,tmp,4);
			mpz_powm(*R,quad_resi,tmp,prime);
			mpz_clear(tmp);
			mpz_clear(Q);
			return R;
		}

		mpz_t z;
		mpz_init_set_ui(z, 0);
		while((mpz_legendre(z, prime)) != -1){
			mpz_add_ui(z,z,1);
		}

		//Step 4
		mpz_t c,t,M,fac2,b,tmp2;
		mpz_init(c);
		mpz_init(t);
		mpz_init(b);
		mpz_init_set_ui(fac2,2);
		mpz_init_set_ui(M, S);
		mpz_powm(c, z, Q, prime);
		mpz_clear(z);

		mpz_add_ui(tmp, Q,1);
		mpz_div_ui(tmp,tmp,2);
		mpz_powm(*R, quad_resi, tmp, prime);

		mpz_powm(t, quad_resi, Q, prime);

		mpz_clear(Q);
		mpz_init_set_ui(tmp2, 0);
		int i;
		while(1){
			//t = 1 mod prime
			if(mpz_cmp_ui(t,1) == 0){
				mpz_clear(c);
				mpz_clear(t);
				mpz_clear(M);
				mpz_clear(b);
				mpz_clear(fac2);
				mpz_clear(tmp2);
				mpz_clear(tmp);
				return R;
			}
			i=0;
			while(1){
				if(mpz_cmp_ui(M,i) == 0)
				{
					mpz_clear(c);
					mpz_clear(t);
					mpz_clear(M);
					mpz_clear(b);
					mpz_clear(fac2);
					mpz_clear(tmp2);
					mpz_clear(tmp);
					mpz_clear(*R);
					mpz_t * fail = malloc(sizeof(mpz_t));
					mpz_init_set(*fail, prime);
					return fail;
				}
				//Check if t^2^i = 1 mod prime
				mpz_pow_ui(tmp2, fac2, i);
				mpz_powm(tmp2, t, tmp2, prime);
				if(mpz_cmp_ui(tmp2,1)==0){
					// b = c^2^(M-tmp-1)
					mpz_sub_ui(tmp2,M,i);
					mpz_sub_ui(tmp2,tmp2,1);
					mpz_powm(tmp2, fac2, tmp2,prime); //MAY INTRODUCE BUG
					mpz_powm(b, c, tmp2, prime);
					// R = R*b
					mpz_mul(*R,*R,b);
					// t = t*b^2
					mpz_pow_ui(tmp2, b,2);
					mpz_mul(t,t,tmp2);
					//M = tmp
					mpz_set(c,tmp2);
					mpz_set_ui(M,i);
					break; //Test now
				}
				i++;
			}
		}
		mpz_clear(c);
		mpz_clear(t);
		mpz_clear(M);
		mpz_clear(b);
		mpz_clear(fac2);
		mpz_clear(tmp2);
		mpz_clear(tmp);
		mpz_clear(*R);
		mpz_t * fail = malloc(sizeof(mpz_t));
		mpz_init_set(*fail, prime);
		return fail;
	}
	else{
		mpz_t * fail = malloc(sizeof(mpz_t));
		mpz_init_set(*fail, prime);
		return fail;
	}

}
