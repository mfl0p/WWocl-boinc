/* 
   prime.h -- Yves Gallot, Bryan Little 12-26-2020

*/


uint64_t invert(uint64_t p)
{
	uint64_t p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}


uint64_t montMul(uint64_t a, uint64_t b, uint64_t p, uint64_t q)
{
	unsigned __int128 res;

	res  = (unsigned __int128)a * b;
	uint64_t ab0 = (uint64_t)res;
	uint64_t ab1 = res >> 64;

	uint64_t m = ab0 * q;

	res = (unsigned __int128)m * p;
	uint64_t mp = res >> 64;

	uint64_t r = ab1 - mp;

	return ( ab1 < mp ) ? r + p : r;
}


uint64_t add(uint64_t a, uint64_t b, uint64_t p)
{
	uint64_t r;

	uint64_t c = (a >= p - b) ? p : 0;

	r = a + b - c;

	return r;
}


void montInit(uint64_t p, uint64_t* qq, uint64_t* o, uint64_t* pmo, uint64_t* r2)
{
	uint64_t q = invert(p);
	uint64_t one = (-p) % p;
	*pmo = p - one;
	
	uint64_t two = add(one, one, p);
	uint64_t t = add(two, two, p);
	for (int i = 0; i < 5; ++i)
		t = montMul(t, t, p, q);	// 4^{2^5} = 2^64
	*r2 = t;
	*o = one;
	*qq = q;
}


/* Used in the prime validator
   Returns 0 only if N is composite.
   Otherwise N is a strong probable prime to base a.
 */
inline int strong_prp(int base, uint64_t N, uint64_t q, uint64_t one, uint64_t mnmo, uint64_t r2)
{
	int retval = 0;
	uint64_t nmo = N-1;
	uint64_t a = base;
	int t = __builtin_ctzll(nmo);
	uint64_t exp = N >> t;
	uint64_t curBit = 0x8000000000000000;
	curBit >>= ( __builtin_clzll(exp) + 1 );

	/* If N is prime and N = d*2^t+1, where d is odd, then either
		1.  a^d = 1 (mod N), or
		2.  a^(d*2^s) = -1 (mod N) for some s in 0 <= s < t    */

	a = montMul(a,r2,N,q);  // convert base to montgomery form
	uint64_t mbase = a;

  	/* r <-- a^d mod N, assuming d odd */
	while( curBit )
	{
		a = montMul(a,a,N,q);

		if(exp & curBit){
			a = montMul(a,mbase,N,q);
		}

		curBit >>= 1;
	}

	/* Clause 1. and s = 0 case for clause 2. */
	if (a == one || a == mnmo){
		retval = 1;
	}

	/* 0 < s < t cases for clause 2. */
	for (int s = 1; !retval && s < t; ++s){

		a = montMul(a,a,N,q);

		if(a == mnmo){
	    		retval = 1;
		}
	}


	return retval;
}


