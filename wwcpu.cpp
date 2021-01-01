/* 

	wwcpu.cpp opencl kernels converted to x86_64 using GCC __int128

	Bryan Little, Yves Gallot, 12/26/2020

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <signal.h>

#include "WW.h"

typedef struct {
	uint64_t s0,s1;
}ulong2;

inline uint64_t mul_hi(uint64_t a, uint64_t b){

	unsigned __int128 res = (unsigned __int128)a * b;
	uint64_t r = res >> 64;
	return r;

}


// r0 + 2^64 * r1 = (a0 + 2^64 * a1) + (b0 + 2^64 * b1)
inline ulong2 add_wide(const ulong2 a, const ulong2 b)
{
	ulong2 r;

	r.s0 = a.s0 + b.s0;
	r.s1 = a.s1 + b.s1;

	if (r.s0 < a.s0) r.s1 += 1;

	return r;
}

// r0 + 2^64 * r1 = a * b
inline ulong2 mul_wide(const uint64_t a, const uint64_t b)
{
	ulong2 r;

	unsigned __int128 res = (unsigned __int128)a * b;

	r.s0 = (uint64_t)res;
	r.s1 = res >> 64;

	return r;
}

// a0 + 2^64 * a1 < b0 + 2^64 * b1
inline bool is_less_than(const ulong2 a, const ulong2 b)
{
	return (a.s1 < b.s1) || ((a.s1 == b.s1) && (a.s0 < b.s0));
}

// p * p_inv = 1 (mod 2^64) (Newton's method)
inline uint64_t invert(const uint64_t p)
{
	uint64_t p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

// r = x + y (mod p) where 0 <= r < p
inline uint64_t add_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t cp = (x >= p - y) ? p : 0;
	return x + y - cp;
}

// r = x + y (mod p) where 0 <= r < p, c is the carry
inline uint64_t add_mod_c(const uint64_t x, const uint64_t y, const uint64_t p, uint64_t * c)
{
	const bool carry = (x >= p - y);
	const uint64_t cp = carry ? p : 0;
	*c = carry ? 1 : 0;
	return x + y - cp;
}

// r = x - y (mod p) where 0 <= r < p
inline uint64_t sub_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t cp = (x < y) ? p : 0;
	return x - y + cp;
}

// r = x - y (mod p) where 0 <= r < p, c is the carry
inline uint64_t sub_mod_c(const uint64_t x, const uint64_t y, const uint64_t p, uint64_t * c)
{
	const bool carry = (x < y);
	const uint64_t cp = carry ? p : 0;
	*c = carry ? 1 : 0;
	return x - y + cp;
}

// "double-precision" variant Montgomery arithmetic. See:
// Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
// Dorais, F. G.; Klyve, D., "A Wieferich Prime Search Up to 6.7x10^15", Journal of Integer Sequences. 14 (9), 2011.

// 2^64 mod p^2 is (2^64, p^2) residue of 1
inline ulong2 m2p_one(const uint64_t p)
{
	ulong2 res;

	if ((p >> 32) == 0)
	{
		const uint64_t p2 = p * p, r_p2 = (-p2) % p2;	// 2^64 mod p^2
		res.s0 = r_p2 % p;
		res.s1 = r_p2 / p;
		return res;
	}
	// 2^64 mod p^2 = 2^64
	const uint64_t mp = -p;	// 2^64 - p
	res.s0 = mp % p;
	res.s1 = mp / p + 1;
	return res;
}

// r0 + p * r1 = 2 * (x0 + p * x1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_dup(const ulong2 x, const uint64_t p)
{
	ulong2 res;
	uint64_t c;
	const uint64_t l = add_mod_c(x.s0, x.s0, p, &c);
	const uint64_t h = add_mod(x.s1 + c, x.s1, p);
	res.s0 = l;
	res.s1 = h;
	return res;
}

// r0 + p * r1 = (x0 + p * x1)^2 (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_square(const ulong2 x, const uint64_t p, const uint64_t q)
{
	ulong2 res;
	const ulong2 t = mul_wide(x.s0, x.s0);
	const uint64_t u0 = q * t.s0;
	const uint64_t t1 = t.s1;
	const uint64_t v1 = mul_hi(p, u0);

	const ulong2 x01 = mul_wide(x.s0, x.s1);
	res.s0 = u0;
	res.s1 = 0;
	const ulong2 x01u = add_wide(x01, res);
	// 0 <= tp < 2p^2: 129 bits
	const ulong2 tp = add_wide(x01u, x01); bool tp_carry = is_less_than(tp, x01);
	// 0 <= tp_h < 2p. tp_h >= p if tp_h >= 2^64 or tp_h >= p
	const uint64_t tp_h = tp.s1, tpc = (tp_carry | (tp_h >= p)) ? p : 0;
	const uint64_t up0 = q * tp.s0;
	const uint64_t t1p = tp_h - tpc;	// 0 <= t1p < p
	const uint64_t v1p = mul_hi(p, up0);

	// 0 <= t1, v1 < p, 0 <= t1p, v1p < p
	uint64_t c;
	const uint64_t z0 = sub_mod_c(t1, v1, p, &c);
	const uint64_t z1 = sub_mod(t1p, v1p + c, p);
	res.s0 = z0;
	res.s1 = z1;
	return res;
}

// r0 + p * r1 = x * y (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_mul_s(const uint64_t x, const uint64_t y, const uint64_t p, const uint64_t q)
{
	ulong2 res;
	const ulong2 t = mul_wide(x, x);
	const uint64_t u0 = q * t.s0;
	const uint64_t t1 = t.s1;
	const uint64_t v1 = mul_hi(p, u0);

	const uint64_t v1p = mul_hi(p, q * u0);

	uint64_t c;
	const uint64_t z0 = sub_mod_c(t1, v1, p, &c);
	const uint64_t z1 = sub_mod(0, v1p + c, p);
	res.s0 = z0;
	res.s1 = z1;
	return res;
}

// To convert a residue to an integer, apply Algorithm REDC
inline ulong2 m2p_get(const ulong2 x, const uint64_t p, const uint64_t q)
{
	ulong2 res;
	const uint64_t u0 = q * x.s0;
	const uint64_t v1 = mul_hi(p, u0);

	const uint64_t tp = x.s1 + u0;
	const uint64_t up0 = q * tp;
	const uint64_t t1p = (tp < x.s1) ? 1 : 0;
	const uint64_t v1p = mul_hi(p, up0);

	uint64_t c;
	const uint64_t z0 = sub_mod_c(0, v1, p, &c);
	const uint64_t z1 = sub_mod(t1p, v1p + c, p);
	res.s0 = z0;
	res.s1 = z1;
	return res;
}


// r0 + p * r1 = (x0 + p * x1) * (y0 + p * y1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_mul(const ulong2 x, const ulong2 y, const uint64_t p, const uint64_t q)
{
	ulong2 res;
	const ulong2 t = mul_wide(x.s0, y.s0);
	const uint64_t u0 = q * t.s0;
	const uint64_t t1 = t.s1;
	const uint64_t v1 = mul_hi(p, u0);

	res.s0 = u0;
	res.s1 = 0;
	const ulong2 x0y1u = add_wide(mul_wide(x.s0, y.s1), res);
	const ulong2 x1y0 = mul_wide(x.s1, y.s0);
	// 0 <= tp < 2p^2: 129 bits
	const ulong2 tp = add_wide(x0y1u, x1y0); bool tp_carry = is_less_than(tp, x1y0);
	// 0 <= tp_h < 2p. tp_h >= p if tp_h >= 2^64 or tp_h >= p
	const uint64_t tp_h = tp.s1, tpc = (tp_carry | (tp_h >= p)) ? p : 0;
	const uint64_t up0 = q * tp.s0;
	const uint64_t t1p = tp_h - tpc;	// 0 <= t1p < p
	const uint64_t v1p = mul_hi(p, up0);

	// 0 <= t1, v1 < p, 0 <= t1p, v1p < p
	uint64_t c;
	const uint64_t z0 = sub_mod_c(t1, v1, p, &c);
	const uint64_t z1 = sub_mod(t1p, v1p + c, p);
	res.s0 = z0;
	res.s1 = z1;
	return res;
}


// r0 + p * r1 = (x0 + p * x1) + (y0 + p * y1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_add(const ulong2 x, const ulong2 y, const uint64_t p)
{
	ulong2 res;
	uint64_t c;
	const uint64_t l = add_mod_c(x.s0, y.s0, p, &c);
	const uint64_t h = add_mod(x.s1 + c, y.s1, p);
	res.s0 = l;
	res.s1 = h;
	return res;
}


int wieferichCPU(uint64_t p, uint64_t & checksum){

	const uint64_t q = invert(p);		// TODO read q from __global
	const ulong2 one = m2p_one(p);	// TODO read one from __global

	const uint64_t e = p >> 1;

	uint64_t curBit = 1ULL << (62 - __builtin_clzll(e));

	// a = 2
	ulong2 a = m2p_dup(one, p);

	// first step: 2^2 = 2 + 2
	a = m2p_dup(a, p);
	if ((e & curBit) != 0) a = m2p_dup(a, p);
	curBit >>= 1;

	while ((curBit != 0) && (a.s1 == 0))
	{
		a = m2p_mul_s(a.s0, a.s0, p, q);
		if ((e & curBit) != 0) a = m2p_dup(a, p);
		curBit >>= 1;
	}

	while (curBit != 0)
	{
		a = m2p_square(a, p, q);
		if ((e & curBit) != 0) a = m2p_dup(a, p);
		curBit >>= 1;
	}

	const ulong2 sp = m2p_get(a, p, q);

	// checksum calc: add sp.s1 mod 2^32 to checksum
	checksum += sp.s1 & UINT32_MAX;

	// store result
	if (((sp.s1 == 0) && (sp.s0 == 1)) || ((sp.s1 == p - 1) && (sp.s0 == p - 1)))
	{
		return 1;
	}
	else if (sp.s1 <= SPECIAL_THRESHOLD)
	{
		return 1;
	}
	else if (sp.s1 >= p - SPECIAL_THRESHOLD)
	{
		return 1;
	}

	return 0;
}

int wallsunsunCPU(uint64_t p, uint64_t & checksum){

	const uint64_t q = invert(p);		// TODO read q from __global
	const ulong2 one = m2p_one(p);	// TODO read one from __global

	const int pMod5 = (int)(p % 5);
	const int kro = (pMod5 == 2 || pMod5 == 3) ? -1 : 1;
	const uint64_t e = p - kro - 1;

	uint64_t curBit = 1ULL << (62 - __builtin_clzll(e));

	// b = F_1, c = F_0
	ulong2 res;
	res.s0 = 0;
	res.s1 = 0;
	ulong2 b = one, c = res;

	while ((curBit != 0) && (b.s1 == 0))
	{
		const ulong2 b2 = m2p_mul_s(b.s0, b.s0, p, q);
		const ulong2 c2 = m2p_mul_s(c.s0, c.s0, p, q);
		// F_{n + 1} = F_n + F_{n - 1}
		const ulong2 a = m2p_add(b, c, p);

		// F_{2n} = F_n (F_{n + 1} + F_{n - 1})
		const ulong2 ac = m2p_add(a, c, p);
		b = (ac.s1 == 0) ? m2p_mul_s(b.s0, ac.s0, p, q) : m2p_mul(b, ac, p, q);
		// F_{2n - 1} = F_n^2 + F_{n - 1}^2
		c = m2p_add(b2, c2, p);
	
		if ((e & curBit) != 0)
		{
			// F_{2n + 1} = F_{2n} + F_{2n - 1}
			const ulong2 a = m2p_add(b, c, p);
			c = b; b = a;
		}

		curBit >>= 1;
	}

	while (curBit != 0)
	{
		const ulong2 b2 = m2p_square(b, p, q);
		const ulong2 c2 = m2p_square(c, p, q);
		// F_{n + 1} = F_n + F_{n - 1}
		const ulong2 a = m2p_add(b, c, p);

		// F_{2n} = F_n (F_{n + 1} + F_{n - 1})
		b = m2p_mul(b, m2p_add(a, c, p), p, q);
		// F_{2n - 1} = F_n^2 + F_{n - 1}^2
		c = m2p_add(b2, c2, p);
	
		if ((e & curBit) != 0)
		{
			// F_{2n + 1} = F_{2n} + F_{2n - 1}
			const ulong2 a = m2p_add(b, c, p);
			c = b; b = a;
		}

		curBit >>= 1;
	}

	// F_{p - (p/5)}
	const ulong2 a = m2p_add(b, c, p);

//		const ulong2 c1 = m2p_get(b, p, q), 
	const ulong2 c2 = m2p_get(a, p, q);

	// If (p/5) = -1 then F_p should be -1 (mod p) and if (p/5) = 1 then F_{p - 2} should be 1 (mod p).
//		const uint64_t c10 = (kro == -1) ? p - c1.s0 : c1.s0, c11 = c1.s1, c20 = c2.s0, c21 = c2.s1;
//	const uint64_t c20 = c2.s0, c21 = c2.s1;
	const uint64_t c21 = c2.s1;

	// checksum calc: add c21 mod 2^32 to checksum
	checksum += c21 & UINT32_MAX;

	if (c21 == 0)
	{
		return 1;
	}
	else if (c21 <= SPECIAL_THRESHOLD)
	{
		return 1;
	}
	else if (c21 >= p - SPECIAL_THRESHOLD)
	{
		return 1;
	}

	return 0;
}






