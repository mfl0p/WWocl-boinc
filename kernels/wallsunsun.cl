/* 
   wallsunsun.cl -- Bryan Little, Yves Gallot, December 2020

   WallSunSun search OpenCL Kernel portion

*/

#define SPECIAL_THRESHOLD 1000
#define UINT32_MAX 0xffffffff
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

// r0 + 2^64 * r1 = (a0 + 2^64 * a1) + (b0 + 2^64 * b1)
inline ulong2 add_wide(const ulong2 a, const ulong2 b)
{
	ulong2 r;

#ifdef __NV_CL_C_VERSION
	const uint a0 = (uint)(a.s0), a1 = (uint)(a.s0 >> 32), a2 = (uint)(a.s1), a3 = (uint)(a.s1 >> 32);
	const uint b0 = (uint)(b.s0), b1 = (uint)(b.s0 >> 32), b2 = (uint)(b.s1), b3 = (uint)(b.s1 >> 32);
	uint c0, c1, c2, c3;

	asm volatile ("add.cc.u32 %0, %1, %2;" : "=r" (c0) : "r" (a0), "r" (b0));
	asm volatile ("addc.cc.u32 %0, %1, %2;" : "=r" (c1) : "r" (a1), "r" (b1));
	asm volatile ("addc.cc.u32 %0, %1, %2;" : "=r" (c2) : "r" (a2), "r" (b2));
	asm volatile ("addc.u32 %0, %1, %2;" : "=r" (c3) : "r" (a3), "r" (b3));

	r.s0 = upsample(c1, c0); r.s1 = upsample(c3, c2);
#else
	r = a + b;
	if (r.s0 < a.s0) r.s1 += 1;
#endif

	return r;
}

// r0 + 2^64 * r1 = a * b
inline ulong2 mul_wide(const ulong a, const ulong b)
{
	ulong2 r;

#ifdef __NV_CL_C_VERSION
	const uint a0 = (uint)(a), a1 = (uint)(a >> 32);
	const uint b0 = (uint)(b), b1 = (uint)(b >> 32);

	uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c3) : "r" (a1), "r" (b1));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c3) : "r" (c3));

	r.s0 = upsample(c1, c0); r.s1 = upsample(c3, c2);
#else
	r.s0 = a * b; r.s1 = mul_hi(a, b);
#endif

	return r;
}

// a0 + 2^64 * a1 < b0 + 2^64 * b1
inline bool is_less_than(const ulong2 a, const ulong2 b)
{
	return (a.s1 < b.s1) || ((a.s1 == b.s1) && (a.s0 < b.s0));
}

// p * p_inv = 1 (mod 2^64) (Newton's method)
inline ulong invert(const ulong p)
{
	ulong p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

// r = x + y (mod p) where 0 <= r < p
inline ulong add_mod(const ulong x, const ulong y, const ulong p)
{
	const ulong cp = (x >= p - y) ? p : 0;
	return x + y - cp;
}

// r = x + y (mod p) where 0 <= r < p, c is the carry
inline ulong add_mod_c(const ulong x, const ulong y, const ulong p, ulong * c)
{
	const bool carry = (x >= p - y);
	const ulong cp = carry ? p : 0;
	*c = carry ? 1 : 0;
	return x + y - cp;
}

// r = x - y (mod p) where 0 <= r < p
inline ulong sub_mod(const ulong x, const ulong y, const ulong p)
{
	const ulong cp = (x < y) ? p : 0;
	return x - y + cp;
}

// r = x - y (mod p) where 0 <= r < p, c is the carry
inline ulong sub_mod_c(const ulong x, const ulong y, const ulong p, ulong * c)
{
	const bool carry = (x < y);
	const ulong cp = carry ? p : 0;
	*c = carry ? 1 : 0;
	return x - y + cp;
}

// "double-precision" variant Montgomery arithmetic. See:
// Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
// Dorais, F. G.; Klyve, D., "A Wieferich Prime Search Up to 6.7x10^15", Journal of Integer Sequences. 14 (9), 2011.

// 2^64 mod p^2 is (2^64, p^2) residue of 1
inline ulong2 m2p_one(const ulong p)
{
	if ((p >> 32) == 0)
	{
		const ulong p2 = p * p, r_p2 = (-p2) % p2;	// 2^64 mod p^2
		return (ulong2)(r_p2 % p, r_p2 / p);
	}
	// 2^64 mod p^2 = 2^64
	const ulong mp = -p;	// 2^64 - p
	return (ulong2)(mp % p, mp / p + 1);
}

// r0 + p * r1 = (x0 + p * x1) + (y0 + p * y1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_add(const ulong2 x, const ulong2 y, const ulong p)
{
	ulong c;
	const ulong l = add_mod_c(x.s0, y.s0, p, &c);
	const ulong h = add_mod(x.s1 + c, y.s1, p);
	return (ulong2)(l, h);
}

// r0 + p * r1 = (x0 + p * x1)^2 (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_square(const ulong2 x, const ulong p, const ulong q)
{
	const ulong2 t = mul_wide(x.s0, x.s0);
	const ulong u0 = q * t.s0;
	const ulong t1 = t.s1;
	const ulong v1 = mul_hi(p, u0);

	const ulong2 x01 = mul_wide(x.s0, x.s1);
	const ulong2 x01u = add_wide(x01, (ulong2)(u0, 0));
	// 0 <= tp < 2p^2: 129 bits
	const ulong2 tp = add_wide(x01u, x01); bool tp_carry = is_less_than(tp, x01);
	// 0 <= tp_h < 2p. tp_h >= p if tp_h >= 2^64 or tp_h >= p
	const ulong tp_h = tp.s1, tpc = (tp_carry | (tp_h >= p)) ? p : 0;
	const ulong up0 = q * tp.s0;
	const ulong t1p = tp_h - tpc;	// 0 <= t1p < p
	const ulong v1p = mul_hi(p, up0);

	// 0 <= t1, v1 < p, 0 <= t1p, v1p < p
	ulong c;
	const ulong z0 = sub_mod_c(t1, v1, p, &c);
	const ulong z1 = sub_mod(t1p, v1p + c, p);
	return (ulong2)(z0, z1);
}

// r0 + p * r1 = (x0 + p * x1) * (y0 + p * y1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_mul(const ulong2 x, const ulong2 y, const ulong p, const ulong q)
{
	const ulong2 t = mul_wide(x.s0, y.s0);
	const ulong u0 = q * t.s0;
	const ulong t1 = t.s1;
	const ulong v1 = mul_hi(p, u0);

	const ulong2 x0y1u = add_wide(mul_wide(x.s0, y.s1), (ulong2)(u0, 0));
	const ulong2 x1y0 = mul_wide(x.s1, y.s0);
	// 0 <= tp < 2p^2: 129 bits
	const ulong2 tp = add_wide(x0y1u, x1y0); bool tp_carry = is_less_than(tp, x1y0);
	// 0 <= tp_h < 2p. tp_h >= p if tp_h >= 2^64 or tp_h >= p
	const ulong tp_h = tp.s1, tpc = (tp_carry | (tp_h >= p)) ? p : 0;
	const ulong up0 = q * tp.s0;
	const ulong t1p = tp_h - tpc;	// 0 <= t1p < p
	const ulong v1p = mul_hi(p, up0);

	// 0 <= t1, v1 < p, 0 <= t1p, v1p < p
	ulong c;
	const ulong z0 = sub_mod_c(t1, v1, p, &c);
	const ulong z1 = sub_mod(t1p, v1p + c, p);
	return (ulong2)(z0, z1);
}

// r0 + p * r1 = x * y (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_mul_s(const ulong x, const ulong y, const ulong p, const ulong q)
{
	const ulong2 t = mul_wide(x, x);
	const ulong u0 = q * t.s0;
	const ulong t1 = t.s1;
	const ulong v1 = mul_hi(p, u0);

	const ulong v1p = mul_hi(p, q * u0);

	ulong c;
	const ulong z0 = sub_mod_c(t1, v1, p, &c);
	const ulong z1 = sub_mod(0, v1p + c, p);
	return (ulong2)(z0, z1);
}

// To convert a residue to an integer, apply Algorithm REDC
inline ulong2 m2p_get(const ulong2 x, const ulong p, const ulong q)
{
	const ulong u0 = q * x.s0;
	const ulong v1 = mul_hi(p, u0);

	const ulong tp = x.s1 + u0;
	const ulong up0 = q * tp;
	const ulong t1p = (tp < x.s1) ? 1 : 0;
	const ulong v1p = mul_hi(p, up0);

	ulong c;
	const ulong z0 = sub_mod_c(0, v1, p, &c);
	const ulong z1 = sub_mod(t1p, v1p + c, p);
	return (ulong2)(z0, z1);
}


__kernel void wallsunsun(
				__global const uint *primeCount,
				__global const ulong *g_prime,
				__global const ulong *g_invert,
				__global const ulong2 *g_one,
				__global ulong *specialPrime,
				__global int *rem,
				__global int *quot,
				__global int *result,
				__global ulong *checksum,
				__global uint *resultcount
				)
{
	const uint gid = get_global_id(0);

	if (gid < primeCount[0])
	{
		const ulong p = g_prime[gid];
		const ulong q = g_invert[gid];
		const ulong2 one = g_one[gid];

		const int pMod5 = (int)(p % 5);
		const int kro = (pMod5 == 2 || pMod5 == 3) ? -1 : 1;
		const ulong e = p - kro - 1;

		ulong curBit = 1ul << (62 - clz(e));

		// b = F_1, c = F_0
		ulong2 b = one, c = (ulong2)(0, 0);

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
//		const ulong c10 = (kro == -1) ? p - c1.s0 : c1.s0, c11 = c1.s1, c20 = c2.s0, c21 = c2.s1;
		const ulong c20 = c2.s0, c21 = c2.s1;

		// Upon return c10 should be 1 and c20 should be 0.
		// If not, then p is not prime or there is a bug in this routine.
/*		if (c10 != 1)
		{
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = -2;
			rem[Ix] = (int)(c10);
			quot[Ix] = (int)(c11);
		}
		else if (c20 != 0)
		{
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = -1;
			rem[Ix] = (int)(c20);
			quot[Ix] = (int)(c21);
		}
*/
		if (c21 == 0)
		{
			// If c21 is 0, then we have a WallSunSun prime.
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = 1;
			rem[Ix] = 0;
			quot[Ix] = 0;
		}
		else if (c21 <= SPECIAL_THRESHOLD)
		{
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = 2;
			rem[Ix] = 0;
			quot[Ix] = (int)(c21);
		}
		else if (c21 >= p - SPECIAL_THRESHOLD)
		{
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = 2;
			rem[Ix] = 0;
			quot[Ix] = (int)(c21 - p);
		}

		// checksum calc: add c21 mod 2^32 to checksum
		atom_add(&checksum[0], c21 & UINT32_MAX);
	}
}
