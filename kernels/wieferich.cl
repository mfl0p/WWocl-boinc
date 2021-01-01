/* 
   wieferich.cl -- Bryan Little, Yves Gallot, December 2020

   Wieferich search OpenCL Kernel portion

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

// r0 + p * r1 = 2 * (x0 + p * x1) (mod p^2) where 0 <= r0, r1 < p
inline ulong2 m2p_dup(const ulong2 x, const ulong p)
{
	ulong c;
	const ulong l = add_mod_c(x.s0, x.s0, p, &c);
	const ulong h = add_mod(x.s1 + c, x.s1, p);
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


__kernel void wieferich(
				__global uint *primeCount,
				__global const ulong *g_prime,
				__global const ulong *g_invert,
				__global ulong2 *g_one,
				__global ulong *specialPrime,
				__global int *rem,
				__global int *quot,
				__global int *result,
				__global ulong *checksum,
				__global uint *resultcount,
				__global ulong *totalprimecount
				)
{
	const uint gid = get_global_id(0);

	if (gid < primeCount[0])
	{
		const ulong p = g_prime[gid];
		const ulong q = g_invert[gid];
		const ulong2 one = m2p_one(p);
		g_one[gid] = one;

		const ulong e = p >> 1;

		ulong curBit = 1ul << (62 - clz(e));

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

		// store result
		if (((sp.s1 == 0) && (sp.s0 == 1)) || ((sp.s1 == p - 1) && (sp.s0 == p - 1)))
		{
			const int Ix = atomic_add(&resultcount[0], 1);
			specialPrime[Ix] = p;
			result[Ix] = 1093;
			rem[Ix] = 1;	// 1^2 or (-1)^2
			quot[Ix] = 0;
		}
		else if (sp.s1 <= SPECIAL_THRESHOLD)
		{
			if (sp.s0 == 1)
			{
				const int Ix = atomic_add(&resultcount[0], 1);
				specialPrime[Ix] = p;
				result[Ix] = 1093;
				rem[Ix] = 1;
				quot[Ix] = (int)(sp.s1);
			}
			else if (sp.s0 == p - 1)
			{
				const int Ix = atomic_add(&resultcount[0], 1);
				specialPrime[Ix] = p;
				result[Ix] = 1093;
				rem[Ix] = -1;
				quot[Ix] = (int)(sp.s1 + 1);
			}
		}
		else if (sp.s1 >= p - SPECIAL_THRESHOLD)
		{
			if (sp.s0 == 1)
			{
				int Ix = atomic_add(&resultcount[0], 1);
				specialPrime[Ix] = p;
				result[Ix] = 1093;
				rem[Ix] = 1;
				quot[Ix] = (int)(sp.s1 - p);
			}
			else if (sp.s0 == p - 1)
			{
				int Ix = atomic_add(&resultcount[0], 1);
				specialPrime[Ix] = p;
				result[Ix] = 1093;
				rem[Ix] = -1;
				quot[Ix] = (int)(sp.s1 + 1 - p);
			}
		}
/*
		ulong w_quot;
		if (sp.s0 == 1)
		{
			// (1 + sp.s1 * p)^2 = 1 + 2 * sp.s1 * p  (mod p^2)
			w_quot = add_mod(sp.s1, sp.s1, p);
		}
		else if (sp.s0 == p - 1)
		{
			// (p - 1 + sp.s1 * p)^2 = 1 - 2 * (sp.s1 + 1) * p  (mod p^2)
			w_quot = add_mod(sp.s1, 1, p);
			w_quot = sub_mod(0, add_mod(w_quot, w_quot, p), p);
		}
		else
		{
			a = m2p_square(a, p, q);
			const ulong2 w = m2p_get(a, p, q);
			w_quot = w.s1;
		}
*/
		// checksum calc: add c21 mod 2^32 to checksum
		atom_add(&checksum[0], sp.s1 & UINT32_MAX);

		if (gid == 0)
		{
			const uint pc = primeCount[0];

			// sum total prime count for host
			totalprimecount[0] += pc;

			// store largest primecount
			if (pc > primeCount[1]) primeCount[1] = pc;
		}
	}
}
