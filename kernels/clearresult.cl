/*

	clearresult kernel

	Bryan Little April 2020

	Clears checksum, result count, total prime count, and largest prime count.

*/


__kernel void clearresult(__global ulong *checksum,
			  __global uint *resultcount,
			  __global ulong *totalprimecount,
			  __global uint *primecount){

	int i = get_global_id(0);

	if(i == 0){
		checksum[0] = 0;
		resultcount[0] = 0;
		totalprimecount[0] = 0;
		primecount[1] = 0;
	}

}

