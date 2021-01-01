/*

	clearn kernel

	Bryan Little April 2020

	Clears prime counter.

*/


__kernel void clearn(__global uint *primecount){

	int i = get_global_id(0);

	if(i==0){
		primecount[0]=0;
	}


}



