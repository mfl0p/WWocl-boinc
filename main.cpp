/*
	WWocl
	Bryan Little, Yves Gallot, 12/26/2020

*/

#define MAJORV 1
#define MINORV 0
#define SUFFIXV "4"

#include <iostream>
#include <cinttypes>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>

# include "boinc_api.h"
# include "filesys.h"
# include "boinc_opencl.h"

#include "simpleCL.h"

#include "clearn.h"
#include "clearresult.h"
#include "getsegprimes.h"
#include "wieferich.h"
#include "wallsunsun.h"

#include "prime.h"
#include "WW.h"

#include "include/primesieve.h"

#define RESULTS_FILENAME "results-WW.txt"
#define STATE_FILENAME_A "stateA-WW.txt"
#define STATE_FILENAME_B "stateB-WW.txt"

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std; 


typedef struct {
	cl_mem gpuResult = NULL;
	cl_mem gpuSpecialPrime = NULL;
	cl_mem gpuChecksum = NULL;
	cl_mem gpuRem = NULL;
	cl_mem gpuQuot = NULL;
	cl_mem gpuResultCount = NULL;
	cl_mem d_segprime = NULL;
	cl_mem d_invert = NULL;
	cl_mem d_one = NULL;
	cl_mem d_segprimecount = NULL;
	cl_mem d_totalprimecount = NULL;

	sclSoft clearn, wieferich, wallsunsun, getsegprimes, clearresult;
	sclHard hardware;
}progData;


void cleanup( progData pd ){

	sclReleaseMemObject(pd.gpuSpecialPrime);
	sclReleaseMemObject(pd.gpuChecksum);
	sclReleaseMemObject(pd.gpuResult);
	sclReleaseMemObject(pd.gpuRem);
	sclReleaseMemObject(pd.gpuQuot);
	sclReleaseMemObject(pd.gpuResultCount);
	sclReleaseMemObject(pd.d_segprime);
	sclReleaseMemObject(pd.d_invert);
	sclReleaseMemObject(pd.d_one);
	sclReleaseMemObject(pd.d_segprimecount);
	sclReleaseMemObject(pd.d_totalprimecount);

        sclReleaseClSoft(pd.clearn);
	sclReleaseClSoft(pd.clearresult);
	sclReleaseClSoft(pd.getsegprimes);
        sclReleaseClSoft(pd.wieferich);
        sclReleaseClSoft(pd.wallsunsun);

        sclReleaseClHard(pd.hardware);
}


FILE *my_fopen(const char * filename, const char * mode)
{
	char resolved_name[512];

	boinc_resolve_filename(filename,resolved_name,sizeof(resolved_name));
	return boinc_fopen(resolved_name,mode);

}


void write_state(uint64_t start, uint64_t stop, uint64_t current, uint64_t cksm, double psum, int pcnt, bool & write_state_a_next)
{
	FILE *out;

        if (write_state_a_next){
		if ((out = my_fopen(STATE_FILENAME_A,"w")) == NULL)
			fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_A);
	}
	else{
                if ((out = my_fopen(STATE_FILENAME_B,"w")) == NULL)
                        fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_B);
        }
	if (fprintf(out,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %lf %d\n",start,stop,current,cksm,psum,pcnt) < 0){
		if (write_state_a_next)
			fprintf(stderr,"Cannot write to %s !!! Continuing...\n",STATE_FILENAME_A);
		else
			fprintf(stderr,"Cannot write to %s !!! Continuing...\n",STATE_FILENAME_B);

		// Attempt to close, even though we failed to write
		fclose(out);
	}
	else{
		// If state file is closed OK, write to the other state file
		// next time round
		if (fclose(out) == 0) 
			write_state_a_next = !write_state_a_next; 
	}
}

/* Return 1 only if a valid checkpoint can be read.
   Attempts to read from both state files,
   uses the most recent one available.
 */
int read_state(uint64_t start, uint64_t stop, uint64_t & current, uint64_t & cksm, double & psum, int & pcnt, bool & write_state_a_next)
{
	FILE *in;
	bool good_state_a = true;
	bool good_state_b = true;
	uint64_t tmp1, tmp2;
	uint64_t current_a, cksm_a, current_b, cksm_b;
	double psum_a, psum_b;
	int pcnt_a, pcnt_b;

        // Attempt to read state file A
	if ((in = my_fopen(STATE_FILENAME_A,"r")) == NULL){
		good_state_a = false;
        }
	else if (fscanf(in,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %lf %d\n",&tmp1,&tmp2,&current_a,&cksm_a,&psum_a,&pcnt_a) != 6){
		fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_A);
		good_state_a = false;
	}
	else{
		fclose(in);
		/* Check that start stop match */
		if (tmp1 != start || tmp2 != stop){
			good_state_a = false;
		}
	}

        // Attempt to read state file B
        if ((in = my_fopen(STATE_FILENAME_B,"r")) == NULL){
                good_state_b = false;
        }
        else if (fscanf(in,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %lf %d\n",&tmp1,&tmp2,&current_b,&cksm_b,&psum_b,&pcnt_b) != 6){
                fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_B);
                good_state_b = false;
        }
        else{
                fclose(in);
		/* Check that start stop match */
		if (tmp1 != start || tmp2 != stop){
                        good_state_b = false;
                }
        }

        // If both state files are OK, check which is the most recent
	if (good_state_a && good_state_b)
	{
		if (current_a > current_b)
			good_state_b = false;
		else
			good_state_a = false;
	}

        // Use data from the most recent state file
	if (good_state_a && !good_state_b)
	{
		current = current_a;
		cksm = cksm_a;
		psum = psum_a;
		pcnt = pcnt_a;
		write_state_a_next = false;
		return 1;
	}
        if (good_state_b && !good_state_a)
        {
                current = current_b;
                cksm = cksm_b;
		psum = psum_b;
		pcnt = pcnt_b;
		write_state_a_next = true;
		return 1;
        }

	// If we got here, neither state file was good
	return 0;
}


void checkpoint(uint64_t start, uint64_t stop, uint64_t current, uint64_t cksm, double psum, int pcnt, bool & write_state_a_next)
{
	write_state(start, stop, current, cksm, psum, pcnt, write_state_a_next);

	if(boinc_is_standalone()){
		printf("Checkpoint, current: %" PRIu64 "\n", current);
	}

	boinc_checkpoint_completed();
}


#ifdef _WIN32
double getSysOpType()
{
    double ret = 0.0;
    NTSTATUS(WINAPI *RtlGetVersion)(LPOSVERSIONINFOEXW);
    OSVERSIONINFOEXW osInfo;

    *(FARPROC*)&RtlGetVersion = GetProcAddress(GetModuleHandleA("ntdll"), "RtlGetVersion");

    if (NULL != RtlGetVersion)
    {
        osInfo.dwOSVersionInfoSize = sizeof(osInfo);
        RtlGetVersion(&osInfo);
        ret = (double)osInfo.dwMajorVersion;
    }
    return ret;
}
#endif


// N is a strong pseudoprime to base 2
int isPRP(uint64_t N)
{
	if (N%2==0)
		return 0;

	// check if the first few primes divide N
	// Must match CPU program to maintain checksum
	if(N % 3)
	if(N % 5)
	if(N % 7)
	if(N % 11)
	if(N % 13)
	if(N % 17)
	if(N % 19)
	if(N % 23)
	if(N % 29)
	if(N % 31)
	if(N % 37)
	if(N % 41)
	if(N % 43)
	if(N % 47)
	if(N % 53)
	if(N % 59)
	if(N % 61)
	if(N % 67)
	if(N % 71)
	if(N % 73)
	if(N % 79)
	if(N % 83)
	if(N % 89)
	if(N % 97)
	if(N % 101)
	if(N % 103)
	if(N % 107)
	if(N % 109)
	if(N % 113){

		uint64_t q, one, r2, nmo;

		montInit(N, &q, &one, &nmo, &r2);

		if( strong_prp(2, N, q, one, nmo, r2) )
			return 1;
	}

	return 0;
}


// N is prime.
// Miller–Rabin primality test
//    if n <  3,825,123,056,546,413,051, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, and 23.
//                  3825123056546413051
//    if n < 18,446,744,073,709,551,616, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
//                 18446744073709551616
// minimum is 29
int isPrime(uint64_t N)
{
	const int base[12] = {2,3,5,7,11,13,17,19,23,29,31,37};

	if (N % 2==0)
		return 0;

	uint64_t q, one, r2, nmo;

	montInit(N, &q, &one, &nmo, &r2);

	if ( N < 3825123056546413051ULL ){
		for (int i = 0; i < 9; ++i)
			if (!strong_prp(base[i], N, q, one, nmo, r2))
				return 0;
	}
	else if ( N <= UINT64_MAX ){
		for (int i = 0; i < 12; ++i)
			if (!strong_prp(base[i], N, q, one, nmo, r2))
				return 0;
	}

	return 1;
}

void sleepCPU(sclHard hardware){

	cl_event kernelsDone;
	cl_int err;
	cl_int info;
	struct timespec sleep_time;
	sleep_time.tv_sec = 0;
	sleep_time.tv_nsec = 1000000;	// 1ms

	err = clEnqueueMarker( hardware.queue, &kernelsDone);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clEnqueueMarker\n");
		fprintf(stderr, "ERROR: clEnqueueMarker\n");
		sclPrintErrorFlags(err); 
	}

	clFlush(hardware.queue);

	err = clGetEventInfo(kernelsDone, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clGetEventInfo\n" );
		fprintf(stderr, "ERROR: clGetEventInfo\n" );
		sclPrintErrorFlags( err );
       	}

	// sleep until event complete
	while(info >= 0 && info != CL_COMPLETE){
		nanosleep(&sleep_time,NULL);
		err = clGetEventInfo(kernelsDone, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
		if ( err != CL_SUCCESS ) {
			printf( "ERROR: clGetEventInfo\n" );
			fprintf(stderr, "ERROR: clGetEventInfo\n" );
			sclPrintErrorFlags( err );
	       	}
	}
	clReleaseEvent(kernelsDone);

}


// test gpu prp generator with cpu
// also returns GPU's PRP count
uint32_t checkPRPS( progData pd, uint64_t range_start, uint64_t range_stop ){

	uint32_t * h_segprimecount = (uint32_t *)malloc(2*sizeof(uint32_t));
	if( h_segprimecount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// get prime count
	sclRead(pd.hardware, 2 * sizeof(uint32_t), pd.d_segprimecount, h_segprimecount);

	// allocate host prime array
	uint64_t * h_segprime = (uint64_t *)malloc(h_segprimecount[0]*sizeof(uint64_t));
	if( h_segprime == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy array
	sclRead(pd.hardware, h_segprimecount[0]*sizeof(uint64_t), pd.d_segprime, h_segprime);

	uint32_t prpcount = h_segprimecount[0];

	// correct the GPU prp array to a prime array
	for(uint32_t y=0; y < prpcount; ++y){
		if( !isPrime(h_segprime[y]) ){
			if(isPRP(h_segprime[y])){
				// GPU found a PRP, set value to max and decrement the prime count
				// set to max will let the prime fall off during sort()
				h_segprime[y] = UINT64_MAX;
				h_segprimecount[0]--;
			}
			else{
				// GPU calculated an invalid result
				fprintf(stderr, "Computation error, GPU found a composite: %" PRIu64 "\n", h_segprime[y]);
				exit(EXIT_FAILURE);
			}
		}
	}

	// sort the GPU prime array
	sort(h_segprime, h_segprime + prpcount);

	// get a CPU generated prime list for the same range
	// -1 because gpu checks # < stop
	// primesieve checks # <= stop
	size_t listsize;
	uint64_t *cprimes = (uint64_t*)primesieve_generate_primes(range_start, range_stop-1, &listsize, UINT64_PRIMES);
	uint32_t cpucount = (uint32_t)listsize;

	// compare prime counts
	if(h_segprimecount[0] != cpucount){
		fprintf(stderr, "ERROR: gpu prime count: %u does not match cpu: %u\n", h_segprimecount[0], cpucount);
		printf("ERROR: gpu prime count: %u does not match cpu: %u\n", h_segprimecount[0], cpucount);
		exit(EXIT_FAILURE);
	}

	// compare the primes
	for(uint32_t y=0; y < cpucount; ++y){
		if(cprimes[y] != h_segprime[y]){
			fprintf(stderr, "ERROR: gpu/cpu prime mismatch %" PRIu64 " vs %" PRIu64 "\n", h_segprime[y], cprimes[y]);
			printf("ERROR: gpu/cpu prime mismatch %" PRIu64 " vs %" PRIu64 "\n",h_segprime[y], cprimes[y]);
			exit(EXIT_FAILURE);
		}
	}

	free(h_segprimecount);
	free(h_segprime);
	free(cprimes);

	return prpcount;

}


// find mod 30 wheel index based on starting N
// this is used by gpu threads to iterate over the number line
void findWheelOffset(uint64_t & start, int32_t & index){

	int32_t wheel[8] = {4, 2, 4, 2, 4, 6, 2, 6};
	int32_t idx = -1;

	// find starting number using mod 6 wheel
	// N=(k*6)-1, N=(k*6)+1 ...
	// where k, k+1, k+2 ...
	uint64_t k = start / 6;
	int32_t i = 1;
	uint64_t N = (k * 6)-1;


	while( N < start || N % 5 == 0 ){
		if(i){
			i = 0;
			N += 2;
		}
		else{
			i = 1;
			N += 4;
		}
	}

	start = N;

	// find mod 30 wheel index by iterating with a mod 6 wheel until finding N divisible by 5
	// forward to find index
	while(idx < 0){

		if(i){
			N += 2;
			i = 0;
			if(N % 5 == 0){
				N -= 2;
				idx = 5;
			}

		}
		else{
			N += 4;
			i = 1;
			if(N % 5 == 0){
				N -= 4;
				idx = 7;
			}
		}
	}

	// reverse to find starting index
	while(N != start){
		--idx;
		if(idx < 0)idx=7;
		N -= wheel[idx];
	}


	index = idx;

}


void testKernels( progData& pd, uint64_t start, uint64_t top, uint64_t calc_range, int debuginfo ){

	uint32_t * hostResultCount = (uint32_t *)malloc(sizeof(uint32_t));
	if( hostResultCount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}
	uint64_t * hostChecksum = (uint64_t *)malloc(sizeof(uint64_t));
	if( hostChecksum == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}
	uint32_t gpuprimecount = 0;
	uint32_t resultcount = 0;
	uint32_t primecount = 0;
	uint64_t checksum = 0;
	uint64_t Prime = start;
	uint64_t stop;

	if(2000000 > calc_range){
		stop = start + calc_range;
	}
	else{
		stop = start + 2000000;
	}

	if(stop < start){
		stop = top;
	}

	fprintf(stderr, "Verifying GPU kernels using range %" PRIu64 " to %" PRIu64 "...\n", start, stop);
	if(debuginfo){
		printf("Verifying GPU kernels using range %" PRIu64 " to %" PRIu64 "...\n", start, stop);
	}

	// clear result arrays
	sclEnqueueKernel(pd.hardware, pd.clearresult);

	// clear segprime count
	sclEnqueueKernel(pd.hardware, pd.clearn);

	int32_t wheelidx;
	uint64_t kernel_start = start;

	findWheelOffset(kernel_start, wheelidx);

	// get primes
	sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &kernel_start);
	sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &stop);
	sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &wheelidx);
	sclEnqueueKernel(pd.hardware, pd.getsegprimes);

	// wieferich test
	sclEnqueueKernel(pd.hardware, pd.wieferich);

	// wallsunsun test
	sclEnqueueKernel(pd.hardware, pd.wallsunsun);

	// run the same tests with CPU
	if(Prime % 2 == 0){
		++Prime;
	}

	for(; Prime < stop; Prime += 2){
		if(Prime % 3)
		if(Prime % 5)
		if(Prime % 7)
		if(Prime % 11)
		if(Prime % 13)
		if(Prime % 17)
		if(Prime % 19)
		if(Prime % 23)
		if(Prime % 29)
		if(Prime % 31)
		if(Prime % 37)
		if(Prime % 41)
		if(Prime % 43)
		if(Prime % 47)
		if(Prime % 53)
		if(Prime % 59)
		if(Prime % 61)
		if(Prime % 67)
		if(Prime % 71)
		if(Prime % 73)
		if(Prime % 79)
		if(Prime % 83)
		if(Prime % 89)
		if(Prime % 97)
		if(Prime % 101)
		if(Prime % 103)
		if(Prime % 107)
		if(Prime % 109)
		if(Prime % 113){
			uint64_t q, one, r2, nmo;

			montInit(Prime, &q, &one, &nmo, &r2);

			if( strong_prp(2, Prime, q, one, nmo, r2) ){
				++primecount;
				if(wieferichCPU(Prime, checksum)){
					++resultcount;
				}
				if(wallsunsunCPU(Prime, checksum)){
					++resultcount;
				}
			}
		}
	}


	// compare GPU primes with primesieve generated list
	gpuprimecount = checkPRPS( pd, start, stop );

	if(gpuprimecount != primecount){
		fprintf(stderr, "ERROR: gpu / cpu prp count mismatch %u / %u\n", gpuprimecount, primecount);
		printf("ERROR: gpu / cpu prp count mismatch %u / %u\n", gpuprimecount, primecount);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "GPU prime generator kernel test passed. (%u PRPs)\n", primecount);
	if(debuginfo){
		printf("GPU prime generator kernel test passed. (%u PRPs)\n", primecount);
	}

	// compare the checksum and result count
	sclRead(pd.hardware, sizeof(uint64_t), pd.gpuChecksum, hostChecksum);	
	sclRead(pd.hardware, sizeof(uint32_t), pd.gpuResultCount, hostResultCount);

	if(*hostChecksum != checksum){
		fprintf(stderr, "ERROR: gpu / cpu checksum mismatch %" PRIu64 " / %" PRIu64 "\n", *hostChecksum, checksum);
		printf("ERROR: gpu / cpu checksum mismatch %" PRIu64 " / %" PRIu64 "\n", *hostChecksum, checksum);
		exit(EXIT_FAILURE);
	}

	if(*hostResultCount != resultcount){
		fprintf(stderr, "ERROR: gpu / cpu result count mismatch %u / %u\n", *hostResultCount, resultcount);
		printf("ERROR: gpu / cpu result count mismatch %u / %u\n", *hostResultCount, resultcount);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "GPU Wieferich/WallSunSun kernel tests passed. (%" PRIu64 " checksum, %u results)\n", checksum, resultcount);
	if(debuginfo){
		printf("GPU Wieferich/WallSunSun kernel tests passed. (%" PRIu64 " checksum, %u results)\n", checksum, resultcount);
	}

	free(hostResultCount);
	free(hostChecksum);

}


void report_solution( char * results ){

	FILE * resfile = my_fopen(RESULTS_FILENAME,"a");

	if( resfile == NULL ){
		fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
		exit(EXIT_FAILURE);
	}

	if( fprintf( resfile, "%s", results ) < 0 ){
		fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
		exit(EXIT_FAILURE);
	}

	fflush(resfile);

#if defined (_WIN32)
	_commit(_fileno(resfile));
#else
	fsync(fileno(resfile));
#endif

	fclose(resfile);

}


void getResults( int ckpt_time, double& p_sum, int& p_cnt, uint64_t mem_size, uint32_t numresults, int debuginfo, uint64_t& checksum, progData pd ){

	uint64_t * h_checksum = (uint64_t *)malloc(sizeof(uint64_t));
	if( h_checksum == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}
	uint32_t * h_resultcount = (uint32_t *)malloc(sizeof(uint32_t));
	if( h_resultcount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}
	uint32_t * h_segprimecount = (uint32_t *)malloc(2*sizeof(uint32_t));
	if( h_segprimecount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}
	uint64_t * h_totalprimecount = (uint64_t *)malloc(sizeof(uint64_t));
	if( h_totalprimecount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy results to host memory
	// blocking read
	sclRead(pd.hardware, sizeof(uint64_t), pd.gpuChecksum, h_checksum);
	sclRead(pd.hardware, sizeof(uint64_t), pd.d_totalprimecount, h_totalprimecount);
	sclRead(pd.hardware, sizeof(uint32_t), pd.gpuResultCount, h_resultcount);
	sclRead(pd.hardware, 2 * sizeof(uint32_t), pd.d_segprimecount, h_segprimecount);

	// add checksum
	checksum += *h_checksum;

	// make sure prime array size was calculated correctly
	if( h_segprimecount[1] > mem_size ){
		fprintf(stderr,"Error: prime array overflow!\n");
		exit(EXIT_FAILURE);
	}

	uint32_t PRPcnt = 0;

	if(*h_resultcount > 0){

		if(*h_resultcount > numresults){
			fprintf(stderr,"Error: number of results (%u) overflowed array.\n", *h_resultcount);
			exit(EXIT_FAILURE);
		}

		uint64_t * hostSpecialPrime = (uint64_t *)malloc(*h_resultcount * sizeof(uint64_t));
		if( hostSpecialPrime == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}
		int32_t * hostResult = (int32_t *)malloc(*h_resultcount * sizeof(int32_t));
		if( hostResult == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}
		int32_t * hostRem = (int32_t *)malloc(*h_resultcount * sizeof(int32_t));
		if( hostRem == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}
		int32_t * hostQuot = (int32_t *)malloc(*h_resultcount * sizeof(int32_t)); 
		if( hostQuot == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}

		sclRead(pd.hardware, *h_resultcount * sizeof(uint64_t), pd.gpuSpecialPrime, hostSpecialPrime);
		sclRead(pd.hardware, *h_resultcount * sizeof(int32_t), pd.gpuRem, hostRem);
		sclRead(pd.hardware, *h_resultcount * sizeof(int32_t), pd.gpuQuot, hostQuot);
		sclRead(pd.hardware, *h_resultcount * sizeof(int32_t), pd.gpuResult, hostResult);

		// sort results by prime size if needed
		if(*h_resultcount > 1){
			for (uint32_t i = 0; i < *h_resultcount-1; i++){    
				for (uint32_t j = 0; j < *h_resultcount-i-1; j++){
					if (hostSpecialPrime[j] > hostSpecialPrime[j+1]){
						swap(hostSpecialPrime[j], hostSpecialPrime[j+1]);
						swap(hostRem[j], hostRem[j+1]);
						swap(hostQuot[j], hostQuot[j+1]);
						swap(hostResult[j], hostResult[j+1]);
					}
				}
			}
		}


		char buffer[256];
		char * resbuff = (char *)malloc( *h_resultcount * sizeof(char) * 256 );
		if( resbuff == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}
		resbuff[0] = '\0';


		// write results
		for (uint32_t  i = 0; i < *h_resultcount; i++){

			if( isPrime(hostSpecialPrime[i]) ){
				// Wieferich Prime
				if(hostResult[i] == 1093){
					if (hostRem[i] == 1 && hostQuot[i] == 0){
						if( sprintf( buffer, "%" PRIu64 " is a Wieferich prime\n", hostSpecialPrime[i] ) < 0 ){
							fprintf(stderr,"error in sprintf()\n");
							exit(EXIT_FAILURE);
						}
						strcat( resbuff, buffer );
					}
					else{
						if( sprintf( buffer, "%" PRIu64 " is a Wieferich special instance (%+d %+d p)\n", hostSpecialPrime[i], hostRem[i], hostQuot[i]) < 0 ){
							fprintf(stderr,"error in sprintf()\n");
							exit(EXIT_FAILURE);
						}
						strcat( resbuff, buffer );
					}
				}
				// WallSunSun Primes
				else if(hostResult[i] == 1){
					if ( sprintf( buffer, "%" PRIu64 " is a WallSunSun prime\n", hostSpecialPrime[i]) < 0 ){
						fprintf(stderr,"error in sprintf()\n");
						exit(EXIT_FAILURE);
					}	
					strcat( resbuff, buffer );			
				}
				else if(hostResult[i] == 2){
					if ( sprintf( buffer, "%" PRIu64 " is a WallSunSun special instance (+0 %+d p)\n", hostSpecialPrime[i], hostQuot[i]) < 0 ){
						fprintf(stderr,"error in sprintf()\n");
						exit(EXIT_FAILURE);
					}
					strcat( resbuff, buffer );				
				}
				else{
					fprintf(stderr,"Error: hostResult returned %d.\n", hostResult[i]);
					exit(EXIT_FAILURE);				
				}
			}
			else if( isPRP(hostSpecialPrime[i]) ){
				// GPU does a base 2 PRP test only.
				// GPU calculations are OK, but we can discard this result
				++PRPcnt;
			}
			else{
				// GPU calculated an invalid result
				fprintf(stderr, "Computation error, GPU found a composite result: %" PRIu64 "\n", hostSpecialPrime[i]);
				exit(EXIT_FAILURE);
			}

		}


		report_solution( resbuff );

		free(hostSpecialPrime);
		free(hostResult);
		free(hostRem);
		free(hostQuot);
		free(resbuff);

	}


	if(debuginfo){
		printf("PRPs tested since last checkpoint: %" PRIu64 "\n", *h_totalprimecount);
		printf("Max number of results:             %u\n",numresults);
		printf("Number of results:                 %u\n",*h_resultcount);
		printf("PRPs removed from results:         %u\n",PRPcnt);
		printf("Largest kernel PRP count:          %u\n",h_segprimecount[1]);
		printf("Max kernel PRP count for malloc:   %" PRIu64 "\n",mem_size);
	}


	if(ckpt_time == 0){
		ckpt_time = 1;
	}

	double p_sec = ((double)(*h_totalprimecount / ckpt_time)) / 1000000.0;

	p_sum += p_sec;
	++p_cnt;

	free(h_checksum);
	free(h_resultcount);
	free(h_segprimecount);
	free(h_totalprimecount);

}


void profileGPU(progData& pd, int computeunits, int COMPUTE, int clworksize, uint64_t current, uint64_t& crange, int debuginfo, uint64_t& msize ){

	// calculate approximate chunk size based on gpu's compute units
	cl_int err = 0;
	uint64_t multiplier;

	if(COMPUTE){
		multiplier = 22000000;
	}
	else{
		multiplier = 2200000;
	}

	uint64_t calc_range = computeunits * multiplier;

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	uint64_t estimated = calc_range;

	uint64_t prof_start = current;

	// don't profile at very low or high N
	if(prof_start < 100000000){
		prof_start = 100000000;
	}
	else if(prof_start > 18446743073709551615ULL){
		prof_start = 18446743073709551615ULL;
	}

	uint64_t prof_stop = prof_start + calc_range;

	// Benchmark the GPU
	double total_time_ms = 0.0;
	double kernel_ms = 0.0;

	sclSetGlobalSize( pd.getsegprimes, (calc_range/60)+1 );

	// get a count of primes in the gpu worksize
	uint64_t prof_range_primes = primesieve_count_primes( prof_start, prof_stop );
	// calculate prime array size based on result
	double prof_mem_multi = 1.5 * ((double)prof_range_primes / (double)calc_range);
	uint64_t prof_mem_size = (uint64_t)(prof_mem_multi * (double)calc_range);

	// wieferich and wallsunsun kernels use uint for global id
	if(prof_mem_size > UINT32_MAX){
		fprintf(stderr, "ERROR: prof_mem_size too large.\n");
                printf( "ERROR: prof_mem_size too large.\n" );
		exit(EXIT_FAILURE);
	}

	// allocate temporary gpu prime array for profiling
	cl_mem d_profileprime = clCreateBuffer( pd.hardware.context, CL_MEM_READ_WRITE, prof_mem_size*sizeof(uint64_t), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
	        printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	cl_mem d_profileinvert = clCreateBuffer( pd.hardware.context, CL_MEM_READ_WRITE, prof_mem_size*sizeof(uint64_t), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
	        printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

	int32_t prof_wheelidx;
	uint64_t prof_kernel_start = prof_start;

	findWheelOffset(prof_kernel_start, prof_wheelidx);

	// set static args
	sclSetKernelArg(pd.clearn, 0, sizeof(cl_mem), &pd.d_segprimecount);

	sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &prof_kernel_start);
	sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &prof_stop);
	sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &prof_wheelidx);
	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &d_profileprime);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &d_profileinvert);
	sclSetKernelArg(pd.getsegprimes, 5, sizeof(cl_mem), &pd.d_segprimecount);

	// first kernel profile is ignored
	sclEnqueueKernel(pd.hardware, pd.clearn);
	kernel_ms = ProfilesclEnqueueKernel(pd.hardware, pd.getsegprimes);

	// spin up some kernels while profiling
	for(int w = 0; w < 3; w++){

		sclEnqueueKernel(pd.hardware, pd.clearn);

		kernel_ms = ProfilesclEnqueueKernel(pd.hardware, pd.getsegprimes);

		total_time_ms += kernel_ms;
	}

	// avg of the 3 profiles
	double prof_avg_ms = total_time_ms / 3.0;
	double prof_multi;

	if(COMPUTE){
		// target kernel time is 100ms
		prof_multi = 100.0 / prof_avg_ms;
	}
	else if(clworksize){
		// target kernel time is specified on command line, between 1 and 100ms
		prof_multi = (double)clworksize / prof_avg_ms;
	}
	else{
		// target kernel time is 11ms
		prof_multi = 11.0 / prof_avg_ms;
	}

	// update chunk size based on the profile
	double new_range = (double)calc_range * prof_multi;
	calc_range = (uint64_t)new_range;

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	fprintf(stderr,"Kernel profile: %0.3f ms. Estimated / Actual worksize: %" PRIu64 " / %" PRIu64 "\n",prof_avg_ms,estimated,calc_range);
	if(debuginfo){
		printf("Kernel profile: %0.3f ms. Estimated / Actual worksize: %" PRIu64 " / %" PRIu64 "\n",prof_avg_ms,estimated,calc_range);
	}

	// get a count of primes in the new gpu worksize
	uint64_t range_primes = primesieve_count_primes( prof_start, prof_start+calc_range );

	// calculate prime array size based on result
	double mem_multi = 1.5 * ((double)range_primes / (double)calc_range);
	uint64_t mem_size = (uint64_t)(mem_multi * (double)calc_range);

	// wieferich and wallsunsun kernels use uint for global id
	if(mem_size > UINT32_MAX){
		fprintf(stderr, "ERROR: mem_size too large.\n");
                printf( "ERROR: mem_size too large.\n" );
		exit(EXIT_FAILURE);
	}

	crange = calc_range;
	msize = mem_size;

	// free temporary arrays
	sclReleaseMemObject(d_profileprime);
	sclReleaseMemObject(d_profileinvert);

}


int main(int argc, char *argv[])
{ 
	uint64_t bottom;
	uint64_t top;
	uint64_t current;
	uint64_t cksm;
	uint64_t calc_range;
	double p_sum = 0.0;
	int p_cnt = 0;
	uint32_t numresults = 15000;
	bool write_state_a_next;
	int computeunits;
	int testPRPgen = 0;
	int debuginfo = 0;
	int COMPUTE = 0;
	int clworksize = 0;
	time_t boinc_last, boinc_curr;
	time_t ckpt_curr, ckpt_last;

	progData pd;

	primesieve_set_num_threads(1);

        // Initialize BOINC
        BOINC_OPTIONS options;
        boinc_options_defaults(options);
        options.normal_thread_priority = true;
        boinc_init_options(&options);

	fprintf(stderr, "\nWWocl version %d.%d%s by Bryan Little and Yves Gallot\n",MAJORV,MINORV,SUFFIXV);
	fprintf(stderr, "Compiled " __DATE__ " with GCC " __VERSION__ "\n");
	fprintf(stderr, "A Wieferich and WallSunSun prime number search program\n");
	fprintf(stderr, "Based on wwwwcl by Mark Rodenkirch and primesieve by Kim Walisch\n");
	if(boinc_is_standalone()){
		printf("WWocl version %d.%d%s by Bryan Little and Yves Gallot\n",MAJORV,MINORV,SUFFIXV);
		printf("Compiled " __DATE__ " with GCC " __VERSION__ "\n");
		printf("A Wieferich and WallSunSun prime number search program\n");
		printf("Based on wwwwcl by Mark Rodenkirch and primesieve by Kim Walisch\n");	
	}

        // Print out cmd line for diagnostics
        fprintf(stderr, "Command line: ");
        for (int i = 0; i < argc; i++)
        	fprintf(stderr, "%s ", argv[i]);
        fprintf(stderr, "\n");

	/* Get search parameters from command line */
	if(argc < 3){
		printf("Usage: %s START END -testprimes -compute -debug -ws #\nWhere 127 <= START < END < 2^64\n",argv[0]);
		printf("-testprimes is optional and will compare all GPU generated primes with CPU generated primes (very slow)\n\n");
		printf("-compute is optional and will fully load the gpu with large kernel and queue sizes.  Enabled by default with Windows 10 and Pascal or newer NVIDIA GPU.\n");
		printf(" Compute flag is NOT recommended for GPUs that have a monitor attached when using Linux, MacOS, or older versions of Windows.\n\n");
		printf("-debug is optional and will display more information in standalone mode.\n\n");
		printf("-ws # is optional worksize for adjusting the GPU load higher for faster computation or lower for better screen refresh rate.\n");
		printf(" If using compute mode, -ws is ignored.  Usage: -ws 11 for default, -ws 1 for lowest GPU load, -ws 100 for highest GPU load.\n\n");
		printf("Example format: %s 10000000 20000000\n",argv[0]);
		printf("            Or: %s 1e7 2e7\n",argv[0]);
		exit(EXIT_FAILURE);
	}


	for(int w=3; w<argc; ++w){
		if( strcmp(argv[w], "-compute") == 0 ){
			COMPUTE = 1;
			fprintf(stderr, "Compute mode enabled.\n");
			printf("Compute mode enabled.\n");
		}
		else if( strcmp(argv[w], "-testprimes") == 0 ){
			if(boinc_is_standalone()){  // only allow this test in stand alone mode
				testPRPgen = 1;
				fprintf(stderr, "Comparing all GPU generated primes with CPU generated primes.\n");
				printf("Comparing all GPU generated primes with CPU generated primes.\n");
			}
		}
		else if( strcmp(argv[w], "-debug") == 0 ){
			if(boinc_is_standalone()){  // only allow this test in stand alone mode
				debuginfo = 1;
				printf("Debug info enabled.\n");
			}
		}
	}


	// get START value
	char * a1 = strstr(argv[1], (char*)"e");
	if(a1 != NULL){
		// exp notation
		uint64_t num = 0;
		int32_t exp = 0;

		num = strtoull(argv[1], &a1, 0);
		exp = strtoul(a1+1, NULL, 0);

		if(num == 0 || exp == 0){
			printf("ERROR: unable to parse START\n");
			fprintf(stderr, "ERROR: unable to parse START\n");
			exit(EXIT_FAILURE);
		}

		for(int i=0; i < exp; ++i){
			if(num > UINT64_MAX/10){
				printf("ERROR: START is out of range\n");
				fprintf(stderr, "ERROR: START is out of range\n");
				exit(EXIT_FAILURE);
			}
			else{
				num *= 10;
			}
		}
		
		bottom = num;

	}
	else{
		// std notation
		sscanf(argv[1],"%" PRIu64 "",&bottom);
	}

	
	// get END value
	char * a2 = strstr(argv[2], (char*)"e");
	if(a2 != NULL){
		// exp notation
		uint64_t num = 0;
		int32_t exp = 0;

		num = strtoull(argv[2], &a2, 0);
		exp = strtoul(a2+1, NULL, 0);

		if(num == 0 || exp == 0){
			printf("ERROR: unable to parse END\n");
			fprintf(stderr, "ERROR: unable to parse END\n");
			exit(EXIT_FAILURE);
		}

		for(int i=0; i < exp; ++i){
			if(num > UINT64_MAX/10){
				printf("ERROR: END is out of range\n");
				fprintf(stderr, "ERROR: END is out of range\n");
				exit(EXIT_FAILURE);
			}
			else{
				num *= 10;
			}
		}
		
		top = num;		
	}
	else{
		// std notation
		sscanf(argv[2],"%" PRIu64 "",&top);
	}


	// check search range
	if(bottom < 127){
		printf("START below minimum. Starting search at 127\n");
		fprintf(stderr, "START below minimum. Starting search at 127\n");
		bottom = 127;
	}
	if(top <= bottom){
		printf("ERROR: END must be > START\n");
		fprintf(stderr, "ERROR: END must be > START\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "Starting search at: %" PRIu64 "\nStopping search at: %" PRIu64 "\n", bottom,top);
	if(boinc_is_standalone()){
		printf("Starting search at: %" PRIu64 "\nStopping search at: %" PRIu64 "\n", bottom,top);
	}

	cl_platform_id platform = 0;
	cl_device_id device = 0;
	cl_context ctx;
	cl_command_queue queue;
	cl_int err = 0;

	int retval = 0;
	retval = boinc_get_opencl_ids(argc, argv, 0, &device, &platform);
	if (retval) {
		if(boinc_is_standalone()){
			printf("init_data.xml not found, using device 0.\n");

			err = clGetPlatformIDs(1, &platform, NULL);
			if (err != CL_SUCCESS) {
				printf( "clGetPlatformIDs() failed with %d\n", err );
				exit(EXIT_FAILURE);
			}
			err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
			if (err != CL_SUCCESS) {
				printf( "clGetDeviceIDs() failed with %d\n", err );
				exit(EXIT_FAILURE);
			}
		}
		else{
			fprintf(stderr, "Error: boinc_get_opencl_ids() failed with error %d\n", retval );
			exit(EXIT_FAILURE);
		}
	}

	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

	ctx = clCreateContext(cps, 1, &device, NULL, NULL, &err);
	if (err != CL_SUCCESS) {
		fprintf(stderr, "Error: clCreateContext() returned %d\n", err);
        	exit(EXIT_FAILURE); 
   	}

	queue = clCreateCommandQueue(ctx, device, CL_QUEUE_PROFILING_ENABLE, &err);	
	//queue = clCreateCommandQueue(ctx, device, 0, &err);
	if(err != CL_SUCCESS) { 
		fprintf(stderr, "Error: Creating Command Queue. (clCreateCommandQueue) returned %d\n", err );
		exit(EXIT_FAILURE);
    	}

	pd.hardware.platform = platform;
	pd.hardware.device = device;
	pd.hardware.queue = queue;
	pd.hardware.context = ctx;

 	char device_string0[1024];
 	char device_string1[1024];
 	char device_string2[1024];
	cl_uint CUs;

	clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_string0), &device_string0, NULL);
	clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(device_string1), &device_string1, NULL);
	clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(device_string2), &device_string2, NULL);
	clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &CUs, NULL);

	fprintf(stderr, "GPU Info:\n  Name: \t\t%s\n  Vendor: \t\t%s\n  Driver: \t\t%s\n  Compute Units: \t%u\n", device_string0, device_string1, device_string2, CUs);
	if(boinc_is_standalone()){
		printf("GPU Info:\n  Name: \t\t%s\n  Vendor: \t\t%s\n  Driver: \t\t%s\n  Compute Units: \t%u\n", device_string0, device_string1, device_string2, CUs);
	}

	// check vendor and normalize compute units
	// kernel size will be determined by profiling so this doesn't have to be accurate.
	computeunits = (int)CUs;
	char amd_s[] = "Advanced";
	char intel_s[] = "Intel";
	char nvidia_s[] = "NVIDIA";
	
	if(strstr((char*)device_string1, (char*)nvidia_s) != NULL){
#ifdef _WIN32
		float winVer = (float)getSysOpType();

		if(winVer >= 10.0f && !COMPUTE){

		 	cl_uint ccmajor;
			err = clGetDeviceInfo(pd.hardware.device, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(ccmajor), &ccmajor, NULL);
			if ( err != CL_SUCCESS ) {
				printf( "Error checking device compute capability\n" );
				fprintf(stderr, "Error checking device compute capability\n");
				exit(EXIT_FAILURE);
			}

			if(ccmajor >= 6){  // pascal or newer gpu on windows 10
				fprintf(stderr, "Detected Windows %0.1f, using compute preemption.\n", winVer);
				if(debuginfo){
					printf( "Detected Windows %0.1f, using compute preemption.\n", winVer );
				}
				COMPUTE = 1;
			}
		}
#endif
	}
	else if(strstr((char*)device_string1, (char*)amd_s) != NULL){
		computeunits /= 2;
	}
	else if(strstr((char*)device_string1, (char*)intel_s) != NULL){
		computeunits /= 30;
	}
        else{
		computeunits /= 2;
        }

	if(computeunits < 1){
		computeunits = 1;
	}


	// check for command line worksize override
	if(!COMPUTE){
		for(int w=3; w<argc; ++w){
			if( strcmp(argv[w], "-ws") == 0 ){
				if(w+1 < argc){
					int wso;
					sscanf(argv[w+1],"%d",&wso);
					if(wso > 0 && wso <= 100){
						fprintf(stderr, "Kernel work size of %d specified on command line.\n", wso);
						if(debuginfo){
							printf( "Kernel work size of %d specified on command line.\n", wso);
						}
						clworksize = wso;
					}
				}
			}
		}
	}


        pd.d_totalprimecount = clCreateBuffer( ctx, CL_MEM_READ_WRITE, sizeof(uint64_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_segprimecount = clCreateBuffer( ctx, CL_MEM_READ_WRITE, 2*sizeof(uint32_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.gpuSpecialPrime = clCreateBuffer( ctx, CL_MEM_READ_WRITE, numresults*sizeof(uint64_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.gpuChecksum = clCreateBuffer( ctx, CL_MEM_READ_WRITE, sizeof(uint64_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.gpuResultCount = clCreateBuffer( ctx, CL_MEM_READ_WRITE, sizeof(uint32_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.gpuResult = clCreateBuffer( ctx, CL_MEM_READ_WRITE,  numresults*sizeof(int32_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.gpuRem = clCreateBuffer( ctx, CL_MEM_READ_WRITE,  numresults*sizeof(int32_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.gpuQuot = clCreateBuffer( ctx, CL_MEM_READ_WRITE,  numresults*sizeof(int32_t), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}


        pd.clearn = sclGetCLSoftware(clearn_cl,"clearn",pd.hardware, 1, debuginfo);

        pd.wieferich = sclGetCLSoftware(wieferich_cl,"wieferich",pd.hardware, 1, debuginfo);

        pd.wallsunsun = sclGetCLSoftware(wallsunsun_cl,"wallsunsun",pd.hardware, 1, debuginfo);

        pd.clearresult = sclGetCLSoftware(clearresult_cl,"clearresult",pd.hardware, 1, debuginfo);

        pd.getsegprimes = sclGetCLSoftware(getsegprimes_cl,"getsegprimes",pd.hardware, 1, debuginfo);


	// kernel has __attribute__ ((reqd_work_group_size(256, 1, 1)))
	// it's still possible the CL complier picked a different size
	if(pd.getsegprimes.local_size[0] != 256){
		pd.getsegprimes.local_size[0] = 256;
		fprintf(stderr, "Set getsegprimes kernel local size to 256\n");
		if(debuginfo){
			printf("Set getsegprimes kernel local size to 256\n");
		}
	}


	// Resume from checkpoint if there is one
	if (read_state(bottom,top,current,cksm,p_sum,p_cnt,write_state_a_next)){
		if(boinc_is_standalone()){
			printf("Resuming search from checkpoint. Current: %" PRIu64 "\n",current);
		}
		fprintf(stderr,"Resuming search from checkpoint. Current: %" PRIu64 "\n",current);
	}
	// starting from beginning
	else{
		if(boinc_is_standalone()){
			printf("Beginning a new search with parameters from the command line\n");
		}
		current = bottom;
		cksm = 0;
                write_state_a_next = true;

		// clear result file
		FILE * temp_file = my_fopen(RESULTS_FILENAME,"w");
		if (temp_file == NULL){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
		fclose(temp_file);
	}

	//trying to resume a finished workunit
	if(current == top){
		if(boinc_is_standalone()){
			printf("Workunit complete.\n");
		}
		fprintf(stderr,"Workunit complete.\n");
		boinc_finish(EXIT_SUCCESS);
		return 0;
	}

	// global size for clearn kernel
	sclSetGlobalSize( pd.clearn, 64 );

	// global size for clear result kernel
	sclSetGlobalSize( pd.clearresult, 64 );

	// get global work size by profiling gpu
	uint64_t mem_size;
	profileGPU( pd, computeunits, COMPUTE, clworksize, current, calc_range, debuginfo, mem_size );

	// allocate gpu prime array
	pd.d_segprime = clCreateBuffer( ctx, CL_MEM_READ_WRITE, mem_size*sizeof(uint64_t), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_invert = clCreateBuffer( ctx, CL_MEM_READ_WRITE, mem_size*sizeof(uint64_t), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_one = clCreateBuffer( ctx, CL_MEM_READ_WRITE, mem_size*sizeof(cl_ulong2), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

	// global size for prime kernel
	sclSetGlobalSize( pd.getsegprimes, (calc_range/60)+1 );

	// global size for wieferich and wallsunsun kernels
	sclSetGlobalSize( pd.wieferich, mem_size );
	sclSetGlobalSize( pd.wallsunsun, mem_size );


	if(debuginfo){
		printf("kernel global sizes:\n  getsegprimes\t%" PRIu64 "\n  wieferich\t%" PRIu64 "\n  wallsunsun\t%" PRIu64 "\n",
			pd.getsegprimes.global_size[0], pd.wieferich.global_size[0], pd.wallsunsun.global_size[0]);
	}

	// set static kernel args
	sclSetKernelArg(pd.clearresult, 0, sizeof(cl_mem), &pd.gpuChecksum);
	sclSetKernelArg(pd.clearresult, 1, sizeof(cl_mem), &pd.gpuResultCount);
	sclSetKernelArg(pd.clearresult, 2, sizeof(cl_mem), &pd.d_totalprimecount);
	sclSetKernelArg(pd.clearresult, 3, sizeof(cl_mem), &pd.d_segprimecount);

	sclSetKernelArg(pd.clearn, 0, sizeof(cl_mem), &pd.d_segprimecount);

	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &pd.d_segprime);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &pd.d_invert);
	sclSetKernelArg(pd.getsegprimes, 5, sizeof(cl_mem), &pd.d_segprimecount);

	sclSetKernelArg(pd.wieferich, 0, sizeof(cl_mem), &pd.d_segprimecount);
	sclSetKernelArg(pd.wieferich, 1, sizeof(cl_mem), &pd.d_segprime);
	sclSetKernelArg(pd.wieferich, 2, sizeof(cl_mem), &pd.d_invert);
	sclSetKernelArg(pd.wieferich, 3, sizeof(cl_mem), &pd.d_one);
	sclSetKernelArg(pd.wieferich, 4, sizeof(cl_mem), &pd.gpuSpecialPrime);
	sclSetKernelArg(pd.wieferich, 5, sizeof(cl_mem), &pd.gpuRem);
	sclSetKernelArg(pd.wieferich, 6, sizeof(cl_mem), &pd.gpuQuot);
	sclSetKernelArg(pd.wieferich, 7, sizeof(cl_mem), &pd.gpuResult);
	sclSetKernelArg(pd.wieferich, 8, sizeof(cl_mem), &pd.gpuChecksum);
	sclSetKernelArg(pd.wieferich, 9, sizeof(cl_mem), &pd.gpuResultCount);
	sclSetKernelArg(pd.wieferich, 10, sizeof(cl_mem), &pd.d_totalprimecount);

	sclSetKernelArg(pd.wallsunsun, 0, sizeof(cl_mem), &pd.d_segprimecount);
	sclSetKernelArg(pd.wallsunsun, 1, sizeof(cl_mem), &pd.d_segprime);
	sclSetKernelArg(pd.wallsunsun, 2, sizeof(cl_mem), &pd.d_invert);
	sclSetKernelArg(pd.wallsunsun, 3, sizeof(cl_mem), &pd.d_one);
	sclSetKernelArg(pd.wallsunsun, 4, sizeof(cl_mem), &pd.gpuSpecialPrime);
	sclSetKernelArg(pd.wallsunsun, 5, sizeof(cl_mem), &pd.gpuRem);
	sclSetKernelArg(pd.wallsunsun, 6, sizeof(cl_mem), &pd.gpuQuot);
	sclSetKernelArg(pd.wallsunsun, 7, sizeof(cl_mem), &pd.gpuResult);
	sclSetKernelArg(pd.wallsunsun, 8, sizeof(cl_mem), &pd.gpuChecksum);
	sclSetKernelArg(pd.wallsunsun, 9, sizeof(cl_mem), &pd.gpuResultCount);


	if(debuginfo){

		printf("Profiling kernels...\n");
		int32_t widx;
		uint64_t k_start = current;
		uint64_t k_stop = current + calc_range;

		if(k_stop > top || k_stop < k_start){  // ck overflow
			k_stop = top;
		}

		findWheelOffset(k_start, widx);
		sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &k_start);
		sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &k_stop);
		sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &widx);

		sclEnqueueKernel(pd.hardware, pd.clearresult);
		sclEnqueueKernel(pd.hardware, pd.clearn);
		double prime_ms = ProfilesclEnqueueKernel(pd.hardware, pd.getsegprimes);
		double wie_ms = ProfilesclEnqueueKernel(pd.hardware, pd.wieferich);
		double wall_ms = ProfilesclEnqueueKernel(pd.hardware, pd.wallsunsun);
		double total_ms = prime_ms + wie_ms + wall_ms;

		int prime_pct = (int)( (prime_ms / total_ms)*100.0 );
		int wie_pct = (int)( (wie_ms / total_ms)*100.0 );
		int wall_pct = (int)( (wall_ms / total_ms)*100.0 );

		printf("Percent of time spent in each kernel:\n");
		printf("\tPrime Generator:\t%d %%\t(%0.1f ms)\n\tWieferich Test:\t\t%d %%\t(%0.1f ms)\n\tWallSunSun Test:\t%d %%\t(%0.1f ms)\n",prime_pct,prime_ms,wie_pct,wie_ms,wall_pct,wall_ms);
	}


	// verify proper operation of the compiled kernels
	testKernels(pd, current, top, calc_range, debuginfo);

	fprintf(stderr,"Starting search...\n");
	if(boinc_is_standalone()){
		printf("Starting search...\n");
	}

	time(&boinc_last);
	time(&ckpt_last);

	// clear result arrays
	sclEnqueueKernel(pd.hardware, pd.clearresult);

	int in_queue = 0;
	uint64_t stop;

	for( uint64_t start = current; start < top; start = stop ){

		stop = start + calc_range;

		if(stop > top || stop < start){  // ck overflow
			stop = top;
		}

		// update BOINC fraction done every 2 sec
		time(&boinc_curr);
		if( ((int)boinc_curr - (int)boinc_last) > 1 ){
    			double fd = (double)(start-bottom)/(double)(top-bottom);
			boinc_fraction_done(fd);
			if(boinc_is_standalone()){
				printf("Tests done: %.4f%%\n",fd*100.0);
			}
			boinc_last = boinc_curr;
		}

		// clear segprime count
		sclEnqueueKernel(pd.hardware, pd.clearn);

		int32_t wheelidx;
		uint64_t kernel_start = start;

		findWheelOffset(kernel_start, wheelidx);

		// get primes
		sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &kernel_start);
		sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &stop);
		sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &wheelidx);
		sclEnqueueKernel(pd.hardware, pd.getsegprimes);

		// command line argument to check all GPU primes with CPU for testing
		if(testPRPgen){
			uint32_t prpcnt = checkPRPS( pd, start, stop );
			printf("GPU primes match CPU (%u PRPs)\n", prpcnt);
		}

		// wieferich test
		sclEnqueueKernel(pd.hardware, pd.wieferich);

		// wallsunsun test
		sclEnqueueKernel(pd.hardware, pd.wallsunsun);

		// sleep cpu thread
		if(COMPUTE){
			if( ++in_queue == 3 ){  // about 1 second of kernels queued
				in_queue = 0;
				sleepCPU(pd.hardware);
			}
		}
		else{
			sleepCPU(pd.hardware);
		}


		// 1 minute checkpoint
		time(&ckpt_curr);
		int ckpt_diff = (int)ckpt_curr - (int)ckpt_last;

		if( ckpt_diff > 60 || stop == top ){

			sleepCPU(pd.hardware);
			in_queue = 0;

			boinc_begin_critical_section();

			// get results
			getResults( ckpt_diff, p_sum, p_cnt, mem_size, numresults, debuginfo, cksm, pd );

			if( stop == top ){
				// print checksum
				char buffer[256];
				if( sprintf( buffer, "%016" PRIX64 "\n",cksm ) < 0 ){
					fprintf(stderr,"error in sprintf()\n");
					exit(EXIT_FAILURE);
				}
				report_solution( buffer );

				double d = 1.0;
				boinc_fraction_done(d);

				// average primes/sec speed
				float p_avg = p_sum / p_cnt;

				if(boinc_is_standalone()){
					printf("Workunit complete. Average speed: %.2f Million primes per second.\n", p_avg);
				}
				fprintf(stderr,"Workunit complete. Average speed: %.2f Million primes per second.\n", p_avg);

				if(boinc_is_standalone()){
					printf("Checksum: %016" PRIX64 "\n", cksm);
				}

			}

			checkpoint(bottom,top,stop,cksm,p_sum,p_cnt,write_state_a_next);

			boinc_end_critical_section();

			ckpt_last = ckpt_curr;

			// clear result arrays
			sclEnqueueKernel(pd.hardware, pd.clearresult);
		}


	}


	cleanup(pd);

	boinc_finish(EXIT_SUCCESS);

	return 0; 
} 

