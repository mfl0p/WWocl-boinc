# WWocl-boinc

WWocl by Bryan Little and Yves Gallot

A BOINC-enabled OpenCL Wieferich and WallSunSun prime search program.

Based on wwwwcl by Mark Rodenkirch and primesieve by Kim Walisch.

All calculations on GPU except some error checking.

A checksum will be printed at the end of the results-WWocl.txt upon completion of the test range.
This checksum will allow comparing results between different hardware.


## How it works

1) A search range is given on the command line.

2) Optimal gpu kernel work size is calculated using profiling.

3) Fast cpu-side verifications of the gpu kernels are next.  A small range of numbers are tested and
   compared to the same tests done on the cpu.  If the cpu detects the kernels are producing invalid
   results an error is generated.  This should catch any kernels that are broken by the video driver's
   compiler.  In addition, the PRP generator algorithm is verified for proper operation by comparing
   results to Kim Walisch's primesieve program.

4) The gpu takes a small chunk of the search range and generates a list of base 2 probable primes.
   These are "industrial grade primes" requiring ~7 times fewer calculations than testing for primality.
   This way we can quickly find candidate primes to Wieferich / WallSunSun test and remove any false 
   positives later on the cpu.  The kernel uses a combined wheel + bitsieve + base 2 prp test.
   This is a fast approach compared to the many memory accesses required for a simple sieve of Eratosthenes.

5) The list of prps are tested for Wieferich primes.  A checksum is recorded for each test.

6) The list of prps are tested for WallSunSun primes.  A checksum is recorded for each test.

7) This repeats until the checkpoint range has been reached, where results are moved to the cpu.
   Results are verified for primality and recorded to the solutions file.

8) The checksum is printed in the result file at the end of the search range.


## Running the program stand-alone

Usage: WWocl START END

where 127 <= START < END < 2^64

Example formats:
WWocl.exe 10000000 20000000
WWocl.exe 1e7 2e7

The program will default to GPU 0 unless an init_data.xml is in the same directory with the format:

<app_init_data>
<gpu_type>NVIDIA</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>

or

<app_init_data>
<gpu_type>ATI</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>


## Related Links

* [primesieve](https://github.com/kimwalisch/primesieve)


