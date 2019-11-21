#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <NTL/ZZ.h>

#include "cycle.h"
#include "rsa.h"

using namespace NTL;

static __inline__ uint64_t rdtsc(void) {
	uint32_t hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return lo | (uint64_t) hi << 32;
}

unsigned long its = 1000;
unsigned long breakIts = 1000;
unsigned long cutoff = 100;
unsigned long benchMulIts = 1000;
unsigned long calibrateIts = 10000; /* Iterations to calculate empirical error in timing */

double median(int n, double x[]) {
    int temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

uint64_t medianInt(int n, uint64_t x[]) {
    int temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

int main() {
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  encryptTotalTime;
  uint64_t encryptTotalCycles;
  uint64_t encryptTotalTicks;
	clock_t  encryptClocks[its];
	uint64_t encryptTimings[its];
	uint64_t encryptTicks[its];

  clock_t  decryptTotalTime;
  uint64_t decryptTotalCycles;
  uint64_t decryptTotalTicks;
	clock_t  decryptClocks[its];
	uint64_t decryptTimings[its];
	uint64_t decryptTicks[its];

  clock_t  decryptCRTTotalTime;
  uint64_t decryptCRTTotalCycles;
  uint64_t decryptCRTTotalTicks;
	clock_t  decryptCRTClocks[its];
	uint64_t decryptCRTTimings[its];
	uint64_t decryptCRTTicks[its];

  clock_t  decryptCRTFaultTotalTime;
  uint64_t decryptCRTFaultTotalCycles;
  uint64_t decryptCRTFaultTotalTicks;
	clock_t  decryptCRTFaultClocks[its];
	uint64_t decryptCRTFaultTimings[its];
	uint64_t decryptCRTFaultTicks[its];

  clock_t  decryptCRTRandExpTotalTime;
  uint64_t decryptCRTRandExpTotalCycles;
  uint64_t decryptCRTRandExpTotalTicks;
	clock_t  decryptCRTRandExpClocks[its];
	uint64_t decryptCRTRandExpTimings[its];
	uint64_t decryptCRTRandExpTicks[its];

  clock_t  decryptCRTBlindTotalTime;
  uint64_t decryptCRTBlindTotalCycles;
  uint64_t decryptCRTBlindTotalTicks;
	clock_t  decryptCRTBlindClocks[its];
	uint64_t decryptCRTBlindTimings[its];
	uint64_t decryptCRTBlindTicks[its];

  clock_t  decryptCRTSecureTotalTime;
  uint64_t decryptCRTSecureTotalCycles;
  uint64_t decryptCRTSecureTotalTicks;
	clock_t  decryptCRTSecureClocks[its];
	uint64_t decryptCRTSecureTimings[its];
	uint64_t decryptCRTSecureTicks[its];

  uint64_t aveTi;
  uint64_t aveti;
  ZZ breakCiphertexts[breakIts];
  uint64_t breakTi[breakIts]; //T_i observed time to decrypt message
  uint64_t breakti[breakIts]; //t_i observed time to perform M_i * M_i^2 mod N
  
  uint64_t decryptModAverage; //Average computation time when there is an additional modulus operation based on d_i = 1
  uint64_t decryptNoModAverage; //Average computation time when there is not an additional modulus operation based on d_i = 1
  uint64_t decryptBreakModAverage; //Average computation time when there is an additional modulus operation based on d_i = 1
  uint64_t decryptBreakNoModAverage; //Average computation time when there is not an additional modulus operation based on d_i = 1
  uint64_t timingError;

  int primeSizes[3] = {512, 1024, 2048}; /* For RSA1024, 2048, 4096 */
  int pubExp = 65537; // Use 2**16 + 1 for efficiency and security

  for (int i = 0; i < 3; i ++)
  {
    
    int primeSize = primeSizes[i]; // 512 bit prime

    /* Construct public and private parameters */
    rsaPrivateParameters_t privParams;
    rsaPublicParameters_t pubParams;

    ZZ m, c, mprime, mprimeCRT, mprimeCRTFault, mprimeCRTRandExp, mprimeCRTBlind, mprimeCRTSecure; //Message m, ciphertext c

	  for (unsigned long iter = 0; iter < its; ++iter) {
      /* Generate random primes of the specified bit length. Use failure rate
       * 2**-30 */
      privParams.p = GenPrime_ZZ(primeSize,30);
      privParams.q = GenPrime_ZZ(primeSize,30);
      privParams.n = privParams.p*privParams.q;
      pubParams.n  = privParams.n;

      /* Generate public and private exponent values */
      privParams.e = ZZ(pubExp);
      pubParams.e  = ZZ(pubExp);

      /* Compute private exponent d such that e*d = 1 mod phi_n, so 
       * d = e^-1 mod phi_n */

      ZZ phi_p, phi_q, phi_n; /* Order of primes p, q, n */
      phi_p = privParams.p-1;
      phi_q = privParams.q-1;
      phi_n = phi_p*phi_q; /* Order of n is the product of its suborders */

      privParams.d = InvMod(pubParams.e,phi_n);

      /* Find private parameters for CRT computations */
      privParams.d_p = privParams.d % phi_p;
      privParams.d_q = privParams.d % phi_q;
      privParams.q_inv = InvMod(privParams.q % privParams.p, privParams.p);

      /* Find a random message to encrypt */
      m = RandomBnd(pubParams.n); //Find a random message in range [0,n-1]
      break;

      /* Encrypt message */
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      c = rsa_encrypt(m,pubParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	encryptTotalTicks = encryptTotalTicks + elapsed(ticks2, ticks1);
	  	encryptTotalCycles += c1 - c0;
	  	encryptTotalTime += t1 - t0;

      encryptTimings[iter] = c1-c0;
      encryptTicks[iter] = elapsed(ticks2,ticks1); 
      encryptClocks[iter] = t1-t0;

      /* Decrypt message using slow method */
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprime = rsa_decrypt(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptTotalTicks = decryptTotalTicks + elapsed(ticks2, ticks1);
	  	decryptTotalCycles += c1 - c0;
	  	decryptTotalTime += t1 - t0;

      decryptTimings[iter] = c1-c0;
      decryptTicks[iter] = elapsed(ticks2,ticks1); 
      decryptClocks[iter] = t1-t0;

      /* Decrypt message using CRT method */
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprimeCRT = rsa_decryptCRT(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptCRTTotalTicks = decryptCRTTotalTicks + elapsed(ticks2, ticks1);
	  	decryptCRTTotalCycles += c1 - c0;
	  	decryptCRTTotalTime += t1 - t0;

      decryptCRTTimings[iter] = c1-c0;
      decryptCRTTicks[iter] = elapsed(ticks2,ticks1); 
      decryptCRTClocks[iter] = t1-t0;

      /* Decrypt message using CRT method and fault validation */
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprimeCRTFault = rsa_decryptCRTFault(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptCRTFaultTotalTicks = decryptCRTFaultTotalTicks + elapsed(ticks2, ticks1);
	  	decryptCRTFaultTotalCycles += c1 - c0;
	  	decryptCRTFaultTotalTime += t1 - t0;

      decryptCRTFaultTimings[iter] = c1-c0;
      decryptCRTFaultTicks[iter] = elapsed(ticks2,ticks1); 
      decryptCRTFaultClocks[iter] = t1-t0;

      /* Decrypt message using CRT method and exponent randomization */
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprimeCRTRandExp = rsa_decryptCRTRandExp(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptCRTRandExpTotalTicks = decryptCRTRandExpTotalTicks + elapsed(ticks2, ticks1);
	  	decryptCRTRandExpTotalCycles += c1 - c0;
	  	decryptCRTRandExpTotalTime += t1 - t0;

      decryptCRTRandExpTimings[iter] = c1-c0;
      decryptCRTRandExpTicks[iter] = elapsed(ticks2,ticks1); 
      decryptCRTRandExpClocks[iter] = t1-t0;

      /* Decrypt message using CRT method and blinding*/
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprimeCRTBlind = rsa_decryptCRTBlind(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptCRTBlindTotalTicks = decryptCRTBlindTotalTicks + elapsed(ticks2, ticks1);
	  	decryptCRTBlindTotalCycles += c1 - c0;
	  	decryptCRTBlindTotalTime += t1 - t0;

      decryptCRTBlindTimings[iter] = c1-c0;
      decryptCRTBlindTicks[iter] = elapsed(ticks2,ticks1); 
      decryptCRTBlindClocks[iter] = t1-t0;

      /* Decrypt message using CRT method with fault validation, exponent
       * randomization, and blinding*/
	  	t0 = clock();
	  	c0 = rdtsc();
	  	ticks1 = getticks();
      mprimeCRTSecure = rsa_decryptCRTSecure(c,privParams);
	  	ticks2 = getticks();
	  	c1 = rdtsc();
	  	t1 = clock();
	  	decryptCRTSecureTotalTicks = decryptCRTSecureTotalTicks + elapsed(ticks2, ticks1);
	  	decryptCRTSecureTotalCycles += c1 - c0;
	  	decryptCRTSecureTotalTime += t1 - t0;

      decryptCRTSecureTimings[iter] = c1-c0;
      decryptCRTSecureTicks[iter] = elapsed(ticks2,ticks1); 
      decryptCRTSecureClocks[iter] = t1-t0;

      /* Ensure message is the same before and after RSA function */
      if (compare(m,mprime) != 0)
      {
        std::cout << "m = " << m << "\n";
        std::cout << "c = " << c << "\n";
        std::cout << "mp= " << mprime << "\n";
        std::cout << "mp= " << rsa_decrypt(c,privParams) << "\n";
        std::cout << "p = " << privParams.p << "\n";
        std::cout << "q = " << privParams.q << "\n";
        std::cout << "d = " << privParams.d << "\n";
        std::cout << "d_p = " << privParams.d_p << "\n";
        std::cout << "d_q = " << privParams.d_q << "\n";
        std::cout << "q_inv = " << privParams.q_inv << "\n";
      }
      assert((compare(m,mprime) == 0));
      assert((compare(m,mprimeCRT) == 0));
      assert((compare(m,mprimeCRTFault) == 0));
      assert((compare(m,mprimeCRTRandExp) == 0));
      assert((compare(m,mprimeCRTBlind) == 0));
      assert((compare(m,mprimeCRTSecure) == 0));

	  }
	  printf("Encrypt Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) encryptTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) encryptTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * encryptTotalTime / CLOCKS_PER_SEC / its);
	  printf("Encrypt Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,encryptTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,encryptTicks));

	  printf("Decrypt Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptTotalTime / CLOCKS_PER_SEC / its);
	  printf("Decrypt Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptTicks));

	  printf("DecryptCRT Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptCRTTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptCRTTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptCRTTotalTime / CLOCKS_PER_SEC / its);
	  printf("DecryptCRT Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptCRTTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptCRTTicks));

	  printf("DecryptCRTFault Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptCRTFaultTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptCRTFaultTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptCRTFaultTotalTime / CLOCKS_PER_SEC / its);
	  printf("DecryptCRTFault Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptCRTFaultTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptCRTFaultTicks));

	  printf("DecryptCRTRandExp Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptCRTRandExpTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptCRTRandExpTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptCRTRandExpTotalTime / CLOCKS_PER_SEC / its);
	  printf("DecryptCRTRandExp Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptCRTRandExpTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptCRTRandExpTicks));

	  printf("DecryptCRTBlind Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptCRTBlindTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptCRTBlindTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptCRTBlindTotalTime / CLOCKS_PER_SEC / its);
	  printf("DecryptCRTBlind Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptCRTBlindTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptCRTBlindTicks));

	  printf("DecryptCRTSecure Averages for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) decryptCRTSecureTotalCycles / its);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) decryptCRTSecureTotalTicks / its);
	  printf("wall clock time: %.31f\n", 1000. * decryptCRTSecureTotalTime / CLOCKS_PER_SEC / its);
	  printf("DecryptCRTSecure Medians for %lu iterations of RSA%d\n", its, 2*primeSize);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(its,decryptCRTSecureTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(its,decryptCRTSecureTicks));

    encryptTotalTicks = 0;
  	encryptTotalCycles = 0;
  	encryptTotalTime = 0;
    decryptTotalTicks = 0;
  	decryptTotalCycles = 0;
  	decryptTotalTime = 0;
    decryptCRTTotalTicks = 0;
  	decryptCRTTotalCycles = 0;
  	decryptCRTTotalTime = 0;
    decryptCRTFaultTotalTicks = 0;
  	decryptCRTFaultTotalCycles = 0;
  	decryptCRTFaultTotalTime = 0;
    decryptCRTRandExpTotalTicks = 0;
  	decryptCRTRandExpTotalCycles = 0;
  	decryptCRTRandExpTotalTime = 0;
    decryptCRTBlindTotalTicks = 0;
  	decryptCRTBlindTotalCycles = 0;
  	decryptCRTBlindTotalTime = 0;
    decryptCRTSecureTotalTicks = 0;
  	decryptCRTSecureTotalCycles = 0;
  	decryptCRTSecureTotalTime = 0;

    /* RSA fault attack */
    /* Attempt to break RSA when there is a random fault on one of the bits
     * of a signature. For instance, for RSA1024, we try flipping a bit
     * of a prime's component in the range[0,511] */
    for (int bitToCorrupt = 0; bitToCorrupt < primeSize; bitToCorrupt++)
    {
      break;
      /* Find a random message to encrypt. This attack won't work for m = 0,1 */
      do {
        m = RandomBnd(pubParams.n); //Find a random message in range [2,n-1]
      } while (compare(m,ZZ(1)) <= 0);

      /* Encrypt the message */
      c = rsa_encrypt(m,pubParams);

      /* Force the RSA decrypt function to produce a faulty decryption */
      ZZ mCorrupt = rsa_decryptCRT_fault(c, privParams, bitToCorrupt);

      /* gcd(M - E^e, N) = q */
      ZZ faultyEncrypt = PowerMod(mCorrupt % pubParams.n, pubParams.e, pubParams.n);
      ZZ diff = (m - faultyEncrypt) % pubParams.n;
      ZZ qSuspect = GCD(diff,pubParams.n);

      /* Check if N | qSuspect. If so RSA private parameters have been
       * retrieved as p*q = N */
      assert(divide(pubParams.n,qSuspect) == 1);
    }

    /* RSA timing attack */
    aveTi = 0; 
    for (int i = 0; i < breakIts; i++)
    {
      /* Find a random message to decrypt */
      c = RandomBnd(pubParams.n); //Find a random message in range [0,n-1]

      /* Encrypt message */
	  	c0 = rdtsc();
      m = rsa_decryptTimingWeakness(c,privParams);
	  	c1 = rdtsc();

      aveTi += c1-c0;
      breakTi[i] = c1-c0;
      breakCiphertexts[i] = c;

    }
    aveTi = aveTi / breakIts;

    /* Start with d = 1 mod 2 */
    ZZ reconstructedD = ZZ(1);
    aveti = 0;
    for (int bitPos = 1; bitPos < primeSize; bitPos++)
    {
      /* Iteratively break RSA */
      std::cout << "reconstructedD " << reconstructedD << "\n"; 
      int powerOfTwo = 1<<bitPos;
      std::cout << "2^x " << powerOfTwo << "\n"; 
      for (int i = 0; i < breakIts; i++)
      {
        /* M_i */
        ZZ valueToTest = PowerMod(breakCiphertexts[i], reconstructedD, pubParams.n);
        /* M_i ^ bitPosition */
        ZZ newValue = PowerMod(breakCiphertexts[i], powerOfTwo, pubParams.n);
        
        /* Time to perform mult and reduction x times */
	    	c0 = rdtsc();
        //for (int j = 0; j < benchMulIts; j++)
        //{
          ZZ tmp = (valueToTest * newValue) % pubParams.n;
        //}
	    	c1 = rdtsc();

        aveti += (c1-c0);
        breakti[i] = (c1-c0); 

      }

      aveti = aveti / breakIts;

      int bothGreaterThanExpect = 0;
      int oneGreaterThanExpect = 0;
      for (int i = 0; i < breakIts; i++)
      {
        if ((breakti[i] > aveti) && (breakTi[i] > aveTi))
        { 
          bothGreaterThanExpect++;
        }
        else if (breakti[i] > aveti)
        {
          oneGreaterThanExpect++;
        }
      }

      std::cout << "Both greater than expect = " << bothGreaterThanExpect << "\n";
      std::cout << "One greater than expect  = " << oneGreaterThanExpect << "\n";

      if (bothGreaterThanExpect > cutoff)
      {
        std::cout << "Attacker guesses d_i = 1" << "\n";
        reconstructedD = reconstructedD + ZZ(1<<bitPos);
      }
      else
      {
        std::cout << "Attacker guesses d_i = 0" << "\n";
      }

      std::cout << "d_i = " << bit(privParams.d, bitPos) << "\n";

      reconstructedD = reconstructedD + ZZ(bit(privParams.d, bitPos)<<bitPos);
    }

    /* Find an error term e in timing based on implementation */
    timingError = 0;
    decryptModAverage = 0;
    decryptNoModAverage = 0;

    /* Pick a random d to do calibration with */
    rsaPrivateParameters_t privParamsMock;
    privParamsMock.n = pubParams.n;
    privParamsMock.d = RandomBnd(pubParams.n);

    /* Set bit d_1 = 1 */
    if (!IsOdd(privParamsMock.d >> 1))
      privParamsMock.d = privParamsMock.d + 2;

    ZZ Y, Z;
    for (int i = 0; i < calibrateIts; i++)
    {
      break;
      /* Find a message Y such that Y^3 < N. This will not cause only one
       * modulus rather than 2 */
      do {
        Y = RandomBits_ZZ(primeSize/2);
        //Y = SqrRoot(pubParams.n) - RandomBits_ZZ(32) + RandomBits_ZZ(32);
        //std::cout << "Y_i = " << Y << "\n";
        //std::cout << "n   = " << pubParams.n << "\n";
        //std::cout << "Y_i2= " << ((Y*Y)%pubParams.n) << "\n";
        //std::cout << "Y_i3= " << ((Y*Y)%pubParams.n)*Y << "\n";
      } while (compare(pubParams.n, ((Y*Y)%pubParams.n)*Y) != 1);

      /* Find a message Z such that Z^2 < N < Z^3. This will force 2
       * modulus operations */
      do {
        Z = RandomBnd(pubParams.n);
      } while ((compare(pubParams.n, (Z*Z)%pubParams.n) != 1) || (compare(pubParams.n, ((Z*Z)%pubParams.n)*Z) != -1));
        //std::cout << "Z_i" << Z << "\n";


	  	c0 = rdtsc();
      c = rsa_decryptTimingWeakness(Y,privParamsMock);
	  	c1 = rdtsc();

      decryptNoModAverage += c1 - c0;

      //std::cout << "decryptNoMod_i" << c1-c0 << "\n";

	  	c0 = rdtsc();
      c = rsa_decryptTimingWeakness(Z,privParamsMock);
	  	c1 = rdtsc();

	  	decryptModAverage += c1 - c0;
     // std::cout << "decryptMod_i  " << c1-c0 << "\n";

    }

    decryptModAverage = decryptModAverage / calibrateIts;
    decryptNoModAverage = decryptNoModAverage / calibrateIts;

    /* Use margin 3/4 */
    timingError = (3*(decryptModAverage - decryptNoModAverage))/4;

    std::cout << "  Mod Ave   = " << decryptModAverage << "\n";
    std::cout << "NoMod Ave   = " << decryptNoModAverage << "\n";
    std::cout << "timingError = " << timingError << "\n";

    decryptBreakModAverage = 0;
    decryptBreakNoModAverage = 0;

    /* Attempt to crack bit d_1 in actual implementation */
    for (int i = 0; i < breakIts; i++)
    {
      break;
      /* Find a message Y such that Y^3 < N. This will not cause only one
       * modulus rather than 2 */
      do {
        Y = RandomBits_ZZ(primeSize/2);
        //Y = SqrRoot(pubParams.n) - RandomBits_ZZ(32) + RandomBits_ZZ(32);
        //std::cout << "Y_i = " << Y << "\n";
        //std::cout << "n   = " << pubParams.n << "\n";
        //std::cout << "Y_i2= " << ((Y*Y)%pubParams.n) << "\n";
        //std::cout << "Y_i3= " << ((Y*Y)%pubParams.n)*Y << "\n";
      } while (compare(pubParams.n, ((Y*Y)%pubParams.n)*Y) != 1);

      /* Find a message Z such that Z^2 < N < Z^3. This will force 2
       * modulus operations */
      do {
        Z = RandomBnd(pubParams.n);
      } while ((compare(pubParams.n, (Z*Z)%pubParams.n) != 1) || (compare(pubParams.n, ((Z*Z)%pubParams.n)*Z) != -1));
        //std::cout << "Z_i" << Z << "\n";


	  	c0 = rdtsc();
      c = rsa_decryptTimingWeakness(Y,privParams);
	  	c1 = rdtsc();

      decryptBreakNoModAverage += c1 - c0;

      //std::cout << "decryptNoMod_i" << c1-c0 << "\n";

	  	c0 = rdtsc();
      c = rsa_decryptTimingWeakness(Z,privParams);
	  	c1 = rdtsc();

	  	decryptBreakModAverage += c1 - c0;
     // std::cout << "decryptMod_i  " << c1-c0 << "\n";

    }

    decryptBreakModAverage = decryptBreakModAverage / breakIts;
    decryptBreakNoModAverage = decryptBreakNoModAverage / breakIts;

    std::cout << "  Mod Aveb = " << decryptBreakModAverage << "\n";
    std::cout << "NoMod Aveb = " << decryptBreakNoModAverage << "\n";

    /* If the decryption time with mod is larger than noMod and timingError,
     * then bit d_1 = 1 */
    if (decryptBreakModAverage > (decryptBreakNoModAverage + timingError))
    {
      std::cout << "Attacker guesses d_1 = 1" << "\n";
    }
    else
    {
      std::cout << "Attacker guesses d_1 = 0" << "\n";
    }

    std::cout << "d_1 = " << bit(privParams.d, 1) << "\n";




  }
}

