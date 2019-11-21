#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <NTL/ZZ.h>

#include "cycle.h"
#include "jpake.h"

using namespace NTL;

static __inline__ uint64_t rdtsc(void) {
	uint32_t hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return lo | (uint64_t) hi << 32;
}

unsigned long its = 1000;
unsigned long paramIts = 100;
unsigned long breakIts = 10;
unsigned long cutoff = 100;
unsigned long benchMulIts = 1000;
unsigned long calibrateIts = 10000; /* Iterations to calculate empirical error in timing */

unsigned long ORDERLEN = 256;
unsigned long PWDLEN = 32;

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

void timeDHParamGen(int primeLen)
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  paramGenTotalTime;
  uint64_t paramGenTotalCycles;
  uint64_t paramGenTotalTicks;
	clock_t  paramGenClocks[paramIts];
	uint64_t paramGenTimings[paramIts];
	uint64_t paramGenTicks[paramIts];

  uint64_t paramGenPTotalLoops;
  uint64_t paramGenGTotalLoops;
  
  int orderLens[5] = {32, 64, 128, 256, 512}; 
  for (int i = 0; i < 5; i++)
  {
    int curOrderLen = orderLens[i];
    if (curOrderLen >= primeLen)
    {
      continue;
    }

    paramGenPTotalLoops = 0;
    paramGenGTotalLoops = 0;
    paramGenTotalTime = 0;
    paramGenTotalCycles = 0;
    paramGenTotalTicks = 0;

    for (int iter = 0; iter < paramIts; iter++)
    {
      t0 = clock();
      c0 = rdtsc();
      ticks1 = getticks();
      /* Time to generate params */
      ZZ q = GenPrime_ZZ(curOrderLen, 10);
      ZZ randMult = RandomBits_ZZ(primeLen - curOrderLen);
      MakeOdd(randMult);
      SetBit(randMult, primeLen - curOrderLen);

      ZZ p = q *randMult + ZZ(1);
      paramGenPTotalLoops++;
      while (!ProbPrime(p,10))
      {
        p = p + q;
        paramGenPTotalLoops++;
      }

      ZZ g = ZZ(1);
      ZZ e = (p - ZZ(1))/q;
      ZZ h = ZZ(2);
      while (IsOne(g))
      {
        g = PowerMod(h,e,p);
        paramGenGTotalLoops++;
      }
      ticks2 = getticks();
      c1 = rdtsc();
      t1 = clock();
      paramGenTotalTicks = paramGenTotalTicks + elapsed(ticks2, ticks1);
      paramGenTotalCycles += c1 - c0;
      paramGenTotalTime += t1 - t0;
      paramGenTimings[iter] = c1-c0;
      paramGenTicks[iter] = elapsed(ticks2,ticks1); 
      paramGenClocks[iter] = t1-t0;
    }
	  printf("ParamGen Averages for %lu iterations of DH%d with order %d\n", paramIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) paramGenTotalCycles / paramIts);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) paramGenTotalTicks / paramIts);
	  printf("wall clock time: %.31f\n", 1000. * paramGenTotalTime / CLOCKS_PER_SEC / paramIts);
	  printf("Iterations for p %lu\n", paramGenPTotalLoops / paramIts);
	  printf("Iterations for g %lu\n", paramGenGTotalLoops / paramIts);
	  printf("ParamGen Medians for %lu iterations of DH%d with order %d\n", paramIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(paramIts,paramGenTimings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(paramIts,paramGenTicks));
  } 
}

void timeDH(int primeLen, int orderLen)
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  roundOneTotalTime;
  uint64_t roundOneTotalCycles = 0;
  uint64_t roundOneTotalTicks = 0;
	clock_t  roundOneClocks[its];
	uint64_t roundOneTimings[its];
	uint64_t roundOneTicks[its];

  clock_t  roundTwoTotalTime;
  uint64_t roundTwoTotalCycles = 0;
  uint64_t roundTwoTotalTicks = 0;
	clock_t  roundTwoClocks[its];
	uint64_t roundTwoTimings[its];
	uint64_t roundTwoTicks[its];

  /* Time round 1 */
  for (int iter = 0; iter < its; iter++)
  {

    /* Find parameters */
    dhParameters_t params = dh_findParams(primeLen, orderLen);

    /* Perform Alice Round 1 and time it */
    t0 = clock();
    c0 = rdtsc();
    ticks1 = getticks();
    dhKeyPair_t aliceDH1 = dh_genKeyPair(params.g, params);
    ticks2 = getticks();
    c1 = rdtsc();
    t1 = clock();
    roundOneTotalTicks = roundOneTotalTicks + elapsed(ticks2, ticks1);
    roundOneTotalCycles += c1 - c0;
    roundOneTotalTime += t1 - t0;
    roundOneTimings[iter] = c1-c0;
    roundOneTicks[iter] = elapsed(ticks2,ticks1); 
    roundOneClocks[iter] = t1-t0;

    dhKeyPair_t bobDH1 = dh_genKeyPair(params.g, params);

    t0 = clock();
    c0 = rdtsc();
    ticks1 = getticks();
    ZZ A = dh_exp(bobDH1.pubKey, aliceDH1.privKey, params);
    ticks2 = getticks();
    c1 = rdtsc();
    t1 = clock();
    roundTwoTotalTicks = roundTwoTotalTicks + elapsed(ticks2, ticks1);
    roundTwoTotalCycles += c1 - c0;
    roundTwoTotalTime += t1 - t0;
    roundTwoTimings[iter] = c1-c0;
    roundTwoTicks[iter] = elapsed(ticks2,ticks1); 
    roundTwoClocks[iter] = t1-t0;

    ZZ B = dh_exp(aliceDH1.pubKey, bobDH1.privKey, params);

    assert((compare(A,B) == 0));
    assert((compare(A,0) == 1));
  }

	printf("R1 Averages for %lu iterations of DH%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",(uint64_t) roundOneTotalCycles / its);
	printf("clock cycles: %lu (getticks)\n",(uint64_t) roundOneTotalTicks / its);
	printf("wall clock time: %.31f\n", 1000. * roundOneTotalTime / CLOCKS_PER_SEC / its);
	printf("R1 Medians for %lu iterations of DH%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundOneTimings));
	printf("clock cycles: %lu (getticks)\n\n",medianInt(its,roundOneTicks));

	printf("R2 Averages for %lu iterations of DH%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",(uint64_t) roundTwoTotalCycles / its);
	printf("clock cycles: %lu (getticks)\n",(uint64_t) roundTwoTotalTicks / its);
	printf("wall clock time: %.31f\n", 1000. * roundTwoTotalTime / CLOCKS_PER_SEC / its);
	printf("R2 Medians for %lu iterations of DH%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundTwoTimings));
	printf("clock cycles: %lu (getticks)\n\n",medianInt(its,roundTwoTicks));

  return;
}


void timeJPAKE(int primeLen, int orderLen, int pwdLen)
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  roundOneTotalTime;
  uint64_t roundOneTotalCycles = 0;
  uint64_t roundOneTotalTicks = 0;
	clock_t  roundOneClocks[its];
	uint64_t roundOneTimings[its];
	uint64_t roundOneTicks[its];

  clock_t  roundTwoTotalTime;
  uint64_t roundTwoTotalCycles = 0;
  uint64_t roundTwoTotalTicks = 0;
	clock_t  roundTwoClocks[its];
	uint64_t roundTwoTimings[its];
	uint64_t roundTwoTicks[its];

  clock_t  roundThreeTotalTime;
  uint64_t roundThreeTotalCycles = 0;
  uint64_t roundThreeTotalTicks = 0;
	clock_t  roundThreeClocks[its];
	uint64_t roundThreeTimings[its];
	uint64_t roundThreeTicks[its];

  /* Time round 1 */
  for (int iter = 0; iter < its; iter++)
  {

    /* Find a random 32-bit password */
    ZZ pwd = RandomBits_ZZ(pwdLen);

    /* Find parameters */
    dhParameters_t params = dh_findParams(primeLen, orderLen);

    /* Perform Alice Round 1 and time it */
    t0 = clock();
    c0 = rdtsc();
    ticks1 = getticks();
    jpakeRound1Struct_t alice1 = jpake_roundOne(params);
    ticks2 = getticks();
    c1 = rdtsc();
    t1 = clock();
    roundOneTotalTicks = roundOneTotalTicks + elapsed(ticks2, ticks1);
    roundOneTotalCycles += c1 - c0;
    roundOneTotalTime += t1 - t0;
    roundOneTimings[iter] = c1-c0;
    roundOneTicks[iter] = elapsed(ticks2,ticks1); 
    roundOneClocks[iter] = t1-t0;

    jpakeRound1Struct_t bob1   = jpake_roundOne(params);

    t0 = clock();
    c0 = rdtsc();
    ticks1 = getticks();
    jpakeRound2Struct_t alice2 = jpake_roundTwo(alice1, bob1, pwd, params);
    ticks2 = getticks();
    c1 = rdtsc();
    t1 = clock();
    roundTwoTotalTicks = roundTwoTotalTicks + elapsed(ticks2, ticks1);
    roundTwoTotalCycles += c1 - c0;
    roundTwoTotalTime += t1 - t0;
    roundTwoTimings[iter] = c1-c0;
    roundTwoTicks[iter] = elapsed(ticks2,ticks1); 
    roundTwoClocks[iter] = t1-t0;

    jpakeRound2Struct_t bob2   = jpake_roundTwo(bob1, alice1, pwd, params);

    t0 = clock();
    c0 = rdtsc();
    ticks1 = getticks();
    ZZ alice3 = jpake_roundThree(alice1, bob1, alice2, bob2, pwd, params);
    ticks2 = getticks();
    c1 = rdtsc();
    t1 = clock();
    roundThreeTotalTicks = roundThreeTotalTicks + elapsed(ticks2, ticks1);
    roundThreeTotalCycles += c1 - c0;
    roundThreeTotalTime += t1 - t0;
    roundThreeTimings[iter] = c1-c0;
    roundThreeTicks[iter] = elapsed(ticks2,ticks1); 
    roundThreeClocks[iter] = t1-t0;
    
    ZZ bob3 =   jpake_roundThree(bob1, alice1, bob2, alice2, pwd, params);

    assert((compare(alice3,bob3) == 0));
    assert((compare(alice3,0) == 1));
  }

	printf("R1 Averages for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",(uint64_t) roundOneTotalCycles / its);
	printf("clock cycles: %lu (getticks)\n",(uint64_t) roundOneTotalTicks / its);
	printf("wall clock time: %.31f\n", 1000. * roundOneTotalTime / CLOCKS_PER_SEC / its);
	printf("R1 Medians for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundOneTimings));
	printf("clock cycles: %lu (getticks)\n\n",medianInt(its,roundOneTicks));

	printf("R2 Averages for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",(uint64_t) roundTwoTotalCycles / its);
	printf("clock cycles: %lu (getticks)\n",(uint64_t) roundTwoTotalTicks / its);
	printf("wall clock time: %.31f\n", 1000. * roundTwoTotalTime / CLOCKS_PER_SEC / its);
	printf("R2 Medians for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundTwoTimings));
	printf("clock cycles: %lu (getticks)\n\n",medianInt(its,roundTwoTicks));

	printf("R3 Averages for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",(uint64_t) roundThreeTotalCycles / its);
	printf("clock cycles: %lu (getticks)\n",(uint64_t) roundThreeTotalTicks / its);
	printf("wall clock time: %.31f\n", 1000. * roundThreeTotalTime / CLOCKS_PER_SEC / its);
	printf("R3 Medians for %lu iterations of JPAKE%d\n", its, primeLen);
	printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundThreeTimings));
	printf("clock cycles: %lu (getticks)\n\n",medianInt(its,roundThreeTicks));

  return;

}

void timeBruteForce(int primeLen)
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  totalTime;
  uint64_t totalCycles;
  uint64_t totalTicks;
	clock_t  clocks[paramIts];
	uint64_t timings[paramIts];
	uint64_t ticks[paramIts];

  uint64_t totalGuesses;
  
  int orderLens[3] = {8, 16, 24}; 
  for (int i = 0; i < 3; i++)
  {
    int curOrderLen = orderLens[i];

    totalGuesses = 0;
    totalTime = 0;
    totalCycles = 0;
    totalTicks = 0;

    for (int iter = 0; iter < breakIts; iter++)
    {
      /* Get parameters set */
      dhParameters_t params = dh_findParams(primeLen, curOrderLen);
      /* Get keyPair */
      dhKeyPair_t keyPair = dh_genKeyPair(params.g, params);

      t0 = clock();
      c0 = rdtsc();
      ticks1 = getticks();

      assert(dl_bruteForce(keyPair, params) == 1);

      ticks2 = getticks();
      c1 = rdtsc();
      t1 = clock();
      totalTicks = totalTicks + elapsed(ticks2, ticks1);
      totalCycles += c1 - c0;
      totalTime += t1 - t0;
      timings[iter] = c1-c0;
      ticks[iter] = elapsed(ticks2,ticks1); 
      clocks[iter] = t1-t0;
    }
    printf("\n");
	  printf("BruteForce Averages for %lu iterations of DH%d with order %d\n", breakIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) totalCycles / breakIts);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) totalTicks / breakIts);
	  printf("wall clock time: %.31f\n", 1000. * totalTime / CLOCKS_PER_SEC / breakIts);
	  printf("Iterations for key %lu\n", totalGuesses / breakIts);
	  printf("BruteForce Medians for %lu iterations of DH%d with order %d\n", breakIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(breakIts,timings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(breakIts,ticks));
  } 
}

void timePollardRho(int primeLen)
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  totalTime;
  uint64_t totalCycles;
  uint64_t totalTicks;
	clock_t  clocks[paramIts];
	uint64_t timings[paramIts];
	uint64_t ticks[paramIts];

  uint64_t totalGuesses;
  
  int orderLens[6] = {8, 16, 24, 32, 40, 48}; 
  for (int i = 0; i < 6; i++)
  {
    int curOrderLen = orderLens[i];

    totalGuesses = 0;
    totalTime = 0;
    totalCycles = 0;
    totalTicks = 0;

    for (int iter = 0; iter < breakIts; iter++)
    {
      /* Get parameters set */
      dhParameters_t params = dh_findParams(primeLen, curOrderLen);
      /* Get keyPair */
      dhKeyPair_t keyPair = dh_genKeyPair(params.g, params);

      t0 = clock();
      c0 = rdtsc();
      ticks1 = getticks();

      //assert(dl_pollardRho(keyPair, params) == 1);
      dl_pollardRho(keyPair, params);

      ticks2 = getticks();
      c1 = rdtsc();
      t1 = clock();
      totalTicks = totalTicks + elapsed(ticks2, ticks1);
      totalCycles += c1 - c0;
      totalTime += t1 - t0;
      timings[iter] = c1-c0;
      ticks[iter] = elapsed(ticks2,ticks1); 
      clocks[iter] = t1-t0;
    }
    printf("\n");
	  printf("PollardRho Averages for %lu iterations of DH%d with order %d\n", breakIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",(uint64_t) totalCycles / breakIts);
	  printf("clock cycles: %lu (getticks)\n",(uint64_t) totalTicks / breakIts);
	  printf("wall clock time: %.31f\n", 1000. * totalTime / CLOCKS_PER_SEC / breakIts);
	  printf("Iterations for key %lu\n", totalGuesses / breakIts);
	  printf("PollardRho Medians for %lu iterations of DH%d with order %d\n", breakIts, primeLen, curOrderLen);
	  printf("clock cycles: %lu (rdtsc)\n",medianInt(breakIts,timings));
	  printf("clock cycles: %lu (getticks)\n\n",medianInt(breakIts,ticks));
  } 
}

void timeJPAKEPwd()
{
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  roundOneTotalTime;
  uint64_t roundOneTotalCycles = 0;
  uint64_t roundOneTotalTicks = 0;
	clock_t  roundOneClocks[its];
	uint64_t roundOneTimings[its];
	uint64_t roundOneTicks[its];

  int orderLens[3] = {128, 256, 512};
  /* Time round 1 */
  for (int i = 0; i < 3; i++)
  {
    int curOrderLen = orderLens[i];

    for (int pwdLen = 64; pwdLen < 256+1; pwdLen = pwdLen + 64)
    {
      roundOneTotalTime = 0;
      roundOneTotalCycles = 0;
      roundOneTotalTicks = 0;
      for (int iter = 0; iter < its; iter++)
      {

        /* Find a random 32-bit password */
        ZZ pwd = RandomBits_ZZ(pwdLen);
        ZZ randOrder = GenPrime_ZZ(curOrderLen,10); 
        ZZ secretKey = RandomBnd(randOrder);

        /* Perform password computation and time it */
        t0 = clock();
        c0 = rdtsc();
        ticks1 = getticks();
        ZZ secretExp = MulMod(secretKey,pwd,randOrder);
        ticks2 = getticks();
        c1 = rdtsc();
        t1 = clock();
        roundOneTotalTicks = roundOneTotalTicks + elapsed(ticks2, ticks1);
        roundOneTotalCycles += c1 - c0;
        roundOneTotalTime += t1 - t0;
        roundOneTimings[iter] = c1-c0;
        roundOneTicks[iter] = elapsed(ticks2,ticks1); 
        roundOneClocks[iter] = t1-t0;

      }
	    printf("PWDGen Medians for %lu iterations of order %d with pwdLen %d\n", its, curOrderLen,pwdLen);
	    printf("clock cycles: %lu (rdtsc)\n",medianInt(its,roundOneTimings));
    }
  }
  return;
}


int main() {
	clock_t t0, t1;
	uint64_t c0, c1;
	ticks ticks1, ticks2;

  clock_t  paramGenTotalTime;
  uint64_t paramGenTotalCycles;
  uint64_t paramGenTotalTicks;
	clock_t  paramGenClocks[its];
	uint64_t paramGenTimings[its];
	uint64_t paramGenTicks[its];

  clock_t  expTotalTime;
  uint64_t expTotalCycles;
  uint64_t expTotalTicks;
	clock_t  expClocks[its];
	uint64_t expTimings[its];
	uint64_t expTicks[its];

  int primeSizes[5] = {512, 1024, 2048, 3072, 4096}; /* For DH512, 1024, 2048, 3072, 4096 */

  for (int i = 0; i < 5; i ++)
  {
    
    int primeLen = primeSizes[i]; // 512 bit prime

    /* DH Parameter generation */
    //timeDHParamGen(primeLen);

    /* Time DH operations */
    //timeDH(primeLen, ORDERLEN);

    /* Time JPAKE operations */
    // timeJPAKE(primeLen, ORDERLEN, PWDLEN);

    /* Time Brute Force Discrete log */
    // timeBruteForce(primeLen);

    /* Time Pohlig-Hellman */
    /* not a feasible attack for large prime orders as we do in our JPAKE */

    /* Time Pollard Rho */
    //timePollardRho(primeLen);
    
  }
  /* Time JPAKE pwd computation */
  timeJPAKEPwd();

  return 1;
}
