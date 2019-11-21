
#include "jpake.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

dhParameters_t dh_findParams(int primeLen, int orderLen)
{
  dhParameters_t params;
  /* Find order of subgroup first */
  params.q = GenPrime_ZZ(orderLen,10);

  /* Find a random multiple of the order */
  ZZ randMult = RandomBits_ZZ(primeLen - orderLen);
  /* Set MSB and LSB to be '1' so that the mult result is odd */
  MakeOdd(randMult);
  SetBit(randMult, primeLen - orderLen);

  params.p = params.q * randMult + 1;
  int cnt = 1;
  /* Increment p until p is a prime */
  while (!ProbPrime(params.p,10))
  {
    params.p = params.p + params.q;
    cnt++;
  }

  /* Find a generator value g that is of the prime order q */
  params.g = 1;
  cnt = 1;
  /* Cancel out orders on other subgroups */
  ZZ e = (params.p - ZZ(1))/params.q;
  ZZ h = ZZ(2);
  while (IsOne(params.g))
  {
    params.g = PowerMod(h,e,params.p);
    cnt++;
  }

  return params;
}

ZZ dh_exp(ZZ base, ZZ exp, dhParameters_t params)
{
  return PowerMod(base % params.p, exp % params.q, params.p);
}

dhKeyPair_t dh_genKeyPair(ZZ gen, dhParameters_t params)
{
  dhKeyPair_t keyPair;
  /* Generate a private key modulo the subgroup */
  RandomBnd(keyPair.privKey, params.q);
  /* Perform the exponentiation to get the public key */
  PowerMod(keyPair.pubKey, gen, keyPair.privKey, params.p);
  return keyPair;
}

jpakePubKey_t jpake_genZKP(dhKeyPair_t keyPair, ZZ gen, dhParameters_t params)
{
  jpakePubKey_t pubKeyOut;
  dhKeyPair_t vKeyPair = dh_genKeyPair(gen, params);
  pubKeyOut.pubKey = keyPair.pubKey;
  pubKeyOut.zkpV = vKeyPair.pubKey;
  /* Use a simple value for digest */
  ZZ digest = pubKeyOut.pubKey - ZZ(1);
  pubKeyOut.zkpR = SubMod(vKeyPair.privKey, MulMod(keyPair.privKey, digest, params.q), params.q);

  return pubKeyOut;
}

int jpake_verifyZKP(jpakePubKey_t zkpIn, ZZ gen, dhParameters_t params)
{
  /* Use a simple value for digest */
  ZZ digest = zkpIn.pubKey - ZZ(1);
  ZZ exp1 = dh_exp(gen, zkpIn.zkpR, params);
  ZZ exp2 = dh_exp(zkpIn.pubKey, digest, params);
  ZZ VPrime = MulMod(exp1, exp2, params.p);

  return (compare(zkpIn.zkpV,VPrime) == 0); 
}

jpakeRound1Struct_t jpake_roundOne(dhParameters_t params)
{
  jpakeRound1Struct_t ret;
  ret.keyPair1 = dh_genKeyPair(params.g, params);
  ret.keyPair2 = dh_genKeyPair(params.g, params);
  ret.pubKey1 = jpake_genZKP(ret.keyPair1, params.g, params);
  ret.pubKey2 = jpake_genZKP(ret.keyPair2, params.g, params);

  return ret;
}

jpakeRound2Struct_t jpake_roundTwo(jpakeRound1Struct_t myRound1, jpakeRound1Struct_t theirRound1, ZZ pwd, dhParameters_t params)
{
  jpakeRound2Struct_t ret;
  /* Verify ZKPs */
  if (!(jpake_verifyZKP(theirRound1.pubKey1, params.g, params) ||
      !(jpake_verifyZKP(theirRound1.pubKey2, params.g, params))))
  {
    printf("Invalid Round 1 jpake ZKP!\n");
    return ret;
  }

  /* Compute new generator point for Round 2 */ 
  ZZ base = MulMod(MulMod(myRound1.pubKey1.pubKey, theirRound1.pubKey1.pubKey, params.p), theirRound1.pubKey2.pubKey, params.p);

  /* Compute new private key mixed with shared password */
  ret.combinedKeyPair.privKey = MulMod(myRound1.keyPair2.privKey, pwd, params.q);
  ret.combinedKeyPair.pubKey = dh_exp(base, ret.combinedKeyPair.privKey, params);

  ret.combinedPubKey = jpake_genZKP(ret.combinedKeyPair, base, params);

  return ret;
}

ZZ jpake_roundThree(jpakeRound1Struct_t myRound1, jpakeRound1Struct_t theirRound1, jpakeRound2Struct_t myRound2, jpakeRound2Struct_t theirRound2, ZZ pwd, dhParameters_t params)
{
  /* find base for ZKP */
  ZZ baseGen = MulMod(MulMod(theirRound1.pubKey1.pubKey, myRound1.pubKey1.pubKey, params.p), myRound1.pubKey2.pubKey, params.p);

  /* Verify ZKP */
  if (!jpake_verifyZKP(theirRound2.combinedPubKey, baseGen, params))
  {
    printf("Invalid Round 2 jpake ZKP!\n");
    return ZZ(0);
  }

  /* Get exp for denominator exponentiation */
  ZZ expDivisor = MulMod(myRound1.keyPair2.privKey, pwd, params.q);

  /* Compute denominator */
  ZZ denom = dh_exp(theirRound1.pubKey2.pubKey, expDivisor, params);

  /* Compute inverse of denominator */
  ZZ denomInv = InvMod(denom, params.p);

  /* Get base #2 */
  ZZ base2 = MulMod(theirRound2.combinedPubKey.pubKey, denomInv, params.p);

  /* Perform last exponentiation */
  return dh_exp(base2, myRound1.keyPair2.privKey, params);
}

int dl_bruteForce(dhKeyPair_t keyPair, dhParameters_t params)
{
  ZZ privKey = ZZ(1);
  while(privKey < params.q)
  {
    if (dh_exp(params.g, privKey, params) == keyPair.pubKey)
    {
      assert(privKey == keyPair.privKey);
      return 1;
    }
    privKey = privKey + ZZ(1);
  }
  return 0;
}

void dl_updatePRIteration(ZZ& x, ZZ& a, ZZ&b, dhKeyPair_t keyPair, dhParameters_t params)
{
  ZZ xMod3 = x % ZZ(3);
  if (xMod3 == ZZ(0))
  {
    x = MulMod(x,x, params.p);
    a = MulMod(a,2, params.q);
    b = MulMod(b,2, params.q);
  } 
  else if (xMod3 == ZZ(1))
  {
    x = MulMod(x,params.g, params.p);
    a = AddMod(a,1, params.q);
  }
  else //xMod3 == 2
  {
    x = MulMod(x,keyPair.pubKey, params.p);
    b = AddMod(b,1, params.q);
  }
  return;
}

int dl_pollardRho(dhKeyPair_t keyPair, dhParameters_t params)
{
  /* Single iteration values */
  ZZ x = ZZ(1);
  ZZ a = ZZ(0);
  ZZ b = ZZ(0);

  /* Double iteration values */
  ZZ X = x;
  ZZ A = a;
  ZZ B = b;

  for (ZZ i = ZZ(1); i < params.q; i = i + ZZ(1))
  {
    dl_updatePRIteration(x,a,b, keyPair, params);
    dl_updatePRIteration(X,A,B, keyPair, params);
    dl_updatePRIteration(X,A,B, keyPair, params);

    if (x == X)
    {
      ZZ r = SubMod(B,b,params.q);
      ZZ rInv = InvMod(r,params.q);
      ZZ aDiff = SubMod(a,A,params.q);
      ZZ privKey = MulMod(rInv,aDiff,params.q);
      assert(privKey == keyPair.privKey);
      return 1;
    }
  }
  return 0;
}
