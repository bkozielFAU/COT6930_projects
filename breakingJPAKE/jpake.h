#ifndef RSA_H
#define RSA_H
  
#include <NTL/ZZ.h>

using namespace NTL;

#define RAND_EXP_LENGTH (32) /* Use a random 32-bit value to scale the exponent to a different representation */
#define RAND_BLIND_LENGTH (256) /* Use a random 256-bit value to blind the base to a different representation */

typedef struct
{
  ZZ g; // Generator value g
  ZZ p; // Modulus p
  ZZ q; // Prime subgroup order q
} dhParameters_t;

typedef struct
{
  ZZ pubKey;  // Public Key
  ZZ privKey; // Private Key
} dhKeyPair_t;

typedef struct
{
  ZZ pubKey; // Public Key
  ZZ zkpV;   // ZKP sig V
  ZZ zkpR;   // ZKP sig R
} jpakePubKey_t;

typedef struct
{
  dhKeyPair_t keyPair1; //g1
  dhKeyPair_t keyPair2; //g2
  jpakePubKey_t pubKey1; //From g1
  jpakePubKey_t pubKey2; //From g2
} jpakeRound1Struct_t;

typedef struct
{
  dhKeyPair_t combinedKeyPair; //A = (g1*g2*g3)^(x2*s)
  jpakePubKey_t combinedPubKey; //From A
} jpakeRound2Struct_t;


dhParameters_t dh_findParams(int primeLen, int orderLen);

ZZ dh_exp(ZZ base, ZZ exp, dhParameters_t params);

dhKeyPair_t dh_genKeyPair(ZZ gen, dhParameters_t params);

jpakePubKey_t jpake_genZKP(dhKeyPair_t keyPair, ZZ gen, dhParameters_t params);

int jpake_verifyZKP(jpakePubKey_t zkpIn, ZZ gen, dhParameters_t params); 

jpakeRound1Struct_t jpake_roundOne(dhParameters_t params);

jpakeRound2Struct_t jpake_roundTwo(jpakeRound1Struct_t myRound1, jpakeRound1Struct_t theirRound1, ZZ pwd, dhParameters_t params);

ZZ jpake_roundThree(jpakeRound1Struct_t myRound1, jpakeRound1Struct_t theirRound1, jpakeRound2Struct_t myRound2, jpakeRound2Struct_t theirRound2, ZZ pwd, dhParameters_t params);

int dl_bruteForce(dhKeyPair_t keyPair, dhParameters_t params);

int dl_pollardRho(dhKeyPair_t keyPair, dhParameters_t params);

#endif
