#ifndef RSA_H
#define RSA_H
  
#include <NTL/ZZ.h>

using namespace NTL;

#define RAND_EXP_LENGTH (32) /* Use a random 32-bit value to scale the exponent to a different representation */
#define RAND_BLIND_LENGTH (256) /* Use a random 256-bit value to blind the base to a different representation */

typedef struct
{
  ZZ n; // Modulus n = p*q
  ZZ e; // Public exponent e
} rsaPublicParameters_t;

typedef struct
{
  ZZ p;    // Prime p
  ZZ q;    // Prime q
  ZZ n;    // Modulus n = p*q
  ZZ e;    // Public exponent e
  ZZ d;    // Private exponent d
  ZZ d_p;  // d mod (p-1)
  ZZ d_q;  // d mod (q-1)
  ZZ q_inv; // q^-1 mod p
} rsaPrivateParameters_t;

ZZ rsa_encrypt(ZZ message, rsaPublicParameters_t pubParams);

ZZ rsa_decrypt(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRT(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRTFault(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRTRandExp(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRTBlind(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRTSecure(ZZ ciphertext, rsaPrivateParameters_t privParams);

ZZ rsa_decryptCRT_fault(ZZ ciphertext, rsaPrivateParameters_t privParams, int bitToCorrupt);

ZZ rsa_decryptTimingWeakness(ZZ ciphertext, rsaPrivateParameters_t privParams);

#endif
