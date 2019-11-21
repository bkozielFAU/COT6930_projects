
#include "rsa.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/* Use a timing-weak binary ladder for modular exponentation. */
static ZZ PowerModBinary(ZZ a, ZZ e, ZZ n)
{
  if (e == 0) return ZZ(1);

  long k = NumBits(e);

  ZZ res, tmp;
  res = 1;
  tmp = a;

  for (long i = 0; i < k; i++) {
    if (bit(e,i) == 1) res = (res*tmp) % n;
    tmp = (tmp*tmp) % n;
  }
  return res;
}

ZZ rsa_encrypt(ZZ message, rsaPublicParameters_t pubParams)
{
  /* c = m^e mod n */
  return PowerMod(message, pubParams.e, pubParams.n);
}

ZZ rsa_decrypt(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  /* m = c^d mod n */
  return PowerMod(ciphertext, privParams.d, privParams.n);
}

ZZ rsa_decryptCRT(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = c^d_p mod p */
  m_1 = PowerMod(ciphertext % privParams.p, privParams.d_p, privParams.p);
  /* m_2 = c^d_q mod q */
  m_2 = PowerMod(ciphertext % privParams.q, privParams.d_q, privParams.q);

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  return (m_2 + (h*privParams.q)); 
}

ZZ rsa_decryptCRTFault(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = c^d_p mod p */
  m_1 = PowerMod(ciphertext % privParams.p, privParams.d_p, privParams.p);
  /* m_2 = c^d_q mod q */
  m_2 = PowerMod(ciphertext % privParams.q, privParams.d_q, privParams.q);

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  ZZ m = m_2 + (h*privParams.q); 
  assert(compare(PowerMod(m, privParams.e, privParams.n),ciphertext) == 0);
  return m;
}

ZZ rsa_decryptCRTRandExp(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  /* Find random representations of the exponents */
  /* d_p = d_p + rand*(p-1) */
  ZZ randExpP = privParams.d_p + (RandomBits_ZZ(RAND_EXP_LENGTH)*(privParams.p-1)); 
  /* d_q = d_q + rand*(q-1) */
  ZZ randExpQ = privParams.d_q + (RandomBits_ZZ(RAND_EXP_LENGTH)*(privParams.q-1)); 

  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = c^d_p mod p */
  m_1 = PowerMod(ciphertext % privParams.p, randExpP, privParams.p);
  /* m_2 = c^d_q mod q */
  m_2 = PowerMod(ciphertext % privParams.q, randExpQ, privParams.q);

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  ZZ m = m_2 + (h*privParams.q); 
  return m;
}

ZZ rsa_decryptCRTBlind(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  ZZ randBlindP = RandomBits_ZZ(RAND_BLIND_LENGTH);
  ZZ randBlindQ = RandomBits_ZZ(RAND_BLIND_LENGTH);

  /* compute r_p^-1 mod p, r_q^-d mod q */ 
  ZZ randBlindPInv = PowerMod(randBlindP % privParams.p, privParams.p - 1 - privParams.d_p, privParams.p); 
  ZZ randBlindQInv = PowerMod(randBlindQ % privParams.q, privParams.q - 1 - privParams.d_q, privParams.q); 

  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = (r*c)^d_p mod p */
  m_1 = PowerMod((ciphertext*randBlindP) % privParams.p, privParams.d_p, privParams.p);
  m_1 = (m_1*randBlindPInv) % privParams.p;
  /* m_2 = (r*c)^d_q mod q */
  m_2 = PowerMod((ciphertext*randBlindQ) % privParams.q, privParams.d_q, privParams.q);
  m_2 = (m_2*randBlindQInv) % privParams.q;

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  ZZ m = m_2 + (h*privParams.q); 
  return m;
}

ZZ rsa_decryptCRTSecure(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  ZZ randBlindP = RandomBits_ZZ(RAND_BLIND_LENGTH);
  ZZ randBlindQ = RandomBits_ZZ(RAND_BLIND_LENGTH);

  /* compute r_p^-1 mod p, r_q^-d mod q */ 
  ZZ randBlindPInv = PowerMod(randBlindP % privParams.p, privParams.p - 1 - privParams.d_p, privParams.p); 
  ZZ randBlindQInv = PowerMod(randBlindQ % privParams.q, privParams.q - 1 - privParams.d_q, privParams.q); 

  /* Find random representations of the exponents */
  /* d_p = d_p + rand*(p-1) */
  ZZ randExpP = privParams.d_p + (RandomBits_ZZ(RAND_EXP_LENGTH)*(privParams.p-1)); 
  /* d_q = d_q + rand*(q-1) */
  ZZ randExpQ = privParams.d_q + (RandomBits_ZZ(RAND_EXP_LENGTH)*(privParams.q-1)); 

  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = c^d_p mod p */
  m_1 = PowerMod((ciphertext*randBlindP) % privParams.p, randExpP, privParams.p);
  m_1 = (m_1*randBlindPInv) % privParams.p;
  /* m_2 = c^d_q mod q */
  m_2 = PowerMod((ciphertext*randBlindQ) % privParams.q, randExpQ, privParams.q);
  m_2 = (m_2*randBlindQInv) % privParams.q;

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  ZZ m = m_2 + (h*privParams.q); 
  assert(compare(PowerMod(m, privParams.e, privParams.n),ciphertext) == 0);
  return m;
}

ZZ rsa_decryptCRT_fault(ZZ ciphertext, rsaPrivateParameters_t privParams, int bitToCorrupt)
{
  /* perform same steps as normal CRT */
  /* Compute message components modulo p and modulo q */
  ZZ m_1, m_2, h;
  /* m_1 = c^d_p mod p */
  m_1 = PowerMod(ciphertext % privParams.p, privParams.d_p, privParams.p);
  /* m_2 = c^d_q mod q */
  m_2 = PowerMod(ciphertext % privParams.q, privParams.d_q, privParams.q);

  /* Corrupt a random bit in m_2 */
  ZZ bitXOR = ZZ(1) << bitToCorrupt;
  m_2 = m_2 ^ bitXOR;

  /* h = (q_inv * (m_1 - m_2)) mod p */
  h = (privParams.q_inv * (m_1 - m_2)) % privParams.p;

  /* m = m_2 + h*q */
  return (m_2 + (h*privParams.q)); 

}

ZZ rsa_decryptTimingWeakness(ZZ ciphertext, rsaPrivateParameters_t privParams)
{
  /* m = c^d mod n */
  return PowerModBinary(ciphertext, privParams.d, privParams.n);
}
