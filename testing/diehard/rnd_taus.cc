/* rnd_taus.cc in "nc" */

/* rng/taus.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>

/* This is a maximally equidistributed combined Tausworthe
   generator. The sequence is,

   x_n = (s1_n ^ s2_n ^ s3_n) 

   s1_{n+1} = (((s1_n & 4294967294) <<12) ^ (((s1_n <<13) ^ s1_n) >>19))
   s2_{n+1} = (((s2_n & 4294967288) << 4) ^ (((s2_n << 2) ^ s2_n) >>25))
   s3_{n+1} = (((s3_n & 4294967280) <<17) ^ (((s3_n << 3) ^ s3_n) >>11))

   computed modulo 2^32. In the three formulas above '^' means
   exclusive-or (C-notation), not exponentiation. Note that the
   algorithm relies on the properties of 32-bit unsigned integers (it
   is formally defined on bit-vectors of length 32). I have added a
   bitmask to make it work on 64 bit machines.

   We initialize the generator with s1_1 .. s3_1 = s_n MOD m, where
   s_n = (69069 * s_{n-1}) mod 2^32, and s_0 = s is the user-supplied
   seed.

   The theoretical value of x_{10007} is 2733957125. The subscript
   10007 means (1) seed the generator with s=1 (2) do six warm-up
   iterations, (3) then do 10000 actual iterations.

   The period of this generator is about 2^88.

   From: P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe
   Generators", Mathematics of Computation, 65, 213 (1996), 203--213.

   This is available on the net from L'Ecuyer's home page,

   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
   ftp://ftp.iro.umontreal.ca/pub/simulation/lecuyer/papers/tausme.ps 

   Update: April 2002

   There is an erratum in the paper "Tables of Maximally
   Equidistributed Combined LFSR Generators", Mathematics of
   Computation, 68, 225 (1999), 261--269:
   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps

	... the k_j most significant bits of z_j must be non-
	zero, for each j. (Note: this restriction also applies to the 
	computer code given in [4], but was mistakenly not mentioned in
	that paper.)
   
   This affects the seeding procedure by imposing the requirement
   s1 > 1, s2 > 7, s3 > 15.

   The generator taus2 has been added to satisfy this requirement.
   The original taus generator is unchanged.

   Update: November 2002

   There was a bug in the correction to the seeding procedure for s2.
   It affected the following seeds 254679140 1264751179 1519430319
   2274823218 2529502358 3284895257 3539574397 (s2 < 8).

*/

/* Modified by RGS 2003/6/7 

   Changes:

  Removed taus_set() (not used because of error)
  Renamed taus2_set() to taus_init(),
    changed order of args and added extra dummy arg to taus_init()
  Added taus_setstate
  Added static global definition of the state variable pointer
    static taus_state_t *taus_state = NULL;
    #define state taus_state

  Removed gsl types at end:

     static const gsl_rng_type taus_type =
     const gsl_rng_type *gsl_rng_taus = &taus_type;
     static const gsl_rng_type taus2_type =
     const gsl_rng_type *gsl_rng_taus2 = &taus2_type;

*/


typedef struct
  {
    unsigned long int s1, s2, s3;
  }
taus_state_t;

/* make global state variable pointer so rng remembers */

static taus_state_t taus_state_buf = {1882473311,238124573,323137};
static taus_state_t *taus_state = &taus_state_buf;

/*------------------------------------------------------*/

void taus_setstate (void *vstate)
{
  taus_state = (taus_state_t *) vstate;
}

/*------------------------------------------------------*/

unsigned long taus_get (void)
{
  taus_state_t *state = (taus_state_t *) taus_state;

#define MASK 0xffffffffUL
#define TAUSWORTHE(s,a,b,c,d) (((s &c) <<d) &MASK) ^ ((((s <<a) &MASK)^s) >>b)

  state->s1 = TAUSWORTHE (state->s1, 13, 19, 4294967294UL, 12);
  state->s2 = TAUSWORTHE (state->s2, 2, 25, 4294967288UL, 4);
  state->s3 = TAUSWORTHE (state->s3, 3, 11, 4294967280UL, 17);

  return (state->s1 ^ state->s2 ^ state->s3);
}

/*------------------------------------------------------*/

double taus_get_double (void)
{
  return taus_get () / 4294967296.0 ;
}

/*------------------------------------------------------*/

void taus_init (unsigned long int s, void *vstate, int siz)
{
  taus_state_t *state;

  taus_setstate(vstate);
  state = (taus_state_t *) taus_state;

  if (s == 0)
    s = 1;	/* default seed is 1 */

#define LCG(n) ((69069 * n) & 0xffffffffUL)
  state->s1 = LCG (s);
  if (state->s1 < 2) state->s1 += 2UL;
  state->s2 = LCG (state->s1);
  if (state->s2 < 8) state->s2 += 8UL;
  state->s3 = LCG (state->s2);
  if (state->s3 < 16) state->s3 += 16UL;

  /* "warm it up" */
  taus_get ();
  taus_get ();
  taus_get ();
  taus_get ();
  taus_get ();
  taus_get ();
  return;
}

