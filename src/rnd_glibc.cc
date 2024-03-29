/* Random number routine from libc source code */
/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *	@(#)random.c	5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 *
 * Modified to compile without GNU include files. RGS 11/94, 9/98.
 */

//#include <limits.h>
//#include <stddef.h>
//#include <stdlib.h>
#include <stdio.h>
#include "ncsub.h"

/* structure modified from glibc version of stdlib.h: */

struct random_data
  {
    int *fptr;			/* Front pointer.  */
    int *rptr;			/* Rear pointer.  */
    int *state;			/* Array of state values.  */
    int rand_type;		/* Type of random number generator.  */
    int rand_deg;		/* Degree of random number generator.  */
    int rand_sep;		/* Distance between front and rear.  */
    int *end_ptr;		/* Pointer behind state table.  */
  };

#define MAXINTR 0x7fffffff

/* An improved random number generation package.  In addition to the standard
   rand()/srand() like interface, this package also has a special state info
   interface.  The initstate() routine is called with a seed, an array of
   bytes, and a count of how many bytes are being passed in; this array is
   then initialized to contain information for random number generation with
   that much state information.  Good sizes for the amount of state
   information are 32, 64, 128, and 256 bytes.  The state can be switched by
   calling the setstate() function with the same array as was initialized
   with initstate().  By default, the package runs with 128 bytes of state
   information and generates far better random numbers than a linear
   congruential generator.  If the amount of state information is less than
   32 bytes, a simple linear congruential R.N.G. is used.  Internally, the
   state information is treated as an array of longs; the zeroth element of
   the array is the type of R.N.G. being used (small integer); the remainder
   of the array is the state information for the R.N.G.  Thus, 32 bytes of
   state information will give 7 longs worth of state information, which will
   allow a degree seven polynomial.  (Note: The zeroth word of state
   information also has some other information stored in it; see setstate
   for details).  The random number generation technique is a linear feedback
   shift register approach, employing trinomials (since there are fewer terms
   to sum up that way).  In this approach, the least significant bit of all
   the numbers in the state table will act as a linear feedback shift register,
   and will have period 2^deg - 1 (where deg is the degree of the polynomial
   being used, assuming that the polynomial is irreducible and primitive).
   The higher order bits will have longer periods, since their values are
   also influenced by pseudo-random carries out of the lower bits.  The
   total period of the generator is approximately deg*(2**deg - 1); thus
   doubling the amount of state information has a vast influence on the
   period of the generator.  Note: The deg*(2**deg - 1) is an approximation
   only good for large deg, when the period of the shift register is the
   dominant factor.  With deg equal to seven, the period is actually much
   longer than the 7*(2**7 - 1) predicted by this formula.  */



/* For each of the currently supported random number generators, we have a
   break value on the amount of state information (you need at least this many
   bytes of state info to support this random number generator), a degree for
   the polynomial (actually a trinomial) that the R.N.G. is based on, and
   separation between the two lower order coefficients of the trinomial.  */

/* Linear congruential.  */
#define	TYPE_0		0
#define	BREAK_0		8
#define	DEG_0		0
#define	SEP_0		0

/* x**7 + x**3 + 1.  */
#define	TYPE_1		1
#define	BREAK_1		32
#define	DEG_1		7
#define	SEP_1		3

/* x**15 + x + 1.  */
#define	TYPE_2		2
#define	BREAK_2		64
#define	DEG_2		15
#define	SEP_2		1

/* x**31 + x**3 + 1.  */
#define	TYPE_3		3
#define	BREAK_3		128
#define	DEG_3		31
#define	SEP_3		3

/* x**63 + x + 1.  */
#define	TYPE_4		4
#define	BREAK_4		256
#define	DEG_4		63
#define	SEP_4		1


/* Array versions of the above information to make code run faster.
   Relies on fact that TYPE_i == i.  */

#define	MAX_TYPES	5	/* Max number of types above.  */

static const int degrees[MAX_TYPES] = { DEG_0, DEG_1, DEG_2, DEG_3, DEG_4 };
static const int seps[MAX_TYPES] = { SEP_0, SEP_1, SEP_2, SEP_3, SEP_4 };

/* Initially, everything is set up as if from:
	initstate(1, randtbl, 128);
   Note that this initialization takes advantage of the fact that srandom
   advances the front and rear pointers 10*rand_deg times, and hence the
   rear pointer which starts at 0 will also end up at zero; thus the zeroth
   element of the state information, which contains info about the current
   position of the rear pointer is just
	(MAX_TYPES * (rptr - state)) + TYPE_3 == TYPE_3.  */

static int randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };


static struct random_data unsafe_state =
  {
/* FPTR and RPTR are two pointers into the state info, a front and a rear
   pointer.  These two pointers are always rand_sep places aparts, as they
   cycle through the state information.  (Yes, this does mean we could get
   away with just one pointer, but the code for random is more efficient
   this way).  The pointers are left positioned as they would be from the call:
	initstate(1, randtbl, 128);
   (The position of the rear pointer, rptr, is really 0 (as explained above
   in the initialization of randtbl) because the state table pointer is set
   to point to randtbl[1] (as explained below).)  */

    /* fptr */ &randtbl[SEP_3 + 1],
    /* rptr */ &randtbl[1],

/* The following things are the pointer to the state information table,
   the type of the current generator, the degree of the current polynomial
   being used, and the separation between the two pointers.
   Note that for efficiency of random, we remember the first location of
   the state information, not the zeroth.  Hence it is valid to access
   state[-1], which is used to store the type of the R.N.G.
   Also, we remember the last location, since this is more efficient than
   indexing every time to find the address of the last element to see if
   the front and rear pointers have wrapped.  */

    /* state */ &randtbl[1],

    /* rand_type */ TYPE_3,
    /* rand_deg  */ DEG_3,
    /* rand_sep  */ SEP_3,

    /* end_ptr */ &randtbl[sizeof (randtbl) / sizeof (randtbl[0])]
};

/* POSIX.1c requires that there is mutual exclusion for the `rand' and
   `srand' functions to prevent concurrent calls from modifying common
   data.  */

/* __libc_lock_define_initialized (static, lock) */

/* Initialize the random number generator based on the given seed.  If the
   type is the trivial no-state-information type, just remember the seed.
   Otherwise, initializes state[] based on the given "seed" via a linear
   congruential generator.  Then, the pointers are set to known locations
   that are exactly rand_sep places apart.  Lastly, it cycles the state
   information a given number of times to get rid of any initial dependencies
   introduced by the L.C.R.N.G.  Note that the initialization of randtbl[]
   for default usage relies on values produced by this routine.  */

static int r_random(void);

int r_srandom_r (unsigned int x, struct random_data *buf)
{

  if (buf == NULL || buf->rand_type < TYPE_0 || buf->rand_type > TYPE_4)
    return -1;

  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  buf->state[0] = x ? x : 1;
  if (buf->rand_type != TYPE_0)
    {
      long int i, word;
      int *dst;
      dst = buf->state;
      word = *dst;
      for (i = 1; i < buf->rand_deg; ++i)
	{
	  /* This does:
	       state[i] = (16807 * state[i - 1]) % 2147483647;
	     but avoids overflowing 31 bits.  */
	  long int hi = word / 127773;
	  long int lo = word % 127773;
	  word = 16807 * lo - 2836 * hi;
	  if (word < 0) word += 2147483647;
	  *++dst = word;
	}
      buf->fptr = &buf->state[buf->rand_sep];
      buf->rptr = &buf->state[0];
      for (i = 0; i < 10 * buf->rand_deg; ++i)
	{
	  (void) r_random ();
	}
    }
  return 0;
}

/* Initialize the state information in the given array of N bytes for
   future random number generation.  Based on the number of bytes we
   are given, and the break values for the different R.N.G.'s, we choose
   the best (largest) one we can and set things up for it.  srandom is
   then called to initialize the state information.  Note that on return
   from srandom, we set state[-1] to be the type multiplexed with the current
   value of the rear pointer; this is so successive calls to initstate won't
   lose this information and will be able to restart with setstate.
   Note: The first thing we do is save the current state, if any, just like
   setstate so that it doesn't matter when initstate is called.
   Returns a pointer to the old state.  */

int r_initstate_r (unsigned long int seed, void *arg_state, 
			int n, struct random_data *buf)
{
  if (buf == NULL)
    return -1;

  if (n < BREAK_1)
    {
      if (n < BREAK_0)
	{
/*	  __set_errno (EINVAL); */
	  return -1;
	}
      buf->rand_type = TYPE_0;
      buf->rand_deg = DEG_0;
      buf->rand_sep = SEP_0;
    }
  else if (n < BREAK_2)
    {
      buf->rand_type = TYPE_1;
      buf->rand_deg = DEG_1;
      buf->rand_sep = SEP_1;
    }
  else if (n < BREAK_3)
    {
      buf->rand_type = TYPE_2;
      buf->rand_deg = DEG_2;
      buf->rand_sep = SEP_2;
    }
  else if (n < BREAK_4)
    {
      buf->rand_type = TYPE_3;
      buf->rand_deg = DEG_3;
      buf->rand_sep = SEP_3;
    }
  else
    {
      buf->rand_type = TYPE_4;
      buf->rand_deg = DEG_4;
      buf->rand_sep = SEP_4;
    }

  buf->state = &((int *) arg_state)[1];	/* First location.  */
  /* Must set END_PTR before srandom.  */
  buf->end_ptr = &buf->state[buf->rand_deg];

  r_srandom_r (seed, buf);

  if (buf->rand_type == TYPE_0)
    buf->state[-1] = buf->rand_type;
  else
    buf->state[-1] = (MAX_TYPES * (buf->rptr - buf->state)) + buf->rand_type;

  return 0;
}

char *r_initstate (unsigned long int seed, void *arg_state, int n)
{
  static char *ostate;

  ostate = (char *) &unsafe_state.state[-1];
  r_initstate_r (seed, arg_state, n, &unsafe_state);
  return ostate;
}

/* Restore the state from the given state array.
   Note: It is important that we also remember the locations of the pointers
   in the current state information, and restore the locations of the pointers
   from the old state information.  This is done by multiplexing the pointer
   location into the zeroth word of the state information. Note that due
   to the order in which things are done, it is OK to call setstate with the
   same state as the current state
   Returns a pointer to the old state information.  */

int r_setstate_r (void *arg_state, struct random_data *buf)
{
  static int *new_state;
  static int type,rear;

  new_state = (int *) arg_state;
  if (buf == NULL)
    return -1;

  if (buf->rand_type == TYPE_0)
    buf->state[-1] = buf->rand_type;
  else
    buf->state[-1] = (MAX_TYPES * (buf->rptr - buf->state)) + buf->rand_type;

  type = new_state[0] % MAX_TYPES;
  rear = new_state[0] / MAX_TYPES;

  switch (type)
    {
    case TYPE_0:
    case TYPE_1:
    case TYPE_2:
    case TYPE_3:
    case TYPE_4:
      buf->rand_type = type;
      buf->rand_deg = degrees[type];
      buf->rand_sep = seps[type];
      break;
    default:
      /* State info munged.  */
      /*  __set_errno (EINVAL); */
      return -1;
    }

  buf->state = &new_state[1];
  if (buf->rand_type != TYPE_0)
    {
      buf->rptr = &buf->state[rear];
      buf->fptr = &buf->state[(rear + buf->rand_sep) % buf->rand_deg];
    }
  /* Set end_ptr too.  */
  buf->end_ptr = &buf->state[buf->rand_deg];

  return 0;
}

char *r_setstate (void *arg_state)

{
  static char *ostate;

  ostate = (char *) &unsafe_state.state[-1];
  if (r_setstate_r (arg_state, &unsafe_state) < 0)
    ostate = (char *)NULL;
  return ostate;
}


/* If we are using the trivial TYPE_0 R.N.G., just do the old linear
   congruential bit.  Otherwise, we do our fancy trinomial stuff, which is the
   same in all the other cases due to all the global variables that have been
   set up.  The basic operation is to add the number at the rear pointer into
   the one at the front pointer.  Then both pointers are advanced to the next
   location cyclically in the table.  The value returned is the sum generated,
   reduced to 31 bits by throwing away the "least random" low bit.
   Note: The code takes advantage of the fact that both the front and
   rear pointers can't wrap on the same call by not testing the rear
   pointer if the front one has wrapped.  Returns a 31-bit random number.  */


int r_random (void)
{
   static int result;
   static int *state;
   static struct random_data *buf=&unsafe_state;

  state = buf->state; 

  //  /* make faster by removing "if" statement -- never use type 0 */
  //if (buf->rand_type == TYPE_0)
  //  {
  //    state[0] = ((state[0] * 1103515245) + 12345) & MAXINTR;
  //    result = state[0];
  //  }
  //else
    {	
	int *fptr = buf->fptr;
	int *rptr = buf->rptr;
	int *end_ptr = buf->end_ptr;
	int val;

      val = *fptr += *rptr;
      /* Chucking least random bit.  */
      result = (val >> 1) & MAXINTR;
      ++fptr;
      if (fptr >= end_ptr)
	{
	  fptr = state;
	  ++rptr;
	}
      else
	{
	  ++rptr;
	  if (rptr >= end_ptr)
	    rptr = state;
	}
      buf->fptr = fptr;
      buf->rptr = rptr;
    }
  return result;
}

#define MAXRAND 2147483647.0
#define INVMAXRAND 4.6566128752457969231e-10            /* 1 / MAXRAND */

double r_random_double(void)
{

  return (r_random() * INVMAXRAND);
}

