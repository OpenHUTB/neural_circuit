/* Tausworthe 258 random number generator */

/* From: https://apps.dtic.mil/sti/pdfs/AD1024837.pdf */
 
/* Modified by R.G.Smith, Nov 2021 */

/* Random number generator of the type: Linear Feedback Shift Register (LFSR) */

/* This version has a sequence length of 2^258 (4.63e77), much longer than taus11. */
/* This taus258 also requires 40 bytes for the state vs 16 for taus113, */
/* but it runs only 10% slower */
 
typedef struct
{
  unsigned long long int z1, z2, z3, z4, z5;	// 64-bit integers
}
taus258_state_t;

/* make global state variable pointer so rng remembers */

static taus258_state_t  taus258_state_buf = {0ULL,0ULL,0ULL,0ULL,0ULL};
static taus258_state_t *taus258_state = &taus258_state_buf;

/*-----------------------------------------------------------*/

void taus258_setstate (void *vstate)
{
  taus258_state = (taus258_state_t *) vstate;
}

/*-----------------------------------------------------------*/

#define c1 0xfffffffffffffffeULL // 18446744073709551614ULL
#define c2 0xfffffffffffffe00ULL // 18446744073709551104ULL
#define c3 0xfffffffffffff000ULL // 18446744073709547520ULL
#define c4 0xfffffffffffe0000ULL // 18446744073709420544ULL
#define c5 0xffffffffff800000ULL // 18446744073701163008ULL

unsigned long long taus258_get (void)
{
  static taus258_state_t *state;

  state = (taus258_state_t *) taus258_state;
  
  state->z1 = ((state->z1 & c1) << 10) ^ (((state->z1 <<  1) ^ state->z1) >> 53);
  state->z2 = ((state->z2 & c2) << 5)  ^ (((state->z2 << 24) ^ state->z2) >> 50);
  state->z3 = ((state->z3 & c3) << 29) ^ (((state->z3 <<  3) ^ state->z3) >> 23);
  state->z4 = ((state->z4 & c4) << 23) ^ (((state->z4 <<  5) ^ state->z4) >> 24);
  state->z5 = ((state->z5 & c5) << 8)  ^ (((state->z5 <<  3) ^ state->z5) >> 33);
  return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4 ^ state->z5);
}

/*-----------------------------------------------------------*/

/* 2^53 is the largest number with properties, as the internal IEEE representation */
/* of double uses 53 bits implicitly. */
 
double taus258_get_double (void)
{
    static long long z;
    const static double d = 1.0/(1ULL<<53);

    do z = taus258_get() >> 11; while(!z);
    return d * z;					// 2^53 double
}

/*-----------------------------------------------------------*/

void taus258_init (unsigned long long int s, void *vstate, int siz)
{
  static taus258_state_t *state;
  static unsigned long long c = 0x27bb2ee687b0b0fd, d = 0x891087b8e3b70cb1;

  taus258_setstate(vstate);
  state = taus258_state;

  do state->z1 = c*s+d;         while(state->z1  < 0x000002);
  do state->z2 = c*state->z1+d; while(state->z2  < 0x000200);
  do state->z3 = c*state->z2+d; while(state->z3  < 0x001000);
  do state->z4 = c*state->z3+d; while(state->z4  < 0x020000);
  do state->z5 = c*state->z4+d; while(state->z5  < 0x800000);

 /* Calling RNG ten times to satify recurrence condition */

  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
  taus258_get ();
 
 return;
}


