/* segment chank11 in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);
void chanrate(chan *cpnt, chanparm *chp);

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);

#ifdef __cplusplus
}
#endif

double ncabs(double x);

/*--------------------------------------------*/

/* From: 

Labro AJ Priest MF, Lacroix JJ, Snyders DJ, Bezanilla F. (2015) Kv3.1 uses a
timely resurgent K(+) current to secure action potential repolarization.  
Nat Commun. 6:10173. doi: 10.1038/ncomms10173. 
https://www.ncbi.nlm.nih.gov/pubmed/26673941
/* */

/*--------------------------------------------*/
/*
 
Supplementary Table 1: 
Parameters for 4-State Markov Model.

Rate constants       Black, pink, green trace  Red trace	 Blue trace
ap / zp			0.05 / 3.5		 0.05 / 3.5	 0.05 / 3.5
bp / zp			0.15 / 3.5 		 0.15 / 3.5	 0.15 / 3.5
al / zl			  6  / 0.4		    6 / 0.4	    6 / 0.4
bl / zl			0.6  / 0.4		  1.8 / 0.4	  0.4 / 0.4
as / zs			  1  / 0.001		    1 / 0.001	    1 / 0.001
bs / zs			0.8  / 0.001		  0.8 / 0.001	  0.8 / 0.001

	a = a0(exp(   z * (V - Vh)/kT)) 
	b = b0(exp( - z * (V - Vh)/kT))
	
	where Vh = 6.2 mV

 - - - - - - - - - - - - - - - - - - - - - - 

Markov diagram for Kv3.1b channel:

               ap (zap)
   Resting  0    ----->  1 Pre-active
   (Closed)    <-----      (Closed)
	       bp (zbp) 


                 bl(zbl)   /|    |
                            |    |  al (zal)
                            |    |/

                                             as (zas)
		  2   Relaxed-pre-active      ----->   3 Relaxed-active
                          (Closed)            <-----      (Open)
	                                     bs (zbs) 
/* */
/*--------------------------------------------*/

#define    K11RTF  (RF*300*MVOLT)                // kT = RT/F = 25.85 mV at 300 deg K

/* For the purpose of normalizing kinetics to 22 deg C, */
/*  we assume Q10=2 */

/* #define K11A27RATE  (exp(log(2.0) * (BASETC- 27.0)/10.0); */
#define    K11A27RATE  0.70710678

/*-------------------------------------------*/

double calck11ab (double v, int func)

/* Calculate K type 11 (Kv3.1b) rate functions alpha, beta.
 
   Rate functions are taken from Labro et al., 2015,
   Nat Commun. 6:10173. doi: 10.1038/ncomms10173. 
   https://www.ncbi.nlm.nih.gov/pubmed/26673941
*/

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

#define Vh 6.2

#define ap 0.05
#define al 6
#define as 1.0

#define bp 0.15
#define bl 0.6
#define bs 0.8

#define zp 3.5
#define zl 0.4
#define zs 0.001

/* 
#define bl   1.8    for red trace
#define bl   0.4    for blue trace
*/

{
   double val;

  if (func > 6) func = 6;
  switch (func) {

    case 1:  val = ap * exp (  zp * (v-Vh)/K11RTF);  break; /* ap */
    case 2:  val = al * exp (  zl * (v-Vh)/K11RTF);  break; /* al */
    case 3:  val = as * exp (  zs * (v-Vh)/K11RTF);  break; /* as */

    case 4:  val = bp * exp ( -zp * (v-Vh)/K11RTF);  break; /* bp */
    case 5:  val = bl * exp ( -zl * (v-Vh)/K11RTF);  break; /* bl */
    case 6:  val = bs * exp ( -zs * (v-Vh)/K11RTF);  break; /* bs */

  }
  return val * K11A27RATE*MSSEC;
}

/*--------------------------------------------*/

double k11ap(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*--------------------------------------------*/

double k11al(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*--------------------------------------------*/

double k11as(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*--------------------------------------------*/

double k11bp(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*--------------------------------------------*/

double k11bl(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[4];
}

/*--------------------------------------------*/

double k11bs(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[5];
}


/*--------------------------------------------*/

chantype *makk11(void)
{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 4;                          /* 4 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value, (or if >1 must be nfval) */
   ch=makchantype(K,11,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dku;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 6;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck11ab;       /* default rate function */
   parm[0].funcname = (char *)"calck11ab"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = k11ap;
   spnt[0].ratemul[0] = 1.0*n;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;              /*   0 <-> 1 <-> 2 <-> 3 */
   spnt[1].trate  [0] = k11al;
   spnt[1].ratemul[0] = 1.0*n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = k11bp;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = k11as;
   spnt[2].ratemul[0] = 1.0*n;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = k11bl;
   spnt[2].ratemul[1] = 1.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 1;             /* state 4 = the open state */
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = k11bs;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 1;
 
   return ch;
}

/*----------------------------------------*/


/*
Taken from Labro et al., 2015,

   Nat Commun. 6:10173. doi: 10.1038/ncomms10173. 
   https://www.ncbi.nlm.nih.gov/pubmed/26673941

Simulations. 

Simulations were performed in MATLAB using a sequential four-state Markov model
of the VSD (Fig. 7b). To produce ionic currents, the fourth state was made an
open state, and its open probability was used to determine the ionic current
based on a reversal potential of -58 mV. Rate parameters at a potential of 0
mV were estimated from Kv3.1b gating and ionic current kinetic data using n =
a0/(a0 + b0) and t0 = 1/(a0 + b0), where n is the proportion of channels in the
final state, a0 is the rate constant going to that state, b0 is the rate
constant going from that state and t0 is the kinetic t of the experimental
state change. For each voltage, rate constants were determined using 
a = a0(exp(z(V - V1/2)/kT)) and b = b0(exp( - z(V - V1/2)/kT)), where V1/2 = 6.2 mV
and zp, zl and zs are 3.5, 0.4 and 0.001 respectively, to sum to a total
apparent valence of 3.9, as taken from the experimental QV curve for the Kv3.1b
channels.
/* */

#undef Vh 

#undef ap
#undef al
#undef as

#undef zp
#undef zl
#undef zs

#undef bp
#undef bl
#undef bs

