/* segment chanach in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);

double rnt(chan *spnt);
double rt(chan *spnt);
int set_perms (chantype *chtyp, iontab *ions, double svrev);

/*--------------------------------------------*/

double calcach1(double v, int func) 

{
return 1.0;
}

/*--------------------------------------------*/

chantype *makach1(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

/* 


/* Taken from:
 
  Keleshian AM, Edeson RO, Liu G-J and Madsen BW. (2000) 
  Biophys J 78:1-12.

   in which a 5-state ACh markov diagram is provided.  
   Several sets of rate constants were derived, with different
   combinations of identicality and independence between neighbor channels. 
   The one defined here is the one assuming identical and independence.

     2 kon          kon    
  0  <->      1     <->     2   
     koff           2 koff  
              |             |         
            ^ |             |  ^ 
           a1 | b1       a2 | b2 
              | \/          | \/     
              |             |

              3    <->      4 

State 0 is unbound
      1 is partially bound
      2 is fully bound
      3 is partially bound and open
      4 is fully bound and open
 
*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k, pca;
     int nstate, nparm, nq;

#define ach_kon  5.0e7 
#define ach_koff 7.6923e3
#define ach_a1   9.6875e3
#define ach_a2   0.96875e3
#define ach_b1   3.1e3
#define ach_b2   31.0e3

   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(ACH,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->unitary = dachu;
   ch->trconc = dstr;
   pca = dpcaampa;
   ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcach1;        /* default rate function */
   parm[0].funcname = (char *)"calcach1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[ACH-GLU]  =  1;		/* response to ACh */
   respamp[MLA-GLU]  = -1;		/* response to MLA */
   respamp[HEX-GLU]  = -1;		/* response to HEX */

   respamp[GLU-GLU]  = 0;		/* response to glutamate */
   respamp[AMPA-GLU] = 0;		/* response to AMPA */
   respamp[NMDA-GLU] = 0;		/* response to NMDA */
   respamp[CNQX-GLU] = 0;		/* response to CNQX */

   k=1;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = ach_kon * 2 * k;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0; 
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = ach_kon * k;
   spnt[1].rateo  [0] = 0;      /* offsm */
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ach_koff * k;
   spnt[1].rateo  [1] = 1;      /* offsm */
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = ach_b1 * k;
   spnt[1].rateo  [2] = 2;      /* offsh */

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;			/* function of time */
   spnt[2].ratemul[0] = ach_koff * 2 * k;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 4;
   spnt[2].trate  [1] = rt;	
   spnt[2].ratemul[1] = ach_b2 * k;
   spnt[2].rateo  [1] = 2;	/* offsh */

   spnt[3].numtrans   = 1;			/* state 3 = open state */
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 1;			
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = ach_a1 * k;
   spnt[3].rateo  [0] = 2;	/* offsh */

   spnt[4].numtrans   = 1;			/* state 3 = open state */
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = ach_a2 * k;	
   spnt[4].rateo  [0] = 2;      /* offsh */

   return ch;
}

#undef ach_kon
#undef ach_koff
#undef ach_a1
#undef ach_a2
#undef ach_b1
#undef ach_b2

