/* segment chanclca in program nc */

/* sets up Ca-activated Chloride channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "control.h"
#include "ncio.h"

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);
void makchanimpl (chanparm *chp);
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

extern int interp;		/* running interpreter */

double ncabs(double x);
double qrate(chanparm *chp);
double rca(chan *cpnt);
double rca2(chan *cpnt);
double rt(chan *cpnt);
void varcopyu(void);
attrib *getattrib(int elnum);
chantype *getchantype(int ctype, int stype);
double calcchaninf(double v, chanparm *chp, double arate, double brate);
double calcchantau(double v, chanparm *chp, double arate, double brate);

/*----------------------------------------*/

double calcclca1 (double v, int func)

{
	  return 1.0;
}

/*--------------------------------------------*/

chantype *makclca1(void)

/* Sequential-state ClCa channel, 
   Reference: Lalonde MR, Kelley ME, Barnes S (2008) 
   	      Calcium-activated chloride channels in the retina.
	      Channels 2:1-9.

   Markov diagram taken from:  SKCa channel (skca4) described in:
       Hirschberg et al., (1998) J. Gen. Physiol. 111: 565-581

  in which a 6-state markov diagram is provided to describe the
cloned rSK2 channel. Two sets of rate constants were provided, one for
normal high-probability gating, and another for simulating
"low-probability" gating behavior (not used here).  This channel
has no gating voltage sensitivity.

        f1      f2        f3        f4
    0  <->  1  <->   2   <->    3  <-> 4
        r1      r2        r3        r4


State 0 is unbound
      1 is bound with 1 Ca++
      2 is bound with 2 Ca++
      3 is bound with 3 Ca++
      4 is open, short duration

From Hirschberg et al (1998): Single channel currents were
recorded from Xenopus oocytes expressing the apamin-sensitive
clone rSK2. Channel activity was detectable in 0.2 uM [Ca]i and
was maximal above 2 uM [Ca]i.  Analysis of stationary currents
revealed two open times and three closed times, with only the
longest closed time being ca dependent, decreasing with
increasing [Ca]i concentrations.  In addition, elevated [Ca]i
concentrations resulted in a larger percentage of long openings
and short closures.  Membrane voltage did not have significant
effects on either open or closed times.  The open probability was
~0.6 in 1 uM free [Ca].

To set greater [Ca] sensitivity, increase "arate" by setting 
to a value less than 1, i.e. "taua = 0.1" will increase the [Ca] 
sensitivity by a factor of 10.

*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define ca4f1 40e6
#define ca4r1 1
#define ca4f2 40e6
#define ca4r2 1
#define ca4f3 40e6
#define ca4r3 1
#define ca4f4 40e6
#define ca4r4 1
#define ca4f5 40e6
#define ca4r5 1

   nstate = 6;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(ClCa,1,nstate,nparm,nq,dbasetc);  /* make chan state info */
   ch->unitary = dclcasu;		/* 14.2 pS @22 makes 22 ps @ 35 deg */
   ch->perm[PCL]  = 1.0;		/* permeability to Cl- ions */
   ch->perm[PNA] = 0;			/* permeability to Na+ ions / perm K */
   ch->perm[PCA] = 0;			/* permeability to Ca++ ions /perm K */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calcclca1;        /* dummy rate func */
   parm[0].funcname = (char *)"calcclca1";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rca rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rca;
   spnt[0].ratemul[0] = ca4f1;
   spnt[0].rateo  [0] = 0;      /* arate */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rca;
   spnt[1].ratemul[0] = ca4f2;
   spnt[1].rateo  [0] = 0;	/* arate */
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ca4r1;
   spnt[1].rateo  [1] = 1;	/* brate */

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rca;             /* function of ca conc, time */
   spnt[2].ratemul[0] = ca4f3;
   spnt[2].rateo  [0] = 0;      /* arate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;             /* function of time */
   spnt[2].ratemul[1] = ca4r2;
   spnt[2].rateo  [1] = 1;      /* brate */
 
   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rca;             /* function of time */
   spnt[3].ratemul[0] = ca4f4;
   spnt[3].rateo  [0] = 0;      /* arate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;             /* function of time */
   spnt[3].ratemul[1] = ca4r3;
   spnt[3].rateo  [1] = 1;      /* brate */
 
   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rca;             /* function of time */
   spnt[4].ratemul[0] = ca4f5;
   spnt[4].rateo  [0] = 0;      /* arate */
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = rt;             /* function of time */
   spnt[4].ratemul[1] = ca4r4;
   spnt[4].rateo  [1] = 1;      /* brate */
 
   spnt[5].numtrans   = 1;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;             /* function of time */
   spnt[5].ratemul[0] = ca4r5;
   spnt[5].rateo  [0] = 1;      /* brate */
 
   return ch;
}

/*--------------------------------------------*/

