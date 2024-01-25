/* stochastic search  */

/* R.G.Smith	Mar, 2020 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "ncfuncs.h"

#define NCOEFF 12 		/* free parameters */

#define SAVAL 0
#define SAMAX 1
#define SAMIN 2
#define SASD  3		/* standard deviation */
#define SARNG 4		/* normalized standard deviation */
#define NSTV  5		/* number of starting values */

#define abs(x)          ((x) < 0 ? -(x) : (x))
#define max(x, y)       (((x) < (y)) ? (y) : (x))
#define min(x, y)       (((x) < (y)) ? (x) : (y))

/*------------------------------------------------------------*/

FILE *prout = NULL;

extern int info;
extern int rseed;
extern const char *p1;
extern const char *p2;
extern const char *p3;
extern const char *p4;
extern const char *p5;
extern const char *p6;
extern const char *p7;
extern const char *p8;
extern const char *p9;
extern const char *p10;
extern const char *p11;
extern const char *p12;

double sastart[NSTV][NCOEFF]; /* Starting, max, min values, set by user */
double satest[NCOEFF];	      /* The test values, used by simulation */
const char *saname[NCOEFF];

int sa_ilim;  		/* maximum # of iterations */
double sa_st;	  	/* starting temperature */
double sa_ftol;  	/* value of match for stopping */
double sa_mthr;		/* val of match to decr stuckcrit*/
double sa_tstk;		/* base level for tstuckcrit */
double sa_nstk;		/* number crit for stuck */
double sa_itb;		/* base for iter to the n power */
double sa_td;		/* temperature decrement */
double sa_tj;		/* temperature jump */
double sa_bj;		/* best match jump */
double sa_ybc;		/* best match criterion */
double sa_sdmul; 	/* multiplier for s.d. */

int sa_sdi;		/* compute sd criterion */
int sa_sanb;		/* trials to compute sd */
int calc_sdi;		/* calculate s.d. */
int sa_wavg;		/* calculate weighted average */

/* Simulated annealing routine to be called from the main script.
This routine calculates new values for the free parameters and
calls the simulation function in the main script.  The simulation
function runs the simulation (from a different set of test values for
the free parameters each time) and computes how close the
simulation is to a template using a matching function.
The simulated annealing procedure decides how to modify the test
values using only the value of the matching function as a guide.
This process is repeated until the difference between the
simulation and the template is minimized.  The best matching
parameters are left in the "bp" array. 

-----------------------------------------------------------

Run the ssrch() procedure, which will run the search, calling
the eval_sim() function to evaluate the match.  At the end, you can
print out the bp[] array for the best match.

The algorithm starts out with a high temperature, which sets the
width of the Gaussian distribution for each parameter (which is
also set by the max,min range values).  After several ("iter")
iterations the best match is selected as the initial test value
for the next set of iterations. The temperature is reduced (by
10-15%, set by sa_td) and the process repeats.  If no better
match (lower value for pe[]) is found, after a few more
temperature reduction steps the temperature is increased by a
factor (10 - 20, set by sa_tj), and the process of running
several iterations at each temperature reduction step is
repeated. If still stuck, the acceptable match value is increased
randomly, to allow the algorithm to skip over sub-optimal peaks.

Several parameters may need to be tuned. The "sa_mthr" parameter
sets the value of match below which tstuckcrit goes down. Normally
sa_mthr should be set about 2 log units below the initial match
value (i.e. for an initial match of 1e-10, set sa_mthr to 1e-12).
The "sa_ftol" parameter sets the ultimate threshold for stopping.
Normally it should be set about 4-5 log units below the initial
match.  The "sa_ilim" parameter sets the maximum mumber of
iterations.  You can set it to a lower value to enforce a time
limit on the computation.

The number of iterations per temperature step is set by an
automatic algorithm, which sets a larger number according to the
number of free parameters.  If a parameter has SAMAX and SAMIN
values within a ratio of 1e-3, that parameter is considered to be
fixed and the number of iterations is reduced by a factor of 2.

The "info" parameter sets the level of printouts.  A level of 2
prints out the current best match, and a level of 4 prints out
all the tested values.

Advice:

To reduce computation time, the most important factor is to
reduce the number of free parameters.  A simulation with 6 free
parameters that takes several thousand iterations to converge may
converge in only several hundred iterations with 4 free
parameters.  Another problem is correlated parameters. If 2 or
more parameters are correlated (i.e. a reduction in one can be
compensated for to some extent by an increase in the other) then
both their correct values may take several times the computation
time required to find one value.  Also it is always helpful to
constrain a parameter by setting a narrow range for SAMAX and
SAMIN.

*/

/*------------------------------------------------------------*/

/* Set by the Stochastic Search algorithm: */

int sa_rnd = 100;		/* sets which random number generator to use */
int sa_itot = 0;		/* total iterations */
int sa_yitit  = 0;		/* best in this iteration */
int sa_yitit2 = 0;		/* 2nd best in this iteration */
int sa_yitit3 = 0;		/* 3rd best in this iteration */
double sstdevthr = 1e-3;	/* stdev thresh for param to be counted free */

double sa_y  = 1e20;		/* latest y value */
double sa_ye = 1e20;		/* latest estimate for best y value */
double sa_yit = 1e20;		/* best y value for this round of iterations */
double sa_yit2 = 1e20;		/* 2nd best y value for this round of iterations */
double sa_yit3 = 1e20;		/* 3rd best y value for this round of iterations */
double sa_bit[NCOEFF];		/* estimate for best from this round */
double sa_bit2[NCOEFF];		/* estimate for 2nd best from this round */
double sa_bit3[NCOEFF];		/* estimate for 3rd best from this round */
double sa_pe[NCOEFF];		/* latest estimate for best test value */
double sa_tsav[NCOEFF];		/* temporary store */

#define SANB 100		/* max number of previous best values to save */
int sannb = 0;			/* number of previous best values saved */
double sa_tb[SANB];		/* best temperature */
double sa_yb[SANB];		/* best y val */
double sa_it[SANB];		/* iteration for best y val */
double sa_pb[SANB][NCOEFF];	/* best test values */
#define RITER 400		/* max number of iters to save */
double sa_ps[RITER][NCOEFF]={0}; /* iter saved values */

double sabest[NCOEFF];		/* The best values, avg of best */
double sabest_sq[NCOEFF];	/* The best values, avg squared of best */
double sabest_sd[NCOEFF];	/* The best values, stdev */

extern int tracelen;

/*------------------------------------------------------------*/

double get_model_data(double x, double *coeff, int n);

double eval_sim(double *coeff, int n, double *xydata, int datasize)
{
    int i;
    double sum, yval, dy;

  if (datasize<=0) return 0;
  for (sum=i=0; i<tracelen; i++) {			/* calculate match */
     yval = get_model_data(i,coeff,n);
     dy = xydata[datasize+i] - yval;
     sum += dy*dy; 
     if (info>=6) fprintf (stderr, "%d %-10.3g %-10.3g  %-10.3g %g\n",i,xydata[datasize+i], yval, dy, sum);
  } 
  return sum/tracelen;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void insert_b (double yb, double temptr, double *sabest, int nparms, int iter)

/* insert a local best run into the save list */

{
    int i, j, k;

  for (k=0; k<SANB; k++) { 
    if (yb < sa_yb[k]) {
       for (i=SANB-1; i>k; i--) {	/* make space by pushing down the previous ones */
          for (j=0; j<nparms; j++) {
	      sa_pb[i][j] = sa_pb[i-1][j];
          };
          sa_tb[i] = sa_tb[i-1];
          sa_yb[i] = sa_yb[i-1];
          sa_it[i] = sa_it[i-1];
       }
       if (++sannb > SANB) sannb = SANB;
       sa_tb[k] = temptr;
       sa_yb[k] = yb;
       sa_it[k] = iter;
      for (j=0; j<nparms; j++) sa_pb[k][j]=sabest[j]; /* save the new one */
      break;
    }  /* if */
  }   /* for (i=0;;)  */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void runstsrch(int iter, int nparms, double temptr, double *coeff, double *xydata, int datasize)

/* Run simulation a specified number of iterations, with a different */
/*  set of free params each time, given a temperature (stddev multiplier). */
/*  Range of each parameter is limited, so reject param values outside range */

{
    int i,ii,j,k,titot;
    double sstd,x;

  titot = sa_itot;
  for (k=1; k<=iter; k++) {
    for (j=0; j<nparms; j++) {
      if (sastart[SARNG][j] > sstdevthr) {
        sstd = temptr * sastart[SASD][j]; 
        satest[j] = sa_pe[j] + (rgasdev(sa_rnd)*sstd);
	// fprintf (prout,"j %d sstd %g %g %g %g\n",j,sstd, sastart[SARNG][j],  abs(sa_pe[j]), satest[j]);
        for (ii=0; (satest[j] > sastart[SAMAX][j] || /* reject if outside rng */
           satest[j] < sastart[SAMIN][j]); ii++) {
           satest[j] = sa_pe[j] + (rgasdev(sa_rnd)*sstd);
	   // fprintf (stderr,"j %d sstd %g %g %g %g\n",j,sstd, sastart[SARNG][j],  abs(sa_pe[j]), satest[j]);
	   // fprintf (stderr,"j %d satest %g ii %d\n",j, satest[j],ii);

           if (ii>1e5) {			/* trouble finding OK val */
             satest[j] = sa_pe[j];
             break;
          }
        }  /* for (ii;;) */
      }
      else satest[j] = (sastart[SAMAX][j]+sastart[SAMIN][j])*.5;

    }  /* for (j;;) */

    if (info>=4) {
	fprintf (prout,"\n");
	fprintf (prout," p %-4d ",++titot);
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",satest[j]);
	    // fprintf (prout,"r %-10.5g",sastart[SARNG][j]);
        };
	fprintf (prout,"\n");
    }
    for (j=0; j<nparms; j++) {
	sa_ps[k][j] = satest[j];	/* save params for all the runs */
    }
    x = get_model_data(0,satest,-k);	/* run models in parallel */

  }  /* for (k;;) */
  
  sa_yit = sa_yit2 = sa_yit3 = 1e20;

  for (k=1; k<=iter; k++) {
    sa_y = eval_sim(satest,k,xydata,datasize);	/* get data after model done */
    sa_itot++;
    if (info>=4) fprintf (prout, " %d m %-9.4g\n",sa_itot,sa_y);

    if (sa_y < sa_yit3) {			/* save  3rd best from this round of iters */
       if (sa_y < sa_yit2) {			/* save  2nd best from this round of iters */
          if (sa_y < sa_yit) {			/* save best from this round of iters */
            for (j=0; j<nparms; j++) sa_bit[j]=sa_ps[k][j]; /* save the new one as best */
            sa_yit = sa_y;
            sa_yitit = sa_itot;
          } else {
            for (j=0; j<nparms; j++) sa_bit2[j]=sa_ps[k][j]; /* save the 2nd best */
            sa_yit2 = sa_y;
            sa_yitit2 = sa_itot;
          }
       } else {
           for (j=0; j<nparms; j++) sa_bit3[j]=sa_ps[k][j]; /* save the 2nd best */
           sa_yit3 = sa_y;
           sa_yitit3 = sa_itot;
       }
    }

    // if (sa_y < sa_yit) {			/* save best, 2nd, and 3rd best from this round of iters */
    //   for (j=0; j<nparms; j++) sa_bit[j]=sa_ps[k][j]; /* save the new one */
    //   sa_yit=sa_y;
    //   sa_yitit = sa_itot;
    // }
        
    if (sa_y < sa_ye) {			/* Save for next iter if better than latest */
      for (j=0; j<nparms; j++) sa_pe[j]=sa_ps[k][j]; /* save the new one */
      sa_ye=sa_y;
    }
  }  /* for (k;;) */

  insert_b (sa_yit,  temptr, sa_bit,  nparms, sa_yitit);  /* save the best one in push-down stack */
  // insert_b (sa_yit2, temptr, sa_bit2, nparms, sa_yitit2);  /* save the 2nd best one in push-down stack */
  // insert_b (sa_yit3, temptr, sa_bit3, nparms, sa_yitit3);  /* save the 3rd best one in push-down stack */

//    if (sa_y < sa_yb[0]) {		/* Save the best one in push-down stack. */
//       for (i=SANB-1; i>0; i--) { 
//          for (j=0; j<nparms; j++) {
// 	    sa_pb[i][j] = sa_pb[i-1][j];
//          };
//          sa_tb[i] = sa_tb[i-1];
//          sa_yb[i] = sa_yb[i-1];
//          sa_it[i] = sa_it[i-1];
//       };
//       if (++sannb > SANB) sannb = SANB;
//       sa_tb[0] = temptr;
//       sa_yb[0] = sa_y;
//       sa_it[0] = sa_itot;
//       for (j=0; j<nparms; j++) sa_pb[0][j]=satest[j]; /* save the new one */
       // if (sannb>0) {
       //     for (j=0; j<nparms; j++) sa_pe[j]= (satest[j]*0.99 + sa_pb[1][j]*0.01); /* save average of new,old */
       // }
//    }

}

/*------------------------------------------------------------*/

double saybest=0, saybest_sq=0, saybest_sd=0;
double satbest=0, satbest_sq=0, satbest_sd=0;

void calc_sd (int nparms, int nbest, const char *header)

/* calculate mean and stdev of the near best matches */

{
      int i,j,k,s;
      double n,val,mean;

  for (j=0; j<nparms; j++) {		/* initialize satest[] */
      sabest[j] = 0;
      sabest_sq[j] = 0;
  }
  satbest=saybest=0;
  satbest_sq=saybest_sq=0;

  for (i=0,s=sa_sanb-1; i<=s && i<nbest; i++) {   /* use 5 best matches */
       k = i;
       if (i==s && sa_wavg>0) k=0;	   /*  make weighted avg, i=0 twice */
       for (j=0; j<nparms; j++) {	   /* calculate sum, sumsq */
	       val = sa_pb[k][j];
	       sabest[j]    += val;
	       sabest_sq[j] += val*val;
       }
       val = sa_tb[k];		/* temperature */
       satbest    += val;
       satbest_sq += val*val;
       val = sa_yb[k];		/* best match */
       saybest    += val;
       saybest_sq += val*val;
  }
  n = i;
  if (n<3) return;

  for (j=0; j<nparms; j++) {
        mean = sabest[j] / n;
        sabest_sd[j] = sqrt((sabest_sq[j] - (sabest[j]*mean)) / (n-1));
  }
  mean = satbest / n;
  satbest_sd = sqrt((satbest_sq - (satbest*mean)) / (n-1));
  mean = saybest / n;
  saybest_sd = sqrt((saybest_sq - (saybest*mean)) / (n-1));

  sa_ye = saybest/n;
  for (j=0; j<nparms; j++) {
	  sa_pe[j]         = sabest[j]/n; 		/* new avg of previous best matches  */
     	  sastart[SASD][j] = sabest_sd[j] * sa_sdmul;
  }
  if (info >= 2) {
      if (sa_wavg) fprintf (prout,"\n%s wavg  ", header);
      else         fprintf (prout,"\n%s avg   ", header);
      for (j=0; j<nparms; j++) {
           fprintf (prout,"%-10.5g ",sabest[j]/n);
      };
      fprintf (prout," m %-9.4g  t %-6.3g\n", saybest/n,satbest/n);

      fprintf (prout,"%s sd ",header);
      fprintf (prout,"%-2g ",n);
      for (j=0; j<nparms; j++) {
         fprintf (prout,"%-10.4g ",sabest_sd[j]);
      }
      fprintf (prout," m %-9.4g  t %-7.3g\n", saybest_sd,satbest_sd);
  }
}

/*------------------------------------------------------------*/

void stsrchfit (double (*user_func) (double user_x_point, double *coeff, int model_num),
                int datasize, double *xydata, int nparms, double *coeff, double *coeffc)

/* find best set of matching parameter values, given arrays:

     sastart[NSTV][nparms]     contains starting values and ranges
            satest[nparms]     contains set of params for each simulation run

*/

 {  int i,j,ii,jj,u; 
    int iter, sditer, nfp, nbest;
    int stuck, prevok;
    double temptr, oldy, x;
    double nstuckcrit, tstuckcrit;

  if (info==6) prout=stderr;	/* set output to stderr when using vid */
  else         prout=stdout;	/*  otherwise use stdout */

  if (info >= 2)      fprintf (prout,"# Stochastic Search\n");

  setptr("sa_ilim", &sa_ilim);		/* set pointers for command line switches */
  setptr("sa_st",   &sa_st);
  setptr("sa_ftol", &sa_ftol);
  setptr("sa_mthr", &sa_mthr);
  setptr("sa_tstk", &sa_tstk);
  setptr("sa_nstk", &sa_nstk);
  setptr("sa_itb",  &sa_itb);
  setptr("sa_td",   &sa_td);
  setptr("sa_tj",   &sa_tj);
  setptr("sa_bj",   &sa_bj);
  setptr("sa_ybc",  &sa_ybc);
  setptr("sa_sdi",  &sa_sdi);
  setptr("sa_sanb", &sa_sanb);
  setptr("sa_sdmul",&sa_sdmul);
  setptr("calc_sdi", &calc_sdi);
  setptr("sa_wavg", &sa_wavg);
  setptr("rseed",   &rseed);

  setvariables();

  /* - - - - - - - - - - - - - - - - - */
  							/* set up initial values */
  for (i=0; i<nparms; i++) {
	  
      j = i*2;
      sastart[SAVAL][i] = coeff[i];
      sastart[SAMIN][i] = coeffc[j];
      sastart[SAMAX][i] = coeffc[j + 1];
  };
      saname[0]  = p1; 
      saname[1]  = p2; 
      saname[2]  = p3; 
      saname[3]  = p4; 
      saname[4]  = p5; 
      saname[5]  = p6; 
      saname[6]  = p7; 
      saname[7]  = p8; 
      saname[8]  = p9; 
      saname[9]  = p10; 
      saname[10] = p11; 
      saname[11] = p12; 


  /* - - - - - - - - - - - - - - - - - */

  for (j=0; j<nparms; j++) {
        double swap;
	double mean, rng;
        double samax, samin;

    if (sastart[SAMAX][j] < sastart[SAMIN][j]) {  /* check max min */
             swap = sastart[SAMAX][j];            /*  swap if max is min */
             sastart[SAMAX][j] = sastart[SAMIN][j];
             sastart[SAMIN][j] = swap;
    }; 
    samax  = sastart[SAMAX][j];
    samin  = sastart[SAMIN][j];
     
    if (samax==0 && samin==0) {
        sastart[SASD][j]    = 0; 
        sastart[SARNG][j]   = 0; 
    }
     else {
         // sastart[SARNG][j] =  (sastart[SAMAX][j] - sastart[SAMIN][j]) /  (sastart[SAMAX][j] + sastart[SAMIN][j]);
         sastart[SASD][j] = sastart[SAMAX][j] - sastart[SAMIN][j];
         mean = (samax + samin) * 0.5;
         rng  = (samax - samin) * 0.5;
        if (mean==0) sastart[SARNG][j] = rng / max(abs(samax),abs(samin)); 
        else         sastart[SARNG][j] = rng / mean; 
     }
   sabest[j] = 0;

   if (info>=4) fprintf (prout,"# %-9s stdev %-8.4g nstdev %-8.4g\n",saname[j],sastart[SASD][j],sastart[SARNG][j]);
  }

  for (nfp=j=0; j<nparms; j++) {	/* find how many params are free */
    nfp += (sastart[SARNG][j] > sstdevthr); 
  }
  if (info >=3) fprintf (prout,"# %d free params\n",nfp);

  for (j=0; j<nparms; j++) {
    if (sastart[SAVAL][j] > sastart[SAMAX][j]) 
		sastart[SAVAL][j] = sastart[SAMAX][j]; 
    if (sastart[SAVAL][j] < sastart[SAMIN][j]) 
		sastart[SAVAL][j] = sastart[SAMIN][j]; 
  }
  for (i=0; i<SANB; i++) {		/* copy starting values */
    for (j=0; j<nparms; j++) {
       sa_pb[i][j] = sastart[SAVAL][j];
    }; 
    sa_yb[i] = 1e20;
    sa_tb[i]  = 0;
    sa_it[i]  = 0;
  }
  for (j=0; j<nparms; j++) {		/* set first test to starting values */
       satest[j] = sastart[SAVAL][j];
       sa_pe[j]  = sastart[SAVAL][j];
  }

  if (notinit(sa_st))   sa_st   = 1.0;	  /* starting temperature */
  if (notinit(sa_ilim)) sa_ilim = 5000;   /* maximum # of iterations */
  if (notinit(sa_ftol)) sa_ftol = 1e-26;  /* value of match for stopping */
  if (notinit(sa_mthr)) sa_mthr = sa_ftol*1e1; /*val of match to decr stuckcrit*/
  // if (notinit(sa_tstk)) sa_tstk = .08; 	  /* base level for tstuckcrit */
  if (notinit(sa_tstk)) sa_tstk = .28; 	  /* base level for tstuckcrit */
  if (notinit(sa_nstk)) sa_nstk = 5; 	  /* number crit for stuck */
  // if (notinit(sa_itb))  sa_itb  = 1.8; 	  /* base for iter to the n power */
  if (notinit(sa_itb))  sa_itb  = 1.8; 	  /* base for iter to the n power */
  // if (notinit(sa_td))   sa_td   = .88;	  /* temperature decrement */
  if (notinit(sa_td))   sa_td   = .88;	  /* temperature decrement */
  if (notinit(sa_tj))   sa_tj   = 10;	  /* temperature jump */
  if (notinit(sa_bj))   sa_bj   = .05;	  /* best match jump */
  if (notinit(sa_ybc))  sa_ybc  = .86;	  /* best match criterion */
  if (notinit(sa_sdi))  sa_sdi  = 10;	  /* compute sd criterion */
  if (notinit(sa_sanb))  sa_sanb = 15;	  /* trials to compute sd */
  if (notinit(sa_sdmul))  sa_sdmul = 4;	  /* multiplier for sd */
  if (notinit(calc_sdi)) calc_sdi = 1;	  /* =1 -> calculate s.d. */
  if (notinit(sa_wavg)) sa_wavg = 0;	  /* =1 -> calc weighted avg. */
  if (notinit(rseed))     rseed = 427713; /* rseed for st_srch */

  initrand(sa_rnd,rseed);
  if (info >= 2)      fprintf (prout,"#      rseed %d\n",rseed);

   runstsrch(iter=1, nparms, temptr=0, coeff, xydata, datasize); /* run first with orig params */

  for (j=0; j<nparms; j++) sa_bit[j]=satest[j]; /* save the first one */
  sa_yit=sa_y;

  if (info>=2) {	
	fprintf (prout,"#  i    ");			/* print header */
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10s ",saname[j]);
        };
	fprintf (prout," match        T\n");

	fprintf (prout,"bi %-5d ",sa_yitit);		/* print first iteration */
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",sa_bit[j]);
        }
	fprintf (prout," m %-9.4g  t %-7.3g %d\n",
			sa_yit, temptr, sannb);

	fprintf (prout,"bp %-5g ",sa_it[0]);		/* print current best match */
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",sa_pb[0][j]);
        }
	fprintf (prout," m %-9.4g t %-7.3g %d\n",
			sa_yb[0], temptr, sannb);
  }
  oldy = sa_ye;

  sditer = 0;				/* count for sa_sdi */
  stuck = 0;
  prevok = 0;
  temptr = sa_st;			  /* starting temperature */

  for (u=1; sa_itot < sa_ilim; u) {
    iter = 4 + int(pow(sa_itb,(nfp-1)))+(stuck);
    if (iter >= RITER) iter = RITER-1;

    // fprintf (stderr, "iter %d sa_itb %g nfp %d\n",iter, sa_itb, nfp);

    if (info>=2) {				/* print current best estimate */
	fprintf (prout,"pe %-4d ",sa_itot);
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",sa_pe[j]);
        }
	fprintf (prout," m %-9.4g t %-7.3g\n",
			sa_ye, temptr);
    }
    runstsrch(iter, nparms, temptr, coeff, xydata, datasize); /* run simulation with different */
							      /*  set of params each time */

    if (info>=2) {	
	fprintf (prout,"bi %-5d ",sa_yitit);		/* print best from last iteration */
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",sa_bit[j]);
        }
	fprintf (prout," m %-9.4g  t %-7.3g %d\n",
			sa_yit, temptr, sannb);

	fprintf (prout,"bp %-5g ",sa_it[0]);		/* print current best match */
        for (j=0; j<nparms; j++) {
	    fprintf (prout,"%-10.5g ",sa_pb[0][j]);
        }
	fprintf (prout," m %-9.4g t %-7.3g %d\n",
			sa_yb[0], temptr, sannb);
    }

    tstuckcrit = sa_tstk / (1+(sa_mthr/sa_ye)); 
    temptr *= pow(sa_td,(1+3*temptr/sa_st));	/* temperature decrement */

    if (sa_ye/oldy >= sa_ybc) {	/* we're stuck here */
				/* try looking harder and vary the range */

      if (++stuck > sa_nstk && temptr < tstuckcrit) 
      {
	  // temptr *= sa_tj; 	/* temperature jump */
	  // if (temptr > sa_st) temptr = sa_st;
          if (stuck > sa_nstk*2) {
             if (info>=3) fprintf (prout,"stuck, raising match crit\n");
             sa_ye *= 1 + sa_bj * -log(rrand(sa_rnd));
             oldy = sa_ye;
	     stuck = sa_nstk/2; /* reset stuck so it doesn't make huge number of iters */
          }
          // temptr = sa_st;
	  if (++prevok >= sannb) prevok = sannb;
          temptr = sa_tb[prevok] * sa_tj;		/* temperature jump */
	  if (temptr > sa_st) temptr = sa_st; 
          if (info>=3) fprintf (prout,"stuck, increasing temp to %g\n",temptr);
          sa_ye = sa_yb[prevok];
          for (j=0; j<nparms; j++) sa_pe[j] = sa_pb[prevok][j]; /* restore the previous best one */
      }

    } else {
      stuck = 0;
      prevok = 0;
      oldy = sa_ye;
      //fprintf (prout,"tch %g match %g\n",temptr, sa_ye);
    };

  if (info>=3) fprintf (prout," stuck %d tstkcrit %g oldyb %g\n",
			stuck, tstuckcrit,oldy);

  for (i=nbest=0; i<sannb; i++) {
      if (sa_yb[i] < sa_yb[0] * 10.0) {	/* find near-best matches */
	  nbest++;
      }
  }
  if (calc_sdi>0) {			/* calculate mean, sd every sa_sdi iterations */
      if (sditer++ >= sa_sdi) {
          sditer = 0;
          if (nbest>1) {
	      fprintf (prout,"bi  \n");
	      calc_sd(nparms,nbest,"bi"); 
              // temptr /= pow(sa_td,(1+3*temptr/sa_st));	/* temperature inecrement */
              temptr = sa_st; 				/* temperature inecrement */
          }
      }
  }
  if (sa_ye < sa_ftol) break;

 }  /* for (;;) */

 if (info>=1) {					 /* When done, print out final answer */
       int npr;

     if (info >=2) {
       npr = sannb;			/* number to print */
       fprintf (prout,"\n");
       fprintf (prout,"Summary of best values:\n");
     }
     else          npr = 1;

     for (i=npr-1; i>=0; i--) {			/* print best values */
        fprintf (prout,"run ");
        fprintf (prout,"%-5g ",sa_it[i]);
        for (j=0; j<nparms; j++) {
           fprintf (prout,"%-10.5g ",sa_pb[i][j]);
        };
        fprintf (prout,"m %-9.4g  t %-6.3g t/m %-9.3g\n", sa_yb[i], sa_tb[i],
					 sa_tb[i]/sa_yb[i]);
     }

     if (nbest>1) {
	calc_sd(nparms,nbest,""); 
     } 

     for (j=0; j<nparms; j++) sa_tsav[j] = sa_pb[1][j]; /* save the next-to-best one */

     /* run the best match once more to leave data file */

     for (j=0; j<nparms; j++) sa_pe[j] = sa_pb[0][j]; /* restore the previous best one */
     runstsrch(iter=1, nparms, temptr=0, coeff, xydata, datasize); /* run first with orig params */

     fprintf (prout,"\n    ");
     fprintf (prout,"%-4d ",sa_itot);
     for (j=0; j<nparms; j++) {
          fprintf (prout,"%-10.5g ",satest[j]);
     };
     fprintf (prout,"m %-9.4g  t %-6.3g t/m %-9.3g\n", sa_y, temptr, temptr/sa_y);

 }
  else if (info==0) {
     fprintf (prout,"%-9.5g ",sa_yb[0]);
 }


}

/*------------------------------------------------------------*/

