/* Experiment gc_cbp_flash */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

double temp_freq;
double ntrials;
double dstim;
double sdia;
double stimtime;
double minten;
double scontrast;
double setploti;

int rec_ct;
int rec_cn;

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("setploti",  &setploti);
  nvalfile = "nval_gc_cbp_am_flash.n";
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  make_cones= 1;        /* make cones, dbp, gc */
  make_ha   = 1;
  make_hb   = 1;
  make_hbat = 0;
  make_dbp1 = 1;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_am   = 1;
  make_gca  = 1;
  make_gcb  = 0;
  make_dsgc = 0;

  make_hb_dbp1 = 0;

  if(notinit(rec_ct)) rec_ct = gca;
  //if (notinit(arrsiz)) arrsiz = 300;
  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */
  //if (arrsiz==100) {
  //  setsv (dbp1,SCOND,1, 25e-10);
  // } 
 }

/*------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, plnum;
    int colr;
    int midcone, midcbp, midcbp2;
    int synin1, synin2;
    double t, fmax,fmin;
    double rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax;

  if (notinit(temp_freq)) temp_freq = 2;
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  dtrial = 1 / temp_freq;
  exptdur = dtrial * ntrials;
  endexp  = exptdur;
  ploti = 1e-4;

  if (notinit(dstim))         dstim = .05;      /* stimulus duration */
  if (notinit(sdia))           sdia = 300;      /* spot diameter */
  if (notinit(stimtime))   stimtime = .10;      /* stimulus time */
  if (notinit(minten))       minten = bg_inten; /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = .5;       /* intensity increment */
  if (!notinit(setploti))     ploti = setploti; /* plot time increment */

  midcone = findmid(xcone,0,0);
  midcbp  = findmid(dbp1,0,0);
  midcbp2 = findmid(dbp1,10,10);
  //midcbp2 = find_gtconn(dbp1, 8);
  synin1 = ncel_in(dbp1,midcbp,xcone);
  synin2 = ncel_in(dbp1,midcbp2,xcone);
  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d ncones %d\n",  midcbp,synin1);
  if (ninfo >=1) fprintf (stderr,"# mid cbp2 # %d ncones %d\n",  midcbp2,synin2);


  plot_v_nod(ct=xcone,cn=midcone,n=soma,Vmin=-.030,Vmax =-.025,colr=cyan,"", -1, -1); /* plot Vcones*/
  plot_synrate_out(ct=xcone,cn=midcone,rmin=0,rmax=400,colr=magenta);	              /* plot rate out */
  plot_v_nod(ct=dbp1,cn=midcbp,n=soma, Vmin=-.045,Vmax =-.035,colr=red,"", -1, -1);     /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=midcbp2,n=soma,Vmin=-.045,Vmax =-.035,colr=green,"", -1, -1);  /* plot Vcbp */
  plot_synrate_out(ct=dbp1,cn=midcbp,rmin=0,rmax=200,colr=magenta);	              /* plot rate out */
  plot_v_nod(ct=gca, cn=1,n=soma,Vmin=-.070,Vmax =-.050,colr=blue,"", -1, -1);	      /* plot Vgc */
  if (make_gca && getn(gca,BIOPHYS)) {plot(CA, 1, ndn(gca,1,soma), fmax=0.5e-6, fmin=0); 
			plot_param("Cai", colr=yellow,plnum=0,plsize=0.3);}

  stim_backgr(minten);
  for (t=0; t<exptdur; t+= dtrial){
     double start, dur;

    stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=dstim);
    step(dtrial);
  }
}
