/* Experiment aii_flash */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"


double temp_freq;
double dstim;
double sdia;
double stimtime;
double minten;
double scontrast;
int ntrials;

int rec_ct;
int rec_cn;

void defparams(void)
{
  setptr("rec_ct",    &rec_ct);
  setptr("rec_cn",    &rec_cn);
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  nvalfile = "nval_aii_flash.n";
}

void setparams(void)

{
  make_rods  = 1;        /* make rods, rbc, aii */
  make_cones = 0;
  make_ha    = 0;
  make_hb    = 0;
  make_hbat  = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_rbp   = 1;
  make_aii   = 1;
  make_gca   = 1;
  make_gcb   = 0;

  make_dbp1_aii  = 1;
  make_aii_aii   = 1;

  if(notinit(rec_ct)) rec_ct = aii;
  if(notinit(rec_cn)) rec_cn = 1;
  if (notinit(arrsiz)) arrsiz = 60;
  if (notinit(bg_inten)) bg_inten = 10;
}

/*------------------------------------------------*/

void runexpt(void)

{
   int c, pl, midrod,plnum;
   double s, t, start, dur, dtrial;
   double max,min;
   double plsize,predur;

  if (notinit(ntrials))     ntrials = 1;	/* number of trials */
  if (notinit(dstim))         dstim = .01;	/* stimulus duration */
  if (notinit(sdia))           sdia = 300;	/* spot diameter */
  if (notinit(stimtime))   stimtime = .05;	/* stimulus time */
  if (notinit(minten))       minten = bg_inten;	/* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0.99;	/* intensity increment */

  if (make_rods) {
    plot (V, ndn(xrod,findmid(xrod,-2,0), soma),max=-0.034,min=-0.039); plot_param ("Rod1", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod, 0,0), soma),max=-0.034,min=-0.039); plot_param ("Rod2", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod, 2,0), soma),max=-0.034,min=-0.039); plot_param ("Rod3", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod, 4,0), soma),max=-0.034,min=-0.039); plot_param ("Rod4", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod, 0,2), soma),max=-0.034,min=-0.039); plot_param ("Rod5", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod, 0,4), soma),max=-0.034,min=-0.039); plot_param ("Rod6", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod,-2,4), soma),max=-0.034,min=-0.039); plot_param ("Rod7", plnum=11, plsize=1);
    plot (V, ndn(xrod,findmid(xrod,-2,-4),soma),max=-0.034,min=-0.039); plot_param ("Rod8", plnum=11, plsize=1);
  }
  if (make_hbat) {plot (V, ndn(hbat,1,soma),max=-.02,min-=0.04); plot_param("hbat",c=2,pl=5,s=1);}
  if (make_rbp)  {plot (V, ndn(rbp,findmid(rbp,0,0),soma),max=-0.02, min=-0.065); plot_param("rbp",c=6,pl=6,s=1);}
  if (make_aii)  {plot (V, ndn(aii,findmid(aii,0,0),soma),max=-0.04, min=-0.07); plot_param("aii",c=14,pl=5,s=1);}
  /* plot (FA9, fwsynat[midrod+2] max 200 min 0);  */
  /* plot (FB1, fwsynrb[midrod] max 10 min 0); */
  /* plot (FC1, fwsynrb[midrod] max 10 min 0); */

  //stim_backgr(minten);
  t = 0;
  dtrial = 0.5; 
  stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=dstim);
  endexp = dtrial*ntrials;

  setxmin = 0;
  predur = 0.02;
  simtime = 0-predur;
  step (predur);

  run();
}
