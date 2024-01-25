/* Experiment gcoff_hbp_flicker */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"

double temp_freq;
double sdia;
double stimdur;
double stimtime;
double inten;
double minten;
double scontrast;
double setploti;
double set_conehbpdiv;
double exptdur;
double postdur;

double set_shoffs;
double set_shgain;
double set_shdur;

int vcmode;
int ntrials;

int rec_ct;
int rec_cn;

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("exptdur",   &exptdur);
  setptr("sdia",      &sdia);
  setptr("stimdur",   &stimdur);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("setploti",  &setploti);
  setptr("ntrials",   &ntrials);
  setptr("vcmode",    &vcmode);
  setptr("set_conehbpdiv",  &set_conehbpdiv);
  setptr("set_shoffs",  &set_shoffs);
  setptr("set_shgain",  &set_shgain);
  setptr("set_shdur",  &set_shdur);
  nvalfile = "nval_gcoff_hbp_flashes.n";
}

/*------------------------------------------------------*/

void setparams(void)
{
   double div;

  make_rods = 0;
  make_cones= 1;        /* make cones, hbp, gcoff */
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 0;
  make_dbp2 = 0;
  make_hbp1 = 1;
  make_hbp2 = 0;
  make_rbp  = 0;
  make_gca  = 0;
  make_gcb  = 0;
  make_gcaoff=1;
  make_gcboff=1;
  make_dsgc = 0;

  if(notinit(rec_ct)) rec_ct = gcaoff;
  //if (notinit(arrsiz)) arrsiz = 300;
  if (notinit(bg_inten)) bg_inten = 1.8e4;      /* background light intensity */
  if (!notinit(set_conehbpdiv)) {
	  div = getn(xcone,CELDIV2); 	      /* cone -> hbp div */
	  if (set_conehbpdiv==0) set_conehbpdiv = div;
	  setn(xcone,CELDIV2,set_conehbpdiv); /* cone -> hbp div */
	  if (set_conehbpdiv==1) setn(xcone,SCOND2,getn(xcone,SCOND2)*1.5); /* cone -> hbp div */
	  else setn(xcone,SCOND2,getn(xcone,SCOND2)*div/set_conehbpdiv); /* cone -> hbp div */
  }
  if (!notinit(set_shoffs)) setn(hbp1,SHOFFS1,set_shoffs);   /* offset for high pass filter */
  if (!notinit(set_shgain)) setn(hbp1,SHGAIN1,set_shgain);   /* gain for high pass filter */
  if (!notinit(set_shdur))  setn(hbp1,SDURH1, set_shdur);    /* time constant for high pass filter */

  srseed =  2000 + 1000 * scontrast;		/* make synaptic/photon/dark noise vary with stimulus */
  if (ninfo>=2) fprintf (stderr,"# expt_gcoff_hbp_flicker: celdiv %g scond %g\n",getn(xcone,CELDIV2),getn(xcone,SCOND2));
}

/*------------------------------------------------------*/

double light_plot(double x, double y)
{
   return (minten*(1+inten));
}

/*------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, ct2, cn2, n, r, plnum;
    int colr,synin1;
    int midcone, midcone2, midcone3, midcone4;
    int midhbp, midhbp2, midhbp3, midhbp4;
    int midhbp5, midhbp6, midhbp7, midhbp8;
    double m, t, fmax, fmin, tim, start, dur, dtrial;
    double lmin, lmax, rmin, rmax, plsize;
    double Vmin, Vmax;
    double Imin, Imax;
    char pbuf[100];
#define GAUSRNG 1

  if (notinit(vcmode))       vcmode = 1;         /* voltage clamp */
  if (notinit(temp_freq)) temp_freq = 1;         /* temporal frequency */
  if (notinit(sdia))           sdia = 500;       /* spot diameter */
  if (notinit(stimdur))     stimdur = 0.1;       /* stimulus duration */
  if (notinit(stimtime))   stimtime = 0.1;       /* stimulus time */
  if (notinit(postdur))     postdur = 0;        /* time after stimulus */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0;      /* intensity increment */
  if (!notinit(setploti))     ploti = setploti;  /* plot time increment */
  if (notinit(ntrials))     ntrials = 1;         /* number of flashes */

  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  dtrial = 1 / temp_freq;
  if (notinit(exptdur)) exptdur  = dtrial * ntrials;
  dtrial = exptdur / ntrials;
  endexp  = exptdur + postdur; 

  midcone  = findmid(xcone,0,0);
  midcone2 = findmid(xcone,20,0);
  midcone3 = findmid(xcone,0,10);
  midcone4 = findmid(xcone,10,10);
  midhbp   = findmid(hbp1,  0,0);
  midhbp2  = findmid(hbp1, 20,10);
  midhbp3  = findmid(hbp1, 20,20);
  midhbp4  = findmid(hbp1,-10,10);
  midhbp5  = findmid(hbp1,-20,10);
  midhbp6  = findmid(hbp1,-20,-10);
  midhbp6  = findmid(hbp1,-10,-20);
  midhbp8  = findmid(hbp1,-20,-30);
  synin1 = ncel_in(hbp1,midhbp, xcone);

  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid hbp  # %d ncones %d\n",  midhbp, synin1);

  sprintf (pbuf,"P cone %d",midcone);
  plot_l_nod(ct=xcone,cn=midcone,n=soma,lmin=0,lmax =1e4,colr=magenta,pbuf, 100, 0.5); /* plot Lcones*/
  // plot_func((double(*)(double, double))light_plot,n,fmax=5e4,fmin=0); plot_param(pbuf,colr=magenta,100,0.5);

  plot_v_nod(ct=xcone,cn=midcone,n=soma,Vmin=-.028,Vmax =-.023,colr=cyan,"", -1, -1); /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone2,n=soma,Vmin=-.028,Vmax =-.023,colr=red,"", -1, -1); /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone3,n=soma,Vmin=-.028,Vmax =-.023,colr=green,"", -1, -1); /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone4,n=soma,Vmin=-.028,Vmax =-.023,colr=blue,"", -1, -1); /* plot Vcones*/
  plot_synrate_out(ct=xcone,cn=midcone,rmin=0,rmax=400,colr=magenta);	              /* plot rate out */
  plot_v_nod(ct=hbp1,cn=midhbp,n=soma, Vmin=-.045,Vmax =-.035,colr=red,"", -1, -1);     /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp2,n=soma,Vmin=-.045,Vmax =-.035,colr=green,"", -1, -1);   /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp3,n=soma,Vmin=-.045,Vmax =-.035,colr=blue,"", -1, -1);    /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp4,n=soma,Vmin=-.045,Vmax =-.035,colr=magenta,"", -1, -1); /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp5,n=soma,Vmin=-.045,Vmax =-.035,colr=ltred,"", -1, -1);   /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp6,n=soma,Vmin=-.045,Vmax =-.035,colr=ltgreen,"", -1, -1); /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp7,n=soma,Vmin=-.045,Vmax =-.035,colr=ltblue,"", -1, -1);  /* plot Vhbp */
  plot_v_nod(ct=hbp1,cn=midhbp8,n=soma,Vmin=-.045,Vmax =-.035,colr=ltmag,"", -1, -1);   /* plot Vhbp */
  plot_synrate_out(ct=hbp1,cn=midhbp,ct2=gcaoff,cn2=1,rmin=0,rmax=200,colr=blue,1);     /* plot rate out */
  plot_synrate_out(ct=hbp1,cn=midhbp,ct2=gcaoff,cn2=2,rmin=0,rmax=200,colr=magenta,0);  /* plot rate out */
  if (vcmode) {
     plot_i_nod(ct=gcaoff, cn=1,n=soma,Imin=-1e-10,Imax =1e-10,colr=blue,"", -1, -1);      /* plot Vgc */
     plot_i_nod(ct=gcaoff, cn=2,n=soma,Imin=-1e-10,Imax =1e-10,colr=red,"", -1, -1);       /* plot Vgc */
  }
  else {
     plot_v_nod(ct=gcaoff, cn=1,n=soma,Vmin=-.070,Vmax =-.050,colr=blue,"", -1, -1);      /* plot Vgc */
     plot_v_nod(ct=gcaoff, cn=2,n=soma,Vmin=-.070,Vmax =-.050,colr=red,"", -1, -1);       /* plot Vgc */
  }
  if (getn(gcaoff,BIOPHYS)) {plot(CA, 1, ndn(gcaoff,1,soma), fmax=0.5e-6, fmin=0); 
			plot_param("Cai", colr=yellow,plnum=0,plsize=0.3);}

  stim_backgr(minten);

  if (vcmode) {
     vclamp (ndn(ct=gcaoff,cn=1,n=soma), -0.058, 0, endexp);
     vclamp (ndn(ct=gcaoff,cn=2,n=soma), -0.058, 0, endexp);
  }

  for (t=r=0; r<ntrials; r++,t+=dtrial){
      stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=stimdur);
      step(dtrial);
  }
}
