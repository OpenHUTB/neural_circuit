/* Experiment gc_cbp_aii_flash */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"

double temp_freq;
double exptdur;
double ntrials;
double stimdur;
double sdia;
double stimtime;
double stimamp;
double minten;
double scontrast;
double setploti;
double set_conecbpdiv;
double set_conecbpcond;

double set_shoffs;
double set_shgain;
double set_shdur;

double set_cbpcbp;
double set_cbpaii;
double set_aiiaii;
int set_sf1h;

int vcmode;

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("exptdur",   &exptdur);
  setptr("stimdur",   &stimdur);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("stimamp",   &stimamp);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("setploti",  &setploti);
  setptr("vcmode",    &vcmode);
  setptr("set_conecbpdiv",  &set_conecbpdiv);
  setptr("set_conecbpcond",  &set_conecbpcond);
  setptr("set_shoffs",  &set_shoffs);
  setptr("set_shgain",  &set_shgain);
  setptr("set_shdur",   &set_shdur);
  setptr("set_cbpcbp",  &set_cbpcbp);
  setptr("set_cbpaii",  &set_cbpaii);
  setptr("set_aiiaii",  &set_aiiaii);
  setptr("set_sf1h",    &set_sf1h);
  nvalfile = "nval_gc_cbp_aii_flash.n";
}

/*------------------------------------------------------*/

void setparams(void)
{
   double div;

  make_rods = 0;
  make_cones= 1;        /* make cones, cbp, gc */
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 1;
  make_dbp2 = 1;
  make_hbp1 = 0;
  make_hbp2 = 0;
  make_rbp  = 0;
  make_aii  = 1;
  make_gca  = 1;
  make_gcb  = 0;
  make_gcaoff=0;
  make_gcboff=0;
  make_dsgc = 0;

  make_dbp1_dbp1 = 1;
  make_dbp1_aii = 1;
  make_aii_aii = 1;

  set_conecbpdiv = 1;

  //if (notinit(arrsiz)) arrsiz = 300;
  if (notinit(bg_inten)) bg_inten = 1.8e4;      /* background light intensity */
  if (!notinit(set_conecbpcond))   setn(xcone,SCOND1,set_conecbpcond); /* cone->cbp cond */
  if (!notinit(set_conecbpdiv)) {
	  div = getn(xcone,CELDIV1); 	      /* cone -> cbp div */
	  if (set_conecbpdiv==0) set_conecbpdiv = div;
	  setn(xcone,CELDIV1,set_conecbpdiv); /* cone -> cbp div */
	  if (set_conecbpdiv==1) setn(xcone,SCOND1,getn(xcone,SCOND1)*1.5); /* cone -> cbp div */
	  else setn(xcone,SCOND1,getn(xcone,SCOND1)*div/set_conecbpdiv); /* cone -> cbp div */
  }
  if (!notinit(set_cbpcbp)) setn(dbp1,SCOND5,set_cbpcbp); /* cbp -> cbp gj cond */
  if (!notinit(set_cbpaii)) setn(dbp1,SCOND6,set_cbpaii); /* cbp -> aii gj cond */
  if (!notinit(set_aiiaii)) setn(aii,SCOND1,set_aiiaii); /* aii -> aii gj cond */
  if (!notinit(set_sf1h))   setn(dbp1,SNFILTH1,set_sf1h); /* synap high-pass filter */

  srseed =  2000 + 1000 * scontrast;		/* make synaptic/photon/dark noise vary with stimulus */
  //if (ninfo>=2) fprintf (stderr,"# expt_gc_cbp_aii_flash: celdiv %g scond %g\n",getn(xcone,CELDIV2),getn(xcone,SCOND2));
}

/*------------------------------------------------------*/

void stim_pulse (int ct, int cn, double i, double ctime, double dur)

{
   node *nd;

  nd = ndn(ct,cn,soma);
  cclamp (nd, i, ctime, dur);
}

/*------------------------------------------------------*/

void runexpt(void)

{
#define NUMCBPS 8
#define NUMAIIS 8
#define STIMRATE 10
    int i, j, ct, ct2, cn, cn2, n, plnum;
    int colr,synin1;
    int midcone, midcone2, midcone3, midcone4;
    int cbps[NUMCBPS];
    int aiis[NUMAIIS];
    double di, t, fmax, fmin, lmin, lmax;
    double rmin, rmax, plsize;
    double start, dur;
    double dtrial;
    double Vmin, Vmax;
    double Imin, Imax;
    double vcdur;
    char pbuf[100];

  if (notinit(vcmode))       vcmode = 1;         /* voltage clamp */
  if (notinit(temp_freq)) temp_freq = 2;
  if (notinit(sdia))           sdia = 300;       /* spot diameter */
  if (notinit(stimdur))     stimdur = .1;        /* stimulus duration */
  if (notinit(stimtime))   stimtime = .10;       /* stimulus time */
  if (notinit(stimamp))     stimamp =  0;        /* current stimulus amplitude */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = .5;        /* intensity increment */
  if (!notinit(setploti))     ploti = setploti;  /* plot time increment */
  if (notinit(ntrials))     ntrials = 10;

  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  dtrial = 1;
  if (notinit(exptdur)) exptdur = dtrial * ntrials;
  dtrial = exptdur / ntrials;
  endexp  = exptdur;

  midcone  = findmid(xcone,0,0);
  midcone2 = findmid(xcone,20,0);
  midcone3 = findmid(xcone,0,10);
  midcone4 = findmid(xcone,10,10);
  for (i=0; i<NUMCBPS; i++) {
     cbps[i] = 0;
  }
  cbps[0]  = findmid(dbp1,  0, 0);
  cbps[1]  = findmid(dbp1, 20, 20);
  cbps[2]  = findmid(dbp1, 20,-20);
  cbps[3]  = findmid(dbp1,-20, 20);
  cbps[4]  = findmid(dbp1,-20,-20);
  cbps[5]  = findmid(dbp1,-50, 20);
  cbps[6]  = findmid(dbp1,-50,-50);
  cbps[7]  = findmid(dbp1, 50,-20);

  aiis[0]  = findmid(aii,  0, 0);
  aiis[1]  = findmid(aii, 20, 20);
  aiis[2]  = findmid(aii, 20,-20);
  aiis[3]  = findmid(aii,-20, 20);
  aiis[4]  = findmid(aii,-20,-20);
  aiis[5]  = findmid(aii,-50, 20);
  aiis[6]  = findmid(aii,-50,-50);
  aiis[7]  = findmid(aii, 50,-20);

  synin1 = ncel_in(dbp1,cbps[0], xcone);

  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d ncones %d\n",  cbps[0], synin1);

  sprintf (pbuf,"P cone %d",midcone);
  plot_l_nod(ct=xcone,cn=midcone,n=soma,lmin=0,lmax=2e4,colr=magenta,pbuf, 100, 0.4); /* plot Lcones*/
    
  plot_v_nod(ct=xcone,cn=midcone, n=soma,Vmin=-.026,Vmax =-.023,colr=cyan,"", -1, 0.6); /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone2,n=soma,Vmin=-.026,Vmax =-.023,colr=red,"", -1, 0.6);  /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone3,n=soma,Vmin=-.026,Vmax =-.023,colr=green,"", -1,0.6); /* plot Vcones*/
  plot_v_nod(ct=xcone,cn=midcone4,n=soma,Vmin=-.026,Vmax =-.023,colr=blue,"", -1, 0.6); /* plot Vcones*/
  plot_synrate_out(ct=xcone,cn=midcone,rmin=0,rmax=300,colr=magenta,0.6);                /* plot rate out */

  plot_v_nod(ct=dbp1,cn=cbps[0],n=soma,Vmin=-.047,Vmax =-.035,colr=red,"", -1, -1);     /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[1],n=soma,Vmin=-.047,Vmax =-.035,colr=green,"", -1, -1);   /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[2],n=soma,Vmin=-.047,Vmax =-.035,colr=blue,"", -1, -1);    /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[3],n=soma,Vmin=-.047,Vmax =-.035,colr=magenta,"", -1, -1); /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[4],n=soma,Vmin=-.047,Vmax =-.035,colr=ltred,"", -1, -1);   /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[5],n=soma,Vmin=-.047,Vmax =-.035,colr=ltgreen,"", -1, -1); /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[6],n=soma,Vmin=-.047,Vmax =-.035,colr=ltblue,"", -1, -1);  /* plot Vcbp */
  plot_v_nod(ct=dbp1,cn=cbps[7],n=soma,Vmin=-.047,Vmax =-.035,colr=ltmag,"", -1, -1);   /* plot Vcbp */

  plot_v_nod(ct=aii,cn=aiis[0],n=soma,Vmin=-.05,Vmax =-.035,colr=red,"",   -1, -1);   /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[1],n=soma,Vmin=-.05,Vmax =-.035,colr=green,"", -1, -1);   /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[2],n=soma,Vmin=-.05,Vmax =-.035,colr=blue,"",  -1, -1);   /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[3],n=soma,Vmin=-.05,Vmax =-.035,colr=magenta,"",-1, -1);  /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[4],n=soma,Vmin=-.05,Vmax =-.035,colr=ltred,"", -1, -1);   /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[5],n=soma,Vmin=-.05,Vmax =-.035,colr=ltgreen,"",-1, -1);  /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[6],n=soma,Vmin=-.05,Vmax =-.035,colr=ltblue,"",-1, -1);   /* plot Vaii */
  plot_v_nod(ct=aii,cn=aiis[7],n=soma,Vmin=-.05,Vmax =-.035,colr=ltmag,"", -1, -1);   /* plot Vaii */

  plot_synrate_out(ct=dbp1,cn=cbps[0],ct2=gca,cn2=1,rmin=0,rmax=200,colr=blue,1);     /* plot rate out */
  plot_synrate_out(ct=dbp1,cn=cbps[0],ct2=gca,cn2=2,rmin=0,rmax=200,colr=magenta,0);  /* plot rate out */


  if (vcmode) {
     if (nde(gca,1,soma))
	plot_i_nod(ct=gca, cn=1,n=soma,Imin=-1e-10,Imax =1e-10,colr=blue,"", -1, -1);      /* plot Vgc */
     if (nde(gca,2,soma))
	plot_i_nod(ct=gca, cn=2,n=soma,Imin=-1e-10,Imax =1e-10,colr=red,"", -1, -1);       /* plot Vgc */
  }
  else {
     if (nde(gca,1,soma))
	plot_v_nod(ct=gca, cn=1,n=soma,Vmin=-.070,Vmax =-.050,colr=blue,"", -1, -1);      /* plot Vgc */
     if (nde(gca,2,soma))
	plot_v_nod(ct=gca, cn=2,n=soma,Vmin=-.070,Vmax =-.050,colr=red,"", -1, -1);       /* plot Vgc */
  }
  if (nde(gca,1,soma)&&getn(gca,BIOPHYS)) {plot(CA, 1, ndn(gca,1,soma), fmax=0.5e-6, fmin=0); 
			plot_param("Cai", colr=yellow,plnum=0,plsize=0.3);}
  if (vcmode) {
     if (nde(gca,1,soma)) vclamp (ndn(ct=gca,cn=1,n=soma), -0.058, 0, endexp);
     if (nde(gca,2,soma)) vclamp (ndn(ct=gca,cn=2,n=soma), -0.058, 0, endexp);
  }
  
  
  stim_backgr(minten);
  stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=stimdur);

  for (t=0; t<exptdur; t+= dtrial){
     double start, dur;

    for (i=0; i<NUMAIIS; i++) {
      for (j=0; j<STIMRATE; j++) {
        stim_pulse (ct=aii, cn=aiis[i], di=stimamp, simtime+0.005+j*0.1+rrange(0,0.095), vcdur=0.005);
      }
    }
    step(dtrial);
  }
}
