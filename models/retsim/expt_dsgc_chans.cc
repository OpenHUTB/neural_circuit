/* Experiment dsgc_sbac for retsim */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"

#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int n_dsgc;
extern int rec_ct;

double theta;
double rectheta;
double iroff;
int light_inhib;
int dsgc_prefdir;
int node_dist;

int no_inhib;
int sbarr;
int rec_ct;
int rec_cn;
int set_vclamp;
int ivplot;
int elnode;
int direction=0;
int outward=0;
int skipv;
int capfilt;
int revdir;
int set_kchan;

int rev_inhib;
int ams1;
int ams2;
int ams3;
int ams4;
int ams5;
int ams6;
int ams7;
int ams8;
int ams9;
int ams10;

int dsgc_n1;    /* CF dsgc nodes selected for plots on command line */
int dsgc_n2;
int dsgc_n3;
int dsgc_n4;

double ampa_cond;
double nmda_cond;
double gaba_cond;

int    stimtype;
int    sqwave;
double barwidth;
double barlength;
double minten;
double econtrast;
double icontrast;
double scontrast;
double eincr;
double iincr;
double sincr;
double velocity;
double stimx;
double stimy;
double stimr;
double stimscale;

double dbp1_anpi;
double dbp1_anpo;
double dbp2_anpi;
double dbp2_anpo;

double vstart;
double vstop;
double vstep;
double vhold;
double tailvolt;
double gvrev;

double disptime;
double stimtime;
double stimdur;
double prestimdur;
double poststimdur;
double tailcurdur;
double spotdur;
double sblur;
double ioffset;
double istim;
double predur;
double set_drm;
double elec_cap;
double elec_rs;
double shunt_res;
double shunt_v;
double dendrm;
double dendcm;
double somacm;
double tipcap;
double dendrs;
double cplam;
double arb_scale;
double kdr_voff;
double setendexp;   // set to shorten length of simulation 
double kexp;        // activation exponent of K chan 
double soma_z;       // soma z dist
double filt_tau;     // filter tau (not used) 
double filt_cfreq;   // filter cutoff freq 
double dsomadia;       // soma dia
double btrend;         // linear slope trend in baseline
double ntrend;         // sine trend in baseline
double nperiod;         // sine phase for trend in baseline
double nphase;         // sine phase for trend in baseline

double speriod;
double sphase;
double orient;
double tfreq;
double drift;

double reczmax;
double reczmin;

double set_timinc;
double set_crit;

char savefile[30] = {0};

void sb_init(void);

/*--------------------------------------------------------*/

void defparams(void)

{
  if (make_movie) {
     defparams_dsgc_movie();
     defparams_onplot_movie();
  }
  setptr("theta", 	 &theta);
  setptr("rectheta", 	 &rectheta);
  setptr("iroff", 	 &iroff);
  setptr("light_inhib",  &light_inhib);
  setptr("dsgc_prefdir", &dsgc_prefdir);
  setptr("node_dist",    &node_dist);
  setptr("no_inhib",	 &no_inhib);
  setptr("ampa_cond",	 &ampa_cond);
  setptr("nmda_cond",	 &nmda_cond);
  setptr("gaba_cond",	 &gaba_cond);
  setptr("rec_ct",	 &rec_ct);
  setptr("rec_cn",	 &rec_cn);
  setptr("sbarr",	 &sbarr);
  setptr("set_vclamp",	 &set_vclamp);
  setptr("ivplot",	 &ivplot);
  setptr("elnode",	 &elnode);
  setptr("skipv",	 &skipv);
  setptr("capfilt",	 &capfilt);
  setptr("revdir",	 &revdir);
  setptr("set_kchan",	 &set_kchan);

  setptr("rev_inhib",    &rev_inhib);
  setptr("ams1",         &ams1);
  setptr("ams2",         &ams2);
  setptr("ams3",         &ams3);
  setptr("ams4",         &ams4);
  setptr("ams5",         &ams5);
  setptr("ams6",         &ams6);
  setptr("ams7",         &ams7);
  setptr("ams8",         &ams8);
  setptr("ams9",         &ams9);
  setptr("ams10",        &ams10);
  setptr("dsgc_n1",      &dsgc_n1);
  setptr("dsgc_n2",      &dsgc_n2);
  setptr("dsgc_n3",      &dsgc_n3);

  setptr("vstart",	 &vstart);
  setptr("vstop",	 &vstop);
  setptr("vstep",	 &vstep);
  setptr("vhold",	 &vhold);
  setptr("tailvolt",	 &tailvolt);
  setptr("gvrev",	 &gvrev);

  setptr("stimtype",	 &stimtype);
  setptr("sqwave",	 &sqwave);
  setptr("barwidth",   &barwidth);
  setptr("barlength",  &barlength);
  setptr("minten",     &minten);
  setptr("econtrast",  &econtrast);
  setptr("icontrast",  &icontrast);
  setptr("scontrast",  &scontrast);
  setptr("eincr",      &eincr);
  setptr("iincr",      &iincr);
  setptr("sincr",      &sincr);
  setptr("velocity",   &velocity);
  setptr("stimx",      &stimx);
  setptr("stimy",      &stimy);
  setptr("stimr",      &stimr);
  setptr("stimscale",  &stimscale);
  
  setptr("dbp1_anpi",  &dbp1_anpi);
  setptr("dbp1_anpo",  &dbp1_anpo);
  setptr("dbp2_anpi",  &dbp2_anpi);
  setptr("dbp2_anpo",  &dbp2_anpo);

  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("stimtime",   &stimtime);
  setptr("stimdur",    &stimdur);
  setptr("prestimdur", &prestimdur);
  setptr("poststimdur",&poststimdur);
  setptr("tailcurdur", &tailcurdur);
  setptr("spotdur",    &spotdur);
  setptr("sblur",      &sblur);
  setptr("predur",     &predur);
  setptr("set_drm",    &set_drm);
  setptr("elec_cap",   &elec_cap);
  setptr("elec_rs",    &elec_rs);
  setptr("shunt_res",  &shunt_res);
  setptr("shunt_v",    &shunt_v);
  setptr("dendrm",     &dendrm);
  setptr("dendcm",     &dendcm);
  setptr("somacm",     &somacm);
  setptr("tipcap",     &tipcap);
  setptr("dendrs",     &dendrs);
  setptr("arb_scale",  &arb_scale);
  setptr("kdr_voff",   &kdr_voff);
  setptr("setendexp",  &setendexp);
  setptr("kexp",       &kexp);
  setptr("soma_z",     &soma_z);
  setptr("filt_tau",   &filt_tau);
  setptr("filt_cfreq", &filt_cfreq);
  setptr("dsomadia",   &dsomadia);
  setptr("btrend",     &btrend);
  setptr("ntrend",     &ntrend);
  setptr("nperiod",    &nperiod);
  setptr("nphase",     &nphase);

  setptr("speriod",    &speriod);
  setptr("sphase",     &sphase);
  setptr("orient",     &orient);
  setptr("tfreq",      &tfreq);
  setptr("drift",      &drift);

  setptr("reczmax",    &reczmax);
  setptr("reczmin",    &reczmin);

  setptr("set_timinc", &set_timinc);
  setptr("set_crit",   &set_crit);

  nvalfile = "nval_dsgc_sbac.n";
  if (notinit(arb_scale)) arb_scale = 1.0;
  chanparamsfile = "chanparams_dsgc_chans";

  if (notinit(dbp1_anpi)) dbp1_anpi = 0;
  if (notinit(dbp1_anpo)) dbp1_anpo = 0;
  if (notinit(dbp2_anpi)) dbp2_anpi = 0;
  if (notinit(dbp2_anpo)) dbp2_anpo = 0;

  //fprintf (stderr,"makestim %d\n",makestim);
  sb_init();
}

/*------------------------------------------------------*/

#include "gprim.h"

void syn_draw2 (synapse *spnt, int color, double vrev, double dscale, double dia,
                           double length, double foreshorten, int hide)
{
    int fill=1;

  dia *= dscale;                        /* draw circle with line */
  if (dia < 0) dia = -dia;
  color = -1;
  if (color < 0) {
     if (vrev < -0.04) color = RED;
     else              color = CYAN;
  }
  gpen (color);
 if (length > 1e-3) {
     gmove (length/2.0,0.0);
     if (dia > 0.001) gcirc (dia/2.0,fill);
     else             gcirc (0.001,fill);
     gmove (0,0);
     gdraw (length,0);
   }
  else                gcirc (0.001,fill);
}




/*--------------------------------------------------------*/

double calck1nx (double v, int func)

/* Calculate K rate functions given voltage in mv.
   All rates are calculated exactly as in HH (1952) paper.
   The "func" parameter defines:

    1	alpha n 
    2	beta  n
*/

#define MSSEC 1000.0
#define Vo (-65.)

{
   double val,x,y;

  switch (func) {
						/* alpha functions */

  case 1:					/* alpha n */
/*    y = 0.1 * (v+10.);                        /* the old way */
    y = -0.1 * (v+55.) * kexp;                  /* the modern way */
    x = exp (y) - 1.;
    if (abs(x) > 1e-5)                          /* singularity when x==0 */
       val = MSSEC * 0.1 * y / x;
    else
       val = MSSEC * 0.1;
    val *= dratehhna;				/* normalize to dbasetc */ 
    break;

  case 2:					/* beta n */
/*    val = 0.125 * exp (v / 80.);              /* the old way */

    val = MSSEC * 0.125 * exp ((v-Vo) / -80. * kexp);  /* the modern way */
    val *= dratehhnb;				/* normalize to dbasetc */ 
    break;

  }
  return val;
}						/*  i.e. 22 deg C */

/*--------------------------------------------------------*/

void setparams(void)

  /*  set up default configuration for sb expts */
  /* cones, cone bipolars, sb, dsgc */

{
   int i, ct;
   double zmax, zmin;

  make_rods  = 0;
  make_cones = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_hbp1  = 1;
  make_hbp2  = 0;
  make_ams   = 1;
  make_amhs  = 1;
  make_sbac  = 0;
  make_dsgc  = 1;

  make_dbp1_ams  = 0;
  make_hbp1_amhs = 0;

  ct=dsgc;

  setn(ct,NCOLOR,RCOLOR);
  set_synapse_dr (syn_draw2);
  if (!notinit(soma_z)) setn(ct,SOMAZ,soma_z);			// -44 for morph_R3RB140523_01

//  DENDD     = R_1;
//  DEND_DIST = R_1;
//  DEND      = R_2;	/* definitions of regions for dens_ file */
//  DENDP     = R_3;
//  DEND_PROX = R_3;
//  SOMA      = R_4;
//  HCK       = R_5;
//  HILLOCK   = R_5;
//  AXONT     = R_6;
//  AXON_THIN = R_6;
//  AXON      = R_7;
//  AXONP     = R_7;
//  AXON_PROX = R_7;
//  AXOND     = R_8;
//  AXON_DIST = R_8;
//  VARIC     = R_9;
//  VARICOS   = R_9;
  
#define SBARR 20

  if (make_movie) {
    onplot_dsgc_movie_init();		/* initialize dsgc movie stuff */
    onplot_movie_init();			/* initialize onplot_movie stuff */
  }

  pickden[dsgc] = 0; //941;       	/* pick one dendrite to allow connections */
  if (notinit(ath_dia)) ath_dia = 0.4;	/* diameter of initial thin axon segment */
  if (notinit(dsomadia)) dsomadia = 12;	/* diameter of dsgc soma */

  if (n_dsgc>0) make_ct(dsgc);		/* make dsgcs if user specifies */

  if (!notinit(no_inhib)) { if (no_inhib==1) setsv(sbac,SCOND,3,0); } 
  if (notinit(revdir)) revdir = 0;	/* reverse pref and null directions */ 

  if (notinit(capfilt)) capfilt = 0;	/* set high res printout and 80 usec filter on isoma */
  if (notinit(btrend)) btrend = 0;	/* linear baseline trend */
  if (notinit(ntrend))  ntrend = 0;	/* sine baseline trend */
  if (notinit(nperiod)) nperiod = 16.6666e-3; /* sine phase baseline trend */
  if (notinit(nphase))  nphase = 0;	/* sine phase baseline trend */

  if (notinit(ivplot)) ivplot = 0;	/* make I/V plot */
  if (notinit(dvrev)) dvrev = -0.06;
  if (notinit(dvst))  dvst  = -0.06;

  if (!notinit(ampa_cond)) setsv(dbp1,SCOND,2,ampa_cond);
  if (!notinit(ampa_cond)) setsv(hbp1,SCOND,4,ampa_cond);
  if (!notinit(nmda_cond)) setsv(dbp1,SCOND,7,nmda_cond);
  if (!notinit(nmda_cond)) setsv(hbp1,SCOND,5,nmda_cond);
  if (!notinit(gaba_cond)) setsv(ams, SCOND,1,gaba_cond);
  if (!notinit(gaba_cond)) setsv(amhs,SCOND,1,gaba_cond);
  if (notinit(dendrm)) dendrm = drm;
  if (notinit(dendcm)) dendcm = dcm;
  if (notinit(somacm)) somacm = dcm;
  if (notinit(tipcap)) tipcap = 0;
  if (notinit(dendrs)) dendrs = 1.0/dri;	// 65000 for dri=200 to predict currents;
  if (notinit(cplam))   cplam = complam;	// size of compartments for density file 

  if (stimtype==5) {				// set linear gain for sine wave stimulus
	   setsv(hbp1,SGAIN,4,5);
	   setsv(hbp1,SVGAIN,4,0);		// set linear gain
  }


  if (notinit(dsgc_prefdir)) dsgc_prefdir=0;

  if (strcmp(sbac_file,"morph_sb1")==0) {
    setsv (sbac,SYNANNI,3,75);                    /* don't connect close to soma */
  }
  else if (strcmp(sbac_file,"morph_sbac3b")==0) {
    setsv (sbac,SYNANNI,3,20);                    /* don't connect close to soma */
  }

  if (!notinit(set_drm)) {                      /* user set default Rm */
       setn(ct,NRM,set_drm);                    /* set default Rm */
       drm = set_drm;
  }
  if (notinit(set_kchan)) set_kchan = 1;
  if (set_kchan>=0) {
       set_chancalc(K,1,0,calck1nx); 		/* sets the 0 (n) param for K type 1 */
       if (notinit(kexp)) kexp = 1;		/* exponent for K activation in calck1nx() above */
  }
  dsintinc = 1e-4;

   //display_z(zmax=-5, zmin=-15);	    /* display sublamina a */
   //display_z(zmax=-15, zmin=-25);	    /* display sublamina b */
}

/*--------------------------------------------------------*/

void setdens (void)
/* Run by retsim after channel densities and properties have been read in. */

{
     int ct;

  if (!notinit(kdr_voff)) set_chan_offset(ct=dsgc, C_K1 , CHOFFM, kdr_voff); /* set KDR voltage offset */
}

/*--------------------------------------------------------*/

void addlabels(void) 

{
     int dsgc_rn;

  if (!notinit(ams1)) { dsgc_rn = findsynlocp(ams,ams1,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams1");  }
  if (!notinit(ams2)) { dsgc_rn = findsynlocp(ams,ams2,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams2");  }
  if (!notinit(ams3)) { dsgc_rn = findsynlocp(ams,ams3,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams3");  }
  if (!notinit(ams4)) { dsgc_rn = findsynlocp(ams,ams4,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams4");  }
  if (!notinit(ams5)) { dsgc_rn = findsynlocp(ams,ams5,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams5");  }
  if (!notinit(ams6)) { dsgc_rn = findsynlocp(ams,ams6,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams6");  }
  if (!notinit(ams7)) { dsgc_rn = findsynlocp(ams,ams7,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams7");  }
  if (!notinit(ams8)) { dsgc_rn = findsynlocp(ams,ams8,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams8");  }
  if (!notinit(ams9)) { dsgc_rn = findsynlocp(ams,ams9,dsgc,1);  label(ndn(dsgc,1,dsgc_rn),blue,"ams9");  }
  if (!notinit(ams10)) { dsgc_rn = findsynlocp(ams,ams10,dsgc,1);label(ndn(dsgc,1,dsgc_rn),blue,"ams10");  }
}

/*--------------------------------------------------------*/

void runonexit (void)
{
       if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/


bool flag = false;
int cellpair=0;
double  ct = dsgc;

double maxCurrent;
double cond,Gmax;
double voltage;
double current, current2;
double idiff;

void onplot(void) {
    current  = i(ndn(ct,  1, elnode));
    if (cellpair) current2 = i(ndn(ct, 2, elnode));
    else          current2 = 0;
    idiff = current - current2;

    voltage  = v(ndn(ct,  1, soma));

  if (flag) {
   if (outward<=0) {
      cond = idiff / (voltage - gvrev);
      if (cond > Gmax) Gmax = cond;
      if (maxCurrent > idiff) {         // inward current
        maxCurrent = idiff;
      //  fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   } else {                             // outward current 
      cond = idiff / (voltage - vk);
      if (cond > Gmax) Gmax = cond;
      if (maxCurrent < idiff) {
        maxCurrent = idiff;
        //fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   }
   //  fprintf(stderr, "i: %g, maxi: %g\n", current, maxCurrent);    
  }
  else cond = 0;  
}

/*------------------------------------------------------*/

double mvbrt1;
double mvbrt2;
double ixoff, iyoff;
double ncycles = 0;

double move_stim(double stimtime,double celldia, double roffset, double barwidth,double theta,double velocity, double ncycles, double econtrast, double icontrast, double eincr, double iincr, int direction, double mask)
{
    static int runyet=0;
    double lbar, rbar, irad, orad;
    double inten, start, dur, wavel;
    double s,send;
    double sdia;

 // stim_spot(100, 100,0, inten=minten,  start=0.02, dur=0.05, wavel=1, 0);

 if (stimtype==1) {			// move bar to & fro
   lbar = -(celldia + barwidth) * 0.6 + roffset;
   rbar =  (celldia + barwidth) * 0.6 + roffset;

   if (revdir) {
	   double temp;
	   temp = lbar;
	   lbar = rbar;
	   rbar = temp;
   }

   // pref: move bar from left to right
 
   mvbrt1 = movebar (stimtime, 0,0,         lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);

   if (light_inhib) {

	 mvbrt1 = movebar (stimtime, ixoff,iyoff, lbar+ioffset, rbar+ioffset, 
		   				barwidth,barlength,theta,velocity,icontrast);

	/* rev_inhib, make bar with different offset, null */

	if (rev_inhib) movebar (stimtime, -ixoff,-iyoff, lbar+ioffset, rbar+ioffset, 
		   				barwidth,barlength,theta,velocity,icontrast+iincr);
   }

   // null: move bar from right to left
  
   mvbrt2 = movebar (mvbrt1+stimtime,0,0, rbar, lbar, barwidth,barlength,theta,velocity,econtrast);

   // make null inhibition stronger
   
   if (light_inhib) {
	mvbrt2 = movebar (mvbrt1+stimtime,ixoff,iyoff, rbar+ioffset, lbar+ioffset, 
						barwidth,barlength,theta,velocity,icontrast+iincr);
	/* rev_inhib, make bar with different offset, pref */

	if (rev_inhib) movebar (mvbrt1+stimtime,-ixoff,-iyoff, rbar+ioffset, lbar+ioffset, 
		   				barwidth,barlength,theta,velocity,icontrast);
   }
 }
 else if (stimtype==2) {		// move annulus to & fro
	 double masktheta=0;

   //lbar =  (celldia/4 + barwidth) * 0.6;
   lbar =  -(celldia + barwidth) * 0.5;
   rbar =  (celldia   + barwidth) * 0.5;

   // mask bar 

   stim_bar(barwidth*1.5,barlength, 0, 0, masktheta, -0.06, start=stimtime, dur=1.0, mask);
   if (light_inhib) stim_bar(barwidth*1.5,barlength, ixoff, iyoff, masktheta, -0.06, start=stimtime, dur=1.0, mask);

   // move bar from left to right

   mvbrt1 = movebar (stimtime, 0,0,         lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);
   if (light_inhib) mvbrt1 = movebar (stimtime, ixoff,iyoff, lbar+ioffset, rbar+ioffset, 
		   					barwidth,barlength,theta,velocity,icontrast);
   // move bar from right to left
  
   mvbrt2 = movebar (mvbrt1+stimtime,0,0, rbar, lbar, barwidth,barlength,theta,velocity,econtrast);

   // make inhibition stronger from right to left
   
   if (light_inhib) mvbrt2 = movebar (mvbrt1+stimtime,ixoff,iyoff, rbar+ioffset, lbar+ioffset, 
		   					   barwidth,barlength,theta,velocity,icontrast+iincr);
 }
 else if (stimtype==3) {		// static annulus, radius = stimx
     //       stim_spot(100, 0, 0, econtrast, stimtime, spotdur);
        mvbrt1 = moveannulus(stimtime,        0,     0,     stimx+50,  stimx+55,  barwidth, velocity, econtrast);
        mvbrt2 = moveannulus(stimtime+0.1,    0,     0,     stimx+120, stimx+125, barwidth, velocity, econtrast);
    if (light_inhib) {
        mvbrt1 = moveannulus(stimtime,       ixoff,  iyoff, stimx+50,  stimx+55,  barwidth, velocity, icontrast);
        mvbrt2 = moveannulus(stimtime+0.1,   ixoff,  iyoff, stimx+120, stimx+125, barwidth, velocity, icontrast);
    }
    mvbrt2 = 0.2;

 }
 else if (stimtype==4) {		// move spot
	  sdia = barwidth;
 	  send = int(2000/sdia + 0.5);
	  for (s=0; s<send; s++) {
            stim_spot(sdia, s*sdia, 0, econtrast, stimtime+s*spotdur, spotdur);
          }
	  mvbrt2 = spotdur * send;
       }
  else if (stimtype==5) {		// sine waves 2,5,10,20,50 hz 
          stim_sine (speriod, sphase, orient, 0, 0, tfreq, drift, 1, scontrast, sqwave, start=stimtime,  stimdur);
	  mvbrt2 = stimdur+stimtime;
 }
 return mvbrt2;
}

/*--------------------------------------------------------*/

double vdiff (double nod1, double nod2, double time)

{
    int ct=dsgc,cn=1;
    double v1, v2;

  v1 = v(ct,cn,int(nod1));
  v2 = v(ct,cn,int(nod2));
  return (v2 - v1)/(dri*dendrs);
}

/*--------------------------------------------------------*/

void runexpt(void)

{
#define NREVINHIB 20
    int c, i, ct, cn, colr, pl, plnum, electrode_node, dist_node;
    int rev_cn, revnod;
    int nsynap, nsynapa, nsynapn, nsynapi;
    int n_rev_inhibs, rev_inhibs[NREVINHIB];
    double start, starttime, dur, dscale, plsize;
    double disp_end, mask, psize, pulsedur;
    double rmax, rmin;
    double Imax, Imin, imax, imin, gmax;
    double vpulse;
    double mvbrt1, sign;
    double celldia, barloc;
    double cmin, cmax;
    node *npnt;
    elem *epnt;
    photorec *p;
    chattrib *a;
    char sbuf[30];

  //ploti = 1e-3;
  ploti = 1e-4;
  timinc = 1e-5;
  if (!notinit(set_timinc)) timinc = set_timinc;
  if (!notinit(set_crit))     crit = set_crit;

  if (capfilt==1) ploti = 20e-6;

  ct = dsgc;
  electrode_node = 5000;

  Vmax  = -0.02;
  Vmaxg = 0.00;
  Vmin = -0.07;
  
   /* add light transducer to each bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }
   for(npnt=nodepnt; npnt=foreach(npnt,hbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(hbp1,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }

   if (notinit(theta))   theta = 0;	/* motion of bar */
   if (notinit(iroff))   iroff = 1000;	/* r offset for inhib transducers */

   ixoff = iroff * cosdeg(theta);
   iyoff = iroff * -sindeg(theta);

   if (notinit(light_inhib)) light_inhib = 0; 
   if (light_inhib) {

     /* Add light transducer to each small-field amacrine cell */
     /*  offset so that inhibition can be controlled separately */

     for(npnt=nodepnt; npnt=foreach(npnt,ams,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       if ((p=(photorec*)make_transducer(ndn(ams,cn,soma))) != NULL) { 
            p->xpos=npnt->xloc + ixoff; 
            p->ypos=npnt->yloc + iyoff;
       }
     }

     /* Test of inhibition with reverse DS at one node */
     /* Make transducer with opposite offset so it can have greater pref inhibition */

     if (notinit(rev_inhib)) rev_inhib = 0;
     if (rev_inhib) {			/* test of inhibition with reverse DS at one node */

	     /* get amacrine cells that will be reversed and their postsynaptic nodes on dsgc */

	 // if ((rev_cn=findcell(ams,110,30)) > 0) {
         i = 0;
	 if (!notinit(ams1)) { rev_inhibs[i++] = ams1; }
	 if (!notinit(ams2)) { rev_inhibs[i++] = ams2; }
	 if (!notinit(ams3)) { rev_inhibs[i++] = ams3; }
	 if (!notinit(ams4)) { rev_inhibs[i++] = ams4; }
	 if (!notinit(ams5)) { rev_inhibs[i++] = ams5; }
	 if (!notinit(ams6)) { rev_inhibs[i++] = ams6; }
	 if (!notinit(ams7)) { rev_inhibs[i++] = ams7; }
	 if (!notinit(ams8)) { rev_inhibs[i++] = ams8; }
	 if (!notinit(ams9)) { rev_inhibs[i++] = ams9; }
	 if (!notinit(ams10)) { rev_inhibs[i++] = ams10; }

	 if ((n_rev_inhibs = i) > NREVINHIB) { fprintf(stderr,"n_rev_inhibs: too many entries %d\n",n_rev_inhibs);
		 			n_rev_inhibs = NREVINHIB-1; }

	 for (i=0; i<n_rev_inhibs; i++) {
	   rev_cn = rev_inhibs[i];
           for(epnt=elempnt; epnt=foreach(epnt,VTRANSDUCER,ams,rev_cn,soma,NULL,NULL,NULL); epnt=epnt->next) {
	     p = (photorec*)epnt;
	     npnt = ndn(ams,rev_cn,soma);
             p->xpos=npnt->xloc - ixoff; 
             p->ypos=npnt->yloc - iyoff;
	     break;
	   }
	}	
     }

     for(npnt=nodepnt; npnt=foreach(npnt,amhs,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       if ((p=(photorec*)make_transducer(ndn(amhs,cn,soma))) != NULL) { 
            p->xpos=npnt->xloc + ixoff; 
            p->ypos=npnt->yloc + iyoff;
       }
     }
   }

  /*  - - - - - - - - - - - - - - - - - - - */

   /* make lists of synaptic inputs to the dsgc */

    nsynapi  = synapse_add (1,ams, -1,-1,dsgc,1,1);       /* make list of inhib synapses from ams's */
    nsynapi += synapse_add (1,amhs,-1,-1,dsgc,1,1);       /* make list of inhib synapses from amhs's */
    nsynapa  = synapse_add (2,dbp1,-1,-1,dsgc,1,2);       /* make list of dbp1 ampa synapses for recording below */
    nsynapa += synapse_add (2,hbp1,-1,-1,dsgc,1,4);       /* make list of hbp1 ampa synapses for recording below */
    nsynapn  = synapse_add (7,dbp1,-1,-1,dsgc,1,7);       /* make list of dbp1 nmda synapses */
    nsynapn += synapse_add (7,hbp1,-1,-1,dsgc,1,5);       /* make list of hbp1 nmda synapses */
    nsynap = nsynapa + nsynapn + nsynapi;
    fprintf (stderr,"# nsynap %d nsynapa %d nsynapn %d nsynapi %d\n",nsynap, nsynapa, nsynapn, nsynapi);

  /*  - - - - - - - - - - - - - - - - - - - */

   if (notinit(barwidth))        barwidth = 100;
   if (notinit(barlength))      barlength = 500;
   if (notinit(minten))            minten = -0.05;
   if (notinit(econtrast))      econtrast =  0.007;
   if (notinit(scontrast))      scontrast =  0.007;
   if (notinit(icontrast))      icontrast =  scontrast;
   if (notinit(eincr))              eincr =  0;
   if (notinit(sincr))              sincr =  0.002;
   if (notinit(iincr))              iincr =  sincr;
   if (notinit(velocity))        velocity =  2000; 
   if (notinit(stimscale))      stimscale =  1; 
   if (notinit(stimdur))          stimdur =  0.5;	   /* used for non-moving stimuli */
   if (notinit(prestimdur))    prestimdur =  0.02;
   if (notinit(poststimdur))  poststimdur =  0.05;
   if (notinit(tailcurdur))    tailcurdur =  0.02;
   if (notinit(stimtime))        stimtime =  0.1;
   if (notinit(disptime))        disptime =  0.15;
   if (notinit(stimtype))        stimtype =  1;
   if (notinit(spotdur))          spotdur =  0.05;

   if (notinit(vhold))              vhold = -0.06;
   if (notinit(vstart))            vstart = -0.120;
   if (notinit(vstop))              vstop =  0.04;
   if (notinit(vstep))              vstep =  0.02;
   if (notinit(tailvolt))        tailvolt =  vhold;
   if (notinit(gvrev))              gvrev =  vna;

   if (notinit(speriod))          speriod = 1000;
   if (notinit(sphase))            sphase = 90;
   if (notinit(orient))            orient = 0;
   if (notinit(tfreq))              tfreq = 10;
   if (notinit(drift))              drift = 0;
   if (notinit(sqwave))            sqwave = 0;
   if (notinit(stimx))              stimx = 0;
   if (notinit(stimy))              stimy = 0;
   if (notinit(stimr))              stimr = 0;

   if (notinit(rectheta))        rectheta = theta;	/* orientation of stimulus recordings */
   if (notinit(reczmax))          reczmax = 0;		/* z pos max (um) for stimulus recordings */
   if (notinit(reczmin))          reczmin = -100;	/* z pos min (um) for stimulus recordings */

   set_run_on_exit(runonexit);                         // set to erase savefile on ^C
   sprintf (savefile,"dsgc_chans%06d",getpid());       // add pid to file name

   if (notinit(ioffset)) ioffset = -barwidth;

   if (notinit(elec_rs))  elec_rs   = 20e6;
   if (notinit(elec_cap)) elec_cap  = 1e-14;
   if (notinit(shunt_res)) shunt_res  = 1e30;
   if (notinit(shunt_v))   shunt_v  = dvrev;
   if (notinit(elnode))     elnode  = electrode_node;

   if (elnode==electrode_node) {
       make_electrode  (nd(ct,cn=1,elnode), ndn(ct,cn=1,soma), elec_rs, elec_cap, tipcap);
       make_shunt  (ndn(ct,cn=1,soma), shunt_res, shunt_v);
   }
       
   cn = 1;
   if (make_movie) {
     if (space_time) {  /* movie */ 
      plot_v_nod(ct=dsgc,cn,soma,Vmin,Vmaxg,c=1,"Vsoma",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn,1336,Vmin,Vmaxg,c=green,"Vtip1",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn,582,Vmin,Vmaxg,c=red,"Vtip2",pl=10,0.35);
      //plot_na_inact(dsgc, 1, soma, red, "Na[soma]", 12, .35);
     };
   }
   else { 
    Vmin = -0.049; Vmax = -0.045;
     if (stimtype==5)
        plot_v_nod(hbp1,findmida(hbp1,50,-150),soma,Vmin, Vmax, 1, "", 22, 0.3);
     // plot_v_nod(hbp1,findmida(hbp1,50,-100),soma,Vmin, Vmax, 2, "", 22, 0.3);
     // plot_v_nod(hbp1,findmida(hbp1,50,-50), soma,Vmin, Vmax, 3, "", 22, 0.3);
     // plot_v_nod(hbp1,findmida(hbp1,50,0),   soma,Vmin, Vmax, 4, "", 22, 0.3);
     // plot_v_nod(hbp1,findmida(hbp1,50,50),  soma,Vmin, Vmax, 5, "", 22, 0.3);
     // plot_v_nod(hbp1,findmida(hbp1,50,100), soma,Vmin, Vmax, 6, "", 22, 0.3);

    if (stimtype==2)
         barloc = 100;
    else barloc = 50;

    plot_syncond(findsynlocra(dbp1,-150,rectheta),  cmin=0,cmax=100e-12, 1, 20,"",0.3);
    plot_syncond(findsynlocra(dbp1,-100,rectheta),  cmin=0,cmax=100e-12, 2, 20,"",0.3);
    plot_syncond(findsynlocra(dbp1,-50, rectheta),  cmin=0,cmax=100e-12, 3, 20,"",0.3);
    plot_syncond(findsynlocra(dbp1,0,   rectheta),  cmin=0,cmax=100e-12, 4, 20,"",0.3);
    plot_syncond(findsynlocra(dbp1,50,  rectheta),  cmin=0,cmax=100e-12, 5, 20,"",0.3);
    plot_syncond(findsynlocra(dbp1,100, rectheta),  cmin=0,cmax=100e-12, 6, 20,"",0.3);

    plot_syncond(findsynlocra(ams,-150, rectheta),  cmin=0,cmax=100e-12, 1, 18,"",0.3);
    plot_syncond(findsynlocra(ams,-100, rectheta),  cmin=0,cmax=100e-12, 2, 18,"",0.3);
    plot_syncond(findsynlocra(ams,-50,  rectheta),  cmin=0,cmax=100e-12, 3, 18,"",0.3);
    plot_syncond(findsynlocra(ams,0,    rectheta),  cmin=0,cmax=100e-12, 4, 18,"",0.3);
    plot_syncond(findsynlocra(ams,50,   rectheta),  cmin=0,cmax=100e-12, 9, 18,"",0.3);
    plot_syncond(findsynlocra(ams,100,  rectheta),  cmin=0,cmax=100e-12, 10, 18,"",0.3);

    if (rev_inhib) {		/* plot node postsyn to rev_inhib ams */
        for (i=0; i<n_rev_inhibs; i++) {
            plot_syncond(findsynlocr(ams,rev_inhibs[i],0,0),cmin=0,cmax=100e-12, 5+i, 18,"",0.3);
	}
    }
    plot_syncond(findsynlocra(hbp1,-150,rectheta),  cmin=0,cmax=100e-12, 1, 20,"",0.3);
    plot_syncond(findsynlocra(hbp1,-100,rectheta),  cmin=0,cmax=100e-12, 2, 20,"",0.3);
    plot_syncond(findsynlocra(hbp1,-50, rectheta),  cmin=0,cmax=100e-12, 3, 20,"",0.3);
    plot_syncond(findsynlocra(hbp1,0,   rectheta),  cmin=0,cmax=100e-12, 4, 20,"",0.3);
    plot_syncond(findsynlocra(hbp1,50,  rectheta),  cmin=0,cmax=100e-12, 5, 20,"",0.3);
    plot_syncond(findsynlocra(hbp1,100, rectheta),  cmin=0,cmax=100e-12, 6, 20,"",0.3);

    plot_syncond(findsynlocra(amhs,-150,rectheta),  cmin=0,cmax=100e-12, 1, 18,"",0.3);
    plot_syncond(findsynlocra(amhs,-100,rectheta),  cmin=0,cmax=100e-12, 2, 18,"",0.3);
    plot_syncond(findsynlocra(amhs,-50, rectheta),  cmin=0,cmax=100e-12, 3, 18,"",0.3);
    plot_syncond(findsynlocra(amhs,0,   rectheta),  cmin=0,cmax=100e-12, 4, 18,"",0.3);
    plot_syncond(findsynlocra(amhs,50,  rectheta),  cmin=0,cmax=100e-12, 5, 18,"",0.3);
    plot_syncond(findsynlocra(amhs,100, rectheta),  cmin=0,cmax=100e-12, 6, 18,"",0.3);


    /* if (notinit(node_dist)) node_dist = 281;

     plot_chan_current(ct=dsgc, cn, node_dist, NMDA, 1, 3e-12, -3e-12);
     sprintf (sbuf,"Inmda  %d",node_dist);
     plot_param (sbuf, blue, plnum=5, plsize=0.5);

     plot_chan_current(ct=dsgc, cn, node_dist, AMPA, 5, 3e-12, -3e-12);
     sprintf (sbuf,"Iampa  %d",node_dist);
     plot_param (sbuf, gray, plnum=5, plsize=0.5);
    */

    if (notinit(set_vclamp)) set_vclamp = 1;

    Vmin = dvrev; Vmax = max(vstop,-0.04);
    if (set_vclamp > 0) {
      //plot_v_nod(ct, cn, findnodlocraz(ct,cn,-40,rectheta,reczmax,reczmin),  Vmin, Vmax, colr=blue,  "", 5, 1);
      plot_v_nod(ct, cn, elnode,                   Vmin, Vmax, colr=white, "Vdsgc_1_elec", 5, 1);
      plot_v_nod(ct, cn, soma,                     Vmin, Vmax, colr=green, "", 5, 1);
      plot_v_nod(ct, cn, findnodlocraz(ct,cn,40,rectheta,reczmax,reczmin),   Vmin, Vmax, colr=cyan,  "", 5, 1);
      plot_v_nod(ct, cn, findnodlocraz(ct,cn,80,rectheta,reczmax,reczmin),   Vmin, Vmax, colr=brown, "", 5, 1);
      plot_v_nod(ct, cn, findnodlocraz(ct,cn,120,rectheta,reczmax,reczmin),  Vmin, Vmax, colr=red,   "", 5, 1);
      plot_i_nod(ct=dsgc,cn=1,elnode,       Imin=-500e-12,Imax=500e-12,green, "Idsgc_1_elec", 3, 1);

      //if (capfilt==1) plotfilt(1,make_filt(79.6e-6));
      if (capfilt==1) { 	// add temporal filter to last plot (Idsgc_elec), takes ploti into account
	  if (notinit (filt_cfreq)) filt_cfreq = 2e3;
          plotbessfilt(filt_cfreq);		// 4th order bessel filter, cutoff = 2000 Hz
         //plotfilt(4,make_filt(filt_tau,filt_tau,filt_tau,filt_tau));
      }
    } else if (rev_inhib) {		/* plot node postsyn to rev_inhib ams */

         plot_v_nod(ct, cn, soma,                                       Vmin, Vmax, colr=red, "", 5, 1);
	 if (!notinit(dsgc_n1)) plot_v_nod(ct, cn, dsgc_n1,             Vmin, Vmax, colr=green, "", 5, 1);
	 if (!notinit(dsgc_n2)) plot_v_nod(ct, cn, dsgc_n2,             Vmin, Vmax, colr=cyan, "", 5, 1);
	 if (!notinit(dsgc_n3)) plot_v_nod(ct, cn, dsgc_n3,             Vmin, Vmax, colr=ltred, "", 5, 1);
	 for (i=0; i<n_rev_inhibs; i++) {
            revnod = findsynlocp(ams,rev_inhibs[i],ct,cn);
            sprintf (sbuf,"ams%d->Vdsgc%d",rev_inhibs[i],revnod);
            plot_v_nod(ct, cn, revnod,                                  Vmin, Vmax, colr=magenta+i, sbuf, 5, 1);
	 }

      } else { 			/* normal display for no set_vclamp */

      plot_v_nod(ct=dsgc,cn=1,soma,        Vmin,Vmax,red, "", 5, 0.35);
      plot_v_nod(ct=dsgc,cn,findnodlocraz(ct,cn,40,rectheta,reczmax,reczmin),    Vmin,Vmax,cyan,"", 5, 0.35);
      plot_v_nod(ct=dsgc,cn,findnodlocraz(ct,cn,80,rectheta,reczmax,reczmin),    Vmin,Vmax,brown,"", 5, 0.35);
      plot_v_nod(ct=dsgc,cn,findnodlocraz(ct,cn,120,rectheta,reczmax,reczmin),   Vmin,Vmax,blue,"", 5, 0.35);
      // plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.3); 
      //

    } 

      imax =  0e-10;
      imin = -3e-10;
      gmax =  5e-9;

      plot_func(isyn_tot,2,imax,imin);       plot_param("Itotbpa",blue,2,0.3);
      plot_func(gsyn_tot,2,gmax,0);          plot_param("Gtotbpa",blue, 1,0.3);
      plot_func(gsyn_tot,7,gmax,0);          plot_param("Gtotbpn",cyan, 1,0.3);
      plot_func(gsyn_tot,1,gmax,0);          plot_param("Gtotami",red,  1,0.3);

      // dist_node = findnodloc(ct=dsgc,cn=1,120,0);
      // sprintf (sbuf,"Vdiff_%d_soma",dist_node);
      // plot_func(vdiff,soma,dist_node, 0.1, -0.05); plot_param (sbuf,red, 4, 0.3);
   }

   celldia = max(xarrsiz,yarrsiz) * 1.5 * stimscale; 
   // fprintf (stderr,"# celldia %g\n",celldia);

   if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
   if (notinit(predur))        predur = 0.05;

   if (!make_movie) {
     if (disp==16) {
	double t;

      simtime = 0;                                    // must be set ahead of stim_backgr()
      stim_backgr(minten,start=simtime);			 /* turn on  background */

      stimdur = move_stim(stimtime, celldia, stimr, barwidth, theta, velocity, ncycles, 
		      		econtrast, icontrast, eincr, iincr, direction, mask=1);
      if (light_inhib) display_size(2500); else display_size(800); 
      disp_end = stimtime+prestimdur+stimdur+tailcurdur+poststimdur;
      for (starttime=0,t=stimtime; t<disp_end; starttime = t, t+= 0.002) {
           if (rev_inhib) display_stim(starttime, t, dscale=4, 2, -0.040, -0.045);
	   else display_stim(starttime, t, dscale=4, -0.045, -0.049);
           //display_stim(t, dscale=4, -0.035, -0.045); 
 	   //display_stim(0+t, dscale=4); 
           simwait(0.10);
      }
      return;
    }
  }

   simtime = -predur;                                    // must be set ahead of stim_backgr()
   stim_backgr(minten,start=simtime);			 /* turn on  background */

  /* set movie plot routine */

  if (make_movie) setonplot(onplot_movie);

  /* run experiment */

  if (notinit(setendexp)) 
    stimdur = 2 * (celldia*1.6  + barwidth) / velocity;
  if (stimtype==3) {
    stimdur = 0.1;
  }	
  if (!notinit(setendexp)) endexp = setendexp;
  else endexp=stimtime+prestimdur+stimdur+tailcurdur+poststimdur;

  if (notinit(istim)) istim = 0; 
  if (istim != 0) {
      cclamp(ndn(ct,cn,elnode), istim, start=simtime, dur=1);
  }
  if (abs(btrend)>0) {
	  double t, tinc;
      for (tinc=1e-4,t=0; t<endexp; t+=tinc) {
          cclamp(ndn(ct,cn,elnode), btrend*t, start=t, dur=tinc);
     }
  }
  if (abs(ntrend)>0) {
	  double t, tinc;
      for (tinc=20e-6,t=0; t<endexp; t+=tinc) {
          cclamp(ndn(ct,cn,elnode), ntrend*sin((t/nperiod+nphase)*2*PI), start=t, dur=tinc);
      }
  }

  if (set_vclamp > 0) {
    vclamp  (ndn(ct,cn, elnode), vhold, simtime,  predur);
  }

  step(predur);

  if (set_vclamp > 0) {
    dst = 0.0001;         // small time to allow voltage clamp to end

    if (ivplot) graph_pen(i+1,i+1,i+1,i+1,i+1);
    if (!ivplot) flag = true;
    if (vstart <= vstop) sign = 1;
    else                 sign = -1;

    simtime = 0;
    savemodel (savefile);

    for (i=0,vpulse=vstart; (vpulse*sign)<=(vstop*sign+1e-6); i++,vpulse += vstep) {

       if (!notinit(skipv)) {
	  switch (skipv) {
            case 0:  if (i!=0) continue; break;		              // do only -80 mV
            case 1:  if (i!=2) continue; break;		              // do only -80 mV
            case 2:  if ((i!=2) && (i!=6)) continue; break;          // do only -80, 0 mV
            case 3:  if ((i!=2) && (i!=4) && (i!=6)) continue; break; // do only -80, -40, 0 mV
            case 4:  if ((i!=0) && (i!=2) && (i!=4) && (i!=6)) continue; break; // do -120, -80, -40, 0 mV
            case 5:  if ((i!=0) && (i!=2) && (i!=4) && (i!=6) && (i!=8)) continue; break; // do -120, -80, -40, 0, 20 mV
            case 12: if ((i!=0) && (i!=6)) continue; break;           // do only -120, 0 mV
            case 14: if ((i!=0) && (i!=2)) continue; break;           // do only -120, -80 mV
            default: break;
          }
       }

       simtime = 0;
       if (!notinit(setendexp)) pulsedur = stimdur;
       else {
            pulsedur = move_stim(stimtime, celldia, stimr, barwidth, theta, velocity, ncycles, econtrast, icontrast, eincr, iincr, direction, mask=1);
       }
       // fprintf (stderr,"prestimdur %g stimdur %g tailcurdur %g poststimdur %g endexp %g\n",prestimdur,stimdur,tailcurdur,poststimdur,endexp);
       // pulsedur += prestimdur;

       vclamp              (ndn(ct,cn, elnode), vhold, simtime,  prestimdur);
       step (prestimdur);

       if (outward) maxCurrent =  0;               // for inward current
       else         maxCurrent =  1000;            // for outward current
       Gmax = 0;
       vclamp              (ndn(ct,cn, elnode), vpulse, simtime,  pulsedur);

       step (dst);
       if (ivplot) flag = true;

       if (capfilt > 0) {
           step (pulsedur);
	   continue;
       }

       if (stimtype==5) {	// if sine wave stimulus
          step (pulsedur);
       }
       else {

       step (pulsedur/2);

       restoremodel (savefile);

       step (pulsedur/2-2*dst);
       if (ivplot) graph(voltage, maxCurrent, Gmax);
       if (ivplot) flag = false;
       step (dst);

       vclamp              (ndn(ct,cn, elnode), tailvolt, simtime,  tailcurdur);
       step (tailcurdur);

       vclamp              (ndn(ct,cn, elnode), vhold, simtime,  poststimdur);
       step (poststimdur);

       restoremodel (savefile);
      }
    }
  }  /* if (set_vclamp) */ 
  else {
	// fprintf(stderr,"endexp %g %g %g\n", endexp, mvbrt1, mvbrt2);

       if (!notinit(setendexp)) pulsedur = stimdur;
       else {
            pulsedur = move_stim(stimtime, celldia, stimr, barwidth, theta, velocity, ncycles, econtrast, icontrast, eincr, iincr, direction,mask=1);
       }
	if (stimtype==1 || stimtype==2) {
	   savemodel (savefile);
	}
	step(mvbrt1+stimtime);

	if (stimtype==1 || stimtype==2) {
	  restoremodel (savefile);
	  step(mvbrt2-mvbrt1+stimtime);
	}

  } /* no vclamp */

  unlink (savefile);
  savefile[0] = 0;
  fflush(stdout);
  printf ("# done\n");		// for the "modelfit" program to see when the simulation is done
}
