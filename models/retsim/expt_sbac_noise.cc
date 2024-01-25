/* Experiment sbac_noise for retsim */
/*   derived from sbac_chans */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "ncio.h"

#include "stimfuncs.h"

double theta;
double iroff;
int light_inhib;
int node_dist;

int no_inhib;
int sbarr;
int rec_ct;
int rec_cn;
int exptrun;
int elnode;
int direction=0;
int revdir;
int sb_cn;

double ampa_cond;
double nmda_cond;
double gaba_cond;

double namid;
double nadist;
double kdist;
double kprox;
double kmid;
double ksoma;
double cadist;
double camid;
double catdist;
double catmid;
double catprox;
double catsoma;

double sbmaxpool;
double sbmaxrate;
double sbsynang;

double sb_vr;
double sb_vs;
double sb_rm;			// soma Rm
double sb_rmp;			// proximal Rm
double sb_rmd;			// distal Rm
double sb_rid;                  // distal Ri
double sb_rii;                  // intermediate Ri
double sb_rip;                  // proximal Ri
double sb_cmp;                  // proximal Cm
double sb_cmd;                  // distal Cm
double sbaclm;
double sdia;
double spdia;
double sddia;
double spdia_ratio;
double sndia;
double smdia;
double orad1;
double irad2;
double orad2;
double irad3;
double orad3;
double irad4;
double orad4;
double irad5;
double orad5;
double irad6;
double orad6;

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
double noise_dur;
double nbase;
double nscale;
double vscale;
double ppdelay;

double stim_theta;
double rstim_theta;

double dbp1_anpo;
double dbp1_anpi;
double dbp2_anpo;
double dbp2_anpi;

double dbp1_sgain;
double dbp1_vgain;
  
double sbac_isynanni;
double sbac_isynanpi;
double sbac_isynanpo;
double sbac_isynrngi;
double sbac_dens;

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
double arb_scale;
double kdr_voff;
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

const char *datfile;

// char savefile[30] = {0};

void sb_init(void);

/*--------------------------------------------------------*/

void defparams(void)

{
  setptr("theta", 	 &theta);
  setptr("iroff", 	 &iroff);
  setptr("light_inhib",  &light_inhib);
  setptr("node_dist",    &node_dist);
  setptr("no_inhib",	 &no_inhib);
  setptr("ampa_cond",	 &ampa_cond);
  setptr("nmda_cond",	 &nmda_cond);
  setptr("gaba_cond",	 &gaba_cond);
  setptr("rec_ct",	 &rec_ct);
  setptr("rec_cn",	 &rec_cn);
  setptr("sbarr",	 &sbarr);
  setptr("exptrun",	 &exptrun);
  setptr("elnode",	 &elnode);
  setptr("revdir",	 &revdir);
  setptr("sb_cn",	 &sb_cn);

  setptr("namid",	 &namid);
  setptr("nadist",	 &nadist);
  setptr("cadist",	 &cadist);
  setptr("camid",	 &camid);
  setptr("catdist",	 &catdist);
  setptr("catmid",	 &catmid);
  setptr("catprox",	 &catprox);
  setptr("catsoma",	 &catsoma);
  setptr("kprox",	 &kprox);
  setptr("kdist",	 &kdist);
  setptr("kmid",	 &kmid);
  setptr("ksoma",	 &ksoma);

  setptr("sbmaxrate",	 &sbmaxrate);
  setptr("sbmaxpool",	 &sbmaxpool);
  setptr("sbsynang",	 &sbsynang);
  
  setptr("sb_rm",	 &sb_rm);
  setptr("sb_rmp",	 &sb_rmp);
  setptr("sb_rmd",	 &sb_rmd);
  setptr("sb_vr",	 &sb_vr);
  setptr("sb_vs",	 &sb_vs);
  setptr("sb_rid",	 &sb_rid);
  setptr("sb_rii",	 &sb_rii);
  setptr("sb_rip",	 &sb_rip);
  setptr("sb_cmp",	 &sb_cmp);
  setptr("sb_cmd",	 &sb_cmd);
  setptr("sbaclm",       &sbaclm);
  setptr("sdia",         &sdia);
  setptr("spdia",        &spdia);
  setptr("sddia",        &sddia);
  setptr("spdia_ratio",  &spdia_ratio);
  setptr("smdia",        &smdia);
  setptr("sndia",        &sndia);

  setptr("orad1",        &orad1);
  setptr("irad2",        &irad2);
  setptr("orad2",        &orad2);
  setptr("irad3",        &irad3);
  setptr("orad3",        &orad3);
  setptr("irad4",        &irad4);
  setptr("orad4",        &orad4);
  setptr("irad5",        &irad5);
  setptr("orad5",        &orad5);
  setptr("irad6",        &irad6);
  setptr("orad6",        &orad6);

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
  setptr("rstim_theta",&rstim_theta);
  setptr("stim_theta", &stim_theta);
  setptr("noise_dur",  &noise_dur);
  setptr("nbase",      &nbase);
  setptr("nscale",     &nscale);
  setptr("vscale",     &vscale);
  setptr("ppdelay",    &ppdelay);

  setptr("dbp1_anpo", &dbp1_anpo);
  setptr("dbp1_anpi", &dbp1_anpi);
  setptr("dbp2_anpo", &dbp2_anpo);
  setptr("dbp2_anpi", &dbp2_anpi);

  setptr("dbp1_sgain",&dbp1_sgain);
  setptr("dbp1_vgain",&dbp1_vgain);

  setptr("sbac_isynanni",&sbac_isynanni);
  setptr("sbac_isynanpi",&sbac_isynanpi);
  setptr("sbac_isynanpo",&sbac_isynanpo);
  setptr("sbac_isynrngi",&sbac_isynrngi);
  setptr("sbac_dens",    &sbac_dens);

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
  setptr("datfile",    &datfile);

  nvalfile      = "nval_dsgc_sbac_noise.n";
  sbac_densfile = "dens_sbac_chans.n";
  if (notinit(arb_scale)) arb_scale = 1.0;
  chanparamsfile = "chanparams_dsgc_chans";

  dbp1_anpi = 20;
  dbp1_anpo = 90;
  dbp2_anpi = 90;
  dbp2_anpo = 110;

  sbmaxpool = 10;
  sbmaxrate = 20;
  sbsynang = 0;		// dendritic orientation for synapses onto dsgc

  dbp1_sgain = 3;  	// for noise (stimtype 7) set = 6
  dbp1_vgain = 0;  	// for noise (stimtype 7) set = 8

  sbac_isynanni = 80;		// inner radius of inhib annulus in presyn cell (for nval file)
  sbac_isynanpi = 0;		// inner radius of inhib annulus in postsyn cell (for nval file)
  sbac_isynanpo = 60;		// outer radius of inhib annulus in postsyn cell (for nval file)
  sbac_isynrngi = -1300;	// range of angles for postsyn cell (opposite, but 300 deg range)

  sbac_dens = 500;	// density of sbacs in nval_dsgc_sbac_noise.n

  _CA_T = _CA6;         /* set type of T-type calcium channel for dens_sbac_chans.n */
			/* then set it to slowly inactivating in chanparams_dsgc_chans */
  sdia = 0.6;
  spdia_ratio = 0.9;
  spdia = sdia * spdia_ratio;
  sddia = sdia;
  sndia = spdia;
  smdia = spdia;

  make_sbac_sbac = 1;

  //fprintf (stderr,"makestim %d\n",makestim);
  sb_init();
}

/*------------------------------------------------------*/

#include "gprim.h"

void syn_draw2 (synapse *spnt, int color, double vrev, double dscale, double dia,
                           double length, double foreshorten, int hide)
{
    int fill=1;
    double tlen;
    char tbuf[10];

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

  // sprintf (tbuf,"%d>%d",spnt->node1b,spnt->node2b);     /* print pre and postsynaptic cell number */
  gpen (blue);
  sprintf (tbuf,"%d %d",spnt->node2b,spnt->node2c);     /* print postsynaptic cell number */
  tlen = strlen(tbuf);
  gmove (length/2.0 -tlen*0.3, -1.0);
  gcwidth (2.5*dscale);
  gtext (tbuf);
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
  make_dbp2  = 1;
  make_hbp1  = 0;
  make_hbp2  = 0;
  make_ams   = 0;
  make_amhs  = 0;
  make_sbac  = 1;
  make_dsgc  = 1;

  make_dbp1_ams  = 0;
  make_hbp1_amhs = 0;

  ct=sbac;

  // setn(ct,NCOLOR,RCOLOR);
  set_synapse_dr (syn_draw2);
  if (!notinit(soma_z)) setn(ct,SOMAZ,soma_z);			// -44 for morph_R3RB140523_01

  SOMA      = R_4;
  if (n_dsgc<0) n_dsgc = 0;

#define SBARR 20

  if (notinit(ath_dia)) ath_dia = 0.4;	/* diameter of initial thin axon segment */
  if (notinit(dsomadia)) dsomadia = 10;	/* diameter of sbac soma */

  if (n_sbac>0) make_ct(sbac);		/* make sbacs if user specifies */

  if (!notinit(no_inhib)) { if (no_inhib==1) setsv(sbac,SCOND,3,0); } 
  if (notinit(revdir)) revdir = 0;	/* reverse pref and null directions */ 

  if (notinit(btrend)) btrend = 0;	/* linear baseline trend */
  if (notinit(ntrend))  ntrend = 0;	/* sine baseline trend */
  if (notinit(nperiod)) nperiod = 16.6666e-3; /* sine phase baseline trend */
  if (notinit(nphase))  nphase = 0;	/* sine phase baseline trend */

  if (notinit(dvrev)) dvrev = -0.06;
  if (notinit(dvst))  dvst  = -0.06;
 
  if (!notinit(ampa_cond)) setsv(dbp1,SCOND,2,ampa_cond);
  if (!notinit(ampa_cond)) setsv(hbp1,SCOND,4,ampa_cond);
  nmda_cond = 0;
  if (!notinit(nmda_cond)) setsv(dbp1,SCOND,7,nmda_cond);
  if (!notinit(nmda_cond)) setsv(hbp1,SCOND,5,nmda_cond);
  if (!notinit(gaba_cond)) setsv(ams, SCOND,1,gaba_cond);
  if (!notinit(gaba_cond)) setsv(amhs,SCOND,1,gaba_cond);
  if (notinit(dendrm)) dendrm = drm;
  if (notinit(dendcm)) dendcm = dcm;
  if (notinit(somacm)) somacm = dcm;
  if (notinit(tipcap)) tipcap = 0;
  if (notinit(dendrs)) dendrs = 1.0/dri;		// 65000 for dri=200 to predict currents;

  if (notinit(sb_rm))   sb_rm   = drm;          /* sbac soma Rm,     used in dens_sbac_chans.n */
  if (notinit(sb_rmp)) sb_rmp   = sb_rm;        /* proximal sbac Rm, used in dens_sbac_chans.n */
  if (notinit(sb_rmd)) sb_rmd   = sb_rm;        /* distal sbac Rm,   used in dens_sbac_chans.n */
  if (notinit(sb_vs)) sb_vs   = -0.07;          /* sbac vstart,      used in dens_sbaca.n */
  if (notinit(sb_vr)) sb_vr   = -0.07;          /* sbac vrev,        used in dens_sbaca.n */
  if (notinit(sb_rid)) sb_rid   = dri;
  if (notinit(sb_rii)) sb_rii   = dri;
  if (notinit(sb_rip)) sb_rip   = dri;

  if (notinit(sb_cmd)) sb_cmd   = dcm;
  if (notinit(sb_cmp)) sb_cmp   = sb_cmd;

  if (notinit(nadist)) nadist   = 0e-3;         /* distal Na density, dens_sbac_chans.n */
  if (notinit(namid))   namid   = 0e-3;         /* middle Na density, dens_sbac_chans.n */
  if (notinit(kdist))   kdist   = 0e-3;         /* dist Kdr density, dens_sbac_chans.n */
  if (notinit(kmid))     kmid   = 0e-3;         /* mid  Kdr density, dens_sbac_chans.n */
  if (notinit(kprox))   kprox   = 0e-3;         /* prox Kdr density, dens_sbac_chans.n */
  if (notinit(ksoma))   ksoma   = 0e-3;         /* soma Kdr density, dens_sbac_chans.n */
  if (notinit(cadist)) cadist   = 0e-3;         /* dist Ca density, dens_sbac_chans.n */
  if (notinit(camid))   camid   = 0e-3;         /* mid  Ca density, dens_sbac_chans.n */
  if (notinit(catdist)) catdist = 0e-3;         /* dist Ca-T density, dens_sbac_chans.n */
  if (notinit(catmid))   catmid = 0e-3;         /* mid  Ca-T density, dens_sbac_chans.n */
  if (notinit(catprox)) catprox = 0e-3;         /* prox Ca-T density, dens_sbac_chans.n */
  if (notinit(catsoma)) catsoma = 0e-3;         /* prox Ca-T density, dens_sbac_chans.n */

  if (notinit(sbaclm)) sbaclm   = 0.1;

  if (strcmp(sbac_file,"morph_sb1")==0) {
    setsv (sbac,SYNANNI,3,75);                    /* don't connect close to soma */
  }
  else if (strcmp(sbac_file,"morph_sbac3b")==0) {
    setsv (sbac,SYNANNI,3,20);                    /* don't connect close to soma */
  }

  set_chancalc(K,1,0,calck1nx);		/* sets the 0 (n) param for K type 1 */
  if (notinit(kexp)) kexp = 1;		/* exponent for K activation in calck1nx() above */
  dsintinc = 1e-4;

  if (notinit(datfile)) datfile = NULL;

   //display_z(zmax=-5, zmin=-15);	    /* display sublamina a */
   //display_z(zmax=-15, zmin=-25);	    /* display sublamina b */
   //
#define NGCS 10

   gcxarr   = (double *)emalloc(NGCS*sizeof(double));
   gcyarr   = (double *)emalloc(NGCS*sizeof(double));
   // gcxarr[0] = -150;
   gcxarr[0] = 0;
   gcyarr[0] = 0;

}

/*--------------------------------------------------------*/

void setdens (void)
/* Run by retsim after channel densities and properties have been read in. */

{
     int ct;

  if (!notinit(kdr_voff)) set_chan_offset(ct=sbac, C_K1 , CHOFFM, kdr_voff); /* set KDR voltage offset */
}

/*--------------------------------------------------------*/

void runonexit (void)
{
//        if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/


bool flag = false;
int cellpair=0;
double  ct = sbac;

double maxCurrent;
double cond;
double voltage;
double current, current2;
double idiff;

void onplot(void) {
    current  = i(ndn(ct,  1, elnode));
    if (cellpair) current2 = i(ndn(ct, 2, elnode));
    else          current2 = 0;
    idiff = current - current2;
    voltage  = v(ndn(ct,  1, soma));
}

/*------------------------------------------------------*/

double time1, time2;
double ixoff, iyoff;
double ncycles = 0;
double *rndarr = NULL;

double move_stim(double stimtime,double stimdur,double celldia, double roffset, double barwidth,double theta,double velocity, double ncycles, double econtrast, double icontrast, double eincr, double iincr, int direction, double mask)

/* Functions to generate light stimuli */
/*  returns the duration of the stimulus */

{
    static int runyet=0;
    double lbar, rbar, irad, orad, idia, odia;
    double inten, start, dur, wavel;
    double s,send;
    double sdia;
    double time2;

 // stim_spot(100, 100,0, inten=minten,  start=0.02, dur=0.05, wavel=1, 0);

 if (stimtype==0) {			// short pulse
     time2 = stimdur;
 } 
 else if (stimtype==1) {			// move bar to & fro
   lbar = -(celldia + barwidth) * 0.6 + roffset;
   rbar =  (celldia + barwidth) * 0.6 + roffset;

   if (revdir) {
	   double temp;
	   temp = lbar;
	   lbar = rbar;
	   rbar = temp;
   }

   // pref: move bar from left to right
 
   time1 = movebar (stimtime, 0,0,         lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);
   if (light_inhib) time1 = movebar (stimtime, ixoff,iyoff, lbar+ioffset, rbar+ioffset, 
		   					barwidth,barlength,theta,velocity,icontrast);
   // null: move bar from right to left
  
   time2 = movebar (time1+stimtime,0,0, rbar, lbar, barwidth,barlength,theta,velocity,econtrast);

   // make null inhibition stronger
   
   if (light_inhib) time2 = movebar (time1+stimtime,ixoff,iyoff, rbar+ioffset, lbar+ioffset, 
		   					   barwidth,barlength,theta,velocity,icontrast+iincr);
 }
 else if (stimtype==2) {			// move bar from left to right
   lbar = -(celldia + barwidth) * 0.6 + roffset;
   rbar =  (celldia + barwidth) * 0.6 + roffset;

   if (revdir) {
	   double temp;
	   temp = lbar;
	   lbar = rbar;
	   rbar = temp;
   }

   // pref: move bar from left to right
 
   time2 = movebar (stimtime, 0,0,         lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);
   // if (light_inhib) time1 = movebar (stimtime, ixoff,iyoff, lbar+ioffset, rbar+ioffset, 
   // 		   					barwidth,barlength,theta,velocity,icontrast);
   // // null: move bar from right to left
   // 
   // time2 = movebar (time1+stimtime,0,0, rbar, lbar, barwidth,barlength,theta,velocity,econtrast);
   //
   // // make null inhibition stronger
   //  
   // if (light_inhib) time2 = movebar (time1+stimtime,ixoff,iyoff, rbar+ioffset, lbar+ioffset, 
   // 		   					   barwidth,barlength,theta,velocity,icontrast+iincr);
  
 }
 else if (stimtype==3) {		// static annulus, radius = stimx
     //       stim_spot(100, 0, 0, econtrast, stimtime, spotdur);
        time1 = moveannulus(stimtime,        0,     0,     stimx+50,  stimx+55,  barwidth, velocity, econtrast);
        time2 = moveannulus(stimtime+0.1,    0,     0,     stimx+120, stimx+125, barwidth, velocity, econtrast);
    if (light_inhib) {
        time1 = moveannulus(stimtime,       ixoff,  iyoff, stimx+50,  stimx+55,  barwidth, velocity, icontrast);
        time2 = moveannulus(stimtime+0.1,   ixoff,  iyoff, stimx+120, stimx+125, barwidth, velocity, icontrast);
    }
    time2 = 0.2;

 }
 else if (stimtype==4) {		// move spot
	  sdia = barwidth;
 	  send = int(2000/sdia + 0.5);
	  for (s=0; s<send; s++) {
            stim_spot(sdia, s*sdia, 0, econtrast, stimtime+s*spotdur, spotdur);
          }
	  time2 = spotdur * send;
       }
  else if (stimtype==5) {		// sine waves 2,5,10,20,50 hz 
          stim_sine (speriod, sphase, orient, 0, 0, tfreq, drift, 1, scontrast, sqwave, 
			  				start=stimtime,  stimdur);
	  time2 = stimtime+stimdur;
  }
  else if (stimtype==6) {                 // checkerboard
           int npixels, nframes;
           double xoffset,yoffset;

       stimdur = 1;
       npixels = 16;
       stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
                          tfreq=15, 0, scontrast*2, stimtime, stimdur, &rndarr, &nframes, rseed);
       time2 = stimtime + stimdur;

       //   write out file containing random stimulus
       //
       // sprintf (checkerboard_file,"rnd_chck_file_%d",rseed);
       // if ((fchk=fopen(checkerboard_file,"w"))==NULL) { /* open file */
       //     ncfprintf (stderr,"Error creating checkerboard file\n");
       // } else {
       //     ncfprintf (stderr,"Creating checkerboard file %s\n",checkerboard_file);
       //     fwrite (rndarr,sizeof(double),npixels*npixels*nframes,fchk);
       // }
       efree(rndarr);
  }
  else if (stimtype==7) {                 // moving bar + checkerboard
           int npixels, nframes;
           double xoffset,yoffset;

    lbar = -(celldia + barwidth) * 0.6 + roffset;
    rbar =  (celldia + barwidth) * 0.6 + roffset;
    // time1 = movebar (stimtime+noise_dur, 0,0, lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);
    // time2 = movebar (time1+stimtime,     0,0, rbar, lbar, barwidth,barlength,theta,velocity,econtrast);
    time2 = movebar (stimtime+noise_dur, 0,0, lbar, rbar, barwidth,barlength,theta,velocity,econtrast+eincr);

    stimdur = time2 + noise_dur * 0.5;
    // stimdur = 1.0;
 
    if (noise_dur > 0) {
        npixels = 16;
        stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
                     tfreq=15, 0, scontrast*0.5, stimtime, stimdur, &rndarr, &nframes, rseed);
        efree(rndarr);
    }
    time2 = stimdur;
  }
  return time2;
}

/*--------------------------------------------------------*/

double vdiff (double nod1, double nod2, double time)

{
    int ct=sbac,cn=1;
    double v1, v2;

  v1 = v(ct,cn,int(nod1));
  v2 = v(ct,cn,int(nod2));
  return (v2 - v1)/(dri*dendrs);
}

/*--------------------------------------------------------*/

void addlabels(void)
{
    int cn;
    int nsynap, nsynapa, nsynapg, nsynapn, nsynapi, nsynapsb;
    node *npnt;
    photorec *p;

   /* add light transducer to each bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }
   for(npnt=nodepnt; npnt=foreach(npnt,dbp2,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp2,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }

   if (notinit(theta))   theta = 0;	/* orientation of bar */
   if (notinit(iroff))   iroff = 1000;	/* r offset for inhib transducers */

   ixoff = iroff * cosdeg(theta);
   iyoff = iroff * -sindeg(theta);

   if (notinit(light_inhib)) light_inhib = 0; 
   if (light_inhib) {

     /* Add light transducer to each small-field amacrine cell */
     /*  offset so that inhibition can be controlled separately */

     for(npnt=nodepnt; npnt=foreach(npnt,ams,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       p = (photorec*)make_transducer(ndn(ams,cn,soma)); 
       p->xpos=npnt->xloc + ixoff; 
       p->ypos=npnt->yloc + iyoff;
     }
     for(npnt=nodepnt; npnt=foreach(npnt,amhs,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       p = (photorec*)make_transducer(ndn(amhs,cn,soma)); 
       p->xpos=npnt->xloc + ixoff; 
       p->ypos=npnt->yloc + iyoff;
     }
   }

  /*  - - - - - - - - - - - - - - - - - - - */

   /* make lists of synaptic inputs to the sbac */

    nsynapa   = synapse_add  (1,dbp1,-1,-1,sbac,-1,3);      /* make list of dbp1 ampa syns to sbac */
    nsynapg   = synapse_add  (2,dbp1,-1,-1,dsgc, 1,2);      /* make list of dbp1 ampa synapses to dsgc */
    nsynapi   = synapse_add  (3,sbac,-1,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
    nsynapsb  = synapse_add  (4,sbac,-1,-1,sbac,-1,3);      /* make list of inhib synapses from sbacs to sbacs */
//    nsynapi += synapse_add (1,amhs,-1,-1,sbac,1,1);       /* make list of inhib synapses from amhs's */
//    nsynapa += synapse_add (2,hbp1,-1,-1,sbac,1,2);       /* make list of hbp1 ampa syns for recording below */
//    nsynapn  = synapse_add (7,dbp1,-1,-1,sbac,1,7);       /* make list of dbp1 nmda synapses */
//    nsynapn += synapse_add (7,hbp1,-1,-1,sbac,1,5);       /* make list of hbp1 nmda synapses */
    nsynapn = 0;
    nsynap = nsynapa + nsynapg + nsynapi + nsynapsb;
    // ncfprintf (stderr,"# nsynap %d nsyn_dbp1->sbac %d nsyn_dbp1->dsgc %d nsyn_sbac_dsgc %d nsyn_sbac->sbac %d\n",
    //		    nsynap, nsynapa, nsynapg, nsynapi, nsynapsb);
    ncfprintf (stderr,"# nsyn_dbp1->sbac %4d\n",nsynapa); 
    ncfprintf (stderr,"# nsyn_dbp1->dsgc %4d\n",nsynapg); 
    ncfprintf (stderr,"# nsyn_sbac->dsgc %4d\n",nsynapi); 
    ncfprintf (stderr,"# nsyn_sbac->sbac %4d\n",nsynapsb); 
    ncfprintf (stderr,"# nsyn total      %4d\n",nsynap); 
    ncfprintf (stderr,"#\n"); 
}

/*--------------------------------------------------------*/

void runexpt(void)

{
#define NUMCBPS 12
    int c, i, ct, cn, colr, pl, plnum, electrode_node, dist_node;
    double start, starttime, dur, dscale, plsize;
    double disp_end, mask; 
    double rmax, rmin;
    double Imax, Imin, imax, imin, gmax, Vmax, Vmin;
    double celldia, barloc;
    double cmin, cmax;
    double r = 20;		/* radius increment for plots */
    elem *e;
    chattrib *a;
    char sbuf[30];

  //ploti = 1e-3;
  ploti = 1e-4;
  timinc = 2e-6;

  gcdistnod = 582;
  electrode_node = 5000;

  Vmax  = -0.02;
  Vmin = -0.07;
  
  /*  - - - - - - - - - - - - - - - - - - - */

   if (notinit(barwidth))        barwidth = 110;
   if (notinit(barlength))      barlength = 500;
   if (notinit(minten))            minten = -0.058;
   if (notinit(econtrast))      econtrast =  0.017;
   if (notinit(scontrast))      scontrast =  econtrast;
   if (notinit(icontrast))      icontrast =  scontrast;
   if (notinit(eincr))              eincr =  0;
   if (notinit(sincr))              sincr =  0.002;
   if (notinit(iincr))              iincr =  sincr;
   if (notinit(velocity))        velocity =  2000; 
   if (notinit(stimscale))      stimscale =  1; 
   if (notinit(stimdur))          stimdur =  0.5;	   /* used for non-moving stimuli */
   if (notinit(prestimdur))    prestimdur =  0.05;
   if (notinit(poststimdur))  poststimdur =  0.005;
   if (notinit(tailcurdur))    tailcurdur =  0.0;
   if (notinit(noise_dur))      noise_dur =  1; 
   if (notinit(stimtime))        stimtime =  prestimdur;
   if (notinit(disptime))        disptime =  0.15;
   if (notinit(stimtype))        stimtype =  2;
   if (notinit(spotdur))          spotdur =  0.05;

   if (notinit(vhold))              vhold = -0.07;
   if (notinit(vstart))            vstart = -0.70;
   if (notinit(vstop))              vstop =  vstart;
   if (notinit(vstep))              vstep =  0.005;
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
   if (notinit(stim_theta))    stim_theta = 0;
   if (notinit(rstim_theta))  rstim_theta = 0;
    
   // set_run_on_exit(runonexit);                         // set to erase savefile on ^C
   // sprintf (savefile,"sbac_noise%06d",getpid());       // add pid to file name

   if (notinit(ioffset)) ioffset = -barwidth;

   if (notinit(elec_rs))  elec_rs   = 20e6;
   if (notinit(elec_cap)) elec_cap  = 1e-14;
   if (notinit(shunt_res)) shunt_res  = 1e30;
   if (notinit(shunt_v))   shunt_v  = dvrev;
   if (notinit(elnode))     elnode  = electrode_node;

   if (notinit(sb_cn)) sb_cn = 1;
   if (sb_cn > n_sbac) sb_cn = 1;

   if (elnode==electrode_node) {
       make_electrode  (nd(sbac,sb_cn,elnode), ndn(sbac,sb_cn,soma), elec_rs, elec_cap, tipcap);
       make_shunt  (ndn(sbac,sb_cn,soma), shunt_res, shunt_v);
   }
       
   cn = 1;

   Vmin = -0.049; Vmax = -0.045;
   plot_v_nod(dbp1,findmida(dbp1,50,-150),soma,-0.045, -0.03, 1, "", 22, 0.3);

   if (stimtype==2)
         barloc = 100;
   else barloc = 50;

	     
    plot_syncond(findsynloca(dbp1,0,0),              cmin=0,cmax=200e-12, 1, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*1,rstim_theta),  cmin=0,cmax=200e-12, 2, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*2,rstim_theta),  cmin=0,cmax=200e-12, 3, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*3,rstim_theta),  cmin=0,cmax=200e-12, 4, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*4,rstim_theta),  cmin=0,cmax=200e-12, 5, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*5,rstim_theta),  cmin=0,cmax=200e-12, 6, 20,"",0.5);
    plot_syncond(findsynloca(dbp1,r*6,rstim_theta),  cmin=0,cmax=200e-12, 7, 20,"",0.5);

    if (notinit(exptrun)) exptrun = 1;

    Vmin = -0.07; Vmax = -0.02;

    //plot_v_nod(sbac, sb_cn, findnodloc(sbac,sb_cn,-40,0),  Vmin, Vmax, colr=blue,  "", 10, 1);
    plot_v_nod(sbac, sb_cn, elnode,                          Vmin, Vmax, colr=white, "Vsbac_elec", 10, 1);
    plot_v_nod(sbac, sb_cn, soma,                            Vmin, Vmax, colr=green, "", 10, 1);
    plot_v_nod(sbac, sb_cn, findnodlocr(sbac,sb_cn,40,0),   Vmin, Vmax, colr=cyan,  "", 10, 1);
    plot_v_nod(sbac, sb_cn, findnodlocr(sbac,sb_cn,80,0),   Vmin, Vmax, colr=brown, "", 10, 1);
    plot_v_nod(sbac, sb_cn, findnodlocr(sbac,sb_cn,120,0),  Vmin, Vmax, colr=red,   "", 10, 1);

    // plot_v_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,0),  Vmin, Vmax, colr=magenta,   "", 10, 1);
    // plot_v_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,45),  Vmin, Vmax, colr=green,   "", 10, 1);
    // plot_v_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,90),  Vmin, Vmax, colr=blue,   "", 10, 1);
    // plot_v_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,135),  Vmin, Vmax, colr=brown,   "", 10, 1);
    // plot_v_nod(sbac, sb_cn, 2840,  Vmin, Vmax, colr=red,   "", 10, 1);


    plot_ca_nod(sbac, sb_cn, findnodlocr(sbac,sb_cn,120,0),     20e-6, colr=cyan,     "", 8, 0.5);
    plot_ca_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,0),    20e-6, colr=magenta,  "", 8, 0.5);
    plot_ca_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,45),   20e-6, colr=green,    "", 8, 0.5);
    plot_ca_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,90),   20e-6, colr=blue,     "", 8, 0.5);
    plot_ca_nod(sbac, sb_cn, findnodlocra(sbac,sb_cn,150,135),  20e-6, colr=brown,    "", 8, 0.5);
    plot_ca_nod(sbac, sb_cn, 2840,                              20e-6, colr=red,      "", 8, 0.5);

    Imin = 0e-12;    Imax = 5e-10; 
    if (exptrun==1 || exptrun==4 || exptrun==5) { 
        plot_i_nod(sbac,sb_cn,elnode,  Imin=-2e-10,Imax=2e-10,green, "Isbac_1_elec", 7, 0.5);
    }
    if (n_dsgc > 0) 
        plot_i_nod(dsgc,cn=1,soma,    Imin=-5e-9,Imax=5e-9,green, "Isoma_dsgc  ", 2, 0.5);

     // Vmin = -0.07; Vmax = -0.02;
     // plot_v_nod(sbac, sb_cn, soma,        	   Vmin, Vmax, green, "", 5, 0.35);
     // plot_v_nod(sbac, sb_cn, findnodloc(sbac,cn,40,0),   Vmin, Vmax, colr=cyan,  "", 5, 1);
     // plot_v_nod(sbac, sb_cn, findnodloc(sbac,cn,80,0),   Vmin, Vmax, colr=brown, "", 5, 1);
     // plot_v_nod(sbac, sb_cn, findnodloc(sbac,cn,120,0),  Vmin, Vmax, colr=red,   "", 5, 1);

    imax =  0e-10;
    imin = -3e-10;
    gmax =  5e-9;
    plot_func(isyn_tot,1,imax,imin);       plot_param("Itot_bpa  ",blue,4,0.3);
    plot_func(gsyn_tot,1,gmax,0);          plot_param("Gtot_bpa  ",blue, 3,0.3);
    // plot_func(gsyn_tot,7,gmax,0);          plot_param("Gtot_bpn  ",cyan, 2,0.3);
    // plot_func(gsyn_tot,1,gmax,0);          plot_param("Gtot_ami  ",red,  2,0.3);
    if (n_dsgc > 0) {
      plot_func(isyn_tot,3,50e-12,0);          plot_param("Iinh_dsgc  ",red,  1,0.5);
    }

   celldia = max(xarrsiz,yarrsiz) * 1.5 * stimscale; 
   // fprintf (stderr,"# celldia %g\n",celldia);

   if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
   if (notinit(predur))        predur = 0.05;

   /* - - - - - - - - - - - - - - */

   if (disp==16) {
	double t;

      simtime = 0;                                    // must be set ahead of stim_backgr()
      stim_backgr(minten,start=simtime);			 /* turn on  background */

      stimdur = move_stim(stimtime, stimdur, celldia, stimr, barwidth, theta, velocity, ncycles, 
		      		econtrast, icontrast, eincr, iincr, direction, mask=1);
      // if (light_inhib) display_size(2500); else display_size(300); 
      display_size(500); 
      disp_end = stimtime+prestimdur+stimdur+poststimdur;
      for (starttime=0,t=stimtime; t<disp_end; starttime = t, t+= 0.002) {
           display_stim(starttime, t, dscale=4, -0.035, -0.045);
           //display_stim(t, dscale=4, -0.035, -0.045); 
 	   //display_stim(0+t, dscale=4); 
           simwait(0.10);
      }
      return;
   }
   /* - - - - - - - - - - - - - - */

   simtime = -predur;                                    // must be set ahead of stim_backgr()
   stim_backgr(minten,start=simtime);			 /* turn on  background */

  /* run experiment */

  if (exptrun == 1) {

     if (stimtype==1) noise_dur = 0;

     stimdur = move_stim(stimtime, stimdur, celldia, stimr, barwidth, theta, velocity, ncycles, 
			 econtrast, icontrast, eincr, iincr, direction, mask=1);

        // fprintf (stderr,"prestimdur %g stimdur %g poststimdur %g endexp %g\n",
       //  		prestimdur,stimdur,poststimdur,endexp);
  
     endexp=prestimdur+stimdur+poststimdur;

     vclamp (ndn(sbac,sb_cn, elnode), vhold, simtime,  predur);
     step(predur);

     vclamp (ndn(sbac,sb_cn, elnode), vhold, simtime,  endexp);
     step (prestimdur);
     // savemodel (savefile);
     step (stimdur + poststimdur);

     // restoremodel (savefile);
     // step (stimdur/2 - noise_dur/2 + noise_dur/4 + poststimdur);

 }    /* run 1 */ 

  /* -  -  -  -  -  -  -  -  -  - */

 else if (exptrun==2) {		/* cclamp, just display responses to moving bar */

	// fprintf(stderr,"endexp %g %g %g\n", endexp, time1, time2);

    if (stimtype==1) noise_dur = 0;

    stimdur = move_stim(stimtime, stimdur, celldia, stimr, barwidth, theta, velocity, ncycles, 
			econtrast, icontrast, eincr, iincr, direction,mask=1);
    // fprintf (stderr,"restore %g\n",noise_dur/4+stimdur/2+poststimdur);

    endexp = prestimdur+stimdur+poststimdur;
    step(predur);

    step (prestimdur);
    // savemodel (savefile);
    step (stimdur + poststimdur);

    // restoremodel (savefile);
    // step (stimdur/2 - noise_dur/2 + noise_dur/4 + poststimdur);

  } 

 /* -  -  -  -  -  -  -  -  -  - */

  else if (exptrun==3) {	/* run noise stimulus */

    stimdur = move_stim(stimtime, stimdur, celldia, stimr, barwidth, theta, velocity, ncycles, 
			econtrast, icontrast, eincr, iincr, direction,mask=1);
    endexp=prestimdur+stimdur+poststimdur;
    step(predur);
    if (n_dsgc > 0) vclamp (ndn(dsgc,1,soma), 0, simtime,  10.0);
    step (prestimdur);
    step (stimdur + poststimdur);

  } 

 /* -  -  -  -  -  -  -  -  -  - */

  else if (exptrun==4) {	/* vclamp, short pulse, noise */
      double t, randv, stimdur2;
      double  tstep = 0.002;
      double pulseampl = 0.02;


    //stimdur2 = 1.0;
    stimdur2 = stimdur;
    endexp = 4*prestimdur + stimdur + stimdur2;
    step(predur);
    if (n_dsgc > 0) vclamp (ndn(dsgc,1,soma), 0, simtime,  10.0);
    step (prestimdur);
    vclamp (ndn(sbac,sb_cn,elnode), vhold+pulseampl, simtime,  stimdur);
    step (stimdur);
    step (prestimdur);
    for (t=0; t<stimdur2; t+=tstep) {
        randv = pulseampl * 1.0 * drand();
        vclamp (ndn(sbac,sb_cn, elnode), vhold + randv, simtime,  tstep);
	step (tstep);
    }
    step (prestimdur);

  } 

 /* -  -  -  -  -  -  -  -  -  - */

  else if (exptrun==5) {	/* read file, vclamp */
	  int i, nlong, nwid;
          double *datarr;
	  double val;
	  double tstep = ploti;
#define TIMESCAL 10

    if (notinit(nbase)) nbase = -0.07;
    if (notinit(nscale)) nscale = 1.0;
    if (notinit(vscale)) vscale = 1.0;
    if (datfile != NULL) {
        if ((datarr=fread (datfile,&nlong,&nwid)) == NULL) {
             ncfprintf (stderr,"expt_sbac_noise: error reading data file '%s'\n",datfile);
        }
        else { 
             endexp = nlong * tstep;
             step(predur);
             if (n_dsgc > 0) vclamp (ndn(dsgc,1,soma), -0.07, simtime,  10.0);
             for (i=0; i<nlong; i+=TIMESCAL) {
	          val = (datarr[i*2+1] * vscale - nbase) * nscale + nbase;
		  // fprintf (stderr,"t %g val %g\n",i*tstep,val);
                  vclamp (ndn(sbac,sb_cn, elnode), val, i*tstep,  tstep*TIMESCAL);
             }
 	     step (endexp);
        }
    }
  } 

 /* -  -  -  -  -  -  -  -  -  - */

  else if (exptrun==6) {	/* paired-pulse depression, vclamp */
      // double pulseampl = 0;
      double pulseampl = 0.06;
      double pulsedur = 0.02;

    if (notinit(ppdelay)) ppdelay = 0.2;
    stimdur = 1.0;
    // endexp = prestimdur + 2 * pulsedur + ppdelay + poststimdur;
    endexp = prestimdur + 1.0;
    vclamp (ndn(sbac,sb_cn,elnode), vhold, simtime, predur);
    if (n_dsgc > 0) vclamp (ndn(dsgc,1,soma), 0, simtime,  10.0);
    step(predur);
    vclamp (ndn(sbac,sb_cn,elnode), vhold, simtime, prestimdur);
    step (prestimdur);
    vclamp (ndn(sbac,sb_cn,elnode), vhold+pulseampl, simtime,  pulsedur);
    step (pulsedur);
    vclamp (ndn(sbac,sb_cn,elnode), vhold, simtime, ppdelay);
    step (ppdelay);
    vclamp (ndn(sbac,sb_cn,elnode), vhold+pulseampl, simtime,  pulsedur);
    step (pulsedur);
    vclamp (ndn(sbac,sb_cn,elnode), vhold, simtime, poststimdur);
    step (poststimdur);

  } 

  // unlink (savefile);
  // savefile[0] = 0;
  fflush(stdout);
  printf ("# done\n");		// for the "modelfit" program to see when the simulation is done
}
