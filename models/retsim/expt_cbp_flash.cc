/* Experiment cbp_flash */
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

double kdr_cond;
double axon_br_dia;
double varicos_dia;

int do_vclamp;
double vclamp_dur;
double clamp_voltage;

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

  setptr("kdr_cond",    &kdr_cond);
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);

  setptr("do_vclamp", &do_vclamp);
  setptr("vclamp_dur",&vclamp_dur);
  setptr("clamp_voltage",&clamp_voltage);

  nvalfile = "nval_cbp_flash.n";

  // dbp1_file = "morph_bp";		// default set in retsim.cc
  // dbp1_file = "morph_DB_111005_12";
  
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  make_cones= 1;        /* make cones, dbp, gc */
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 1;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 0;
  make_gcb  = 0;
  make_dsgc = 0;

  setn(dbp1,MORPH,0);		/* set dbp1 morphology from file, def = "morph_bp" */
  setn(dbp1,BIOPHYS,1);		/* set dbp1 biophys from file, def = "dens_dbp1" */
  setn(dbp1,NCOLOR,RCOLOR);       /* set cell display color from region */

  //if (notinit(arrsiz)) arrsiz = 300;
  if (notinit(dvrev)) dvrev = -0.06;            /* Vrev for Rm */
  if (notinit(dvst))   dvst = -0.06;            /* Vstart */

  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5; 	/* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5; 	/* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2; 		/* default dia for dendrites */
  
  if (notinit(cbplam)) cbplam = 0.1; 		/* default compartment size */
  // if (notinit(dispsize)) dispsize = 100; 	/* default display size */
  if (notinit(node_scale)) node_scale = -3.05;  /* 3: nodenum, 0.05: small font */

  // Set default channel types in density file.
  // To change, check manual, and uncomment and modify these:
  //
  //   _NA  = _NA2;
  //   _KA  = _K3;
  //   _KH  = _K4;
  //   _KIR = _K5;

  pnoise = 0;  		// =1 -> photon noise
  dnoise = 0;  		// =1 -> continuous dark noise
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{
    if (!notinit(kdr_cond))    celdens[dbp1][0][_KDR][SOMA] = kdr_cond;  
    //setsv (hb,SCOND,1, 0);
    // if (notinit(arrsiz)) arrsiz = 20;
    //if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 

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
    double Imin, Imax;

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
  if (notinit(scontrast)) scontrast = 3;       /* intensity increment */
  if (!notinit(setploti))     ploti = setploti; /* plot time increment */
  
  if (notinit(do_vclamp))        do_vclamp = 0;     /* =1 -> turn on vclamp */
  if (notinit(vclamp_dur))      vclamp_dur = 1.0;   /* duration of vclamp */
  if (notinit(clamp_voltage))clamp_voltage = -0.06; /* clamp voltage */

  midcone = findmid(xcone,0,0);
  midcbp  = findmid(dbp1,0,0);
  midcbp2 = findmid(dbp1,10,10);
  //midcbp2 = find_gtconn(dbp1, 8);
  synin1 = ncel_in(dbp1,midcbp,xcone);
  synin2 = ncel_in(dbp1,midcbp2,xcone);
  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d ncones %d\n",  midcbp,synin1);
  if (ninfo >=1) fprintf (stderr,"# mid cbp2 # %d ncones %d\n",  midcbp2,synin2);


  plot_v_nod(ct=xcone,cn=midcone,n=soma,Vmin=-.030,Vmax =-.025,colr=cyan,"", -1, -1); /* Vcones*/
  plot_synrate_out(ct=xcone,cn=midcone,rmin=0,rmax=400,colr=magenta);	              /* rate out */
  plot_v_nod(ct=dbp1,cn=midcbp,n=soma, Vmin=-.06,Vmax =-.035,colr=red,"", -1, -1);   /* Vcbp */
  if (do_vclamp)
    plot_i_nod(ct,cn=midcbp,n=soma, Imin= 1.0e-11, Imax = -5.0e-11,colr=magenta,"", 1, -1);   /* Icbp */

  //plot_v_nod(ct=dbp1,cn=midcbp2,n=soma,Vmin=-.045,Vmax =-.035,colr=green,"", -1, -1); /* Vcbp */
  // plot_synrate_out(ct=dbp1,cn=midcbp,rmin=0,rmax=200,colr=magenta);	              /* rate out */
  //plot_v_nod(ct=gca, cn=1,n=soma,Vmin=-.070,Vmax =-.050,colr=blue,"", -1, -1);      /* Vgc */
  // if (make_gca && getn(gca,BIOPHYS)) {plot(CA, 1, ndn(gca,1,soma), fmax=0.5e-6, fmin=0); 
  //			plot_param("Cai", colr=yellow,plnum=0,plsize=0.3);}


  if (do_vclamp) vclamp(ndn(ct=dbp1,midcbp,soma), clamp_voltage, simtime,  vclamp_dur);

  stim_backgr(minten);
  for (t=0; t<exptdur; t+= dtrial){
     double start, dur;

    stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=dstim);
    step(dtrial);
  }
}
