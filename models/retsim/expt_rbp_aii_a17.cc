/* Experiment rbp_aii_a17 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"


int ntrials;
int rodarr;
int rbparr;
int aiiarr;
int a17arr;

double temp_freq;
double dstim;
double sdia;
double stimtime;
double prestimdur;
double stimdur;
double poststimdur;
double minten;
double scontrast;
double setploti;
double dvrev;
double dvst;
double cbplam;
double istart;
double istop;
double istep;
double vstart;
double vstop;
double vstep;
double cstart;
double cstop;
double cstep;
double iaii;
double drm;

void vinit(void);

int rec_ct;
int rec_cn;

char savefile[30] = {0};


void defparams(void)
{
  setptr("rodarr",    &rodarr);
  setptr("rbparr",    &rbparr);
  setptr("aiiarr",     &aiiarr);
  setptr("a17arr",    &a17arr);

  setptr("rec_ct",     &rec_ct);
  setptr("rec_cn",     &rec_cn);
  setptr("temp_freq",  &temp_freq);
  setptr("ntrials",    &ntrials);
  setptr("dstim",      &dstim);
  setptr("sdia",       &sdia);
  setptr("stimtime",   &stimtime);
  setptr("prestimdur", &prestimdur);
  setptr("stimdur",    &stimdur);
  setptr("poststimdur",&poststimdur);
  setptr("minten",     &minten);
  setptr("scontrast",  &scontrast);
  setptr("setploti",   &setploti);
  setptr("dvst",       &dvst);
  setptr("dvrev",      &dvrev);
  setptr("cbplam",     &cbplam);
  setptr("istart",     &istart);
  setptr("istop",      &istop);
  setptr("istep",      &istep);
  setptr("vstart",     &vstart);
  setptr("vstop",      &vstop);
  setptr("vstep",      &vstep);
  setptr("cstart",     &cstart);
  setptr("cstop",      &cstop);
  setptr("cstep",      &cstep);
  setptr("iaii",       &iaii);
  setptr("drm",       &drm);
  
  nvalfile	    = "nval_rbp_aii_a17.n";
  rod_densfile	    = "dens_rod.n";
  rbp_densfile	    = "dens_rbp.n";
  aii_densfile	    = "dens_aii.n";
  a17_densfile	    = "dens_a17.n";
  gca_densfile	    = "dens_gca_aii_rbp.n";

  if (notinit(rbp_file))   rbp_file    = (char *)"morph_rbp_yoshi_cell14_normalized";     /* real rbp morphology from Tsukamoto & Omi, 2013 */
  if (notinit(aii_file))   aii_file    = (char *)"morph_aii_yoshi_cell2_normalized";     /* real aii morphology from Tsukamoto & Omi, 2013 */

  if (notinit(rbp_thetax)) rbp_thetax=55;
  if (notinit(rbp_thetay)) rbp_thetay=14;
  if (notinit(aii_thetax)) aii_thetax=48;
  if (notinit(aii_thetay)) aii_thetay=-150;	/* optimal x,y rotations for rbp and aii cells from Yoshi - determined by trial and error */

  if (notinit(rbp_nscale)) rbp_nscale=-3.1;
  if (notinit(rod_nscale)) rod_nscale=-3.1;
  if (notinit(a17_nscale)) a17_nscale=-3.1;
  if (notinit(aii_nscale)) aii_nscale=-3.1;	/* node display params */

}

/*------------------------------------------------------*/
#include "gprim.h"

void syn_draw2 (synapse *spnt, int color, double vrev, double dscale, double dia,
		                                double length, double foreshorten, int hide)
	
	/* draw synapse within small oriented frame */
	
{
    int fill=1;
    double tlen;
    char tbuf[10];
		    
    dia *= dscale;                        /* draw circle with line */
    if (dia < 0) dia = -dia;
	  color = -1;
    	if (color < 0) {
        	if (vrev < -0.04) color = RED;
	        else              color = BLUE;
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
			    gpen (black);
	 sprintf (tbuf,"%d %d %d %d",spnt->node1b,spnt->node1c, spnt->node2b, spnt->node2c);     /* print pre and postsynaptic cell number */
	 tlen = strlen(tbuf);
	 gmove (length/2.0 -tlen*0.3, 0);
	 gcwidth (0.5);
	 gtext (tbuf);
}

/*------------------------------------------------------*/
void setparams(void)

{
  make_rods  = 1;        /* make rods, rbc, aii, a17 */
  make_cones = 0;
  make_ha    = 0;
  make_hb    = 0;
  make_hbat  = 0;
  make_dbp2  = 0;
  make_rbp   = 1;
  make_aii   = 1;
  make_a17   = 1;
  make_gca   = 0;
  make_gcb   = 0;

  make_rod_rbp  = 1;
  make_rbp_aii  = 1;
  make_aii_aii  = 1;
  make_rbp_a17  = 1;
  make_a17_rbp  = 1;

  if (notinit (n_rods)) n_rods=1;
  if (notinit (n_rbp)) n_rbp=1;
  if (notinit (n_aii)) n_aii=1;
  if (notinit (n_a17)) n_a17=1;

  set_synapse_dr (syn_draw2);

  vna =  0.04;
  vcl = -0.060;
  drm = 30e3;
  if (notinit(drm)) drm = 30e3;
  dicafrac = 0;				// set ca pump current = 0

  _CA_T = C_CA7;

  dnoise = 0;
  pnoise = 1;

  if(notinit(rec_ct)) rec_ct = aii;
  if(notinit(rec_cn)) rec_cn = 1;
//  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 1000;

  if (notinit(dvst))   dvst = -0.065;
  if (notinit(dvrev)) dvrev = -0.065;
  if (notinit(cbplam)) cbplam = 0.005;

  if (notinit(rodarr)) rodarr = 1;
  if (notinit(rbparr)) rbparr = 1;
  if (notinit(aiiarr)) aiiarr = 1;
  if (notinit(a17arr)) a17arr = 1;

#define RODARR 30
  if (!notinit(rodarr)) {
    rodxarr  = (double *)emalloc(RODARR*sizeof(double));
    rodyarr  = (double *)emalloc(RODARR*sizeof(double));
    rodtharr = (double *)emalloc(RODARR*sizeof(double));
    rodnarr  = (int *) emalloc(RODARR*sizeof(int));
    if (rodarr==0) {
	n_rods = 0;
    }
    if (rodarr==1) {
      rodxarr[0] = 0;    rodyarr[0] = 0;    rodtharr[0] = 0;    rodnarr[0] = 1;
      rodarrsiz=8;
      n_rods = 4;
    }
    else if (rodarr==2) {
      rodxarr[0] =  0;   rodyarr[0] = 0;    rodtharr[0] = 0;    rodnarr[0] = 1;
      rodxarr[1] =  5;   rodyarr[1] = 0;    rodtharr[1] = 0;    rodnarr[1] = 2;
      n_rods = 4;
    }
  }
#define RBPARR 30
  if (!notinit(rbparr)) {
    rbpxarr  = (double *)emalloc(RBPARR*sizeof(double));
    rbpyarr  = (double *)emalloc(RBPARR*sizeof(double));
    rbptharr = (double *)emalloc(RBPARR*sizeof(double));
    rbpnarr  = (int *) emalloc(RBPARR*sizeof(int));
    if (rbparr==0) {
	n_rbp = 0;
    }
    if (rbparr==1) {
      rbpxarr[0] =    0;    rbpyarr[0] = 0;    rbptharr[0] = 0;    rbpnarr[0] = 1;
      n_rbp = 1;
    }
    else if (rbparr==2) {
      rbpxarr[0] =  0;   rbpyarr[0] = 0;    rbptharr[0] = 0;    rbpnarr[0] = 1;
      rbpxarr[1] =  5;   rbpyarr[1] = 0;    rbptharr[1] = 0;    rbpnarr[1] = 2;
      n_rbp = 2;
    }
  }
#define AIIARR 30
  if (!notinit(aiiarr)) {
    aiixarr  = (double *)emalloc(AIIARR*sizeof(double));
    aiiyarr  = (double *)emalloc(AIIARR*sizeof(double));
    aiitharr = (double *)emalloc(AIIARR*sizeof(double));
    aiinarr  = (int *) emalloc(AIIARR*sizeof(int));
    if (aiiarr==0) {
	n_aii = 0;
    }
    if (aiiarr==1) {
      aiixarr[0] =  3;      aiiyarr[0] = 0;       aiitharr[0] = 0;    aiinarr[0] = 1;
      n_aii = 1;
    }
    else if (aiiarr==2) {
      aiixarr[0] =  3;      aiiyarr[0] = 0;       aiitharr[0] = 0;    aiinarr[0] = 1;
      aiixarr[1] = -5;      aiiyarr[1] = 0;       aiitharr[1] = 0;    aiinarr[1] = 2;
      n_aii = 2;
    }
  }
#define A17ARR 30
  if (!notinit(a17arr)) {
    a17xarr  = (double *)emalloc(A17ARR*sizeof(double));
    a17yarr  = (double *)emalloc(A17ARR*sizeof(double));
    a17tharr = (double *)emalloc(A17ARR*sizeof(double));
    a17narr  = (int *) emalloc(A17ARR*sizeof(int));
    if (a17arr==1) {
      a17xarr[0] =  50;      a17yarr[0] = 0;       a17tharr[0] = 0;    a17narr[0] = 1;
      n_a17 = 1;
    }
    else if (a17arr==2) {
      a17xarr[0] =  40;      a17yarr[0] = 0;       a17tharr[0] = 0;    a17narr[0] = 1;
      a17xarr[1] = -40;      a17yarr[1] = 0;       a17tharr[1] = 0;    a17narr[1] = 2;
      n_a17 = 2;
    }
  }
}

/*------------------------------------------------------*/

void runonexit (void)
{
   if (savefile[0]!=0)  unlink (savefile);
}

/*-----------------------------------------------*/

void runexpt(void)

{
   int ct, cn, n, c, r, colr, i, v, pl, plnum;
   int aiin, dbpn, lob_append=1;
   int midrod, midrbp, mida17, midaii, middbp ;
   int rbp_recnode, aii_recnode, a17_recnode1, a17_recnode2;
   char rbp_n[20], aii_n[20], a17_n1[20], a17_n2[20];
   double s, t, start, dur, dtrial, exptdur;
   double max,min, rmax, rmin;
   double plsize,predur;
   double ipulse, cvolt, vhold, cont, sign;
   double cmax, cmin, Vmin, Vmax;

   node *npt;

  if (notinit(temp_freq)) temp_freq = 2;
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  dtrial = 1 / temp_freq;
  exptdur = dtrial * ntrials;
//  endexp  = exptdur;
  ploti = 1e-4;
  timinc = 1e-5;

  if (notinit(ntrials))         ntrials = 1;	/* number of trials */
  if (notinit(stimtime))       stimtime = .10;  /* stimulus time */
  if (notinit(dstim))             dstim = .05;	/* stimulus duration */
  if (notinit(sdia))               sdia = 30;	/* spot diameter */
  if (notinit(prestimdur))   prestimdur = .05;	/* before stimulus */
  if (notinit(stimdur))         stimdur = .002;	/* stimulus duration */
  if (notinit(poststimdur)) poststimdur = .2;	/* after stimulus */

  if (notinit(minten))           minten = bg_inten;	/* background intensity (for make cone)*/
  if (notinit(scontrast))     scontrast = 9;		/* intensity increment */
  if (!notinit(setploti))         ploti = setploti;	/* plot time increment */
  if (notinit(istart))           istart = 0e-12;
  if (notinit(istep))             istep = 5e-12;
  if (notinit(istop))             istop = 25e-12;
  if (notinit(vstart))           vstart = -0.04;
  if (notinit(vstep))             vstep = 0.002;
  if (notinit(vstop))             vstop = -0.05;
  if (notinit(cstart))           cstart = 2;
  if (notinit(cstep))             cstep = 2;
  if (notinit(cstop))             cstop = 3;
  if (notinit(iaii))               iaii = 2e-12;

  set_run_on_exit(runonexit);                         // set to erase savefile on ^C

  /*
  if (make_rbp) {
      for (npt=nodepnt; npt=foreach (npt, rbp, -1, soma, NULL, &dbpn, NULL); npt=npt->next)  {
          cclamp (ndn(rbp,dbpn,soma), dbp_hstim, 0, 0.05);
          cclamp (ndn(rbp,dbpn,soma), dbp_hstim2, 0.05, 0.15);
	  plot (V, ndn(rbp,dbpn,soma),max=-0.0, min=-0.07); plot_param("rbp",c=dbpn,pl=6,s=1);
      }

  }
  if (make_aii)  {
      for (npt=nodepnt; npt=foreach (npt, aii, -1, soma, NULL, &aiin, NULL); npt=npt->next)  {
          cclamp (ndn(aii,aiin,soma), aii_hstim, 0.2, 1);
          cclamp (ndn(aii,aiin,soma), aii_hstim2, 0, 0.01);
	  plot (V, ndn(aii,aiin,soma),max=-0.0, min=-0.07); plot_param("aii",c=aiin,pl=5,s=1);
      }
      //plot (V, ndn(aii,midaii,lob_append),max=-0.0, min=-0.07); plot_param("aii_lob",c=13,pl=5,s=1);
  }
  */

  /* plot (FA9, fwsynat[midrod+2] max 200 min 0);  */
  /* plot (FB1, fwsynrb[midrod] max 10 min 0); */
  /* plot (FC1, fwsynrb[midrod] max 10 min 0); */



//  midrod = findmid(xrod,0,0);
//  plot_v_nod(ct=xrod,cn=1,n=soma,Vmin=-.050,Vmax =-.025,colr=cyan,"", -1, -1); /* Vrod*/
  plot_synrate_out(xrod,1,rmin=0,rmax=400,colr=magenta);	              /* rate out */
//  plot_syncond(xrod,1,rmin=0,rmax=1e-9,colr=red,89,"",1);	              /* conductance out */
//  plot_v_nod(ct=xrod,cn=midrod,n=soma,Vmin=-.060,Vmax =-.045,colr=cyan,"", -1, -1); /* Vrod*/
//  plot_synrate_out(ct=xrod,cn=midrod,rmin=0,rmax=400,colr=magenta);	              /* rate out */

  endexp=prestimdur+stimdur+poststimdur;

  sprintf (savefile,"aii_dbp%06d",getpid());       // add pid to file name


  plot (V, ndn(xrod,1,soma),max=-0.03, min=-0.045); plot_param("rod soma",c=white,pl=170,1);
//  plot (V, ndn(xrod,1,axtrm),max=-0.04, min=-0.07); plot_param("rod axon",c=magenta,pl=70,1);

  plot (V, ndn(rbp,1,soma),max=-0.04, min=-0.06); plot_param("rbp soma",c=cyan,pl=60,1);
//  plot (V, ndn(rbp,1,120),max=-0.04, min=-0.06); plot_param("120",c=white,pl=60,1);
//  plot (V, ndn(rbp,1,130),max=-0.03, min=-0.05); plot_param("130",c=blue,pl=60,1);

  /* find a specific node/dendrite on rbp and record from it */
//  rbp_recnode = 170;
//  sprintf (rbp_n, "rbp_%d\n", rbp_recnode);
//  plot (V,ndn(rbp,1,rbp_recnode),max=-0.02, min=-0.06); plot_param(rbp_n,c=white,pl=50,1);

//  plot_syncond(findsynloc(rbp,0,0), cmax=200e-12, cmin=0, red, 45,"",1);
//  plot_syncond(findsynloc(rbp,0,-40), cmax=200e-12, cmin=0, red, 40,"",1);

//  plot_syncondp(findsynloc(rbp,1,a17,1,  0,40,0.0), cmax=100e-12, cmin=0, red, 50,"",1);
//  plot_syncondp(findsynloc(rbp,1,a17,1,-10,40,0.0), cmax=100e-12, cmin=0, red, 35,"",1);
//  plot_syncondp(findsynloc(rbp,1,aii,1,  0,-40,0.0), cmax=100e-12, cmin=0, red, 25,"",1);
//  plot (V, ndn(rbp,1,axtrm),max=-0.02, min=-0.06); plot_param("rbp axon term",c=white,pl=50,1);

//  plot_syncond(findsynloc(rbp,1,0), cmax=100e-12, cmin=0, white, 30, "",1);
//  plot_syncond(findsynloc(a17,1,0), cmax=10e-12, cmin=0, white, 30, "",1);
//  plot_syncondp(findsynloc(a17,1,rbp,1,0,0,-0.06), cmax=100e-12, cmin=0, white, 30, "",1);


//  plot (V, ndn(aii,1,soma),max=-0.04, min=-0.05); plot_param("aii soma",c=brown,pl=25,1);
//  plot (V, ndn(a17,1,soma),max=-0.04, min=-0.05); plot_param("a17 soma",c=magenta,pl=20,1);

//  plot_syncond(findsynloc(a17,0,40), cmax=10e-12, cmin=0, red, 10,"",1);

  /* find a specific node/dendrite on aii and record from it */
//  aii_recnode = findnodlocr(aii,1,-20,10);

//  aii_recnode = 401;
//  sprintf (aii_n,"aii_%d\n",aii_recnode);
//  plot (V, ndn(aii,1,aii_recnode),max=-0.04, min=-0.05); plot_param(aii_n,c=magenta,pl=15,1);

//  a17_recnode1 = 178;
//  sprintf (a17_n1,"a17_%d\n",a17_recnode1);
//  plot (V, ndn(a17,1,a17_recnode1),max=-0.04, min=-0.05); plot_param(a17_n1,c=blue,pl=10,1);

//  a17_recnode2 = 234;
//  sprintf (a17_n2,"a17_%d\n",a17_recnode2);
//  plot (V, ndn(a17,1,a17_recnode2),max=-0.04, min=-0.05); plot_param(a17_n2,c=blue,pl=0,1);

  //  plot (V, ndn(aii,1,lob_append),max=-0.045, min=-0.05); plot_param("aii lob append",c=blue,pl=0,1);

//  plot (V, ndn(rbp,1,2),max=-0.02, min=-0.06); plot_param("rbp",c=magenta,pl=30,1);
//  plot (FA9,findsynlocr(rbp,1,0,0)->elnum,max=100, min=0); plot_param("FA9",c=blue,pl=24,1);
//  plot (V, ndn(gca, 1,soma),max=-0.04, min=-0.06); plot_param("gca",c=blue,pl=2,1);
//  plot_ca_syn(findsynlocr(rbp,1,0,0), max=2e-6, min=0, green, "",30,1);

  setxmin = 0;
  predur = 0.1;
  // predur = 0;
  simtime = 0-predur; /* start simulation here for calibration before the expt */

  stim_backgr(bg_inten);
  step (predur);

  vhold = -0.035;

//     cclamp (ndn(rbp,1,soma), ipulse=7e-12, simtime, 10); /* run this for 10 sec which is much larger than the total expt duration */
//     vclamp (ndn(rbp,1,soma), vhold, 0, 10); /* run this for 10 sec which is much larger than the total expt duration */

//     step (predur);

//  if (istart < istop) sign = 1;
//  else                sign = -1;

  if (vstart < vstop) sign = 1;
  else                sign = -1;

  savemodel (savefile);

//  for (i=0,ipulse=istart; (ipulse*sign)<=(istop*sign); i++, ipulse += istep*sign) {
//
//     simtime = 0;
//     step (prestimdur);
//     cclamp (ndn(rbp,1,soma), ipulse, prestimdur, stimdur);
//     step (stimdur);
//     step (poststimdur);
//     restoremodel (savefile);
//  }


//  for (v=0,cvolt=vstart; (cvolt*sign)<=(vstop*sign); v++, cvolt += vstep*sign) {
//     simtime = 0;
//      vclamp (ndn(rbp,1,soma), vhold, simtime, prestimdur);
//     step (prestimdur);
//     vclamp (ndn(rbp,1,soma), cvolt, prestimdur, stimdur);
//     step (stimdur);
//      vclamp (ndn(rbp,1,soma), vhold, prestimdur+stimdur, poststimdur);
//     step (poststimdur);
//     restoremodel (savefile);
//  }

/*
  for (v=0,cvolt=vstart; (cvolt*sign)<=(vstop*sign); v++, cvolt += vstep*sign) {
//  for (r=1; r<= n_rods; r++) {
     simtime = 0;
//    fprintf (stderr,"# rod number %d\n", r);
//      vclamp (ndn(xrod,1,soma), vhold, simtime, prestimdur);
     step (prestimdur);
      vclamp (ndn(xrod,1,soma), cvolt, prestimdur, stimdur);
      vclamp (ndn(xrod,2,soma), cvolt, prestimdur, stimdur);
      vclamp (ndn(xrod,3,soma), cvolt, prestimdur, stimdur);
      vclamp (ndn(xrod,4,soma), cvolt, prestimdur, stimdur);
     step (stimdur);
//      vclamp (ndn(xrod,1,soma), vhold, prestimdur+stimdur, poststimdur);
     step (poststimdur);
//  }  
     restoremodel (savefile);
  }
*/

  for (cont=cstart; cont<=cstop; cont+=cstep){
    simtime = 0;
    step (prestimdur);
    stim_spot(sdia, 0, 0, minten*cont, prestimdur, stimdur);
    //stim_spot(sdia, 0, 0, minten, prestimdur, stimdur);
    step(stimdur);
    step(poststimdur);
    restoremodel (savefile);
  }

  unlink (savefile);
  savefile[0] = 0;

//  for (t=0; t<exptdur; t+= dtrial){
//     double start, dur;
//    stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=dstim);
//    step(dtrial);
//  }

}
