/* Experiment rod_cone_hz_rbp */

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
int rodstimchan;
int conestimchan;
int hbarr;
int hbatarr;

double temp_freq;
double dstim;
double sdia;
double predur;
double stimtime;
double prestimdur;
double stimdur;
double poststimdur;
double minten;
double csinten;
double rsinten;
double cbg_inten;
double rbg_inten;
double scontrast;
double setploti;
double dvrev;
double dvst;
double dvsr;
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
double conerm;
double rodrm;
double rbrm;
double spotx;
double spoty;
double hbaxdia;
double hbaxlen;

double roddens;
double conedens;

double g_cone_cone;
double g_cone_rod;
double g_cone_hb;
double g_rod_hb;
double g_rod_rbp;

double disp_stim_min;
double disp_stim_max;
double disp_stim_incr;

void vinit(void);

int rec_ct;
int rec_cn;
node **alt_hb_somas;

char savefile[30] = {0};


void defparams(void)
{
  setptr("rodarr",     &rodarr);
  setptr("rbparr",     &rbparr);
  setptr("aiiarr",     &aiiarr);
  setptr("a17arr",     &a17arr);

  setptr("rec_ct",     &rec_ct);
  setptr("rec_cn",     &rec_cn);
  setptr("temp_freq",  &temp_freq);
  setptr("ntrials",    &ntrials);
  setptr("dstim",      &dstim);
  setptr("sdia",       &sdia);
  setptr("predur",     &predur);
  setptr("stimtime",   &stimtime);
  setptr("prestimdur", &prestimdur);
  setptr("stimdur",    &stimdur);
  setptr("poststimdur",&poststimdur);
  setptr("minten",     &minten);
  setptr("csinten",    &csinten);
  setptr("rsinten",    &rsinten);
  setptr("rbg_inten",  &rbg_inten);
  setptr("cbg_inten",  &cbg_inten);
  setptr("scontrast",  &scontrast);
  setptr("setploti",   &setploti);
  setptr("dvst",       &dvst);
  setptr("dvsr",       &dvsr);
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
  setptr("conerm",     &conerm);
  setptr("rodrm",      &rodrm);
  setptr("rbrm",       &rbrm);
  setptr("roddens",    &roddens);
  setptr("conedens",   &conedens);
  setptr("g_cone_cone",&g_cone_cone);
  setptr("g_cone_rod", &g_cone_rod);
  setptr("g_cone_hb",  &g_cone_hb);
  setptr("g_rod_hb",   &g_rod_hb);
  setptr("g_rod_rbp",  &g_rod_rbp);
  setptr("spotx",      &spotx);
  setptr("spoty",      &spoty);
  setptr("hbaxdia",    &hbaxdia);
  setptr("hbaxlen",    &hbaxlen);

  setptr("disp_stim_min",  &disp_stim_min);
  setptr("disp_stim_max",  &disp_stim_max);
  setptr("disp_stim_incr", &disp_stim_incr);

  setptr("hbarr",	&hbarr);
  setptr("hbatarr",	&hbatarr);
 
  chanparamsfile    = "chanparams_rod_cone_hz_rbp"; 
  nvalfile	    = "nval_rbp_aii_a17.n";
  cone_densfile	    = "dens_cone_fb.n";
  rod_densfile	    = "dens_rod_fb.n";
  rbp_densfile	    = "dens_rbp.n";
  aii_densfile	    = "dens_aii.n";
  a17_densfile	    = "dens_a17.n";
  gca_densfile	    = "dens_gca_aii_rbp.n";

  if (notinit(rbp_file)) rbp_file = (char *)"morph_rbp_yoshi_cell14_normalized"; /* real rbp morphology from Tsukamoto & Omi, 2013 */
  if (notinit(aii_file)) aii_file = (char *)"morph_aii_yoshi_cell2_normalized";  /* real aii morphology from Tsukamoto & Omi, 2013 */

  if (notinit(rbp_thetax)) rbp_thetax=55;
  if (notinit(rbp_thetay)) rbp_thetay=14;
  if (notinit(aii_thetax)) aii_thetax=48;
  if (notinit(aii_thetay)) aii_thetay=-150;	/* optimal x,y rotations for rbp and aii cells from Yoshi - determined by trial and error */

  if (notinit(rbp_nscale)) rbp_nscale=-3.1;
  if (notinit(rod_nscale)) rod_nscale=-3.1;
  if (notinit(a17_nscale)) a17_nscale=-3.1;
  if (notinit(aii_nscale)) aii_nscale=-3.1;	/* node display params */

  if (notinit(hbaxdia)) hbaxdia = 0.4;		/* diameter of axon connecting hb soma to its axon field */
  if (notinit(hbaxlen)) hbaxlen = -50;		/* length of axon connecting hb soma to its axon field */

  // arrsiz = 100;

  make_rod_rbp  = 1;
  make_rbp_aii  = 1;
  make_aii_aii  = 1;
  make_rbp_a17  = 1;
  make_a17_rbp  = 1;

  make_cone_cone  = 1;
  make_cone_rod  = 1;
   
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
	 sprintf (tbuf,"%d %d %d %d", 
		spnt->node1b,spnt->node1c, spnt->node2b, spnt->node2c);  /* print pre and postsyn cell number */
	 tlen = strlen(tbuf);
	 gmove (length/2.0 -tlen*0.3, 0);
	 gcwidth (0.5);
	 gtext (tbuf);
}

/*------------------------------------------------------*/

void rayphot_draw2 ( int ctype, double dia, double xr, double yr, double zr,
		                double tx1, double ty1, double tz1, int color)
{ }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void  phot_draw2 (int ctype, int pigm, int color, double dscale, double dia, double dist, int hide)
{ }

/*------------------------------------------------------*/

void setparams(void)

{
  make_rods  = 1;        /* make rods, rbc, aii, a17 */
  make_cones = 1;
  make_ha    = 0;
  make_hb    = 1;
  make_hbat  = 1;
  make_dbp2  = 0;
  make_rbp   = 1;
  make_aii   = 0;
  make_a17   = 0;
  make_gca   = 0;
  make_gcb   = 0;

  set_synapse_dr (syn_draw2);
  set_phot_dr (phot_draw2);
  set_rayphot_dr (rayphot_draw2);

  vna =  0.04;
  vcl = -0.060;
  drm = 30e3;
  dicafrac = 0;				// set ca pump current = 0

  _CA_T = C_CA7;

 dnoise = 0;
 // pnoise = 0;

  if(notinit(rec_ct)) rec_ct = aii;
  if(notinit(rec_cn)) rec_cn = 1;
//  if (notinit(arrsiz)) arrsiz = 50;

  if (notinit(dvsr))   dvsr = -0.050;   // rod start voltage
  if (notinit(dvst))   dvst = -0.065;
  if (notinit(dvrev)) dvrev = -0.065;
  if (notinit(cbplam)) cbplam = 0.005;

  if (!notinit(roddens)) setn(xrod,DENS,roddens);
  if (!notinit(conedens)) setn(xcone,DENS,conedens);

  if (notinit(rod_stimchan))   rod_stimchan = 1;
  if (notinit(cone_stimchan)) cone_stimchan = 2;

  				// colors defined in ccols2.txt
  setn(xrod,NCMAP,9);
  setn(xrod,NCOLOR,117);	// 117 magenta4 
  
  setn(xcone,NCMAP,9);
  setn(xcone,NCOLOR,49);	// 49 OliveDrab
  
  setn(hb,NCMAP,9);
  // setn(hb,NCOLOR,37);
  setn(hb,NCOLOR,48);		// 48 ForestGreen
 
  setn(hbat,NCMAP,9);
  setn(hbat,NCOLOR,14);		// 14 MidnightBlue

  setn(rbp,NCMAP,9);
  setn(rbp,NCOLOR,76);		// 76 DarkOrange

  // setsv(xrod,SCOND,2,20e-10);
  // setsv(xcone,SCOND,6,50e-10);
 
  if (notinit(rodrm))   rodrm = 10e3;	// in dens_rod_fb.n
  if (notinit(conerm)) conerm = 10e3;	// in dens_cone_fb.n
  if (notinit(rbrm))     rbrm = 10e3;	// in dens_rbp.n

  if (!notinit(g_cone_cone)) setsv(xcone,SCOND,7,g_cone_cone);
  else                       setsv(xcone,SCOND,7,28*100e-12);  // 27 connexins of 100e-12 S

  if (!notinit(g_cone_rod)) setsv(xcone,SCOND,8,g_cone_rod);
  else                      setsv(xcone,SCOND,8,40*100e-12);   // 40 connexins of 100e-12 S

  if (!notinit(g_cone_hb))   setsv(xcone,SCOND,6,g_cone_hb);
  else                       setsv(xcone,SCOND,6,50e-10); 

  if (!notinit(g_rod_hb))    setsv(xrod,SCOND,4,g_rod_hb);    // rod -> hbat 
  else                       setsv(xrod,SCOND,4,10e-10); 

  if (!notinit(g_rod_rbp))   setsv(xrod,SCOND,1,g_rod_rbp);    // rod -> rbp 
  else                       setsv(xrod,SCOND,1,1e-10); 

  // setsv(hbat,SCOND,1,1e-10);
  if (notinit(dvfbg)) dvfbg = 0.2;

  // rod_timec = 0.3;
  // rod_maxcond = 0e-12;

  dscavg = 7e4;
  _CA_L = _CA1;                 /* set type of L-type calcium channel for dens_rod_fb.n */
  // _CA_T = _CA7;                 /* set type of T-type calcium channel for dens_rod_fb.n */

   if (notinit(disp_stim_incr))  disp_stim_incr = 0.002;         /* time step for stimulus display */
   if (notinit(disp_stim_max))   disp_stim_max = 1e5;            /* max level for stimulus display */
   if (notinit(disp_stim_min))   disp_stim_min = 1e3;            /* min level for stimulus display */

   dscavg = 1e5;

#define HBARR 10
#define HBATARR 20

   if (notinit(hbarr)) hbarr = 0;

   if (hbarr>0) {         // array of hb locations, out of the way of rbp
      n_hb = hbarr;
      if (n_hb >=HBARR) n_hb = HBARR-1;

      hbxarr  = (double *)emalloc(HBARR*sizeof(double));
      hbyarr  = (double *)emalloc(HBARR*sizeof(double));
      hbtharr = (double *)emalloc(HBARR*sizeof(double));  // rotation array
   
     if (hbarr==1) { 
       //hbxarr[0] = 30;
       hbxarr[0] = arrsiz * 0.4;
       hbyarr[0] = 0;
      hbtharr[0] = 0;
     }
     else { 
      hbxarr[0] = 60;
      hbyarr[0] = 20;
     hbtharr[0] = -20;

      hbxarr[1] = -60;
      hbyarr[1] = -20;
     hbtharr[1] = 160;

      hbxarr[2] = 20;
      hbyarr[2] = -20;
     hbtharr[2] = 0;
    }
   }

   if (notinit(hbatarr)) hbatarr = 0;

   if (hbatarr>0) {         // array of hb locations, out of the way of rbp
      n_hbat = hbatarr;
      if (n_hbat >=HBATARR) n_hbat = HBATARR-1;

      hbatxarr  = (double *)emalloc(HBATARR*sizeof(double));
      hbatyarr  = (double *)emalloc(HBATARR*sizeof(double));
      hbattharr = (double *)emalloc(HBATARR*sizeof(double));  // rotation array
   
      hbatxarr[0] = -10;
      hbatyarr[0] = 10;
     hbattharr[0] = 0;

      hbatxarr[1] = 10;
      hbatyarr[1] = -10;
     hbattharr[1] = 0;

      hbatxarr[2] = 20;
      hbatyarr[2] = -20;
     hbattharr[2] = 0;

      hbatxarr[3] = -20;
      hbatyarr[3] = 20;
     hbattharr[3] = 0;
   }
   make_hb_hb = 0;
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
   int midrbp, mida17, midaii, middbp ;
   int midrod,  rod2, rod3, rod4, rod5;
   int midcone, cone2, cone3, cone4, cone5;
   int rbp_recnode, aii_recnode, a17_recnode1, a17_recnode2;
   char rbp_n[20], aii_n[20], a17_n1[20], a17_n2[20];
   double s, t, start, dur, dtrial, exptdur;
   double max,min, rmax, rmin;
   double dscale, plsize;
   double ipulse, cvolt, vhold, cont, sign;
   double cmax, cmin, Vmin, Vmax;
   double starttime, disp_end;

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
  if (notinit(stimtime))       stimtime = 0.10;  /* stimulus time */
  if (notinit(dstim))             dstim = 0.05;	/* stimulus duration */
  if (notinit(sdia))               sdia = 30;	/* spot diameter */
  if (notinit(prestimdur))   prestimdur = 0.05;	/* before stimulus */
  if (notinit(stimdur))         stimdur = 0.2;	/* stimulus duration */
  if (notinit(poststimdur)) poststimdur = 0.2;	/* after stimulus */

  if (notinit(bg_inten))       bg_inten = 5000;
  if (notinit(minten))           minten = bg_inten;	/* background intensity (for rods and cones)*/
  // if (notinit(rsinten))         rsinten = bg_inten;	/* background intensity (for rods )*/
  // if (notinit(csinten))         csinten = bg_inten;	/* background intensity (for cones )*/
  if (notinit(rbg_inten))     rbg_inten = bg_inten;	/* background intensity (for rods)*/
  if (notinit(cbg_inten))     cbg_inten = bg_inten;	/* background intensity (for cones)*/
  if (notinit(scontrast))     scontrast = 1.0;
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
          
  
midrod = findmid(xrod,0,0);
rod2   = findmid(xrod,10,0);
rod3   = findmid(xrod,20,0);
rod4   = findmid(xrod,-20,0);
rod5   = findmid(xrod,-30,0);

midcone = findmid(xcone,0,0);
cone2   = findmid(xcone,10,0);
cone3   = findmid(xcone,20,0);
cone4   = findmid(xcone,30,0);
cone5   = findmid(xcone,40,0);

//  plot_v_nod(ct=xrod,cn=1,n=soma,Vmin=-.050,Vmax =-.025,colr=cyan,"", -1, -1); /* Vrod*/
//  plot_synrate_out(xrod,1,rmin=0,rmax=400,colr=magenta);	              /* rate out */
//  plot_syncond(xrod,1,rmin=0,rmax=1e-9,colr=red,89,"",1);	              /* conductance out */
//  plot_v_nod(ct=xrod,cn=midrod,n=soma,Vmin=-.060,Vmax =-.045,colr=cyan,"", -1, -1); /* Vrod*/
//  plot_synrate_out(ct=xrod,cn=midrod,rmin=0,rmax=400,colr=magenta);	              /* rate out */

  endexp=prestimdur+stimdur+poststimdur;

  rodstimchan  = 1;
  conestimchan = 2;

  setxmin = 0;
  if (notinit(predur)) predur = 0.3;
  simtime = 0-predur; /* start simulation here for calibration before the expt */

  if (notinit(spotx)) spotx = 0;
  if (notinit(spoty)) spoty = 0;

  stim_backgr(rbg_inten,rodstimchan);
  stim_backgr(cbg_inten,conestimchan);
  if (notinit(rsinten))  rsinten = rbg_inten * scontrast;  /* background intensity (for rods)*/
  if (notinit(csinten))  csinten = cbg_inten * scontrast;  /* background intensity (for rods)*/
  stim_spotc(sdia, spotx, spoty, csinten, 0.25, 0.010,conestimchan);
  stim_spotc(sdia, spotx, spoty, rsinten, 0.05, 0.010,rodstimchan);

  if (disp & DSTIM) {
         double t;
     if (makestim) return;
     if (notinit(dispsize)) {
         maxsize = max(xarrsiz,yarrsiz);
         dispsize=maxsize*(1+disp_margin);
     }
     display_size(dispsize);
     if (notinit(disp_calib_len)) disp_calib_len = set_int_val(dispsize/7);
     if (disp_calib_len > 0) drcalib(0.97,0.05,disp_calib_len,dispsize,white);

     disp_end = stimtime+stimdur+poststimdur;
     for (starttime=simtime,t=0; t<disp_end; starttime = t, t+= disp_stim_incr) {
        display_stim(starttime, t, dscale=2, disp_stim_max, disp_stim_min);
        // display_stim(starttime, t, dscale=4, -0.025, -0.035); 
	//         // display_stim(t, dscale=4, -0.035, -0.045); 
	simwait(0.10);
	if (disp & DMOVIE) display_page();
     }
     return;
  }

  step (predur);

  sprintf (savefile,"rod_cone_aii_dbp%06d",getpid());       // add pid to file name

  r = 10;

  plot_v_nod(xrod,midrod,soma, min= -0.06, max= -0.03,white,"",pl=88,0.5);
  plot_v_nod(xrod,midrod,axtrm,min= -0.06, max= -0.03,magenta,"",pl=88,0.5);
  plot_v_nod(xrod,rod3,    axtrm,min= -0.06, max= -0.03,cyan,"",pl=88,0.5);
  plot_v_nod(xrod,rod5,    axtrm,min= -0.06, max= -0.03,brown,"",pl=88,0.5);
  
  plot_ca_nod(xrod,midrod,axtrm, 2e-6,green,    "",  pl=85, 0.5);
  plot_ca_nod(xrod,rod3,  axtrm, 2e-6,brown,    "",  pl=85, 0.5);

  plot_synrate(findsynloc(xrod,midrod,hb,-1,0,0),rmin=0,rmax=200, blue,    80,"",0.5);
  plot_synrate(findsynloc(xrod,rod3,  hb,-1,0,0),rmin=0,rmax=200, magenta, 80,"",0.5);
  plot_synrate(findsynloc(xrod,rod5,  hb,-1,0,0),rmin=0,rmax=200, green,   80,"",0.5);

  // plot_v_nod(hbat,1,1,       min= -0.07, max= -0.03,green,"hbat axon",pl=70,0.5);
  plot_v_nod(alt_hb_somas[0],min= -0.06, max= -0.04,blue,"hbat_axon",pl=55,0.5);

  /* - - - - - - - - - - - - - - */


  plot_v_nod(xcone,midcone,soma, min= -0.06, max= -0.03,white,  "",pl=68,0.5);
  plot_v_nod(xcone,midcone,axtrm,min= -0.06, max= -0.03,magenta,"",pl=68,0.5);
  plot_v_nod(xcone,cone3,  axtrm,min= -0.06, max= -0.03,cyan,"",pl=68,0.5);
  plot_v_nod(xcone,cone5,  axtrm,min= -0.06, max= -0.03,brown,"",pl=68,0.5);

  plot_ca_nod(xcone,midcone,axtrm, 5e-6,green,     "",  pl=65, 0.5);
  plot_ca_nod(xcone,cone3,  axtrm, 5e-6,brown,      "",  pl=65, 0.5);

  plot_synrate(findsynloc(xcone,midcone,hb,-1,0,0), rmin=0,rmax=200, blue,  60,"",0.5);
  plot_synrate(findsynloc(xcone,cone3,hb,-1,0,0),rmin=0,rmax=200, magenta,  60,"",0.5);
  plot_synrate(findsynloc(xcone,cone5,hb,-1,0,0),rmin=0,rmax=200, green,    60,"",0.5);

  plot_v_nod(hb,1,1,min= -0.06, max= -0.04,green,"hb",pl=50,0.5);
    
  /* - - - - - - - - - - - - - - */

  plot_v_nod(rbp,1,soma, min= -0.07, max= -0.04,cyan,"rbp axon",pl=20,0.5);

//  plot (V, ndn(rbp,1,120),max=-0.04, min=-0.06); plot_param("120",c=white,pl=60,1);
//  plot (V, ndn(rbp,1,130),max=-0.03, min=-0.05); plot_param("130",c=blue,pl=60,1);

  /* find a specific node/dendrite on rbp and record from it */

//  rbp_recnode = 170;
//  sprintf (rbp_n, "rbp_%d\n", rbp_recnode);
//  plot (V,ndn(rbp,1,rbp_recnode),max=-0.02, min=-0.06); plot_param(rbp_n,c=white,pl=50,1);

//  plot_syncond(findsynloc(rbp,0,0), cmax=200e-12, cmin=0, red, 45,"",1);
//  plot_syncond(findsynloc(rbp,0,-40), cmax=200e-12, cmin=0, red, 40,"",1);

//  plot_syncond(findsynloc(rbp,1,0), cmax=100e-12, cmin=0, white, 30, "",1);


//  plot (V, ndn(rbp,1,2),max=-0.02, min=-0.06); plot_param("rbp",c=magenta,pl=30,1);
//  plot (FA9,findsynlocr(rbp,1,0,0)->elnum,max=100, min=0); plot_param("FA9",c=blue,pl=24,1);
//  plot (V, ndn(gca, 1,soma),max=-0.04, min=-0.06); plot_param("gca",c=blue,pl=2,1);
//  plot_ca_syn(findsynlocr(rbp,1,0,0), max=2e-6, min=0, green, "",30,1);

  // for (r=1; r <= nrods; r++) {
  //    cclamp (ndn(xrod,r,soma), 0e-12, simtime, 1);
  // }
  // vclamp (ndn(hbat,1,2), -0.040, 0.05, 0.02);
  // vclamp (ndn(hbat,1,2), -0.050, 0.07, 0.02);
  
  // for (r=1; r <= 20; r++) {
  //   vclamp (nde(xrod,r,soma), -0.030, 0.05, 0.02);
  // }

  /*
  for (c=1; c <= 20; c++) {
    // vclamp (nde(xcone,c,soma), -0.080, 0.05, 0.02);
    // vclamp (nde(xcone,c,soma), -0.060, 0.01, 1);
  } /* */

  // vclamp (nde(xrod,1,soma), -0.040, 0.01, 0.02);

  // vclamp (nde(hb,1,soma), -0.040, 0.01, 0.02);

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
    // stim_spot(sdia, 0, 0, minten*cont, prestimdur, stimdur);
    //  stim_spot(sdia, 0, 0, minten, prestimdur, stimdur);
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
