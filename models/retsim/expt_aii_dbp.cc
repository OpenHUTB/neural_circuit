/* Experiment aii_dbp */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"

//#include "onplot_dsgc_movie.cc"

int ntrials;
int aiiarr;
int dbp1arr;

double temp_freq;
double dstim;
double sdia;
double stimtime;
double prestimdur;
double stimdur;
double poststimdur;
double minten;
double scontrast;
double dvrev;
double dvst;
double cbplam;
double istart;
double istop;
double istep;
double iaii;
double drma;

void vinit(void);

int rec_ct;
int rec_cn;

char savefile[30] = {0};

void defparams(void)
{
  setptr("aiiarr",     &aiiarr);
  setptr("dbp1arr",    &dbp1arr);

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
  setptr("dvst",       &dvst);
  setptr("dvrev",      &dvrev);
  setptr("cbplam",     &cbplam);
  setptr("istart",     &istart);
  setptr("istop",      &istop);
  setptr("istep",      &istep);
  setptr("iaii",       &iaii);
  setptr("drma",       &drma);

  nvalfile	    = "nval_aii_dbp.n";
  aii_densfile	    = "dens_aii_dbp.n";
  dbp1_densfile	    = "dens_aii_dbp1.n";
  gca_densfile	    = "dens_gca_aii_dbp1.n";

}

void setparams(void)

{
  make_rods  = 0;        /* make rods, rbc, aii */
  make_cones = 0;
  make_ha    = 0;
  make_hb    = 0;
  make_hbat  = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_rbp   = 0;
  make_aii   = 1;
  make_gca   = 1;
  make_gcb   = 0;

  make_dbp1_aii  = 1;
  make_aii_aii   = 1;


  DENDD     = R_1;
  DEND_DIST = R_1;
  DEND      = R_2;      /* definitions of regions for dens_ file */
  DENDP     = R_3;
  DEND_PROX = R_3;
  SOMA      = R_4;
  HCK       = R_5;
  HILLOCK   = R_5;
  AXONT     = R_6;
  AXON_THIN = R_6;
  AXON      = R_7;
  AXONP     = R_7;
  AXON_PROX = R_7;
  AXOND     = R_8;
  AXON_DIST = R_8;
  VARIC     = R_9;
  VARICOS   = R_9;

  vna =  0.04;
  vcl = -0.060;
  drm = 30e3;
  if (notinit(drma)) drma = 30e3;
  dicafrac = 0;				// set ca pump current = 0

  _CA_T = C_CA7;

  if(notinit(rec_ct)) rec_ct = aii;
  if(notinit(rec_cn)) rec_cn = 1;
//  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 10;

  if (notinit(dvst))   dvst = -0.065;
  if (notinit(dvrev)) dvrev = -0.065;
  if (notinit(cbplam)) cbplam = 0.005;

  //if (notinit(n_dbp1)) n_dbp1 = 1;
  if (notinit(aiiarr)) aiiarr = 1;
  if (notinit(dbp1arr)) dbp1arr = 1;

#define DBP1ARR 30
  if (!notinit(dbp1arr)) {
    dbp1xarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1yarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1tharr = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1narr  = (int *) emalloc(DBP1ARR*sizeof(int));
    if (dbp1arr==1) {
      dbp1xarr[0] =  10;    dbp1yarr[0] = -5;    dbp1tharr[0] = 0;    dbp1narr[0] = 1;
      n_dbp1 = 1;
    }
    else if (aiiarr==2) {
      dbp1xarr[0] =  10;    dbp1yarr[0] = 0;    dbp1tharr[0] = 0;    dbp1narr[0] = 1;
      dbp1xarr[1] = -10;    dbp1yarr[1] = 0;    dbp1tharr[1] = 0;    dbp1narr[1] = 2;
      n_dbp1 = 2;
    }
  }
#define AIIARR 30
  if (!notinit(aiiarr)) {
    aiixarr  = (double *)emalloc(AIIARR*sizeof(double));
    aiiyarr  = (double *)emalloc(AIIARR*sizeof(double));
    aiitharr = (double *)emalloc(AIIARR*sizeof(double));
    aiinarr  = (int *) emalloc(AIIARR*sizeof(int));
    if (aiiarr==1) {
      aiixarr[0] =  -10;      aiiyarr[0] = 0;       aiitharr[0] = 0;    aiinarr[0] = 1;
      n_aii = 1;
    }
    else if (aiiarr==2) {
      aiixarr[0] =  10;      aiiyarr[0] = 0;       aiitharr[0] = 0;    aiinarr[0] = 1;
      aiixarr[1] = -10;      aiiyarr[1] = 0;       aiitharr[1] = 0;    aiinarr[1] = 2;
      n_aii = 2;
    }
  }
}

/*------------------------------------------------*/

void addcells(void)

// Make a ganglion cell compartment postsynaptic to the bipolar cell.
// To use, set make_gca = 0 above

{      
     node *npt;
     int dbpn, gcn, gcnod;
     node *dbpnod;
     sphere *s;
     cable *c;

  /*
  if (make_dbp1) {
      for (npt=nodepnt; npt=foreach (npt, dbp1, -1, axtrm, NULL, &dbpn, NULL); npt=npt->next)  {
	  dbpnod = ndn(dbp1,dbpn,axtrm);
          loc (nd(gca,gcn=1,soma), dbpnod->xloc+20, dbpnod->yloc, dbpnod->zloc - 20);
          loc (nd(gca,gcn=1,gcnod=10), dbpnod->xloc, dbpnod->yloc, dbpnod->zloc - 2);
          s = make_sphere (ndn(gca,gcn,soma), 5, 20e3);
          c = make_cable (ndn(gca,gcn,soma), ndn(gca,gcn,gcnod), 0.5);
	  c->elabl = "dend_med";
	  ngca = 1;
	  make_gca = 1;
	  setn (gca,MAKE,make_gca);
	  setn (gca,NMADE,ngca);

	  //make_synapse (ndn(dbp1,dbpn,axtrm), ndn(gca,gcn,gcnod));
	  //setcelconn (dbp1,1,gca,1);  // remember the synaptic connection
	  //connected  = ncel_out(dbp1,1,gca);
	  //fprintf (stderr,"addcells: connected %d\n",connected);
      }
   }
   */

}

/*------------------------------------------------------*/

void runonexit (void)
{
   if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------*/

void runexpt(void)

{
   int c, i, middbp, pl, midrod, plnum;
   int aiin, dbpn, lob_append=1;
   int midaii, midaii2, midaii3;
   double s, t, start, dur, dtrial;
   double max,min;
   double plsize,predur;
   double ipulse, sign;
   double cmax, cmin;

   node *npt;

  if (notinit(ntrials))         ntrials = 1;	/* number of trials */
  if (notinit(dstim))             dstim = .01;	/* stimulus duration */
  if (notinit(sdia))               sdia = 300;	/* spot diameter */
  if (notinit(prestimdur))   prestimdur = .01;	/* before stimulus */
  if (notinit(stimdur))         stimdur = .02;	/* stimulus duration */
  if (notinit(poststimdur)) poststimdur = .1;	/* after stimulus */

  if (notinit(minten))           minten = bg_inten;	/* background intensity (for make cone)*/
  if (notinit(scontrast))     scontrast = 0.99;	/* intensity increment */
  if (notinit(istart))           istart = 0e-12;
  if (notinit(istep))             istep = 0.5e-12;
  if (notinit(istop))             istop = 5e-12;
  if (notinit(iaii))               iaii = 2e-12;

  set_run_on_exit(runonexit);                         // set to erase savefile on ^C

  /*
  if (make_dbp1) {
      for (npt=nodepnt; npt=foreach (npt, dbp1, -1, soma, NULL, &dbpn, NULL); npt=npt->next)  {
          cclamp (ndn(dbp1,dbpn,soma), dbp_hstim, 0, 0.05);
          cclamp (ndn(dbp1,dbpn,soma), dbp_hstim2, 0.05, 0.15);
	  plot (V, ndn(dbp1,dbpn,soma),max=-0.0, min=-0.07); plot_param("dbp1",c=dbpn,pl=6,s=1);
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

  //stim_backgr(minten);
  //
  endexp=prestimdur+stimdur+poststimdur;

  sprintf (savefile,"aii_dbp%06d",getpid());       // add pid to file name


  plot (V, ndn(aii, 1,soma),max=-0.0, min=-0.07); plot_param("aii",c=brown,pl=40,1);
  plot (V, ndn(dbp1,1,soma),max=-0.02, min=-0.06); plot_param("dbp",c=cyan,pl=30,1);
  // plot (V, ndn(dbp1,1,2),max=-0.02, min=-0.06); plot_param("dbp",c=magenta,pl=30,1);
  plot_ca_syn(findsynlocr(dbp1,1,0,0),   max=2e-6, green, "", 25, 1);
  //plot (FA9,findsynlocr(dbp1,1,0,0)->elnum,max=100, min=0); plot_param("FA9",c=blue,pl=24,1);
  plot_syncond(findsynloc(dbp1,0,20),   cmin=0,cmax=100e-12, red, 20,"",1);
  plot (V, ndn(gca, 1,soma),max=-0.04, min=-0.06); plot_param("gca",c=blue,pl=2,1);
 

  setxmin = 0;
  predur = 0.2;
  simtime = 0-predur;

  cclamp (ndn(aii,1,soma), iaii, simtime, 1);

  step (predur);

  if (istart < istop) sign = 1;
  else                sign = -1;

  savemodel (savefile);

  for (i=0,ipulse=istart; (ipulse*sign)<=(istop*sign); i++, ipulse += istep) {

     simtime = 0;
     step (prestimdur);
     cclamp (ndn(dbp1,1,soma), ipulse, prestimdur, stimdur);
     step (stimdur);
     step (poststimdur);
     restoremodel (savefile);
  }
  unlink (savefile);
  savefile[0] = 0;

}
