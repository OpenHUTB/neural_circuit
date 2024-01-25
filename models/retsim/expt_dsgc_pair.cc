/* Experiment dsgc_pair */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int n_dsgc;

int rec_ct;
int rec_cn;
int dsgcarr;
int excl_cn;

double minten;
double sinten;
double stimtime;
double stimdur;
double endwait;
double predur;
double currinjd;
double currinja;
double sdia;
double sloc1;
double sloc2;
double iloc;
double set_drm;
double gj_cond;
double dbp1_cond;
double sq_wave;
double fscale;
double theta;
double iroff;
double hroff;
double velocity;
double barwidth;
double barsep;
double dsgc_spacing;
double set_vclamp1;
double set_vclamp2;
double set_endexp;

int qx314;
int spike_pref;
int ninhib;
int light_inhib;
int spot_stim;
int bar_stim;
int move_left;
int ntrials;
int dsgc1, dsgc2;

/*--------------------------------------------------------*/

void defparams(void)
{
  defparams_dsgc_movie();
  defparams_onplot_movie();

  setptr("rec_cn",     &rec_cn);
  setptr("rec_ct",     &rec_ct);
  setptr("dsgcarr",    &dsgcarr);
  setptr("excl_cn",    &excl_cn);

  setptr("minten",     &minten);
  setptr("sinten",     &sinten);
  setptr("stimtime",   &stimtime);
  setptr("stimdur",    &stimdur);
  setptr("endwait",    &endwait);
  setptr("predur",     &predur);
  setptr("currinjd",   &currinjd);
  setptr("currinja",   &currinja);
  setptr("sdia",       &sdia);
  setptr("sloc1",      &sloc1);
  setptr("sloc2",      &sloc2);
  setptr("iloc",       &iloc);
  setptr("set_drm",    &set_drm);
  setptr("gj_cond",    &gj_cond);
  setptr("dbp1_cond",  &dbp1_cond);
  setptr("ninhib",     &ninhib);
  setptr("spike_pref", &spike_pref);
  setptr("sq_wave",    &sq_wave);
  setptr("fscale",     &fscale);
  setptr("spot_stim",  &spot_stim);
  setptr("bar_stim",   &bar_stim);
  setptr("theta",      &theta);
  setptr("iroff",      &iroff);
  setptr("hroff",      &hroff);
  setptr("light_inhib",&light_inhib);
  setptr("velocity",   &velocity);
  setptr("barwidth",   &barwidth);
  setptr("barsep",     &barsep);
  setptr("dsgc_spacing", &dsgc_spacing);
  setptr("set_vclamp1", &set_vclamp1);
  setptr("set_vclamp2", &set_vclamp2);
  setptr("set_endexp", &set_endexp);

  setptr("qx314",      &qx314);
  setptr("move_left",  &move_left);
  setptr("ntrials",    &ntrials);
  if (notinit(dsgc_file)) dsgc_file = "morph_03-17-2010c4";
  nvalfile = "nval_dsgc_pair.n";
}

/*--------------------------------------------------------*/

void setparams(void)

        /*  set up default configuration for sb expts */
        /* cones, cone bipolars, sb, dsgc */
{
    int i,ct;
    double zmax, zmin;

  make_rods  = 0;
  make_cones = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_hbp1  = 1;
  make_hbp2  = 0;
  // make_ams   = 0;
  make_sbac  = 0;
  make_dsgc  = 1;
  if (!notinit(light_inhib) && light_inhib==1) make_ams = 1;   /* make amacrines for inhibition */

  ct = dsgc;
  make_dsgc_dsgc = 1;
 
  onplot_dsgc_movie_init();             /* initialize dsgc movie stuff */
  onplot_movie_init();                  /* initialize onplot_movie stuff */

  if(notinit(rec_ct)) rec_ct = dsgc;    /* type of cell to record from */
  if(notinit(rec_cn)) rec_cn= -1;         /* cellnum to record from */

  if(notinit(dsgc_spacing)) dsgc_spacing=80;         /* dsgc spacing */
  setn(ct,NCOLOR,RCOLOR);

  if (notinit(dend_dia_factor))   dend_dia_factor = 1;
  if (notinit(dendd_dia_factor)) dendd_dia_factor = 0.3;
  if (notinit(dendp_dia_factor)) dendp_dia_factor = 0.5;
  if (notinit(ax_dia_factor))       ax_dia_factor = 0.5;

  if (notinit(dvrev)) dvrev = -0.06;
  if (notinit(dvst)) dvst = dvrev;
  if (notinit(qx314)) qx314 = 0;
    
  if (!notinit(set_drm)) {             /* user set default Rm */
     setn(ct,NRM,set_drm); /* set default Rm */
     drm = set_drm;
  }
  if (!notinit(gj_cond)) setn(ct,SCOND1,gj_cond); /* set gap junction conductance */
  if (gj_cond==0) make_dsgc_dsgc = 0;		  /* if no conductance, make no connections */


  if (!notinit(dbp1_cond)) setn(dbp1,SCOND2,dbp1_cond); /* set bipolar synapse conductance */

  dpcanmda = 0;					  /* turn off calcium flux for NMDA receptors */

  // display_z(zmin=-20, zmax=-40);            /* exclude dsgc Off-arborization layer */

#define GCARR 20

  if (!notinit(dsgcarr)) {         	// arrays of dsgc locations and angles
     gcxarr  = (double *)emalloc(GCARR*sizeof(double));
     gcyarr  = (double *)emalloc(GCARR*sizeof(double));
     gctharr = (double *)emalloc(GCARR*sizeof(double));

     if (dsgcarr==0) {				// aligned in the same direction				
       gcxarr[0] =  0;  gcyarr[0] = 0; gctharr[0] = 55;
       n_dsgc = 1;
     }
     if (dsgcarr==1) {				// aligned in the same direction				
       gcxarr[0] =  0;  gcyarr[0] = 0; gctharr[0] = 55;
       gcxarr[1] =  95; gcyarr[1] = 0; gctharr[1] = 55;
       n_dsgc = 2;
     }
     if (dsgcarr==2 || dsgcarr==102) {		// aligned in the same direction				
       gcxarr[0] =  0;               gcyarr[0] = 0; gctharr[0] = 55;
       gcxarr[1] =  1*dsgc_spacing;  gcyarr[1] = 0; gctharr[1] = 55;
       n_dsgc = 2;
     }
     if (dsgcarr==3) {				// aligned in the same direction				
       gcxarr[0] =  0;   gcyarr[0] = 0; gctharr[0] = 55;
       gcxarr[1] =  50;  gcyarr[1] = 0; gctharr[1] = 55;
       n_dsgc = 2;
     }
     if (dsgcarr==4) {				// aligned in the same direction				
       gcxarr[0] =  0;              gcyarr[0] = 0; gctharr[0] = 55;
       gcxarr[1] =  1*dsgc_spacing; gcyarr[1] = 0; gctharr[1] = 55;
       gcxarr[2] =  2*dsgc_spacing; gcyarr[2] = 0; gctharr[2] = 55;
       gcxarr[3] =  3*dsgc_spacing; gcyarr[3] = 0; gctharr[3] = 55;
       n_dsgc = 4;
     }
     if (dsgcarr==6) {				// aligned in the same direction				
       gcxarr[0] =  0;              gcyarr[0] = 0; gctharr[0] = 55;
       gcxarr[1] =  1*dsgc_spacing; gcyarr[1] = 0; gctharr[1] = 55;
       gcxarr[2] =  2*dsgc_spacing; gcyarr[2] = 0; gctharr[2] = 55;
       gcxarr[3] =  3*dsgc_spacing; gcyarr[3] = 0; gctharr[3] = 55;
       gcxarr[4] =  4*dsgc_spacing; gcyarr[4] = 0; gctharr[4] = 55;
       gcxarr[5] =  5*dsgc_spacing; gcyarr[5] = 0; gctharr[5] = 55;
       n_dsgc = 6;
     }
     switch (dsgcarr) {
       case 0: dsgc1 = 1; dsgc2 = 2; break;
       case 1:
       case 2:
       case 102:
       case 3: dsgc1 = 1; dsgc2 = 2; break;
       case 4: dsgc1 = 2; dsgc2 = 3; break;
       case 6: dsgc1 = 3; dsgc2 = 4; break;
     }
  }
  if (notinit(bar_stim)) bar_stim = 0;
  if (notinit(spike_pref)) spike_pref = 1;
  if (notinit(n_dsgc)) n_dsgc = 2;
  if (notinit(dsgc_densfile))  dsgc_densfile = "dens_dsgc_pair.n";
  if (notinit(dsgc_densfile2)) dsgc_densfile = "dens_dsgc_pairx2.n";
}

/*--------------------------------------------------------*/
void setdens(void)
	/* set density and conductance parameters */

{
    // if (!notinit(kdr_cond))    celdens[dbp1][0][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    // setsv (dbp1,SCOND,1, 25e-10);
    // } 
    // ndens[dsgc][cn=1] = 0;              // set cn 1 to use dsgc_densfile
    ndens[dsgc][dsgc1] = qx314;         // qx314 == 1 -> set dsgc1 to use dsgc_densfile2

}

/*--------------------------------------------------------*/

#define CSIZ 80
#define NTIMEC 10

void runexpt(void) 

{
    int c, n, ct, cn, cbp_cn, pl;
    int stimnod, inod, distnod1,distnod2;
    double t, start, dia, dur, trialdur, dscale, bpsize;
    double ixoff, iyoff, hxoff, hyoff;
    double Vminp, Vmaxp, Imin, Imax;
    double rmin, rmax, fmin, fmax;
    double bar_xmax, bar_xmin, bar_dist;
    double fymin, fymax, fgain, foffs;
    node *npnt;
    elem *epnt;
    synapse *spnt;
    photorec *p;
    chattrib *a;
    nattrib *napnt;
    char gbuf[CSIZ];
    double timec[NTIMEC]={0};
    int ntimec = 2;
    double smin, smax;

  timinc = 5e-6;
  if (!strcmp(dsgc_file,"morph_03-17-2010c4")) {
    cbp_cn    = 68;

    if (dsgcarr==1) {
      distnod1 = 970;
      distnod2 = 1140;
    } else
    if (dsgcarr==2 || dsgcarr==102) {
      distnod1 = 179;
      distnod2 =  179;
    }
    if (dsgcarr==3) {
      distnod1 = 790;
      distnod2 =  790;
    }
 }
  if (notinit(sq_wave)) sq_wave = 0;
  if (notinit(fscale)) fscale = 0.1;

  if (notinit(currinjd)) currinjd = 50e-12;
  if (notinit(stimtime)) stimtime = 0.01;
  if (notinit(stimdur))   stimdur = 0.05;	/* used for non-moving stimuli */
  //if (notinit(endwait))   endwait = 0.005;
  if (notinit(endwait))   endwait = stimdur;
  if (notinit(ntrials))   ntrials = 1;

   if (dsgcarr<4) bpsize = 1;			/* for bp plots below */
   else           bpsize = 1;

   /* Remove excitatory connections from one dsgc */

   if (!notinit(excl_cn)) {
      for(epnt=elempnt; epnt=foreach(epnt,SYNAPSE,dbp1,-1,-1,NULL,&cn,NULL); epnt=epnt->next) {
        if (epnt->node2b==excl_cn) {
 	  spnt = (synapse *)epnt;
	  spnt->maxcond = 0;			/* remove synapse to excluded cell */
        }
      }
      for(epnt=elempnt; epnt=foreach(epnt,SYNAPSE,hbp1,-1,-1,NULL,&cn,NULL); epnt=epnt->next) {
        if (epnt->node2b==excl_cn) {
 	  spnt = (synapse *)epnt;
	  spnt->maxcond = 0;			/* remove synapse to excluded cell */
        }
      }
   }



   smin = -0.05;
   smax = -0.03;

   Imax =  1e-10;
   Imin = -1e-10;

   Vmaxg =  0.04;
   Vming = -0.075;
   Vmaxp = -0.060;
   Vminp = -0.070;
   stimnod = distnod1;
   fgain = 0.01/fscale;
   foffs = -0.068;
   fymin = foffs;
   fymax = foffs+fgain;

   // fprintf (stderr,"fymax %g fymin %g\n",fymax,fymin);

   if (make_movie) {
       setonplot(onplot_movie); /* set movie plot routine */
       if (space_time) {  /* movie */
         plot_v_nod(ct=dsgc,cn=1,soma,Vming,Vmaxg,c=1,"Vsoma1",pl=10,0.35);
         plot_v_nod(ct=dsgc,cn=2,soma,Vming,Vmaxg,c=4,"Vsoma2",pl=10,0.35);
         //plot_v_nod(ct=dsgc,cn,1336,Vming,Vmaxg,c=green,"Vtip1",pl=10,0.35); // for morph_ds1e
         //plot_v_nod(ct=dsgc,cn,582,Vming,Vmaxg,c=red,"Vtip2",pl=10,0.35);
         //plot_na_inact(dsgc, 1, soma, red, "Na[soma]", 12, .35);
	 //
       };
   } else if (spike_pref) {
    timec[0]=0.0002; timec[1]= 0.008;
    if (!notinit(currinja) || !notinit(spot_stim)) Vmaxp = -0.02;
    if (bar_stim>0) Vmaxp = 0.02;
    plot_v_nod(ct=dsgc,cn=1,soma,    Vminp,Vmaxp,red,"", 1, -1);	/* plot V at soma */
    plot_v_nod(ct=dsgc,cn=1,distnod1,Vminp,Vmaxp,brown,"", 1, -1);/* plot V at dendritic tip */
    plot_v_nod(ct=dsgc,cn=2,soma,    Vminp,Vmaxp,blue,"", 2, -1);	/* plot V at soma */
    plot_v_nod(ct=dsgc,cn=2,distnod2,Vminp,Vmaxp,green,"", 2, -1);/* plot V at dendritic tip */
    //plot_v_nod(ct=dsgc,cn=1,distnod1,fymin,fymax,magenta,"Fdsgc_1_970", 2, -1);/* filtered distal APs */
    //plotfilt (ntimec,timec,-0.065);
    if (sq_wave) {
	square_wave_i(ndn(dsgc,dsgc1,soma),sq_wave,currinjd,stimtime,stimdur);
    } else {
	if (!notinit(currinja && currinja != 0)) {
		cclamp(ndn(dsgc,dsgc1,soma), currinjd, start=stimtime, dur=stimdur);
		//cclamp(ndn(dsgc,dsgc1,soma), currinjd*1.5, start=stimtime+2*stimdur, dur=stimdur);
		cclamp(ndn(dsgc,dsgc2,soma), currinja, start=0.01,               dur=stimdur);
		//cclamp(ndn(dsgc,dsgc2,soma), currinja*1.5, start=stimtime+2*stimdur, dur=stimdur);
	} else if (!notinit(currinjd) && currinjd != 0) {
		cclamp(ndn(dsgc,dsgc1,soma), currinjd, start=stimtime, dur=stimdur);
	}
    }
    // vclamp(ndn(dsgc,dsgc1,stimnod), -0.06, start=stimtime, dur=stimdur);
   } else {	// spike_pref==0, null dir
    timec[0]=0.0005; timec[1]= 0.020;
    if (!notinit(currinja) || !notinit(spot_stim)) Vmaxp = -0.00;
    if (bar_stim>0) Vmaxp = 0.00;
    if (dsgcarr==2) {
      plot_v_nod(ct=dsgc,cn=1,soma,Vminp,Vmaxp,blue,"",  30, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,soma,Vminp,Vmaxp,green,"", 45, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,100,0, -20,-32),Vminp,Vmaxp,ltred,"",   35, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,130,80,-20,-32),Vminp,Vmaxp,ltgreen,"", 36, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,100,0, -32,-48),Vminp,Vmaxp,ltred,"",   37, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,130,80,-32,-48),Vminp,Vmaxp,ltgreen,"", 38, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=2,findnodlocz(dsgc,2,100,0, -20,-32),Vminp,Vmaxp,gray,"",    40, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,findnodlocz(dsgc,2,130,80,-20,-32),Vminp,Vmaxp,ltblue,"",  41, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,findnodlocz(dsgc,2,100,0, -32,-48),Vminp,Vmaxp,gray,"",    42, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,findnodlocz(dsgc,2,130,80,-32,-48),Vminp,Vmaxp,ltblue,"",  43, -1);	/* plot V at soma */
    }
    else if (dsgcarr==102) {
      plot_v_nod(ct=dsgc,cn=1,soma,Vminp,Vmaxp,blue,"",  30, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,soma,Vminp,Vmaxp,green,"", 85, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,25,0,-32,-48),Vminp,Vmaxp,blue,"",  40, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,50,0,-32,-48),Vminp,Vmaxp,blue,"",  41, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,75,0,-32,-48),Vminp,Vmaxp,blue,"",  42, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,100,0,-32,-48),Vminp,Vmaxp,blue,"", 43, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,125,0,-32,-48),Vminp,Vmaxp,blue,"", 44, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,150,0,-32,-48),Vminp,Vmaxp,blue,"", 45, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,25,75,-32,-48),Vminp,Vmaxp,blue,"", 46, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,50,75,-32,-48),Vminp,Vmaxp,blue,"", 47, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,75,75,-32,-48),Vminp,Vmaxp,blue,"", 49, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,100,75,-32,-48),Vminp,Vmaxp,blue,"",50, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,125,75,-32,-48),Vminp,Vmaxp,blue,"",51, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=1,findnodlocz(dsgc,1,150,75,-32,-48),Vminp,Vmaxp,blue,"",52, -1);	/* plot V at soma */
    }
    else if (dsgcarr<4) {
      plot_v_nod(ct=dsgc,cn=2,soma,Vminp,Vmaxp,red,"", 1, -1);	/* plot V at soma */
      //plot_v_nod(ct=dsgc,cn=2,distnod2,Vminp,Vmaxp,brown,"", 1, -1);/* plot V at dendritic tip */
      plot_v_nod(ct=dsgc,cn=1,soma,Vminp,Vmaxp,blue,"", 1, -1);	/* plot V at soma */
      //plot_v_nod(ct=dsgc,cn=1,distnod1,Vminp,Vmaxp,green,"", 2, -1);/* plot V at dendritic tip */
      //plot_v_nod(ct=dsgc,cn=2,distnod2,fymin,fymax,magenta,"Fdsgc_1_1140", 2, -1);/* filtered distal APs */
      //plotfilt (ntimec,timec,-0.065);
    } else if (dsgcarr==4) {
      plot_v_nod(ct=dsgc,cn=1,soma,Vminp,Vmaxp,blue,"",  10, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,soma,Vminp,Vmaxp,green,"", 11, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=3,soma,Vminp,Vmaxp,cyan,"",  12, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=4,soma,Vminp,Vmaxp,red,"",   13, -1);	/* plot V at soma */
    } else if (dsgcarr==6) {
      plot_v_nod(ct=dsgc,cn=1,soma,Vminp,Vmaxp,blue,"",  10, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=2,soma,Vminp,Vmaxp,green,"", 11, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=3,soma,Vminp,Vmaxp,cyan,"",  30, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=3,findnodlocz(dsgc,3,260,20, -20,-32),Vminp,Vmaxp,ltred,"",   35, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=3,findnodlocz(dsgc,3,280,130,-20,-32),Vminp,Vmaxp,ltgreen,"", 36, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=3,findnodlocz(dsgc,3,260,20, -32,-48),Vminp,Vmaxp,ltred,"",   37, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=3,findnodlocz(dsgc,3,280,130,-32,-48),Vminp,Vmaxp,ltgreen,"", 38, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=4,findnodlocz(dsgc,4,260,20, -20,-32),Vminp,Vmaxp,gray,"",    40, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=4,findnodlocz(dsgc,4,280,130,-20,-32),Vminp,Vmaxp,ltblue,"",  41, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=4,findnodlocz(dsgc,4,260,20, -32,-48),Vminp,Vmaxp,gray,"",    42, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=4,findnodlocz(dsgc,4,280,130,-32,-48),Vminp,Vmaxp,ltblue,"",  43, -1);	/* plot V at soma */

      plot_v_nod(ct=dsgc,cn=4,soma,Vminp,Vmaxp,red,"",     45, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=5,soma,Vminp,Vmaxp,magenta,"", 50, -1);	/* plot V at soma */
      plot_v_nod(ct=dsgc,cn=6,soma,Vminp,Vmaxp,brown,"",   60, -1);	/* plot V at soma */

    }
    if (sq_wave) {
	square_wave_i(ndn(dsgc,dsgc2,soma),sq_wave,currinjd,stimtime,stimdur);
    } else {
        if (!notinit(currinja) && currinja != 0) {
                cclamp(ndn(dsgc,dsgc1,soma), currinjd, start=stimtime, dur=stimdur);
                //cclamp(ndn(dsgc,dsgc1,soma), currinjd*1.5, start=stimtime+2*stimdur, dur=stimdur);
		cclamp(ndn(dsgc,dsgc2,soma), currinja, start=0.01,               dur=stimdur);
		//cclamp(ndn(dsgc,dsgc2,soma), currinja*1.5, start=stimtime+2*stimdur, dur=stimdur);
	} else if (!notinit(currinjd) && currinjd != 0) {
                cclamp(ndn(dsgc,dsgc1,soma), currinjd, start=stimtime, dur=stimdur);
	}
    }
    if (!notinit(set_vclamp1)) {
	    vclamp(ndn(dsgc,dsgc1,soma), set_vclamp1, start=stimtime, dur=stimdur);
            plot_i_nod(ct=dsgc,dsgc1,soma,Imin,Imax,green,"",     31, -1);	/* plot I at soma */
    }
    if (!notinit(set_vclamp2)) {
	    vclamp(ndn(dsgc,dsgc2,soma), set_vclamp2, start=0, dur=10);
            plot_i_nod(ct=dsgc,dsgc2,soma,Imin,Imax,red,"",       46, -1);	/* plot I at soma */
    }
   }  // spike_pref == 0; null dir

  if (notinit(minten))     minten = -0.048;
  if (notinit(sinten))     sinten =  0.005;

  if (notinit(theta))   theta = 0;     /* orientation of bar */
  if (notinit(iroff))   iroff = 150;   /* r offset for inhib transducers */
  if (notinit(hroff))   hroff = 2000;  /* r offset for inhib transducers */

  if (!notinit(spot_stim) || !notinit(bar_stim)) {

   /* add light transducer to each On-bipolar cell */

    for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
      p = (photorec*)make_transducer(ndn(dbp1,cn,soma));
      p->xpos=npnt->xloc;
      p->ypos=npnt->yloc;
    }

    /* add light transducer with spatial offset to each Off-cell */

    hxoff = hroff *  cosdeg(theta);
    hyoff = hroff * -sindeg(theta);
    for(npnt=nodepnt; npnt=foreach(npnt,hbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
      p = (photorec*)make_transducer(ndn(hbp1,cn,soma));
      p->xpos=npnt->xloc + hxoff;
      p->ypos=npnt->yloc + hyoff;
    }
  }

   if (notinit(light_inhib)) light_inhib = 1;
   if (light_inhib) {

     /* Add light transducer to each small-field amacrine cell */
     /*  offset so that inhibition can be controlled separately */

     ixoff = iroff *  cosdeg(theta);
     iyoff = iroff * -sindeg(theta);
     for(npnt=nodepnt; npnt=foreach(npnt,ams,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       p = (photorec*)make_transducer(ndn(ams,cn,soma));
       p->xpos=npnt->xloc + ixoff;
       p->ypos=npnt->yloc + iyoff;
    }
  }

  if (notinit(velocity))   velocity = 500;
  if (notinit(barwidth))   barwidth = 50;
  if (notinit(barsep))       barsep = 200;
  if (notinit(sdia))           sdia = 100;
  if (notinit(move_left)) move_left = 1;

  if (make_movie==0) {
   if (bar_stim>0) { 

     if (dsgcarr<4 || dsgcarr==102) {
       bar_xmax = 180;
       bar_xmin = -80;
     } else if (dsgcarr==4) {
       bar_xmax = 360;
       bar_xmin = -80;
     } else if (dsgcarr==6) {
       bar_xmax = 600;
       bar_xmin = -200;
     }
     if (light_inhib==1) bar_xmax += iroff;
     bar_dist = bar_xmax - bar_xmin;

     /* Bar On-edge */
     if (make_dbp1) {
        stimdur = movebar (stimtime,0,0,bar_xmin,bar_xmax,barwidth,theta,velocity,sinten);  /* pref */
     }
     /* Bar Off-edge */
     if (make_hbp1) {
        stimdur = movebar (stimtime+barsep/velocity,0,0,bar_xmin+hroff,bar_xmax+hroff,barwidth,theta,velocity,sinten); 
     }
   //  stimdur = movebar (stimtime+0.1+bar_dist/velocity,0,0,bar_xmax,bar_xmin,barwidth,theta,velocity,sinten);  /* null */

     if (light_inhib==1) {
      plot_v_nod(ct=ams,cn=31,soma,smin,smax,blue,   "", 68, bpsize);/* plot V at dendritic tip */
      plot_v_nod(ct=ams,cn=79,soma,smin,smax,green,  "", 68, bpsize);/* plot V at dendritic tip */
      plot_v_nod(ct=ams,cn=88,soma,smin,smax,cyan,   "", 68, bpsize);/* plot V at dendritic tip */
      plot_v_nod(ct=ams,cn=134,soma,smin,smax,red,   "", 68, bpsize);/* plot V at dendritic tip */
      plot_v_nod(ct=ams,cn=32,soma,smin,smax,magenta,"", 68, bpsize);/* plot V at dendritic tip */
      plot_v_nod(ct=ams,cn=10,soma,smin,smax,brown,  "", 68, bpsize);/* plot V at dendritic tip */
    }

  }  // bar stim

    if (dsgcarr==2 || dsgcarr==102) {
       if (spot_stim==2) {
         plot_synrate(ct=dbp1,cn=findmid(dbp1,60,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   65,"",bpsize); /* plot rate out */
         plot_synrate(ct=dbp1,cn=findmid(hbp1,60,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   65,"",bpsize); /* plot rate out */
       } else {
     bar_xmin = 0;
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*0,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,green,  65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1.2,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,cyan,   65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1.4,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,red,    65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1.6,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,magenta,65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1.8,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,brown,  65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1.4,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,gray,   65,"",bpsize); /* plot rate out */


     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*1.4,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   62,"",bpsize); /* plot rate out */
       }
    }
    else {
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*0,0), soma,smin,smax,blue,   "", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1,0), soma,smin,smax,green,  "", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*2,0), soma,smin,smax,cyan,   "", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*3,0), soma,smin,smax,red,    "", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*4,0), soma,smin,smax,magenta,"", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*5,0), soma,smin,smax,brown,  "", 67, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*6,0), soma,smin,smax,gray,   "", 67, bpsize);/* plot V at dendritic tip */

     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*0,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*1,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,green,  65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*2,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,cyan,   65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*3,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,red,    65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*4,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,magenta,65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*5,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,brown,  65,"",bpsize); /* plot rate out */
     plot_synrate(ct=dbp1,cn=findmid(dbp1,bar_xmin+dsgc_spacing*6,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,gray,   65,"",bpsize); /* plot rate out */


//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*0,0), soma,smin,smax,blue,   "", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*1,0), soma,smin,smax,green,  "", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*2,0), soma,smin,smax,cyan,   "", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*3,0), soma,smin,smax,red,    "", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*4,0), soma,smin,smax,magenta,"", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*5,0), soma,smin,smax,brown,  "", 63, bpsize);/* plot V at dendritic tip */
//     plot_v_nod(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*6,0), soma,smin,smax,gray,   "", 63, bpsize);/* plot V at dendritic tip */

     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*0,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,blue,   62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*1,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,green,  62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*2,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,cyan,   62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*3,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,red,    62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*4,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,magenta,62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*5,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,brown,  62,"",bpsize); /* plot rate out */
     plot_synrate(ct=hbp1,cn=findmid(hbp1,bar_xmin+dsgc_spacing*6,0),1,rmin=0,rmax=200,1,fmin=0,fmax=1.5,gray,   62,"",bpsize); /* plot rate out */

  }
 } /* make movie == 0 */

  if (!notinit(bar_stim) && bar_stim>0) {
       endexp=stimtime+0.1+bar_dist/velocity+barwidth/velocity;
       trialdur = endexp;
  } else {
       trialdur = stimtime+3*stimdur+endwait;
       endexp = trialdur * ntrials;
  }
  if (!notinit(set_endexp)) endexp = set_endexp;

  stim_backgr(minten,start=0);		/* background */

  if (notinit(predur)) predur=0.0;
  setxmin=0;
  simtime=0-predur;
  step (predur);		// equilibrate model

  if (!notinit(spot_stim)) {			// set the spot locations
     if (dsgcarr==2 || dsgcarr==102) {
        if (spot_stim==2)  	{ if (notinit(sloc1))   sloc1 =  60; 
        			  if (notinit(sloc2))   sloc2 =  185; }
        if (spot_stim==3)  	{ if (notinit(sloc1))   sloc1 = 120; }
     }
     else if (dsgcarr==4) {
        if (spot_stim<3)        { if (notinit(sloc1))   sloc1 = 90;  }
        else if (spot_stim==3)  { if (notinit(sloc1))   sloc1 = 200; }
        if (spot_stim==2)       { if (notinit(sloc2))   sloc2 = 280; }
     }
     else if (dsgcarr==6) {
        if (spot_stim<3)        { if (notinit(sloc1))   sloc1 = 190;  }
        else if (spot_stim==3)  { if (notinit(sloc1))   sloc1 = 280; }
        if (spot_stim==2)       { if (notinit(sloc2))   sloc2 = 340; }
     }
  }


  for (n=t=0; n<ntrials; n++,t+=trialdur) {
 
    // simtime = 0; 
    if (!notinit(spot_stim)) {
        if (make_dbp1) stim_spot(sdia, sloc1,0,          sinten,stimtime+t,stimdur);
        if (make_hbp1) stim_spot(sdia, sloc1+hxoff,0,    sinten,stimtime+t+2*stimdur,stimdur);
        if (spot_stim==2) {
            if (make_dbp1) stim_spot(sdia, sloc2,0,      sinten,stimtime+t,stimdur);
            if (make_hbp1) stim_spot(sdia, sloc2+hxoff,0,sinten,stimtime+t+2*stimdur,stimdur);
	}
    //   plot_v_nod(ct=dbp1,cn=25,soma,Vminp,Vmaxp,brown,"", 5, -1);/* plot V at dendritic tip */

    } /* spot_stim */

    step (trialdur);
  }
}
