/* Experiment gc_cbp_snr */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "stimfuncs.h"

int flseed;

int stimtype;
int run_vclamp;
int node_dist;
int zerotime;
int itransducer;
double ampar;

double temp_freq;
double ntrials;
double sdur;
double spotdia;
double flashdur;
double stimtime;
double minten;
double inten_mul;
double scontrast;
double scontrast2;
double setploti;
double predur;
double somadia;
double gca_arbscale;
double vhold;
double cbplam;
double dbp1_mp;
double dbp1_mr;
double dbp1_stau;
double dbp1_k1_offm;
double dbp1_gca_reg;
double g_dbp1_dbp1;
double g_dbp1_gca;
double trialdur;
double db1_axdia;
double vbias;

int rec_ct;
int rec_cn;

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("flseed", &flseed);
  setptr("run_vclamp",&run_vclamp);
  setptr("node_dist", &node_dist);
  setptr("zerotime",  &zerotime);
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("sdur",      &sdur);
  setptr("spotdia",   &spotdia);
  setptr("flashdur",  &flashdur);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("inten_mul", &inten_mul);
  setptr("scontrast", &scontrast);
  setptr("scontrast2", &scontrast2);
  setptr("setploti",  &setploti);
  setptr("predur",    &predur);
  setptr("stimtype",  &stimtype);
  setptr("gca_arbscale",&gca_arbscale);
  setptr("somadia",   &somadia);
  setptr("vhold",     &vhold);
  setptr("cbplam",    &cbplam);
  setptr("dbp1_mp",   &dbp1_mp);
  setptr("dbp1_mr",   &dbp1_mr);
  setptr("dbp1_stau", &dbp1_stau);
  setptr("db1_axdia", &db1_axdia);
  setptr("dbp1_gca_reg",&dbp1_gca_reg);
  setptr("g_dbp1_gca",&g_dbp1_gca);
  setptr("g_dbp1_dbp1",&g_dbp1_dbp1);
  setptr("ampar",     &ampar);
  setptr("trialdur",  &trialdur);
  setptr("itransducer",&itransducer);
  setptr("dbp1_k1_offm",&dbp1_k1_offm);
  setptr("vbias",     &vbias);

  if (notinit(dbp1_gca_reg)) dbp1_gca_reg = R2;
  if (notinit(g_dbp1_gca)) g_dbp1_gca = 50e-12;
  if (notinit(gca_arbscale)) gca_arbscale = 1.0;

  nvalfile = "nval_gc_cbp_snr.n";
  dbp1_densfile = "dens_dbp1_gc_cbp_snr.n";
  // gca_densfile = "dens_gc_cbp_snr.n";	/* density file for cclamp, spiking */
  gca_densfile = "dens_gca_zero.n";
  chanparamsfile = "chanparams_gc_cbp_snr";

  if (notinit(dbp1_mp))     dbp1_mp = 8;	/* maximum rr pool size */
  if (notinit(dbp1_mr))     dbp1_mr = 50;	/* maximum syn release rate */
  if (notinit(dbp1_stau)) dbp1_stau = 2;	/* mini quanta tau */
  if (notinit(dbp1_k1_offm)) dbp1_k1_offm = 0.025; /* dbp1 K1 offsetm */
  if (notinit(ampar))         ampar = xampa5;  /* dbp1->gca postsyn receptor, non-desens */
  if (notinit(g_dbp1_dbp1)) g_dbp1_dbp1 = 2e-10; /* dbp1 gap juncition conductance */

  _CA_L = _CA1;                 /* set type of L-type calcium channel for dens_dbp1_gc_cbp_snr.n */
  _CA_T = _CA7;                 /* set type of T-type calcium channel for dens_dbp1_gc_cbp_snr.n */
  dscavg = 1e5;                /* Ca sensitivity for vesicle release */
  dscaeg = 2;			/* Ca sensitivity exp gain for ves release */

  gc_vs = -0.065;		/* Vstart for gca */
  // if (notinit(ampar))         ampar = xampa1;	/* dbp1->gca postsyn receptor, desensitizing */
 
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  make_cones= 0;        /* make cones, dbp, gc */
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 1;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 1;
  make_gcb  = 0;
  make_dsgc = 0;

  if (notinit(somadia)) somadia = 15;

  gca_file = "morph_beta8c";

  if(notinit(rec_ct)) rec_ct = gca;

  DENDD     = R_1;
  DEND_DIST = R_1;
  DEND      = R_2;  /* definitions of regions for dens_ file */
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

  if (notinit(ath_dia)) ath_dia = 0.25;
  if (notinit(cbplam)) cbplam = 0.5;		/* make cbp with 2 compartments */
  setn(gca,NCOLOR,1);
  if (notinit(db1_axdia)) db1_axdia = 0.5;      //  dia mult.for dbp1 axon in dens_dbp1_gc_cbp_snr.n
  if (notinit(vbias)) vbias = 0;      		// bias voltage for high-s.d. stimulus 

  //if (notinit(arrsiz)) arrsiz = 300;
  //if (arrsiz==100) {
  //  setsv (dbp1,SCOND,1, 25e-10);
  //} 
}

/*------------------------------------------------------*/
#define SPOT_ON_OFF 2

double move_stim(int stimtype, double spotdia, double scontrast, double stimtime, double spotdur)
{
   int s, send;
   int npixels, nframes;
   double x,y,stimend,orient;
   double *rndarr = NULL;

 if (spotdia < 0) spotdia = max(xarrsiz,yarrsiz) * 1.05; 

 x = y = 0;
 if (stimtype==1) {			// flashed spot 
    stim_spot(spotdia, x, y, scontrast, stimtime, spotdur);
    stimend = stimtime + spotdur;
 }
 else if (stimtype==SPOT_ON_OFF) {	// flashed spot, on, off
    stim_spot(spotdia, x, y, scontrast, stimtime,         spotdur);
    stim_spot(spotdia, x, y, -scontrast, stimtime+spotdur, spotdur);
    stimend = stimtime + 2*spotdur;

 } else if (stimtype==3) {		// spot random flicker modulated
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast,       1.0,scontrast, stimtime,spotdur);
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast,       1.0,scontrast, stimend,spotdur);
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast2+vbias,1.0,scontrast2,stimend,spotdur);
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast2+vbias,1.0,scontrast2,stimend,spotdur);
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast,       1.0,scontrast, stimend,spotdur);
   stimend = spot_flicker(spotdia,x,y,temp_freq,flseed,scontrast,       1.0,scontrast, stimend,spotdur);

 } else if (stimtype==4) {		// random checkerboard 
   npixels = 16;
   stim_checkerboard(spotdia, spotdia, npixels, npixels, orient=0, x=0, y=0,
                     temp_freq, 0, scontrast, stimtime, spotdur, &rndarr, &nframes, flseed);
   efree(rndarr);
   stimend = stimtime + spotdur;

 } else if (stimtype==5) {		// spot sine modulated
    stimend = spot_sine(spotdia, x, y, temp_freq, 1.0, scontrast, stimtime, spotdur);

 } else if (stimtype==6) {		// move spot
	   int cellrad;

   cellrad = xarrsiz / 2;
   send = int(cellrad/spotdia + 0.5) * 2;
   for (s=0; s<send; s++) {
           stim_spot(spotdia, s*spotdia - cellrad, y, scontrast, stimtime+s*spotdur, spotdur);
   }
   stimend = spotdur * send;

 } else if (stimtype==7) {		// flashed annulus, proximal, distal

 }

   return stimend;
}

/*------------------------------------------------------*/

void runexpt(void)

{
    int i, ct, cn, n, plnum;
    int colr;
    int midcone,midcbp,midcbp2,midcbp3;
    int synin1, synin2;
    double t, fmax,fmin;
    double rmin, rmax, plsize;
    double stimend;
    double Vmin, Vmax;
    double Imin, Imax;
    char plbuf[20];
    node *npnt;
    photorec *p;

  if (notinit(temp_freq)) temp_freq = 2;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  ploti = 1e-4;

  /* add light transducer to each bipolar cell */

   if (notinit(itransducer)) itransducer = 1;

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     if (itransducer) p = (photorec*)make_itransducer(ndn(dbp1,cn,soma));
     else             p = (photorec*)make_transducer (ndn(dbp1,cn,soma));
     p->xpos=npnt->xloc;
     p->ypos=npnt->yloc;
   }

  if (notinit(sdur))           sdur = 0.1;      /* stimulus duration */
  if (notinit(flashdur))   flashdur = 0;        /* duration of initial flash */
  if (notinit(spotdia))     spotdia = 300;      /* spot diameter */
  if (notinit(stimtime))   stimtime = 0.05;    /* stimulus time */
  if (itransducer) {
    if (notinit(minten))         minten = 0e-12;  /* background intensity */
    if (notinit(scontrast))  scontrast  = 1e-12;  /* contrast */
    if (notinit(scontrast2)) scontrast2 = 2e-12;  /* contrast 2 */
  } else {
    if (notinit(minten))       minten = -0.05;    /* background intensity */
    if (notinit(scontrast)) scontrast = 0.005;    /* contrast */
    if (notinit(scontrast2)) scontrast = 0.010;    /* contrast 2 */
  }
  if (!notinit(setploti))     ploti = setploti; /* plot time increment */
  if (notinit(predur))       predur = 0.05;      /* equilib time */
  if (notinit(flseed))       flseed = 3751;     /* seed for random flicker */
  if (notinit(stimtype))   stimtype = 1;        /* stimulus type */
  if (notinit(ntrials))     ntrials = 1;
  if (notinit(zerotime))   zerotime = 1;	/* return to zero time after each trial */

  if (notinit(trialdur)) {
	  if (flashdur>0) trialdur = 2*stimtime + 2*flashdur + sdur;
	  else	          trialdur = stimtime + 6 * sdur + 0.05;
  }

  if (zerotime>0)		/* if reset time to zero between trials */
       endexp  = trialdur;
  else endexp  = trialdur * ntrials;

  midcone = findcell(xcone,0,0);
  midcbp  = findcell(dbp1,100,20);
  midcbp2 = findcell(dbp1,20,0);
  midcbp3 = findcell(dbp1,20,10);
  if (getn(xcone,NMADE)>0) {
    synin1 = ncel_in(dbp1,midcbp,xcone);
    synin2 = ncel_in(dbp1,midcbp2,xcone);
    if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
    if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d ncones %d\n",  midcbp,synin1);
    if (ninfo >=1) fprintf (stderr,"# mid cbp2 # %d ncones %d\n",  midcbp2,synin2);
  }

  simtime = -predur;		/* start model before time 0, must be ahead of stim_backgr() */
  setxmin = 0;			/* start plotting at time 0 */

  stim_backgr(minten,simtime);                                     /* background */

  if (disp & DSTIM) {			/* display stimulus */
         double t, dscale, disp_end, starttime;

    if (flashdur==0) stimend = 0;
    //else
    //  stimend = move_stim(SPOT_ON_OFF, spotdia, scontrast, simtime+predur+stimtime, flashdur);
    stimend   = move_stim(stimtype,    spotdia, scontrast, stimend+stimtime, sdur);

    if (makestim) return;
    if (notinit(dispsize)) {
        maxsize = max(xarrsiz,yarrsiz);
        dispsize=maxsize*(1+disp_margin);
    }
    display_size(dispsize);
    if (notinit(disp_calib_len)) disp_calib_len = set_int_val(dispsize/5);
    if (disp_calib_len > 0) drcalib(0.97,0.05,disp_calib_len,dispsize,white);

    disp_end = stimtime+stimend+0.05;
    for (starttime=simtime,t=0; t<disp_end; starttime = t, t+= 0.001) {
         if (itransducer) display_stim(starttime, t, dscale=4, 5e-12, -5e-12);
	 else             display_stim(starttime, t, dscale=4, -0.04, -0.05);
         // display_stim(starttime, t, dscale=4, -0.025, -0.035); 
         //         // display_stim(t, dscale=4, -0.035, -0.045); 
         simwait(0.10);
    }
    return;
  }

  plot_v_nod(xcone,midcone,soma,Vmin=-.040,Vmax =-.035,colr=cyan,"", 45, -1);     /* plot Vcone */
  //plot_synrate(findsynlocr(ct=xcone,cn=midcone),rmin=0,rmax=400,colr=magenta,-1,"",1);  /* cone release rate */
  //plot_synves(findsynlocr(ct=xcone,cn=midcone),rmin=0,rmax=2e-4,colr=red,-1,"",1);       /* cone vesicles */
  //plot_syncond(findsynlocr(ct=xcone,cn=midcone),rmin=0,rmax=10e-12,colr=green,-1,"",1); /* cone syn cond */
  plot_synrate(xcone,midcone,rmin=0,rmax=200,colr=magenta,40,"",1);

  plot_v_nod(dbp1,midcbp,  soma,Vmin=-0.05,Vmax =-0.03, colr=green,"",   18, 0.5);   /* plot Vcbp */
  plot_v_nod(dbp1,midcbp,  7,Vmin=-0.05,Vmax =-0.03, colr=blue,"",    18, 0.5);   /* plot Vcbp axterm */
  
  plot(CA, 1, ndn(dbp1,midcbp,7), fmax=10e-6, fmin=0); 	/* display Ca in dbp1 */
       sprintf (plbuf,"Cai_dbp1_%d",midcbp);
 	plot_param(plbuf, colr=brown,plnum=16,plsize=0.5);

  plot_synrate(dbp1,midcbp, rmin=0,rmax=200,colr=magenta,15,"",0.5);     /* dbp1 release rate */
  plot_synrate(dbp1,midcbp2,rmin=0,rmax=200,colr=blue,15,"",0.5);        /* dbp1 release rate */
  plot_synrate(dbp1,midcbp3,rmin=0,rmax=200,colr=green,15,"",0.5);        /* dbp1 release rate */
  //plot_synves (dbp1,midcbp,rmin=0,rmax=2e-3,colr=red,12,"",0.5);       /* dbp1 vesicles */
  plot_syncond(dbp1,midcbp, rmin=0,rmax=1e-11,colr=magenta,12,"",0.5);   /* dbp1 syn cond */
  plot_syncond(dbp1,midcbp2,rmin=0,rmax=1e-11,colr=blue,12,"",0.5);      /* dbp1 syn cond */
  plot_syncond(dbp1,midcbp3,rmin=0,rmax=1e-11,colr=green,12,"",0.5);      /* dbp1 syn cond */

  if (notinit(node_dist)) node_dist = 197;
  if (notinit(run_vclamp)) run_vclamp = 0;
  if (run_vclamp > 0) {
      if (notinit(vhold)) vhold = -0.068;
      vclamp (ndn(ct=gca,cn=1,soma),    vhold, simtime,  predur+endexp);	/* for first trial */
      plot_i_nod(ct=gca,cn=1,soma,      Imin= -2e-10,Imax=0e-12,blue,"", 1, 0.5);
      //plot_v_nod(ct=gca,cn=1,soma,      Vmin=-.07,Vmax=-.04,green,"", 2, 0.35);
      //plot_v_nod(ct=gca,cn=1,node_dist, Vmin=-.07,Vmax=-.04,brown,"", 2, 0.35);
  } else {
      plot_v_nod(ct=gca,cn=1,soma,Vmin=-.070,Vmax =-.03,colr=blue,"", 1, -1);	      /* plot Vgc */
      // plot_v_nod(ct=gca,cn=1,1514,Vmin=-.070,Vmax =-.040,colr=brown,"", 1, -1);      /* plot Vgc axon */
      // if (make_gca && getn(gca,BIOPHYS)) {plot(CA, 1, ndn(gca,1,soma), fmax=0.5e-6, fmin=0); 
      // 		plot_param("Cai_gca", colr=yellow,plnum=0,plsize=0.3);}
  }

//  synaptau = 0.001;             // set synapses to run with shorter time constant

  step (predur);

//  synaptau = 1.0;               // set synapses to run with normal time constant

  stimend = 0;

  for (i=0; i<ntrials; i++) {
       double start, dur;
//     if (zerotime) simtime = 0;
//     if (i>0) {
//         if (flashdur==0) stimend = simtime;
//         if (run_vclamp) vclamp (ndn(ct=gca,cn=1,soma),    vhold, simtime,  predur+endexp);
//        else
//	   stimend = move_stim(SPOT_ON_OFF, spotdia, scontrast, simtime+stimtime, flashdur);
//	 stimend   = move_stim(stimtype,    spotdia, scontrast, stimend+stimtime, sdur);
//     }
	 stimend   = move_stim(stimtype,    spotdia, scontrast, stimend+stimtime, sdur);
     step(trialdur);
  }

}
