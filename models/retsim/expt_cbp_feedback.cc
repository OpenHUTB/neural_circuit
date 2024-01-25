/* Experiment cbp_feedback */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "ncio.h"
#include "stimfuncs.h"

int ct;
int ivplot;
int elnode;
int bp_morph;
int nod_r10;
int nod_r10b;
int nod_r8;
int nod_am1;
int nod_am2;

const char *celltype;

double kdr_cond;

double drmab;
double driab;
double dria;
double dvreva;
double dcmab;
double axdia;
double naax;
double naab;
double nahd;
double cars;
double krs;
double kax;
double set_drm;

double varicos_dia;

double predur;
double prestimdur;
double stimdur;
double tailcurdur;
double poststimdur;
double stimtime;
double stimfreq;

double istart;
double istop;
double istep;
double ipre;
double itail;
double iinj;
double dci;
double am1v;
double am2v;
double am1_time;
double am2_time;
double am1_dur;
double am2_dur;

double voltage;
double current;
double maxCurrent;
double pscal;
double elec_rs;
double elec_cap;

double set_timinc;
double set_ploti;

/*------------------------------------------------------*/

void defparams(void) { 

  setptr("celltype",    &celltype);
  setptr("ivplot",	&ivplot);
  setptr("bp_morph",	&bp_morph);
  setptr("nod_r8",	&nod_r8);
  setptr("nod_r10",	&nod_r10);
  setptr("nod_r10b",	&nod_r10b);
  setptr("nod_am1",	&nod_am1);
  setptr("nod_am2",	&nod_am2);

  setptr("set_drm",     &set_drm);
  setptr("drmab",       &drmab);
  setptr("driab",       &driab);
  setptr("dria",        &dria);
  setptr("dvreva",      &dvreva);
  setptr("dcmab",       &dcmab);
  setptr("axdia",       &axdia);
  setptr("naax",        &naax);
  setptr("naab",        &naab);
  setptr("nahd",        &nahd);
  setptr("cars",        &cars);
  setptr("krs",         &krs);
  setptr("kax",         &kax);
  setptr("kdr_cond",    &kdr_cond);
  setptr("varicos_dia", &varicos_dia);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("stimdur",	&stimdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststimdur",	&poststimdur);

  setptr("stimtime",	&stimtime);
  setptr("stimfreq",	&stimfreq);
  
  setptr("istart",&istart);
  setptr("istop",&istop);
  setptr("istep",&istep);
  setptr("ipre",&ipre);
  setptr("itail",&itail);
  setptr("iinj",&iinj);
  setptr("dci",&dci);
  setptr("am1v",&am1v);
  setptr("am2v",&am2v);
  setptr("am1_time",&am1_time);
  setptr("am2_time",&am2_time);
  setptr("am1_dur",&am1_dur);
  setptr("am2_dur",&am2_dur);

  setptr("maxCurrent", &maxCurrent);
  setptr("pscal", &pscal);

  setptr("set_timinc", &set_timinc);
  setptr("set_ploti", &set_ploti);

  nvalfile = "nval_cbp_flash.n";

  // dbp1_file = "morph_bp";		// default set in retsim.cc 
  dbp1_file = "morph_DB_111005_12_db4";	// morphology file in retsim.cc 
  dbp1_densfile = "dens_db4.n";
  // dbp1_thetay = -6;
  // dbp1_thetaz = -6;
  //
  if (!notinit(celltype)) {
       ct = find_ct(celltype);
  } else {
       ct = dbp1;
  }
  make_ct(ct);          /* make the cell type */
  set_ncel(ct,1);
}

/*------------------------------------------------------*/
   
void setparams(void)
{
#define AMARR 20

  if (n_am > 0) {
	  make_ct(am);
          setn(ct,AXARBT,BRANCHED);	/* connect to bipolar cell branched axonal arbor */
	  amxarr = (double *)emalloc(AMARR*sizeof(double));
	  amyarr = (double *)emalloc(AMARR*sizeof(double));
	  amxarr[0] = 0;
	  amyarr[0] = 0;
	  amxarr[1] = 50;
	  amyarr[1] = 0;
  }

  if (notinit(bp_morph)) setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_bp" */
  else                   setn(ct,MORPH,bp_morph);
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */
  setn(ct,NCOLOR,RCOLOR);	/* set cell display color from region */

  if (notinit(cbplam)) cbplam = 0.005; 		/* default complam for regions in density file */
  // if (notinit(dispsize)) dispsize = 80; 	/* default display size */
  if (notinit(node_scale)) node_scale = -3.05;  /* 3: nodenum, 0.05: small font */

  if (notinit(axdia))  axdia  = 1; 		/* multiplier for axon diameter */
  if (notinit(dvrev))  dvrev  = -0.065;		/* Vrev for dens_dbp1.n */
  if (notinit(dvst))   dvst   = dvrev;		/* Vstart for dens_dbp1.n */
  if (notinit(dvreva)) dvreva = -0.07;		/* Vrev for axon terminal in dens_dbp1.n */

  if (notinit(ivplot)) ivplot = 0; 		/* make I/V plot */
  if (notinit(pscal)) pscal   = 50e-12;	/* plot scale */

  if (!notinit(set_drm)) {             /* user set default Rm */
       setn(ct,NRM,set_drm); /* set default Rm */
       drm = set_drm;
  }
   if (notinit(drmab)) drmab = drm;              /* user set default Rm for axon branches */
   if (notinit(driab)) driab = dri;              /* user set default Ri for axon branches */
   if (notinit(dria))  dria  = dri;              /* user set default Ri for axon branches */
   if (notinit(dcmab)) dcmab = dcm;              /* user set default Cm for axon branches */
   if (notinit(naax))   naax = 0;                /* user set Na density in axon */
   if (notinit(naab))   naab = 0;                /* user set Na density in axon branches */
   if (notinit(nahd))   nahd = 900e-3;           /* user set Na high density region in axon */
   if (notinit(cars))   cars = 1e-3;             /* user set Ca at release site (R10) in axon */
   if (notinit(krs))    krs  = 0e-3;             /* user set K at in axon (R7-R9) in axon */
   if (notinit(kax))    kax  = 1e-3;             /* user set K at release site (R10) in axon */

  vna = 0.04;

  // dicafrac = 1;				/* remove Ca flux from ICa */

  // Set channel types in density file.
  // To change, check manual, and uncomment and modify these:  
  //
  //  _NA  = _NA5;		// NaV1.1 channel from Clancy & Kass (2004)
  //  _KA  = _K3;
  //  _KH  = _K4;
  //  _KIR = _K5;
  
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{
    // if (!notinit(kdr_cond))    celdens[dbp1][0][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    setsv (am,SCOND,1, 20e-10);

}

/*------------------------------------------------------*/

void onplot(void) {
    current = i(ndn(ct, 1, soma));
    voltage = v(ndn(ct, 1, soma));
    
    if (maxCurrent > current) {
      maxCurrent = current;
      //fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
    }
    
     // fprintf(stderr, "i: %g, maxi: %g\n", current, maxCurrent);    
}


void runexpt(void) {

    int cn, i, n, plnum, electrode_node;
    int colr;
    double dst, t, fmax,fmin;
    double rmin, rmax, plsize;
    double dtrial;
    double Vmin, Vmax;
    double Imin, Imax;
    double ipulse, mult, sdia, rm;
    node *npnt;
    photorec *p;
    char sbuf[30];

  if (!notinit(set_timinc)) timinc = set_timinc;
  else                      timinc = 2e-6;
  if (!notinit(set_ploti)) ploti = set_ploti;
  else 			   ploti = 1e-4;
  crit = 1e-10;

  cn = 1;
  electrode_node = 5000;
  setonplot(onplot);

  if (notinit(prestimdur))   prestimdur = 0.01;
  if (notinit(stimdur))      stimdur    = 0.01;
  if (notinit(poststimdur)) poststimdur = 0.05;

  if (notinit(stimtime))       stimtime = 0.05;
  if (notinit(stimfreq))       stimfreq = 10;

  dtrial = prestimdur+stimdur+poststimdur;
  endexp  = dtrial;

  if (notinit(itail)) itail  = 0e-12;
  if (notinit(iinj))   iinj  = 0e-12;
  if (notinit(dci))    dci   = 0e-12;
  if (notinit(ipre))   ipre  = dci;
  if (notinit(am1v))    am1v   = -0.07;
  if (notinit(am2v))    am2v   = -0.07;
  if (notinit(am1_time)) am1_time = 0.005;
  if (notinit(am2_time)) am2_time = 0.005;
  if (notinit(am1_dur))  am1_dur  = 0.005;
  if (notinit(am2_dur))  am2_dur  = 0.005;

  if (notinit(istart))      istart =   40e-12;
  if (notinit(istop))        istop =   240e-12;
  if (notinit(istep))        istep =   40e-12;

  if (notinit(elec_rs))  elec_rs   = 20e6;
  if (notinit(elec_cap)) elec_cap  = 1e-14;
  if (notinit(elnode))     elnode  = soma;

  if (elnode==electrode_node) {
      make_electrode  (nd(ct,cn=1,elnode), ndn(ct,cn=1,soma), elec_rs, elec_cap);
  }

  // make an amacrine cell that contacts the bipolar cell
  //  Use nval file to define synapse 
  
  if (!notinit(nod_am1) && nod_am1 > 0) {
     npnt = loc (nd(am,1,soma),100,0,-20);  
     make_sphere (npnt,sdia=2,rm=10e3);
     // p = (photorec*)make_transducer(npnt);
     // p->xpos=npnt->xloc;
     // p->ypos=npnt->yloc;
     connect_synapse (npnt,ndn(ct,1,nod_am1));
     vclamp(npnt, am1v, am1_time, am1_dur);
  }
  if (!notinit(nod_am2) && nod_am2 > 0) {
     npnt = loc (nd(am,2,soma),100,0,-20);  
     make_sphere (npnt,sdia=2,rm=10e3);
     // p = (photorec*)make_transducer(npnt);
     // p->xpos=npnt->xloc;
     // p->ypos=npnt->yloc;
     connect_synapse (npnt,ndn(ct,1,nod_am2));
     vclamp(npnt, am2v, am2_time, am2_dur);
  }

  // midcbp  = findmid(ct,0,0);

 
  if (ivplot) {
     graph_x(istop, istart);
     graph_y(pscal, -pscal);
     graph_init(); 
  } 
  else { 
     mult = 1;					// number of boutons
     if (notinit(nod_r8))  nod_r8  = 14;
     if (notinit(nod_r10)) nod_r10 = 104;
     plot_v_nod(ct, cn,nod_r10,Vmin = -0.07,   Vmax = -0.0,    colr=white,   "", 10, 0.75);
     if (!notinit(nod_r10b)) 
      plot_v_nod(ct, cn,nod_r10b,Vmin = -0.07,   Vmax = -0.0,    colr=red,   "", 10, 0.75);
     //plot_v_nod(ct, cn, 12,   Vmin = -0.07,   Vmax = -0.0,    colr=red,   "", 10, 0.75);
     //plot_v_nod(ct, cn, 14,   Vmin = -0.07,   Vmax = -0.0,    colr=green, "", 10, 0.75);
     plot_v_nod(ct, cn, soma, Vmin = -0.07,   Vmax = -0.0,    colr=blue,  "", 10, 0.75);
     plot_chan_current(ct, cn, nod_r8, NA, 5, mult=1.0, mult*3.0*10e-12, -mult*3.0*10e-12);
        sprintf (sbuf,"INa %d",nod_r8);
        plot_param (sbuf, colr=green,plnum=5,plsize=0.5);
     plot_chan_current(ct, cn, nod_r10, CA, 1, mult=1.0, mult*3.0*10e-12, -mult*3.0*10e-12);
        sprintf (sbuf,"ICa %d",nod_r10);
        plot_param (sbuf, colr=white,plnum=5,plsize=0.5);
     plot_chan_current(ct, cn, nod_r10, K, 1, mult=1.0, mult*3.0*10e-12, -mult*3.0*10e-12);
        sprintf (sbuf,"IK  %d",nod_r10);
        plot_param (sbuf, colr=red,plnum=5,plsize=0.5);
     plot_ca_nod(ct, cn, nod_r10, 1,      2e-5, brown, "Ca axterm", 2, 0.5); 
     plot_ca_nod(ct, cn, nod_r10, 15,     2e-5, cyan,   "Ca 10 axterm", 2, 0.5); 
     plot_cabufb_nod(ct, cn, nod_r10, 15, 2e-5, colr=gray, "Cabuf axterm", 2, 0.5);

  if (!notinit(nod_am1) && nod_am1 > 0) {
     plot_v_nod(am, 1, soma, Vmin = -0.07,   Vmax = -0.03,    colr=blue,  "", 9, 0.2);
     plot_chan_cond(ct, cn, nod_am1, GABA, 1, mult=1.0, mult*5*100e-12, -mult*1.5*10e-12);
        sprintf (sbuf,"GGABA %d",nod_am1);
        plot_param (sbuf, colr=magenta,plnum=8,plsize=0.2);
  }

     plot_i_nod(ct, cn, soma, Imin = 0e-12, Imax = 100e-12, colr=blue, "", 1, 0.2); 
  }

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 0.05;
  simtime = 0 - predur;
  cclamp(ndn(ct,cn,soma), ipre, simtime, predur+stimtime);
  step (predur);

  if (iinj > 0 || dci>0) {
     sine_wave_i(ndn(ct,cn,elnode),stimfreq,iinj,stimtime,stimdur,0.001);
     cclamp (ndn(ct,cn,elnode),dci,stimtime,stimdur);
     endexp = stimtime+stimdur+poststimdur;
     step (endexp);
  } 
  else {

    for (i=0,ipulse=istart; ipulse <= istop; i++,ipulse += istep) {

     simtime = 0;
     //graph_pen(i+1,i+1,i+1);
     cclamp(ndn(ct,cn,elnode), ipre,   simtime,  prestimdur);
     step (prestimdur);

     cclamp(ndn(ct,cn,elnode), ipulse, simtime,  stimdur);
     step (stimdur);

     cclamp(ndn(ct,cn,elnode), ipre, simtime,  poststimdur);
     step (poststimdur);
     
    }
  }
  
}




