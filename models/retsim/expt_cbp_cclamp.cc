/* Experiment cbp_cclamp */
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

double voltage;
double current;
double maxCurrent;
double pscal;
double elec_rs;
double elec_cap;

/*------------------------------------------------------*/

void defparams(void) { 

  setptr("celltype",    &celltype);
  setptr("ivplot",	&ivplot);
  setptr("bp_morph",	&bp_morph);

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

  setptr("maxCurrent", &maxCurrent);
  setptr("pscal", &pscal);

  nvalfile = "nval_cbp_flash.n";

  // dbp1_file = "morph_bp";		// default set in retsim.cc 
  dbp1_file = "morph_DB_111005_12_db4";	// morphology file in retsim.cc 
  dbp1_densfile = "dens_db4.n";
  // dbp1_thetay = -6;
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
  DEND = R1;
  DENDP = R2;

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
    double ipulse;

  timinc = 10e-6;
  ploti = 1e-3;
  crit = 1e-10;

  cn = 1;
  electrode_node = 5000;
  setonplot(onplot);

  if (notinit(prestimdur))   prestimdur = 0.02;
  if (notinit(stimdur))      stimdur    = 0.1;
  if (notinit(poststimdur)) poststimdur = 0.05;

  if (notinit(stimtime))       stimtime = 0.05;
  if (notinit(stimfreq))       stimfreq = 10;

  dtrial = prestimdur+stimdur+poststimdur;
  endexp  = dtrial;

  if (notinit(itail)) itail  = 0e-12;
  if (notinit(iinj))   iinj  = 0e-12;
  if (notinit(dci))    dci   = 0e-12;
  if (notinit(ipre))   ipre  = dci;

  if (notinit(istart))      istart =   40e-12;
  if (notinit(istop))        istop =   240e-12;
  if (notinit(istep))        istep =   40e-12;

  if (notinit(elec_rs))  elec_rs   = 20e6;
  if (notinit(elec_cap)) elec_cap  = 1e-14;
  if (notinit(elnode))     elnode  = soma;

  if (elnode==electrode_node) {
      make_electrode  (nd(ct,cn=1,elnode), ndn(ct,cn=1,soma), elec_rs, elec_cap);
  }
	
	
  // midcbp  = findmid(ct,0,0);

 
  if (ivplot) {
     graph_x(istop, istart);
     graph_y(pscal, -pscal);
     graph_init(); 
  } 
  else { 
     plot_v_nod(ct, cn, 12,   Vmin = -0.08,   Vmax = -0.01,    colr=red,   "", 2, 1.0);
     plot_v_nod(ct, cn, 6,    Vmin = -0.08,   Vmax = -0.01,    colr=green, "", 2, 1.0);
     plot_v_nod(ct, cn, soma, Vmin = -0.08,   Vmax = -0.01,    colr=blue,  "", 2, 1.0);
     plot_i_nod(ct, cn, soma, Imin = -iinj*2, Imax = iinj*2, colr=magenta, "", 1, 0.5); 
     // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 
  }

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 0.15;
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




