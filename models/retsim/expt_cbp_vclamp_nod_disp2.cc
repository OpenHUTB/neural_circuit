/* Experiment cbp_vclamp */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "ncio.h"

int ct;
int ivplot;
int elnode;
int dbpair;
int bp_morph;

int stimloc;
int tip1;
int recpnt1;
int recpnt2;
int recpnt3;

const char *celltype;

double kdr_cond;

double drmab;
double dria;
double driab;
double dcmab;
double naax;
double naab;
double nahd;
double axdia;
double set_drm;

double varicos_dia;

double predur;
double prestepdur;
double stepdur;
double tailcurdur;
double poststepdur;

double vstart;
double vstop;
double vstep;
double vleak;
double vhold;
double tailvolt;
double gvrev;
double dvreva;

double voltage;
double current;
double current2;
double idiff;
double iscal;
double gscal;
double outward;
double set_timinc;
double set_ploti;
double elec_rs;
double elec_cap;


char savefile[30] = {0};

/*------------------------------------------------------*/

void defparams(void) {

  setptr("celltype",    &celltype);
  setptr("ivplot",	&ivplot);
  setptr("outward",	&outward);
  setptr("elnode",	&elnode);
  setptr("bp_morph",	&bp_morph);

  setptr("set_drm",     &set_drm);
  setptr("drmab",	&drmab);
  setptr("dria",	&dria);
  setptr("driab",	&driab);
  setptr("dcmab",	&dcmab);
  setptr("naax",	&naax);
  setptr("naab",	&naab);
  setptr("nahd",	&nahd);
  setptr("axdia",	&axdia);
  setptr("kdr_cond",    &kdr_cond);
  setptr("varicos_dia", &varicos_dia);

  setptr("predur",	&predur);
  setptr("prestepdur",	&prestepdur);
  setptr("stepdur",	&stepdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststepdur",	&poststepdur);
  
  setptr("vstart",&vstart);
  setptr("vstop",&vstop);
  setptr("vstep",&vstep);
  setptr("vleak",&vleak);
  setptr("vhold",&vhold);
  setptr("tailvolt",&tailvolt);
  setptr("gvrev",&gvrev);
  setptr("dvreva",&dvreva);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  setptr("set_timinc", &set_timinc);
  setptr("set_ploti",  &set_ploti);
  setptr("elec_rs",    &elec_rs);
  setptr("elec_cap",   &elec_cap);

  setptr("stimloc",&stimloc);
  setptr("tip1",&tip1);
  
  setptr("recpnt1",&recpnt1);
  setptr("recpnt2",&recpnt2);
  setptr("recpnt3",&recpnt3);
  
  nvalfile = "nval_cbp_flash.n";

  // dbp1_file = "morph_bp";		// default set in retsim.cc 
   dbp1_file = "morph_DB_111005_12";	// morphology file in retsim.cc 
  dbp1_densfile = "dens_db3a.n";
  dbp1_densfile2 = "dens_db3x.n";
  // dbp1_thetay = -6;
}

/*------------------------------------------------------*/
   
void setparams(void)
{

  if (!notinit(celltype)) {
       ct = find_ct(celltype);
  } else {
       ct  = dbp1;
  }
  make_ct(ct);          /* make the cell type */

  if (n_dbp1<=0) n_dbp1 = 2;
  set_ncel(ct,n_dbp1);
  dbpair = (n_dbp1 == 2);

  if (notinit(bp_morph)) setn(ct,MORPH,0);	/* set cell morphology from file, default = "morph_bp" */
  else                   setn(ct,MORPH,bp_morph);
  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_bp" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */
  setn(ct,NCOLOR,RCOLOR);	/* set cell display color from region */

  if (notinit(cbplam)) cbplam = 0.005; 		/* default complam for regions in density file */
  // if (notinit(dispsize)) dispsize = 80; 	/* default display size */
  if (notinit(node_scale)) node_scale = -3.05;  /* 3: nodenum, 0.05: small font */

  if (notinit(axdia)) axdia = 1; 		/* multiplier for axon diameter in dens_xxx.n */
  if (notinit(dvrev)) dvrev = -0.065;		/* Vrev for dens_dbp1.n */
  if (notinit(dvst))   dvst = -0.068;		/* Vstart for dens_dbp1.n */
  if (notinit(dvreva)) dvreva = -0.070;		/* Vrev for dens_dbp1.n, axon */

  if (notinit(ivplot)) ivplot = 0; 		/* make I/V plot */
  if (notinit(outward)) outward = 1; 		/* >= 0 => calc K cond, otherwise Na cond */
  if (notinit(iscal)) iscal   = 1e-9;		/* plot scale */
  if (notinit(gscal)) gscal   = 2e-9;		/* plot scale */

  if (!notinit(set_drm)) {             		/* user set default Rm */
       setn(ct,NRM,set_drm); 			/* set default Rm */
       drm = set_drm;
  }
  if (notinit(drmab)) drmab = drm;		/* user set default Rm for axon branches */
  if (notinit(dria))  dria  = dri;		/* user set default Ri for axon hi-res region */
  if (notinit(driab)) driab = dri;		/* user set default Ri for axon branches */
  if (notinit(dcmab)) dcmab = dcm;		/* user set default Cm for axon branches */
  if (notinit(axdia)) axdia = 1;		/* user set default Dia for axon branches */
  if (notinit(naax))   naax = 0;		/* user set Na density in axon */
  if (notinit(naab))   naab = 0;		/* user set Na density in axon branches */
  if (notinit(nahd))   nahd = 900e-3;		/* user set Na high density region in axon */
 
  /* set recording points */
 
  //if (!notinit(tip1)) {label(ndn(dbp1,1,dendn_node(dbp1,tip1)),red);}
  if (!notinit(recpnt1)) recpnt1 = 0;           /*record from node through commandline*/ 
  if (!notinit(recpnt2)) recpnt2 = 802;  
  if (!notinit(recpnt3)) recpnt3 = 803;  
 
  
  //vna = 0.04;
  dicafrac = 0;		/* remove Ca pump from ICa */

  // Set channel types in density file.
  // To change, check manual, and uncomment and modify these:  
  //
  //  _NA  = _NA5;
  //  _KA  = _K3;
  //  _KH  = _K4;
  //  _KIR = _K5;
  
}

/*------------------------------------------------------*/

void addlabels(void)
{
   if (!notinit(stimloc)) {label(ndn(dbp1,1,dendn_node(dbp1,stimloc)),yellow, "stim");}
   // if (!notinit(tip1)) {label(ndn(dbp1,1,dendn_node(dbp1,tip1)),yellow);}
   //label(ndn(dbp1,1,dendn_node(dbp1,0)),blue); 
   //label(ndn(dbp1,1,0),red);
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{
    int cn,dn;

    // if (!notinit(kdr_cond))    celdens[dbp1][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    ndens[dbp1][cn=1] = 0;              // set cn 1 to use dbp1_densfile
    ndens[dbp1][cn=2] = 1;              // set cn 2 to use dbp1_densfile2

}

/*------------------------------------------------------*/

void runonexit (void)
{
  if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/

bool flag = false;
double maxCurrent;
double cond,gmax;

void onplot(void) {
    current              = i(ndn(ct, 1, elnode));
    if (dbpair) current2 = i(ndn(ct, 2, elnode));
    else        current2 = 0;
    idiff = current - current2;

    voltage  = v(ndn(ct,  1, soma));
   
  if (flag) {
   if (outward<=0) { 
      cond = idiff / (voltage - gvrev);
      if (cond > gmax) gmax = cond;
      if (maxCurrent > idiff) {		// inward current
        maxCurrent = idiff;
        //fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   } else {  				// outward current 
      cond = idiff / (voltage - vk);
      if (cond > gmax) gmax = cond;
      if (maxCurrent < idiff) {	
        maxCurrent = idiff;
        //fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   } 
     // fprintf(stderr, "i: %g, maxi: %g\n", current, maxCurrent);    
  }
  else cond = 0;
}

/*------------------------------------------------------*/

void runexpt(void) {

    int cn, i, n, plnum,electrode_node;
    int colr,pl;
    double dst, t, fmax,fmin;
    double rmin, rmax, plsize;
    double dtrial;
    double Vmin, Vmax, vstopx;
    double imax,imin,Imin,Imax,Gmin,Gmax;
    double vpulse,sign;
  
  if (!notinit(set_timinc)) timinc = set_timinc; 
  else                      timinc = 5e-6;
  
  if (!notinit(set_ploti))  ploti = set_ploti; 
  else 			    ploti = 1e-4;
  // crit = 1e-12;
  crit = 1e-10;

  cn = 1;
  electrode_node = 5000;
  vstopx = 0.02;

  setonplot(onplot);

  if (notinit(prestepdur))   prestepdur = 0.02;
  if (notinit(stepdur))      stepdur    = 0.10;
  if (notinit(tailcurdur))   tailcurdur = 0.02;
  if (notinit(poststepdur)) poststepdur = 0.05;

  //dtrial = prestepdur+stepdur+tailcurdur+poststepdur;
  dtrial = prestepdur+stepdur+tailcurdur;
  endexp  = dtrial;

  if (notinit(vhold))        vhold = -0.07;
  if (notinit(vstart))      vstart = -0.06;
  if (notinit(vstop))        vstop =  0.02;
  if (notinit(vstep))        vstep =  0.005;
  if (notinit(vleak))        vleak =  -0.005;
  if (notinit(tailvolt)) tailvolt  = vhold;
  if (notinit(gvrev)) 	    gvrev  = vna;
  
  if (notinit(elec_rs))  elec_rs   = 20e6;
  if (notinit(elec_cap)) elec_cap  = 1e-14;
  if (notinit(elnode))   elnode  = electrode_node;

  // midcbp  = findmid(ct,0,0);

  if (elnode==electrode_node) {
      make_electrode             (nd(ct,cn=1,elnode), ndn(ct,cn=1, dendn_node(ct,recpnt1)), elec_rs, elec_cap);
      if (dbpair) make_electrode (nd(ct,cn=2,elnode), ndn(ct,cn=2, dendn_node(ct,recpnt1)), elec_rs, elec_cap);
  }

  if (ivplot) {
     graph_x(vstop, vstart);
     graph_y(iscal, -iscal);
     graph_y(gscal, 0);
     graph_init(); 
  } 
  else { 

     if      (outward>=1)		{ Imax = iscal; Imin = 0;      Gmax=gscal; Gmin=0;}
     else if (outward>0 && outward<1)   { Imax = iscal; Imin = -iscal; Gmax=gscal; Gmin= -gscal;}
     if      (outward==0)		{ Imax = 0;     Imin = -iscal; Gmax=gscal; Gmin=0;}
     else if (outward< 0)		{ Imax = iscal; Imin = -iscal; Gmax=gscal; Gmin= -gscal;}
     
     
     if (!notinit(recpnt1)) 
	plot_v_nod(ct,cn=1,dendn_node(ct,recpnt1),Vmin = min(vhold,vstop),Vmax=max(vhold,vstopx), colr=magenta,"", 8, 0.5);
     if (!notinit(recpnt2)) 
	plot_v_nod(ct,cn=1,dendn_node(ct,recpnt2),Vmin = min(vhold,vstop),Vmax=max(vhold,vstopx),colr=blue,"", 7, 0.5);
     if (!notinit(recpnt3)) 
	plot_v_nod(ct,cn=1,dendn_node(ct,recpnt3),Vmin = min(vhold,vstop),Vmax = max(vhold,vstopx),colr=cyan,"", 6, 0.5);
     plot_v_nod(ct,cn=1, soma, Vmin = min(vhold,vstop),     Vmax = max(vhold,vstopx),    colr=brown,    "", 5, 0.5);
     plot_v_nod(ct,cn=2, soma, Vmin = min(vhold,vstop),     Vmax = 0.05,    colr=blue,    "", 4, 0.5);
     plot_i_nod(ct,cn=1, soma, imin = -iscal, imax = iscal, colr=magenta, "", 3, 1); 
     // if (dbpair) plot_i_nod(ct,cn=2, soma, imin = -iscal, imax = iscal, colr=brown, "", 2, 1); 
    // plot_var(&idiff,1,Imin,Imax);  				/* plot current minus cap transient */
       // plot_param ("Idbp_soma", colr=brown, pl=3,plsize=1);
    // plot_var(&current,1,Imin,Imax);  				/* plot current minus cap transient */
     //   plot_param ("Idbp_nonsub", colr=brown, pl=2,plsize=1);
     //plot_var(&cond,1,Gmin,Gmax);                               /* plot cond=current/driving force */
        //plot_param ("Gdbp_soma", colr=green, pl=1,plsize=1);

     // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 
  }

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 0.2;
  simtime = 0 - predur;
  if (predur>0) {
    vclamp            (ndn(ct,cn=1,elnode), vhold, simtime,  predur);
    if (dbpair) vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  predur);
    step (predur);
  }

  dst = 0.0001;		// small time to allow voltage clamp to end

   if (!ivplot) flag = true;
  
  if (vstart < vstop) sign = 1;
  else		      sign = -1;

  // save model so it can be restored after each clamp voltage
  
  set_run_on_exit(runonexit);		     	      // set to erase savefile on ^C
  sprintf (savefile,"cbp_vclamp%06d",getpid());       // add pid to file name
  savemodel (savefile);                        
  
  for (i=0,vpulse=vstart; (vpulse*sign)<=(vstop*sign+1e-6); i++,vpulse += vstep) {

     simtime = 0;
     graph_pen(i+1,i+1,i+1,i+1,i+1);
     vclamp            (ndn(ct,cn=1,elnode), vhold, simtime,  prestepdur);
     if (dbpair) vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  prestepdur);
     step (prestepdur);
     
     if (outward>0) maxCurrent = 0; 		// for outward current
     else           maxCurrent = 1000; 		// for inward current
     gmax = 0;
     vclamp            (ndn(ct,cn=1,elnode), vpulse, simtime,  stepdur);
     if (dbpair) vclamp(ndn(ct,cn=2,elnode), vpulse, simtime,  stepdur);
     step (dst);
     if (ivplot) flag = true;
     step (stepdur-2*dst);
     if (ivplot) graph(voltage, maxCurrent, gmax);
     if (ivplot) flag = false;
     step (dst);

     vclamp            (ndn(ct,cn=1,elnode), tailvolt, simtime,  tailcurdur);
     if (dbpair) vclamp(ndn(ct,cn=2,elnode), tailvolt, simtime,  tailcurdur);
     step (tailcurdur);

     vclamp            (ndn(ct,cn=1,elnode), vhold, simtime,  poststepdur);
     if (dbpair) vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  poststepdur);
     step (poststepdur);

     restoremodel (savefile);                        
  }
  unlink (savefile);                        
  savefile[0] = 0; 
}
