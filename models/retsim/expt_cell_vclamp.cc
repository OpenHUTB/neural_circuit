/* Experiment cell_vclamp */
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

const char *celltype;

double kdr_cond;

double set_drm;

double axon_br_dia;
double varicos_dia;

double predur;
double prestimdur;
double stimdur;
double tailcurdur;
double poststimdur;

double vstart;
double vstop;
double vstep;
double vhold;
double tailvolt;
double dcrm;		// default cell Rm
double gvrev;

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

void defparams(void) 
{
  setptr("celltype",	&celltype);
  setptr("ivplot",	&ivplot);
  setptr("outward",	&outward);
  setptr("elnode",	&elnode);

  setptr("set_drm",     &set_drm);
  setptr("kdr_cond",    &kdr_cond);
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("stimdur",	&stimdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststimdur",	&poststimdur);
  
  setptr("vstart",&vstart);
  setptr("vstop",&vstop);
  setptr("vstep",&vstep);
  setptr("vhold",&vhold);
  setptr("tailvolt",&tailvolt);
  setptr("gvrev",&gvrev);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  setptr("set_timinc", &set_timinc);
  setptr("set_ploti", &set_ploti);
  setptr("elec_rs", &elec_rs);
  setptr("elec_cap", &elec_cap);


  setptr("dcrm",&dcrm);

  nvalfile = "nval_cbp_flash.n";

  // dbp1_file = "morph_bp";		// default set in retsim.cc 
  // dbp1_densfile = "dens_dbp1_full.n";
  // dbp2_densfile = "dens_dbp1_fullx.n";  // zero, no chans

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
  set_ncel(ct,1);

  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_bp" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */
  setn(ct,NCOLOR,RCOLOR);       /* set cell display color from region */

  if (notinit(dispsize)) dispsize = 80;      /* default // display size */
  if (notinit(node_scale)) node_scale = -3.05;/* 3: nodenum, 0.05: small font */

  // if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5; 	/* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5; 	/* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2; 		/* default dia for dendrites */
  

  if (notinit(dcrm)) dcrm = 1e4;  	/* default cell Rm */
  if (notinit(cbplam)) cbplam = 0.005;  	/* default dbp1 complam (dens_dbp1.n) */
  if (notinit(dvrev)) dvrev = -0.05;  	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst)) dvst =   dvrev;  	/* default dbp1 vstart (dens_dbp1.n) */

  if (notinit(cone_maxcond)) cone_maxcond = 1000e-12;  	/* default cone OS cond */

  if (notinit(ax_dia_factor)) ax_dia_factor = 1; /* multiplier for axon diameter */
  if (notinit(dvrev)) dvrev = -0.045;           /* Vrev for dens_dbp1.n */
  if (notinit(dvst))   dvst = -0.07;            /* Vstart for dens_dbp1.n */

  if (notinit(ivplot)) ivplot = 0;              /* make I/V plot */
  if (notinit(outward)) outward = 1;            /* >= 0 => calc K cond, otherwise Na */
  if (notinit(iscal)) iscal   = 1e-9;           /* plot scale */
  if (notinit(gscal)) gscal   = 2e-9;           /* plot scale */

  if (!notinit(set_drm)) {                      /* user set default Rm */
       setn(ct,NRM,set_drm);                    /* set default Rm */
       drm = set_drm;
  }
  vna = 0.04;
  dicafrac = 0;         /* remove Ca pump from ICa */

  // Set channel types in density file.
  // To change, check manual, and uncomment and modify these:  
  //
  //  _NA  = _NA2;
  //  _KA  = _K3;
  //  _KH  = _K4;
  //  _KIR = _K5;
  
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{
   int cn;

    // if (!notinit(kdr_cond))    celdens[dbp1][0][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    ndens[dbp1][cn=1] = 0;		// set cn 1 to use dbp1_densfile
    ndens[dbp1][cn=2] = 1;		// set cn 2 to use dbp1_densfile2
}

/*------------------------------------------------------*/

void runonexit (void)
{
  if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/


bool flag = false;
double maxCurrent;
double cond,Gmax;

void onplot(void) {
    current  = i(ndn(ct, 1, elnode));
    current2 = i(ndn(ct, 2, elnode));
    idiff = current - current2;

    voltage  = v(ndn(ct,  1, soma));

  if (flag) {
   if (outward<=0) {
      cond = idiff / (voltage - gvrev);
      if (cond > Gmax) Gmax = cond;
      if (maxCurrent > idiff) {         // inward current
        maxCurrent = idiff;
      //  fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   } else {                             // outward current 
      cond = idiff / (voltage - vk);
      if (cond > Gmax) Gmax = cond;
      if (maxCurrent < idiff) {
        maxCurrent = idiff;
        //fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
      }
   }
   //  fprintf(stderr, "i: %g, maxi: %g\n", current, maxCurrent);    
  }
  else cond = 0;
}


void runexpt(void)

{
    int cn, i, n, plnum,electrode_node;
    int colr,pl;
    double dst, t, fmax,fmin;
    double rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax;
    double Imin, Imax, Gmin, Gmax;
    double vpulse, sign;


  if (!notinit(set_timinc)) timinc = set_timinc;
  else                      timinc = 5e-6;

   if (!notinit(set_ploti))  ploti = set_ploti;
   else                      ploti = 1e-4;
	     
   // crit = 1e-12;

  cn = 1;
  electrode_node = 5000;

  setonplot(onplot);

  if (notinit(prestimdur))   prestimdur = 0.02;
  if (notinit(stimdur))      stimdur    = 0.1;
  if (notinit(tailcurdur))   tailcurdur = 0.03;
  if (notinit(poststimdur)) poststimdur = 0.05;

  // dtrial = prestimdur+stimdur+tailcurdur+poststimdur;
  dtrial = prestimdur+stimdur+tailcurdur;
  endexp  = dtrial;
  ploti = 1e-4;

  if (notinit(vhold))       vhold  = -0.07;
  if (notinit(vstart))     vstart  = -0.06;
  if (notinit(vstop))       vstop  =  0.02;
  if (notinit(vstep))       vstep  =  0.005;
  if (notinit(tailvolt)) tailvolt  =  vhold;
  if (notinit(gvrev))       gvrev  = vna;

  if (notinit(elec_rs))  elec_rs   = 20e6;
  if (notinit(elec_cap)) elec_cap  = 1e-14;
  if (notinit(elnode))     elnode  = electrode_node;


  // midcbp  = findmid(ct,0,0);

  if (elnode==electrode_node) {
      make_electrode (nd(ct,cn=1,elnode), ndn(ct,cn=1,soma), elec_rs, elec_cap);
      make_electrode (nd(ct,cn=2,elnode), ndn(ct,cn=2,soma), elec_rs, elec_cap);
  }
   
  if (ivplot) {
     graph_x(vstop, vstart);
     graph_y(iscal, -iscal);
     graph_y(gscal, 0);
     graph_init();
  }
  else {

     if      (outward>=1)               { Imax = iscal; Imin = 0;      Gmax=gscal; Gmin=0;}
     else if (outward>0 && outward<1)   { Imax = iscal; Imin = -iscal; Gmax=gscal; Gmin= -gscal;}
     if      (outward==0)               { Imax = 0;     Imin = -iscal; Gmax=gscal; Gmin=0;}
     else if (outward< 0)               { Imax = iscal; Imin = -iscal; Gmax=gscal; Gmin= -gscal;}

     plot_v_nod(ct, cn=1, 12, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 7, 0.5);
     //plot_v_nod(ct, cn=1, 6, Vmin = min(vstart,vstop),     Vmax = //0.05,      colr=cyan,    "", 6, 0.5);
     plot_v_nod(ct, cn=1, soma, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 5, 0.5);
     //plot_v_nod(ct,cn=2, soma, Vmin = -0.1,     Vmax = 0.05, colr=blue,    "", 4, 0.5);
     // plot_i_nod(ct, cn=1, soma, Imin = -iscal, Imax = iscal, colr=magenta, "", 3, 1); 
     //plot_i_nod(ct, cn=2, soma, Imin = -iscal, Imax = iscal, colr=brown, "", 2, 1); 
     plot_var(&idiff,1,Imax,Imin);				/* plot current minus cap transient */
        plot_param ("Idbp_soma", colr=brown, pl=2,plsize=1);
     plot_var(&current,1,Imax,Imin);				/* plot current minus cap transient */
        plot_param ("Idbp_nonsub", colr=brown, pl=2,plsize=1);
     plot_var(&cond,1,Gmax,Gmin);				/* plot cond=current/driving force */
        plot_param ("Gdbp_soma", colr=green, pl=1,plsize=1);

     // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 
  }

  if (notinit(setxmin)) setxmin = 0;                    // set plot to start at 0
  if (notinit(predur)) predur = 0.2;
  simtime = 0 - predur;
  vclamp(ndn(ct,cn=1,elnode), vhold, simtime,  predur);
  vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  predur);
  step (predur);

  dst = 0.0001;         // small time to allow voltage clamp to end

  if (!ivplot) flag = true;
  if (vstart < vstop) sign = 1;
  else                sign = -1;

  // save model so it can be restored after each clamp voltage

  set_run_on_exit(runonexit);                         // set to erase savefile on ^C
  sprintf (savefile,"cell_vclamp%06d",getpid());       // add pid to file name
  savemodel (savefile);

  for (i=0,vpulse=vstart; (vpulse*sign)<=(vstop*sign); i++,vpulse += vstep) {

     simtime = 0;
     graph_pen(i+1,i+1,i+1,i+1,i+1);
     vclamp(ndn(ct,cn=1,elnode), vhold, simtime,  prestimdur);
     vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  prestimdur);
     step (prestimdur);

     if (outward>0) maxCurrent = 0;		// for outward current
     else           maxCurrent = 1000;		// for inward current
     Gmax = 0;
     vclamp(ndn(ct,cn=1,elnode), vpulse, simtime,  stimdur);
     vclamp(ndn(ct,cn=2,elnode), vpulse, simtime,  stimdur);
     step (dst);
     flag = true;
     step (stimdur-2*dst);
     if (ivplot) graph(voltage, maxCurrent, Gmax);
     flag = false;
     step (dst);

     vclamp(ndn(ct,cn=1,elnode), tailvolt, simtime,  tailcurdur);
     vclamp(ndn(ct,cn=2,elnode), tailvolt, simtime,  tailcurdur);
     step (tailcurdur);

     vclamp(ndn(ct,cn=1,elnode), vhold, simtime,  poststimdur);
     vclamp(ndn(ct,cn=2,elnode), vhold, simtime,  poststimdur);
     step (poststimdur);

     restoremodel (savefile);
  }
  unlink (savefile);
  savefile[0] = 0;
}
