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
#include "stimfuncs.h"

int ct;
int ivplot;
int outward;
int sbacpair;
int stimtype;
int sbarr;
int no_inhib;
int no_excit;
int direction;
int waveshape;
int sbac_color;

const char *celltype;

double kdr_cond;

double set_drm;
double nadist;
double namid;
double kdist;
double kmid;
double kprox;
double ksoma;
double cadist;
double camid;
double set_excit;
double set_inhib;
double sbspac;
double sbac_synrng;
double sbac_synanp;
double sbac_isynrngi;
double sbac_isynanpi;
double sbac_isynanpo;

double axon_br_dia;
double varicos_dia;
double sbaclm;

double predur;
double prestimdur;
double stimdur;
double tailcurdur;
double poststimdur;
double spotdur;
double ncycles;
double minten;
double scontrast;
double stimtime;
double dstim;
double sdia;

double barwidth;
double barlength;
double theta;
double velocity;

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
int set_vclamp;
int set_tonic;

void sb_init(void);

double sb_vr;
double sb_vs;
double sb_rm;
double sb_rid;			// distal Ri
double sb_rii;			// intermediate Ri
double sb_rip;			// proximal Ri

char savefile[30] = {0};


/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("celltype",    &celltype);
  setptr("ivplot",      &ivplot);
  setptr("outward",     &outward);
  setptr("sbacpair",    &sbacpair);
  setptr("stimtype",    &stimtype);
  setptr("sbarr",       &sbarr);
  setptr("no_inhib",    &no_inhib);
  setptr("no_excit",    &no_excit);
  setptr("set_excit",   &set_excit);
  setptr("set_inhib",   &set_inhib);
  setptr("direction",   &direction);
  setptr("waveshape",   &waveshape);
  setptr("sbspac",      &sbspac);
  setptr("sbac_synrng", &sbac_synrng);
  setptr("sbac_synanp", &sbac_synanp);
  setptr("sbac_isynrngi", &sbac_isynrngi);
  setptr("sbac_isynanpi", &sbac_isynanpi);
  setptr("sbac_isynanpo", &sbac_isynanpo);
  setptr("sbac_color", &sbac_color);

  setptr("set_drm",     &set_drm);
  setptr("kdr_cond",    &kdr_cond);
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);
  setptr("sbaclm",      &sbaclm);
  setptr("nadist",      &nadist);
  setptr("namid",       &namid);
  setptr("kdist",       &kdist);
  setptr("kmid",        &kmid);
  setptr("kprox",       &kprox);
  setptr("ksoma",       &ksoma);
  setptr("cadist",      &cadist);
  setptr("camid",       &camid);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("stimdur",	&stimdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststimdur",	&poststimdur);
  setptr("spotdur",	&spotdur);
  setptr("ncycles",	&ncycles);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("stimtime",  &stimtime);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  
  setptr("vstart",&vstart);
  setptr("vstop",&vstop);
  setptr("vstep",&vstep);
  setptr("vhold",&vhold);
  setptr("tailvolt",&tailvolt);
  setptr("gvrev",&gvrev);
  setptr("sb_vs",&sb_vs);
  setptr("sb_vr",&sb_vr);
  setptr("sb_rm",&sb_rm);

  setptr("sb_rid",&sb_rid);
  setptr("sb_rii",&sb_rii);
  setptr("sb_rip",&sb_rip);
	
  setptr("velocity",&velocity);
  setptr("barwidth",&barwidth);
  setptr("barlength",&barlength);
  setptr("theta",&theta);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  
  setptr("set_vclamp", &set_vclamp);
  setptr("set_tonic", &set_tonic);

  setptr("dcrm",&dcrm);

  nvalfile = "nval_dsgc_sbac2.n";

  // for excitatory interconnections:
  //
  if (notinit(sbac_synrng)) { sbac_synrng = 0; }  // for nval file
  if (notinit(sbac_synanp)) { sbac_synanp = 0; }  // for nval file

  // for inhibitory interconnections:
  //
  if (notinit(sbac_isynrngi)) { sbac_isynrngi = 0; }  // (for nval file)
  if (notinit(sbac_isynanpi)) { sbac_isynanpi = 0; }  // inner radius of annulus in postsyn cell (for nval file)
  if (notinit(sbac_isynanpo)) { sbac_isynanpo = 0; }  // outer radius of annulus in postsyn cell (for nval file)

  make_sbac_sbac = 1;

  chanparamsfile = "chanparams2";
  sb_init();

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
     else              color = CYAN;
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
  sprintf (tbuf,"%d %d",spnt->node1b,spnt->node2b);	/* print pre and postsynaptic cell number */
  tlen = strlen(tbuf);
  gmove (length/2.0 -tlen*0.3, -1.0);
  gcwidth (1.2);
  gtext (tbuf);
}


/*------------------------------------------------------*/
   
void setparams(void)
{
   int i;

  if (!notinit(celltype)) {
       ct = find_ct(celltype);
  } else {
       ct  = sbac;
  }
  make_ct(ct);          /* make the cell type */
  // if (!notinit(set_vclamp) && set_vclamp>0) set_ncel(ct,2);
  set_ncel(ct,1);

#define SBARR 30

 if (!notinit(sbarr)) {
   sbxarr = (double *)emalloc(SBARR*sizeof(double));
   sbyarr = (double *)emalloc(SBARR*sizeof(double));
   sbtharr = (double *)emalloc(SBARR*sizeof(double));
   sbnarr = (int *)emalloc(SBARR*sizeof(int));

   for (i=0; i<SBARR; i++) sbnarr[i] = i+1;

   // if (notinit (sbtheta)) sbtheta = 0;
   if (notinit (sbspac)) sbspac = 50;

   if (strcmp(sbac_file,"morph_sbac3c")==0) {
     if (sbarr==0) {
           sbxarr[0] =  0; sbyarr[0] = 0; sbtharr[0] = 0;  /* only 1 SBAC */
           n_sbac = 1;
     }
     if (sbarr==1) {                    /* sbac3c morphology, opposing dendrites */
        sbxarr[0] =  85; sbyarr[0] = 0; sbtharr[0] = 0;
        sbxarr[1] = -85; sbyarr[1] = 0; sbtharr[1] = 0;
        n_sbac = 2;
     }
     if (sbarr==2) {                    /* sbac3c morphology, aligned dendrites */
         sbxarr[0] =   0; sbyarr[0] = 0; sbtharr[0] = 0;
         sbxarr[1] =   1; sbyarr[1] = 0; sbtharr[1] = 30;
         n_sbac = 2;
     }
     if (sbarr==3) {                    /* sbac3c morphology, 3 aligned */
        sbxarr[0] =    0; sbyarr[0] = 0; sbtharr[0] = 0; sbnarr[0]=1;
        sbxarr[1] =  -sbspac; sbyarr[1] = 0; sbtharr[1] = 30; sbnarr[1]=3;
        sbxarr[2] =   sbspac; sbyarr[2] = 0; sbtharr[2] = 60; sbnarr[2]=13;
        n_sbac = 3;
     }
     if (sbarr==4) {                    /* sbac3c morphology, 3 aligned, 1 opposing */
        sbxarr[0] =   sbspac; sbyarr[0] = 0; sbtharr[0] = 0;
        sbxarr[1] =   sbspac; sbyarr[1] = 0; sbtharr[1] = 0;
        sbxarr[2] =        0; sbyarr[2] = 0; sbtharr[2] = 0;
        sbxarr[3] =  -sbspac; sbyarr[3] = 0; sbtharr[3] = 0;
        n_sbac = 4;
     }
     if (sbarr==7) {                    /* sbac3c morphology, 7 overlapping */
        sbxarr[0] =  3*sbspac; sbyarr[0] = 0; sbtharr[0] = 0;
        sbxarr[1] =  2*sbspac; sbyarr[1] = 0; sbtharr[1] = 0;
        sbxarr[2] =    sbspac; sbyarr[2] = 0; sbtharr[2] = 0;
        sbxarr[3] =         0; sbyarr[3] = 0; sbtharr[3] = 0;
        sbxarr[4] =   -sbspac; sbyarr[4] = 0; sbtharr[4] = 0;
        sbxarr[5] = -2*sbspac; sbyarr[5] = 0; sbtharr[5] = 0;
        sbxarr[6] = -3*sbspac; sbyarr[6] = 0; sbtharr[6] = 0;
        n_sbac = 7;
      }
     if (sbarr==102) {
       sbxarr[0] =  145; sbyarr[0] =  50; sbtharr[0] = 0;
       sbxarr[1] =  145; sbyarr[1] = -50; sbtharr[1] = 0;
       n_sbac = 2;
     }
     if (sbarr==103) {			// 5 sbacs
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0;
       sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;
       sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;
       sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0;
       sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90;
       n_sbac = 5;
     }
     if (sbarr==104) {			// 9 sbacs
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0;
       sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;
       sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;
       sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0;
       sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90;
       sbxarr[5] = sbspac;  sbyarr[5] =  sbspac; sbtharr[5] =  0;
       sbxarr[6] = -sbspac; sbyarr[6] =  sbspac; sbtharr[6] = 120;
       sbxarr[7] = sbspac;  sbyarr[7] =  -sbspac; sbtharr[7] = 180;
       sbxarr[8] = -sbspac; sbyarr[8] =  -sbspac; sbtharr[8] = 210;
       n_sbac = 9;
     }
     if (sbarr==105) {			// 13 sbacs
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0;
       sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;
       sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;
       sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0;
       sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90;

       sbxarr[5] = sbspac;  sbyarr[5] =  sbspac; sbtharr[5] =  0;
       sbxarr[6] = -sbspac; sbyarr[6] =  sbspac; sbtharr[6] = 120;
       sbxarr[7] = sbspac;  sbyarr[7] =  -sbspac; sbtharr[7] = 180;
       sbxarr[8] = -sbspac; sbyarr[8] =  -sbspac; sbtharr[8] = 210;

       sbxarr[9] = -2*sbspac; sbyarr[9] =  0;         sbtharr[9] = 120;
       sbxarr[10] = 2*sbspac; sbyarr[10] = 0;         sbtharr[10] = 90;
       sbxarr[11] = 0;        sbyarr[11] = -2*sbspac; sbtharr[11] = 60;
       sbxarr[12] = 0;        sbyarr[12] =  2*sbspac; sbtharr[12] = 30;
       n_sbac = 13;
     }
     if (sbarr==1052) {			// 2 sbacs: 1, 2
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0;   sbnarr[0] = 1;
       sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;  sbnarr[1] = 2;
       //sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;sbnarr[2] = 3;
       //sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0; sbnarr[3] = 4;
       //sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90;sbnarr[4] = 5;

       //sbxarr[5] = sbspac;  sbyarr[5] =  sbspac; sbtharr[5] =  0;   sbnarr[5] = 6;
       //sbxarr[1] = -sbspac; sbyarr[1] =  sbspac; sbtharr[1] = 120;  sbnarr[6] = 7;
       //sbxarr[7] = sbspac;  sbyarr[7] =  -sbspac; sbtharr[7] = 180; sbnarr[7] = 8;
       //sbxarr[8] = -sbspac; sbyarr[8] =  -sbspac; sbtharr[8] = 210; sbnarr[8] = 9;

       //sbxarr[9] = -2*sbspac; sbyarr[9] =  0;         sbtharr[9] = 120; sbnarr[9] = 10;
       //sbxarr[10] = 2*sbspac; sbyarr[10] = 0;         sbtharr[10] = 90; sbnarr[10] = 11;
       //sbxarr[11] = 0;        sbyarr[11] = -2*sbspac; sbtharr[11] = 60; sbnarr[11] = 12;
       //sbxarr[12] = 0;        sbyarr[12] =  2*sbspac; sbtharr[12] = 30; sbnarr[12] = 13;
       n_sbac = 2;
     }
     if (sbarr==1053) {			// 3 sbacs: 1, 2, 3
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0;    sbnarr[0] = 1;
       sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;   sbnarr[1] = 2;
       sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;   sbnarr[2] = 3;
       //sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0;  sbnarr[2] = 4;
       //sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90; sbnarr[4] = 5;

       //sbxarr[5] = sbspac;  sbyarr[5] =  sbspac; sbtharr[5] =  0; sbnarr[5] = 6;
       //sbxarr[1] = -sbspac; sbyarr[1] =  sbspac; sbtharr[1] = 120; sbnarr[6] = 7;
       //sbxarr[7] = sbspac;  sbyarr[7] =  -sbspac; sbtharr[7] = 180; sbnarr[7] = 8;
       //sbxarr[8] = -sbspac; sbyarr[8] =  -sbspac; sbtharr[8] = 210; sbnarr[8] = 9;

       //sbxarr[9] = -2*sbspac; sbyarr[9] =  0;         sbtharr[9] = 120; sbnarr[9] = 10;
       //sbxarr[10] = 2*sbspac; sbyarr[10] = 0;         sbtharr[10] = 90; sbnarr[10] = 11;
       //sbxarr[11] = 0;        sbyarr[11] = -2*sbspac; sbtharr[11] = 60; sbnarr[11] = 12;
       //sbxarr[12] = 0;        sbyarr[12] =  2*sbspac; sbtharr[12] = 30; sbnarr[12] = 13;
       n_sbac = 3;
     }
     if (sbarr==1057) {			// 2 sbacs: 1, 7
       sbxarr[0] =  0;      sbyarr[0] = 0;       sbtharr[0] = 0; sbnarr[0] = 1;
       //sbxarr[1] = -sbspac; sbyarr[1] = 0;       sbtharr[1] = 30;
       //sbxarr[2] =  sbspac; sbyarr[2] = 0;       sbtharr[2] = 60;
       //sbxarr[3] = 0;       sbyarr[3] = -sbspac; sbtharr[3] = 0;
       //sbxarr[4] = 0;       sbyarr[4] =  sbspac; sbtharr[4] = 90;

       //sbxarr[5] = sbspac;  sbyarr[5] =  sbspac; sbtharr[5] =  0;
       sbxarr[1] = -sbspac; sbyarr[1] =  sbspac; sbtharr[1] = 120; sbnarr[1] = 27;
       //sbxarr[7] = sbspac;  sbyarr[7] =  -sbspac; sbtharr[7] = 180;
       //sbxarr[8] = -sbspac; sbyarr[8] =  -sbspac; sbtharr[8] = 210;

       //sbxarr[9] = -2*sbspac; sbyarr[9] =  0;         sbtharr[9] = 120;
       //sbxarr[10] = 2*sbspac; sbyarr[10] = 0;         sbtharr[10] = 90;
       //sbxarr[11] = 0;        sbyarr[11] = -2*sbspac; sbtharr[11] = 60;
       //sbxarr[12] = 0;        sbyarr[12] =  2*sbspac; sbtharr[12] = 30;
       n_sbac = 2;
     }
   }
 }

 if (notinit(sbacpair)) sbacpair = 0;
 make_dbp1 = 1;

//  dbp2_file = dbp1_file;        /* Use the same morphology for dbp1 and dbp2 */

  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_xx" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */

//  setn(ct2,MORPH,0);            /* set cell morphology from file, default = "morph_bp" */
//  setn(ct2,BIOPHYS,1);          /* set cell biophys from file, default = "dens_default" */


  setn(ct,NCOLOR,RCOLOR);       /* set cell display color from region */
if (!notinit(sbac_color)) setn(ct,NCOLOR,sbac_color);   /* set cell display color from region */

  set_synapse_dr (syn_draw2);

  DENDD     = R_1;
  DEND_DIST = R_1;
  DEND      = R_2;
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

  if (notinit(dispsize)) dispsize = 150;      /* default // display size */
  if (notinit(node_scale)) node_scale = -3.05;/* 3: nodenum, 0.05: small font */

  // if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5; 	/* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5; 	/* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2; 		/* default dia for dendrites */
  

  if (notinit(dcrm)) dcrm = 1e4;  	/* default cell Rm */
  if (notinit(sbaclm)) sbaclm = 0.1;  	/* default sbac complam (dens_sbaca.n) */
  if (notinit(dvrev)) dvrev = -0.05;  	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst)) dvst =   dvrev;  	/* default dbp1 vstart (dens_dbp1.n) */

  if (notinit(cone_maxcond)) cone_maxcond = 1000e-12;  	/* default cone OS cond */

  if (!notinit(set_excit)) { setsv(sbac,SCOND,3,set_excit); }
  if (!notinit(no_excit)) { if (no_excit==1) { setsv(sbac,CELPRE,3,-1); setsv(sbac,CELPOST,3,-1); }}

  if (!notinit(set_inhib)) { setsv(sbac,SCOND,4,set_inhib); }
  if (!notinit(no_inhib)) { if (no_inhib==1) { setsv(sbac,CELPRE,4,-1); setsv(sbac,CELPOST,4,-1); }}

  if (notinit(ax_dia_factor)) ax_dia_factor = 1; /* multiplier for axon diameter */
  if (notinit(dvrev)) dvrev = -0.045;           /* Vrev for dens_dbp1.n */
  if (notinit(dvst))   dvst = -0.07;            /* Vstart for dens_dbp1.n */

  if (notinit(ivplot)) ivplot = 0;              /* make I/V plot */
  if (notinit(outward)) outward = 0;            /* >0 => calc K cond, otherwise Na */
  if (notinit(iscal)) iscal   = 3e-9;           /* plot scale */
  if (notinit(gscal)) gscal   = 1.5e-8;         /* plot scale */
  if (notinit(sb_rm)) sb_rm   = 10e-3;          /* sbac RM, used in dens_sbaca.n */
  if (notinit(sb_vs)) sb_vs   = -0.07;          /* sbac vstart, used in dens_sbaca.n */
  if (notinit(sb_vr)) sb_vr   = -0.07;          /* sbac vrev,   used in dens_sbaca.n */

  if (notinit(nadist)) nadist   = 10e-3;        /* distal Na density, dens_sbaca.n */
  if (notinit(namid))   namid   = 5e-3;         /* middle Na density, dens_sbaca.n */
  if (notinit(kdist))   kdist   = 3e-3;         /* dist Kdr density, dens_sbaca.n */
  if (notinit(kmid))     kmid   = 2e-3;         /* mid  Kdr density, dens_sbaca.n */
  if (notinit(kprox))   kprox   = 2e-3;         /* prox Kdr density, dens_sbaca.n */
  if (notinit(ksoma))   ksoma   = 2e-3;         /* soma Kdr density, dens_sbaca.n */
  if (notinit(cadist)) cadist   = 3e-3;         /* dist Ca density, dens_sbaca2.n */
  if (notinit(camid))   camid   = 0e-3;         /* mid  Ca density, dens_sbaca2.n */

  if (!notinit(set_drm)) {                      /* user set default Rm */
       setn(ct,NRM,set_drm);                    /* set default Rm */
       drm = set_drm;
  }
  if (notinit(sb_rid)) sb_rid = dri;
  if (notinit(sb_rii)) sb_rii = dri;
  if (notinit(sb_rip)) sb_rip = dri;

  vna = 0.04;
  dicafrac = 0;         /* remove Ca pump from ICa */

  if (notinit(set_tonic)) set_tonic = 0;
  if (set_tonic >= 1) setsv (dbp1,SNFILTH,3, 0);

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
    // if (!notinit(kdr_cond))    celdens[dbp1][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    
    // int cn;
    // ndens[sbac][cn=1] = 0;           // set cn 1 to use sbac_densfile
    // ndens[sbac][cn=2] = 1;  		// set cn 2 to use sbac_densfile2
    // ndens[sbac][cn=3] = 0;  		// set cn 3 to use sbac_densfile
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
    current  = i(ndn(ct,  1, soma));
    if (sbacpair) current2 = i(ndn(ct, 2, soma));
    else          current2 = 0;
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

/*------------------------------------------------------*/

double mvbrt1;
double mvbrt2;

double move_stim(double stimtime,double barwidth,double theta,double velocity, double ncycles, double scontrast, int direction, double mask)
{
    static int runyet=0;
    double lbar, rbar, irad, orad;
    double inten, start, dur, wavel;
    double cellrad;
    double s,send;
    double tfreq;

 cellrad = max(xarrsiz,yarrsiz) * 0.525; 

 // stim_spot(100, 100,0, inten=minten,  start=0.02, dur=0.05, wavel=1, 0);

 if (stimtype==1) {			// move bar to & fro
   lbar = -barwidth / 2;
   rbar = cellrad + barwidth/2;
   mvbrt1 = movebar (stimtime,            0,0, lbar, rbar, barwidth,barlength,theta,velocity,scontrast);
   mvbrt2 = movebar (mvbrt1+stimtime+0.05,0,0, rbar, lbar, barwidth,barlength,theta,velocity,scontrast);
   if (runyet==0 && mask>0) {
       stim_bar(200,1000, -100,0, theta, inten=minten, start=0, dur=1, wavel=1, mask);
       runyet = 1;
   }
 }
 else if (stimtype==2) {		// move annulus to & fro

   irad = 20;
   orad = cellrad + barwidth;
   mvbrt1 = moveannulus (stimtime,            0,0, irad, orad, barwidth,velocity,scontrast);
   mvbrt2 = moveannulus (mvbrt1+stimtime+0.05,0,0, orad, irad, barwidth,velocity,scontrast);
   if (runyet==0 && mask>0) {
  //     stim_spot(idia, 0,0, inten=minten,  start=0, dur=1, wavel=1, mask);
       stim_bar(200,1000, -100,0, theta, inten=minten, start=0, dur=1, wavel=1, mask);
       runyet = 1;
   }
 } else if (stimtype==3) {		// move sineann

  //void movesineann (x,y,direction,ann_gaussenv,centdia, phase,speriod,stfreq,0,1.0,sinten,contrast, sq, stimtime,sdur)
  //
  tfreq = velocity/(2*barwidth);
  dur = ncycles / tfreq;
  if (direction==1) 
      movesineann (0, 0, direction ,1e6,0, 180, barwidth*2, tfreq, 0, 1.0,  minten, scontrast, waveshape, stimtime,dur);
  else
      movesineann (0, 0, direction ,1e6,0, 250, barwidth*2, tfreq, 0, 1.0, minten, scontrast, waveshape, stimtime,dur);
  mvbrt1 = mvbrt2 = dur;
 }
 else if (stimtype==4) {		// move spot
 	  send = int(200/sdia + 0.5);
	  for (s=0; s<send; s++) {
            stim_spot(sdia, s*sdia, 0, scontrast, stimtime+s*spotdur, spotdur);
          }
	  mvbrt2 = spotdur * send;
       }

   return mvbrt2;
}

/*------------------------------------------------------*/

void addlabels(void)

{
        //label (ndn(sbac,2,96), red);
        //label (ndn(sbac,3,2726), red);
	label (findsynloc(sbac,2,sbac,1,100, 0, 0.0, 50), red);
	label (findsynloc(sbac,3,sbac,1,100, 0, 0.0, 50), red);
}

/*------------------------------------------------------*/

void runexpt(void)

{
    int cn, i, n, s, plnum;
    int colr,pl,nsynap, nsynape, nsynapi;
    double dst, t, fmax,fmin;
    double cmin, cmax, rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax, Vmino;
    double Imin, Imax;
    double imin, imax, gmax;
    double vpulse, sign, pulsedur;
    node *npnt;
    photorec *p;
    double disp_end, starttime;
    double dscale, mask;

  timinc = 2e-6;
  ploti = 1e-4;
  crit = 1e-12;

  for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
        p = (photorec*)make_transducer(ndn(dbp1,cn,soma));
        p->xpos=npnt->xloc;
        p->ypos=npnt->yloc;
	//	fprintf(stderr, "%g, %g\n", npnt->xloc, npnt->yloc);
   }
  cn = 1;
  setonplot(onplot);

  if (notinit(prestimdur))   prestimdur = 0.05;
  if (notinit(stimdur))      stimdur    = 0.1;
  if (notinit(tailcurdur))   tailcurdur = 0.05;
  if (notinit(poststimdur)) poststimdur = 0.00;
  if (notinit(spotdur))      spotdur    = 0.1;
  if (notinit(minten))       minten = -0.050;
  if (notinit(scontrast)) scontrast = 0.005; 
  if (notinit(stimtime))   stimtime = .1;
  if (notinit(dstim))         dstim = .05;  
  if (notinit(sdia))           sdia = 30;      
  if (notinit(ncycles))       ncycles = 2;      
  if (notinit(stimtype))   stimtype = 1;      

  // dtrial = prestimdur+stimdur+tailcurdur+poststimdur;
//  dtrial = prestimdur+dstim+poststimdur+0.1;
//  endexp  = dtrial;
  ploti = 1e-4;

 
  if (notinit(vhold))       vhold  = -0.060;
  if (notinit(vstart))     vstart  = -0.120;
  if (notinit(vstop))       vstop  = -0.00;
  if (notinit(vstep))       vstep  =  0.03;
  if (notinit(tailvolt)) tailvolt  = vhold;
  if (notinit(gvrev))       gvrev  = vna;

  // midcbp  = findmid(ct,0,0);

  nsynap = synapse_add (1,dbp1,-1,-1,sbac,1);		/* make list of bipolar synapses onto sbac 1 for recording below */
  nsynape = nsynapi = 0;
  for (i=2; i<=n_sbac; i++) {
    nsynape = synapse_add (3,sbac,i,-1,sbac,1,3);	/* make list of excit synapses from sbacs onto sbac1 */
    nsynap += nsynape;
    fprintf (stderr,"# nsynape %d\n",nsynape);
  }
  for (i=2; i<=n_sbac; i++) {
    nsynapi = synapse_add (4,sbac,i,-1,sbac,1,4);	/* make list of inhib synapses from sbacs onto sbac1 */
    nsynap += nsynapi;
    fprintf (stderr,"# nsynapi %d\n",nsynapi);
  }
  fprintf (stderr,"# nsynap %d\n",nsynap);

  if (ivplot) {
     graph_x(vstop, vstart);
     graph_y(iscal, -iscal);
     graph_y(gscal, 0);
     graph_init();
  }
  else {


        plot_v_nod(dbp1,findmida(dbp1,0,0),  soma,Vmin=-0.05, Vmax= -0.03, 1, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,20,0), soma,Vmin=-0.05, Vmax= -0.03, 2, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,40,0), soma,Vmin=-0.05, Vmax= -0.03, 3, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,60,0), soma,Vmin=-0.05, Vmax= -0.03, 4, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,80,0), soma,Vmin=-0.05, Vmax= -0.03, 5, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,100,0),soma,Vmin=-0.05, Vmax= -0.03, 6, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,120,0),soma,Vmin=-0.05, Vmax= -0.03, 7, "", 22, 0.3);

/*
	plot_synrate(findsynloc(dbp1,0,0),   cmin=0,cmax=500, 1, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,20,0),  cmin=0,cmax=500, 2, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,40,0),  cmin=0,cmax=500, 3, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,60,0),  cmin=0,cmax=500, 4, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,80,0),  cmin=0,cmax=500, 5, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,100,0), cmin=0,cmax=500, 6, 21,"",0.3); 
  	plot_synrate(findsynloc(dbp1,120,0), cmin=0,cmax=500, 7, 21,"",0.3); 
*/

	plot_syncond(findsynloc(dbp1,0,0),   cmin=0,cmax=100e-12, 1, 20,"",0.3); 
	plot_syncond(findsynloc(dbp1,20,0),  cmin=0,cmax=100e-12, 2, 20,"",0.3); 
	plot_syncond(findsynloc(dbp1,40,0),  cmin=0,cmax=100e-12, 3, 20,"",0.3); 
	plot_syncond(findsynloc(dbp1,60,0),  cmin=0,cmax=100e-12, 4, 20,"",0.3); 
	plot_syncond(findsynloc(dbp1,80,0),  cmin=0,cmax=100e-12, 5, 20,"",0.3); 
	plot_syncond(findsynloc(dbp1,100,0), cmin=0,cmax=100e-12, 6, 20,"",0.3); 
  	plot_syncond(findsynloc(dbp1,120,0), cmin=0,cmax=100e-12, 7, 20,"",0.3); 

	// plot_syncond(findsynlocr(sbac,1,80, 0, 0.0), cmin=0,cmax=100e-12, 3, 18,"",0.3); 
	// plot_syncond(findsynlocr(sbac,1,120,0, 0.0), cmin=0,cmax=100e-12, 6, 18,"",0.3); 

	plot_syncond(findsynlocr(sbac,2,80, 0, 0.0), cmin=0,cmax=100e-12, 3, 18,"",0.5); 
	plot_syncond(findsynlocr(sbac,2,120,0, 0.0), cmin=0,cmax=100e-12, 6, 18,"",0.5); 
	
	plot_syncond(findsynloc(sbac,2,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 1, 17,"",0.5); 
	plot_syncond(findsynloc(sbac,3,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 2, 17,"",0.5); 
	plot_syncond(findsynloc(sbac,4,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 3, 17,"",0.5); 
	plot_syncond(findsynloc(sbac,5,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 4, 17,"",0.5); 
	plot_syncond(findsynloc(sbac,6,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 5, 17,"",0.5); 
	plot_syncond(findsynloc(sbac,7,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 6, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,8,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 7, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,9,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 8, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,10,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 9, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,11,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 10, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,12,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 11, 17,"",0.5); 
  	plot_syncond(findsynloc(sbac,13,sbac,1,100, 0, 0.0, 50), cmin=0,cmax=100e-12, 12, 17,"",0.5); 

	//plot_synrate(findsynloc(sbac,2,0,0),   cmin=0,cmax=500, 1, 16,"",0.3); 
	//plot_synrate(findsynloc(sbac,2,20,0),  cmin=0,cmax=500, 2, 16,"",0.3); 
	//plot_synrate(findsynloc(sbac,2,40,0),  cmin=0,cmax=500, 3, 16,"",0.3); 
	//plot_synrate(findsynloc(sbac,2,60,0),  cmin=0,cmax=500, 4, 16,"",0.3); 
	//plot_synrate(findsynloc(sbac,2,80,0),  cmin=0,cmax=500, 5, 16,"",0.3); 
	//plot_synrate(findsynloc(sbac,2,100,0), cmin=0,cmax=500, 6, 16,"",0.3); 
  	//plot_synrate(findsynloc(sbac,2,120,0), cmin=0,cmax=500, 7, 16,"",0.3); 

	// plot_ca_nod(ct, 1, 304, 		    20e-6, colr=blue, "", 15, 0.5);
	plot_ca_nod(ct, 3, 2726, 		    20e-6, colr=blue, "", 15, 0.5);
	plot_ca_nod(ct, 3, 2727, 		    20e-6, colr=blue, "", 15, 0.5);
	plot_ca_nod(ct, 3, 2728, 		    20e-6, colr=blue, "", 15, 0.5);
	plot_ca_nod(ct, 3, 2729, 		    20e-6, colr=blue, "", 15, 0.5);
	plot_ca_nod(ct, 3, 2730, 		    20e-6, colr=blue, "", 15, 0.5);
	// plot_ca_nod(ct, 3, findnodloc(ct,3,120,0,20), 10e-6, colr=cyan, "", 14, 0.5);
	// plot_ca_nod(ct, 3, findnodloc(ct,3,100,0,20), 10e-6, colr=brown, "", 14, 0.5);
        // vclamp(ndn(ct,1,304), -0.06, 0.050,  1);
        //vclamp(ndn(ct,3,2735), -0.06, 0.130,  1);

	plot_syncond(findsynloc(sbac,2,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 1, 14,"Gi_2",0.5); 
	plot_syncond(findsynloc(sbac,3,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 2, 14,"Gi_3",0.5); 
	plot_syncond(findsynloc(sbac,4,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 3, 14,"Gi_4",0.5); 
	plot_syncond(findsynloc(sbac,5,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 4, 14,"Gi_5",0.5); 
	plot_syncond(findsynloc(sbac,6,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 5, 14,"Gi_6",0.5); 
	plot_syncond(findsynloc(sbac,7,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 6, 14,"Gi_7",0.5); 
  	plot_syncond(findsynloc(sbac,8,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 7, 14,"Gi_8",0.5); 
  	plot_syncond(findsynloc(sbac,9,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 8, 14,"Gi_9",0.5); 
  	plot_syncond(findsynloc(sbac,10,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 9, 14,"Gi_10",0.5); 
  	plot_syncond(findsynloc(sbac,11,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 10, 14,"Gi_11",0.5); 
  	plot_syncond(findsynloc(sbac,12,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 11, 14,"Gi_12",0.5); 
  	plot_syncond(findsynloc(sbac,13,sbac,1,60, 0, -0.06, 50), cmin=0,cmax=100e-12, 12, 14,"Gi_13",0.5); 


	if (set_vclamp > 0) { Vmin=vstart; if (stimtype==4) Vmax=vstart+0.01; else Vmax = 0; }
	else                { Vmin = -0.07; Vmax = -0.00; }

	Vmino = -0.06;
	if (set_vclamp > 0) {
	  // plot_v_nod(ct, 2,  soma,		       Vmino, Vmax, colr=blue,   "", 11, 0.5);
	  // plot_v_nod(ct, 2,  findnodlocr(ct,2,30,0),   Vmino, Vmax, colr=green,"", 11, 0.5);
	  // plot_v_nod(ct, 2,  findnodlocr(ct,2,80,0),   Vmino, Vmax, colr=brown,"", 11, 0.5);
	  // plot_v_nod(ct, 2,  findnodlocr(ct,2,120,0),  Vmino, Vmax, colr=red,  "", 11, 0.5);

	  // plot_v_nod(ct, 3,  findnodlocr(ct,3,120,0),  Vmino, Vmax, colr=12,"", 11, 0.5);
	  // plot_v_nod(ct, 3,  findnodlocr(ct,3,-120,0), Vmino, Vmax, colr=13,"", 11, 0.5);
	  // plot_v_nod(ct, 4,  findnodlocr(ct,3,0,-120), Vmino, Vmax, colr=14,"", 11, 0.5);
	  // plot_v_nod(ct, 5,  findnodlocr(ct,3,0,120),  Vmino, Vmax, colr=15,"", 11, 0.5);
	  plot_v_nod(ct, 2,  findnodloc(ct,2,50,0,20),      Vmino, Vmax, colr=1,   "", 11, 0.5);
	  plot_v_nod(ct, 3,  findnodloc(ct,3,100,0,20),      Vmino, Vmax, colr=2,  "", 11, 0.5);
	  plot_v_nod(ct, 4,  findnodloc(ct,4,100,0,20),      Vmino, Vmax, colr=3,   "", 11, 0.5);
	  plot_v_nod(ct, 5,  findnodloc(ct,5,100,0,20),      Vmino, Vmax, colr=4,    "", 11, 0.5);
	  plot_v_nod(ct, 6,  findnodloc(ct,6,100,0,20),      Vmino, Vmax, colr=5,"", 11, 0.5);
	  plot_v_nod(ct, 7,  findnodloc(ct,7,100,0,20),      Vmino, Vmax, colr=6,  "", 11, 0.5);
	  plot_v_nod(ct, 8,  findnodloc(ct,8,100,0,20),      Vmino, Vmax, colr=7,  "", 11, 0.5);
	  plot_v_nod(ct, 9,  findnodloc(ct,9,100,0,20),      Vmino, Vmax, colr=8, "", 11, 0.5);
	  plot_v_nod(ct, 10,  findnodloc(ct,10,100,0,20),      Vmino, Vmax, colr=9,"", 11, 0.5);
	  plot_v_nod(ct, 11,  findnodloc(ct,11,100,0,20),      Vmino, Vmax, colr=10,"", 11, 0.5);
	  plot_v_nod(ct, 12,  findnodloc(ct,12,100,0,20),      Vmino, Vmax, colr=11, "", 11, 0.5);
	  plot_v_nod(ct, 13,  findnodloc(ct,13,100,0,20),      Vmino, Vmax, colr=12, "", 11, 0.5);
	}
	else {
	  plot_v_nod(ct, 2,  findnodlocr(ct,2,100,0),  Vmino, Vmax, colr=green,"", 12, 0.5);
	  plot_v_nod(ct, 3,  findnodlocr(ct,3,100,0),  Vmino, Vmax, colr=magenta,"", 12, 0.5);
	}

	plot_v_nod(ct, cn, soma, 		    Vmin, Vmax, colr=blue, "", 12, 0.5);
	plot_v_nod(ct, cn, findnodloc(ct,cn,30,0),  Vmin, Vmax, colr=green,"", 12, 0.5);
	plot_v_nod(ct, cn, findnodloc(ct,cn,80,0),  Vmin, Vmax, colr=brown,"", 12, 0.5);
	plot_v_nod(ct, cn, findnodloc(ct,cn,120,0), Vmin, Vmax, colr=red,  "", 12, 0.5);

	plot_v_nod(ct, 1, 304,  		    Vmin, Vmax, colr=brown, "", 13, 0.5);
	plot_v_nod(ct, 3, 2726, 		    Vmin, Vmax, colr=magenta, "", 13, 0.5);


//	plot_chan_current(ct,cn, findnodloc(ct,cn,soma,0), K, 6, Imin = 0, Imax = iscal/40); plot_param("IKv3S", colr=red, plnum=9, plsize=0.3);
//	plot_chan_current(ct,cn, findnodloc(ct,cn,120,0), K, 6, Imin = 0, Imax = iscal/40); plot_param("IKv3D", colr=blue, plnum=9, plsize=0.3);

//	plot_chan_current(ct,cn, soma,                    NA, 8, Imin = -iscal/40, Imax = 0); plot_param("INaS", colr=red, plnum=8, plsize=0.3);
//	plot_chan_current(ct,cn, findnodloc(ct,cn,30,0),  NA, 8, Imin = -iscal/40, Imax = 0); plot_param("INaP", colr=brown, plnum=8, plsize=0.3);
//	plot_chan_current(ct,cn, findnodloc(ct,cn,69,0),  NA, 8, Imin = -iscal/40, Imax = 0); plot_param("INaI", colr=green, plnum=8, plsize=0.3);
//	plot_chan_current(ct,cn, findnodloc(ct,cn,120,0), NA, 8, Imin = -iscal/40, Imax = 0); plot_param("INaD", colr=blue, plnum=8, plsize=0.3);	
 

 //   plot_v_nod(ct, cn, 28, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 9, 0.5);
//    plot_v_nod(ct, cn, 18, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 8, 0.5);
//    plot_v_nod(ct, cn, 8, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 7, 0.5);
//    plot_v_nod(ct, cn, 6, Vmin = min(vstart,vstop),     Vmax = //0.05,      colr=cyan,    "", 6, 0.5);

////	plot_v_nod(ct, cn+1, soma, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=green,    "", 6, 0.5);
//////	plot_i_nod(ct, cn, 28, Imin = -iscal, Imax = iscal, colr=blue,    "", 5, 0.5);
//	plot_i_nod(ct, cn, findnodloc(ct,cn,100,-120), Imin = -iscal, Imax = 0, colr=blue,    "", 7, 0.5);
//	plot_i_nod(ct, cn, soma, Imin = -iscal, Imax = 0, colr=green,    "", 7, 0.5);
////	plot_i_nod(ct, cn+1, soma, Imin = -iscal, Imax = iscal, colr=green,    "", 4, 0.5);
////    plot_v_nod(ct, cn+1, 28, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 5, 0.5);
////    plot_v_nod(ct, cn+1, soma, Vmin = min(vstart,vstop),     Vmax = 0.05,    colr=blue,    "", 4, 0.5);
    //  plot_v_nod(ct2,cn, soma, Vmin = -0.1,     Vmax = 0.05, colr=blue,    "", 4, 0.5);
     // plot_i_nod(ct, cn, soma, Imin = -iscal, Imax = iscal, colr=magenta, "", 3, 1); 
     // plot_i_nod(ct2,cn, soma, Imin = -iscal, Imax = iscal, colr=brown, "", 2, 1); 
    if (set_vclamp > 0) {	
     // plot_var(&idiff,1,Imax= (outward?iscal:0),Imin = (outward?0:-iscal));  // plot current minus cap transient 
     plot_var(&idiff,1,Imax=0,Imin= -iscal);  // plot current minus cap transient 
        plot_param ("Isbac_soma", colr=brown, pl=3,plsize=0.6);

     switch (stimtype) {
        default: imax=0e-12; imin=-2000e-12; gmax=20e-9; break;
 	case 4:  imax=10e-12; imin=-10e-12; gmax=500e-12; break;
     }
     plot_func(isyn_tot,1,imax,imin);  	    plot_param("Itotsyn",green,2,0.3);
     plot_func(gsyn_tot,1,gmax,0);     	    plot_param("Gtotbp",blue, 1,0.3);
     plot_func(gsyn_tot,3,gmax*0.25,0);     plot_param("Gtotsbe",cyan, 1,0.3);
     plot_func(gsyn_tot,4,gmax*0.25,0);     plot_param("Gtotsbi",red, 1,0.3);

     if (sbacpair) {
        plot_var(&current,1,Imax= (outward?iscal:0),Imin = (outward?0:-iscal));  // plot current minus cap transient 
          plot_param ("Isbac_somax", colr=red, pl=2,plsize=0.5);
	plot_var(&current2,1,Imax= (outward?iscal:0),Imin = (outward?0:-iscal));  
          plot_param ("Isbac2_soma", colr=green, pl=2,plsize=0.5);
     }
//     plot_var(&cond,1,Imax= gscal,Imin = 0);				  // plot cond=current/driving force 
//        plot_param ("Gsbac_soma", colr=green, pl=1,plsize=0.2);
    }

     // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 

  }

    if (notinit(velocity))   velocity  = 1000;
    if (notinit(barwidth))   barwidth  = 100;
    if (notinit(barlength))  barlength = 200;
    if (notinit(theta))         theta  = 0;
    if (notinit(direction))  direction = 1;	
    if (notinit(waveshape))  waveshape = 1;		  // 1 -> square wave

    if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
    if (notinit(predur))        predur = 0.02;

    simtime = -predur;					  // must be set ahead of stim_backgr()
    setxmin = 0;

    stim_backgr(minten);
    if (disp) {			// display the stimulus
	double t;

      stimdur = move_stim(stimtime, barwidth, theta, velocity, ncycles, scontrast, direction, mask=1);
      display_size(500);
      disp_end = stimdur+0.05;
      for (starttime=simtime,t=stimtime; t<disp_end; starttime = t, t+= 0.002) {
	   display_stim(starttime, t, dscale=4, -0.035, -0.045); 
	   //display_stim(t, dscale=4, -0.035, -0.045); 
	   //display_stim(0+t, dscale=4); 
	   simwait(0.10);
      }

      return;
    }


//    if (notinit(set_vclamp)) set_vclamp = 0;
//    if (set_vclamp > 0) {	
//      vclamp              (ndn(ct, cn,soma),   vhold, simtime,  2);
//      if (sbacpair) vclamp(ndn(ct, cn+1,soma), vhold, simtime,  2);
//    }
	

    Gmax = 0;
//	stim_spot(100, 100, 0, scontrast, 0.001, 0.2);
//	stim_spot(1000, 50, 0, scontrast, stimtime, stimdur);

  if (velocity>0) {
      stimdur = ncycles * 2 * barwidth / velocity;
  }

  endexp=stimtime+stimdur+tailcurdur+poststimdur;
  step(predur);

  set_run_on_exit(runonexit);                         // set to erase savefile on ^C
  sprintf (savefile,"sbac_vclamp%06d",getpid());       // add pid to file name


  if (set_vclamp > 0) {
    dst = 0.0001;         // small time to allow voltage clamp to end

    if (ivplot) graph_pen(i+1,i+1,i+1,i+1,i+1);
    if (!ivplot) flag = true;
    if (vstart < vstop) sign = 1;
    else                sign = -1;

    savemodel (savefile);

    for (i=0,vpulse=vstart; (vpulse*sign)<=(vstop*sign+1e-6); i++,vpulse += vstep) {

       simtime = 0;
       pulsedur = move_stim(stimtime, barwidth, theta, velocity, ncycles, scontrast, direction, mask=1);
       pulsedur += prestimdur * 2;

       vclamp              (ndn(ct,cn,  soma), vhold, simtime,  prestimdur);
       if (sbacpair) vclamp(ndn(ct,cn+1,soma), vhold, simtime,  prestimdur);
       step (prestimdur);

       if (outward) maxCurrent =  0;               // for inward current
       else         maxCurrent =  1000;            // for outward current
       Gmax = 0;
       vclamp              (ndn(ct,cn,  soma), vpulse, simtime,  pulsedur);
       if (sbacpair) vclamp(ndn(ct,cn+1,soma), vpulse, simtime,  pulsedur);
       step (dst);
       if (ivplot) flag = true;
       step (pulsedur-2*dst);
       if (ivplot) graph(voltage, maxCurrent, Gmax);
       if (ivplot) flag = false;
       step (dst);

       vclamp              (ndn(ct,cn,  soma), tailvolt, simtime,  tailcurdur);
       if (sbacpair) vclamp(ndn(ct,cn+1,soma), tailvolt, simtime,  tailcurdur);
       step (tailcurdur);

       vclamp              (ndn(ct,cn,  soma), vhold, simtime,  poststimdur);
       if (sbacpair) vclamp(ndn(ct,cn+1,soma), vhold, simtime,  poststimdur);
       step (poststimdur);

       restoremodel (savefile);
    }
  }  /* if (set_vclamp) */ 
  else {
	// fprintf(stderr,"endexp %g %g %g\n", endexp, mvbrt1, mvbrt2);

       pulsedur = move_stim(stimtime, barwidth, theta, velocity, ncycles, scontrast, direction, mask=1);

	if (stimtype==1 || stimtype==2) {
	   savemodel (savefile);
	}
	step(mvbrt1+stimtime);

	if (stimtype==1 || stimtype==2) {
	  restoremodel (savefile);
	  step(mvbrt2-mvbrt1+stimtime);
	}

  } /* no vclamp */

  unlink (savefile);
  savefile[0] = 0;
}
