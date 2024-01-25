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
#include "onplot_dsgc_movie.cc"

int ct;
int ct2;
int rec_cn;
int rec_ct;
int ivplot;
int outward;
int sbacpair;
int stimtype;
int sbarr;
int n_dsgc;
int no_inhib;
int no_excit;
int direction;
int waveshape;
int sbac_color;
int sbr; 

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
double sbac_synanni;
double sbac_isynrngi;
double sbac_isynanpi;
double sbac_isynanpo;
double sbac_ithr;
double dbp1_cond;

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
double spotdia;
double sdia;

double barwidth;
double barlength;
double stim_theta;
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
double soma_z_sbac;
double soma_z_dsgc;
double radincr;

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
  defparams_dsgc_movie();
  defparams_onplot_movie();

  setptr("celltype",    &celltype);
  setptr("ivplot",      &ivplot);
  setptr("outward",     &outward);
  setptr("sbacpair",    &sbacpair);
  setptr("stimtype",    &stimtype);
  setptr("sbarr",       &sbarr);
  setptr("n_dsgc",	&n_dsgc);
  setptr("no_inhib",    &no_inhib);
  setptr("no_excit",    &no_excit);
  setptr("set_excit",   &set_excit);
  setptr("set_inhib",   &set_inhib);
  setptr("dbp1_cond",   &dbp1_cond);
  setptr("direction",   &direction);
  setptr("waveshape",   &waveshape);
  setptr("sbspac",      &sbspac);
  setptr("sbac_synrng", &sbac_synrng);
  setptr("sbac_synanp", &sbac_synanp);
  setptr("sbac_synanni", &sbac_synanni);
  setptr("sbac_isynrngi", &sbac_isynrngi);
  setptr("sbac_isynanpi", &sbac_isynanpi);
  setptr("sbac_isynanpo", &sbac_isynanpo);
  setptr("sbac_color",  &sbac_color);
  setptr("sbac_ithr",   &sbac_ithr);
  setptr("sbr",         &sbr);

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
  setptr("spotdia",   &spotdia);
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
  setptr("stim_theta",&stim_theta);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  setptr("soma_z_sbac", &soma_z_sbac);
  setptr("soma_z_dsgc", &soma_z_dsgc);
  setptr("radincr", &radincr);
  
  setptr("set_vclamp", &set_vclamp);
  setptr("set_tonic", &set_tonic);

  setptr("dcrm",&dcrm);

  nvalfile = "nval_sbac_stim.n";

  // inner radius of presynaptic output zone
  if (notinit(sbac_synanp)) { sbac_synanp = 60; }  // for nval file
  // for excitatory sbac-sbac interconnections:
  //
  if (notinit(sbac_synrng)) { sbac_synrng = 0; }  // for nval file
  if (notinit(sbac_synanp)) { sbac_synanp = 0; }  // for nval file

  // for inhibitory sbac-sbac interconnections:
  //
  if (notinit(sbac_isynrngi)) { sbac_isynrngi = 0; }  // (for nval file)
  if (notinit(sbac_isynanpi)) { sbac_isynanpi = 0; }  // inner radius of annulus in postsyn cell (for nval file)
  if (notinit(sbac_isynanpo)) { sbac_isynanpo = 0; }  // outer radius of annulus in postsyn cell (for nval file)
  if (notinit(sbac_synanni))  { sbac_synanni = 0;  }  // inner radius of annulus in presyn cell (for nval file)
  if (notinit(sbac_ithr))     { sbac_ithr = -0.05; }  // sbac -> sbac synapse threshold

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
  sprintf (tbuf,"%d>%d",spnt->node1b,spnt->node2b);	/* print pre and postsynaptic cell number */
  tlen = strlen(tbuf);
  gmove (length/2.0 -tlen*0.3, -1.0);
  gcwidth (1.2);
  gtext (tbuf);
}


/*------------------------------------------------------*/
   
void setparams(void)
{
   int i;

  make_dsgc = 1;
  ct2 = dsgc;
  if (notinit(n_dsgc)) n_dsgc = 1;
  make_sbac = 1;
  ct = sbac;
  make_sbac_dsgc = 1;
  make_dbp1_sbac = 1;
  make_dbp1_dsgc = 1;
  //make_cones = 1;
  // if (!notinit(set_vclamp) && set_vclamp>0) set_ncel(ct,2);
 // set_ncel(ct,1);

#define SBARR 30

 if (strcmp(sbac_file,"morph_sbac_168")==0) { 
    if (notinit (sbspac)) sbspac = 65;
 }
 if (notinit(sbr)) sbr = 2;			// other sbac to record from besides 1;

 if (!notinit(sbarr)) {
   sbxarr = (double *)emalloc(SBARR*sizeof(double));
   sbyarr = (double *)emalloc(SBARR*sizeof(double));
   sbtharr = (double *)emalloc(SBARR*sizeof(double));
   sbnarr = (int *)emalloc(SBARR*sizeof(int));

   for (i=0; i<SBARR; i++) sbnarr[i] = i+1;

   // if (notinit (sbtheta)) sbtheta = 0;
   if (notinit (sbspac)) sbspac = 85;

     if (sbarr==0) {
           sbxarr[0] =  0; sbyarr[0] = 0; sbtharr[0] = 0;  /* only 1 SBAC */
           n_sbac = 1;
     }
     if (sbarr==1) {                    /* opposing dendrites */
         sbxarr[0] =   0; sbyarr[0] = 0; sbtharr[0] = 0;
         sbxarr[1] =   1; sbyarr[1] = 0; sbtharr[1] = 30;
        n_sbac = 2;
     }
     if (sbarr==2) {                    /* aligned dendrites */
         sbxarr[0] = 0;      sbyarr[0] = 0; sbtharr[0] = 0;
         sbxarr[1] = sbspac; sbyarr[1] = 0; sbtharr[1] = 180;
         n_sbac = 2;
     }
     if (sbarr==3) {                    /* sbac3c morphology, 3 aligned */
        sbxarr[0] =    0; sbyarr[0] = 0; sbtharr[0] = 0; sbnarr[0]=1;
        sbxarr[1] =  -sbspac; sbyarr[1] = 0; sbtharr[1] = 30; sbnarr[1]=3;
        sbxarr[2] =   sbspac; sbyarr[2] = 0; sbtharr[2] = 60; sbnarr[2]=13;
        n_sbac = 3;
     }
     if (sbarr==4) {                    /* 3 aligned, 1 opposing */
        sbxarr[0] =   sbspac; sbyarr[0] = 0; sbtharr[0] = 0;
        sbxarr[1] =   sbspac; sbyarr[1] = 0; sbtharr[1] = 0;
        sbxarr[2] =        0; sbyarr[2] = 0; sbtharr[2] = 0;
        sbxarr[3] =  -sbspac; sbyarr[3] = 0; sbtharr[3] = 0;
        n_sbac = 4;
     }
     if (sbarr==7) {                    /* 7 overlapping */
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
       sbr = 3;
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

 onplot_dsgc_movie_init();             /* initialize dsgc movie stuff */
 onplot_movie_init();                  /* initialize onplot_movie stuff */
 rec_ct = sbac;
 rec_cn = -1;


 if (notinit(sbacpair)) sbacpair = 0;
 make_dbp1 = 1;

//  dbp2_file = dbp1_file;        /* Use the same morphology for dbp1 and dbp2 */

  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_xx" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */

  setn(ct2,MORPH,0);            /* set cell morphology from file, default = "morph_bp" */
  setn(ct2,BIOPHYS,1);          /* set cell biophys from file, default = "dens_default" */


  setn(sbac,NCOLOR,RCOLOR);       /* set cell display color from region */
  if (!notinit(sbac_color)) setn(ct,NCOLOR,sbac_color);   /* set cell display color from region */
  if (!notinit(soma_z_sbac)) setn(ct,SOMAZ,soma_z_sbac);            /* set cell soma z location */
  if (!notinit(soma_z_dsgc)) setn(dsgc,SOMAZ,soma_z_dsgc);

  set_synapse_dr (syn_draw2);

  R4 = SOMA;
  
/*
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
  */

  if (notinit(dispsize)) dispsize = 150;      /* default // display size */
  if (notinit(node_scale)) node_scale = -3.05;/* 3: nodenum, 0.05: small font */
  if (notinit(dend_dia_factor)) dend_dia_factor = 1.0;

  // if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5; 	/* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5; 	/* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2; 		/* default dia for dendrites */
  if (notinit(radincr)) radincr = 20;		/* radius increment for regions in sbac */  

  if (notinit(dcrm)) dcrm = 1e4;  	/* default cell Rm */
  if (notinit(sbaclm)) sbaclm = 0.1;  	/* default sbac complam (dens_sbaca.n) */
  if (notinit(dvrev)) dvrev = -0.05;  	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst)) dvst =   dvrev;  	/* default dbp1 vstart (dens_dbp1.n) */

  if (notinit(cone_maxcond)) cone_maxcond = 1000e-12;  	/* default cone OS cond */

  if (!notinit(dbp1_cond)) { setsv(dbp1,SCOND,3,dbp1_cond); }

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
       setn(sbac,NRM,set_drm);                    /* set default Rm */
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
   lbar = -cellrad - barwidth / 2;
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
      movesineann (0, 0, direction ,1e6,0, 180, barwidth*2, tfreq, 0,1.0, minten, scontrast, waveshape, stimtime,dur);
  else
      movesineann (0, 0, direction ,1e6,0, 250, barwidth*2, tfreq, 0,1.0, minten, scontrast, waveshape, stimtime,dur);
  mvbrt1 = mvbrt2 = dur;
 }
 else if (stimtype==4) {		// move spot
 	  send = int(200/spotdia + 0.5);
	  for (s=0; s<send; s++) {
            stim_spot(spotdia, s*spotdia, 0, scontrast, stimtime+s*spotdur, spotdur);
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
    int cn, i, j, n, s, plnum;
    int colr,pl,nsynap, nsynape, nsynapi;
    double dst, t, fmax,fmin;
    double cmin, cmax, rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax, Vminb, Vmaxb, Vmino;
    double Imin, Imax;
    double imin, imax, gmax;
    double vpulse, sign, pulsedur, maxdist;
    node *npnt;
    photorec *p;
    double disp_end, starttime;
    double dscale, mask, r;

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
  if (notinit(stimtime))   stimtime = .01;
  if (notinit(dstim))         dstim = .05;  
  if (notinit(spotdia))     spotdia = 30;      
  if (notinit(ncycles))     ncycles = 2;      
  if (notinit(stimtype))   stimtype = 1;      

//  endexp  = dtrial;
  ploti = 1e-4;

 
  if (notinit(vhold))       vhold  = -0.060;
  if (notinit(vstart))     vstart  = -0.120;
  if (notinit(vstop))       vstop  = -0.00;
  if (notinit(vstep))       vstep  =  0.03;
  if (notinit(tailvolt)) tailvolt  = vhold;
  if (notinit(gvrev))       gvrev  = vna;

  // midcbp  = findmid(ct,0,0);

  nsynap = synapse_add (10,dbp1,-1,-1,sbac,1);		/* make list of bipolar synapses onto sbac 1 for recording below */
  nsynape = nsynapi = 0;
  //for (i=2; i<=n_sbac; i++) {
  //  nsynape = synapse_add (3,dbp1,-1,-1,sbac,1,3);	/* make list of excit synapses from sbacs onto sbac1 */
  //  nsynap += nsynape;
  //  fprintf (stderr,"# nsynape %d\n",nsynape);
  //}
  for (j=1; j<=n_sbac; j++) {
    for (i=1; i<=n_sbac; i++) {
      if (i==j) continue;
      nsynapi = synapse_add (j,sbac,i,-1,sbac,j,4);	/* make list of inhib synapses from sbacs onto sbac1 */
      nsynap += nsynapi;
      fprintf (stderr,"# nsynapi %d->%d %d\n",i,j,nsynapi);
     }
  }
  fprintf (stderr,"# nsynap %d\n",nsynap);

  Vmin = -0.06;
  Vmax = -0.04;
  maxdist = 15;
  r = radincr;
  
  if (make_movie) {
    setonplot(onplot_movie); /* set movie plot routine */
    if (space_time) {  /* movie */
     plot_v_nod(sbac,1,findnodlocra(sbac,1,0,   stim_theta, maxdist), Vmin,Vmax, blue,"Vsoma",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,1*r, stim_theta, maxdist), Vmin,Vmax, green,   "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,2*r, stim_theta, maxdist), Vmin,Vmax, cyan,    "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,3*r, stim_theta, maxdist), Vmin,Vmax, red,     "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,4*r, stim_theta, maxdist), Vmin,Vmax, magenta, "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,5*r, stim_theta, maxdist), Vmin,Vmax, brown,  "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, stim_theta, maxdist), Vmin,Vmax, gray,   "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,7*r, stim_theta, maxdist), Vmin,Vmax, ltgreen,"",10,0.35); 
     
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,0,   stim_theta, maxdist), Vmin,Vmax, blue,"Vsoma",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,1*r, stim_theta, maxdist), Vmin,Vmax, green,   "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,2*r, stim_theta, maxdist), Vmin,Vmax, cyan,    "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,3*r, stim_theta, maxdist), Vmin,Vmax, red,     "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,4*r, stim_theta, maxdist), Vmin,Vmax, magenta, "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,5*r, stim_theta, maxdist), Vmin,Vmax, brown,  "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,6*r, stim_theta, maxdist), Vmin,Vmax, gray,   "",11,1.35); 
     plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,7*r, stim_theta, maxdist), Vmin,Vmax, ltgreen,"",11,1.35); 
    }
  }
  else {

/*
        plot_v_nod(dbp1,findmida(dbp1,0,0),     soma,Vminb=-0.05, Vmaxb= -0.03, blue, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,1*r,-1*r),  soma,Vminb=-0.05, Vmaxb= -0.03, green, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,2*r,-2*r),  soma,Vminb=-0.05, Vmaxb= -0.03, cyan, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,3*r,-3*r),  soma,Vminb=-0.05, Vmaxb= -0.03, red, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,4*r,-4*r),  soma,Vminb=-0.05, Vmaxb= -0.03, magenta, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,5*r,-5*r),soma,Vminb=-0.05, Vmaxb= -0.03, brown, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,6*r,-6*r),soma,Vminb=-0.05, Vmaxb= -0.03, gray, "", 22, 0.3);
	plot_v_nod(dbp1,findmida(dbp1,7*r,-7*r),soma,Vminb=-0.05, Vmaxb= -0.03, ltgreen, "", 22, 0.3);
*/

/*
	plot_synrate(findsynloc(dbp1,0,0),   cmin=0,cmax=500, 1, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,1*r,0),  cmin=0,cmax=500, 2, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,2*r,0),  cmin=0,cmax=500, 3, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,3*r,0),  cmin=0,cmax=500, 4, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,4*r,0),  cmin=0,cmax=500, 5, 21,"",0.3); 
	plot_synrate(findsynloc(dbp1,5*r,0), cmin=0,cmax=500, 6, 21,"",0.3); 
  	plot_synrate(findsynloc(dbp1,6*r,0), cmin=0,cmax=500, 7, 21,"",0.3); 
*/

/*
	plot_syncond(findsynloca(dbp1,-2*r,stim_theta),cmin=0,cmax=50e-12, ltred,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,-1*r,stim_theta),cmin=0,cmax=50e-12, ltmag,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,0,0),            cmin=0,cmax=50e-12, blue,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,1*r,stim_theta), cmin=0,cmax=50e-12, green,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,2*r,stim_theta), cmin=0,cmax=50e-12, cyan,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,3*r,stim_theta), cmin=0,cmax=50e-12, red,     20,"",0.3); 
	plot_syncond(findsynloca(dbp1,4*r,stim_theta), cmin=0,cmax=50e-12, magenta, 20,"",0.3); 
	plot_syncond(findsynloca(dbp1,5*r,stim_theta), cmin=0,cmax=50e-12, brown,   20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,6*r,stim_theta), cmin=0,cmax=50e-12, gray,    20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,7*r,stim_theta), cmin=0,cmax=50e-12, ltgreen, 20,"",0.3); 
*/
	/*
	 if (!notinit(sbarr) && sbarr > 0) {
	  plot_syncond(findsynlocra(sbac,1,4*r,stim_theta), cmin=0,cmax=50e-12, magenta, 19,"",0.3); 
	  plot_syncond(findsynlocra(sbac,1,5*r,stim_theta), cmin=0,cmax=50e-12, brown,   19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,1,6*r,stim_theta), cmin=0,cmax=50e-12, gray,    19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,1,7*r,stim_theta), cmin=0,cmax=50e-12, ltgreen, 19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,sbr,-7*r,stim_theta),cmin=0,cmax=50e-12, red,     19,"",0.3); 
  	  plot_syncond(findsynloc(dbp1,7*r,-7*r), cmin=0,cmax=50e-12, ltgreen, 19,"",0.3); 
	} */

        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,   0,     0, maxdist),      Vmin,Vmax, blue,"Vsoma_2",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-1*r, stim_theta, maxdist), Vmin,Vmax, green,    "",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-2*r, stim_theta, maxdist), Vmin,Vmax, cyan,     "",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-3*r, stim_theta, maxdist), Vmin,Vmax, red,     "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-4*r, stim_theta, maxdist), Vmin,Vmax, magenta, "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-5*r, stim_theta, maxdist), Vmin,Vmax, brown,   "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-6*r, stim_theta, maxdist), Vmin,Vmax, gray,    "", 12,0.5);	
        // plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-7*r, stim_theta, maxdist), Vmin,Vmax, ltgreen, "", 12,0.5); 

	// orthogonal nodes
        //plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, stim_theta+90, maxdist),  Vmin,Vmax, blue,      "",11,0.5); 
        //plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, stim_theta-90, maxdist),  Vmin,Vmax, green,     "",11,0.5); 

        //plot_v_nod(sbac,1,findnodlocra(sbac,1,-3*r, stim_theta, maxdist), Vmin,Vmax,red,       "",10,0.5); 
        //plot_v_nod(sbac,1,findnodlocra(sbac,1,-2*r, stim_theta, maxdist), Vmin,Vmax, cyan,     "",10,0.5); 
        //plot_v_nod(sbac,1,findnodlocra(sbac,1,-1*r, stim_theta, maxdist), Vmin,Vmax, green,    "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,  0,     0, maxdist),  Vmin,Vmax, blue, "Vsoma_1",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,1*r, stim_theta, maxdist),  Vmin,Vmax, green,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,2*r, stim_theta, maxdist),  Vmin,Vmax, cyan,      "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,3*r, stim_theta, maxdist),  Vmin,Vmax, red,       "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,4*r, stim_theta, maxdist),  Vmin,Vmax, magenta,   "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,5*r, stim_theta, maxdist),  Vmin,Vmax, brown,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, stim_theta, maxdist),  Vmin,Vmax, gray,      "",10,0.5); 
        // plot_v_nod(sbac,1,findnodlocra(sbac,1,7*r, stim_theta, maxdist),  Vmin,Vmax, ltgreen,   "",10,0.5); 

	

	if (set_vclamp > 0) { Vmin=vstart; if (stimtype==4) Vmax=vstart+0.01; else Vmax = 0; }
	else                { Vmin = -0.07; Vmax = -0.00; }

	Vmino = -0.06;

     switch (stimtype) {
        default: imax=0e-12; imin=-100e-12; gmax=5e-9; break;
 	case 4:  imax=10e-12; imin=-10e-12; gmax=500e-12; break;
     }
     /*
     plot_func(isyn_tot,10,imax,imin);  	    plot_param("Itotbpsyn",green,2,0.3);
     plot_func(gsyn_tot,10,gmax,0);     	    plot_param("Gtotbp",blue, 1,0.3);
  //   plot_func(gsyn_tot,3,gmax*0.25,0);     plot_param("Gtote_sb",cyan, 1,0.3);
     plot_func(gsyn_tot,1,gmax*0.25,0);       plot_param("Gtoti_sb1",red, 1,0.3);
     plot_func(gsyn_tot,sbr,gmax*0.25,0);     plot_param("Gtoti_sb2",magenta, 1,0.3);
     */

  }

    if (notinit(velocity))     velocity  = 1000;
    if (notinit(barwidth))     barwidth  = 100;
    if (notinit(barlength))    barlength = 200;
    if (notinit(stim_theta)) stim_theta  = 0;
    if (notinit(direction))  direction = 1;	
    if (notinit(waveshape))  waveshape = 1;		  // 1 -> square wave

    if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
    if (notinit(predur))        predur = 0.02;

    simtime = -predur;					  // must be set ahead of stim_backgr()
    setxmin = 0;

    stim_backgr(minten);
    if (!make_movie) {
     if (disp) {			// display the stimulus
	double t;

      stimdur = move_stim(stimtime, barwidth, stim_theta, velocity, ncycles, scontrast, direction, mask=0);
      // stim_spot(200, 0, 0, scontrast, 0.001, 0.2);
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
   }


    Gmax = 0;
//	stim_spot(100, 100, 0, scontrast, 0.001, 0.2);
//	stim_spot(1000, 50, 0, scontrast, stimtime, stimdur);

  if (velocity>0) {
      stimdur = ncycles * 2 * barwidth / velocity;
  }

   stimdur = move_stim(stimtime, barwidth, stim_theta, velocity, ncycles, scontrast, direction, mask=0);
   // stim_spot(200, 0, 0, scontrast, 0.001, 0.2);

  endexp=stimtime+stimdur+tailcurdur+poststimdur;
  ploti = (endexp > 1.0 ? 1e-3: ploti);

  step(predur);

   // fprintf(stderr,"endexp %g %g %g\n", endexp, mvbrt1, mvbrt2);

   step(stimdur+stimtime);
}
