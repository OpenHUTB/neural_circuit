/* Experiment sbac_stim */
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
int rec_cn;
int rec_ct;
int ivplot;
int elnode;
int outward;
int sbacpair;
int stimtype;
int sbarr;
int no_inhib;
int no_excit;
int direction;
int waveshape;
int sbac_color;
int sbr; 
int draw_synapse; 
int run_vclamp; 
int set_vclamp; 
int makenv; 
int syn_rad_dist; 
int use_stimfile; 

const char *celltype;

double kdr_cond;
double wac_g;
double wac_rad;
int    wac_nsyn;

double mglur2; 

double set_drm;
double elec_cap;
double elec_rs;
double nadist;
double namid;
double kdist;
double kmid;
double kprox;
double ksoma;
double camid;
double cadist;
double catmid;
double catdist;
double set_excit;
double set_inhib;
double set_synspac;
double sbspac;
double sbac_synrng;
double sbac_synanp;
double sbac_synanni;
double sbac_isynrngi;
double dbp1_anpi;
double dbp1_anpo;
double dbp2_anpi;
double dbp2_anpo;
double sbac_isynanpi;
double sbac_isynanpo;
double sbac_isynanni;
double sbac_ithr;
double sbac_maxsdist;
double dbp1_cond;
double dSomaDia;
double tipcap;
double sbsynang;

double axon_br_dia;
double varicos_dia;
double sbaclm;
double sb1mul;

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
double spdia;
double spdia_ratio;
double orad1;
double irad2;
double orad2;
double edist;
double idist;

double barwidth;
double barlength;
double stim_theta;
double rstim_theta;
double velocity;
double stimx;
double stimy;
double annrad;

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
double soma_z;
double radincr;

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
//  defparams_dsgc_movie();
//  defparams_onplot_movie();

  setptr("celltype",    &celltype);
  setptr("ivplot",      &ivplot);
  setptr("elnode",      &elnode);
  setptr("tipcap",      &tipcap);
  setptr("outward",     &outward);
  setptr("sbacpair",    &sbacpair);
  setptr("stimtype",    &stimtype);
  setptr("use_stimfile",&use_stimfile);
  setptr("sbarr",       &sbarr);
  setptr("no_inhib",    &no_inhib);
  setptr("no_excit",    &no_excit);
  setptr("set_excit",   &set_excit);
  setptr("set_inhib",   &set_inhib);
  setptr("set_synspac", &set_synspac);
  setptr("dbp1_cond",   &dbp1_cond);
  setptr("direction",   &direction);
  setptr("waveshape",   &waveshape);
  setptr("makenv",      &makenv);
  setptr("sbspac",      &sbspac);
  setptr("sbac_synrng", &sbac_synrng);
  setptr("sbac_synanp", &sbac_synanp);
  setptr("sbac_synanni", &sbac_synanni);
  setptr("sbac_isynrngi", &sbac_isynrngi);
  setptr("dbp1_anpi",   &dbp1_anpi);
  setptr("dbp1_anpo",   &dbp1_anpo);
  setptr("dbp2_anpi",   &dbp2_anpi);
  setptr("dbp2_anpo",   &dbp2_anpo);
  setptr("sbac_isynanpi", &sbac_isynanpi);
  setptr("sbac_isynanpo", &sbac_isynanpo);
  setptr("sbac_isynanni", &sbac_isynanni);
  setptr("sbac_color",  &sbac_color);
  setptr("sbac_ithr",   &sbac_ithr);
  setptr("sbac_maxsdist",   &sbac_maxsdist);
  setptr("sbr",         &sbr);
  setptr("draw_synapse",&draw_synapse);
  setptr("run_vclamp",  &run_vclamp);
  setptr("set_vclamp",  &set_vclamp);
  setptr("syn_rad_dist",  &syn_rad_dist);
  setptr("mglur2",  	&mglur2);
  setptr("sbsynang",  	&sbsynang);

  setptr("set_drm",     &set_drm);
  setptr("elec_cap",    &elec_cap);
  setptr("elec_rs",      &elec_rs);
  setptr("kdr_cond",    &kdr_cond);
  setptr("wac_g",       &wac_g);
  setptr("wac_rad",     &wac_rad);
  setptr("wac_nsyn",    &wac_nsyn);
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);
  setptr("sbaclm",      &sbaclm);
  setptr("sb1mul",      &sb1mul);
  setptr("nadist",      &nadist);
  setptr("namid",       &namid);
  setptr("kdist",       &kdist);
  setptr("kmid",        &kmid);
  setptr("kprox",       &kprox);
  setptr("ksoma",       &ksoma);
  setptr("cadist",      &cadist);
  setptr("camid",       &camid);
  setptr("catdist",     &catdist);
  setptr("catmid",      &catmid);

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
  setptr("spdia",     &spdia);
  setptr("spdia_ratio",&spdia_ratio);
  setptr("annrad",    &annrad);
  setptr("dSomaDia",  &dSomaDia);
  
  setptr("orad1",     &orad1);
  setptr("irad2",     &irad2);
  setptr("orad2",     &orad2);
  
  setptr("edist",     &edist);
  setptr("idist",     &idist);
  
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
  setptr("rstim_theta",&rstim_theta);
  setptr("stim_theta",&stim_theta);
  setptr("velocity",&velocity);
  setptr("stimx",&stimx);
  setptr("stimy",&stimy);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  setptr("soma_z", &soma_z);
  setptr("radincr", &radincr);
  
  setptr("set_tonic", &set_tonic);

  setptr("dcrm",&dcrm);

  nvalfile = "nval_sbac_stim.n";

  if (notinit(dSomaDia)) { dSomaDia = 10; }  // for nval file

  // for excitatory sbac-sbac interconnections:
  //
  sbac_synrng = 0;	// for nval file
  sbac_synanp = 0;	// for nval file
  sbsynang = 180;	// sbac -> dsgc synang2

  // for inhibitory sbac-sbac interconnections:
  //
  sbac_isynrngi = 0;	// (for nval file)
  dbp1_anpi = 0;	// inner radius of excit annulus in postsyn cell (for nval file)
  dbp1_anpo = 90;	// outer radius of excit annulus in postsyn cell (for nval file)
  dbp2_anpi = 0;	// inner radius of excit annulus in postsyn cell (for nval file)
  dbp2_anpo = 90;	// outer radius of excit annulus in postsyn cell (for nval file)
  sbac_isynanpi = 0;	// inner radius of inhib annulus in postsyn cell (for nval file)
  sbac_isynanpo = 0;	// outer radius of inhib annulus in postsyn cell (for nval file)
  sbac_synanni = 0;	// inner radius of annulus in presyn cell (for nval file)
  sbac_isynanni = 0;	// inner radius of inhib annulus in presyn cell (for nval file)
  sbac_ithr = -0.05;	// sbac -> sbac synapse threshold in nval_sbac_stim.n
  sbac_maxsdist = 20;	// sbac -> sbac synapse max dist in nval_sbac_stim.n
  wac_g = 0;		/* conductance of synapse onto sbac */

  make_sbac_sbac = 1;
  //dscavg = 2e7;				// sets Ca gain for release
  dscavg = 1e6;				// sets Ca gain for release
  dcapkm = 10e-6;			// sets ka for capump

  _CA_T = _CA6;         /* set type of T-type calcium channel for dens_sbac_chans.n */
                        /* then set it to slowly inactivating in chanparams_dsgc_chans */
  sdia = 0.3;
  spdia = 0.6;
  spdia_ratio = 0.9;
  // spdia = sdia * spdia_ratio;

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

  if (!draw_synapse) return;
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
  if (spnt->node1a==am) gpen (brown);
  else                  gpen (black);
  sprintf (tbuf,"%d>%d",spnt->node1b,spnt->node2b);	/* print pre and postsynaptic cell number */
  tlen = strlen(tbuf);
  gmove (length/2.0 -tlen*0.3, -1.0);
  gcwidth (2.5*dscale);
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
     if (sbarr==1) {                    /* aligned dendrites */
         sbxarr[0] =   0; sbyarr[0] = 0; sbtharr[0] = 0;
         sbxarr[1] =   1; sbyarr[1] = 0; sbtharr[1] = 0;
        n_sbac = 2;
     }
     if (sbarr==2) {                    /* opposing dendrites */
         sbxarr[0] = 0;       sbyarr[0] = 0; sbtharr[0] = 90;
         sbxarr[1] = -sbspac; sbyarr[1] = 0; sbtharr[1] = 330;
         n_sbac = 2;
     }
     if (sbarr==3) {                    /* sbac3c morphology, 3 aligned */
        sbxarr[0] =    0; sbyarr[0] = 0; sbtharr[0] = 0; sbnarr[0]=1;
        sbxarr[1] =  -sbspac; sbyarr[1] = 0; sbtharr[1] = 30; sbnarr[1]=2;
        sbxarr[2] =   sbspac; sbyarr[2] = 0; sbtharr[2] = 60; sbnarr[2]=3;
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
     if (sbarr==107) {			// 7 sbacs
	     double r;
       r = 20;
       sbxarr[0] =  0;             sbyarr[0] = 0;             sbtharr[0] = 0   + r;
       sbxarr[1] =  sbspac;        sbyarr[1] = 0;             sbtharr[1] = 180 + r;
       sbxarr[2] =  0.5*sbspac;    sbyarr[2] = -0.866*sbspac; sbtharr[2] = 60  + r;
       sbxarr[3] = -0.5*sbspac;    sbyarr[3] = -0.866*sbspac; sbtharr[3] = 0   + r;
       sbxarr[4] = -sbspac;        sbyarr[4] = 0;             sbtharr[4] = 90  + r;
       sbxarr[5] = -0.5*sbspac;    sbyarr[5] = 0.866*sbspac;  sbtharr[5] = 60  + r;
       sbxarr[6] =  0.5*sbspac;    sbyarr[6] = 0.866*sbspac;  sbtharr[6] = 0   + r;
       n_sbac = 7;
       sbr = 2;
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

 // onplot_dsgc_movie_init();             /* initialize dsgc movie stuff */
 // onplot_movie_init();                  /* initialize onplot_movie stuff */
 rec_ct = sbac;
 rec_cn = -1;


 if (notinit(sbacpair)) sbacpair = 0;
 make_dbp1 = 1;

 if (n_dsgc > 0) make_dsgc = 1;

//  dbp2_file = dbp1_file;        /* Use the same morphology for dbp1 and dbp2 */

  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_xx" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */

//  setn(ct2,MORPH,0);            /* set cell morphology from file, default = "morph_bp" */
//  setn(ct2,BIOPHYS,1);          /* set cell biophys from file, default = "dens_default" */


  setn(ct,NCOLOR,RCOLOR);       /* set cell display color from region */
  if (!notinit(sbac_color)) setn(ct,NCOLOR,sbac_color);   /* set cell display color */
  if (!notinit(soma_z)) setn(ct,SOMAZ,soma_z);            /* set cell soma z location */

  if (notinit(sbac_densfile2)) sbac_densfile2 = sbac_densfile;

  set_synapse_dr (syn_draw2);

  // SOMA      = R_4;
  
  DENDD     = DEND_DIST = R_1;
  DEND        		= R_2;
  DENDP     = DEND_PROX = R_3;
  SOMA       		= R_4;
  HCK       = HILLOCK   = R_5;
  AXONT     = AXON_THIN = R_6;
  AXON      		= R_7;
  AXONP     = AXON_PROX = R_7;
  AXOND     = AXON_DIST = R_8;
  VARIC     = VARICOS   = R_9;

  // if (notinit(dispsize)) dispsize = 150;      /* default // display size */
  if (notinit(node_scale)) node_scale = -3.05;/* 3: nodenum, 0.05: small font */
  if (notinit(dend_dia_factor)) dend_dia_factor = 1.0;

  // if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 2.0e4;      /* background light intensity */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5; 	/* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5; 	/* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2; 		/* default dia for dendrites */
  if (notinit(radincr)) radincr = 20;		/* radius increment for regions in sbac */  
  if (notinit(draw_synapse)) draw_synapse = 1;	/* enable drawing synapses */  

  if (notinit(dcrm)) dcrm = 1e4;  	/* default cell Rm */
  if (notinit(sbaclm)) sbaclm = 0.1;  	/* default sbac complam (dens_sbaca.n) */
  if (notinit(sb1mul)) sb1mul = 0.5;  	/* default sbac1 cmul factor in dens_sbacd.n (inhib input mul) */
  if (notinit(dvrev)) dvrev = -0.05;  	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst)) dvst =   dvrev;  	/* default dbp1 vstart (dens_dbp1.n) */

  if (notinit(cone_maxcond)) cone_maxcond = 1000e-12;  	/* default cone OS cond */

  if (!notinit(dbp1_cond)) { setsv(dbp1,SCOND,3,dbp1_cond); }
  

  if (!notinit(set_excit)) { setsv(sbac,SCOND,3,set_excit); }
  if (!notinit(no_excit)) { if (no_excit==1) { setsv(sbac,CELPRE,3,-1); setsv(sbac,CELPOST,3,-1); }}

  if (!notinit(set_synspac)) { setsv(sbac,SYNSPAC,4,set_synspac); }
  if (!notinit(set_inhib)) { setsv(sbac,SCOND,4,set_inhib); }
  if (!notinit(no_inhib)) { if (no_inhib==1) { setsv(sbac,CELPRE,4,-1); setsv(sbac,CELPOST,4,-1); }}

  if (notinit(syn_rad_dist)) { syn_rad_dist = 0; }

  if (notinit(ax_dia_factor)) ax_dia_factor = 1; /* multiplier for axon diameter */
  if (notinit(dvrev))   dvrev = -0.045;           /* Vrev for dens_dbp1.n */
  if (notinit(dvst))     dvst = -0.07;            /* Vstart for dens_dbp1.n */
  if (notinit(tipcap)) tipcap = 0;  

  if (notinit(ivplot)) ivplot = 0;              /* make I/V plot */
  if (notinit(outward)) outward = 0;            /* >0 => calc K cond, otherwise Na */
  if (notinit(iscal)) iscal   = 3e-9;           /* plot scale */
  if (notinit(gscal)) gscal   = 1.5e-8;         /* plot scale */
  if (notinit(sb_rm)) sb_rm   = 10e3;           /* sbac RM, used in dens_sbaca.n */
  if (notinit(sb_vs)) sb_vs   = -0.07;          /* sbac vstart, used in dens_sbaca.n */
  if (notinit(sb_vr)) sb_vr   = -0.07;          /* sbac vrev,   used in dens_sbaca.n */

  if (notinit(nadist)) nadist   = 3e-3;         /* distal Na density, dens_sbaca.n */
  if (notinit(namid))   namid   = 3e-3;         /* middle Na density, dens_sbaca.n */
  if (notinit(kdist))   kdist   = 2e-3;         /* dist Kdr density, dens_sbaca.n */
  if (notinit(kmid))     kmid   = 2e-3;         /* mid  Kdr density, dens_sbaca.n */
  if (notinit(kprox))   kprox   = 2e-3;         /* prox Kdr density, dens_sbaca.n */
  if (notinit(ksoma))   ksoma   = 3e-3;         /* soma Kdr density, dens_sbaca.n */
  if (notinit(cadist)) cadist   = 1e-3;         /* dist CaL density, dens_sbaca2.n */
  if (notinit(camid))   camid   = 1e-3;         /* mid  CaL density, dens_sbaca2.n */
  if (notinit(catdist)) catdist = 0e-3;         /* dist CaT density, dens_sbaca2.n */
  if (notinit(catmid))  catmid  = 0e-3;         /* mid  CaT density, dens_sbaca2.n */

  if (!notinit(set_drm)) {                      /* user set default Rm */
       setn(ct,NRM,set_drm);                    /* set default Rm */
       drm = set_drm;
  }
  if (notinit(sb_rid)) sb_rid = dri;
  if (notinit(sb_rii)) sb_rii = dri;
  if (notinit(sb_rip)) sb_rip = dri;

  if (notinit(use_stimfile)) use_stimfile = 0;

  if (notinit(mglur2)) mglur2 = 0;

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
    
     int cn;
    // ndens[sbac][cn=1] = 0;           // set cn 1 to use sbac_densfile
     ndens[sbac][cn=1] = 1;  		// set cn 1 to use sbac_densfile2, 50% cond
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
    double lbar, rbar, irad, orad, idia, odia;
    double inten, start, dur, wavel;
    double cellrad;
    double s,send;
    double tfreq;

 cellrad = max(xarrsiz,yarrsiz) * 0.525; 

 // stim_spot(100, 100,0, inten=minten,  start=0.02, dur=0.05, wavel=1, 0);

 if (stimtype==1) {			// move bar to & fro
   lbar = -cellrad - barwidth / 2;
   rbar = cellrad + barwidth/2;
   mvbrt1 = movebar (stimtime,            stimx,stimy, lbar, rbar, barwidth,barlength,theta,velocity,scontrast);
   mvbrt2 = movebar (mvbrt1+stimtime+0.05,stimx,stimy, rbar, lbar, barwidth,barlength,theta,velocity,scontrast);
   stim_spot(spotdia=1,2000,0,scontrast,-predur,2);	/* for wacs */
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

  tfreq = velocity/barwidth;
  dur = ncycles / tfreq;
  movesineann (0, 0, direction=1, annrad, 0,   0, barwidth, tfreq, 0,1.0,minten, scontrast, makenv, waveshape, stimtime,dur);
  movesineann (0, 0, direction=2, annrad, 0, 180, barwidth, tfreq, 0,1.0,minten, scontrast, makenv, waveshape, 2*stimtime+dur,dur);
  mvbrt1 = dur;
  mvbrt2 = 2 * dur + stimtime;
 }
 else if (stimtype==4) {		// move spot
 	  send = int(200/spotdia + 0.5);
	  for (s=0; s<send; s++) {
            stim_spot(spotdia, s*spotdia, 0, scontrast, stimtime+s*spotdur, spotdur);
          }
	  mvbrt2 = spotdur * send;

 } else if (stimtype==5) {		// move grating
      double sphase;
      int drift;

  tfreq = velocity/barwidth;  // speriod = barwidth
  dur = ncycles / tfreq;
  movegrating (0, 0, theta, sphase=0, barwidth, tfreq, drift=1, minten, scontrast, waveshape, stimtime,dur);
  movegrating (0, 0, theta, sphase=0, barwidth, tfreq, drift=-1, minten, scontrast, waveshape, 2*stimtime+dur,dur);
  mvbrt1 = dur;
  mvbrt2 = 2 * dur + stimtime;

 } else if (stimtype==6) {		// flashed annulus, proximal, distal

   if (notinit(orad1)) orad1 = 60;
   if (notinit(irad2)) irad2 = 75;
   if (notinit(orad2)) orad2 = 150;

   stim_spot(odia=2*orad1, 0,  0, scontrast, stimtime, spotdur);		// flashed proximal annulus

   stim_spot(odia=2*orad2, 0, 0, scontrast, 2*stimtime+prestimdur+spotdur, spotdur);		// flashed distal annulus
   stim_spot(idia=2*irad2, 0, 0, -scontrast, 2*stimtime+prestimdur+spotdur, spotdur);
   mvbrt2 = 2*stimtime+prestimdur+2*spotdur;
 }

   return mvbrt2;
}

/*------------------------------------------------------*/

void addlabels(void)

{
    int cn,i,j,nssyn,snode;
    sphere *srepnt = NULL;
    photorec *p;
    node *npnt;
#define CLOSESYNMAX 100
    int close_syn[CLOSESYNMAX];

        // label (ndn(sbac,1,96), red);
        //label (ndn(sbac,3,2726), red);
	// label (findsynloc(sbac,2,sbac,1,100, 0, 0.0, 50), red);
	// label (findsynloc(sbac,3,sbac,1,100, 0, 0.0, 50), red);

  if (notinit(wac_g))         wac_g = 0;	/* conductance of synapse onto sbac */
  if (notinit(wac_rad))    wac_rad  = 30;	/* outer limit radius from some of wac synapses onto sbac */
  if (notinit(wac_nsyn))   wac_nsyn = 10;	/* number of synapses per post-syn sbac */

  if (wac_g>0) { 			/* if wide-field amacrine cond */
      srepnt = (sphere *)at (am,1,soma,SPHERE);
      srepnt->nodp1->xloc = 0; 
      srepnt->nodp1->yloc = 0; 
      srepnt->nodp1->zloc = 0; 
      srepnt->dia = 1;
      initrand(1,5371);
      p = (photorec*)make_transducer(ndn(am,1,soma));
      p->xpos=2000;
      p->ypos=0;
      for (cn=1; cn<=nsbac; cn++) {
	   nssyn = 0;
           for (npnt=nodepnt; 
		npnt=foreach (npnt,sbac,cn,-1,NULL,NULL,&snode,wac_rad,dist2d,ndn(sbac,cn,soma)); 
		npnt=npnt->next) {
             close_syn[nssyn++] = snode;
	     if (nssyn>=CLOSESYNMAX) nssyn=CLOSESYNMAX-1;
           }
           fprintf(stderr,"# wfamac to sbac %d nssyn %d\n",cn,nssyn);
	   if (nssyn > 0) {
	      for (i=0; i<wac_nsyn; i++) { 
		   j = rrand(1) * nssyn; 
		   if (close_syn[j] == 0) { i--; continue; }
                   make_synapse_type (ndn(am,1,soma), ndn(sbac,cn,close_syn[j]), 0);
		   fprintf(stderr,"# wfamac to sbac %d %d\n",cn,close_syn[j]);
		   close_syn[j] = 0;
	      }
           }
      }
  }
}

/*------------------------------------------------------*/

void runexpt(void)

{
    int cn, tcn, snode, i, j, n, s, plnum, dbplist1, dbplist2, pa, pb;
    int electrode_node;
    int colr,pl,nsynap, nsynape, nsynapi, nsynapj, nssyn;
    double dst, t, fmax,fmin;
    double cmin, cmax, rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax, Vminb, Vmaxb;
    double Imin, Imax;
    double imin, imax, gmax;
    double vpulse, sign, pulsedur, maxdist;
    elem *epnt = NULL;
    synapse *sepnt = NULL;
    node *npnt;
    photorec *p;
    double disp_end, starttime;
    double dscale, mask, r;
    double dendang1, dendang2;

  timinc = 2e-6;
  ploti = 1e-4;
  crit = 1e-12;
  setonplot(onplot);
  electrode_node = 5000;

  for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
        p = (photorec*)make_transducer(ndn(dbp1,cn,soma), npnt->xloc, npnt->yloc);
	//	fprintf(stderr, "%g, %g\n", npnt->xloc, npnt->yloc);
   }
  cn = 1;

  if (notinit(run_vclamp))   run_vclamp = 0;
  if (notinit(set_vclamp))   set_vclamp = 0;
  if (notinit(prestimdur))   prestimdur = 0.02;
  if (notinit(stimdur))      stimdur    = 0.05;
  if (notinit(tailcurdur))   tailcurdur = 0;
  if (notinit(poststimdur)) poststimdur = 0.2;
  if (notinit(spotdur))      spotdur    = 0.1;
  if (notinit(minten))         minten = -0.050;
  if (notinit(scontrast))   scontrast = 0.005; 
  if (notinit(stimtime))     stimtime = .02;
  if (notinit(dstim))           dstim = .05;  
  if (notinit(spotdia))       spotdia = 30;      
  if (notinit(ncycles))       ncycles = 2;      
  if (notinit(stimtype))     stimtype = 1;      
  if (notinit(stim_theta))   stim_theta = 0;
  if (notinit(rstim_theta)) rstim_theta = 0;
  if (notinit(stimx))            stimx  = 0;
  if (notinit(stimy))            stimy  = 0;
  if (notinit(annrad))          annrad  = 150;

//  endexp  = dtrial;
  ploti = 1e-4;

 
  if (notinit(vhold))       vhold  = -0.060;
  if (notinit(vstart))     vstart  = -0.120;
  if (notinit(vstop))       vstop  = -0.00;
  if (notinit(vstep))       vstep  =  0.03;
  if (notinit(tailvolt)) tailvolt  = vhold;
  if (notinit(gvrev))       gvrev  = vna;

  if (run_vclamp!=0 || stimtype==5) {
      set_run_on_exit(runonexit);                         // set to erase savefile on ^C
      sprintf (savefile,"sbac_stim%06d",getpid());       // add pid to file name
  }

  // midcbp  = findmid(ct,0,0);

  /* compute statistics on synaptic connections: */

  if (notinit(edist))       edist  = 1000;	// outer limit on distance from soma to excit synapses
  if (notinit(idist))       idist  = 1000;	// outer limit on distance from soma to inhib synapses

  dbplist1 = 20;
  dbplist2 = 21;
  nsynap = synapse_add (dbplist1,dbp1,-1,-1,sbac,1,0,edist);  /* make list of bipolar syns onto sbac 1 for rec below */
  nsynap = synapse_add (dbplist2,dbp1,-1,-1,sbac,2,0,edist);  /* make list of bipolar syns onto sbac 2 for rec below */
  nsynape = nsynapi = 0;
  //for (i=2; i<=n_sbac; i++) {
  //  nsynape = synapse_add (3,dbp1,-1,-1,sbac,1,3);	/* make list of dbp1 synapses from sbacs onto sbac1 */
  //  nsynap += nsynape;
  //  fprintf (stderr,"# nsynape %d\n",nsynape);
  //}
  for (j=1; j<=n_sbac; j++) {
    nsynapj = 0;
    for (i=1; i<=n_sbac; i++) {
      if (i==j) continue;
      nsynapi = synapse_add (j,sbac,i,-1,sbac,j,4,0,idist);	/* make list of inhib synapses between sbacs */
      nsynapj += nsynapi;
      nsynap  += nsynapi;
      fprintf (stderr,"# nsynapi %d->%d %d\n",i,j,nsynapi);
     }
     fprintf (stderr,"# nsynapj  ->%d %d\n#\n",j,nsynapj);
  }
  // fprintf (stderr,"# nsynap %d\n",nsynap);

  /* print out the relative angles for SBAC-SBAC inhibition */

  tcn = 1;
  for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, sbac, -1, &pa, &pb); epnt=epnt->next){
	  double dendang;
     sepnt = (synapse *)epnt;
     if (epnt->node2a==sbac && epnt->node2b==tcn) {
       if (sepnt->vrev < -0.05) {			// look at inhibitory sbac -> sbac synapses
          dendang1 = node_angle(ndn(sbac,epnt->node1b,epnt->node1c), ndn(sbac,epnt->node1b,soma)); // presyn angle
          dendang2 = node_angle(ndn(sbac,epnt->node2b,epnt->node2c), ndn(sbac,epnt->node2b,soma)); // postsyn angle
          dendang =  (dendang1 - dendang2) * DEG;
          if (dendang < 0) dendang += 360;
          dendang = 180 - abs(dendang - 180);
          fprintf (stderr,"# relangle %g\n",dendang);
	  if (syn_rad_dist > 0) {
             fprintf (stderr,"# radial_dist %g\n",
		 dist2d(ndn(sbac,epnt->node2b,epnt->node2c),ndn(sbac,epnt->node2b,soma)));
	  }
       }
     }
  }

   if (notinit(elec_rs))  elec_rs   = 20e6;
   if (notinit(elec_cap)) elec_cap  = 1e-14;
//   if (notinit(shunt_res)) shunt_res  = 1e30;
//   if (notinit(shunt_v))   shunt_v  = dvrev;
   if (notinit(elnode))     elnode  = electrode_node;

   if (elnode==electrode_node) {
       make_electrode  (nd(ct,cn=1,elnode), ndn(ct,cn=1,soma), elec_rs, elec_cap, tipcap);
       // make_shunt  (ndn(ct,cn=1,soma), shunt_res, shunt_v);
   }

  Vmin = -0.06;
  Vmax = -0.02;
  maxdist = 15;
  r = radincr;
  
  if (make_movie) {
    setonplot(onplot_movie); /* set movie plot routine */
    if (space_time) {  /* movie */
     plot_v_nod(sbac,1,findnodlocra(sbac,1,0,   rstim_theta, maxdist), Vmin,Vmax, blue,"Vsoma",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,1*r, rstim_theta, maxdist), Vmin,Vmax, green,   "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,2*r, rstim_theta, maxdist), Vmin,Vmax, cyan,    "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,3*r, rstim_theta, maxdist), Vmin,Vmax, red,     "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,4*r, rstim_theta, maxdist), Vmin,Vmax, magenta, "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,5*r, rstim_theta, maxdist), Vmin,Vmax, brown,  "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, rstim_theta, maxdist), Vmin,Vmax, gray,   "",10,0.35); 
     plot_v_nod(sbac,1,findnodlocra(sbac,1,7*r, rstim_theta, maxdist), Vmin,Vmax, ltgreen,"",10,0.35); 
    }
  }
  else if (run_vclamp) {
	  double Imin, Imax;
      
	if (stimtype==6) { 
            Vmax = -0.02;
            Vmin = -0.07;
        } else {
            Vmax = 0;
            Vmin = -0.08;
	}
        plot_v_nod(sbac,1,findnodlocra(sbac,1,  0,     0, maxdist),        Vmin,Vmax, blue, "Vsoma",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,1*r, rstim_theta, maxdist),  Vmin,Vmax, green,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,2*r, rstim_theta, maxdist),  Vmin,Vmax, cyan,      "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,4*r, rstim_theta, maxdist),  Vmin,Vmax, magenta,   "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,8*r, rstim_theta, maxdist),  Vmin,Vmax, gray,      "",10,0.5); 

        if (stimtype==6) {Imin = -4e-10; Imax = 0e-10; }
	else             {Imin = -1e-9; Imax = 1e-9; }
        plot_i_nod(sbac,1,findnodlocra(sbac,1,  0,     0, maxdist),  Imin,Imax, blue, "Isoma",13,0.5); 
        plot_funci(isyn_tot,dbplist1,imax=0,imin=-500e-12);  plot_param("Itotbpsyn",green,2,0.3);

        if (stimtype==6) {
	plot_syncond(findsynloca(dbp1,-2*r,rstim_theta),cmin=0,cmax=200e-12, ltred,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,-1*r,rstim_theta),cmin=0,cmax=200e-12, ltmag,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,0,0),            cmin=0,cmax=200e-12, blue,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,1*r,rstim_theta), cmin=0,cmax=200e-12, green,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,2*r,rstim_theta), cmin=0,cmax=200e-12, cyan,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,3*r,rstim_theta), cmin=0,cmax=200e-12, red,     20,"",0.3); 
	plot_syncond(findsynloca(dbp1,4*r,rstim_theta), cmin=0,cmax=200e-12, magenta, 20,"",0.3); 
	plot_syncond(findsynloca(dbp1,5*r,rstim_theta), cmin=0,cmax=200e-12, brown,   20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,6*r,rstim_theta), cmin=0,cmax=200e-12, gray,    20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,7*r,rstim_theta), cmin=0,cmax=200e-12, ltgreen, 20,"",0.3); 
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

	plot_syncond(findsynloca(dbp1,-2*r,rstim_theta),cmin=0,cmax=200e-12, ltred,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,-1*r,rstim_theta),cmin=0,cmax=200e-12, ltmag,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,0,0),            cmin=0,cmax=200e-12, blue,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,1*r,rstim_theta), cmin=0,cmax=200e-12, green,   20,"",0.3); 
	plot_syncond(findsynloca(dbp1,2*r,rstim_theta), cmin=0,cmax=200e-12, cyan,    20,"",0.3); 
	plot_syncond(findsynloca(dbp1,3*r,rstim_theta), cmin=0,cmax=200e-12, red,     20,"",0.3); 
	plot_syncond(findsynloca(dbp1,4*r,rstim_theta), cmin=0,cmax=200e-12, magenta, 20,"",0.3); 
	plot_syncond(findsynloca(dbp1,5*r,rstim_theta), cmin=0,cmax=200e-12, brown,   20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,6*r,rstim_theta), cmin=0,cmax=200e-12, gray,    20,"",0.3); 
  	plot_syncond(findsynloca(dbp1,7*r,rstim_theta), cmin=0,cmax=200e-12, ltgreen, 20,"",0.3); 

	/*
	 if (!notinit(sbarr) && sbarr > 0) {
	  plot_syncond(findsynlocra(sbac,1,4*r,rstim_theta), cmin=0,cmax=50e-12, magenta, 19,"",0.3); 
	  plot_syncond(findsynlocra(sbac,1,5*r,rstim_theta), cmin=0,cmax=50e-12, brown,   19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,1,6*r,rstim_theta), cmin=0,cmax=50e-12, gray,    19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,1,7*r,rstim_theta), cmin=0,cmax=50e-12, ltgreen, 19,"",0.3); 
  	  plot_syncond(findsynlocra(sbac,sbr,-7*r,rstim_theta),cmin=0,cmax=50e-12, red,     19,"",0.3); 
  	  plot_syncond(findsynloc(dbp1,7*r,-7*r), cmin=0,cmax=50e-12, ltgreen, 19,"",0.3); 
	} */

        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,   0,     0, maxdist),      Vmin,Vmax, blue,"Vsoma",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-1*r, rstim_theta, maxdist), Vmin,Vmax, green,    "",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-2*r, rstim_theta, maxdist), Vmin,Vmax, cyan,     "",12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-3*r, rstim_theta, maxdist), Vmin,Vmax, red,     "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-4*r, rstim_theta, maxdist), Vmin,Vmax, magenta, "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-5*r, rstim_theta, maxdist), Vmin,Vmax, brown,   "", 12,0.5); 
        plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-6*r, rstim_theta, maxdist), Vmin,Vmax, gray,    "", 12,0.5); 
        // plot_v_nod(sbac,sbr,findnodlocra(sbac,sbr,-7*r, rstim_theta, maxdist), Vmin,Vmax, ltgreen, "", 12,0.5); 

	// orthogonal nodes
        plot_v_nod(sbac,1,findnodlocra(sbac,1, 6*r, rstim_theta+90, maxdist),  Vmin,Vmax, blue,      "",11,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1, 6*r, rstim_theta-90, maxdist),  Vmin,Vmax, green,     "",11,0.5); 

        plot_v_nod(sbac,1,findnodlocra(sbac,1,-3*r, rstim_theta, maxdist), Vmin,Vmax,red,       "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,-2*r, rstim_theta, maxdist), Vmin,Vmax, cyan,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,-1*r, rstim_theta, maxdist), Vmin,Vmax, green,    "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,  0,     0, maxdist),  Vmin,Vmax, blue, "Vsoma",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,1*r, rstim_theta, maxdist),  Vmin,Vmax, green,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,2*r, rstim_theta, maxdist),  Vmin,Vmax, cyan,      "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,3*r, rstim_theta, maxdist),  Vmin,Vmax, red,       "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,4*r, rstim_theta, maxdist),  Vmin,Vmax, magenta,   "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,5*r, rstim_theta, maxdist),  Vmin,Vmax, brown,     "",10,0.5); 
        plot_v_nod(sbac,1,findnodlocra(sbac,1,6*r, rstim_theta, maxdist),  Vmin,Vmax, gray,      "",10,0.5); 
        // plot_v_nod(sbac,1,findnodlocra(sbac,1,7*r, rstim_theta, maxdist),  Vmin,Vmax, ltgreen,   "",10,0.5); 
	
        // fprintf (stderr,"nodloc %d\n",findnodlocra(sbac,1,6*r, rstim_theta, maxdist));

	// Ca at inhibitory output synapses
	plot_ca_nod(sbac,1,findnodlocra(sbac,1, 8*r,rstim_theta), 1, 5e-6, red,    "",8,0.3);
	plot_ca_nod(sbac,2,findnodlocra(sbac,2,-8*r,rstim_theta), 1, 5e-6, magenta,"",8,0.3);

	// Ca at mid-radius, near dbp1 input synapses
	plot_ca_nod(sbac,1,findnodlocra(sbac,1, 3*r,rstim_theta), 1, 5e-6, blue,    "",8,0.3);
	plot_ca_nod(sbac,1,findnodlocra(sbac,1, 4.2*r,rstim_theta), 1, 5e-6, green,    "",8,0.3);
	plot_ca_nod(sbac,1,findnodlocra(sbac,1,-4*r,rstim_theta), 1, 5e-6, magenta,"",8,0.3);

	// conductance of inhibitory output synapses
	// plot_syncond(findsynlocra(sbac,1, 6*r, rstim_theta, -0.07), cmin=0,cmax=1.0e-9, red,     8,"",0.3); 
	// plot_syncond(findsynlocra(sbac,2,-6*r, rstim_theta, -0.07), cmin=0,cmax=1.0e-9, magenta, 8,"",0.3); 

     switch (stimtype) {
        default: imax=0e-12; imin=-200e-12; gmax=5e-9; break;
 	case 4:  imax=10e-12; imin=-10e-12; gmax=500e-12; break;
     }
     plot_funci(isyn_tot,dbplist1,imax,imin);  plot_param("Itotbpsyn",green,2,0.3);
     plot_funcg(gsyn_tot,dbplist1,gmax,0);     plot_param("Gtotbp_sb1",blue, 1,0.3);
     plot_funcg(gsyn_tot,dbplist2,gmax,0);     plot_param("Gtotbp_sb2",green, 1,0.3);
  //   plot_funcg(gsyn_tot,3,gmax,0);            plot_param("Gtote_sb",cyan, 1,0.3);
     plot_funcg(gsyn_tot,1,  gmax,0);          plot_param("Gtoti_sb1",red,     1,0.3);
     plot_funcg(gsyn_tot,sbr,gmax,0);          plot_param("Gtoti_sb2",magenta, 1,0.3);

  }
    if (notinit(predur))        predur = 0.02;
    if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
    if (notinit(velocity))     velocity  = 1000;
    if (notinit(barwidth))     barwidth  = 100;
    if (notinit(barlength))    barlength = 200;
    if (notinit(direction))    direction = 1;	
    if (notinit(waveshape))    waveshape = 1;		  // 1 -> square wave
    if (notinit(makenv))          makenv = 1;		  // 1 -> gaussian envelope, 0 = hard limit

    simtime = -predur;					  // must be set ahead of stim_backgr()
    setxmin = 0;
    // setxmin = simtime;

    if (use_stimfile && (stimtype==3 || stimtype==5)) {		// if sineann stimulus (and --use_stimfile 1 )
	     const char *fname;
	     char stimfile[100];
        if      (stimtype==3) fname = "stim_sbac_sineann_%g_%g";	// set the stimulus file name
	else if (stimtype==5) fname = "stim_sbac_grating_%g_%g";	
        sprintf (stimfile,fname,scontrast,barwidth);
        stim_file (stimfile);					// set the stimulus file name
        stim_blur (0,0,0,0,1);					// set the blur and stimulus arrays
    }
    stim_backgr(minten);
    if (!make_movie) {
     if (disp) {			// display the stimulus
	double t;

      stimdur = move_stim(stimtime, barwidth, stim_theta, velocity, ncycles, scontrast, direction, mask=0);
      // stim_spot(200, 0, 0, scontrast, 0.001, 0.2);
      // display_size(500);
      display_size(1000);
      disp_end = stimtime+stimdur+0.05;
      for (starttime=simtime,t=stimtime; t<disp_end; starttime = t, t+= 0.001) {
	   display_stim(starttime, t, dscale=4, -0.035, -0.045); 
	   // display_stim(starttime, t, dscale=4, -0.025, -0.035); 
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

  if (run_vclamp==0) {

     if (velocity>0) {
         stimdur = ncycles * 2 * barwidth / velocity;
     }

     stimdur = move_stim(stimtime, barwidth, stim_theta, velocity, ncycles, scontrast, direction, mask=0);
     // stim_spot(200, 0, 0, scontrast, 0.001, 0.2);

     endexp=stimtime+stimdur+tailcurdur+poststimdur;

     if (stimtype==3 || stimtype==5) {
        ploti = 1e-3;
     } else {
        ploti = (endexp > 1.0 ? 1e-3: ploti);
     }

     if (set_vclamp==1) { 			/* vclamp at soma, only vhold */
	  double Imin, Imax;
	Imin = -5e-10; Imax = 5e-10; 
        // plot_i_nod(sbac,1,findnodlocra(sbac,1,  0,     0, maxdist),  Imin,Imax, blue, "Isoma",13,0.5); 
        plot_i_nod(sbac,1,elnode,  Imin,Imax, blue, "Ielec",13,0.5); 
	if (elnode==electrode_node) 
                plot_v_nod(sbac,1,elnode,  Vmin,Vmax, blue, "Velec",10,0.5); 
        vclamp (ndn(sbac,1,elnode), vhold, simtime, endexp+predur);
     }
     step(predur);
       // fprintf(stderr,"simtime %g endexp %g %g %g\n", simtime,endexp, mvbrt1, mvbrt2);
     
     if (mglur2 > 0) set2ndmsg(CA,-1,dbp1,sbac,AMPA,10.0,mglur2);	// set glutamate to modulate Ca chans

     if (stimtype==5) {
        savemodel (savefile);
        step(stimtime+mvbrt1);
        restoremodel (savefile);
        step(stimtime+mvbrt1+poststimdur);
	unlink (savefile);
	savefile[0] = 0;
     }
     else step(endexp);
  }
  else {  /* vclamp mode */

     predur = 0.05;
     simtime = -predur;					  // must be set ahead of stim_backgr()
     setxmin = 0;

     if (stimtype==6) {
         stimdur = move_stim(stimtime, barwidth, stim_theta, velocity, ncycles, scontrast, direction, mask=0);
     }
     endexp=prestimdur+(stimdur+poststimdur);

     vclamp (ndn(sbac,1,elnode), vhold,         simtime,  predur);
     step (predur);
     if (mglur2 > 0) set2ndmsg(CA,-1,dbp1,sbac,10.0,mglur2,1);	// set glutamate to modulate Ca chans

     if (stimtype==6) {
        vclamp (ndn(sbac,1,elnode), vhold,         simtime,  endexp);
        step (endexp);
     } else {

       savemodel (savefile);
       vclamp (ndn(sbac,1,elnode), vhold,         simtime,  prestimdur);
       step (prestimdur);
       vclamp (ndn(sbac,1,elnode), vpulse=0,      simtime,  stimdur);
       step (stimdur);
       vclamp (ndn(sbac,1,elnode), vhold,         simtime,  poststimdur);
       step (poststimdur);
       restoremodel (savefile);

       simtime = 0;

       vclamp (ndn(sbac,1,elnode), vhold,         simtime,  prestimdur);
       step (prestimdur);
       vclamp (ndn(sbac,1,elnode), vpulse= -0.07, simtime,  stimdur);
       step (stimdur);
       vclamp (ndn(sbac,1,elnode), vhold,         simtime,  poststimdur);
       step (poststimdur);

     }
     unlink (savefile);
     savefile[0] = 0;
  }
}
