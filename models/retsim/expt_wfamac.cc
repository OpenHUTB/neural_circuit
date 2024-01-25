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
int ivplot;
int outward;
int amacpair;
int stimtype;
int amarr;
int no_inhib;
int no_excit;
int direction;
int waveshape;
int amac_color;
int stimtyp;
int rec_cn;
int rec_ct;

const char *celltype;

double kdr_cond;
double ka_cond;
double na_cond;
double na_condp;
double kdr_condp;
double ka_condp;

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
double dbp1_scond;
double am_am_scond;
double dbp1_svnoise;
double set_inhib;
double amspac;
double amac_synrng;
double amac_synanp;
double amac_isynrngi;
double amac_isynanpi;
double amac_isynanpo;

double axon_br_dia;
double varicos_dia;
double amaclm;
double wf_dia;

double speriod;
double sphase;
double orient;
double tfreq;
int    drift;

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
double xloc;
double yloc;
double sdia2;
double xloc2;
double yloc2;
double sdia3;
double xloc3;
double yloc3;
double sdia4;
double xloc4;
double yloc4;

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

//void am_init(void);

double am_vr;
double am_vs;
double am_rm;
double am_rid;			// distal Ri
double am_rii;			// intermediate Ri
double am_rip;			// proximal Ri

char savefile[30] = {0};


/*------------------------------------------------------*/

void defparams(void) 
{
  defparams_dsgc_movie();
  defparams_onplot_movie();
  
  setptr("celltype",    &celltype);
  setptr("ivplot",      &ivplot);
  setptr("outward",     &outward);
  setptr("amacpair",    &amacpair);
  setptr("stimtype",    &stimtype);
  setptr("amarr",       &amarr);
  setptr("no_inhib",    &no_inhib);
  setptr("no_excit",    &no_excit);
  setptr("set_excit",   &set_excit);
  setptr("dbp1_scond",   &dbp1_scond);
  setptr("am_am_scond", &am_am_scond);
  setptr("dbp1_svnoise", &dbp1_svnoise);
  setptr("set_inhib",   &set_inhib);
  setptr("direction",   &direction);
  setptr("waveshape",   &waveshape);
  setptr("amspac",      &amspac);
  setptr("amac_synrng", &amac_synrng);
  setptr("amac_synanp", &amac_synanp);
  setptr("amac_isynrngi", &amac_isynrngi);
  setptr("amac_isynanpi", &amac_isynanpi);
  setptr("amac_isynanpo", &amac_isynanpo);
  setptr("amac_color",  &amac_color);
  setptr("wf_dia",      &wf_dia);
  setptr("stimtyp",	 &stimtyp);
  
  setptr("set_drm",     &set_drm);
  setptr("kdr_cond",    &kdr_cond);
  setptr("ka_cond",    &ka_cond);
  setptr("na_cond",    &na_cond);
  setptr("kdr_condp", &kdr_condp);
  setptr("ka_condp", &ka_condp);
  setptr("na_condp", &na_condp);
  
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);
  setptr("amaclm",      &amaclm);
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
  setptr("poststimdur",	&poststimdur);
   setptr("ncycles",	&ncycles);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("stimtime",  &stimtime);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  setptr("xloc",      &xloc);
  setptr("yloc",      &yloc);
  setptr("sdia2",      &sdia2);
  setptr("xloc2",      &xloc2);
  setptr("yloc2",      &yloc2);
  setptr("sdia3",      &sdia3);
  setptr("xloc3",      &xloc3);
  setptr("yloc3",      &yloc3);
  setptr("sdia4",      &sdia4);
  setptr("xloc4",      &xloc4);
  setptr("yloc4",      &yloc4);
  
  setptr("speriod",   &speriod);
  setptr("sphase",    &sphase);
  setptr("orient",    &orient);
  setptr("tfreq",     &tfreq);
  setptr("drift",     &drift);
  
  
  setptr("vstart",&vstart);
  setptr("vstop",&vstop);
  setptr("vstep",&vstep);
  setptr("vhold",&vhold);
  setptr("tailvolt",&tailvolt);
  setptr("gvrev",&gvrev);
  setptr("am_vs",&am_vs);
  setptr("am_vr",&am_vr);
  setptr("am_rm",&am_rm);

  setptr("am_rid",&am_rid);
  setptr("am_rii",&am_rii);
  setptr("am_rip",&am_rip);
	
  setptr("velocity",&velocity);
  setptr("barwidth",&barwidth);
  setptr("barlength",&barlength);
  setptr("theta",&theta);

  setptr("iscal", &iscal);
  setptr("gscal", &gscal);
  
  setptr("set_vclamp", &set_vclamp);
  setptr("set_tonic", &set_tonic);

  setptr("dcrm",&dcrm);

  nvalfile = "nval_wfamac.n";

  // for excitatory interconnections:
  //
  if (notinit(amac_synrng)) { amac_synrng = 0; }  // for nval file
  if (notinit(amac_synanp)) { amac_synanp = 0; }  // for nval file

  // for inhibitory interconnections:
  //
  if (notinit(amac_isynrngi)) { amac_isynrngi = 0; }  // (for nval file)
  if (notinit(amac_isynanpi)) { amac_isynanpi = 0; }  // inner radius of annulus in syn cell (for nval file)
  if (notinit(amac_isynanpo)) { amac_isynanpo = 0; }  // outer radius of annulus in postsyn cell (for nval file)

  //make_am_am = 0;

  chanparamsfile = "chanparams_wfa";
  //am_init();

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
       ct  = am;
  }
  make_ct(ct);          /* make the cell type */
  // if (!notinit(set_vclamp) && set_vclamp>0) set_ncel(ct,2);
  //set_ncel(ct,1);

#define AMARR 30

  onplot_dsgc_movie_init();             /* initialize dsgc movie stuff */
  onplot_movie_init();                  /* initialize onplot_movie stuff */
  rec_cn=-1;
  rec_ct=-1;
 if (!notinit(amarr)) {
   amxarr = (double *)emalloc(AMARR*sizeof(double));
   amyarr = (double *)emalloc(AMARR*sizeof(double));
   amtharr = (double *)emalloc(AMARR*sizeof(double));
   amnarr = (int *)emalloc(AMARR*sizeof(int));

   for (i=0; i<AMARR; i++) amnarr[i] = i+1;

   // if (notinit (amtheta)) amtheta = 0;
   if (notinit (amspac)) amspac = 200;

   if (strcmp(am_file,"morph_wfamac1")==0 || strcmp(am_file,"morph_wfamac2")==0) {
     if (amarr==0) {
           amxarr[0] =  0; amyarr[0] = 0; amtharr[0] = 0;  /* only 1 AM AC */
           n_am = 1;
     }
     if (amarr==1) {                    /* amac3c morphology, opposing dendrites */
        amxarr[0] =  85; amyarr[0] = 0; amtharr[0] = 0;
        amxarr[1] = -85; amyarr[1] = 0; amtharr[1] = 0;
        n_am = 2;
     }
     if (amarr==2) {                    /* amac3c morphology, aligned dendrites */
         amxarr[0] =   0; amyarr[0] = 0; amtharr[0] = 0;
         amxarr[1] =   1; amyarr[1] = 0; amtharr[1] = 30;
         n_am = 2;
     }
     if (amarr==3) {                    /* amac3c morphology, 3 aligned */
        amxarr[0] =    0; amyarr[0] = 0; amtharr[0] = 0; amnarr[0]=1;
        amxarr[1] =  -amspac; amyarr[1] = 0; amtharr[1] = 30; amnarr[1]=3;
        amxarr[2] =   amspac; amyarr[2] = 0; amtharr[2] = 60; amnarr[2]=13;
        n_am = 3;
     }
     if (amarr==4) {                    /* amac3c morphology, 3 aligned, 1 opposing */
        amxarr[0] =   amspac; amyarr[0] = 0;      amtharr[0] = 0;
        amxarr[1] =   amspac; amyarr[1] = amspac; amtharr[1] = 90;
        amxarr[2] =        0; amyarr[2] = 0;      amtharr[2] = 180;
        amxarr[3] =  -amspac; amyarr[3] = 0;      amtharr[3] = 270;
        n_am = 4;
     }
     if (amarr==5) {			// 5 amacs
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;
       amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;
       amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 0;
       amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] = 90;
       n_am = 5;
     }
     if (amarr==7) {                    /* amac3c morphology, 7 overlapping */
        amxarr[0] =  3*amspac; amyarr[0] = 0; amtharr[0] = 0;
        amxarr[1] =  2*amspac; amyarr[1] = 0; amtharr[1] = 0;
        amxarr[2] =    amspac; amyarr[2] = 0; amtharr[2] = 0;
        amxarr[3] =         0; amyarr[3] = 0; amtharr[3] = 0;
        amxarr[4] =   -amspac; amyarr[4] = 0; amtharr[4] = 0;
        amxarr[5] = -2*amspac; amyarr[5] = 0; amtharr[5] = 0;
        amxarr[6] = -3*amspac; amyarr[6] = 0; amtharr[6] = 0;
        n_am = 7;
      }
     if (amarr==9) {			// 9 amacs
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;
       amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;
       amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 270;
       amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] =  30;
       amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  80;
       amxarr[6] = -amspac; amyarr[6] =  amspac; amtharr[6] = 270;
       amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] =  90;
       amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 230;
       n_am = 9;
     }
     if (amarr==102) {
       amxarr[0] =  145; amyarr[0] =  50; amtharr[0] = 0;
       amxarr[1] =  145; amyarr[1] = -50; amtharr[1] = 0;
       n_am = 2;
     }
     if (amarr==12){
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;
       amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;
       amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 270;
       amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] =  30;
       amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  80;
       amxarr[6] = -amspac; amyarr[6] =  amspac; amtharr[6] = 270;
       amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] =  90;
       amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 230;
       amxarr[9]= 2*amspac; amyarr[9] =  amspac; amtharr[9] = 70;
       amxarr[10]=2*amspac; amyarr[10]= 0;       amtharr[10]=  50;
       amxarr[11]=2*amspac; amyarr[11]=-amspac;  amtharr[11]= 45;
       n_am = 12;
     }
     if (amarr==105) {			// 13 amacs
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;
       amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;
       amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 0;
       amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] = 90;

       amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  0;
       amxarr[6] = -amspac; amyarr[6] =  amspac; amtharr[6] = 120;
       amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] = 180;
       amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 210;

       amxarr[9] = -2*amspac; amyarr[9] =  0;         amtharr[9] = 120;
       amxarr[10] = 2*amspac; amyarr[10] = 0;         amtharr[10] = 90;
       amxarr[11] = 0;        amyarr[11] = -2*amspac; amtharr[11] = 60;
       amxarr[12] = 0;        amyarr[12] =  2*amspac; amtharr[12] = 30;
       n_am = 13;
     }
     if (amarr==1052) {			// 2 amacs: 1, 2
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;   amnarr[0] = 1;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;  amnarr[1] = 2;
       //amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;amnarr[2] = 3;
       //amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 0; amnarr[3] = 4;
       //amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] = 90;amnarr[4] = 5;

       //amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  0;   amnarr[5] = 6;
       //amxarr[1] = -amspac; amyarr[1] =  amspac; amtharr[1] = 120;  amnarr[6] = 7;
       //amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] = 180; amnarr[7] = 8;
       //amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 210; amnarr[8] = 9;

       //amxarr[9] = -2*amspac; amyarr[9] =  0;         amtharr[9] = 120; amnarr[9] = 10;
       //amxarr[10] = 2*amspac; amyarr[10] = 0;         amtharr[10] = 90; amnarr[10] = 11;
       //amxarr[11] = 0;        amyarr[11] = -2*amspac; amtharr[11] = 60; amnarr[11] = 12;
       //amxarr[12] = 0;        amyarr[12] =  2*amspac; amtharr[12] = 30; amnarr[12] = 13;
       n_am = 2;
     }
     if (amarr==1053) {			// 3 amacs: 1, 2, 3
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0;    amnarr[0] = 1;
       amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;   amnarr[1] = 2;
       amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;   amnarr[2] = 3;
       //amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 0;  amnarr[2] = 4;
       //amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] = 90; amnarr[4] = 5;

       //amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  0; amnarr[5] = 6;
       //amxarr[1] = -amspac; amyarr[1] =  amspac; amtharr[1] = 120; amnarr[6] = 7;
       //amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] = 180; amnarr[7] = 8;
       //amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 210; amnarr[8] = 9;

       //amxarr[9] = -2*amspac; amyarr[9] =  0;         amtharr[9] = 120; amnarr[9] = 10;
       //amxarr[10] = 2*amspac; amyarr[10] = 0;         amtharr[10] = 90; amnarr[10] = 11;
       //amxarr[11] = 0;        amyarr[11] = -2*amspac; amtharr[11] = 60; amnarr[11] = 12;
       //amxarr[12] = 0;        amyarr[12] =  2*amspac; amtharr[12] = 30; amnarr[12] = 13;
       n_am = 3;
     }
     if (amarr==1057) {			// 2 amacs: 1, 7
       amxarr[0] =  0;      amyarr[0] = 0;       amtharr[0] = 0; amnarr[0] = 1;
       //amxarr[1] = -amspac; amyarr[1] = 0;       amtharr[1] = 30;
       //amxarr[2] =  amspac; amyarr[2] = 0;       amtharr[2] = 60;
       //amxarr[3] = 0;       amyarr[3] = -amspac; amtharr[3] = 0;
       //amxarr[4] = 0;       amyarr[4] =  amspac; amtharr[4] = 90;

       //amxarr[5] = amspac;  amyarr[5] =  amspac; amtharr[5] =  0;
       amxarr[1] = -amspac; amyarr[1] =  amspac; amtharr[1] = 120; amnarr[1] = 27;
       //amxarr[7] = amspac;  amyarr[7] =  -amspac; amtharr[7] = 180;
       //amxarr[8] = -amspac; amyarr[8] =  -amspac; amtharr[8] = 210;

       //amxarr[9] = -2*amspac; amyarr[9] =  0;         amtharr[9] = 120;
       //amxarr[10] = 2*amspac; amyarr[10] = 0;         amtharr[10] = 90;
       //amxarr[11] = 0;        amyarr[11] = -2*amspac; amtharr[11] = 60;
       //amxarr[12] = 0;        amyarr[12] =  2*amspac; amtharr[12] = 30;
       n_am = 2;
     }
      if(amarr==33) {		//2 amacs one loop
       amxarr[0] = -260;	amyarr[0] = 200;	amtharr[0] =-115;	amnarr[0] = 1;
       amxarr[1] = 0;	amyarr[1] = 0;   amtharr[1] =0;  amnarr[1] = 1;
       n_am = 2;
      }
      if(amarr==20) {//4amacs spread out
	amxarr[0] = 0;			 amyarr[0] = 0; 		amtharr[0] = 25;
	amxarr[1] = -1.5*amspac;	 amyarr[1] = 0;			amtharr[1] = -25;
	amxarr[2] = 0;		  	 amyarr[2] = -1.5*amspac;	amtharr[2] = 155;
	amxarr[3] = -1.5*amspac;	 amyarr[3] = -1.5*amspac;	amtharr[3] = 205;
	n_am = 4;
      }
   }
 }
 am_densfile = "dens_wfamac.n";

 if (notinit(amacpair)) amacpair = 0;
 make_dbp1 = 1;

//  dbp2_file = dbp1_file;        /* Use the same morphology for dbp1 and dbp2 */

  setn(ct,MORPH,0);		/* set cell morphology from file, default = "morph_xx" */
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */

//  setn(ct2,MORPH,0);            /* set cell morphology from file, default = "morph_bp" */
//  setn(ct2,BIOPHYS,1);          /* set cell biophys from file, default = "dens_default" */


  setn(ct,NCOLOR,RCOLOR);       /* set cell display color from region */
if (!notinit(amac_color)) setn(ct,NCOLOR,amac_color);   /* set cell display color from region */

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
  if (notinit(wf_dia)) wf_dia = 0.15; 		/* default dia for dendrites */
  

  if (notinit(dcrm)) dcrm = 1e4;  	/* default cell Rm */
  if (notinit(amaclm)) amaclm = 0.1;  	/* default amac complam (dens_amaca.n) */
  if (notinit(dvrev)) dvrev = -0.05;  	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst)) dvst =   dvrev;  	/* default dbp1 vstart (dens_dbp1.n) */

  if (notinit(cone_maxcond)) cone_maxcond = 1000e-12;  	/* default cone OS cond */

  if (!notinit(set_excit)) { setsv(am,SCOND,3,set_excit); }
  if (!notinit(no_excit)) { if (no_excit==1) { setsv(am,CELPRE,3,-1); setsv(am,CELPOST,3,-1); }}
  
  if (!notinit(am_am_scond)){setsv(am,SCOND,3,am_am_scond); }

  if (!notinit(set_inhib)) { setsv(am,SCOND,4,set_inhib); }
  if (!notinit(no_inhib)) { if (no_inhib==1) { setsv(am,CELPRE,4,-1); setsv(am,CELPOST,4,-1); }}

  if (!notinit(dbp1_scond)) { setsv(dbp1,SCOND,4,dbp1_scond); }
  
   if (!notinit(dbp1_svnoise)) { setsv(dbp1,SVNOISE,4,dbp1_svnoise); }
  
  if (notinit(ax_dia_factor)) ax_dia_factor = 1; /* multiplier for axon diameter */
  if (notinit(dvrev)) dvrev = -0.045;           /* Vrev for dens_dbp1.n */
  if (notinit(dvst))   dvst = -0.07;            /* Vstart for dens_dbp1.n */

  if (notinit(ivplot)) ivplot = 0;              /* make I/V plot */
  if (notinit(outward)) outward = 0;            /* >0 => calc K cond, otherwise Na */
  if (notinit(iscal)) iscal   = 3e-9;           /* plot scale */
  if (notinit(gscal)) gscal   = 1.5e-8;         /* plot scale */
  if (notinit(am_rm)) am_rm   = 10e-3;          /* amac RM, used in dens_amaca.n */
  if (notinit(am_vs)) am_vs   = -0.07;          /* amac vstart, used in dens_amaca.n */
  if (notinit(am_vr)) am_vr   = -0.07;          /* amac vrev,   used in dens_amaca.n */

  if (notinit(nadist)) nadist   = 10e-3;        /* distal Na density, dens_amaca.n */
  if (notinit(namid))   namid   = 5e-3;         /* middle Na density, dens_amaca.n */
  if (notinit(kdist))   kdist   = 3e-3;         /* dist Kdr density, dens_amaca.n */
  if (notinit(kmid))     kmid   = 2e-3;         /* mid  Kdr density, dens_amaca.n */
  if (notinit(kprox))   kprox   = 2e-3;         /* prox Kdr density, dens_amaca.n */
  if (notinit(ksoma))   ksoma   = 2e-3;         /* soma Kdr density, dens_amaca.n */
  if (notinit(cadist)) cadist   = 3e-3;         /* dist Ca density, dens_amaca2.n */
  if (notinit(camid))   camid   = 0e-3;         /* mid  Ca density, dens_amaca2.n */

  if (!notinit(set_drm)) {                      /* user set default Rm */
       setn(ct,NRM,set_drm);                    /* set default Rm */
       drm = set_drm;
  }
  if (notinit(am_rid)) am_rid = dri;
  if (notinit(am_rii)) am_rii = dri;
  if (notinit(am_rip)) am_rip = dri;
  if (notinit(kdr_cond)) kdr_cond=10e-3;
  if (notinit(ka_cond)) ka_cond=20e-3;
  if (notinit(na_cond)) na_cond=50e-3;
  if (notinit(kdr_condp)) kdr_condp=1e-3;
  if (notinit(ka_condp)) ka_condp=0;
  if (notinit(na_condp)) na_condp=1e-3;
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
    // if (!notinit(setkdr_cond))    celdens[am][_KDR][DENDD] = setkdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    
    // int cn;
    // ndens[amac][cn=1] = 0;           // set cn 1 to use amac_densfile
    // ndens[amac][cn=2] = 1;  		// set cn 2 to use amac_densfile2
    // ndens[amac][cn=3] = 0;  		// set cn 3 to use amac_densfile
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
    if (amacpair) current2 = i(ndn(ct, 2, soma));
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


/*------------------------------------------------------*/

void addlabels(void)

{
        //label (ndn(amac,2,96), red);
        //label (ndn(amac,3,2726), red);
	label (findsynloc(am,2,am,1,100, 0, 0.0, 50), red);
	label (findsynloc(am,3,am,1,100, 0, 0.0, 50), red);
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
    double start;
    //double xloc,yloc;
    
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

  if (notinit(stimdur))      stimdur    = 0.05;
  if (notinit(poststimdur))      poststimdur    = 0.05;
  if (notinit(minten))       minten = -0.050; //mean intensity
  if (notinit(scontrast)) scontrast = 0.010; 
  if (notinit(stimtime))   stimtime = .01;
  if (notinit(sdia))           sdia = 1000;   
  if (notinit(xloc))           xloc=0;
  if (notinit(yloc))	       yloc=0;
  if (notinit(sdia2))           sdia2=0;   
  if (notinit(xloc2))           xloc2=0;
  if (notinit(yloc2))	       yloc2=0;
  if (notinit(sdia3))           sdia3=0;   
  if (notinit(xloc3))           xloc3=0;
  if (notinit(yloc3))	       yloc3=0;
  if (notinit(sdia4))           sdia4=0;   
  if (notinit(xloc4))           xloc4=0;
  if (notinit(yloc4))	       yloc4=0;
if (notinit(stimtype))   stimtype = 1;      

  // dtrial = prestimdur+stimdur+tailcurdur+poststimdur;
//  dtrial = prestimdur+dstim+poststimdur+0.1;
//  endexp  = dtrial;
  ploti = 1e-4;


  // midcbp  = findmid(ct,0,0);

 // nsynap = synapse_add (1,dbp1,-1,-1,am,1);		/* make list of bipolar synapses onto amac 1 for recording below */
 // nsynape = nsynapi = 0;
 // for (i=2; i<=n_am; i++) {
 //   nsynape = synapse_add (3,am,i,-1,am,1,3);	/* make list of excit synapses from amacs onto amac1 */
 //   nsynap += nsynape;
//    fprintf (stderr,"# nsynape %d\n",nsynape);
 // }
 // for (i=2; i<=n_am; i++) {
 //   nsynapi = synapse_add (4,am,i,-1,am,1,4);	/* make list of inhib synapses from amacs onto amac1 */
  //  nsynap += nsynapi;
   // fprintf (stderr,"# nsynapi %d\n",nsynapi);
  //}
 // fprintf (stderr,"# nsynap %d\n",nsynap);

  if (make_movie) {
       setonplot(onplot_movie); /* set movie plot routine */
       if (space_time) {  /* movie */
         plot_v_nod(ct=am,cn=1,soma,Vming,Vmaxg,1,"Vsoma1",pl=10,0.35);
         plot_v_nod(ct=am,cn=2,soma,Vming,Vmaxg,2,"Vsoma2",pl=10,0.35);
	 plot_v_nod(ct=am,cn=3,soma,Vming,Vmaxg,3,"Vsoma3",pl=10,0.35);
	 plot_v_nod(ct=am,cn=4,soma,Vming,Vmaxg,4,"Vsoma4",pl=10,0.35);
// 	 plot_v_nod(ct=am,cn=5,soma,Vming,Vmaxg,5,"Vsoma5",pl=10,0.35);
// 	 plot_v_nod(ct=am,cn=6,soma,Vming,Vmaxg,6,"Vsoma6",pl=10,0.35);
// 	 plot_v_nod(ct=am,cn=7,soma,Vming,Vmaxg,7,"Vsoma7",pl=10,0.35);
// 	 plot_v_nod(ct=am,cn=8,soma,Vming,Vmaxg,8,"Vsoma8",pl=10,0.35);
// 	 plot_v_nod(ct=am,cn=9,soma,Vming,Vmaxg,9,"Vsoma9",pl=10,0.35);

         //plot_v_nod(ct=dsgc,cn,1336,Vming,Vmaxg,c=green,"Vtip1",pl=10,0.35); // for morph_ds1e
         //plot_v_nod(ct=dsgc,cn,582,Vming,Vmaxg,c=red,"Vtip2",pl=10,0.35);
         //plot_na_inact(dsgc, 1, soma, red, "Na[soma]", 12, .35);
	 //
       };
  }
  else {

//          plot_v_nod(dbp1,findmida(dbp1,0,0),  soma,Vmin=-0.05, Vmax= -0.03, 1, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,20,0), soma,Vmin=-0.05, Vmax= -0.03, 2, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,40,0), soma,Vmin=-0.05, Vmax= -0.03, 3, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,60,0), soma,Vmin=-0.05, Vmax= -0.03, 4, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,80,0), soma,Vmin=-0.05, Vmax= -0.03, 5, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,100,0),soma,Vmin=-0.05, Vmax= -0.03, 6, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,120,0),soma,Vmin=-0.05, Vmax= -0.03, 7, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,0,-155),soma,Vmin=-0.05, Vmax= -0.03, 8, "", 22, 0.3);
//  	plot_v_nod(dbp1,findmida(dbp1,30,-10),soma,Vmin=-0.05, Vmax= -0.03, 9, "", 22, 0.3);
//  
//  
// 	plot_syncond(findsynloc(dbp1,0,0),   cmin=0,cmax=100e-12, 1, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,20,0),  cmin=0,cmax=100e-12, 2, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,40,0),  cmin=0,cmax=100e-12, 3, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,60,0),  cmin=0,cmax=100e-12, 4, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,80,0),  cmin=0,cmax=100e-12, 5, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,100,0), cmin=0,cmax=100e-12, 6, 20,"",0.3); 
//    	plot_syncond(findsynloc(dbp1,120,0), cmin=0,cmax=100e-12, 7, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,0,-155), cmin=0,cmax=100e-12, 8, 20,"",0.3); 
//  	plot_syncond(findsynloc(dbp1,30,-10), cmin=0,cmax=100e-12, 9, 20,"",0.3); 
//  	
	// plot_syncond(findsynlocr(am,1,80, 0, 0.0), cmin=0,cmax=100e-12, 3, 18,"",0.3); 
	// plot_syncond(findsynlocr(am,1,120,0, 0.0), cmin=0,cmax=100e-12, 6, 18,"",0.3); 



	if (set_vclamp > 0) { Vmin=vstart; if (stimtype==4) Vmax=vstart+0.01; else Vmax = 0; }
	else                { Vmin = -0.07; Vmax = -0.00; }


		
	
	//plot_v_nod(ct, cn, soma, 		    Vmin, Vmax, colr=blue, "", 10, 0.5);
	
 	plot_v_nod(ct, cn, soma, 		    Vmin, Vmax, colr=blue, "", 15, 0.5);
 	plot_v_nod(ct, cn+1, soma,		    Vmin, Vmax, colr=green,"", 15, 0.5);
// 	
	
//  	//1
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,30,-57),  Vmin, Vmax, colr=green,"", 14, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,50,-105),  Vmin, Vmax, colr=brown,"", 14, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,75,-155), Vmin, Vmax, colr=red,  "", 14, 0.5);
//  	
//  	//2
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,55,30),  Vmin, Vmax, colr=green,"", 13, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,105,45),  Vmin, Vmax, colr=brown,"", 13, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,160,60), Vmin, Vmax, colr=red,  "", 13, 0.5);
//  	
//  	//3
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,55,0),  Vmin, Vmax, colr=green,"", 12, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,118,0),  Vmin, Vmax, colr=brown,"", 12, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,175,0), Vmin, Vmax, colr=red,  "", 12, 0.5);
//  	
//  	//4
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,0,-53.5),  Vmin, Vmax, colr=green,"", 11, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,0,-103),  Vmin, Vmax, colr=brown,"", 11, 0.5);
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,0,-155), Vmin, Vmax, colr=red,  "", 11, 0.5);
// 
//  	//5
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,45,-45),  Vmin, Vmax, colr=green,"", 10, 0.5);
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,82,-82),  Vmin, Vmax, colr=brown,"", 10, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,122,-122), Vmin, Vmax, colr=red,  "", 10, 0.5);
//  	
//  	//6
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-63,0),  Vmin, Vmax, colr=green,"", 9, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-113,0),  Vmin, Vmax, colr=brown,"", 9, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-167,-30), Vmin, Vmax, colr=red,  "", 9, 0.5);
//  	
//  	//7
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-45,-45),  Vmin, Vmax, colr=green,"", 8, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-85,-85),  Vmin, Vmax, colr=brown,"", 8, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-145,-90), Vmin, Vmax, colr=red,  "", 8, 0.5);
//  	
// 	//8
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,-42,52),  Vmin, Vmax, colr=green,"", 7, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-102,52),  Vmin, Vmax, colr=brown,"", 7, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,-154,74), Vmin, Vmax, colr=red,  "", 7, 0.5);
//  
//  	//9
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,30,-10),  Vmin, Vmax, colr=green,"", 6, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,80,-30),  Vmin, Vmax, colr=brown,"", 6, 0.5);
//  	plot_v_nod(ct, cn, findnodloc(ct,cn,132,-52), Vmin, Vmax, colr=red,  "", 6, 0.5);
// 	
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,30,-10),  Vmin, Vmax, colr=green,"", 6, 0.5);
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,-60, 0),  Vmin, Vmax, colr=red,"", 6, 0.5);
// 	
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,80, -55),  Vmin, Vmax, colr=green,"", 7, 0.5);
// 	plot_v_nod(ct, cn, findnodloc(ct,cn,80, -35),  Vmin, Vmax, colr=red,"", 7, 0.5);
// 	plot_v_nod(ct, cn, 27,  Vmin, Vmax, colr=brown,"", 7, 0.5);
// 	plot_v_nod(ct, cn, 276,  Vmin, Vmax, colr=magenta,"", 7, 0.5);
// 	
// 	plot_v_nod(ct, 2, 282, Vmin, Vmax, colr=green,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 219, Vmin, Vmax, colr=brown,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 114, Vmin, Vmax, colr=magenta,"", 5, 0.5);
// 	
// 	plot_v_nod(ct, 2, findnodloc(ct,cn,-280, 0),  Vmin, Vmax, colr=red,"", 5, 0.5);
	
	
	//plot_v_nod(ct, 1, 304,  		    Vmin, Vmax, colr=brown, "", 13, 0.5);
	//plot_v_nod(ct, 3, 2726, 		    Vmin, Vmax, colr=magenta, "", 13, 0.5);
	
//each dendrite end	
// 	plot_v_nod(ct, 1, 297, Vmin, Vmax, colr=green,"", 6, 0.5); 
// 	plot_v_nod(ct, 1, 256, Vmin, Vmax, colr=red,  "", 6, 0.5); 
// 	
// 	plot_v_nod(ct, 1, 207, Vmin, Vmax, colr=blue, "", 7, 0.5); 
// 	plot_v_nod(ct, 1, 159, Vmin, Vmax, colr=brown,"", 7, 0.5); 
// 	plot_v_nod(ct, 1, 102, Vmin, Vmax, colr=cyan, "", 7, 0.5); 
// 	plot_v_nod(ct, 1, 52,  Vmin, Vmax, colr=8, "", 7, 0.5); 
	
// 	plot_v_nod(ct, 1, 351, Vmin, Vmax, colr=yellow,"", 8, 0.5); 
// 	plot_v_nod(ct, 1, 409, Vmin, Vmax, colr=white,"", 8, 0.5); 
//	plot_v_nod(ct, 1, 461, Vmin, Vmax, colr=magenta,  "", 8, 0.5); 
	
	
// 	plot_v_nod(ct, 2, 297, Vmin, Vmax, colr=green,"", 3, 0.5); 
// 	plot_v_nod(ct, 2, 256, Vmin, Vmax, colr=red,  "", 3, 0.5); 
// 	
// 	plot_v_nod(ct, 2, 207, Vmin, Vmax, colr=blue, "", 4, 0.5); 
// 	plot_v_nod(ct, 2, 159, Vmin, Vmax, colr=brown,"", 4, 0.5); 
// 	plot_v_nod(ct, 2, 102, Vmin, Vmax, colr=cyan, "", 4, 0.5); 
// 	plot_v_nod(ct, 2, 52,  Vmin, Vmax, colr=8, "", 4, 0.5); 
// 	
// 	plot_v_nod(ct, 2, 351, Vmin, Vmax, colr=yellow,"", 5, 0.5); 
// 	plot_v_nod(ct, 2, 409, Vmin, Vmax, colr=white,"", 5, 0.5); 
// 	plot_v_nod(ct, 2, 461, Vmin, Vmax, colr=magenta,  "", 5, 0.5); 
// 	
//along dend under stim 
	/*plot_v_nod(ct, 1, 159, Vmin, Vmax, colr=green,"", 6, 0.5);
	plot_v_nod(ct, 1, 118, Vmin, Vmax, colr=red  ,"", 6, 0.5);
	plot_v_nod(ct, 1, 107, Vmin, Vmax, colr=blue,"", 6, 0.5);
	plot_v_nod(ct, 1, 104, Vmin, Vmax, colr=yellow,"", 6, 0.5);
	
	plot_v_nod(ct, 1, 102, Vmin, Vmax, colr=green,"", 7, 0.5);
	plot_v_nod(ct, 1, 64 , Vmin, Vmax, colr=red,"", 7, 0.5);
	plot_v_nod(ct, 1, 58, Vmin, Vmax, colr= blue,"", 7, 0.5);
	plot_v_nod(ct, 1, 53, Vmin, Vmax, colr=yellow,"", 7, 0.5);
	
	plot_v_nod(ct, 1, 52, Vmin, Vmax, colr=green,"", 8, 0.5);
	plot_v_nod(ct, 1, 15, Vmin, Vmax, colr=red,"", 8, 0.5);
	plot_v_nod(ct, 1, 6, Vmin, Vmax, colr=blue,"", 8, 0.5);
	plot_v_nod(ct, 1, 2, Vmin, Vmax, colr=yellow,"", 8, 0.5);
	
	
	
	
	plot_v_nod(ct, 2, 256, Vmin, Vmax, colr=green,"", 5, 0.5);
	plot_v_nod(ct, 2, 235, Vmin, Vmax, colr=red  ,"", 5, 0.5);
	plot_v_nod(ct, 2, 216, Vmin, Vmax, colr=blue,"", 5, 0.5);
	plot_v_nod(ct, 2, 210, Vmin, Vmax, colr=yellow,"",5, 0.5);

	plot_v_nod(ct, 2, 207, Vmin, Vmax, colr=green,"", 4, 0.5);
	plot_v_nod(ct, 2, 184 , Vmin, Vmax, colr=red,"", 4, 0.5);
	plot_v_nod(ct, 2, 165, Vmin, Vmax, colr= blue,"", 4, 0.5);
	plot_v_nod(ct, 2, 161, Vmin, Vmax, colr=yellow,"", 4, 0.5);
	
	plot_v_nod(ct, 2, 159, Vmin, Vmax, colr=green,"", 3, 0.5);
	plot_v_nod(ct, 2, 129, Vmin, Vmax, colr=red,"", 3, 0.5);
	plot_v_nod(ct, 2, 108, Vmin, Vmax, colr=blue,"", 3, 0.5);
	plot_v_nod(ct, 2, 104, Vmin, Vmax, colr=yellow,"", 3, 0.5);*/
	
//on dend between stim area and tip
// 	plot_v_nod(ct, 2, 256, Vmin, Vmax, colr=green,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 250 , Vmin, Vmax, colr=blue ,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 245 , Vmin, Vmax, colr=yellow,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 230 , Vmin, Vmax, colr=cyan ,"", 5, 0.5);
// 	plot_v_nod(ct, 2, 235, Vmin, Vmax, colr=red  ,"", 5, 0.5);
// 	
// 	plot_v_nod(ct, 2, 207, Vmin, Vmax, colr=green,"", 4, 0.5);
// 	plot_v_nod(ct, 2, 201, Vmin, Vmax, colr=blue ,"", 4, 0.5);
// 	plot_v_nod(ct, 2, 195 , Vmin, Vmax, colr=yellow,"", 4, 0.5);
// 	plot_v_nod(ct, 2, 189 , Vmin, Vmax, colr=cyan ,"", 4, 0.5);
// 	plot_v_nod(ct, 2, 184 , Vmin, Vmax, colr=red,"", 4, 0.5);	
// 	
//  	plot_v_nod(ct, 2, 159, Vmin, Vmax, colr=green,"", 3, 0.5);
//  	plot_v_nod(ct, 2, 152 , Vmin, Vmax, colr=blue ,"", 3, 0.50);
//  	plot_v_nod(ct, 2, 144 , Vmin, Vmax, colr=yellow ,"", 3, 0.5);
//  	plot_v_nod(ct, 2, 136 , Vmin, Vmax, colr=cyan ,"", 3, 0.5);	
//  	plot_v_nod(ct, 2, 129, Vmin, Vmax, colr=red,"", 3, 0.5);
// 	
// 	
 	//plot_v_nod(ct, 1, 153, Vmin, Vmax, colr=green,"", 6, 0.5);
 	//plot_v_nod(ct, 1, 150 , Vmin, Vmax, colr=blue ,"", 6, 0.5);
 	//plot_v_nod(ct, 1, 145 , Vmin, Vmax, colr=yellow ,"", 6, 0.5);
	//plot_v_nod(ct, 1, 142 , Vmin, Vmax, colr=magenta ,"", 6, 0.5);
// 	plot_v_nod(ct, 1, 139 , Vmin, Vmax, colr=8 ,"", 6, 0.5);
// 	plot_v_nod(ct, 1, 136 , Vmin, Vmax, colr=brown ,"", 6, 0.5);
// 	plot_v_nod(ct, 1, 133 , Vmin, Vmax, colr=white ,"", 6, 0.5);
//  	plot_v_nod(ct, 1, 130 , Vmin, Vmax, colr=cyan ,"", 6, 0.5);
//  	plot_v_nod(ct, 1, 122, Vmin, Vmax, colr=red  ,"", 6, 0.5);
// 	
// 	plot_v_nod(ct, 1, 102, Vmin, Vmax, colr=green,"", 7, 0.5);


	
//  	plot_v_nod(ct, 1, 72 , Vmin, Vmax, colr=red,"", 1, 0.5);
// 	plot_v_nod(ct, 1, 75 , Vmin, Vmax, colr=7 ,"", 2, 0.5);
// 	plot_v_nod(ct, 1, 78, Vmin, Vmax, colr=yellow ,"", 3, 0.5);
// 	plot_v_nod(ct, 1, 81 , Vmin, Vmax, colr=green ,"", 4, 0.5);
// 	plot_v_nod(ct, 1, 83 , Vmin, Vmax, colr=blue ,"",5, 0.5);
// 	plot_v_nod(ct, 1, 85 , Vmin, Vmax, colr=magenta ,"", 6, 0.5);
// 	plot_v_nod(ct, 1, 87 , Vmin, Vmax, colr=brown ,"", 7, 0.5);
// 	plot_v_nod(ct, 1, 89 , Vmin, Vmax, colr=white ,"", 8, 0.5);
// 	plot_v_nod(ct, 1, 93 , Vmin, Vmax, colr=cyan ,"", 9, 0.5);
// 	plot_v_nod(ct, 1, 95 , Vmin, Vmax, colr=8 ,"", 10, 0.5);
// 	plot_v_nod(ct, 1, 102 , Vmin, Vmax, colr=9 ,"", 11, 0.5);
// 	
// 	plot_v_nod(ct, 1, 52, Vmin, Vmax, colr=green,"", 8, 0.5);
// 	plot_v_nod(ct, 1, 46, Vmin, Vmax, colr=blue ,"", 8, 0.5);
// 	plot_v_nod(ct, 1, 35, Vmin, Vmax, colr=yellow ,"", 8, 0.5);
// 	plot_v_nod(ct, 1, 25, Vmin, Vmax, colr=cyan ,"", 8, 0.5);
// 	plot_v_nod(ct, 1, 15, Vmin, Vmax, colr=red,"", 8, 0.5);
	
	
	
    // plot_func(isyn_tot,1,imax,imin);  	    plot_param("Itotsyn",green,2,0.3);
   //  plot_func(gsyn_tot,1,gmax,0);     	    plot_param("Gtotbp",blue, 1,0.3);
    // plot_func(gsyn_tot,3,gmax*0.25,0);     plot_param("Gtotame",cyan, 1,0.3);
    // plot_func(gsyn_tot,4,gmax*0.25,0);     plot_param("Gtotami",red, 1,0.3);

	//loop with one stim side
   	plot_v_nod(ct, 1, 414, Vmin, Vmax, colr=cyan,"", 4, 0.5);
   	plot_v_nod(ct, 1, 421, Vmin, Vmax, colr=red,"", 4, 0.5);
  	plot_v_nod(ct, 1, 425, Vmin, Vmax, colr=white,"", 4, 0.5);
  	plot_v_nod(ct, 1, 440, Vmin, Vmax, colr=magenta,"", 4, 0.5);
  	plot_v_nod(ct, 1, 454, Vmin, Vmax, colr=blue,"", 4, 0.5);
  	plot_v_nod(ct, 1, 461, Vmin, Vmax, colr=green,"", 4, 0.5);

	plot_v_nod(ct, 1, 359, Vmin, Vmax, colr=cyan,"", 3, 0.5);
	plot_v_nod(ct, 1, 380, Vmin, Vmax, colr=red,"", 3, 0.5);
	plot_v_nod(ct, 1, 387, Vmin, Vmax, colr=white,"", 3, 0.5);
	plot_v_nod(ct, 1, 398, Vmin, Vmax, colr=magenta,"", 3, 0.5);
	plot_v_nod(ct, 1, 405, Vmin, Vmax, colr=blue ,"", 3, 0.5);
	plot_v_nod(ct, 1, 409, Vmin, Vmax, colr=green,"", 3, 0.5);
	
	plot_v_nod(ct, 2, 414, Vmin, Vmax, colr=cyan,"", 2, 0.5);
	plot_v_nod(ct, 2, 433, Vmin, Vmax, colr=red,"", 2, 0.5);
	plot_v_nod(ct, 2, 442, Vmin, Vmax, colr=white,"", 2, 0.5);
	plot_v_nod(ct, 2, 449, Vmin, Vmax, colr=magenta,"", 2, 0.5);
	plot_v_nod(ct, 2, 455, Vmin, Vmax, colr=blue ,"", 2, 0.5);
	plot_v_nod(ct, 2, 461, Vmin, Vmax, colr=green,"", 2, 0.5);

 	plot_v_nod(ct, 2, 359, Vmin, Vmax, colr=cyan,"", 1, 0.5);
 	plot_v_nod(ct, 2, 363, Vmin, Vmax, colr=red,"", 1, 0.5);
 	plot_v_nod(ct, 2, 367, Vmin, Vmax, colr=white,"", 1, 0.5);
 	plot_v_nod(ct, 2, 385, Vmin, Vmax, colr=magenta,"", 1, 0.5);
 	plot_v_nod(ct, 2, 398, Vmin, Vmax, colr=blue ,"", 1, 0.5);
 	plot_v_nod(ct, 2, 409, Vmin, Vmax, colr=green,"", 1, 0.5);
	
  }


  
    if (notinit(stimtyp))    stimtyp=1;
    if (notinit(velocity))   velocity  = 1000;
    if (notinit(barwidth))   barwidth  = 100;
    if (notinit(barlength))  barlength = 200;
    if (notinit(theta))         theta  = 0;
    if (notinit(direction))  direction = 1;	
    if (notinit(waveshape))  waveshape = 1;		  // 1 -> square wave

    if (notinit(setxmin))      setxmin = 0;               // set plot to start at 0
    if (notinit(predur))        predur = 0.05;
    if (notinit(speriod))	speriod=50;
    if (notinit(sphase))	sphase=0;
    if (notinit(orient))	orient=0;
    if (notinit(tfreq))		tfreq=4;
    if (notinit(drift))		drift=1;
      

    simtime = -predur;					  // must be set ahead of stim_backgr()
   // setxmin = 0;

    stim_backgr(minten);
    
    if(!make_movie){
      if (disp==16) {			// display the stimulus
	double t;

    //  stimdur = move_stim(stimtime, barwidth, theta, velocity, ncycles, scontrast, direction, mask=1);
      if (stimtyp==1){stim_spot (sdia, xloc, yloc, scontrast, start=stimtime, stimdur);}
      
      else if(stimtyp==2){ stim_sine (speriod, sphase, orient, xloc, yloc, tfreq, drift,
		1, scontrast, 0, start=stimtime,  stimdur);}
      else if(stimtyp==3){stim_spot(sdia, xloc, yloc, scontrast, start=stimtime, stimdur);
			  stim_spot(sdia2, xloc2, yloc2, scontrast, start=stimtime, stimdur);
			  stim_spot(sdia3, xloc3, yloc3, scontrast, start=stimtime, stimdur);
			  stim_spot(sdia4, xloc4, yloc4, scontrast, start=stimtime, stimdur);
	
      } 
      //stim node 73;
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

//    if (notinit(set_vclamp)) set_vclamp = 0;
//    if (set_vclamp > 0) {	
//      vclamp              (ndn(ct, cn,soma),   vhold, simtime,  2);
//      if (amacpair) vclamp(ndn(ct, cn+1,soma), vhold, simtime,  2);
//    }
	

    Gmax = 0;
//	stim_spot(100, 100, 0, scontrast, 0.001, 0.2);
//	stim_spot(1000, 50, 0, scontrast, stimtime, stimdur);



  endexp=stimtime+stimdur+poststimdur;
		 /* turn on  background */
  step(predur);

  set_run_on_exit(runonexit);                         // set to erase savefile on ^C
  sprintf (savefile,"amac_vclamp%06d",getpid());       // add pid to file name


 if (stimtyp==1){stim_spot (sdia, xloc, yloc, scontrast, start=stimtime, stimdur);}
 else if(stimtyp==2){ stim_sine (speriod, sphase, orient, xloc, yloc, tfreq, drift,
		1, scontrast,0, start=stimtime,  stimdur);}
 else if (stimtyp==3){stim_spot(sdia, xloc, yloc, scontrast, start=stimtime, stimdur);
		      stim_spot(sdia2, xloc2, yloc2, scontrast, start=stimtime, stimdur);
		      stim_spot(sdia3, xloc3, yloc3, scontrast, start=stimtime, stimdur);
		      stim_spot(sdia4, xloc4, yloc4, scontrast, start=stimtime, stimdur);
   
}
	step(stimdur+stimtime+poststimdur);

	
   /* no vclamp */

  //unlink (savefile);
  //savefile[0] = 0;
}
