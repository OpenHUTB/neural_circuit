/* Experiment dsgc_sbac for retsim */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include "ncfuncs.h"
#include "ncio.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"

#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int n_dsgc;
extern int rec_ct;
// extern long int total_lookups;
// extern long int total_finds;


double theta;
double iroff;
int light_inhib;
int dsgc_prefdir;
int node_dist;
int set_vel;

int sbarr;
int rec_ct;
int rec_cn;
int movein;
int stimtype;
int stim_theta;
int rstim_theta;
int direction;
int swaveshape;
int twaveshape;
int makenv;
int use_stimfile;
int run_vclamp;
int run_vclamp_sbac;
int sbac_species;
int revdir;
int c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
int c11, c12, c13, c14, c15, c16, c17, c18;
int c19, c20, c21, c22, c23, c24, c25, c50, c53;

int sbr;
int sbac_color;
int dsgc_color;
int draw_synapse;
int dbp1_dsgc_reg;
int ivplot;
int outward;
int ncycles;

int sreg1;
int sreg2;
int sreg3;
int sreg4;
int sreg5;
int sreg6;
int stonic;
int nstim;
int db1_morph;
int db1_biophys;
int db2_morph;
int db2_biophys;
int t2bipolar;

int plotnod1;
int plotnod2;
int plotnod3;
int plotnod4;

double set_ploti;
double Vmin_c1;
double idiff;
double current;
double current2;
double voltage;
double set_sr;

double iscal;
double gscal;
double vstart;
double vstop;
double vstep;
double tailvolt;
double istart;
double istop;
double istep;
double disp_stim_incr;
double disp_stim_max;
double disp_stim_min;

double min_sp;
double max_sp;
double incr_sp;

double sbac_dia;
double sbac_soma_z;
double sbac_soma_z2;
double sbac_dens;
double sbspac;
double sbac_regu;
double sb_mglur_maxdist;
double sbsynang;
double idist;
double sbac_r;
double sbac_rnd;
double dsgc_x;
double dsgc_y;
double rloc;
double dbpthr;
double dbp2thr;
double sb_denddia;
double sb_denddia2;
double sb_arbscale;
double sb_arbscale2;
double sb_dia1;
double sb_ca6_offm;
double sb_ca6_tauc;
double sb_ca6_taud;
double dbp1_dens;
double dbp2_dens;
double dbp1_maxsdist;
double dbp2_maxsdist;
double db1_r5;
double db2_r5;
double db1_r6;
double db2_r6;
double dsgc_denddia;
double db1_lam;
double db2_lam;
double db1_lax;
double db2_lax;
double db1_sca;
double db2_sca;
double db_vs;

double center_weight;
double surr_weight;
double blur_csize;
double blur_ssize;
double surr_delay;
double surr_tau;
double stim_tstep;

double g_am_dsgc;
double g_am2_dsgc;
double vhold;
double sbac_vhold;
double sbac_vpulse;
double sbac_vpulse_dur;
double sbac_vpulse_time;
double soma_clamp_time;
double g_sbac_dsgc;
double g_sbac_sbac;
double n_sbac_sbac;
double m_sbac_sbac;
double r_sbac_sbac;

double na5offs;
double k1offm;
double k6offm;

double dbp1_mrate;
double dbp1_mpool;
double dbp1_hdur;
double dbp1_hfilt;
double dbp1_hgain;
double dbp1_hoffs;
double dbp2_mrate;
double dbp2_mpool;
double g_dbp1_sbac;
double g_dbp1_ams;
double g_dbp1_dsgc;
double n_dbp1_sbac;
double v_dbp1_sbac;
double g_dbp2_dsgc;
double g_dbp2_sbac;
double n_dbp2_sbac;
double v_dbp2_sbac;
double mglur2;
double g_ams_sbac;
double f_ams_sbac;
double n_ams_sbac;
double ams_synanpi;
double ams_synanpo;
double g_sbac_dbp1;
double n_dbp1_ams;
double sb_db_anni;
double sb_db_dur;
double sb_db_casens;
double sb_db_synr;

double naax;
double nahd;
double naab;
double kax;
double ca5ax;
double ca7ax;
double ca6ax;
double ca1ax;
double db1_ca1ax;
double db1_ca6ax;
double db2_ca1ax;
double db2_ca6ax;
double db1_cap;
double db2_cap;
double db1_axdia;
double db2_axdia;
double db1_ca1_offm;
double db1_ca6_offm;
double db2_ca1_offm;
double db2_ca6_offm;
double db1_ca7_offm;
double db2_ca7_offm;
double db2_ca6_tauc;
double db2_ca6_taud;

double barwidth;
double barlength;
double tfreq;
double minten;
double scontrast;
double scontrastn;
double stimtime;
double velocity;
double stimx;
double stimy;
double annrad;
double disptime;
double stimdur;
double prestimdur;
double poststimdur;
double tailcurdur;
double noise_dur;
double spotdur;
double spotint;

double sblur;
double ioffset;
double istim;
double sbac_istim;
double predur;
double spotdia;
double spotdia_incr;
double spotdia_max;
double annulus_width;
double annulus_odia;
double barwidth_incr;
double barwidth_max;
double sdia;
double spdia;
double smdia;
double sndia;
double sddia;
double orad1;
double irad1;
double irad2;
double orad2;
double mask_dia;
double mask_x;
double mask_y;
double c1_somax;
double c1_somay;

double sbaclm;
double sb1mul;
double dsgclm;
double dsomadia;
double dtreedia;

double wac_g;

double sbac_synrng;
double sbac_synanp;

double dbp1_anpi;
double dbp1_anpo;
double dbp2_anpi;
double dbp2_anpo;

double sbac_maxsdist;
double sbac_synanni;
double sbac_isynanpi;
double sbac_isynanpo;
double sbac_isynanni;
double sbac_isynrngi;
double sbac_ithr;
double sbac_synspac;

double sb_rii;
double sb_rip;
double sb_rid;
double sb_rm;
double sb_rmp;
double sb_rmd;
double sb_vr;
double sb_vs;
double sb_cap;
double sb_capm;
double sb_capp;

double dendrm;

double amna;
double amk;
double amca;
double am2ca;
double amsca;
double amvst;
double amvrev;
double amrm;

double nasoma;
double naprox;
double namid;
double nadist;

double kprox;
double kdist;
double kmid;
double ksoma;

double kdrs;
double kdrp;
double kdrm;
double kdrd;

double camid;
double cadist;
double catsoma;
double catprox;
double catmid;
double catdist;

double cadiap;
double cadiam;
double cadiad;

const char *stim_fnum;

char savefile[30] = {0};

void sb_init(void);

/*--------------------------------------------------------*/

void defparams(void)

{
  defparams_dsgc_movie();
  defparams_onplot_movie();

  setptr("set_ploti", 	 &set_ploti);
  setptr("ivplot", 	 &ivplot);
  setptr("outward", 	 &outward);
  setptr("set_vel", 	 &set_vel);
  setptr("nstim", 	 &nstim);
  setptr("db1_morph", 	 &db1_morph);
  setptr("db1_biophys",	 &db1_biophys);
  setptr("db1_sca",	 &db1_sca);
  setptr("db2_morph", 	 &db2_morph);
  setptr("db2_biophys",	 &db2_biophys);
  setptr("db2_sca",	 &db2_sca);
  setptr("db1_ca1_offm", &db1_ca1_offm);
  setptr("db1_ca6_offm", &db1_ca6_offm);
  setptr("db2_ca1_offm", &db2_ca1_offm);
  setptr("db2_ca6_offm", &db2_ca6_offm);
  setptr("db1_ca7_offm", &db1_ca7_offm);
  setptr("db2_ca7_offm", &db2_ca7_offm);
  setptr("db2_ca6_tauc", &db2_ca6_tauc);
  setptr("db2_ca6_taud", &db2_ca6_taud);
  setptr("t2bipolar",    &t2bipolar);
  setptr("set_sr",       &set_sr);

  setptr("plotnod1",     &plotnod1);
  setptr("plotnod2",     &plotnod2);
  setptr("plotnod3",     &plotnod3);
  setptr("plotnod4",     &plotnod4);

  setptr("iscal", 	 &iscal);
  setptr("gscal", 	 &gscal);
  setptr("vstart", 	 &vstart);
  setptr("vstop", 	 &vstop);
  setptr("vstep", 	 &vstep);
  setptr("tailvolt", 	 &tailvolt);

  setptr("istart", 	 &istart);
  setptr("istop", 	 &istop);
  setptr("istep", 	 &istep);

  setptr("theta", 	 &theta);
  setptr("iroff", 	 &iroff);
  setptr("light_inhib",  &light_inhib);
  setptr("dsgc_prefdir", &dsgc_prefdir);
  setptr("node_dist",    &node_dist);
  setptr("makenv",	 &makenv);
  setptr("direction",	 &direction);
  setptr("swaveshape",	 &swaveshape);
  setptr("twaveshape",	 &twaveshape);
  setptr("revdir",	 &revdir);

  setptr("min_sp",	 &min_sp);
  setptr("max_sp",	 &max_sp);
  setptr("incr_sp",	 &incr_sp);

  setptr("sbspac",	 &sbspac);
  setptr("rec_ct",	 &rec_ct);
  setptr("rec_cn",	 &rec_cn);
  setptr("sbarr",	 &sbarr);
  setptr("stimtype",	 &stimtype);
  setptr("use_stimfile", &use_stimfile);
  setptr("run_vclamp",   &run_vclamp);
  setptr("run_vclamp_sbac",&run_vclamp_sbac);
  setptr("sbac_vhold",   &sbac_vhold);
  setptr("sbac_vpulse",  &sbac_vpulse);
  setptr("sbac_vpulse_dur",&sbac_vpulse_dur);
  setptr("sbac_vpulse_time",&sbac_vpulse_time);
  setptr("soma_clamp_time",&soma_clamp_time);
  setptr("sbac_species", &sbac_species);
  setptr("sbac_r",	 &sbac_r);
  setptr("sbac_rnd",	 &sbac_rnd);
  setptr("dsgc_x",	 &dsgc_x);
  setptr("dsgc_y",	 &dsgc_y);
  setptr("rloc",	 &rloc);
  setptr("dbpthr",	 &dbpthr);
  setptr("dbp2thr",	 &dbp2thr);
  setptr("sb_denddia",	 &sb_denddia);
  setptr("sb_denddia2",	 &sb_denddia2);
  setptr("sb_arbscale",	 &sb_arbscale);
  setptr("sb_arbscale2", &sb_arbscale2);
  setptr("sb_dia1",      &sb_dia1);
  setptr("sb_ca6_offm",  &sb_ca6_offm);
  setptr("sb_ca6_tauc",  &sb_ca6_tauc);
  setptr("sb_ca6_taud",  &sb_ca6_taud);
  setptr("disp_stim_incr", &disp_stim_incr);
  setptr("disp_stim_max", &disp_stim_max);
  setptr("disp_stim_min", &disp_stim_min);
  setptr("db1_lam",	 &db1_lam);
  setptr("db2_lam",	 &db2_lam);
  setptr("db1_lax",	 &db1_lax);
  setptr("db2_lax",	 &db2_lax);
  setptr("db_vs",	 &db_vs);

  setptr("na5offs",	 &na5offs);
  setptr("k1offm",	 &k1offm);
  setptr("k6offm",	 &k6offm);
  
  setptr("g_am_dsgc",	 &g_am_dsgc);
  setptr("g_am2_dsgc",	 &g_am2_dsgc);
  setptr("vhold",	 &vhold);
  setptr("dbp1_mrate",	 &dbp1_mrate);
  setptr("dbp1_mpool",	 &dbp1_mpool);
  setptr("dbp1_hdur",	 &dbp1_hdur);
  setptr("dbp1_hfilt",	 &dbp1_hfilt);
  setptr("dbp1_hgain",	 &dbp1_hgain);
  setptr("dbp1_hoffs",	 &dbp1_hoffs);
  setptr("dbp2_mrate",	 &dbp2_mrate);
  setptr("dbp2_mpool",	 &dbp2_mpool);
  setptr("g_sbac_dsgc",	 &g_sbac_dsgc);
  setptr("g_sbac_sbac",	 &g_sbac_sbac);
  setptr("n_sbac_sbac",	 &n_sbac_sbac);
  setptr("m_sbac_sbac",	 &m_sbac_sbac);
  setptr("r_sbac_sbac",	 &r_sbac_sbac);
  setptr("g_sbac_dbp1",	 &g_sbac_dbp1);
  setptr("g_dbp1_sbac",	 &g_dbp1_sbac);
  setptr("n_dbp1_sbac",	 &n_dbp1_sbac);
  setptr("v_dbp1_sbac",	 &v_dbp1_sbac);
  setptr("g_dbp1_dsgc",	 &g_dbp1_dsgc);
  setptr("g_dbp1_ams",	 &g_dbp1_ams);
  setptr("n_dbp1_ams",	 &n_dbp1_ams);
  setptr("g_ams_sbac",	 &g_ams_sbac);
  setptr("n_ams_sbac",	 &n_ams_sbac);
  setptr("ams_synanpi",	 &ams_synanpi);
  setptr("ams_synanpo",	 &ams_synanpo);
  setptr("g_dbp2_dsgc",	 &g_dbp2_dsgc);
  setptr("n_dbp2_sbac",	 &n_dbp2_sbac);
  setptr("v_dbp2_sbac",	 &v_dbp2_sbac);
  setptr("g_dbp2_sbac",	 &g_dbp2_sbac);
  setptr("idist",	 &idist);
  setptr("mglur2",	 &mglur2);

  setptr("sbac_dia",	 &sbac_dia);
  setptr("sbac_soma_z",	 &sbac_soma_z);
  setptr("sbac_soma_z2", &sbac_soma_z2);
  setptr("sbac_color",	 &sbac_color);
  setptr("sbac_dens",	 &sbac_dens);
  setptr("sbac_regu",	 &sbac_regu);
  setptr("sbsynang",	 &sbsynang);
  setptr("dsgc_color",	 &dsgc_color);
  setptr("draw_synapse", &draw_synapse);
  setptr("sbr",		 &sbr);
  setptr("dbp1_dsgc_reg",&dbp1_dsgc_reg);
  setptr("dbp1_dens",	 &dbp1_dens);
  setptr("dbp2_dens",	 &dbp2_dens);
  setptr("dbp1_maxsdist",&dbp1_maxsdist);
  setptr("dbp2_maxsdist",&dbp2_maxsdist);
  
  setptr("dsgc_denddia",&dsgc_denddia);

  setptr("blur_csize",	 &blur_csize);
  setptr("blur_ssize",	 &blur_ssize);
  setptr("surr_delay",	 &surr_delay);
  setptr("center_weight",&center_weight);
  setptr("surr_weight",	 &surr_weight);
  setptr("surr_tau",	 &surr_tau);
  setptr("stim_tstep",	 &stim_tstep);
  setptr("stim_fnum",	 &stim_fnum);

  setptr("sbac_synrng", &sbac_synrng);
  setptr("sbac_synanp", &sbac_synanp);
  setptr("sb_mglur_maxdist", &sb_mglur_maxdist);

  setptr("sreg1",   &sreg1);
  setptr("sreg2",   &sreg2);
  setptr("sreg3",   &sreg3);
  setptr("sreg4",   &sreg4);
  setptr("sreg5",   &sreg5);
  setptr("sreg6",   &sreg6);
  setptr("stonic",  &stonic);
  setptr("db1_r5",  &db1_r5);
  setptr("db2_r5",  &db2_r5);
  setptr("db1_r6",  &db1_r6);
  setptr("db2_r6",  &db2_r6);

  setptr("dbp1_anpi",   &dbp1_anpi);
  setptr("dbp1_anpo",   &dbp1_anpo);
  setptr("dbp2_anpi",   &dbp2_anpi);
  setptr("dbp2_anpo",   &dbp2_anpo);

  setptr("sbac_maxsdist", &sbac_maxsdist);
  setptr("sbac_synanni",  &sbac_synanni);
  setptr("sbac_isynanpi", &sbac_isynanpi);
  setptr("sbac_isynanpo", &sbac_isynanpo);
  setptr("sbac_isynanni", &sbac_isynanni);
  setptr("sbac_isynrngi", &sbac_isynrngi);
  setptr("sbac_ithr",     &sbac_ithr);
  setptr("sbac_synspac",  &sbac_synspac);
  setptr("sb_db_anni",    &sb_db_anni);
  setptr("sb_db_dur",     &sb_db_dur);
  setptr("sb_db_casens",  &sb_db_casens);
  setptr("sb_db_synr",    &sb_db_synr);

  setptr("naax",          &naax);
  setptr("nahd",          &nahd);
  setptr("naab",          &naab);
  setptr("kax",           &kax);
  setptr("ca5ax",         &ca5ax);
  setptr("ca7ax",         &ca7ax);
  setptr("ca6ax",         &ca6ax);
  setptr("ca1ax",         &ca1ax);
  setptr("db1_ca1ax",     &db1_ca1ax);
  setptr("db1_ca6ax",     &db1_ca6ax);
  setptr("db1_cap",       &db1_cap);
  setptr("db2_cap",       &db2_cap);
  setptr("db2_ca1ax",     &db2_ca1ax);
  setptr("db2_ca6ax",     &db2_ca6ax);
  setptr("db1_axdia",     &db1_axdia);
  setptr("db2_axdia",     &db2_axdia);

  setptr("wac_g",       &wac_g);

  setptr("barwidth",   &barwidth);
  setptr("barlength",  &barlength);
  setptr("tfreq",      &tfreq);
  setptr("stim_theta", &stim_theta);
  setptr("rstim_theta", &rstim_theta);
  setptr("minten",     &minten);
  setptr("scontrast",  &scontrast);
  setptr("scontrastn", &scontrastn);
  setptr("stimtime",   &stimtime);
  setptr("velocity",   &velocity);
  setptr("stimx",      &stimx);
  setptr("stimy",      &stimy);
  setptr("annrad",     &annrad);
  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("sbac_istim", &sbac_istim);
  setptr("stimdur",    &stimdur);
  setptr("prestimdur", &prestimdur);
  setptr("poststimdur",&poststimdur);
  setptr("tailcurdur", &tailcurdur);
  setptr("noise_dur",  &noise_dur);
  setptr("spotdur",    &spotdur);
  setptr("spotint",    &spotint);
  setptr("ncycles",    &ncycles);
  setptr("sblur",      &sblur);
  setptr("movein",     &movein);
  setptr("predur",     &predur);
  setptr("spotdia",    &spotdia);
  setptr("spotdia_incr",&spotdia_incr);
  setptr("spotdia_max",&spotdia_max);
  setptr("annulus_width",&annulus_width);
  setptr("annulus_odia",&annulus_odia);
  setptr("barwidth_max",&barwidth_max);
  setptr("barwidth_incr",&barwidth_incr);
  setptr("sdia",       &sdia);
  setptr("spdia",      &spdia);
  setptr("sndia",      &sndia);
  setptr("smdia",      &smdia);
  setptr("sddia",      &sddia);
  setptr("orad1",      &orad1);
  setptr("irad1",      &irad1);
  setptr("irad2",      &irad2);
  setptr("orad2",      &orad2);
  setptr("mask_dia",   &mask_dia);
  setptr("mask_x",     &mask_x);
  setptr("mask_y",     &mask_y);


  setptr("sbaclm",     &sbaclm);
  setptr("sb1mul",     &sb1mul);
  setptr("dsgclm",     &dsgclm);
  setptr("sb_rid",     &sb_rid);
  setptr("sb_rii",     &sb_rii);
  setptr("sb_rip",     &sb_rip);
  setptr("sb_rm",      &sb_rm);
  setptr("sb_rmp",     &sb_rmp);
  setptr("sb_rmd",     &sb_rmd);
  setptr("sb_vs",      &sb_vs);
  setptr("sb_vr",      &sb_vr);
  setptr("sb_cap",     &sb_cap);
  setptr("sb_capm",    &sb_capm);
  setptr("sb_capp",    &sb_capp);

  setptr("dendrm",     &dendrm);

  setptr("amna",	&amna);
  setptr("amk",		&amk);
  setptr("amca",	&amca);
  setptr("am2ca",	&am2ca);
  setptr("amsca",	&amsca);
  setptr("amvst",	&amvst);
  setptr("amvrev",	&amvrev);
  setptr("amrm",	&amrm);

  setptr("nasoma",     &nasoma);
  setptr("naprox",     &naprox);
  setptr("namid",      &namid);
  setptr("nadist",     &nadist);

  setptr("kdist",      &kdist);
  setptr("kmid",       &kmid);
  setptr("kprox",      &kprox);
  setptr("ksoma",      &ksoma);
  
  setptr("kdrs",       &kdrs);
  setptr("kdrp",       &kdrp);
  setptr("kdrm",       &kdrm);
  setptr("kdrd",       &kdrd);

  setptr("camid",      &camid);
  setptr("cadist",     &cadist);
  setptr("catsoma",    &catsoma);
  setptr("catprox",    &catprox);
  setptr("catmid",     &catmid);
  setptr("catdist",    &catdist);

  setptr("cadiap",    &cadiap);
  setptr("cadiam",    &cadiam);
  setptr("cadiad",    &cadiad);

  setptr("dsomadia",   &dsomadia);
  setptr("dtreedia",   &dtreedia);

  nvalfile       = "nval_dsgc_sbac.n";
  dsgc_densfile  = "dens_dsgc_sbac.n";
  dbp1_densfile  = "dens_dbp1_sbac.n" ;
  dbp2_densfile  = "dens_dbp2_sbac.n" ;
  chanparamsfile = "chanparams_dsgc_sbac";

  setvar();		// set variables from command line so sbac_species is correctly set
			//
  if (notinit(dsomadia)) dsomadia = 10;
  if (notinit(sbac_species)) sbac_species = 1;


  sreg1 = 1001;		// _SREG1, overall bp input to sbac,  sets of regions in dens_sbac_bar.n
  sreg2 = 1002;		// _SREG2, GABA output from sbac  
  sreg3 = 1003;		// _SREG3, inhibition to sbac
  sreg5 = 1005;		// _SREG5, transient bp input to sbac
  sreg6 = 1006;		// _SREG6, tonic bp input to sbac

  if (sbac_species==1) {		 // parameter values in nval file specific to mouse

      if (notinit(sbac_densfile))  sbac_densfile = "dens_sbac_bar.n";
      if (notinit(sbac_densfile2)) sbac_densfile2 = "dens_sbac_bar2.n";

      if (notinit(sbac_dens))     sbac_dens = 500;       // set sbac density, avg dist = 1/sqrt(dens)
      if (notinit(sbac_regu))     sbac_regu = 8;         // set sbac regularity (random locations)
      if (notinit(sbsynang))      sbsynang = 180;        // set dendritic orient for sbac->dsgc synapse
      if (notinit(sbac_isynanni)) sbac_isynanni = 0;     // inner radius of inhib annulus in presyn cell 
      if (notinit(sbac_isynrngi)) sbac_isynrngi = -1300; // range of angles for inhib synapses 
      if (notinit(sbac_isynanpi)) sbac_isynanpi = 0;     // inner radius of inhib annulus in postsyn cell 
      if (notinit(sbac_isynanpo)) sbac_isynanpo = 40;    // outer radius of inhib annulus in postsyn cell 
      if (notinit(sb_db_anni))    sb_db_anni = 100;      // inner radius of presyn annulus in sbac -> dbp1 
      if (notinit(sb_db_dur))     sb_db_dur = 2;         // time const of feedback sbac -> dbp1 
      if (notinit(sb_db_casens))  sb_db_casens = 1;      // Ca sens of feedback sbac -> dbp1 
      if (notinit(sb_db_synr))    sb_db_synr = 0;        // synrng: dendritic angle for sbac->dbp1 syn 
      if (notinit(sbac_ithr))     sbac_ithr = -0.05;     // sbac -> sbac synapse threshold in nval_sbac_stim.n
      if (notinit (sbspac)) {
	    // if (sbarr==100)       sbspac = 0.07 * sbac_dens; // sbac area*dens/mm2, gives coverage factor 
	    if (sbarr==100)       sbspac = 1/sqrt(sbac_dens) * 1e3; // nearest-neighbor dist 
	    else 		  sbspac = 100;
      };
      if (notinit (sbac_synspac)) sbac_synspac = 5;
      if (notinit (db1_morph))    db1_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db2_morph))    db2_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db1_ca1ax))    db1_ca1ax = 0;
      if (notinit (db1_ca6ax))    db1_ca6ax = 1e-3;
      if (notinit (db2_ca1ax))    db2_ca1ax = 0;
      if (notinit (db2_ca6ax))    db2_ca6ax = 1e-3;
      if (notinit (db1_cap))      db1_cap = 1e-6;
      if (notinit (db2_cap))      db2_cap = 4e-6;
      if (notinit (dbp1_dens))    dbp1_dens = 9500;
      if (notinit (dbp2_dens))    dbp2_dens = 9500;
      if (notinit (dbp1_maxsdist)) dbp1_maxsdist = 10;
      if (notinit (db1_biophys))  db1_biophys = 1;	// use density file for dbp1
      if (notinit (db2_biophys))  db2_biophys = 1;	// use density file for dbp2 
      if (notinit (db1_axdia))    db1_axdia = 0.6;      // dia mult.for dbp1 axon  
      if (notinit (db2_axdia))    db2_axdia = 0.6;      // dia mult.for dbp1 axon  
      if (notinit (db1_sca))      db1_sca = 1;		// use calcium for release for dbp1
      if (notinit (db2_sca))      db2_sca = 1;		// use calcium for release for dbp2 
      if (notinit (db1_r5))       db1_r5 = 0.5;		// SREG6 for dens_sbac_bar.n regs 5-6 (tonic BP )
      if (notinit (db2_r5))       db2_r5 = 0.5;		// SREG5 for dens_sbac_bar.n regs 5-6 (transient BP )
      if (notinit (db1_ca1_offm)) db1_ca1_offm = 0;	// dbp1 Ca1 offsetm 
      if (notinit (db1_ca6_offm)) db1_ca6_offm = 0.028;	// dbp1 Ca6 offsetm 
      if (notinit (db2_ca1_offm)) db2_ca1_offm = 0;	// dbp2 Ca1 offsetm 
      if (notinit (db2_ca6_offm)) db2_ca6_offm = 0.028;	// dbp2 Ca6 offsetm 
      if (notinit (db2_ca6_tauc)) db2_ca6_tauc = 0.5;	// dbp2 Ca6 tauc 
      if (notinit (db2_ca6_taud)) db2_ca6_taud = 0.5;	// dbp2 Ca6 taud 
      if (notinit (db1_lam))      db1_lam = 0.05;	// dbp1: make 2 comps
      if (notinit (db1_lax))      db1_lax = 0.05;	// dbp1: make 2 comps
      if (notinit (db2_lam))      db2_lam = 0.05;	// dbp2: make 2 comps
      if (notinit (db2_lax))      db2_lax = 0.05;	// dbp2: make 2 comps
      if (notinit(na5offs))       na5offs = 0;          // Na5 offm, offh for dbp2 chanparams_dsgc_sbac_db2  
      if (notinit(naax))          naax = 0;             //  density of Na chans in axon
      if (notinit(nahd))          nahd = 0e-3;          //  density of Na chans in high-density region
      if (notinit(naab))          naab = 200e-3;        //  density of Na chans in axon base
      if (notinit(k6offm))        k6offm = 0.005;       //  K6 (Kv3b) offm for dbp2 chanparams_dsgc_sbac  
      if (notinit(k1offm))        k1offm = 0.01;        //  K1 (Kdr) offm for dbp2 chanparams_dsgc_sbac

  } else if (sbac_species==2) {			// parameter values specific to rabbit 

      if (notinit(sbac_densfile))  sbac_densfile = "dens_sbac_bar.n";
      if (notinit(sbac_densfile2)) sbac_densfile2 = "dens_sbac_bar2.n";

      if (notinit(sbac_dens))     sbac_dens = 200;       // set sbac density, avg dist = 1/sqrt(dens)
      if (notinit(sbac_regu))     sbac_regu = 8;         // set sbac regularity (random locations)
      if (notinit(sbsynang))      sbsynang = 180;        // set dendritic orient for sbac->dsgc synapse
      if (notinit(sbac_isynanni)) sbac_isynanni = 100;   // inner radius of inhib annulus in presyn cell 
      if (notinit(sbac_isynrngi)) sbac_isynrngi = -1300; // range of angles for inhib synapses 
      if (notinit(sbac_isynanpi)) sbac_isynanpi = 0;     // inner radius of inhib annulus in postsyn cell 
      if (notinit(sbac_isynanpo)) sbac_isynanpo = 40;    // outer radius of inhib annulus in postsyn cell 
      if (notinit(sb_db_anni))    sb_db_anni = 100;      // inner radius of annulus in sbac -> dbp1 
      if (notinit(sb_db_dur))     sb_db_dur = 2;         // time const of feedback sbac -> dbp1 
      if (notinit(sb_db_casens))  sb_db_casens = 1;      // Ca sens of feedback sbac -> dbp1 
      if (notinit(sb_db_synr))    sb_db_synr = 0;        // synrng: dendritic angle for sbac->dbp1 syn 
      if (notinit(sbac_ithr))     sbac_ithr = -0.05;     // sbac -> sbac synapse threshold in nval_sbac_stim.n
      if (notinit (sbspac)) {
	    // if (sbarr==100)       sbspac = 0.07 * sbac_dens; // sbac area*dens/mm2, gives coverage factor 
	    if (sbarr==100)       sbspac = 1/sqrt(sbac_dens) * 1e3; // nearest-neighbor dist 
	    else 		  sbspac = 250;
      };
      if (notinit (sbac_synspac)) sbac_synspac = 5;
      if (notinit (db1_morph))    db1_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db2_morph))    db2_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db1_sca))      db1_sca = 1;		// use calcium for release for dbp1
      if (notinit (db2_sca))      db2_sca = 1;		// use calcium for release for dbp2 
      if (notinit (db1_ca1ax))    db1_ca1ax = 0;
      if (notinit (db1_ca6ax))    db1_ca6ax = 1e-3;
      if (notinit (db2_ca1ax))    db2_ca1ax = 0;
      if (notinit (db2_ca6ax))    db2_ca6ax = 1e-3;
      if (notinit (db1_cap))      db1_cap = 1e-6;
      if (notinit (db2_cap))      db2_cap = 1e-6;
      if (notinit (db1_ca1_offm)) db1_ca1_offm = 0;	// dbp1 Ca1 offsetm 
      if (notinit (db1_ca6_offm)) db1_ca6_offm = 0.028;	// dbp1 Ca6 offsetm 
      if (notinit (db2_ca1_offm)) db2_ca1_offm = 0;	// dbp1 Ca1 offsetm 
      if (notinit (db2_ca6_offm)) db2_ca6_offm = 0.028;	// dbp1 Ca6 offsetm 
      if (notinit (db2_ca6_tauc)) db2_ca6_tauc = 0.5;	// dbp2 Ca6 tauc 
      if (notinit (db2_ca6_taud)) db2_ca6_taud = 0.5;	// dbp2 Ca6 taud 
      if (notinit (dbp1_dens))    dbp1_dens = 7500;
      if (notinit (dbp2_dens))    dbp2_dens = 7500;
      if (notinit (dbp1_maxsdist)) dbp1_maxsdist = 10;
      if (notinit (db1_biophys))  db1_biophys = 1;	// use density file for dbp1
      if (notinit (db2_biophys))  db2_biophys = 1;	// use density file for dbp2
      if (notinit (db1_axdia))    db1_axdia = 0.6;      // dia mult.for dbp1 axon  
      if (notinit (db2_axdia))    db2_axdia = 0.6;      // dia mult.for dbp1 axon  
      if (notinit (db1_lam))      db1_lam = 0.05;	// dbp1: make 2 comps
      if (notinit (db1_lax))      db1_lax = 0.05;	// dbp1: make 2 comps
      if (notinit (db2_lam))      db2_lam = 0.05;	// dbp2: make 2 comps
      if (notinit (db2_lax))      db2_lax = 0.05;	// dbp2: make 2 comps
      if (notinit(na5offs))       na5offs = 0;          // Na5 offm, offh for dbp2 chanparams_dsgc_sbac_db2  
      if (notinit(naax))          naax = 0;             //  density of Na chans in axon
      if (notinit(nahd))          nahd = 0e-3;          //  density of Na chans in high-density region
      if (notinit(naab))          naab = 200e-3;        //  density of Na chans in axon base
      if (notinit(k6offm))        k6offm = 0.005;       //  K6 (Kv3b) offm for dbp2 chanparams_dsgc_sbac  
      if (notinit(k1offm))        k1offm = 0.01;        //  K1 (Kdr) offm for dbp2 chanparams_dsgc_sbac

  } else if (sbac_species==3) {		// parameter values specific to primate

      dbp1_densfile  = "dens_midg_primate.n" ;
      dbp2_densfile  = "dens_db4_primate.n" ;
      if (notinit(sbac_densfile))  sbac_densfile  = "dens_sbac_bar_db4.n";
      if (notinit(sbac_densfile2)) sbac_densfile2 = "dens_sbac_bar2_db4.n";
      chanparamsfile = "chanparams_dsgc_sbac_db4";

      make_dbp2_sbac = 1;

      if (notinit(sbac_dens))     sbac_dens = 500;       // set sbac density, avg dist = 1/sqrt(dens)
      if (notinit(sbac_regu))     sbac_regu = 8;         // set sbac regularity (random locations)
      if (notinit(sbsynang))      sbsynang = 180;        // set dendritic orient for sbac->dsgc synapse
      if (notinit(sbac_isynanni)) sbac_isynanni = 100;   // inner radius of inhib annulus in presyn cell 
      if (notinit(sbac_isynrngi)) sbac_isynrngi = -1300; // range of angles for inhib synapses 
      if (notinit(sbac_isynanpi)) sbac_isynanpi = 0;     // inner radius of inhib annulus in postsyn cell 
      if (notinit(sbac_isynanpo)) sbac_isynanpo = 40;    // outer radius of inhib annulus in postsyn cell 
      if (notinit(sb_db_anni))    sb_db_anni = 100;      // inner radius of annulus in sbac -> dbp1 
      if (notinit(sb_db_dur))     sb_db_dur = 2;         // time const of feedback sbac -> dbp1 
      if (notinit(sb_db_casens))  sb_db_casens = 1;      // Ca sens of feedback sbac -> dbp1 
      if (notinit(sb_db_synr))    sb_db_synr = 0;        // synrng: dendritic angle for sbac->dbp1 syn 
      if (notinit(sbac_ithr))     sbac_ithr = -0.05;     // sbac -> sbac synapse threshold in nval_sbac_stim.n
      if (notinit (sbspac)) {
	    // if (sbarr==100)       sbspac = 0.07 * sbac_dens; // sbac area*dens/mm2, gives coverage factor 
	    if (sbarr==100)       sbspac = 1/sqrt(sbac_dens) * 1e3; // nearest-neighbor dist 
	    else 		  sbspac = 160;
      };
      if (notinit (sbac_synspac)) sbac_synspac = 5;

      if (notinit (dbp1_dens))    dbp1_dens = 6000;	// set in dsgc_sbac_r
      if (notinit (db1_morph))    db1_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db1_biophys))  db1_biophys = 1;	// use density file for dbp1
      if (notinit (db1_axdia))    db1_axdia = 0.6;      // dia mult.for dbp1 axon  
      if (notinit (db1_sca))      db1_sca = 1;		// dbp1->sbac sensca
      if (notinit (db1_lam))      db1_lam = 0.05;	// dbp1: make 2 comps
      if (notinit (db1_lax))      db1_lax = 0.05;	// dbp1: make 2 comps
      if (notinit (db1_ca1ax))    db1_ca1ax = 0;
      if (notinit (db1_ca6ax))    db1_ca6ax = 1e-3;     //  density of Ca chans in bp axon 
      if (notinit (db1_ca1_offm)) db1_ca1_offm = 0;     // dbp1 Ca1 offsetm 
      if (notinit (db1_ca6_offm)) db1_ca6_offm = 0.028;	// dbp1 Ca6 offsetm 
      if (notinit (db1_ca7_offm)) db1_ca7_offm = 0.028;	// dbp1 Ca7 offsetm 
      if (notinit (dbp1_maxsdist)) dbp1_maxsdist = 10;
      if (notinit (dbp1_mpool))	  dbp1_mpool = 15;	// dbp1 smrrpool 
      if (notinit (dbp1_mrate))	  dbp1_mrate = 100;	// dbp1 smaxrate 
      
      if (notinit (dbp2_dens))    dbp2_dens = 6000;	// set in dsgc_sbac_r
      if (notinit (db2_morph))    db2_morph = 4;	// A4 morphology, no dendrites, 4 axon regions 
      if (notinit (db2_biophys))  db2_biophys = 1;	// use density file 
      if (notinit (db2_axdia))    db2_axdia = 0.6;      // dia mult.for dbp2 axon  
      if (notinit (db2_sca))      db2_sca = 1;		// dbp2->sbac sensca
      if (notinit (db2_ca1ax))    db2_ca1ax = 0;        //  density of Ca1 chans in bp axon with Na chans
      if (notinit (db2_ca6ax))    db2_ca6ax = 1e-3;     //  density of Ca6 chans in bp axon with Na chans
      if (notinit (db2_ca1_offm)) db2_ca1_offm = 0;	// dbp2 Ca1 offsetm 
      if (notinit (db2_ca6_offm)) db2_ca6_offm = 0.028;	// dbp2 Ca6 offsetm 
      if (notinit (db2_ca7_offm)) db2_ca7_offm = 0.028;	// dbp1 Ca7 offsetm 
      if (notinit (db2_ca6_tauc)) db2_ca6_tauc = 0.5;	// dbp2 Ca6 tauc 
      if (notinit (db2_ca6_taud)) db2_ca6_taud = 0.5;	// dbp2 Ca6 taud 
      if (notinit (db2_lam))      db2_lam = 0.05;	// dbp2: make 2 comps
      if (notinit (db2_lax))      db2_lax = 0.05;	// dbp2: make 2 comps
      if (notinit (dbp2_maxsdist)) dbp2_maxsdist = 10;
      if (notinit (dbp2_mpool))	  dbp2_mpool = 8;	// dbp2 smrrpool 
      if (notinit (dbp2_mrate))	  dbp2_mrate = 100;	// dbp2 smaxrate 

      if (notinit (db1_r6))       db1_r6 = 0.5;		// SREG6 for dens_sbac_bar_db4.n reg 6-9 (tonic BP )
      if (notinit (db2_r6))       db2_r6 = 0.5;		// SREG5 for dens_sbac_bar_db4.n reg 3 (transient BP )

      if (notinit(dbp1_hfilt))    dbp1_hfilt = 0;       //  1 => dbp1 high-pass filt
      if (notinit(dbp1_hdur))     dbp1_hdur  = 30;      //  time const of dbp1 high-pass filt (msec)
      if (notinit(dbp1_hgain))    dbp1_hgain = 1.0;     //  rel gain of dbp1 high-pass filt (0 - 1) 
      if (notinit(stonic))	  stonic = sreg6;	//  place dbp1s throughout

      if (notinit(na5offs))       na5offs = 0;          //  Na5 offm, offh for dbp2 chanparams_dsgc_sbac_db4  
      if (notinit(k1offm))        k1offm = 0.01;        //  K1 (Kdr) offm for dbp2 chanparams_dsgc_sbac_db4  
      if (notinit(k6offm))        k6offm = 0.005;       //  K6 (Kv3b) offm for dbp2 chanparams_dsgc_sbac_db4  
      if (notinit(naax))          naax = 0;             //  density of Na chans in axon
      if (notinit(nahd))          nahd = 0e-3;          //  density of Na chans in high-density region
      if (notinit(naab))          naab = 200e-3;        //  density of Na chans in axon base
      if (notinit(kax))           kax = 4e-3;           //  density of K chans in axon
      if (notinit(ca1ax))         ca1ax = 0e-3;         //  density of Ca chans in dbp1 (db5) bp axon (sustained)
      if (notinit(ca5ax))         ca5ax = 0e-3;         //  density of Ca chans in bp axon
      if (notinit(db1_cap))       db1_cap = 3e-6;       //  Ca pump in db1 terminal
      if (notinit(db2_cap))       db2_cap = 3e-6;       //  Ca pump in db2 terminal 
      // DEND  = R1;					//  defs for db4 bipolar in morph_DB4_120605_02.na1b
      // DENDP = R2;					//    but conflicts with dsgc defs
      if (notinit(disp_stim_max)) disp_stim_max = -0.045; /* max level for stimulus display */
      if (notinit(disp_stim_min)) disp_stim_min = -0.052; /* min level for stimulus display */
    
  };  /* if (sbac_species==3) */

  if (notinit(center_weight))	 center_weight = 1;
  if (notinit(surr_weight))	 surr_weight = -0.8;
  if (notinit(surr_delay))	 surr_delay = 0.004;
  if (notinit(surr_tau))	 surr_tau = 0.02;
  if (notinit (blur_csize))	 blur_csize = 30;
  if (notinit (blur_ssize))	 blur_ssize = 120;

  setvarstr("stim_fnum");
  if (notinit(stim_fnum))	 stim_fnum = "1";
  if (notinit(draw_synapse))	 draw_synapse = 1;
  if (notinit(dsgc_denddia))	 dsgc_denddia = 1.0;
  if (notinit(dbp2_maxsdist))	 dbp2_maxsdist = dbp1_maxsdist;
  if (notinit (db2_sca))         db2_sca = 0;		// default: dbp2->sbac, don't use sensca
  if (notinit (db1_morph))       db1_morph = 2;		// default: use simple artificial morphology 
  if (notinit (db2_morph))       db2_morph = 2;		// default: use simple artificial morphology 
  if (notinit (db2_biophys))     db2_biophys = 0;	// default: don't use density file
  
  if (notinit(cadiap))		cadiap = 1;
  if (notinit(cadiam))		cadiam = 1;
  if (notinit(cadiad))		cadiad = 0.5;
  dcashd = 0.1;					// Ca shell depth (thickness, um)

  	/* excitatory synapses, for nval file */

  if (notinit(dbp1_mpool))  dbp1_mpool = 100;   // smrrpool 
  if (notinit(dbp1_mrate))  dbp1_mrate = 25;    // smaxrate 
  if (notinit(dbp2_mpool))  dbp2_mpool = 100;   // smrrpool 
  if (notinit(dbp2_mrate))  dbp2_mrate = 25;    // smaxrate 

  if (notinit(dbp1_hfilt))  { dbp1_hfilt = 0; }     //  1 => dbp1 high-pass filt
  if (notinit(dbp1_hdur))   { dbp1_hdur  = 50; }    //  time const of dbp1 high-pass filt (msec)
  if (notinit(dbp1_hgain))  { dbp1_hgain = 0.7; }   //  rel gain of dbp1 high-pass filt (0 - 1) 
  if (notinit(dbp1_hoffs))  { dbp1_hoffs = 0; }     //  offset of dbp1 high-pass filt (V) 

  if (notinit(dbp1_anpi))  { dbp1_anpi = 0; }       // inner radius of excit annulus in postsyn cell 
  if (notinit(dbp1_anpo))  { dbp1_anpo = 0; }       // outer radius of excit annulus in postsyn cell
  if (notinit(dbp2_anpi))  { dbp2_anpi = 0; }       // inner radius of excit annulus in postsyn cell
  if (notinit(dbp2_anpo))  { dbp2_anpo = 0; }       // outer radius of excit annulus in postsyn cell
  if (notinit(dbpthr))     { dbpthr    = -0.050; }  // dbp1->sbac thresh (not needed with dbp1 sensca 1)
  if (notinit(dbp2thr))    { dbp2thr   = -0.050; }  // dbp2->sbac thresh (not needed with dbp2 sensca 1) 
  if (notinit(n_dbp1_sbac))   n_dbp1_sbac = 0;	    //  synaptic noise dbp1 -> sbac
  if (notinit(v_dbp1_sbac))   v_dbp1_sbac = 1;	    //  synaptic gain dbp1 -> sbac (from 4)
  if (notinit(n_dbp2_sbac))   n_dbp2_sbac = 0;	    //  synaptic noise dbp1 -> sbac
  if (notinit(v_dbp2_sbac))   v_dbp2_sbac = 1;	    //  synaptic gain dbp1 -> sbac (from 4)

  if (notinit(sbac_synrng))  { sbac_synrng = -50; }    // range of angles for excitatory synapses 
  if (notinit(sbac_synanni)) { sbac_synanni = 0; }     // inner radius of annulus in presyn cell
  if (notinit(sbac_synanp))  { sbac_synanp = 0; }      // inner radius of excit annulus in postsyn cell
  if (notinit(sb_denddia))   sb_denddia  = 1;          // thickness factor for sbac #1 dendrites
  if (notinit(sb_denddia2))  sb_denddia2 = 1;          // thickness factor for sbac #2 dendrites
  if (notinit(sb_arbscale))  sb_arbscale  = 1;         // scale factor (x,y) for sbac #1 dendrites
  if (notinit(sb_arbscale2)) sb_arbscale2 = 1;         // scale factor (x,y) for sbac #2 dendrites
  if (notinit(sb_ca6_offm))  sb_ca6_offm = 0.026;      // offset for Ca6 chans in sbac dendrites
  if (notinit(sb_ca6_tauc))  sb_ca6_tauc = 8;          // tauc (inactivation) for Ca6 chans in sbac dendrites
  if (notinit(sb_ca6_taud))  sb_ca6_taud = 1;          // taud (reactivation) for Ca6 chans in sbac dendrites
  if (notinit(sb_dia1))      sb_dia1 = 1.0;	       // iniitial dendr thickness in morph_11_17_18_cell_1

  if (notinit(sbac_maxsdist)) { sbac_maxsdist = 7; }  // sbac -> sbac synapse max dist in nval_dsgc_sbac.n
  if (notinit(wac_g))         wac_g = 0;              /* conductance of synapse onto sbac */

  if (notinit(ivplot)) ivplot   = 0;       /* ivplot */
  if (notinit(outward)) outward   = 0;     /* outward */

  if (notinit(iscal)) iscal   = 2e-10;      /* ivplot scale */
  if (notinit(gscal)) gscal   = 1.5e-8;    /* ivplot scale */

  if (notinit(mask_dia)) mask_dia = 100;    /* mask dia for stimtype 9, masked moving bar */
  if (notinit(mask_x)) mask_x = 0;         /* mask dia for stimtype 9, masked moving bar */
  if (notinit(mask_y)) mask_y = 0;         /* mask dia for stimtype 9, masked moving bar */


//  if using dbp2s, set stonic = sreg6 to limit tonic input to proximal in mouse (see dens_sbac_bar.n) 

  if (make_dbp2_sbac > 0) {
	 if (n_dbp2 != 0) make_dbp2 = 1;	// when dbp2 inputs are included,
         if (notinit(stonic)) stonic = sreg6;	//  move dbp1 inputs proxomal (in nval_dsgc_sbac.n) 
  } else if (notinit(stonic)) stonic = sreg1;

  make_sbac_sbac = 1;

  make_am_dbp1 = 0;			// no feedback to dbp1s
  make_am2_dbp1 = 0;
  make_ams_dbp1 = 0;
  make_am_dbp2 = 0;
  make_am2_dbp2 = 0;
  make_ams_dbp2 = 0;
  make_ams_dsgc = 0;

  sb_init();
}

/* - - - - - - - - - - - - - - - - - - - */

int make_am_cell (int ct, double theta, double xloc, double yloc, int i, double somadist)
{
   int x;
   double *xarr,*yarr,*tharr;
#define AMARR  50

  make_ct(ct);
  switch (ct) {
     case AM:  xarr = amxarr;  yarr = amyarr;  tharr = amtharr;  break;
     case AM2: xarr = am2xarr; yarr = am2yarr; tharr = am2tharr; break;
     case AM3: xarr = am3xarr; yarr = am3yarr; tharr = am3tharr; break;
     case AM4: xarr = am4xarr; yarr = am4yarr; tharr = am4tharr; break;
												    }
     if (xarr==NULL) {
          xarr = (double *)emalloc(AMARR*sizeof(double));
	  yarr = (double *)emalloc(AMARR*sizeof(double));
	  tharr = (double *)emalloc(AMARR*sizeof(double));
	  switch (ct) {
	     case AM:   amxarr = xarr;  amyarr = yarr;  amtharr = tharr; break;
	     case AM2: am2xarr = xarr; am2yarr = yarr; am2tharr = tharr; break;
	     case AM3: am3xarr = xarr; am3yarr = yarr; am3tharr = tharr; break;
	     case AM4: am4xarr = xarr; am4yarr = yarr; am4tharr = tharr; break;
	  }
     }
     if (xarr!=NULL && i<AMARR) {
         xarr[i]  = -cosdeg(theta) * somadist + xloc;
	 yarr[i]  =  sindeg(theta) * somadist + yloc;
	 tharr[i++] = theta;
     } else
         fprintf (stderr,"expt_cbp_inh: make_am_cell, can't allocate memory for %s\n",cname[ct]);
  return (i);
}

/* - - - - - - - - - - - - - - - - - - - */

int make_am_cell (int ct, double theta, int i, double somadist)

{
	  return make_am_cell (ct, theta, 0, 0, i, somadist);
}


/*--------------------------------------------------------*/

#include "gprim.h"

void syn_draw2 (synapse *spnt, int color, double vrev, double dscale, double dia,
                                double length, double foreshorten, int hide)

/* draw synapse within small oriented frame */
/*   draw_synzpse==0 -> don't draw */
/*   draw_synzpse==1 -> draw, default */
/*   draw_synzpse==3 -> draw line but not circle or text*/

{
    int fill=1;
    double tlen;
    static char tbuf[15];

  if (!draw_synapse) return;
  dia *= dscale;                        /* draw circle with line */
  if (dia < 0) dia = -dia;
  color = -1;
  if (color < 0) {
      // if (vrev < -0.04) color = RED;
      // else              color = CYAN;
      if      (spnt->node1a==dbp1 && spnt->node2a==sbac) color = green;
      else if (spnt->node1a==dbp2 && spnt->node2a==sbac) color = blue;
      else if (spnt->node1a==dbp1 && spnt->node2a==dsgc) color = cyan;
      else if (spnt->node1a==sbac && spnt->node2a==sbac) color = blue;
      // else if (spnt->node1a==sbac && spnt->node2a==sbac) color = brown;
      else if (spnt->node1a==sbac && spnt->node2a==dbp1) color = red;
      else if (spnt->node1a==sbac && spnt->node2a==dsgc) color = ltblue;
  }
  if (draw_synapse==2) color = NOCOLOR;
  if (draw_synapse==3) color = GRAY;
  gpen (color);
  if (length > 1e-3) {
      if (draw_synapse!=3) {
         gmove (length/2.0,0.0);
         if (dia > 0.001) gcirc (dia/2.0,fill);
         else             gcirc (0.001,fill);
      }
      gmove (0,0);
      gdraw (length,0);
  }
  else if (draw_synapse!=3) gcirc (0.001,fill);
  // if      (spnt->node2a==sbac) gpen (green);
  // else if (spnt->node1a==dbp1) gpen (cyan);
  // else                         gpen (yellow);
  if (draw_synapse!=3) { 
     sprintf (tbuf,"%d>%d",spnt->node1b,spnt->node2b);     /* print pre and postsynaptic cell number */
     tlen = strlen(tbuf);
     gmove (length/2.0 -tlen*0.5, -2.0);
     gcwidth (1.5*dscale);
     gtext (tbuf);
  }
}

/*--------------------------------------------------------*/

void makcyl(double x1, double y1, double z1, double x2, double y2, double z2, double dia, int color);
const char *getcol (int color);

void raysyn_draw2 (synapse *spnt, double dia, double tx1, double ty1, double tz1,
                                             double tx2, double ty2, double tz2, int color)
/* draw synapse in 3d image */

{
    int fill=1;
    double cx, cy, cz, dx, dy, dz, tx, ty, tz;
    double tlen, length, theta;
    static char tbuf[15];
        
  if (!draw_synapse) return;
  if (dia < 0) dia = -dia;
  color = -1;
  if (color < 0) {
      // if (vrev < -0.04) color = RED;
      // else              color = CYAN;
      if      (spnt->node1a==dbp1 && spnt->node2a==sbac) color = green;
      else if (spnt->node1a==dbp2 && spnt->node2a==sbac) color = blue;
      else if (spnt->node1a==dbp1 && spnt->node2a==dsgc) color = cyan;
      else if (spnt->node1a==sbac && spnt->node2a==sbac) color = blue;
      //else if (spnt->node1a==sbac && spnt->node2a==sbac) color = brown;
      else if (spnt->node1a==sbac && spnt->node2a==dbp1) color = red;
      else if (spnt->node1a==sbac && spnt->node2a==dsgc) color = ltblue;
  }
  if (draw_synapse==2) color = NOCOLOR;

  cx = (tx1+tx2)*0.5;
  cy = (ty1+ty2)*0.5;
  cz = (tz1+tz2)*0.5;
  if (tx1!=tx2 || ty1!=ty2 || tz1!=tz2) {
     //makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia,color);
     makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia*2,color);
       
    ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
                  cx,cy,cz,dia/2.0,getcol(color));
  
     dx = tx2-tx1;
     dy = ty2-ty1;
     dz = tz2-tz1;
     length = sqrt(dx*dx + dy*dy);
     if (length > 1e-6) {
	     if (dx < 0) {
		theta = -asin (dy/length)*DEG; 
		tx = tx2; ty = ty2; tz = tz2;
	     } else      {
	      	 theta =  asin (dy/length)*DEG;
		 tx = tx1; ty = ty1; tz = tz1;
            }
     }  else theta = 0;
     sprintf (tbuf,"%d>%d",spnt->node1b,spnt->node2b);		// pre and postsyn cell number
     //rprintf (tbuf, tx, ty, tz, 3, theta, color);
  }


}

/*--------------------------------------------------------*/

void setparams(void)

  /*  set up default configuration for sb expts */
  /* cones, cone bipolars, sb, dsgc */

{
   int i, n, ct;
   double zmax, zmin;
   double t, theta_incr;
   double amx=0, amy=0;
   double somadist;
   double *sbrarr; 		// sbac rnd rotation array

  make_rods  = 0;
  make_cones = 0;
  make_dbp1  = 1;
  // make_dbp2  = 0;
  make_hbp1  = 0;
  make_hbp2  = 0;
  // make_ams   = 0;		// allow ams to be included from command line (default=0)
  make_sbac  = 1;
  make_dsgc  = 1;

  if (sb_biophys) setn(sbac,NCOLOR,RCOLOR);		/* set sbac color by region from dens_sbac_bar.n */
  // setn(sbac,NCOLOR,CNCOLOR);				/* set sbac color by region from dens_sbac_bar.n */
  //if (notinit(sbac_color)) sbac_color = CNCOLOR;	/* set sbac color by cell number */
  //if (notinit(sbac_color)) sbac_color = 8;
  if (!notinit(sbac_color)) setn(sbac,NCOLOR,sbac_color);   /* set sbac display color */
  if (!notinit(sbac_soma_z)) setn(sbac,SOMAZ,sbac_soma_z);  /* set sbac soma z location */
  if (!notinit(sbac_soma_z2)) setn(sbac,SOMAZ2,sbac_soma_z2);  /* set sbac soma z location */

  if (notinit(dbp1_dsgc_reg)) dbp1_dsgc_reg = R2;
  setsv(dbp1,SYNREGP,2,dbp1_dsgc_reg);

  if (notinit(dsgc_color)) dsgc_color = RCOLOR;		/* set dsgc color by cell region */
  setn(dsgc,NCOLOR,dsgc_color);  /* set dsgc display color */
  if (notinit(dsgclm)) dsgclm = 0.03;
  if (db2_morph>2) setn(dbp2,NCOLOR,RCOLOR);  /* set dbp2 display color */

  if (notinit(disp_stim_incr))  disp_stim_incr = 0.01;		/* time step for stimulus display */
  if (notinit(disp_stim_max))   disp_stim_max = -0.03;		/* max level for stimulus display */
  if (notinit(disp_stim_min))   disp_stim_min = -0.06;		/* min level for stimulus display */

  if (notinit(dendrm))         dendrm = drm;

  set_synapse_dr (syn_draw2);
  set_raysyn_dr  (raysyn_draw2);

#define GCARR 5
#define SBARR 100

   if (notinit(use_stimfile)) use_stimfile = 0;

   // if (notinit (sbtheta)) sbtheta = 0;
   if (notinit(spdia))        spdia = 0.54;	/* proximal dend dia factor in dens_sbac_bar.n */
   if (notinit(sndia))        sndia = spdia;	/* nearest to prox dend dia factor in dens_sbac_bar.n */
   if (notinit(smdia))        smdia = spdia;	/* medial dend dia factor in dens_sbac_bar.n */
   if (notinit(sdia))          sdia = 0.28;	/* dendrite dia factor in dens_sbac_bar.n */ 
   if (notinit(sddia))        sddia = 0.65;	/* distal dend dia factor in dens_sbac_bar.n */
   if (notinit(sbaclm))      sbaclm = 0.05;
   if (notinit(sb1mul))      sb1mul = 1.0;	/* sbac1 cmul factor in dens_sbac_bar2.n (inhib input mul) */
   if (notinit(sbac_r))      sbac_r = 0;	/* default sbac1 rotation */
   if (notinit(sbac_rnd))  sbac_rnd = 0;	/* default sbac1 random rotation multiplier for specified arrays */
   if (notinit(dsgc_x))      dsgc_x = -200;	/* default dsgc x location */
   if (notinit(dsgc_y))      dsgc_y = 50;	/* default dsgc y location */
   if (notinit(rloc))          rloc = 0.75;	/* default relative loc for bar starting position */

   if (notinit(sb_mglur_maxdist)) sb_mglur_maxdist = 10; /* default sbac mglur -> ca chan dist */

   if (notinit(sb_rid)) sb_rid = dri;
   if (notinit(sb_rii)) sb_rii = dri;
   if (notinit(sb_rip)) sb_rip = dri;
   if (notinit(sb_rm))   sb_rm = drm;
   if (notinit(sb_rmp)) sb_rmp = sb_rm;
   if (notinit(sb_rmd)) sb_rmd = sb_rmp;
   if (notinit(sb_vr))   sb_vr = -0.07;	        /* sbac vrev,   used in dens_sbac_bar.n */
   if (notinit(sb_vs))   sb_vs = -0.063;        /* sbac vstart, used in dens_sbac_bar.n */
   if (notinit(sb_cap))   sb_cap = 1e-6;        /* sbac ca pump, used in dens_sbac_bar.n */
   if (notinit(sb_capm)) sb_capm = 4e-7;        /* sbac ca pump, medial dendrites, used in dens_sbac_bar.n */
   if (notinit(sb_capp)) sb_capp = 4e-7;        /* sbac ca pump, prox dendrites, used in dens_sbac_bar.n */
   dcapkm = 30e-6;

   if (notinit(nasoma))   nasoma = 0e-3;	/* soma Na density, dens_sbac_bar.n */
   if (notinit(naprox))   naprox = 0e-3;	/* prox Na density, dens_sbac_bar.n */
   if (notinit(namid))     namid = 0e-3;	/* mid Na density, dens_sbac_bar.n */
   if (notinit(nadist))   nadist = 0e-3;	/* distal Na density, dens_sbac_bar.n */
   if (notinit(ksoma))     ksoma = 3e-3;	/* soma Kv3 density, dens_sbac_bar.n */
   if (notinit(kprox))     kprox = 2e-3;	/* prox Kv3 density, dens_sbac_bar.n */
   if (notinit(kmid))       kmid = 0e-3;	/* mid  Kv3 density, dens_sbac_bar.n */
   if (notinit(kdist))     kdist = 0e-3;	/* dist Kv3 density, dens_sbac_bar.n */
   if (notinit(kdrs))       kdrs = 0e-3;	/* soma Kdr density, dens_sbac_bar.n */
   if (notinit(kdrp))       kdrp = 0e-3;	/* prox Kdr density, dens_sbac_bar.n */
   if (notinit(kdrm))       kdrm = 0e-3;	/* mid  Kdr density, dens_sbac_bar.n */
   if (notinit(kdrd))       kdrd = 0e-3;	/* dist Kdr density, dens_sbac_bar.n */
   if (notinit(camid))     camid = 0e-3;	/* mid  Ca density, dens_sbac_bar.n */
   if (notinit(cadist))   cadist = 0e-3;	/* dist Ca density, dens_sbac_bar.n */
   if (notinit(catsoma)) catsoma = 0.2e-3;	/* soma Ca-T density, dens_sbac_bar.n */
   if (notinit(catprox)) catprox = 0.2e-3;	/* prox Ca-T density, dens_sbac_bar.n */
   if (notinit(catmid))   catmid = 3e-3;	/* mid  Ca-T density, dens_sbac_bar.n */
   if (notinit(catdist)) catdist = 3e-3;	/* dist Ca-T density, dens_sbac_bar.n */

   if (notinit(amna))      amna = 0e-3;		/* Na chan density, dens_dsgc_sbac_am.n */
   if (notinit(amk))        amk = 0e-3;		/* K  chan density, dens_dsgc_sbac_am.n */
   if (notinit(amca))      amca = 1e-3;		/* Ca chan density, dens_dsgc_sbac_am.n */
   if (notinit(am2ca))    am2ca = 0.5e-3;	/* Ca chan density, dens_dsgc_sbac_am2.n */
   if (notinit(amsca))    amsca = 0e-3;		/* Ca chan density, dens_dsgc_sbac_am.n */
   if (notinit(amvst))    amvst = -0.07;	/* Vstart, dens_dsgc_sbac_am.n */
   if (notinit(amvrev))  amvrev = -0.07;	/* Vrev, dens_dsgc_sbac_am.n */
   if (notinit(amrm))      amrm = 20e3;		/* Rm, dens_dsgc_sbac_am.n */

   make_am_dbp1 = 0;
   make_am2_dbp1 = 0;

   if (!notinit(g_sbac_dbp1)) {                 // make positive feedback to dbp1
           make_sbac_dbp1 = 1;
           setsv(sbac,SCOND,1,g_sbac_dbp1);
           if (notinit(db1_lam)) db1_lam = 0.1; // complam for dbp1 (needed for sbac -> dbp1 feedback)
   }
   if (notinit(db1_lam)) db1_lam = 0.5;
   setn(dbp1,COMPLAM,db1_lam);
   if (notinit(db1_lax)) db1_lax = 0.5;
   setn(dbp1,COMPLAM,db1_lax);

   if (notinit(db2_lam)) db2_lam = 0.5;
   setn(dbp2,COMPLAM,db2_lam);
   if (notinit(db2_lax)) db2_lax = 0.5;
   setn(dbp2,COMPLAM,db2_lax);

   if (notinit(mglur2))      mglur2 = 0;		/* if > 0 -> mglur2 modulation of Ca chans */

	/* definitions of regions for dens_dsgc_sbac.n dens_ file */

    DENDD     = DEND_DIST = R_1;
    DEND      		  = R_2;
    DENDP     = DEND_PROX = R_3;
    SOMA      	          = R_4;
    HCK       = HILLOCK   = R_5;
    AXONT     = AXON_THIN = R_6;
    AXON      	          = R_7;
    AXONP     = AXON_PROX = R_7;
    AXOND     = AXON_DIST = R_8;
    VARIC     = VARICOS   = R_8;

   if (n_dsgc > 0) {					// array of dsgc locations
     gcxarr = (double *)emalloc(GCARR*sizeof(double));
     gcyarr = (double *)emalloc(GCARR*sizeof(double));
     gctharr = (double *)emalloc(GCARR*sizeof(double));
     // gcxarr[0] = -150;
     // gcxarr[0] = -100;
     gcxarr[0] = dsgc_x;
     gcyarr[0] = dsgc_y;
     for (i=0; i<GCARR; i++) gctharr[i] = 0;		/* set up default cell numbering */
   }
   if (strcmp(sbac_file,"morph_sbac3b")==0 || strcmp(sbac_file,"morph_sbac3c")==0) {
        if (notinit(sbac_dia)) sbac_dia = 315;
   }
   else if (strcmp(sbac_file,"morph_sb1")==0) {
        if (notinit(sbac_dia)) sbac_dia = 315;
   }
   else if (strncmp(sbac_file,"morph_sbac_168",14)==0) {
        if (notinit(sbac_dia)) sbac_dia = 315;
        if (strncmp(sbac_file,"morph_sbac_168c5",14)==0) {  // 168c4/5 has flat arborization
           if (notinit(sbac_thetax)) sbac_thetax = 0;
           if (notinit(sbac_thetay)) sbac_thetay = 0;
	} else {
           if (notinit(sbac_thetax)) sbac_thetax = 2;
           if (notinit(sbac_thetay)) sbac_thetay = -2;
	};
   }

 if (!notinit(sbarr)) {		// array of sbac locations
   sbxarr = (double *)emalloc(SBARR*sizeof(double));
   sbyarr = (double *)emalloc(SBARR*sizeof(double));
   sbtharr = (double *)emalloc(SBARR*sizeof(double));
   sbnarr = (int *)emalloc(SBARR*sizeof(int));
   sbrarr = (double *)emalloc(SBARR*sizeof(double));

   for (i=0; i<SBARR; i++) {
	   sbnarr[i] = i+1;		/* set up default cell numbering */
	   sbrarr[i] = rrange(0,359);	/* set up random cell rotation */
   }

   if (sbarr==0) {
       sbxarr[0] =  0; sbyarr[0] = 0;	sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;  /* only 1 SBAC */
       n_sbac = 1;
   }
   if (sbarr== -2) {			/* second sbac for subtraction of capacitive transient */
       sbxarr[0] = 0; sbyarr[0] = 0; sbtharr[0] = 0;
       sbxarr[1] = 1; sbyarr[1] = 0; sbtharr[1] = 0;
       n_sbac = 2;
       make_sbac_sbac = 0;		/* don't connect these sbacs */
   }
   if (sbarr==1) {			/* aligned dendrites */
       sbxarr[0] = 0; sbyarr[0] = 0; sbtharr[0] = 0;
       sbxarr[1] = 1; sbyarr[1] = 0; sbtharr[1] = 0;
       n_sbac = 2;
   }
   if (sbarr==2) {			/* opposing dendrites */
       sbxarr[0] =  0;      sbyarr[0] = 0; sbtharr[0] = sbac_r;
       sbxarr[1] = -sbspac; sbyarr[1] = 0; sbtharr[1] = sbac_r;
       n_sbac = 2;
   }
   if (sbarr==3) {                    	/* 3 aligned */
       sbxarr[0] =  0;       sbyarr[0] = 0; sbtharr[0] = sbac_r; sbnarr[0]=1;
       sbxarr[1] =  -sbspac; sbyarr[1] = 0; sbtharr[1] = sbac_r; sbnarr[1]=2;
       sbxarr[2] =   sbspac; sbyarr[2] = 0; sbtharr[2] = sbac_r; sbnarr[2]=3;
       n_sbac = 3;
   }
   if (sbarr==4) {			/* 3 aligned, 1 opposing */
       sbxarr[0] =  sbspac; sbyarr[0] = 0; sbtharr[0] = 0;
       sbxarr[1] =  sbspac; sbyarr[1] = 0; sbtharr[1] = 0;
       sbxarr[2] =       0; sbyarr[2] = 0; sbtharr[2] = 0;
       sbxarr[3] = -sbspac; sbyarr[3] = 0; sbtharr[3] = 0;
       n_sbac = 4;
   }
   if (sbarr==7) {			/* 7 overlapping */
       sbxarr[0] =  3*sbspac; sbyarr[0] = 0; sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  2*sbspac; sbyarr[1] = 0; sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =    sbspac; sbyarr[2] = 0; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] =         0; sbyarr[3] = 0; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] =   -sbspac; sbyarr[4] = 0; sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -2*sbspac; sbyarr[5] = 0; sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] = -3*sbspac; sbyarr[6] = 0; sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;
        n_sbac = 7;
   }
   if (sbarr==9) {			/* 9 overlapping */
       sbxarr[0] =  4*sbspac; sbyarr[0] = 0; sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  3*sbspac; sbyarr[1] = 0; sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  2*sbspac; sbyarr[2] = 0; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] =    sbspac; sbyarr[3] = 0; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] =         0; sbyarr[4] = 0; sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] =   -sbspac; sbyarr[5] = 0; sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] = -2*sbspac; sbyarr[6] = 0; sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;
       sbxarr[7] = -3*sbspac; sbyarr[7] = 0; sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] = -4*sbspac; sbyarr[8] = 0; sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       n_sbac = 9;
   }
   if (sbarr==11) {			/* 11 overlapping */
       sbxarr[0] =  4*sbspac; sbyarr[0] = 0; sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  3*sbspac; sbyarr[1] = 0; sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  2*sbspac; sbyarr[2] = 0; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] =    sbspac; sbyarr[3] = 0; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] =         0; sbyarr[4] = 0; sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] =   -sbspac; sbyarr[5] = 0; sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] = -2*sbspac; sbyarr[6] = 0; sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;
       sbxarr[7] = -3*sbspac; sbyarr[7] = 0; sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] = -4*sbspac; sbyarr[8] = 0; sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] = -5*sbspac; sbyarr[9] = 0; sbtharr[9] = sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10]= -6*sbspac; sbyarr[10]= 0; sbtharr[10]= sbac_r + sbrarr[10]* sbac_rnd;
       n_sbac = 11;
   }
   if (sbarr==100) {			/* random placement, sbxarrsiz, sbyzrrsiz */
        if (notinit(sbarrsiz) && notinit(sbxarrsiz)) sbarrsiz = 100;
   }
   if (sbarr==102) {
       sbxarr[0] = 145; sbyarr[0] =  50; sbtharr[0] = 0;
       sbxarr[1] = 145; sbyarr[1] = -50; sbtharr[1] = 0;
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
       // r = 60;
       // r = 90;
       // sbac_r = 240;
       sbxarr[0] =  0;             sbyarr[0] = 0;             sbtharr[0] = 0   + sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;        sbyarr[1] = 0;             sbtharr[1] = 180 + sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;    sbyarr[2] = -0.866*sbspac; sbtharr[2] = 60  + sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;    sbyarr[3] = -0.866*sbspac; sbtharr[3] = 0   + sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;        sbyarr[4] = 0;             sbtharr[4] = 90  + sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;    sbyarr[5] = 0.866*sbspac;  sbtharr[5] = 240 + sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;    sbyarr[6] = 0.866*sbspac;  sbtharr[6] = 300 + sbac_r + sbrarr[6] * sbac_rnd;
       n_sbac = 7;
       sbr = 2;
   }
   if (sbarr==112) {			// 13 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;             sbyarr[0] = 0;              sbtharr[0] = 0   + sbac_r;
       // sbxarr[1] =  0;             sbyarr[1] = 0.2887*sbspac; sbtharr[1] = 90  + sbac_r;
       sbxarr[1] =  0;             sbyarr[1] = sbspac; sbtharr[1] = 90  + sbac_r;

       n_sbac = 2;
       sbr = 2;
   }
   if (sbarr==113) {			// 13 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;             sbyarr[0] = 0;             sbtharr[0] = 0   + sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;        sbyarr[1] = 0;             sbtharr[1] = 180 + sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;    sbyarr[2] = -0.866*sbspac; sbtharr[2] = 60  + sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;    sbyarr[3] = -0.866*sbspac; sbtharr[3] = 0   + sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;        sbyarr[4] = 0;             sbtharr[4] = 90  + sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;    sbyarr[5] = 0.866*sbspac;  sbtharr[5] = 60  + sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;    sbyarr[6] = 0.866*sbspac;  sbtharr[6] = 0   + sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;    sbyarr[7] = 0.2887*sbspac;  sbtharr[7] = 0   + sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;             sbyarr[8] = 0.5773*sbspac;  sbtharr[8] = 180 + sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;   sbyarr[9] = 0.2887*sbspac;  sbtharr[9] = 0   + sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;   sbyarr[10] = -0.2887*sbspac; sbtharr[10] = 0  + sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;            sbyarr[11] = -0.5773*sbspac; sbtharr[11] = 90  + sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;   sbyarr[12] = -0.2887*sbspac; sbtharr[12] = 60  + sbac_r + sbrarr[12] * sbac_rnd;

       n_sbac = 13;
       sbr = 2;
   }
   if (sbarr==117) {			// 17 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = 0   + sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = 180 + sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = 60  + sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = 0   + sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = 90  + sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = 60  + sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = 0   + sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = 0   + sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = 180 + sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = 0   + sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = 0  + sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = 90 + sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = 60 + sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = 180 + sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = 60  + sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = 90  + sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = 0   + sbac_r + sbrarr[16] * sbac_rnd;

       n_sbac = 17;
       sbr = 2;
   }
   if (sbarr==121) {			// 21 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = 0   + sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = 180 + sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = 60  + sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = 0   + sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = 90  + sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = 60  + sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = 0   + sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = 0   + sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = 180 + sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = 0   + sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = 0  + sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = 90 + sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = 60 + sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = 180 + sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = 60  + sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = 90  + sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = 0   + sbac_r + sbrarr[16] * sbac_rnd;

       sbxarr[17] =  1.5*sbspac;  sbyarr[17] =  0.2887*sbspac; sbtharr[17] = 30   + sbac_r + sbrarr[17] * sbac_rnd;
       sbxarr[18] =  1.5*sbspac;  sbyarr[18] = -0.2887*sbspac; sbtharr[18] = 60   + sbac_r + sbrarr[18] * sbac_rnd;
       sbxarr[19] = -1.5*sbspac;  sbyarr[19] =  0.2887*sbspac; sbtharr[19] = 90   + sbac_r + sbrarr[19] * sbac_rnd;
       sbxarr[20] = -1.5*sbspac;  sbyarr[20] = -0.2887*sbspac; sbtharr[20] = 120  + sbac_r + sbrarr[20] * sbac_rnd;

       n_sbac = 21;
       sbr = 2;
   }
   if (sbarr==123) {			// 23 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = sbac_r + sbrarr[16] * sbac_rnd;

       sbxarr[17] =  1.5*sbspac;  sbyarr[17] =  0.2887*sbspac; sbtharr[17] = sbac_r + sbrarr[17] * sbac_rnd;
       sbxarr[18] =  1.5*sbspac;  sbyarr[18] = -0.2887*sbspac; sbtharr[18] = sbac_r + sbrarr[18] * sbac_rnd;
       sbxarr[19] = -1.5*sbspac;  sbyarr[19] =  0.2887*sbspac; sbtharr[19] = sbac_r + sbrarr[19] * sbac_rnd;
       sbxarr[20] = -1.5*sbspac;  sbyarr[20] = -0.2887*sbspac; sbtharr[20] = sbac_r + sbrarr[20] * sbac_rnd;

       sbxarr[21] =  2*sbspac;    sbyarr[21] = 0;              sbtharr[21] = sbac_r + sbrarr[21] * sbac_rnd;
       sbxarr[22] = -2*sbspac;    sbyarr[22] = 0;              sbtharr[22] = sbac_r + sbrarr[22] * sbac_rnd;

       n_sbac = 23;
       sbr = 2;
   }
   if (sbarr==137) {			// 37 sbacs
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = sbac_r + sbrarr[16] * sbac_rnd;

       sbxarr[17] =  1.5*sbspac;  sbyarr[17] =  0.2887*sbspac; sbtharr[17] = sbac_r + sbrarr[17] * sbac_rnd;
       sbxarr[18] =  1.5*sbspac;  sbyarr[18] = -0.2887*sbspac; sbtharr[18] = sbac_r + sbrarr[18] * sbac_rnd;
       sbxarr[19] = -1.5*sbspac;  sbyarr[19] =  0.2887*sbspac; sbtharr[19] = sbac_r + sbrarr[19] * sbac_rnd;
       sbxarr[20] = -1.5*sbspac;  sbyarr[20] = -0.2887*sbspac; sbtharr[20] = sbac_r + sbrarr[20] * sbac_rnd;

       sbxarr[21] =  1.5*sbspac;  sbyarr[21] =  0.866*sbspac; sbtharr[21] = sbac_r + sbrarr[21] * sbac_rnd;
       sbxarr[22] =  1.5*sbspac;  sbyarr[22] = -0.866*sbspac; sbtharr[22] = sbac_r + sbrarr[22] * sbac_rnd;
       sbxarr[23] = -1.5*sbspac;  sbyarr[23] = -0.866*sbspac; sbtharr[23] = sbac_r + sbrarr[23] * sbac_rnd;
       sbxarr[24] = -1.5*sbspac;  sbyarr[24] =  0.866*sbspac; sbtharr[24] = sbac_r + sbrarr[24] * sbac_rnd;

       					// 2 / sqrt(3) -> 2 * 0.5773 -> 1.1547
					
       sbxarr[25] = 0;            sbyarr[25] =  1.1547*sbspac; sbtharr[25] = sbac_r + sbrarr[25] * sbac_rnd;
       sbxarr[26] = 0;            sbyarr[26] = -1.1547*sbspac; sbtharr[26] = sbac_r + sbrarr[26] * sbac_rnd;

       sbxarr[27] =  2 * sbspac;  sbyarr[27] = 0;              sbtharr[27] = sbac_r + sbrarr[27] * sbac_rnd;
       sbxarr[28] =  2 * sbspac;  sbyarr[28] = 0.5773*sbspac;  sbtharr[28] = sbac_r + sbrarr[28] * sbac_rnd;
       sbxarr[29] =  2 * sbspac;  sbyarr[29] = -0.5773*sbspac; sbtharr[29] = sbac_r + sbrarr[29] * sbac_rnd;
       
       sbxarr[30] = -2 * sbspac;  sbyarr[30] = 0;              sbtharr[30] = sbac_r + sbrarr[30] * sbac_rnd;
       sbxarr[31] = -2 * sbspac;  sbyarr[31] = 0.5773*sbspac;  sbtharr[31] = sbac_r + sbrarr[31] * sbac_rnd;
       sbxarr[32] = -2 * sbspac;  sbyarr[32] = -0.5773*sbspac; sbtharr[32] = sbac_r + sbrarr[32] * sbac_rnd;

       sbxarr[33] =  sbspac;      sbyarr[33] =  1.1547*sbspac; sbtharr[33] = sbac_r + sbrarr[33] * sbac_rnd;
       sbxarr[34] = -sbspac;      sbyarr[34] =  1.1547*sbspac; sbtharr[34] = sbac_r + sbrarr[34] * sbac_rnd;
       sbxarr[35] =  sbspac;      sbyarr[35] = -1.1547*sbspac; sbtharr[35] = sbac_r + sbrarr[35] * sbac_rnd;
       sbxarr[36] = -sbspac;      sbyarr[36] = -1.1547*sbspac; sbtharr[36] = sbac_r + sbrarr[36] * sbac_rnd;

       n_sbac = 37;
       sbr = 2;
   }
   if (sbarr==145) {			// 45 sbacs, all same rotation
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = sbac_r + sbrarr[16] * sbac_rnd;

       sbxarr[17] =  1.5*sbspac;  sbyarr[17] =  0.2887*sbspac; sbtharr[17] = sbac_r + sbrarr[17] * sbac_rnd;
       sbxarr[18] =  1.5*sbspac;  sbyarr[18] = -0.2887*sbspac; sbtharr[18] = sbac_r + sbrarr[18] * sbac_rnd;
       sbxarr[19] = -1.5*sbspac;  sbyarr[19] =  0.2887*sbspac; sbtharr[19] = sbac_r + sbrarr[19] * sbac_rnd;
       sbxarr[20] = -1.5*sbspac;  sbyarr[20] = -0.2887*sbspac; sbtharr[20] = sbac_r + sbrarr[20] * sbac_rnd;

       sbxarr[21] =  1.5*sbspac;  sbyarr[21] =  0.866*sbspac; sbtharr[21] = sbac_r + sbrarr[21] * sbac_rnd;
       sbxarr[22] =  1.5*sbspac;  sbyarr[22] = -0.866*sbspac; sbtharr[22] = sbac_r + sbrarr[22] * sbac_rnd;
       sbxarr[23] = -1.5*sbspac;  sbyarr[23] = -0.866*sbspac; sbtharr[23] = sbac_r + sbrarr[23] * sbac_rnd;
       sbxarr[24] = -1.5*sbspac;  sbyarr[24] =  0.866*sbspac; sbtharr[24] = sbac_r + sbrarr[24] * sbac_rnd;

       					// 2 / sqrt(3) -> 2 * 0.5773 -> 1.1547
					
       sbxarr[25] = 0;            sbyarr[25] =  1.1547*sbspac; sbtharr[25] = sbac_r + sbrarr[25] * sbac_rnd;
       sbxarr[26] = 0;            sbyarr[26] = -1.1547*sbspac; sbtharr[26] = sbac_r + sbrarr[26] * sbac_rnd;

       sbxarr[27] =  2 * sbspac;  sbyarr[27] = 0;              sbtharr[27] = sbac_r + sbrarr[27] * sbac_rnd;
       sbxarr[28] =  2 * sbspac;  sbyarr[28] = 0.5773*sbspac;  sbtharr[28] = sbac_r + sbrarr[28] * sbac_rnd;
       sbxarr[29] =  2 * sbspac;  sbyarr[29] = -0.5773*sbspac; sbtharr[29] = sbac_r + sbrarr[29] * sbac_rnd;
       
       sbxarr[30] = -2 * sbspac;  sbyarr[30] = 0;              sbtharr[30] = sbac_r + sbrarr[30] * sbac_rnd;
       sbxarr[31] = -2 * sbspac;  sbyarr[31] = 0.5773*sbspac;  sbtharr[31] = sbac_r + sbrarr[31] * sbac_rnd;
       sbxarr[32] = -2 * sbspac;  sbyarr[32] = -0.5773*sbspac; sbtharr[32] = sbac_r + sbrarr[32] * sbac_rnd;

       sbxarr[33] =  sbspac;      sbyarr[33] =  1.1547*sbspac; sbtharr[33] = sbac_r + sbrarr[33] * sbac_rnd;
       sbxarr[34] = -sbspac;      sbyarr[34] =  1.1547*sbspac; sbtharr[34] = sbac_r + sbrarr[34] * sbac_rnd;
       sbxarr[35] =  sbspac;      sbyarr[35] = -1.1547*sbspac; sbtharr[35] = sbac_r + sbrarr[35] * sbac_rnd;
       sbxarr[36] = -sbspac;      sbyarr[36] = -1.1547*sbspac; sbtharr[36] = sbac_r + sbrarr[36] * sbac_rnd;

                                   // .2887+1.1547 -> 1.4434

       sbxarr[37] =  0.5*sbspac;  sbyarr[37] =  1.4434*sbspac; sbtharr[37] = sbac_r + sbrarr[37] * sbac_rnd;
       sbxarr[38] = -0.5*sbspac;  sbyarr[38] =  1.4434*sbspac; sbtharr[38] = sbac_r + sbrarr[38] * sbac_rnd;
       sbxarr[39] =  0.5*sbspac;  sbyarr[39] = -1.4434*sbspac; sbtharr[39] = sbac_r + sbrarr[39] * sbac_rnd;
       sbxarr[40] = -0.5*sbspac;  sbyarr[40] = -1.4434*sbspac; sbtharr[40] = sbac_r + sbrarr[40] * sbac_rnd;

       sbxarr[41] =  2.5*sbspac;  sbyarr[41] =  0.2887*sbspac; sbtharr[41] = sbac_r + sbrarr[41] * sbac_rnd;
       sbxarr[42] =  2.5*sbspac;  sbyarr[42] = -0.2887*sbspac; sbtharr[42] = sbac_r + sbrarr[42] * sbac_rnd;
       sbxarr[43] = -2.5*sbspac;  sbyarr[43] =  0.2887*sbspac; sbtharr[43] = sbac_r + sbrarr[43] * sbac_rnd;
       sbxarr[44] = -2.5*sbspac;  sbyarr[44] = -0.2887*sbspac; sbtharr[44] = sbac_r + sbrarr[44] * sbac_rnd;

       n_sbac = 45;
       sbr = 2;
   }
   if (sbarr==157) {			// 57 sbacs, all same rotation
	     double r;
       // r = 60;
       // r = 90;
       // sbac_r = 120;
       sbxarr[0] =  0;            sbyarr[0] = 0;             sbtharr[0] = sbac_r + sbrarr[0] * sbac_rnd;
       sbxarr[1] =  sbspac;       sbyarr[1] = 0;             sbtharr[1] = sbac_r + sbrarr[1] * sbac_rnd;
       sbxarr[2] =  0.5*sbspac;   sbyarr[2] = -0.866*sbspac; sbtharr[2] = sbac_r + sbrarr[2] * sbac_rnd;
       sbxarr[3] = -0.5*sbspac;   sbyarr[3] = -0.866*sbspac; sbtharr[3] = sbac_r + sbrarr[3] * sbac_rnd;
       sbxarr[4] = -sbspac;       sbyarr[4] = 0;             sbtharr[4] = sbac_r + sbrarr[4] * sbac_rnd;
       sbxarr[5] = -0.5*sbspac;   sbyarr[5] = 0.866*sbspac;  sbtharr[5] = sbac_r + sbrarr[5] * sbac_rnd;
       sbxarr[6] =  0.5*sbspac;   sbyarr[6] = 0.866*sbspac;  sbtharr[6] = sbac_r + sbrarr[6] * sbac_rnd;

       					//  0.5 / sqrt(3) => 0.2887
				 	//  1   / sqrt(3) -> 0.5773
					
       sbxarr[7] =  0.5*sbspac;   sbyarr[7] = 0.2887*sbspac;   sbtharr[7] = sbac_r + sbrarr[7] * sbac_rnd;
       sbxarr[8] =  0;            sbyarr[8] = 0.5773*sbspac;   sbtharr[8] = sbac_r + sbrarr[8] * sbac_rnd;
       sbxarr[9] =  -0.5*sbspac;  sbyarr[9] = 0.2887*sbspac;   sbtharr[9] = sbac_r + sbrarr[9] * sbac_rnd;
       sbxarr[10] =  0.5*sbspac;  sbyarr[10] = -0.2887*sbspac; sbtharr[10] = sbac_r + sbrarr[10] * sbac_rnd;
       sbxarr[11] =  0;           sbyarr[11] = -0.5773*sbspac; sbtharr[11] = sbac_r + sbrarr[11] * sbac_rnd;
       sbxarr[12] = -0.5*sbspac;  sbyarr[12] = -0.2887*sbspac; sbtharr[12] = sbac_r + sbrarr[12] * sbac_rnd;

       sbxarr[13] =  sbspac;      sbyarr[13] = 0.5773*sbspac;  sbtharr[13] = sbac_r + sbrarr[13] * sbac_rnd;
       sbxarr[14] =  sbspac;      sbyarr[14] = -0.5773*sbspac; sbtharr[14] = sbac_r + sbrarr[14] * sbac_rnd;
       sbxarr[15] = -sbspac;      sbyarr[15] = 0.5773*sbspac;  sbtharr[15] = sbac_r + sbrarr[15] * sbac_rnd;
       sbxarr[16] = -sbspac;      sbyarr[16] = -0.5773*sbspac; sbtharr[16] = sbac_r + sbrarr[16] * sbac_rnd;

       sbxarr[17] =  1.5*sbspac;  sbyarr[17] =  0.2887*sbspac; sbtharr[17] = sbac_r + sbrarr[17] * sbac_rnd;
       sbxarr[18] =  1.5*sbspac;  sbyarr[18] = -0.2887*sbspac; sbtharr[18] = sbac_r + sbrarr[18] * sbac_rnd;
       sbxarr[19] = -1.5*sbspac;  sbyarr[19] =  0.2887*sbspac; sbtharr[19] = sbac_r + sbrarr[19] * sbac_rnd;
       sbxarr[20] = -1.5*sbspac;  sbyarr[20] = -0.2887*sbspac; sbtharr[20] = sbac_r + sbrarr[20] * sbac_rnd;

       sbxarr[21] =  1.5*sbspac;  sbyarr[21] =  0.866*sbspac; sbtharr[21] = sbac_r + sbrarr[21] * sbac_rnd;
       sbxarr[22] =  1.5*sbspac;  sbyarr[22] = -0.866*sbspac; sbtharr[22] = sbac_r + sbrarr[22] * sbac_rnd;
       sbxarr[23] = -1.5*sbspac;  sbyarr[23] = -0.866*sbspac; sbtharr[23] = sbac_r + sbrarr[23] * sbac_rnd;
       sbxarr[24] = -1.5*sbspac;  sbyarr[24] =  0.866*sbspac; sbtharr[24] = sbac_r + sbrarr[24] * sbac_rnd;

       					// 2 / sqrt(3) -> 2 * 0.5773 -> 1.1547
					
       sbxarr[25] = 0;            sbyarr[25] =  1.1547*sbspac; sbtharr[25] = sbac_r + sbrarr[25] * sbac_rnd;
       sbxarr[26] = 0;            sbyarr[26] = -1.1547*sbspac; sbtharr[26] = sbac_r + sbrarr[26] * sbac_rnd;

       sbxarr[27] =  2 * sbspac;  sbyarr[27] = 0;              sbtharr[27] = sbac_r + sbrarr[27] * sbac_rnd;
       sbxarr[28] =  2 * sbspac;  sbyarr[28] = 0.5773*sbspac;  sbtharr[28] = sbac_r + sbrarr[28] * sbac_rnd;
       sbxarr[29] =  2 * sbspac;  sbyarr[29] = -0.5773*sbspac; sbtharr[29] = sbac_r + sbrarr[29] * sbac_rnd;
       
       sbxarr[30] = -2 * sbspac;  sbyarr[30] = 0;              sbtharr[30] = sbac_r + sbrarr[30] * sbac_rnd;
       sbxarr[31] = -2 * sbspac;  sbyarr[31] = 0.5773*sbspac;  sbtharr[31] = sbac_r + sbrarr[31] * sbac_rnd;
       sbxarr[32] = -2 * sbspac;  sbyarr[32] = -0.5773*sbspac; sbtharr[32] = sbac_r + sbrarr[32] * sbac_rnd;

       sbxarr[33] =  sbspac;      sbyarr[33] =  1.1547*sbspac; sbtharr[33] = sbac_r + sbrarr[33] * sbac_rnd;
       sbxarr[34] = -sbspac;      sbyarr[34] =  1.1547*sbspac; sbtharr[34] = sbac_r + sbrarr[34] * sbac_rnd;
       sbxarr[35] =  sbspac;      sbyarr[35] = -1.1547*sbspac; sbtharr[35] = sbac_r + sbrarr[35] * sbac_rnd;
       sbxarr[36] = -sbspac;      sbyarr[36] = -1.1547*sbspac; sbtharr[36] = sbac_r + sbrarr[36] * sbac_rnd;

                                   // .2887+1.1547 -> 1.4434

       sbxarr[37] =  0.5*sbspac;  sbyarr[37] =  1.4434*sbspac; sbtharr[37] = sbac_r + sbrarr[37] * sbac_rnd;
       sbxarr[38] = -0.5*sbspac;  sbyarr[38] =  1.4434*sbspac; sbtharr[38] = sbac_r + sbrarr[38] * sbac_rnd;
       sbxarr[39] =  0.5*sbspac;  sbyarr[39] = -1.4434*sbspac; sbtharr[39] = sbac_r + sbrarr[39] * sbac_rnd;
       sbxarr[40] = -0.5*sbspac;  sbyarr[40] = -1.4434*sbspac; sbtharr[40] = sbac_r + sbrarr[40] * sbac_rnd;

       sbxarr[41] =  2.5*sbspac;  sbyarr[41] =  0.2887*sbspac; sbtharr[41] = sbac_r + sbrarr[41] * sbac_rnd;
       sbxarr[42] =  2.5*sbspac;  sbyarr[42] = -0.2887*sbspac; sbtharr[42] = sbac_r + sbrarr[42] * sbac_rnd;
       sbxarr[43] = -2.5*sbspac;  sbyarr[43] =  0.2887*sbspac; sbtharr[43] = sbac_r + sbrarr[43] * sbac_rnd;
       sbxarr[44] = -2.5*sbspac;  sbyarr[44] = -0.2887*sbspac; sbtharr[44] = sbac_r + sbrarr[44] * sbac_rnd;

                                   // 12 more 
				   
       sbxarr[45] =  1.5*sbspac;  sbyarr[45] =  1.4434*sbspac; sbtharr[45] = sbac_r + sbrarr[45] * sbac_rnd;
       sbxarr[46] =  2.0*sbspac;  sbyarr[46] =  1.1547*sbspac; sbtharr[46] = sbac_r + sbrarr[46] * sbac_rnd;
       sbxarr[47] =  2.5*sbspac;  sbyarr[47] =  0.866*sbspac;  sbtharr[47] = sbac_r + sbrarr[47] * sbac_rnd;
       sbxarr[48] =  1.5*sbspac;  sbyarr[48] = -1.4434*sbspac; sbtharr[48] = sbac_r + sbrarr[48] * sbac_rnd;
       sbxarr[49] =  2.0*sbspac;  sbyarr[49] = -1.1547*sbspac; sbtharr[49] = sbac_r + sbrarr[49] * sbac_rnd;
       sbxarr[50] =  2.5*sbspac;  sbyarr[50] = -0.866*sbspac;  sbtharr[50] = sbac_r + sbrarr[50] * sbac_rnd;

       sbxarr[51] = -1.5*sbspac;  sbyarr[51] =  1.4434*sbspac; sbtharr[51] = sbac_r + sbrarr[51] * sbac_rnd;
       sbxarr[52] = -2.0*sbspac;  sbyarr[52] =  1.1547*sbspac; sbtharr[52] = sbac_r + sbrarr[52] * sbac_rnd;
       sbxarr[53] = -2.5*sbspac;  sbyarr[53] =  0.866*sbspac;  sbtharr[53] = sbac_r + sbrarr[53] * sbac_rnd;
       sbxarr[54] = -1.5*sbspac;  sbyarr[54] = -1.4434*sbspac; sbtharr[54] = sbac_r + sbrarr[54] * sbac_rnd;
       sbxarr[55] = -2.0*sbspac;  sbyarr[55] = -1.1547*sbspac; sbtharr[55] = sbac_r + sbrarr[55] * sbac_rnd;
       sbxarr[56] = -2.5*sbspac;  sbyarr[56] = -0.866*sbspac;  sbtharr[56] = sbac_r + sbrarr[56] * sbac_rnd;

       n_sbac = 57;
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
 }

  if (notinit(sbac_dia)) sbac_dia = 300;
  if (ninfo >= 2) ncfprintf (stderr,"#\n# sbac dia %-5.3g dens %-5.3g nnd %-5.3g regularity %-5.3g coverage %-5.3g\n",
		  sbac_dia, sbac_dens, 1/sqrt(sbac_dens*1e-6), getn(sbac,REGU),(sbac_dia*sbac_dia)/4*PI*sbac_dens*1e-6);

  // onplot_dsgc_movie_init();		/* initialize dsgc movie stuff */
  // onplot_movie_init();		/* initialize onplot_movie stuff */

  gcdistnod  = 582;
  pickden[dsgc] = 0; //941;       	/* pick one dendrite to allow connections */
  if (notinit(ath_dia)) ath_dia = 0.4;	/* diameter of initial thin axon segment */

  if (n_dsgc>0) make_ct(dsgc);		/* make dsgcs if user specifies */
   
  if (notinit(dvrev)) dvrev = -0.06;	/* default dbp1 vrev (dens_dbp1.n) */
  if (notinit(dvst))  dvst  = -0.06;	/* default vstart */
  if (notinit(db_vs))  db_vs  = -0.056;	/* default dbp1 vstart */
  
  if (!notinit(g_dbp1_dsgc)) setsv(dbp1,SCOND,2,g_dbp1_dsgc);
  if (!notinit(g_dbp1_sbac)) setsv(dbp1,SCOND,3,g_dbp1_sbac);
  if (!notinit(n_dbp1_sbac)) setsv(dbp1,SVNOISE,3,n_dbp1_sbac); /* dbp1 -> sbac noise */
  if (!notinit(v_dbp1_sbac)) setsv(dbp1,SGAIN,3,v_dbp1_sbac);   /* dbp1 -> sbac ves release gain */

  if (!notinit(n_dbp1_ams))  setsv(dbp1,SVNOISE,7,n_dbp1_ams);  /* dbp1 -> ams noise */
  if (!notinit(g_dbp1_ams))  setsv(dbp1,SCOND,  7,g_dbp1_ams);

  if (!notinit(g_dbp2_dsgc)) setsv(dbp2,SCOND,  2,g_dbp2_dsgc);

  if (make_dbp2_sbac) {
     if      (!notinit(g_dbp2_sbac)) setsv(dbp2,SCOND,3,g_dbp2_sbac);
     else if (!notinit(g_dbp1_sbac)) setsv(dbp2,SCOND,3,g_dbp1_sbac);
     if (!notinit(n_dbp2_sbac)) setsv(dbp2,SVNOISE,3,n_dbp2_sbac); /* dbp2 -> sbac noise */
     if (!notinit(v_dbp2_sbac)) setsv(dbp2,SGAIN,  3,v_dbp2_sbac);   /* dbp2 -> sbac ves release gain */
  }

  	/* inhibitory synapses, for nval file */

  if (!notinit(g_am_dsgc))   setsv(am,  SCOND, 2, g_am_dsgc);
  if (!notinit(g_am2_dsgc))  setsv(am2, SCOND, 2, g_am2_dsgc);
  if (!notinit(g_sbac_dsgc)) setsv(sbac,SCOND, 2, g_sbac_dsgc);
  if (!notinit(g_sbac_sbac)) setsv(sbac,SCOND, 4, g_sbac_sbac);
  if (!notinit(n_sbac_sbac)) setsv(sbac,SVNOISE,4,n_sbac_sbac);

  if (!notinit(m_sbac_sbac)) setsv(sbac,SMRRPOOL,4,m_sbac_sbac); /* sbac - sbac synaptic depression */
  if (!notinit(r_sbac_sbac)) setsv(sbac,SMAXRATE,4,r_sbac_sbac); /* sbac - sbac synaptic depression */

  if (!notinit(g_ams_sbac))   setsv(ams,  SCOND, 2, g_ams_sbac);
  if (!notinit(g_ams_sbac))   setsv(ams,  SCOND, 3, g_ams_sbac);
  if (!notinit(f_ams_sbac))   setsv(ams,  SDUR,  2, f_ams_sbac);
  if (!notinit(f_ams_sbac))   setsv(ams,  SDUR,  3, f_ams_sbac);
  if (!notinit(n_ams_sbac))   setsv(ams,  SVNOISE,2, n_ams_sbac);
  if (!notinit(n_ams_sbac))   setsv(ams,  SVNOISE,3, n_ams_sbac);
  if (!notinit(ams_synanpi))  setsv(ams,  SYNANPI,3, ams_synanpi);      // slow diffusive connection
  if (!notinit(ams_synanpo))  setsv(ams,  SYNANPO,3, ams_synanpo);


  if(notinit(dsgc_prefdir)) dsgc_prefdir=0;
  if(notinit(revdir)) revdir=0;

  if (strcmp(sbac_file,"morph_sb1")==0) {
    setsv (sbac,SYNANPI,3,75);                    /* don't connect close to soma */
  }
  else if (strcmp(sbac_file,"morph_sbac3b")==0) {
    setsv (sbac,SYNANPI,3,20);                    /* don't connect close to soma */
  }

  if (notinit(dtreedia)) dtreedia = 500;
  somadist  = dtreedia * 0.3;

  /* 					// don't make am's
   
  if (n_dsgc > 0) {
    amx = gcxarr[0] * 0.5;
    amy = gcyarr[0] * 0.5;
  } 
  setn(am, DTREEDIA,dtreedia);
  setn(am2,DTREEDIA,dtreedia);
  n = 8;
  theta_incr = 360.0/n;
  i = make_am_cell (am,  theta=0,     amx,  amy, i=0, somadist);
  for (n=1,t=theta_incr; t<180; t+= theta_incr,n++) {
       i = make_am_cell (am, theta=t, amx,  amy, i,   somadist);
  }
  if (n_am < 0) n_am = n;
  i = make_am_cell (am2, theta=theta_incr/2,    amx,  amy, i=0, somadist);
  for (n=1,t=theta_incr*1.5; t<180; t+= theta_incr,n++) {
       i = make_am_cell (am2, theta=t, amx,  amy, i,   somadist);
  }
  if (n_am2 < 0) n_am2 = n;
  /* */

  // display_z(zmax=-24, zmin=-15);	    /* exclude dsgc Off-arborization layer */
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

    ndens[sbac][cn=1] = 0;          // ndens array initialized to 0
    				    // set cn 1 to use sbac_densfile
				    //
    if (sbarr== -2) {		    // set cn 2 to use sbac_densfile2
	    ndens[sbac][cn=2] = 1;        
    } else
    				    //  set sb1mul > 1 to increase inhibition to other sbacs
    if (sbarr==3) {
       ndens[sbac][cn=2] = 1;         // set cells 2,3 to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=3] = 1;   
    } 
    if (sbarr==7) {
       ndens[sbac][cn=1] = 1;         // set cells 1,7 to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=7] = 1;         
    } 
    if (sbarr==9) {
       ndens[sbac][cn=1] = 1;         // set cells 1,9 to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=9] = 1; 
    } 
    if (sbarr==11) {
       ndens[sbac][cn=1] = 1;         // set cells 1,11 to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=11] = 1; 
    } 
    else if (sbarr==107) {
       ndens[sbac][cn=2] = 1;         // set cell 2-7 to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=3] = 1;   
       ndens[sbac][cn=4] = 1;   
       ndens[sbac][cn=5] = 1;   
       ndens[sbac][cn=6] = 1;   
       ndens[sbac][cn=7] = 1;   
    }
    else if (sbarr==123) {
       ndens[sbac][cn=3] = 1;         // set cells to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=4] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=6] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=7] = 1;         // set cells to use sbac_densfile2
       
       ndens[sbac][cn=14] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=15] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=16] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=17] = 1;         // set cells to use sbac_densfile2

       ndens[sbac][cn=18] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=18] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=20] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=21] = 1;         // set cells to use sbac_densfile2

       ndens[sbac][cn=22] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=23] = 1;         // set cells to use sbac_densfile2

    } else if (sbarr==157) {

       ndens[sbac][cn=38] = 1;         // set cells to use sbac_densfile2, inhib input cond set by sb1mul
       ndens[sbac][cn=39] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=40] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=41] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=42] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=43] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=44] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=45] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=46] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=47] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=48] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=49] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=50] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=51] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=52] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=53] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=54] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=55] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=56] = 1;         // set cells to use sbac_densfile2
       ndens[sbac][cn=57] = 1;         // set cells to use sbac_densfile2
    }

    // ndens[sbac][cn=3] = 0;          // set cn 3 to use sbac_densfile
}
	    
/*--------------------------------------------------------*/

void runonexit (void)
{
       if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/

double ixoff, iyoff;
double mvbrt1;
double mvbrt2;
double *rndarr = NULL;

double move_stim(double stimtime,double barwidth,double theta, double scontrast, int direction, double mask)
{
    static int runyet=0;
    double lbar, rbar, tbar, irad, orad, idia, odia;
    double inten, start, dur, wavel;
    double cellrad, celldia;
    double s,send;
    double tf;

 celldia = max(xarrsiz,yarrsiz) * 1.5;
 cellrad = max(xarrsiz,yarrsiz) * 0.5; 
 // cellrad 446 for array 500 x 300
 
 // stim_spot(100, 100,0, inten=minten,  start=0.02, dur=0.05, wavel=1, 0);

 if (stimtype==1) {			// move bar to & fro
   lbar = -cellrad*rloc - barwidth - blur_ssize;	// blur_ssize added to be compat with stimtype==13
   rbar =  cellrad*rloc + barwidth + blur_ssize;
   if (revdir) { tbar = lbar; lbar = rbar; rbar = tbar; }
   if (notinit(velocity)) velocity = 2000;
   mvbrt1 = movebar (stimtime,                   stimx,stimy, lbar, rbar, barwidth,barlength,theta,velocity,scontrast);
   mvbrt2 = movebar (mvbrt1+poststimdur+stimtime,stimx,stimy, rbar, lbar, barwidth,barlength,theta,velocity,scontrast);

   // stim_spot(spotdia=1,2000,0,scontrast,-predur,2);	/* for wacs */
   if (runyet==0 && mask>0) {
       stim_bar(200,1000, -100,0, theta, inten=minten, start=0, dur=1, wavel=1, mask);
       runyet = 1;
   }
 }
 else if (stimtype==2) {		// move annulus out & in  = CF, then CP

   irad = 20;
   orad = cellrad + barwidth*0.5;
   if (revdir) { tbar = irad; irad = orad; orad = tbar; }
   if (notinit(velocity)) velocity = 2000;
   if (notinit(ncycles)) ncycles = 4;
   mvbrt1 = moveannulus (stimtime,                   0,0, irad, orad, barwidth,velocity,scontrast);
   mvbrt2 = moveannulus (mvbrt1+poststimdur+stimtime,0,0, orad, irad, barwidth,velocity,scontrast);
   if (runyet==0 && mask>0) {
  //     stim_spot(idia, 0,0, inten=minten,  start=0, dur=1, wavel=1, mask);
       stim_bar(200,1000, -100,0, theta, inten=minten, start=0, dur=1, wavel=1, mask);
       runyet = 1;
   }
 } else if (stimtype==3) {		// move sineann

  if (notinit(velocity)) velocity = 2000;
  tf = velocity/barwidth;
  dur = ncycles / tf;
  movesineann (0, 0, direction=1, annrad, 0,   0,barwidth,tf,0,1.0,minten,scontrast, makenv, swaveshape, stimtime,dur);
  movesineann (0, 0, direction=2, annrad, 0, 180,barwidth,tf,0,1.0,minten,scontrast, makenv, swaveshape, 2*stimtime+dur,dur);
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

  if (notinit(velocity)) velocity = 2000;
  if (notinit(ncycles)) ncycles = 4;
  tf = velocity/barwidth;  // speriod = barwidth
  dur = ncycles / tf;
  movegrating (0, 0, theta, sphase=0, barwidth, tf, drift=1, minten, scontrast, swaveshape, stimtime,dur);
  movegrating (0, 0, theta, sphase=0, barwidth, tf, drift=-1, minten, scontrast, swaveshape, 2*stimtime+dur,dur);
  mvbrt1 = dur;
  mvbrt2 = 2 * dur + stimtime;

 } else if (stimtype==6) {		// flashed annulus, proximal, distal

   if (notinit(orad1)) orad1 = 100;
   if (notinit(irad1)) irad1 = 0;
   if (notinit(irad2)) irad2 = 0;
   if (notinit(orad2)) orad2 = 0;

   if (orad1 > 0) {
     stim_spot (odia=2*orad1, 0,  0,  scontrast, stimtime, spotdur);		// flashed proximal annulus
     if (irad1 > 0)
       stim_spot(idia=2*irad1, 0,  0, -scontrast, stimtime, spotdur);		// flashed proximal annulus
   }

   mvbrt1 = stimtime+prestimdur;
   if (orad1 > 0) mvbrt1 += spotdur;
   mvbrt2 = mvbrt1;

   if (orad2 > 0) {
      stim_spot(odia=2*orad2, 0, 0,  scontrast, 2*stimtime+prestimdur+spotdur, spotdur);  // flashed distal annulus
      stim_spot(idia=2*irad2, 0, 0, -scontrast, 2*stimtime+prestimdur+spotdur, spotdur);
      mvbrt2 += stimtime+spotdur;
   }
 } else if (stimtype==7) {		// just vclamp series, no stim
					//
      mvbrt2 = spotdur; 

 } else if (stimtype==8) {		// just iclamp series, no stim
					//
      mvbrt2 = spotdur; 

 } else if (stimtype==9) {		// move bar to & fro inside spot mask
	 				// like stimtype 1, but with mask
	int stimchan, invert;
        double mask, stimtime2;
					// allow background on stimulus channel 0
   stimchan = 1;			// set stimulus channel 1 for masked moving bar
  // lbar = -cellrad*rloc - barwidth/2;
  // rbar =  cellrad*rloc + barwidth/2;
   lbar = -mask_dia/2*rloc - barwidth/2 + c1_somax + mask_x;
   rbar =  mask_dia/2*rloc + barwidth/2 + c1_somax + mask_x;
   if (revdir) { tbar = lbar; lbar = rbar; rbar = tbar; }
   if (notinit(velocity)) velocity = 2000;
   mvbrt1 = movebar (stimtime, stimx,stimy, lbar, rbar, barwidth,barlength,theta,velocity,scontrast,stimchan);

   stimtime2 = mvbrt1+poststimdur;
   mvbrt2 = movebar (stimtime2,stimx,stimy, rbar, lbar, barwidth,barlength,theta,velocity,scontrast,stimchan);
   // stim_spot(spotdia=1,2000,0,scontrast,-predur,2);	/* for wacs */
   if (runyet==0 && mask_dia>0) {
       stim_spot(mask_dia, c1_somax+mask_x, c1_somay+mask_y, inten=0, start=0, dur=10, mask=1.0, stimchan, invert=1);
       runyet = 1;
   }
 } else if (stimtype==10) {                 // checkerboard
          int npixels, nframes;
          double xoffset,yoffset;
          double orient=0;

    stimdur = 1;
    npixels = 16;
    stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
    tf=15, 0, scontrast*2, stimtime, stimdur, &rndarr, &nframes, rseed);
    mvbrt2 = stimtime + stimdur;

    //   write out file containing random stimulus
    //
    // sprintf (checkerboard_file,"rnd_chck_file_%d",rseed);
    // if ((fchk=fopen(checkerboard_file,"w"))==NULL) { /* open file */
    //     ncfprintf (stderr,"Error creating checkerboard file\n");
    // } else {
    //     ncfprintf (stderr,"Creating checkerboard file %s\n",checkerboard_file);
    //     fwrite (rndarr,sizeof(double),npixels*npixels*nframes,fchk);
    // }
    efree(rndarr);
 }

 /* **************************************************************************************** */

 else if (stimtype==11) {                 // moving bar + checkerboard
         int npixels, nframes;
         double xoffset,yoffset;
         double orient=0, roffset=0;
         double time2;

    lbar = -(celldia + barwidth) * 0.6 + roffset;
    rbar =  (celldia + barwidth) * 0.6 + roffset;
    if (notinit(velocity)) velocity = 2000;
    // time1 = movebar (stimtime+noise_dur, 0,0, lbar, rbar, barwidth,barlength,theta,velocity,scontrast);
    // time2 = movebar (time1+stimtime,     0,0, rbar, lbar, barwidth,barlength,theta,velocity,scontrast);
    time2 = movebar (stimtime+noise_dur, 0,0, lbar, rbar, barwidth,barlength,theta,velocity,scontrast);

    stimdur = time2 + noise_dur * 0.5;
    // stimdur = 1.0;

    if (noise_dur > 0) {
        npixels = 16;
        stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
                     tf=15, 0, scontrast*0.5, stimtime, stimdur, &rndarr, &nframes, rseed);
        efree(rndarr);
     }
     mvbrt2 = stimdur;
 }

 /* **************************************************************************************** */

 else if (stimtype==12) {                 // flashed spot with blur d-o-g, outer and inner dia  (annulus)
	 	//  retsim --expt_dsgc_sbac_bar --use_stimfile 1 --makestim 1 ..  // make the stim file
		//  retsim --expt_dsgc_sbac_bar --use_stimfile 1 --makestim 0 ... // run model with the stim file
 	  int stimchan;
	  double incr;
          const char *fname;
          static char stimfil[150];
	  double spottime, postspotdur=0.05;
	  int s;
          double sd, spdia, annulus_idia;

     //   if tfreq==0, to set spot interval, set spotint
     //   to set number of cycles, set ncycles
     //   to set dia increment, set spotdia_incr
     //   to make annulus, set annulus_width 
     //   to make annulus outer dia constant, set annulus_odia > 0
     //   to make repeating sine/square modulation, set tfreq > 0
     //   if tfreq==0, to add negative spot after pos spot is done, set scontrastn > 0

     if (notinit (blur_csize)) blur_csize = 30;
     if (notinit (blur_ssize)) blur_ssize = 120;
     if (notinit (spotdia_incr)) spotdia_incr = 50;
     if (notinit (spotdia_max)) spotdia_max = 300;
     if (notinit (annulus_odia)) annulus_odia = 0;	// if > 0, sets annulus outer diameter
     if (notinit (annulus_width)) annulus_width = 0;
     if (notinit(ncycles)) ncycles = 2;
     if (use_stimfile) {    				     // if d-o-g stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_spot_%d_%g_%g_%g_%g_%g_%g_%g_%g_%g_%g";	     // set the stimulus file name
       sprintf (stimfil,fname,stimtype,sbac_r,spotdia,scontrast,minten,blur_csize,blur_ssize,surr_weight,
		       spotdia_incr,spotdia_max,annulus_width);
       sprintf (stimfil,"%s_%s",stimfil,stim_fnum);
       stim_file (stimfil);                                  // set the stimulus file 
     }

     if (tfreq>0) spotint = ncycles / tfreq;
     if (spotdia_incr <= 10) spotdia_incr = 10; 
     if (notinit(nstim)) {
	     nstim = spotdia_max / spotdia_incr;
             if (annulus_width > 0) nstim++;
     }
     for (s=0; s<nstim; s++) {
         spottime = stimtime + s*(spotint);
	 spdia = (s+1) * spotdia_incr; 
         if (use_stimfile) {    				    // if d-o-g stimulus (and --use_stimfile 1 )
             stim_blur (blur_csize,center_weight);                  // set the blur array
	 }
         stim_backgr(minten,start=simtime);	  	     /* background */
	 if (annulus_odia > 0) sd = annulus_odia;
	 else                  sd = spdia;
         if (tfreq==0) {
	         stim_spot(sd, stimx, stimy, inten=scontrast, start=spottime, dur=spotdur, mask=0, stimchan=0);
	     if (scontrastn > 0) 
                 stim_spot(sd, stimx, stimy, inten= -scontrastn, start=spottime+spotdur, dur=spotdur, mask=0,stimchan=0);
	 }
	 else    spot_sine_ncycle(sd, stimx, stimy, tfreq, 0.0, twaveshape, 1.0, scontrast,spottime,ncycles);

         if (annulus_width > 0) {
	     annulus_idia = spdia - annulus_width;
             if (annulus_idia > 0) {
                 if (tfreq==0) {
		    stim_spot(annulus_idia, stimx, stimy, inten= -scontrast, start=spottime, 
			     				dur=spotdur, mask=0, stimchan=0);
	            if (scontrastn > 0) 
		     stim_spot(annulus_idia, stimx, stimy, inten= scontrastn, start=spottime+spotdur, 
			     				dur=spotdur, mask=0, stimchan=0);
	         } else spot_sine_ncycle(annulus_idia, stimx, stimy, tfreq, 0.0, twaveshape, 1.0, -scontrast,
				 			spottime,ncycles);
	     }
	 }
         if (surr_weight != 0) {
            if (use_stimfile) {    				     // if d-o-g stimulus (and --use_stimfile 1 )
               stim_blur (blur_ssize,surr_weight);                   // set the blur array
               if (tfreq==0) {
                   spot_ramp(sd, stimx, stimy, 1.0, 0, scontrast, start=spottime, spotdur, surr_tau, incr= 0.1);
                   if (annulus_width > 0) {
                      spot_ramp(annulus_idia, stimx, stimy, 1.0, 0, -scontrast, start=spottime, spotdur,
				  			surr_tau, incr= 0.1);
	           }
	           if (scontrastn > 0) {
                      spot_ramp(sd, stimx, stimy, 1.0, 0, -scontrastn, start=spottime+spotdur, spotdur,
				  surr_tau, incr= 0.1);
                      if (annulus_width > 0) {
                         spot_ramp(annulus_idia, stimx, stimy, 1.0, 0, scontrastn, start=spottime+spotdur, spotdur,
				  			surr_tau, incr= 0.1);
	             }
	          }
	       } else {
		    					// digital low-pass filter for twaveshape == 1
	            spot_sine_ncycle(sd, stimx, stimy, tfreq, surr_tau, twaveshape, 1.0, scontrast,
				 			spottime+surr_delay,ncycles);
                    if (annulus_idia > 0) {
	               spot_sine_ncycle(annulus_idia, stimx, stimy, tfreq, surr_tau, twaveshape, 1.0, -scontrast,
				 			spottime+surr_delay,ncycles);
	            } 
	       }
            }
         }
     }
      mvbrt2 = s * spotint;
 }

 /* **************************************************************************************** */

 else if (stimtype==13) {                 // moving bar with diff-of-gaussians blur
	  double offsetx, delay, rstep, tstep;
          const char *fname;
          static char stimfil[150];

     dur = 2;
     if (notinit (blur_csize)) blur_csize = 30;
     if (notinit (blur_ssize)) blur_ssize = 120;
     if (notinit(velocity)) velocity = 2000;
     if (notinit(stim_tstep)) {
	 if      (velocity >= 5000) stim_tstep = 0.0005;
	 else if (velocity >= 1000) stim_tstep = 0.001;
	 else                  stim_tstep = 0.002;
     }
     // delay = surr_delay/2;
     rstep = velocity * stim_tstep;
     delay = surr_delay;
     if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_bar_%d_%g_%g_%g_%g_%g_%g";	     // set the stimulus file name
       sprintf (stimfil,fname,stimtype,scontrast,barwidth,barlength,velocity,blur_csize,blur_ssize);
       sprintf (stimfil,"%s_%s_%g_%d",stimfil,stim_fnum,sbac_rnd,drseed);
       stim_file (stimfil);                                  // set the stimulus file 
       stim_blur (blur_csize,center_weight);                 // set the blur array
     }
     stim_backgr(minten,start=simtime);	  	     /* background */

     lbar = -(cellrad*rloc + barwidth) - blur_ssize;
     rbar =  (cellrad*rloc + barwidth) + blur_ssize;
     if (revdir) { tbar = lbar; lbar = rbar; rbar = tbar; }
     mvbrt1 = movebar (stimtime,                   stimx,stimy, lbar, rbar, barwidth,barlength,
		     						theta,velocity,scontrast,rstep);
     mvbrt2 = movebar (mvbrt1+poststimdur+stimtime,stimx,stimy, rbar, lbar, barwidth,barlength,
		     						theta,velocity,scontrast,rstep);
     if (surr_weight != 0) { 
        if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
           stim_blur (blur_ssize,surr_weight);               // set the blur array
           offsetx = velocity * delay;
           mvbrt1 = movebar (stimtime,                   stimx-offsetx,stimy, lbar, rbar, barwidth,barlength,
			   					theta,velocity,scontrast,rstep);
           mvbrt2 = movebar (mvbrt1+poststimdur+stimtime,stimx+offsetx,stimy, rbar, lbar, barwidth,barlength,
			   					theta,velocity,scontrast,rstep);
        }
    }
 }
 
 /* **************************************************************************************** */

 else if (stimtype==14) {		// move annulus out CF, then in CP, with diff-of-gaussian blur
	  double offsetx, delay;
          const char *fname;
          static char stimfil[150];

     barwidth = 20;
     if (notinit (blur_csize)) blur_csize = 30;
     if (notinit (blur_ssize)) blur_ssize = 120;
     if (notinit(velocity)) velocity = 2000;
     // delay = surr_delay/2;
     delay = surr_delay;
     if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_ann_%d_%g_%g_%g_%g_%g";	     // set the stimulus file name
       sprintf (stimfil,fname,stimtype,scontrast,barwidth,velocity,blur_csize,blur_ssize);
       sprintf (stimfil,"%s_%s_%d_%g",stimfil,stim_fnum,drseed,sbac_rnd);
       stim_file (stimfil);                                  // set the stimulus file 
       stim_blur (blur_csize,center_weight);                 // set the blur array
     }
     stim_backgr(minten,start=simtime);	  	     /* background */

   irad = -20;
   orad = cellrad + 2*barwidth + blur_ssize/2;
   if (revdir) { tbar = irad; irad = orad; orad = tbar; }
   mvbrt1 = moveannulus (stimtime,                   0,0, irad, orad, barwidth,velocity,scontrast);
   mvbrt2 = moveannulus (mvbrt1+poststimdur+stimtime,0,0, orad, irad, barwidth,velocity,scontrast);

   if (runyet==0 && mask_dia>0) {
       stim_spot(mask_dia, 0,0, inten=minten,  start=0, dur=1, wavel=1, mask=1.0);
  //     stim_bar(200,1000, -100,0, theta, inten=minten, start=0, dur=1, wavel=1, mask=1.0);
       runyet = 1;
   }
   if (surr_weight != 0) {
      if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
          stim_blur (blur_ssize,surr_weight);                  // set the blur array
      }
      // offsetx = velocity * delay;
      offsetx = 0;
      mvbrt1 = moveannulus (stimtime+delay,                   0,0, irad-offsetx, orad, barwidth,velocity,scontrast);
      mvbrt2 = moveannulus (mvbrt1+poststimdur+stimtime+delay,0,0, orad+offsetx, irad, barwidth,velocity,scontrast);
   }
 }
 
 /* **************************************************************************************** */

 else if (stimtype==15) {		// move concentric sine- or square-wave grating out & in: CF, then CP, 
	 				//  with diff-of-gaussian blur.
 	  int stimchan, invert;		// This stimulus generates a larger input for the central bipolar cells
	  double offsetx, delay;	//  because when the sine wave peak leaves the center, its blurred
          const char *fname;		//  result gets energy from both sides of the center.
          static char stimfil[150];

     if (notinit (blur_csize)) blur_csize = 30;
     if (notinit (blur_ssize)) blur_ssize = 100;
     if (notinit(velocity)) velocity = 2000;
     if (notinit(ncycles)) ncycles = 4;
     // delay = surr_delay/2;
     delay = surr_delay;
     if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_sineann_%d_%d_%g_%g_%g_%g_%g_%g_%g";	     // set the stimulus file name
       sprintf (stimfil,fname,stimtype,swaveshape,barwidth,ncycles,velocity,
		       scontrast,minten,blur_csize,blur_ssize);
       sprintf (stimfil,"%s_%s",stimfil,stim_fnum);
       stim_file (stimfil);                                  // set the stimulus file 
       stim_blur (blur_csize,center_weight);                 // set the blur array
     }
     stim_backgr(minten,start=simtime);	  	     /* background */

  tf = velocity/barwidth;
  dur = ncycles / tf;
  movesineann (0, 0, direction=1, annrad, 0,   0,barwidth,tf,0,1.0,
		  	minten,scontrast, makenv, swaveshape, stimtime,dur);
  movesineann (0, 0, direction=2, annrad, 0, 180,barwidth,tf,0,1.0,
		  	minten,scontrast, makenv, swaveshape, 2*stimtime+dur,dur);
  mvbrt1 = dur;
  mvbrt2 = stimtime + dur;

  if (surr_weight != 0) {
     if (use_stimfile) {    				      // if moving bar stimulus (and --use_stimfile 1 )
         stim_blur (blur_ssize,surr_weight);                  // set the blur array
         movesineann (0, 0, direction=1, annrad, 0,   0,barwidth,tf,0,1.0,
		  	minten,scontrast, makenv, swaveshape, stimtime+delay,dur);
         movesineann (0, 0, direction=2, annrad, 0, 180,barwidth,tf,0,1.0,
		  	minten,scontrast, makenv, swaveshape, 2*stimtime+dur+delay,dur);
     }
  }

  if (runyet==0 && mask_dia>0) {
       if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
            stim_blur (0,1);                                 // remove blur for central mask
       }
       unmaskthr2 = -0.001;
       // stim_backgr(minten-0.01, mask, start=simtime);	  	     /* background for mask */
       stim_spot(mask_dia, 0,0, inten=minten,  start=stimtime, stimtime+2*dur+delay, mask=1.0);
  //     stim_bar(200,1000, -100,0, theta, inten=minten, start=stimtime, dur, wavel=1, mask=1.0);
       runyet = 1;
   }

  mvbrt1 = dur+delay;			    // time for each pass
  mvbrt2 = stimtime + 2*(dur + delay);    // total time for display of stimulus
 }

 /* **************************************************************************************** */

 else if (stimtype==16) {	// drifting linear grating, sine or square, with diff-of-gaussian blur
	 			// can set velocity or tfreq, and ncycles
          int s, dir;
	  int set_vel = 0, set_ncycles = 0;
	  double bwidth, tf, sphase;
	  double celldia = 150; 
	  double offsetx, delay, dstart;
          const char *fname;
          static char stimfil[150];
	  double stim_incr;


    if (notinit (blur_csize)) blur_csize = 30;
    if (notinit (blur_ssize)) blur_ssize = 120;
    if (!notinit(velocity)) set_vel = 1;
    if (!notinit(ncycles)) set_ncycles = ncycles;
    if (notinit (min_sp)) min_sp = 75;
    if (notinit (max_sp)) max_sp = 500;
    if (notinit(incr_sp)) incr_sp = 25;

    // delay = surr_delay/2;
    delay = surr_delay;

    if (ninfo>=2) {
	if (set_vel) 
 	     ncfprintf(stderr,"# drifting grating, velocity %g, vary temporal, spatial freq\n", velocity);
       else  ncfprintf(stderr,"# drifting grating, temporal freq %g, vary velocity,spatial freq\n", tfreq);
    }
    if (notinit(nstim))                      nstim = 1 + (max_sp - min_sp) / incr_sp;
    if (min_sp + (nstim-1)*incr_sp > max_sp) nstim = 1 + (max_sp - min_sp) / incr_sp;
    max_sp = min_sp + (nstim-1) * incr_sp;
    fprintf (stderr,"# min_sp %g max_sp %g incr_sp %g nstim %d\n",min_sp,max_sp,incr_sp,nstim);

    if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_drifting_grating_%d_%d_%d_%g_%g_%g_%g_%g_%g_%g_%g";	     // set the stimulus file name
       if (set_vel) sprintf (stimfil,fname,stimtype,swaveshape,set_ncycles,min_sp,max_sp,velocity,
		       	scontrast,minten,blur_csize,blur_ssize,surr_weight);
       else         sprintf (stimfil,fname,stimtype,swaveshape,set_ncycles,min_sp,max_sp,tfreq,
		       	scontrast,minten,blur_csize,blur_ssize,surr_weight);
       sprintf (stimfil,"%s_%s",stimfil,stim_fnum);
       stim_file (stimfil);                                  // set the stimulus file 
    }
    dstart = stimtime;
    for (s=0; s<nstim; s++) {
      bwidth = min_sp + s * incr_sp;
      if (use_stimfile) {
          stim_blur (blur_csize,center_weight);                 // set the blur array
      }
      stim_backgr(minten,start=simtime);	  	     /* background */
      if (set_ncycles==0) {
          ncycles = celldia / bwidth;
          if (ncycles < 2) ncycles = 2;
      }
      if (set_vel) tf = velocity / bwidth; 
      else         tf = tfreq;
      dur = ncycles / tf;
      if (direction==1) movegrating (0, 0, theta, sphase=0,   bwidth, tf, dir=1,0,1.0, 
		  	scontrast, swaveshape, twaveshape, dstart, dur);
      if (direction==2) movegrating (0, 0, theta, sphase=180, bwidth, tf, dir=2,0,1.0, 
			scontrast, swaveshape, twaveshape, dstart, dur);

     if (surr_weight != 0) { 
       if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
         stim_blur (blur_ssize,surr_weight);                  // set the blur array
         if (direction==1) movegrating (0, 0, theta, sphase=0,   bwidth, tf, dir=1,0,1.0, 
		  	scontrast, swaveshape, twaveshape, dstart+delay, dur);
         if (direction==2) movegrating (0, 0, theta, sphase=180, bwidth, tf, dir=2,0,1.0, 
			scontrast, swaveshape, twaveshape, dstart+delay, dur);
        }
     }
     dstart += dur;
   }

   mvbrt1 = dstart + dur + delay;
   mvbrt2 = dstart + 2*dur + delay;

 }

/* **************************************************************************************** */

 else if (stimtype==17) {	// contrast-reversing stationary grating, 
	 			//   with variable spatial and temporal freq
	  int s;
	  double tf, barw, grtime;
	  double incr;
	  double sphase;
	  double offsetx, delay;
          const char *fname;
          static char stimfil[250];

     if (notinit (blur_csize)) blur_csize = 30;
     if (notinit (blur_ssize)) blur_ssize = 100;
     if (notinit(barwidth_max)) barwidth_max = 400;
     if (notinit(barwidth_incr)) barwidth_incr = 25;

     // delay = surr_delay/2;
     delay = surr_delay;
     if (notinit(ncycles)) ncycles = 2;
     if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_contrev_grating_%d_%d_%d_%g_%g_%g_%g_%g_%g_%g";     // set the stimulus file name
       sprintf (stimfil,fname,stimtype,swaveshape,twaveshape,barwidth_max,tfreq,ncycles,
		       	scontrast,minten,blur_csize,blur_ssize);
       sprintf (stimfil,"%s_%s",stimfil,stim_fnum);
       stim_file (stimfil);                                  // set the stimulus file 
     }

  if (tfreq <= 0) tfreq = 1;
  if (tfreq>0) dur = ncycles / tfreq;
  else         dur = 0.1;

  if (barwidth_incr > 0) { 
     if (notinit(nstim)) nstim = barwidth_max / barwidth_incr;
  }
  for (s=0; s<nstim; s++) {
     grtime = stimtime + s*(dur);
     barw = (s+1) * barwidth_incr;
     if (use_stimfile) {
       stim_blur (blur_csize,center_weight);                 // set the blur array
     }
     stim_backgr(minten,start=simtime);	  	     /* background */
     counterphase_grating_ncycle(stimx, stimy, theta, sphase=0, barw, tfreq, 0, 1.0, 
		                  scontrast, swaveshape, twaveshape, grtime, ncycles);

//     movegrating (stimx, stimy, theta, sphase=0, barw, tfreq, direction=0,0,1.0, 
//		  	scontrast, swaveshape, twaveshape, grtime, dur);
    if (surr_weight != 0) { 
       if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
           stim_blur (blur_ssize,surr_weight);                  // set the blur array
           counterphase_grating_ncycle(stimx, stimy, theta, sphase=0, barw, tfreq, surr_tau, 1.0, 
                                       scontrast, swaveshape, twaveshape, grtime+surr_delay, ncycles);

//           grating_ramp(stimx, stimy, theta, sphase=0, barw, tfreq, direction=0, 1.0, 0, scontrast,
//			 swaveshape,twaveshape,start=grtime+surr_delay, dur, surr_tau, incr= 0.1);
       }
     }
  }
  mvbrt2 = s * (dur + delay);
 }

 /* **************************************************************************************** */

 else if (stimtype==18) {	// drifting linear grating, to & fro: first rightward, then leftward,
	 			//      sine or square, with diff-of-gaussian blur
          int s, dir;
	  int set_vel = 0, set_ncycles = 0;
	  double tf, sphase;
	  double delay, dstart;
          const char *fname;
          static char stimfil[150];
	  double stim_incr;


    if (notinit (blur_csize)) blur_csize = 30;
    if (notinit (blur_ssize)) blur_ssize = 120;
    if (notinit(barwidth)) barwidth = 400;
    if (notinit(velocity)) velocity = 200;
    if (notinit(ncycles)) ncycles = 4;

    delay = surr_delay;

    if (ninfo>=2) {
      ncfprintf(stderr,"# drifting grating, velocity %g, spatial period %g\n", velocity,barwidth);
    }

    if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
       fname = "dsgc_sbac_dbp1_drifting_grating_%d_%d_%d_%g_%g_%g_%g_%g_%g_%g_%g";  // set the stimulus file name
       sprintf (stimfil,fname,stimtype,swaveshape,ncycles,barwidth,velocity,
		       	scontrast,minten,delay,blur_csize,blur_ssize,surr_weight);
       sprintf (stimfil,"%s_%s",stimfil,stim_fnum);
       stim_file (stimfil);                                  // set the stimulus file 
    }
    dstart = stimtime;
      if (use_stimfile) {
          stim_blur (blur_csize,center_weight);                 // set the blur array
      }
      stim_backgr(minten,start=simtime);	  	     /* background */
      tf = velocity / barwidth; 
      dur = ncycles / tf;
      movegrating (0, 0, theta, sphase=0,   barwidth, tf, dir=1,0,1.0, 
		  	scontrast, swaveshape, twaveshape, dstart, dur);
      movegrating (0, 0, theta, sphase=180, barwidth, tf, dir=2,0,1.0, 
			scontrast, swaveshape, twaveshape, dstart+dur, dur);

     if (surr_weight != 0) { 
       if (use_stimfile) {    				     // if moving bar stimulus (and --use_stimfile 1 )
         stim_blur (blur_ssize,surr_weight);                  // set the blur array
         movegrating (0, 0, theta, sphase=0,   barwidth, tf, dir=1,0,1.0, 
		  	scontrast, swaveshape, twaveshape, dstart+delay, dur);
         movegrating (0, 0, theta, sphase=180, barwidth, tf, dir=2,0,1.0, 
			scontrast, swaveshape, twaveshape, dstart+dur+delay, dur);
        }
     }
     dstart += dur;

   mvbrt1 = dstart + dur + delay;
   mvbrt2 = dstart + 2*dur + delay;

  }

 runsort();			/* sort the stimfile if makestim==1 */

  return mvbrt2;
}

/*--------------------------------------------------------*/

void addlabels(void)
{
    int i,j,k,cn;
    int nsynapi, nsynapj;
    int nsynap, nsynapa, nsynapb, nsynapc, nsynapd, nsynape, nsynapg, nsynapds, nsynapsb;
    int nsynapam, nsynapam2, nsynapams, nsynapn;
    int nsynapdx, nsynapsb1, nsynapsb2, nsynapsb5, nsynapsb16, nsynapsb25, nsynapsb53;
    node *npnt;
    photorec *p;

    /* find starburst cells by soma r,theta offset from 0,0 */

 if (sbarr >= 0) {
    c1  = findcella(sbac,0,0);
    c2  = findcella(sbac,sbspac,0);
    c3  = findcella(sbac,sbspac,60);
    c4  = findcella(sbac,sbspac,120);
    c5  = findcella(sbac,sbspac,180);
    c6  = findcella(sbac,sbspac,240);
    c7  = findcella(sbac,sbspac,300);
    c8  = findcella(sbac,sbspac/2,-30);
    c9  = findcella(sbac,sbspac/2,-90);
    c10 = findcella(sbac,sbspac/2,-150);
    c11 = findcella(sbac,sbspac/2,30);
    c12 = findcella(sbac,sbspac/2,90);
    c13 = findcella(sbac,sbspac/2,150);
    c14 = findcella(sbac,sbspac,-30);
    c15 = findcella(sbac,sbspac,30);
    c16 = findcella(sbac,sbspac*1.15,-150);
    c17 = findcella(sbac,sbspac,150);
    c18 = findcella(sbac,sbspac*1.5,-15);
    c19 = findcella(sbac,sbspac*1.5,15);

    c20 = findcella(sbac,sbspac*2,180);
    c21 = findcella(sbac,sbspac*3,180);
    c22 = findcella(sbac,sbspac*4,180);
    c23 = findcella(sbac,sbspac*5,180);
    c24 = findcella(sbac,sbspac*6,180);

    c25 = findcella(sbac,sbspac*1.7,-150);
    c50 = findcella(sbac,sbspac*2,0);
    c53 = findcella(sbac,sbspac*2.3,-150);
  }
  else if (sbarr == -2) {
    c1 = 1;
    c2 = 2;
  }
   if (c1 < 0 && nodepnt==NULL) {
	fprintf (stderr,"# retsim --expt dsgc_sbac_bar:  no cells left\n"); 
        exit(0);
   }
   fprintf (stderr, "# c1 %d c2 %d c5 %d c16 %d c25 %d c53 %d\n",c1,c2,c5,c16,c25,c53);

   /* add light transducer to each dbp1 bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }

   /* add light transducer to each dbp2 bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp2,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp2,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }

   if (notinit(theta))   theta = 0;	/* orientation of bar */
   if (notinit(iroff))   iroff = 1000;	/* r offset for inhib transducers */

   ixoff = iroff * cosdeg(theta);
   iyoff = iroff * -sindeg(theta);

   if (notinit(light_inhib)) light_inhib = 0; 
   if (light_inhib) {

     /* Add light transducer to each small-field amacrine cell */
     /*  offset so that inhibition can be controlled separately */

     /*
     for(npnt=nodepnt; npnt=foreach(npnt,ams,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       p = (photorec*)make_transducer(ndn(ams,cn,soma)); 
       p->xpos=npnt->xloc + ixoff; 
       p->ypos=npnt->yloc + iyoff;
     }
     /* */
   }

  /*  - - - - - - - - - - - - - - - - - - - */

   /* make lists of synaptic inputs to the sbac */

   if (notinit(idist)) idist = 1000;

   nsynapa   = synapse_add  (k=1,dbp1,-1,-1,sbac,-1,3);      /* make list of dbp1 ampa syns to sbac */
   nsynapb   = synapse_add  (k=2,dbp1,-1,-1,am,  -1,4);      /* make list of dbp1 ampa syns to am */
   nsynapc   = synapse_add  (k=3,dbp1,-1,-1,am2, -1,7);      /* make list of dbp1 ampa syns to am2 */
   nsynapg   = synapse_add  (k=4,dbp1,-1,-1,dsgc, 1,2);      /* make list of dbp1 ampa synapses to dsgc */
   nsynapd   = synapse_add  (k=5,dbp2,-1,-1,sbac,-1,3);      /* make list of dbp2 ampa syns to sbac */
   nsynape   = synapse_add  (k=6,sbac,-1,-1,dbp1,-1,1);      /* make list of sbac nACh syns to dbp1 */

   nsynapds  = synapse_add  (k=10,sbac,-1,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */

   nsynapdx  = synapse_add  (k=11,sbac, c1,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=12,sbac, c2,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=13,sbac, c3,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=14,sbac, c8,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=15,sbac, c9,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=16,sbac,c10,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=17,sbac,c11,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=18,sbac,c12,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=19,sbac,c13,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=20,sbac,c14,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=21,sbac,c15,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */
   nsynapdx += synapse_add  (k=22,sbac,c18,-1,dsgc, 1,2);      /* make list of inhib synapses from sbacs to dsgc */

   nsynapsb  = synapse_add  (k=30,sbac,-1,-1,sbac,-1,4);      /* make list of inhib synapses from sbacs to sbacs */
   nsynapsb2 = synapse_add  (k=31,sbac,-1,-1,sbac,c2,4);     /* make list of inhib synapses from sbacs to sbac 2 */
   nsynapsb5 = synapse_add  (k=32,sbac,-1,-1,sbac,c5,4);     /* make list of inhib synapses from sbacs to sbac 5 */

   nsynapsb1 = synapse_add  (k=33,sbac,-1,-1,sbac,c1,4);     /* make list of inhib synapses from sbacs to sbac 1 */
   nsynapsb16 = synapse_add (k=34,sbac,-1,-1,sbac,c2,4);     /* make list of inhib synapses from sbacs to sbac 2 */
   nsynapsb16 = synapse_add (k=35,sbac,-1,-1,sbac,c5,4);     /* make list of inhib synapses from sbacs to sbac 5 */
   nsynapsb25 = synapse_add (k=36,sbac,-1,-1,sbac,c16,4);     /* make list of inhib synapses from sbacs to sbac 16 */
   nsynapsb53 = synapse_add (k=37,sbac,-1,-1,sbac,c25,4);     /* make list of inhib synapses from sbacs to sbac 25 */

   nsynapams  = synapse_add  (k=40,ams,  -1,-1,sbac,-1,2);     /* make list of inhib synapses from ams's to sbacs */
   nsynapams += synapse_add  (k=41,ams,  -1,-1,sbac,-1,3);     /* make list of inhib synapses from ams's to sbacs */
   nsynapam  = synapse_add  (k=42,am,  -1,-1,dsgc,-1,2);      /* make list of inhib synapses from am's to dsgcs */
   nsynapam2 = synapse_add  (k=43,am2, -1,-1,dsgc,-1,2);      /* make list of inhib synapses from am2's to dsgcs */

//    nsynapi += synapse_add (1,amhs,-1,-1,sbac,c1,1);      /* make list of inhib synapses from amhs's */
//    //    nsynapa += synapse_add (2,hbp1,-1,-1,sbac,c1,2);       /* make list of hbp1 ampa syns for recording below */
//    //    nsynapn  = synapse_add (7,dbp1,-1,-1,sbac,c1,7);       /* make list of dbp1 nmda synapses */
//    //    nsynapn += synapse_add (7,hbp1,-1,-1,sbac,c1,5);       /* make list of hbp1 nmda synapses */

//   k = number of synapse lists up to here

   if (!(disp && nsbac>60) && make_sbac_sbac > 0 && syn_savefile != NULL) {
     for (j=1; j<=nsbac; j++) {
       nsynapj = 0;
       for (i=1; i<=nsbac; i++) {
          if (i==j) continue;
          nsynapi = synapse_add (j+k,sbac,i,-1,sbac,j,4,0,idist); /* make list of inhib synapses between sbacs */
          nsynapj += nsynapi;
          nsynap  += nsynapi;
	  ncfprintf (stderr,"# nsynapi %d->%d %d\n",i,j,nsynapi);
       }
       ncfprintf (stderr,"# nsynapj  ->%d %d\n#\n",j,nsynapj);
     }
   }

   nsynapn = 0;
   nsynap = nsynapa + nsynapb + nsynapc + nsynapg + nsynapd + nsynape + nsynapds + nsynapsb + nsynapam + nsynapam2;
   // ncfprintf (stderr,"# nsynap %d nsyn_dbp1->sbac %d nsyn_dbp1->dsgc %d nsyn_sbac_dsgc %d nsyn_sbac->sbac %d\n",
   //              nsynap, nsynapa, nsynapg, nsynapi, nsynapsb);
   ncfprintf (stderr,"# nsyn_dbp1->sbac  %4d\n",nsynapa);
   ncfprintf (stderr,"# nsyn_dbp1->am    %4d\n",nsynapb);
   ncfprintf (stderr,"# nsyn_dbp1->am2   %4d\n",nsynapc);
   ncfprintf (stderr,"# nsyn_dbp1->dsgc  %4d\n",nsynapg);
   ncfprintf (stderr,"# nsyn_dbp2->sbac  %4d\n",nsynapd);
   ncfprintf (stderr,"# nsyn_sbac->dbp1  %4d\n",nsynape);
   ncfprintf (stderr,"# nsyn_sbac->dsgc  %4d\n",nsynapds);
   ncfprintf (stderr,"# nsyn_sbac->sbac  %4d\n",nsynapsb);
   ncfprintf (stderr,"# nsyn_sbac->sbac1 %4d\n",nsynapsb1);
   ncfprintf (stderr,"# nsyn_ams->sbac   %4d\n",nsynapams);
   // ncfprintf (stderr,"# nsyn_am->dsgc   %4d\n",nsynapam);
   // ncfprintf (stderr,"# nsyn_am2->dsgc  %4d\n",nsynapam2);
   ncfprintf (stderr,"# nsyn total      %4d\n",nsynap);
   ncfprintf (stderr,"#\n");
}
		//
/*--------------------------------------------------------*/

#define MAXSTIME 10000
double isbac[MAXSTIME];
int stimtim = 0;
double tailcur = 0;

bool flag = false;
double maxCurrent;
double cond,Gmax;
double isub;

void onplot (void) 

{
    int stime;
   
   stime = int(simtime * 1000);		// msec resolution
   if (stime >= MAXSTIME) stime = MAXSTIME-1;
   if (stime >= 0 && stime < (stimtim + 200)) {			// for stimtype==6
      isbac[stime] =  i(ndn(sbac,1,soma));
      if (stime >= stimtim) tailcur = isbac[stime] - isbac[stime-stimtim];
      else tailcur = 0;
   } else tailcur = 0;
   if (sbarr== -2) {
	if (simtime>=prestimdur && simtime <= prestimdur + soma_clamp_time) {
	  isub = i(ndn(sbac,c1,soma)) - i(ndn(sbac,c2,soma));
	} else isub = 0;
   }
   if (flag) {
      if (sbarr==-2) current = isub;
      else           current  = i(ndn(sbac,  1, soma));
      if (outward<=0) {
         if (maxCurrent > current) { maxCurrent = current;
          	 //  fprintf(stderr, "v: %g, i: %g\n", vpulse, maxCurrent);
         }
         //  fprintf(stderr, "i: %g, maxi: %g\n", vpulse, maxCurrent);    
      } else {
        if (maxCurrent < current) { maxCurrent = current;
	        //  fprintf(stderr, "v: %g, i: %g\n", voltage, maxCurrent);
        }
      }
  }
}

/*--------------------------------------------------------*/

void runexpt(void)

{
#define NUMCBPS 15
    int c, i, ct, cn, cbp_cn, pl, plnum, dbp2x;
    double start, dur, dscale, plsize;
    double ixoff, iyoff, disp_end, psize;
    double rmax, rmin;
    double Imax, Imin;
    double a, r, sr, cmax, cmin;
    double mask, c1_somax, c1_somay;
    double starttime, stim_end;
    node *npnt;
    elem *e;
    photorec *p;
    // chattrib *a;
    int dbp1s[NUMCBPS];
    int dbp2s[NUMCBPS];
    char sbuf[30];

  timinc = 10e-6;
  // ploti = 1e-3;
  // ploti = 1e-4;
  if (!notinit (set_ploti)) ploti = set_ploti;

// printsyns(dbp1,-1, dsgc,-1); 	// print out synapses 
// fprintf (stderr,"lookups %ld finds %ld\n",total_lookups,total_finds);

  cbp_cn    = 68;
  gcdistnod = 582;

  // e = at (ndn(dsgc,1,297), CHAN);
  // a = make_chan (e,NA,2);
  // chset(e);
  // xxy = e->elnum;
  // /* at [dsgc][1][297] chan Na type 2 chset ename xxy; */

  if (notinit(sblur)) sblur = 0;
//  if (notinit(stimdur))   stimdur=0.45;	/* used for non-moving stimuli */
  if (notinit(stimdur))     stimdur=0;	/* used for non-moving stimuli */
  
  Vmax  = 0.04;
  Vmaxg = 0.00;
  Vmin = -0.07;
  if (stimtype==7) {
	  Vmaxg = vstop;
          Vmin = vstart;
  }

  /*  - - - - - - - - - - - - - - - - - - - */

   if (notinit(barwidth))        barwidth = 50;
   if (notinit(barlength))      barlength = 500;
   if (notinit(tfreq))              tfreq = 2;
   if (notinit(direction))      direction = 1;
   if (notinit(swaveshape))    swaveshape = 0;
   if (notinit(twaveshape))    twaveshape = 0;
   if (notinit(makenv))            makenv = 1;

   if (notinit(minten))            minten = -0.05;
   if (notinit(scontrast))      scontrast =  0.01;
   if (notinit(scontrastn))    scontrastn =  0;
   // if (notinit(velocity))     velocity =  2000;
   if (notinit(stimx))              stimx =  0; 
   if (notinit(stimy))              stimy =  0; 
   if (notinit(annrad))            annrad =  300; 

   if (notinit(stimtime))        stimtime =  0.02;
   if (notinit(prestimdur))    prestimdur =  0.05;
   if (notinit(poststimdur))  poststimdur =  0.1;
   if (notinit(tailcurdur))    tailcurdur =  0.05;
   if (notinit(noise_dur))      noise_dur =  1;
   if (notinit(disptime))        disptime =  0.15;
   if (notinit(stimtype))        stimtype =  1;
   if (notinit(stim_theta))    stim_theta =  0;
   if (notinit(rstim_theta))  rstim_theta =  0;
   if (notinit(spotdur))          spotdur =  0.1;
   if (notinit(spotdia))          spotdia =  50;
   // if (notinit(ncycles))          ncycles =  3;
   if (notinit(spotint))          spotint =  0.2;

   if (notinit(vstart))            vstart = -0.06;
   if (notinit(vstop))              vstop =  0.00;
   if (notinit(vstep))              vstep =  0.01;
   if (notinit(sbac_vhold))    sbac_vhold = -0.08; 
   if (notinit(tailvolt))        tailvolt = -0.08;

   if (notinit(istart))            istart = 0e-12;
   if (notinit(istop))              istop = 120e-12;
   if (notinit(istep))              istep = 20e-12;

   stimtim = int(stimtime * 1000);		// integer time for stimulus for calculating tail currents
   if (run_vclamp_sbac > 0) setonplot(onplot);

   // if (run_vclamp!=0 || stimtype==5) {
   //     set_run_on_exit(runonexit);                         // set to erase savefile on ^C
   //     sprintf (savefile,"dsgc_sbac_bar%06d",getpid());       // add pid to file name
   // }

   for (i=0; i<NUMCBPS; i++) {
      dbp1s[i] = 0;
      dbp2s[i] = 0;
   }

   dbp1s[0]  = findmida(dbp1,-150,0);
   dbp1s[1]  = findmida(dbp1,-125,0);
   dbp1s[2]  = findmida(dbp1,-100,0);
   dbp1s[3]  = findmida(dbp1,-75,0);
   dbp1s[4]  = findmida(dbp1,-50,0);
   dbp1s[5]  = findmida(dbp1,-25,0);
   dbp1s[6]  = findmida(dbp1, 0,0);
   dbp1s[7]  = findmida(dbp1,25,0);
   dbp1s[8]  = findmida(dbp1,50,0);
   dbp1s[9]  = findmida(dbp1,75,0);
   dbp1s[10]  = findmida(dbp1,100,0);
   dbp1s[11]  = findmida(dbp1,125,0);
   dbp1s[12]  = findmida(dbp1,150,0);

   if ((ndbp2-ncell_erased[dbp2])==0) dbp2x = dbp1;		// if no dbp2s, make dbp2s[]=dbp1s[]
   else                               dbp2x = dbp2;

   dbp2s[0]  = findmida(dbp2x,-150,0);
   dbp2s[1]  = findmida(dbp2x,-125,0);
   dbp2s[2]  = findmida(dbp2x,-100,0);
   dbp2s[3]  = findmida(dbp2x,-75,0);
   dbp2s[4]  = findmida(dbp2x,-50,0);
   dbp2s[5]  = findmida(dbp2x,-25,0);
   dbp2s[6]  = findmida(dbp2x, 0,0);
   dbp2s[7]  = findmida(dbp2x,25,0);
   dbp2s[8]  = findmida(dbp2x,50,0);
   dbp2s[9]  = findmida(dbp2x,75,0);
   dbp2s[10]  = findmida(dbp2x,100,0);
   dbp2s[11]  = findmida(dbp2x,125,0);
   dbp2s[12]  = findmida(dbp2x,150,0);

  if (strcmp(sbac_file,"morph_sbac3b")==0) {
    if (sbarr==1) { 					// opposing sbacs
      synapse_add (1, sbac, c1, 1215, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 1332, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 1810, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 3589, sbac, c2);		// add synapse for rate plot

      synapse_add (2, sbac, c2,  632, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2,  293, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2,   76, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2,  271, sbac, c1);		// add synapse for rate plot
    }
    else if (sbarr==2) { 				// opposing sbacs
      synapse_add (1, sbac, c1, 1215, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 1589, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 1110, sbac, c2);		// add synapse for rate plot
      //synapse_add (1, sbac, c1, 1934, sbac, c2);		// add synapse for rate plot

      synapse_add (2, sbac, c2, 1718, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2, 2023, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2, 1157, sbac, c1);		// add synapse for rate plot
      //synapse_add (2, sbac, c2, 1166, sbac, c1);		// add synapse for rate plot
    }
  }
  if (strcmp(sbac_file,"morph_sb1")==0) {
  }

   if ((npnt=findnode (sbac,c1,soma))==NULL) {		// find central cell
         ncfprintf (stderr,"retsim: no cells, exiting...\n");
	 exit(0);
   }
   c1_somax = npnt->xloc;				// get soma (x,y) of central cell
   c1_somay = npnt->yloc;


   set_run_on_exit(runonexit);                         // set to erase savefile on ^C
   sprintf (savefile,"dsgc_sbac_bar%06d",getpid());       // add pid to file name

   if (notinit(predur)) predur=0.05;

   simtime = -predur;                                    // must be set ahead of stim_backgr()
   if (!ivplot) setxmin = 0;
   // setxmin = simtime;
  
  if (use_stimfile && (stimtype==3 || stimtype==5)) {         // if sineann stimulus (and --use_stimfile 1 )
           const char *fname;
           char stimfile[100];
      if      (stimtype==3) fname = "dsgc_sbac_sineann_%g_%g";   // set the stimulus file name
      else if (stimtype==5) fname = "dsgc_sbac_grating_%g_%g";
      sprintf (stimfile,fname,scontrast,barwidth);
      stim_file (stimfile);                                   // set the stimulus file
      stim_blur (sblur,0,0,0,1);                              // set the blur and stimulus arrays
  }

  if (stimtype < 12) stim_backgr(minten,start=simtime);				  	 /* background */
  stimdur = move_stim(stimtime, barwidth, stim_theta, scontrast, direction, mask=0);
 
   if (disp & DSTIM) {
	double t;

     if (makestim) return;

     if (notinit(dispsize)) {
	maxsize = max(xarrsiz,yarrsiz);
	dispsize=maxsize*(1+disp_margin);
     }
     display_size(dispsize);
     if (notinit(disp_calib_len)) disp_calib_len = set_int_val(dispsize/7);
     if (disp_calib_len > 0) drcalib(0.97,0.05,disp_calib_len,dispsize,white);

     disp_end = stimtime+stimdur+poststimdur;
     for (starttime=simtime,t=0; t<disp_end; starttime = t, t+= disp_stim_incr) {
        display_stim(starttime, t, dscale=4, disp_stim_max, disp_stim_min);
        // display_stim(starttime, t, dscale=4, -0.025, -0.035); 
        // display_stim(t, dscale=4, -0.035, -0.045); 
        simwait(0.10);
	if (disp & DMOVIE) display_page();
     }
     return;
   }

   if (notinit(sbac_istim)) sbac_istim = 0; 
   if (sbac_istim != 0) cclamp(ndn(sbac,1,soma), sbac_istim, start=prestimdur, dur=1);

   if (notinit(istim)) istim = 0; 
   if (istim != 0) cclamp(ndn(dsgc,1,soma), istim, start=prestimdur, dur=1);

   if (!notinit(make_movie) && make_movie) {
       if (space_time) {  /* movie */ 
         plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmaxg,c=1,"Vdsgc",pl=10,0.35);
         plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmaxg,c=green,"Vtip1",pl=10,0.35);
         plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmaxg,c=red,"Vtip2",pl=10,0.35);
         //plot_na_inact(dsgc, 1, soma, red, "Na[soma]", 12, .35);
       }
   }
   else { 


   if (sbarr==0) {

     plot_v_nod(ct=dbp1,dbp1s[6], 1, -0.08,-0.03,green,   "", pl=50,0.5);  // plot stimulus in central bipolar
     plot_v_nod(ct=dbp1,dbp1s[10],1, -0.08,-0.03,blue,    "", pl=50,0.5); // plot stimulus in distal bipolar
     plot_v_nod(ct=dbp1,dbp1s[2], 1, -0.08,-0.03,red,     "", pl=50,0.5); // plot stimulus in distal bipolar

     plot_v_nod(ct=dbp2x,dbp2s[4],0, -0.08,-0.03,cyan,    "", pl=50,0.5);  // plot stimulus in transient bipolar
     plot_v_nod(ct=dbp2x,dbp2s[4],7, -0.08,-0.03,brown,   "", pl=50,0.5);  // plot stimulus in transient bipolar
     plot_v_nod(ct=dbp2x,dbp2s[4],8, -0.08,-0.03,magenta, "", pl=50,0.5);  // plot stimulus in transient bipolar
     plot_v_nod(ct=dbp2x,dbp2s[4],9, -0.08,-0.03,yellow,  "", pl=50,0.5);  // plot stimulus in transient bipolar

     plot_ca_nod(ct=dbp1,dbp1s[6], 9, 5e-6,brown,         "", pl=49,0.5);  // plot Ca in tonic db5 bipolar
     plot_ca_nod(ct=dbp1,dbp1s[10],9, 5e-6,blue,          "", pl=49,0.5);  // plot Ca in tonic db5 bipolar
     plot_ca_nod(ct=dbp1,dbp1s[2], 9, 5e-6,red,           "", pl=49,0.5);  // plot Ca in tonic db5 bipolar
     plot_ca_nod(ct=dbp2x,dbp2s[4], 9, 5e-6,green,        "", pl=49,0.5);  // plot Ca in transient db4 bipolar
   }

   if (sbarr!=-2) {
    r = 15;
    plot_syncond(findsynloca(dbp1,dsgc,1,-6*r,rstim_theta),cmin=0,cmax=2e-10, blue,    48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,-4*r,rstim_theta),cmin=0,cmax=2e-10, green,   48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,-2*r,rstim_theta),cmin=0,cmax=2e-10, cyan,    48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,0,0),             cmin=0,cmax=2e-10, red,     48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,2*r,rstim_theta), cmin=0,cmax=2e-10, magenta, 48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,4*r,rstim_theta), cmin=0,cmax=2e-10, brown,   48,"",0.5);
    plot_syncond(findsynloca(dbp1,dsgc,1,6*r,rstim_theta), cmin=0,cmax=2e-10, white,   48,"",0.5);

    if (g_dbp1_sbac < 1e-14) cmax = 1e-14;
    else                     cmax = 0.5e-10;

    plot_syncond(findsynloc(dbp1,dbp1s[6],0,0),  cmin=0,cmax, blue,   46,"",0.5);
    plot_syncond(findsynloc(dbp1,dbp1s[10],0,0), cmin=0,cmax, blue,   46,"",0.5);
    plot_syncond(findsynloc(dbp1,dbp1s[2],0,0),  cmin=0,cmax, red,    46,"",0.5);

    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1,-6*r,rstim_theta),cmin=0,cmax, blue,    46,"",0.5);
    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1,-4*r,rstim_theta),cmin=0,cmax, green,   46,"",0.5);
    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1,-2*r,rstim_theta),cmin=0,cmax, cyan,    46,"",0.5);
    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1, 0*r,rstim_theta),cmin=0,cmax, red,     46,"",0.5);
    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1, 2*r,rstim_theta), cmin=0,cmax, magenta,46,"",0.5);
    plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1, 4*r,rstim_theta), cmin=0,cmax, brown,  46,"",0.5);
    if (!notinit(t2bipolar) && (t2bipolar>0)) 
      plot_syncond(findsynlocat(dbp2,dbp2,sbac,-1, 6*r,rstim_theta), cmin=0,cmax, white,  46,"",0.5);
    else
      plot_syncond(findsynlocat(dbp1,dbp2,sbac,-1, 6*r,rstim_theta), cmin=0,cmax, white,  46,"",0.5);

    // plot_synrate(findsynloc(dbp1,dbp1s[10],0,0),rmin=0,rmax=200, blue,   45,"",0.5);
    // plot_synrate(findsynloc(dbp1,dbp1s[2],0,0), rmin=0,rmax=200, red,    45,"",0.5);

    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1,-6*r,rstim_theta),rmin=0,rmax=200, blue,   45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1,-4*r,rstim_theta),rmin=0,rmax=200, green,  45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1,-2*r,rstim_theta),rmin=0,rmax=200, cyan,   45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1, 0*r,rstim_theta),rmin=0,rmax=200, red,    45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1, 2*r,rstim_theta),rmin=0,rmax=200, magenta,45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1, 4*r,rstim_theta),rmin=0,rmax=200, brown,  45,"",0.5);
    plot_synrate(findsynlocat(dbp1,dbp2,sbac,-1, 6*r,rstim_theta),rmin=0,rmax=200, white,  45,"",0.5);

    a = 30;
    if (notinit(run_vclamp_sbac)) run_vclamp_sbac = 0; 
    if (run_vclamp_sbac>0) Vmin_c1 = sbac_vhold;
    else                   Vmin_c1 = Vmin;

    // plot specific nodes of c1
    
       if (!notinit (plotnod1)) 
         plot_v_nod(ct=sbac,c1,plotnod1, Vmin_c1,Vmaxg,blue,   "", pl=43,0.8);
       if (!notinit (plotnod2)) 
         plot_v_nod(ct=sbac,c1,plotnod2, Vmin_c1,Vmaxg,gray,   "", pl=43,0.8);
       if (!notinit (plotnod3)) 
         plot_v_nod(ct=sbac,c1,plotnod3, Vmin_c1,Vmaxg,cyan,   "", pl=43,0.8);
       if (!notinit (plotnod4)) 
         plot_v_nod(ct=sbac,c1,plotnod4, Vmin_c1,Vmaxg,ltblue, "", pl=43,0.8);

    // look along axis of left-right pointing dendrites
    
       if (notinit(set_sr)) set_sr = 15;
       sr = set_sr; 
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,    0,0),  Vmin_c1,Vmaxg,blue,   "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,  1*sr,0), Vmin_c1,Vmaxg,gray,   "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,  2*sr,0), Vmin_c1,Vmaxg,cyan,   "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,  3*sr,0), Vmin_c1,Vmaxg,ltblue, "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,  4*sr,0), Vmin_c1,Vmaxg,yellow, "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1,  5*sr,0), Vmin_c1,Vmaxg,magenta,"", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, -1*sr,0), Vmin_c1,Vmaxg,brown,  "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, -2*sr,0), Vmin_c1,Vmaxg,white,  "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, -3*sr,0), Vmin_c1,Vmaxg,ltgreen,"", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, -4*sr,0), Vmin_c1,Vmaxg,green,  "", pl=40,0.8);
       plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, -5*sr,0), Vmin_c1,Vmaxg,red,    "", pl=40,0.8);

       if (mglur2 > 0) {
          plot_mglur2(findsynloca(ct=dbp1,40,0),0,1,gray, pl=39, "", 0.3);
          plot_mglur2(findsynloca(ct=dbp1,70,0),0,1,cyan, pl=39, "", 0.3);
          plot_mglur2(findsynloca(ct=dbp1,90,0),0,1,blue, pl=39, "", 0.3);
       }

       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,    0,0),  30e-6,blue,   "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,  1*sr,0), 30e-6,gray,   "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,  2*sr,0), 30e-6,cyan,   "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,  3*sr,0), 30e-6,ltblue, "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,  4*sr,0), 30e-6,yellow, "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1,  5*sr,0), 30e-6,magenta,"",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1, -1*sr,0), 30e-6,brown,  "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1, -2*sr,0), 30e-6,white,  "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1, -3*sr,0), 30e-6,ltgreen,"",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1, -4*sr,0), 30e-6,green,  "",pl=38,0.8);
       plot_ca_nod(ct=sbac,c1,findnodlocra(sbac,c1, -5*sr,0), 30e-6,red,    "",pl=38,0.8);

    //plot_ca_nod(dsgc,1,soma,cyan,1e-6,"Ca_soma", -1, -1);

    // look at distal dendrites at increments of "a" (30) deg
  
   /*  
   if (stimtype !=8) { 
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,0*a), Vmin_c1,Vmaxg,blue,  "",pl=37,0.8);
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,2*a), Vmin_c1,Vmaxg,green, "",pl=37,0.8);
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,4*a), Vmin_c1,Vmaxg,cyan,  "",pl=37,0.8);
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,6*a), Vmin_c1,Vmaxg,red,   "",pl=37,0.8);
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,8*a), Vmin_c1,Vmaxg,magenta,"",pl=37,0.8);
    plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,10*a), Vmin_c1,Vmaxg,brown, "",pl=37,0.8);
    // plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,9*a), Vmin_c1,Vmaxg,gray,  "",pl=37,0.8);
    // plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,12*a), Vmin_c1,Vmaxg,white, "",pl=37,0.8);
   }
   /* */

      if (sbarr==9 || sbarr==11) {   // look along axis of left-right pointing dendrites
			//
       plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2,    0,0),   Vmin,Vmaxg,blue,   "",pl=37,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 3*sr,0),   Vmin,Vmaxg,green,  "",pl=37,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 5*sr,0),   Vmin,Vmaxg,magenta,"",pl=37,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2,-3*sr,0),   Vmin,Vmaxg,ltgreen,"",pl=37,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2,-5*sr,0),   Vmin,Vmaxg,red,    "",pl=37,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5,    0,0),   Vmin,Vmaxg,blue,   "",pl=36,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 3*sr,0),   Vmin,Vmaxg,green,  "",pl=36,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 5*sr,0),   Vmin,Vmaxg,magenta,"",pl=36,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5,-3*sr,0),   Vmin,Vmaxg,ltgreen,"",pl=36,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5,-5*sr,0),   Vmin,Vmaxg,red,    "",pl=36,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c20,findnodlocra(sbac,c20,    0,0), Vmin,Vmaxg,blue,   "",pl=35,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c20,findnodlocra(sbac,c20, 3*sr,0), Vmin,Vmaxg,green,  "",pl=35,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c20,findnodlocra(sbac,c20, 5*sr,0), Vmin,Vmaxg,magenta,"",pl=35,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c20,findnodlocra(sbac,c20,-3*sr,0), Vmin,Vmaxg,ltgreen,"",pl=35,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c20,findnodlocra(sbac,c20,-5*sr,0), Vmin,Vmaxg,red,    "",pl=35,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c21,findnodlocra(sbac,c21,    0,0), Vmin,Vmaxg,blue,   "",pl=34,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c21,findnodlocra(sbac,c21, 3*sr,0), Vmin,Vmaxg,green,  "",pl=34,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c21,findnodlocra(sbac,c21, 5*sr,0), Vmin,Vmaxg,magenta,"",pl=34,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c21,findnodlocra(sbac,c21,-3*sr,0), Vmin,Vmaxg,ltgreen,"",pl=34,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c21,findnodlocra(sbac,c21,-5*sr,0), Vmin,Vmaxg,red,    "",pl=34,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c22,findnodlocra(sbac,c22,    0,0), Vmin,Vmaxg,blue,   "",pl=33,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c22,findnodlocra(sbac,c22, 3*sr,0), Vmin,Vmaxg,green,  "",pl=33,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c22,findnodlocra(sbac,c22, 5*sr,0), Vmin,Vmaxg,magenta,"",pl=33,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c22,findnodlocra(sbac,c22,-3*sr,0), Vmin,Vmaxg,ltgreen,"",pl=33,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c22,findnodlocra(sbac,c22,-5*sr,0), Vmin,Vmaxg,red,    "",pl=33,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c23,findnodlocra(sbac,c23,    0,0), Vmin,Vmaxg,blue,   "",pl=32,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c23,findnodlocra(sbac,c23, 3*sr,0), Vmin,Vmaxg,green,  "",pl=32,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c23,findnodlocra(sbac,c23, 5*sr,0), Vmin,Vmaxg,magenta,"",pl=32,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c23,findnodlocra(sbac,c23,-3*sr,0), Vmin,Vmaxg,ltgreen,"",pl=32,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c23,findnodlocra(sbac,c23,-5*sr,0), Vmin,Vmaxg,red,    "",pl=32,0.8);/* V at dist nod */

       plot_v_nod(ct=sbac,c24,findnodlocra(sbac,c24,    0,0), Vmin,Vmaxg,blue,   "",pl=31,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c24,findnodlocra(sbac,c24, 3*sr,0), Vmin,Vmaxg,green,  "",pl=31,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c24,findnodlocra(sbac,c24, 5*sr,0), Vmin,Vmaxg,magenta,"",pl=31,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c24,findnodlocra(sbac,c24,-3*sr,0), Vmin,Vmaxg,ltgreen,"",pl=31,0.8);/* V at dist nod */
       plot_v_nod(ct=sbac,c24,findnodlocra(sbac,c24,-5*sr,0), Vmin,Vmaxg,red,    "",pl=31,0.8);/* V at dist nod */
      }
     else if (nsbac > 1) {
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,0*a), Vmin,Vmaxg,blue,  "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,3*a), Vmin,Vmaxg,green, "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,4*a), Vmin,Vmaxg,cyan,  "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,5*a), Vmin,Vmaxg,red,   "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,6*a), Vmin,Vmaxg,magenta,"",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,7*a), Vmin,Vmaxg,brown, "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,8*a), Vmin,Vmaxg,white, "",pl=37,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,9*a), Vmin,Vmaxg,gray,  "",pl=37,0.8);/* V at dist nod */

      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,0*a), Vmin,Vmaxg,blue,  "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,3*a), Vmin,Vmaxg,green, "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,4*a), Vmin,Vmaxg,cyan,  "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,5*a), Vmin,Vmaxg,red,   "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,6*a), Vmin,Vmaxg,magenta,"",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,7*a), Vmin,Vmaxg,brown, "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,8*a), Vmin,Vmaxg,white, "",pl=36,0.8);/* V at dist nod */
      plot_v_nod(ct=sbac,c5,findnodlocra(sbac,c5, 110,9*a), Vmin,Vmaxg,gray,  "",pl=36,0.8);/* V at dist nod */
     }
    }
    if (ivplot) {
	graph_x(vstop, vstart);
	graph_y(1e-10, -iscal);
	// graph_y(gscal, 0);
	graph_init();
    }


    // plot_v_nod(ct=sbac,c1,findnodlocra(sbac,c1, 110,0*a), Vmin,Vmaxg,blue,  "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c2,findnodlocra(sbac,c2, 110,3*a), Vmin,Vmaxg,green, "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c3,findnodlocra(sbac,c3, 110,4*a), Vmin,Vmaxg,cyan,  "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c8,findnodlocra(sbac,c8, 110,5*a), Vmin,Vmaxg,red,   "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c9,findnodlocra(sbac,c9, 110,6*a), Vmin,Vmaxg,magenta,"",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c10,findnodlocra(sbac,c10, 110,7*a), Vmin,Vmaxg,brown, "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c11,findnodlocra(sbac,c11, 110,8*a), Vmin,Vmaxg,white, "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c12,findnodlocra(sbac,c12, 110,9*a), Vmin,Vmaxg,gray,  "",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c13,findnodlocra(sbac,c13, 110,9*a), Vmin,Vmaxg,ltblue,"",pl=23,0.5);/* V at dist nod */
    // plot_v_nod(ct=sbac,c14,findnodlocra(sbac,c14, 110,9*a), Vmin,Vmaxg,ltmag,"",pl=23,0.5);/* V at dist nod */

    if (n_sbac >= 23) {
      plot_v_nod(ct=sbac,c1,soma,Vmin,Vmaxg,blue,    "",pl=22,0.5);/* V at soma */
      plot_v_nod(ct=sbac,c10,soma,Vmin,Vmaxg,magenta,"",pl=22,0.5);/* V at soma */
      plot_v_nod(ct=sbac,c16,soma,Vmin,Vmaxg,green,  "",pl=22,0.5);/* V at soma */
      if (n_sbac >= 45)
          plot_v_nod(ct=sbac,c25,soma,Vmin,Vmaxg,ltblue, "",pl=22,0.5);/* V at soma */
      if (n_sbac >= 57) 
          plot_v_nod(ct=sbac,c53,soma,Vmin,Vmaxg,red,  "",pl=22,0.5);/* V at soma */
    }
    // else if (nsbac > 1 && sbarr != -2) {
    //  plot_v_nod(ct=sbac,c1,soma,Vmin,Vmaxg,blue,   "",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c2,soma,Vmin,Vmaxg,green,  "",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c3,soma,Vmin,Vmaxg,cyan,   "",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c4,soma,Vmin,Vmaxg,red,    "",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c5,soma,Vmin,Vmaxg,magenta,"",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c6,soma,Vmin,Vmaxg,yellow, "",pl=22,0.5);/* V at soma */
    //  plot_v_nod(ct=sbac,c7,soma,Vmin,Vmaxg,gray,   "",pl=22,0.5);/* V at soma */
    // }

    // plot_func(gsyn_tot,34,rmax=1e-8,rmin=0); plot_param("Gsbac_sbac2",blue,  21,0.5);
    // plot_func(gsyn_tot,35,rmax=1e-8,rmin=0); plot_param("Gsbac_sbac5",red,   21,0.5);

    if (sbarr==100) {
      sprintf (sbuf, "Gsbac_sbac%d",c1);
      plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,blue,   21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c2);
      plot_func(gsyn_tot,34,rmax=1e-9,rmin=0); plot_param(sbuf,magenta,21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c5);
      plot_func(gsyn_tot,35,rmax=1e-9,rmin=0); plot_param(sbuf,green, 21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c16);
      plot_func(gsyn_tot,36,rmax=1e-9,rmin=0); plot_param(sbuf,ltblue,21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c25);
      plot_func(gsyn_tot,37,rmax=1e-9,rmin=0); plot_param(sbuf,red,   21,0.5);

    } else if (sbarr == 107) {
      sprintf (sbuf, "Gsbac_sbac%d",c1);
      plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,blue,   21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c2);
      plot_func(gsyn_tot,31,rmax=1e-9,rmin=0); plot_param(sbuf,green, 21,0.5);
      // sprintf (sbuf, "Gsbac_sbac%d",c3);
      // plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,cyan,   21,0.5);
      // sprintf (sbuf, "Gsbac_sbac%d",c4);
      // plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,red,   21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c5);
      plot_func(gsyn_tot,32,rmax=1e-9,rmin=0); plot_param(sbuf,magenta,21,0.5);
      // sprintf (sbuf, "Gsbac_sbac%d",c6);
      // plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,gray,   21,0.5);

    } else if (nsbac > 1 && make_sbac_sbac==1) {
      sprintf (sbuf, "Gsbac_sbac%d",c1);
      plot_func(gsyn_tot,33,rmax=1e-9,rmin=0); plot_param(sbuf,blue,   21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c2);
      plot_func(gsyn_tot,34,rmax=1e-9,rmin=0); plot_param(sbuf,magenta,21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c5);
      plot_func(gsyn_tot,35,rmax=1e-9,rmin=0); plot_param(sbuf,green, 21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c16);
      plot_func(gsyn_tot,36,rmax=1e-9,rmin=0); plot_param(sbuf,ltblue,21,0.5);
      sprintf (sbuf, "Gsbac_sbac%d",c25);
      plot_func(gsyn_tot,37,rmax=1e-9,rmin=0); plot_param(sbuf,red,   21,0.5);
    }

    if (ndsgc > 0) {
       sprintf (sbuf, "Ca_sbac%d_dsgc",c1);
       // plot_func(casyn_avg,11,rmax=10e-6,rmin=0); plot_param("Ca_sbac1_dsgc", magenta, 10,0.5);
       plot_func(casyn_avg,11,rmax=10e-6,rmin=0); plot_param(sbuf, magenta, 10,0.5);

       sprintf (sbuf, "R_sbac%d_dsgc",c1);
       // plot_func(rsyn_avg,11,rmax=400,rmin=0);  plot_param("R_sbac1_dsgc", brown, 9,0.5);
       plot_func(rsyn_avg,11,rmax=400,rmin=0);  plot_param(sbuf, brown, 9,0.5);

      plot_func(gsyn_tot,10,rmax=2e-8,rmin=0); plot_param("Gsbac_dsgc",  magenta, 8,0.5);
      if (nsbac > 1 && ndsgc > 0) {
         plot_func(gsyn_tot,11,rmax=2e-8,rmin=0); plot_param("Gsbac1_dsgc", blue,    8,0.5);
         plot_func(gsyn_tot,12,rmax=2e-8,rmin=0); plot_param("Gsbac2_dsgc", green,   8,0.5);
         plot_func(gsyn_tot,13,rmax=2e-8,rmin=0); plot_param("Gsbac3_dsgc", cyan,    8,0.5);
         plot_func(gsyn_tot,14,rmax=2e-8,rmin=0); plot_param("Gsbac8_dsgc", red,     8,0.5);
         plot_func(gsyn_tot,15,rmax=2e-8,rmin=0); plot_param("Gsbac9_dsgc", brown,   8,0.5);
         plot_func(gsyn_tot,16,rmax=2e-8,rmin=0); plot_param("Gsbac10_dsgc",white,   8,0.5);
         plot_func(gsyn_tot,17,rmax=2e-8,rmin=0); plot_param("Gsbac11_dsgc",gray,    8,0.5);
         plot_func(gsyn_tot,18,rmax=2e-8,rmin=0); plot_param("Gsbac12_dsgc",ltblue,  8,0.5);
         plot_func(gsyn_tot,19,rmax=2e-8,rmin=0); plot_param("Gsbac13_dsgc",ltgreen, 8,0.5);
         plot_func(gsyn_tot,20,rmax=2e-8,rmin=0); plot_param("Gsbac14_dsgc",ltcyan,  8,0.5);
         plot_func(gsyn_tot,21,rmax=2e-8,rmin=0); plot_param("Gsbac15_dsgc",yellow,  8,0.5);
         plot_func(gsyn_tot,22,rmax=2e-8,rmin=0); plot_param("Gsbac18_dsgc",ltmag,   8,0.5);
      }
      plot_func(gsyn_tot,4, rmax=2e-8,rmin=0); plot_param("Gdbp1_dsgc",blue,    7,0.5);
    }

    //plot_ca_nod(dsgc,1,soma,cyan,1e-6,"Ca_soma", -1, -1);
    //plot_ca_nod(dsgc,1,gcdistnod,1e-6,yellow,"", -1, -1);
    //plot_v_nod(ct=cbp,cbp_cn,axtrm,   Vmin,Vmax,red,"", -1, 0.35);
    //plot_v_nod(ct=cbp,41,axtrm,      Vmin,Vmax,blue,"", -1, 0.35);
    //plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmax,c=green,"Vtip1",pl=10,0.35);
    //plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmax,c=red,"Vtip2",pl=10,0.35);
    //plot_v_nod(ct=dsgc,cn=1,422,      Vmin,Vmax,red,"", 10, 0.35);
    //plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,c=red,"Vtip2", pl=10, .35);
    //plot_v_nod(ct=dsgc,cn=1,1336,      Vmin,Vmaxg,c=green,"Vtip1", pl=10, .35);
    //plot_v_nod(ct=dsgc,cn=1,1464,     Vmin,Vmaxg,c=blue,"", -1, -1);
    //plot_v_nod(ct=dsgc,cn=1,2328,     Vmin,Vmaxg,c=magenta,"",pl=10,1);
    //plot_synrate_out(cbp,cbp_cn,0,500,green);
    //plot_synrate_out(cbp,241,0,500,blue);
    //plot_currents(ct=dsgc,plgain=200e-12);

     // if (notinit(node_dist)) node_dist = 281;
     // plot_chan_current(ct=dsgc, cn=1, node_dist, NMDA, 1, 3e-12, -3e-12);
     // sprintf (sbuf,"Inmda  %d",node_dist);
     // plot_param (sbuf, blue, plnum=5, plsize=0.5);
     // plot_chan_current(ct=dsgc, cn=1, node_dist, AMPA, 5, 3e-12, -3e-12);
     // sprintf (sbuf,"Iampa  %d",node_dist);
     // plot_param (sbuf, gray, plnum=5, plsize=0.5);


    // /* display amacrine cell responses */
    //
    // plot_v_nod(ct=am, cn=3,8,  Vmin,Vmax,green,"", 6, 0.3);
    // plot_v_nod(ct=am2,cn=3,8,  Vmin,Vmax,blue, "", 6, 0.3);
    // if (g_am_dsgc > 1e-12) {
    //   r =20;
    //   plot_syncond(findsynloca(am, 3,dsgc,1,1*r,rstim_theta), cmin=0,cmax=4e-10, green, 5,"",0.3);
    //   plot_syncond(findsynloca(am2,3,dsgc,1,1*r,rstim_theta), cmin=0,cmax=4e-10, blue,  5,"",0.3);
    // }

    if (notinit(run_vclamp))           run_vclamp = 0;
    if (notinit(run_vclamp_sbac))      run_vclamp_sbac = 0; 
    if (notinit(sbac_vpulse))          sbac_vpulse = -0.00; 
    if (notinit(sbac_vpulse_dur))      sbac_vpulse_dur = 0.1; 
    if (notinit(sbac_vpulse_time))     sbac_vpulse_time = 0.02; 
    if (notinit(soma_clamp_time))      soma_clamp_time = 0.01; 

    if (run_vclamp > 0 && ndsgc > 0) {
      if (notinit(vhold)) vhold = -0.07;
      vclamp (ndn(ct=dsgc,cn=1,soma),   vhold, simtime,  10);
      plot_i_nod(ct=dsgc,cn=1,soma,     Imin=-2000e-12,Imax=1000e-12,magenta,"", 1, 0.5);
      plot_v_nod(ct=dsgc,cn=1,soma,   Vmin,Vmaxg,red,"", 2, 0.5);
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,0,DEND),   Vmin,Vmaxg,blue,"",   2, 0.5);
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,-30,DEND), Vmin,Vmaxg,green,"",  2, 0.5);
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,30,DEND), Vmin,Vmaxg,magenta,"",2, 0.5);/* V at dist nod */
    } else {
      plot_v_nod(ct=dsgc,cn=1,soma,   Vmin,Vmaxg,red,"", 1, 0.5);
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,0,DEND), Vmin,Vmaxg,blue,    "",1,0.5);/* V at dist nod */
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,-30,DEND), Vmin,Vmaxg,green, "",1,0.5);/* V at dist nod */
      plot_v_nod(ct=dsgc,cn=1,findnodlocra(dsgc,1, 110,30,DEND), Vmin,Vmaxg,magenta,"",1,0.5);/* V at dist nod */
      // plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.3); 
    } 
   }

  if (stimtype==7) {
     stimdur = sbac_vpulse_dur;
     vclamp                 (ndn(sbac,c1, soma), sbac_vhold, simtime,  predur);
     if (sbarr==-2) vclamp  (ndn(sbac,c2, soma), sbac_vhold, simtime,  predur);
     if (!ivplot) {
	plot_i_nod   (ct=sbac,c1,soma,   Imax=2e-9,Imin=-2e-9,red,"", 5, 0.8);
	if (sbarr== -2) {
	  plot_i_nod (ct=sbac,c2,soma,   Imax=2e-9,Imin=-2e-9,blue,"", 5, 0.5);
	  plot_var (&isub,               1, Imax=0e-12,Imin= -200e-12);
	  plot_param ("Isbac sub",magenta, 4, 0.5);
	}
        endexp=prestimdur+sbac_vpulse_dur+tailcurdur+poststimdur;
     }
  }
  else endexp=stimtime+stimdur+poststimdur;

  step (1e-4);		// short step before run to create model with channels
  if (mglur2 > 0) set2ndmsg(CA,-1,dbp1,sbac,AMPA,sb_mglur_maxdist,mglur2,50,0.15);    // set glutamate to modulate Ca chans

  synaptau = 0.001;	// set synapses to run with shorter time constant

  step(predur);

  synaptau = 1.0;	// set synapses to run with normal time constant

  /* set movie plot routine */
  // setonplot(onplot_movie);

  /* run experiment */

   if (run_vclamp_sbac > 0) {
       if (stimtype==6) { 
	       vclamp (ndn(ct=sbac,c1,soma), sbac_vhold, simtime, 10);
	       vclamp (ndn(ct=sbac,c1,soma), sbac_vpulse, sbac_vpulse_time, sbac_vpulse_dur);
	       vclamp (ndn(ct=sbac,c1,soma), sbac_vhold, sbac_vpulse_time+sbac_vpulse_dur, sbac_vpulse_dur);
	       vclamp (ndn(ct=sbac,c1,soma), sbac_vpulse, stimtime + sbac_vpulse_time, sbac_vpulse_dur);
	       vclamp (ndn(ct=sbac,c1,soma), sbac_vhold, stimtime+sbac_vpulse_time+sbac_vpulse_dur, 10);
	       plot_i_nod (ct=sbac,c1,soma,   Imin=-100e-12,Imax=100e-12,magenta,"", 3, 0.5);
	       plot_var (&tailcur,1,Imin,Imax);
	       plot_param ("Itail-Ctrl", green, 3, 0.5);

	} else if (stimtype !=7 && stimtype !=8) {
		vclamp (ndn(ct=sbac,c1,soma), sbac_vhold, stimtime,  stimdur);
        }
    }

    if (stimtype==7) {			// vclamp
             double dst, sign, vpulse;

	 dst = 0.0001;         // small time to allow voltage clamp to end

	 // if (ivplot) graph_pen(i+1,i+1,i+1,i+1,i+1);
	 if (!ivplot) flag = true;
	 if (vstart < vstop) sign = 1;
	 else                sign = -1;

	 savemodel (savefile);
	 ct = sbac;
	 c1 = 1;

	 for (i=0,vpulse=vstart; (vpulse*sign)<=(vstop*sign+1e-6); i++,vpulse += vstep) {

	        if (!ivplot) graph_pen(red,blue,i+1,i+1);
		simtime = 0;
                stimdur = move_stim(stimtime, barwidth, stim_theta, scontrast, direction, mask=0);

		vclamp              (ndn(ct,c1, soma), sbac_vhold, simtime,  prestimdur);
		if (sbarr==-2)
		  vclamp            (ndn(ct,c2, soma), sbac_vhold, simtime,  prestimdur);
		step (prestimdur);

		if (outward) maxCurrent =  0;
		else         maxCurrent =  1000;
		Gmax = 0;
		vclamp              (ndn(ct,c1, soma), vpulse, simtime,  sbac_vpulse_dur);
		if (sbarr==-2)
		  vclamp            (ndn(ct,c2, soma), vpulse, simtime,  sbac_vpulse_dur);
		step (dst);
		if (ivplot) flag = true;
		step (sbac_vpulse_dur-2*dst);
		//if (ivplot) graph(vpulse, maxCurrent, Gmax);
		if (ivplot) graph(vpulse, maxCurrent);
		if (ivplot) flag = false;
		step (dst);
		
		vclamp              (ndn(ct,c1, soma), tailvolt, simtime,  tailcurdur);
		if (sbarr==-2)
		  vclamp            (ndn(ct,c2, soma), tailvolt, simtime,  tailcurdur);
		step (tailcurdur);
		
		vclamp              (ndn(ct,c1, soma), sbac_vhold, simtime,  poststimdur);
		if (sbarr==-2)
		  vclamp            (ndn(ct,c2, soma), sbac_vhold, simtime,  poststimdur);
		step (poststimdur);
		
		restoremodel (savefile);
    	   }
   	   unlink (savefile);

   } else if (stimtype==8) {		// cclamp
             double sign, ipulse;

	 if (istart < istop) sign = 1;
	 else                sign = -1;

	 savemodel (savefile);
	 ct = sbac;

	 for (i=0,ipulse=istart; (ipulse*sign)<=(istop*sign+1e-14); i++,ipulse += istep) {

		simtime = 0;
                stimdur = move_stim(stimtime, barwidth, stim_theta, scontrast, direction, mask=0);

		step (prestimdur);

		if (outward) maxCurrent =  0;
		else         maxCurrent =  1000;
		Gmax = 0;
		cclamp   (ndn(ct,c1, soma), ipulse, simtime,  sbac_vpulse_dur);

		step (sbac_vpulse_dur);
		
		// cclamp   (ndn(ct,c1, soma), 0, simtime,  tailcurdur);
		step (tailcurdur);
		
		// cclamp   (ndn(ct,c1, soma), 0, simtime,  poststimdur);
		step (poststimdur);
		
		restoremodel (savefile);
    	   }
   	   unlink (savefile);

   } else {

    if (stimtype != 6) savemodel (savefile);
    step(stimtime+mvbrt1+poststimdur);
 
    if (stimtype != 6) {
      restoremodel (savefile);
      step(stimtime+mvbrt1+poststimdur);
      unlink (savefile);
    }

  }
  savefile[0] = 0;
}
