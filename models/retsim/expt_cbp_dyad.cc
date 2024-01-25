/* Experiment cbp_dyad */
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
int ct2 = -1;
int cn2 = 1;
int ivplot;
int elnode;
int cone_type;
int bp_morph;
int amarr;
int no_tip_labels;
int no_rpnt_labels;
int disp_c2;
int nspots;
int s;
int stimtype;
int draw_synapse;

int stimloc;
int axon_base;
int vcloc;
int itransducer;
int npixels;
int nsines;

#define R8 R_8

int tip1 = R8*100+1;
int tip2 = R8*100+2;
int tip3 = R8*100+3;
int tip4 = R8*100+4;
int tip5 = R8*100+5;
int tip6 = R8*100+6;
int tip7 = R8*100+7;
int tip8 = R8*100+8;
int tip9 = R8*100+9;
int tip10 = R8*100+10;
int tip11 = R8*100+11;
int tip12 = R8*100+12;
int tip13 = R8*100+13;
int tip14 = R8*100+14;
int tip15 = R8*100+15;
int tip16 = R8*100+16;
int tip17 = R8*100+17;
int tip18 = R8*100+18;
int tip19 = R8*100+19;
int tip20 = R8*100+20;
int tip21 = R8*100+21;
int tip22 = R8*100+22;
int tip23 = R8*100+23;
int tip24 = R8*100+24;
int tip25 = R8*100+25;
int tip26 = R8*100+26;
int tip27 = R8*100+27;
int tip28 = R8*100+28;
int tip29 = R8*100+29;
int tip30 = R8*100+30;
int recpnt1;
int recpnt2;
int recpnt3;
int recpnt4;
int recpnt5;
int recpnt6;
int recpnt7;
int recpnt8;
int recpnt9;
int presyn1;
int presyn2;
int presyn3;
int ampar;
int ampar2;
int gabac;
int dbp1_am;
int dbp2_am;
int dbp1_morph;
double gabactrc;

const char *cbptype;
const char *cbptype2;

double kdr_cond;
double dbp1_gtau;
double dbp2_gtau;
double dbp1_g4tau;
double dbp2_g4tau;
double am_atau;
double am2_atau;

double am_sdur;
double am_sdur2;
double am2_sdur;
double am2_sdur2;
double am_sfall;
double am_sfall2;
double am2_sfall;
double am2_sfall2;

double g_am_dbp1;
double g_am_dbp2;
double g_am2_dbp1;
double g_am2_dbp2;
double g_dbp1_am;
double g_dbp1_am2;
double g_dbp2_am;
double g_dbp2_am2;
double g_dbp1_gca;
double dbp1_mp;
double dbp1_mr;
double dbp1pca4;
double dbp1pca7;
double dbp2pca5;
double am2dbp1ca;
double dbp1_dyad;

double drmab;
double driab;
double dria;
double cbp_rm;
double dvreva;
double dcmab;
double axdia;
double axdiap;
double ddia;
double naax;
double naab;
double nahd;
double ksoma;
double kr6;
double kr7;
double kasoma;
double kar6;
double kar7;
double dbp1_k1o;
double khsoma;
double casoma;
double cadist;
double cadist2;
double catdist;
double catsoma;
double capsoma;
double cone_soma_z;
double dbp1_soma_z;
double dbp2_soma_z;
double spacing;
double dtreedia;
double dtreedia2;
double c2rot;                   // z rot for 2nd cell of stereo pair (sets disparity for side view)
double c2yrot;                  // y rot for 2nd cell of stereo pair (sets disparity for bottom view)

double amvst;                   // for dens_cbp_inh_am.n
double amvrev;                  // for dens_cbp_inh_am.n
double amrm;
double amca;
double amsca;
double am2ca;
double amna;
double amk;

double varicos_dia;

double predur;
double prestimdur;
double spotdur;
double spotduri;
double stimdur;
double tailcurdur;
double poststimdur;
double stimtime;
double tfreq;
double speriod;
double spotdia;
double mask_center;

double fstart;
double fincr;
double c1;
double c2;
double cmult;
double minten;

double istart;
double istop;
double istep;
double ipre;
double itail;
double dci;
double scontrast;
double vcontrast;

double voltage;
double current;
double maxCurrent;
double pscal;
double elec_rs;
double elec_cap;
double somadia;

char savefile[30] = {0};

/*------------------------------------------------------*/

void defparams(void) { 

  setptr("cbptype",     &cbptype);
  setptr("cbptype2",    &cbptype2);
  setptr("ivplot",      &ivplot);
  setptr("cone_type",   &cone_type);
  setptr("bp_morph",    &bp_morph);
  setptr("amarr",       &amarr);
  setptr("dbp1_am",     &dbp1_am);
  setptr("dbp2_am",     &dbp2_am);
  setptr("no_tip_labels",&no_tip_labels);
  setptr("no_rpnt_labels",&no_rpnt_labels);
  setptr("disp_c2",     &disp_c2);
  setptr("nspots",      &nspots);
  setptr("stimtype",    &stimtype);
  setptr("mask_center", &mask_center);
  setptr("draw_synapse",&draw_synapse);
    
  setptrn("tip1",       &tip1);         /* Use predefined values (801,802, etc) from above, */
  setptrn("tip2",       &tip2);         /*  but allow user to change on command line.       */
  setptrn("tip3",       &tip3);
  setptrn("tip4",       &tip4);
  setptrn("tip5",       &tip5);
  setptrn("tip6",       &tip6);
  setptrn("tip7",       &tip7);
  setptrn("tip8",       &tip8);
  setptrn("tip9",       &tip9);
  setptrn("tip10",      &tip10);
  setptrn("tip11",      &tip11);
  setptrn("tip12",      &tip12);
  setptrn("tip13",      &tip13);
  setptrn("tip14",      &tip14);
  setptrn("tip15",      &tip15);
  setptrn("tip16",      &tip16);
  setptrn("tip17",      &tip17);
  setptrn("tip18",      &tip18);
  setptrn("tip19",      &tip19);
  setptrn("tip20",      &tip20);
  setptrn("tip21",      &tip21);
  setptrn("tip22",      &tip22);
  setptrn("tip23",      &tip23);
  setptrn("tip24",      &tip24);
  setptrn("tip25",      &tip25);
  setptrn("tip26",      &tip26);
  setptrn("tip27",      &tip27);
  setptrn("tip28",      &tip28);
  setptrn("tip29",      &tip29);
  setptrn("tip30",      &tip30);
    
  setptr("stimloc",     &stimloc);
  setptr("axon_base",   &axon_base);
  setptr("vcloc",       &vcloc);
  setptr("recpnt1",     &recpnt1);
  setptr("recpnt2",     &recpnt2);
  setptr("recpnt3",     &recpnt3);
  setptr("recpnt4",     &recpnt4);
  setptr("recpnt5",     &recpnt5);
  setptr("recpnt6",     &recpnt6);
  setptr("recpnt7",     &recpnt7);
  setptr("recpnt8",     &recpnt8);
  setptr("recpnt9",     &recpnt9);
  setptr("presyn1",     &presyn1);
  setptr("presyn2",     &presyn2);
  setptr("presyn3",     &presyn3);
  setptr("itransducer", &itransducer);
  setptr("dbp1_morph",  &dbp1_morph);
  setptr("npixels",     &npixels);
  setptr("nsines",      &nsines);
    
  setptr("drmab",       &drmab);
  setptr("driab",       &driab);
  setptr("dria",        &dria);
  setptr("cbp_rm",      &cbp_rm);
  setptr("dvreva",      &dvreva);
  setptr("amvst",       &amvst);
  setptr("amvrev",      &amvrev);
  setptr("amrm",        &amrm);
  setptr("dcmab",       &dcmab);
  setptr("axdia",       &axdia);
  setptr("axdiap",      &axdiap);
  setptr("ddia",        &ddia);
  setptr("naax",        &naax);
  setptr("naab",        &naab);
  setptr("nahd",        &nahd);
  setptr("ksoma",       &ksoma);
  setptr("kr6",         &kr6);
  setptr("kr7",         &kr7);
  setptr("kasoma",      &kasoma);
  setptr("kar6",        &kar6);
  setptr("kar7",        &kar7);
  setptr("dbp1_k1o",    &dbp1_k1o);
  setptr("khsoma",      &khsoma);
  setptr("amca",        &amca);
  setptr("amsca",       &amsca);
  setptr("am2ca",       &am2ca);
  setptr("amna",        &amna);
  setptr("amk",         &amk);
  setptr("casoma",      &casoma);
  setptr("cadist",      &cadist);
  setptr("cadist2",     &cadist2);
  setptr("catdist",     &catdist);
  setptr("catsoma",     &catsoma);
  setptr("capsoma",     &capsoma);
  setptr("kdr_cond",    &kdr_cond);
  setptr("g_am_dbp1",   &g_am_dbp1);
  setptr("g_am_dbp2",   &g_am_dbp2);
  setptr("g_am2_dbp1",  &g_am2_dbp1);
  setptr("g_am2_dbp2",  &g_am2_dbp2);
  setptr("g_dbp1_am",   &g_dbp1_am);
  setptr("g_dbp2_am",   &g_dbp2_am);
  setptr("g_dbp1_am2",  &g_dbp1_am2);
  setptr("g_dbp2_am2",  &g_dbp2_am2);
  setptr("g_dbp1_gca",  &g_dbp1_gca);
  setptr("elec_rs",     &elec_rs);
  setptr("elec_cap",    &elec_cap);
    
  setptr("varicos_dia", &varicos_dia);
  setptr("am_sdur",     &am_sdur);
  setptr("am_sdur2",    &am_sdur2);
  setptr("am2_sdur",    &am2_sdur);
  setptr("am2_sdur2",   &am2_sdur2);
  setptr("am_sfall",    &am_sfall);
  setptr("am_sfall2",   &am_sfall2);
  setptr("am2_sfall",   &am2_sfall);
  setptr("am2_sfall2",  &am2_sfall2);
  setptr("am_atau",     &am_atau);
  setptr("am2_atau",    &am2_atau);
  setptr("dbp1_gtau",   &dbp1_gtau);
  setptr("dbp2_gtau",   &dbp2_gtau);
  setptr("dbp1_g4tau",  &dbp1_g4tau);
  setptr("dbp2_g4tau",  &dbp2_g4tau);
  setptr("ampar",       &ampar);
  setptr("ampar2",      &ampar2);
  setptr("gabac",       &gabac);
  setptr("gabactrc",    &gabactrc);
  setptr("dbp1_mp",     &dbp1_mp);
  setptr("dbp1_mr",     &dbp1_mr);
  setptr("dbp1pca4",    &dbp1pca4);
  setptr("dbp1pca7",    &dbp1pca7);
  setptr("dbp2pca5",    &dbp2pca5);
  setptr("am2dbp1ca",   &am2dbp1ca);
  setptr("dbp1_dyad",   &dbp1_dyad);
    
  setptr("cone_soma_z", &cone_soma_z);
  setptr("dbp1_soma_z", &dbp1_soma_z);
  setptr("dbp2_soma_z", &dbp2_soma_z);
  setptr("spacing",     &spacing);
  setptr("dtreedia",    &dtreedia);
  setptr("dtreedia2",   &dtreedia2);
  setptr("c2rot",       &c2rot);
  setptr("c2yrot",      &c2yrot);
    
  setptr("predur",      &predur);
  setptr("prestimdur",  &prestimdur);
  setptr("spotdur",     &spotdur);
  setptr("spotduri",    &spotduri);
  setptr("stimdur",     &stimdur);
  setptr("tailcurdur",  &tailcurdur);
  setptr("poststimdur", &poststimdur);
  setptr("spotdia",     &spotdia);
    
  setptr("stimtime",    &stimtime);
  setptr("tfreq",       &tfreq);
  setptr("speriod",     &speriod);
    
  setptr("fstart",      &fstart);
  setptr("fincr",       &fincr);
  setptr("c1",          &c1);
  setptr("c2",          &c2);
  setptr("cmult",       &cmult);
  setptr("minten",      &minten);
    
  setptr("istart",      &istart);
  setptr("istop",       &istop);
  setptr("istep",       &istep);
  setptr("ipre",        &ipre);
  setptr("itail",       &itail);
  setptr("dci",         &dci);
  setptr("scontrast",   &scontrast);
  setptr("vcontrast",   &vcontrast);
    
  setptr("somadia",     &somadia);
     
  setptr("maxCurrent", &maxCurrent);
  setptr("pscal", &pscal);
    
  dbp1_file = "cbp_0576_t5os";
  am_file   = "morph_simple_am";
  gca_file  = "morph_simple_gca";

  nvalfile = "nval_cbp_inh.n";
  chanparamsfile = "chanparams_cbp_inh";
  dbp1_densfile  = "dens_cbp_inh.n";            // cbplam, axdia, dria
  am_densfile    = "dens_cbp_inh_am.n";
  // gca_densfile   = "dens_gca.n";
  gca_densfile   = "dens_gca_zero.n";


  g_dbp1_am  = 1e-10;   /* dbp1 cond, output to am, nval_cbp_inh.n */
  g_dbp1_am2 = 1e-10;   /* dbp1 cond, output to am2, nval_cbp_inh.n */
  g_dbp1_gca = 4e-10;   /* dbp1 cond, output to gca, nval_cbp_inh.n */
  g_dbp2_am  = 1e-10;   /* dbp2 cond, output to am, nval_cbp_inh.n */
  g_dbp2_am2 = 1e-10;   /* dbp2 cond, output to am2, nval_cbp_inh.n */
  g_am_dbp1  = 1e-14;   /* am cond, feedback to dbp1, nval_cbp_inh.n */
  g_am_dbp2  = 1e-10;   /* am cond, feedback to dbp2, nval_cbp_inh.n */
  g_am2_dbp1 = 1e-10;   /* am2 cond, feedback to dbp1, nval_cbp_inh.n */
  g_am2_dbp2 = 1e-10;   /* am2 cond, feedback to dbp2, nval_cbp_inh.n */
 
  if (notinit(dbp1_morph)) dbp1_morph = MORPH_REAL;
  if (dbp1_morph==MORPH_REAL)
     presyn1 = 823;        /* presynaptic locus for synaptic output (label in nval.n) */
  else
     presyn1 = 0;       /* presynaptic locus for synaptic output (label in nval.n) */
  presyn2 = 0;          /* presynaptic locus for dbp2 synaptic output (label in nval.n) */
  presyn3 = 0;          /* presynaptic locus for dbp1 2nd synaptic output (label in nval.n) */
  // recpnt1 = 0;          /* recording point (label) */
  // recpnt2 = 0;          /* recording point (label) */
  dbp1_gtau = 1;        /* tau multiplier for GABA channel in dbp1 */
  dbp2_gtau = 1;        /* tau multiplier for GABA channel in dbp2 */
  dbp1_g4tau = 1;       /* tau multiplier for GABAC channel in dbp1, also changes ligand sens */
  dbp2_g4tau = 1;       /* tau multiplier for GABAC channel in dbp2, also changes ligand sens */
  am_atau  = 1;         /* tau multiplier for AMPA channel in am */
  am2_atau = 1;         /* tau multiplier for AMPA channel in am2 */
  am_sdur  = 2;         /* sdur (ms) for am->dbp1 synapse */
  am_sdur2 = 2;         /* sdur (ms) for am->dbp2 synapse */
  am2_sdur = 2;         /* sdur (ms) for am2->dbp1 synapse */
  am2_sdur2 = 2;        /* sdur (ms) for am2->dbp2 synapse */
  am_sfall  = 0;        /* fall time const for am->dbp1 synapse */
  am_sfall2 = 0;        /* fall time const for am->dbp2 synapse */
  am2_sfall = 0;        /* fall time const for am2->dbp1 synapse */
  am2_sfall2 = 0;       /* fall time const for am2->dbp2 synapse */
  ampar  = xampa5;      /* AMPA response for dbp1->am synapse */
  ampar2 = xampa5;      /* AMPA response for dbp2->am synapse */
  gabac  = xgaba4;      /* GABA-c response for am2->dbp1 synapse */
  gabactrc = 4e-5;      /* trconc for gabac (try to set to give approx same fb as gaba-a) */
  
  dbp1_mp = 500;         /* maximum rr pool size */
  dbp1_mr = 500;        /* rr pool replenishment rate */
 
  //dbp1pca4 = 0.05;      /* dbp1->am  pca, sets gain of dbp1->am synapse */
  //dbp1pca7 = 0.05;      /* dbp1->am2 pca, sets gain of dbp1->am2 synapse */
  //dbp2pca5 = 0.05;      /* dbp2->am2 pca, sets gain of dbp2->am2 synapse */
  dbp1pca4 = 0.0;       /* dbp1->am  pca, sets gain of dbp1->am synapse */
  dbp1pca7 = 0.0;       /* dbp1->am2 pca, sets gain of dbp1->am2 synapse */
  dbp2pca5 = 0.0;       /* dbp2->am2 pca, sets gain of dbp2->am2 synapse */
  am2dbp1ca = 1;        /* am2->dbp1 sensca, sets Ca sens of am2->dbp1 synapse, Ca chans or AMPA perm */
  dbp1_dyad = 1;        /* sets dyad in dbp1 */

  cbp_rm = 5e3;
  gc_rm = 10e3;

  make_am_dbp1  = 1;
  make_dbp1_am2 = 1;
  make_dbp2_am  = 1;
  make_dbp2_gca = 0;
  make_am2_dbp1 = 1;            /* lateral feedback from second amacrine type to the first bp type */
  make_am3_dbp1 = 1;
  make_am_gca   = 0;
  make_am2_gca  = 0;

  _CA_L = _CA1;                 /* set type of L-type calcium channel for dens_cbp_inh.n */
  _CA_T = _CA7;                 /* set type of T-type calcium channel for dens_cbp_inh.n */

  dscavg = 1e6;                 /* Ca sensitivity for vesicle release */
  dscaeg = 1;			/* Ca release exponential gain (def 1) */
}

/*------------------------------------------------------*/

void setparams(void)
{
   make_ct (xcone);
   make_ct (dbp1);
   make_ct (am);
   make_ct (gca);
   n_cones = 1;
   n_dbp1 = 1;
   n_am   = 1;
   n_gca  = 1;

   if (make_cones) {
      if (notinit(cone_type)) cone_type=2;        /* 1 => cone phototransduction */
      setn(xcone,MORPH,cone_type);
      setn(xcone,SOMAZ,25);			  /* cone connects to dbp1 dendrites */
      setn(xcone,SCOND1,7e-9);	
      if (cone_type == 2) itransducer = 0;        /* cone -> vtransducer */
      else
      if (cone_type  > 2) itransducer = 1;        /* cone -> itransducer */
   }

   SOMA = R_3;                   /* defined in retsim_var.cc, retsim.h */

#define NCELLXY 10

   conexarr = (double *)emalloc(NCELLXY*sizeof(double));
   coneyarr = (double *)emalloc(NCELLXY*sizeof(double));
   dbp1xarr = (double *)emalloc(NCELLXY*sizeof(double));
   dbp1yarr = (double *)emalloc(NCELLXY*sizeof(double));
   amxarr   = (double *)emalloc(NCELLXY*sizeof(double));
   amyarr   = (double *)emalloc(NCELLXY*sizeof(double));
   gcxarr   = (double *)emalloc(NCELLXY*sizeof(double));
   gcyarr   = (double *)emalloc(NCELLXY*sizeof(double));

   if (dbp1_morph==MORPH_REAL) {

     double presyn1x = 0, presyn1y = -19;	/* loc of presyn1, node 264 */

     conexarr[0] = 0;
     coneyarr[0] = 0;
     dbp1xarr[0] = 0;
     dbp1yarr[0] = 0;
     amxarr[0] = 5;
     amyarr[0] = presyn1y;
     gcxarr[0] = 0;
     gcyarr[0] = presyn1y-5;
   } else {
     dbp1xarr[0] = 0;
     dbp1yarr[0] = 0;
     amxarr[0] = 5;
     amyarr[0] = 0;
     gcxarr[0] = 0;
     gcyarr[0] = -5;
   }

   remove_nconns = 0;

   setn(dbp1,MORPH,dbp1_morph); 
   setn(am,  MORPH,MORPH_SOMA);
   setn(gca, MORPH,MORPH_SOMA);

   if (dbp1_morph==MORPH_REAL) {
     setn(dbp1, SOMAZ,0);
     setn(am,   SOMAZ,-38);
     setn(gca,  SOMAZ,-38);
   } else {
     setn(dbp1, SOMAZ,0);
     setn(am,   SOMAZ,0);
     setn(gca,  SOMAZ,0);
   } 
   setn(dbp1, SOMADIA, 5);
   setn(am,   SOMADIA, 3);
   setn(gca,  SOMADIA, 5);

   if (dbp1_morph==MORPH_REAL) {
      setsv(dbp1, SYNREG,1,presyn1);
      setsv(dbp1, SYNREG,4,presyn1);
   } else {
      dispsize = 20;
      setsv(dbp1, SYNREG,1,-1);
      setsv(dbp1, SYNREG,4,-1);
   }
   setsv(am,   SYNREG,1,-1);
   setsv(am,   SYNREG,2,-1);

   setsv(dbp1, SYNREGP,1,-1);
   setsv(dbp1, SYNREGP,4,-1);
   setsv(am,   SYNREGP,1,-1);
   setsv(am,   SYNREGP,2,-1);

  if (notinit(cbplam))   cbplam = 0.02;         /* default complam for regions in density file */
  if (notinit(cbplam2)) cbplam2 = 0.2;          /* default complam for regions in density file 2 (large comps) */
  if (notinit(node_scale)) node_scale = -3.15;  /* 3: nodenum, 0.2 medium font, 0.05: small font */
  if (notinit(dbp1_nscale)) dbp1_nscale = -2.09;
  if (notinit(dbp2_nscale)) dbp2_nscale = -2.09;

  if (notinit(ddia))     ddia  = 1;             /* multiplier for dendrite diameter */
  if (notinit(axdia))   axdia  = 1;             /* multiplier for axon diameter */
  if (notinit(axdiap)) axdiap  = 1;             /* multiplier for proximal axon diameter */
  if (notinit(dvrev))   dvrev  = -0.065;        /* Vrev for dens_dbp1.n */
  if (notinit(dvst))    dvst   = -0.04;         /* Vstart for dens_dbp1.n */
  if (notinit(dvreva))  dvreva = -0.07;         /* Vrev for axon terminal in dens_dbp1.n */
  if (notinit(casoma))  casoma = 2e-3;          /* Ca chan density in soma, dens_cbp_inh.n */
  if (notinit(cadist))  cadist = 1.0e-3;        /* Ca chan density in axonal tips, dens_cbp_inh.n */
  if (notinit(cadist2)) cadist2 = 1.0e-3;       /* Ca chan density in axonal tips, dens_cbp2_inh.n */
  if (notinit(catdist)) catdist = 0e-3;         /* Ca T-type chan density in axonal tips, dens_cbp_inh.n */
  if (notinit(catsoma)) catsoma = 0e-3;         /* Ca T-type chan density in soma, dens_cbp_inh.n */
  if (notinit(capsoma)) capsoma = 2.0e-6;       /* Ca pump density in soma, dens_cbp_inh.n */
    
  if (notinit(amvst))     amvst = -0.065;       /* Vst  for dens_cbp_inh_am.n */
  if (notinit(amvrev))   amvrev = -0.065;       /* Vrev for dens_cbp_inh_am.n */
  if (notinit(amrm))       amrm = 20e3;         /* Rm for dens_cbp_inh_am.n */
  if (notinit(amca))       amca = 0e-3;         /* Ca chan density in dend for dens_cbp_inh_am.n, use AMPA Ca perm */
  if (notinit(amsca))     amsca = 0.2e-3;       /* Ca chan density in soma for dens_cbp_inh_am.n, use AMPA Ca perm */
  if (notinit(am2ca))     am2ca = 0.2e-3;       /* Ca chan density for dens_cbp_inh_am2.n */
  if (notinit(amna))       amna = 0e-3;         /* Na chan density for dens_cbp_inh_am.n */
  if (notinit(amk))         amk = 0e-3;         /* K  chan density for dens_cbp_inh_am.n */
    
  if (notinit(ivplot)) ivplot = 0;              /* make I/V plot */
  if (notinit(pscal)) pscal   = 50e-12;         /* plot scale */
    
  if (notinit(drmab)) drmab = drm;              /* user set default Rm for axon branches */
  if (notinit(driab)) driab = dri;              /* user set default Ri for axon branches */
  if (notinit(dria))  dria  = dri;              /* user set default Ri for axon branches */
  if (notinit(dcmab)) dcmab = dcm;              /* user set default Cm for axon branches */
  if (notinit(naax))   naax = 0;                /* user set Na density in axon */
  if (notinit(naab))   naab = 0;                /* user set Na density in axon branches */
  if (notinit(nahd))   nahd = 0e-3;             /* user set Na high density region in axon */
  if (notinit(ksoma)) ksoma = 0e-3;             /* Kdr chan density in soma, dens_cbp_inh.n */
  if (notinit(kr6))     kr6 = 0e-3;             /* Kdr chan density in R6, dens_cbp_inh.n */
  if (notinit(kr7))     kr7 = 0e-3;             /* Kdr chan density in R7, dens_cbp_inh.n */
  if (notinit(dbp1_k1o)) dbp1_k1o = 0.010;      /* Kdr chan voltage offset, chanparams_cbp_inh */
  if (notinit(kasoma)) kasoma = 0e-3;           /* KA chan density in soma, dens_cbp_inh.n */
  if (notinit(kar6))     kar6 = 0e-3;           /* KA chan density in R6, dens_cbp_inh.n */
  if (notinit(kar7))     kar7 = 0e-3;           /* KA chan density in R7, dens_cbp_inh.n */
  if (notinit(khsoma)) khsoma = 0e-3;           /* Kh chan density in soma, dens_cbp_inh.n */
 
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{

}

/*------------------------------------------------------*/

void addlabels(void)
{ 
}

/*------------------------------------------------------*/

void addcells() 
{
     double dia, rm;


}

/*------------------------------------------------------*/


void runonexit (void)
{
  if (savefile[0]!=0)  unlink (savefile);
}

/*------------------------------------------------------*/

void onplot(void) {
}

/*------------------------------------------------------*/

void runexpt(void) {

    int ct,cn;
    int transdloc, cashell, canode;
    double msize,mask,orient,time1,time2;
    double Vmin, Vmax, Imin, Imax;
    double cmin, cmax, fmin, fmax;
    double *rndarr = NULL;
    node *npnt;
    photorec *p;

    timinc = 2e-5;
    ploti = 1e-3;
    crit = 1e-9;

  ct = dbp1;
  if (notinit(itransducer))       itransducer = 0;

  if (!make_cones) {
    for(npnt=nodepnt; npnt=foreach(npnt,ct,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {  	// add transducer to dbp1  
        transdloc = soma;
	if (itransducer) {
	   if ((p=(photorec*)make_itransducer(ndn(ct,cn,transdloc))) != NULL) {
             p->xpos=npnt->xloc;                                                                     
             p->ypos=npnt->yloc;
	   }
	} else {
	   if ((p=(photorec*)make_transducer(ndn(ct,cn,transdloc))) != NULL) {
             p->xpos=npnt->xloc;                                                                     
             p->ypos=npnt->yloc;
	   }
	}
	if (ninfo >= 3) fprintf(stderr,"# ct  %s cn %-2d transducer node %d\n",cname[ct],cn, transdloc);
    }
  }

  if (notinit(prestimdur))   prestimdur = 1.0;
  if (notinit(spotdur))      spotdur    = 0.5;
  if (notinit(spotduri))    spotduri    = 0.5;
  if (notinit(stimdur))      stimdur    = 4.0;
  if (notinit(poststimdur)) poststimdur = 0.1;
  if (notinit(nspots))           nspots = 3;

  if (notinit(stimtime))       stimtime = prestimdur;
  if (notinit(tfreq))             tfreq = 2;
  if (notinit(speriod))         speriod = 10;
  if (notinit(fstart))           fstart = 0.5;
  if (notinit(fincr))             fincr = 0.00003;
  if (notinit(c1))                   c1 = 0;
  if (notinit(c2))                   c2 = 1;
  if (notinit(vcontrast))     vcontrast = 0;
  if (notinit(scontrast))     scontrast = 1;
  if (notinit(spotdia))         spotdia = 300;

  if (itransducer==1) {
      if (notinit(minten))           minten = 2e-12;
      if (notinit(cmult))             cmult = 1.0e-12;
  } else {
     if (notinit(minten))           minten = -0.040;
     if (notinit(cmult))             cmult = 0.005;
  }

/* - - - - - - - - - - - - - - - - - - - - - - - - */
/* plots */


   Vmax = -0.02;
   Vmin = -0.07;

   canode = dendn_node(dbp1,presyn1);

   if (make_cones) 
      plot_v_nod(xcone,1, 0,    -0.045, -0.035, magenta,"Vcone    ", 15, 0.5);
   plot_v_nod(dbp1,1, 0,    Vmin, Vmax, cyan,   "Vdbp_soma ", 12, 0.5);
   plot_v_nod(dbp1,1,canode,Vmin, Vmax, red,    "Vdbp1_presyn",12,0.5);
   plot_v_nod(am,  1, 0,    Vmin, Vmax, magenta,"Vam_soma  ", 4,  0.5);
//   plot_v_nod(gca, 1, 0,    Vmin, Vmax, brown,  "Vgc_soma  ", 2,  0.5);
   plot_i_nod(gca, 1, 0,    Imin= -2e-12, Imax=0e-12, brown,  "Igc     ", 2,  0.5);

   if (canode!=soma) Imin = -5e-15; else Imin = -2e-12;
   plot_chan_current(dbp1,1,canode,CA,1,Imin,Imax=0);  plot_param("dbp_Ica   ", blue, 11, 0.5);
   if (catdist > 0) {
	   plot_chan_current(dbp1,1,canode,CA,7,Imin= -2e-12,Imax=0);plot_param("dbp_Icat  ",magenta,11,0.5);
   }
   plot_ca_nod(am,   1, 0,       cashell=1, 2.0e-6, magenta, "am_Cai      ",    10,0.5);
   plot_ca_nod(dbp1, 1, canode,  cashell=1, 2.0e-6, green,   "dbp_Cai     ",    10,0.5);
   plot_syncond(findsyn(dbp1,1,canode, gca), cmin=0,cmax=50e-12,  green, 6,"g_dbp_gca  ",0.5);
   plot_syncond(findsyn(dbp1,1,canode, am),  cmin=0,cmax=50e-12,  cyan,  6,"g_dbp_am  ",0.5);
   plot_syncond(findsyn(am,1,0, dbp1) , cmin=0,cmax=50e-12,  red,        6,"g_am_dbp  ",0.5);

   plot_synrate(findsyn(dbp1,1,canode, gca), fmin=0,fmax=600,    blue,  8,"rate_dbp_gca",0.5);



/* - - - - - - - - - - - - - - - - - - - - - - - - */
/* stimuli */

   if (notinit(stimtype))       stimtype = 1;      

   if (stimtype == 1) { 		// flashed spot and chirp
           double x,y;
       s = 0;
       if (nspots==3) {
         stim_spot (spotdia, x=0, y=0, -cmult*scontrast, stimtime, spotduri);
       }
       stim_spot (spotdia, x=0, y=0,  cmult*scontrast, stimtime+spotduri+s*spotdur, spotdur); s++;
       stim_spot (spotdia, x=0, y=0, -cmult*scontrast, stimtime+spotduri+s*spotdur, spotdur); s++;
  			// mask out stim to central cones
       if (!notinit(mask_center)) stim_spot (30, x=0, y=0, mask_center, stimtime, spotduri+s*spotdur+stimdur, mask=1); 

       if (notinit (nsines)) nsines = 0;
       if (nsines>0) {			// discrete steps in frequency
            double f,tend;
 	    int i;
	 stimdur = 4.0;
         // tend = stimtime+spotduri+s*spotdur+prestimdur;
         tend = stimdur+prestimdur;	// first segment with steps
         for (f=fstart,i=0; i<nsines; i++, f*=fincr) {
              // tend = spot_sine_ncycle (spotdia, x=0, y=0, f, cmult, scontrast, tend,2);
	      tend = spot_sine(spotdia, x=0, y=0, f, cmult, scontrast, tend+prestimdur, stimdur);
         }
         time2 = tend;
       }
       else {				// continuous chirp
          time2 = spot_chirp (spotdia, x=0, y=0, fstart, fincr, cmult, scontrast, 
		       				stimtime+spotduri+s*spotdur+prestimdur, stimdur);
       }
       // fprintf (stderr,"time2+prestimdur %g\n",time2+prestimdur);
       // if (vcontrast>0) time2 = spot_vcontrast (spotdia, x=0, y=0, tfreq, 
       //					cmult, c1, c2, time2+prestimdur, stimdur);
       // if (vcontrast>0) time2 = spot_vdcontrast (spotdia, x=0, y=0, tfreq, 
       //  					cmult, c1, c2, time2+prestimdur, stimdur);
       if (vcontrast>0) time2 = spot_vfcontrast (spotdia, x=0, y=0, tfreq, 
         					cmult, c1, c2, time2, stimdur, prestimdur);
       poststimdur = 0.2;
       // fprintf (stderr,"end %g\n",time2);
   }

   // mask out left half of the field with a sector mask
   //
   if (stimtype == 2) { 		// annulus
       stimdur = 0.1;
       sector_mask(msize=180, orient=90, -0.05, stimtime, 3*spotdur+2*stimdur);
       stim_annulus (40,80,0,0,0.01,stimtime,3*spotdur);
       time2 = stimtime + stimdur;
   }
   
   // make random checkerboard 
   //
   if (stimtype == 3) { 		// checkerboard
 	   int nframes;
           double xoffset,yoffset;
       if (notinit(npixels)) npixels = 1;
       stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
         	          tfreq=100, 0, cmult*scontrast*2, stimtime, stimdur, &rndarr, &nframes, rseed);
       time2 = stimtime + stimdur;

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
   else if (stimtype==4) {               // sine waves 2,5,10,20,50 hz 
	   int drift;
	   double sqwave, sphase, orient;
      stim_sine (speriod, sphase=0, orient=0, 0, 0, tfreq, drift=1, 1, cmult*scontrast, sqwave=0, stimtime,  stimdur);
      time2 = stimdur+stimtime;
   }
   else if (stimtype==5) {               // sine waves 2,5,10,20,50 hz 
           double f,x,y,tend;
	   int i, nsines;
      nsines = 7;
      tend = stimtime+spotduri+s*spotdur+prestimdur;
      for (f=fstart,i=0; i<nsines; i++, f*=fincr) {
           tend = spot_sine_ncycle (spotdia, x=0, y=0, f, cmult, scontrast, tend,2);
      }
      time2 = tend;
      stimdur = tend - stimtime;
   }


/* - - - - - - - - - - - - - - - - - - - - - - - - */

   if (disp) {                        // display the stimulus
         double t, dscale, starttime, disp_end;

      stim_backgr(minten);		// set background for display at simtime=0
      display_size(500);
      disp_end = time2+0.05;
      for (starttime=simtime,t=stimtime; t<disp_end; starttime = t, t+= 0.01) {
           display_stim(starttime, t, dscale=4, -0.037, -0.043);
           //display_stim(t, dscale=4, -0.035, -0.045); 
           //display_stim(0+t, dscale=4); 
           simwait(0.10);
      }
      return;
   }
/* - - - - - - - - - - - - - - - - - - - - - - - - */

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 1.0;
  simtime = 0 - predur;
  if (stimtype==1) {
      // if (vcontrast>0) {endexp = stimtime+2*stimdur+nspots*spotdur+2*prestimdur+poststimdur;}
      // else             {endexp = stimtime+1*stimdur+nspots*spotdur+1*prestimdur+poststimdur;}
      if (vcontrast>0) {endexp = time2 + poststimdur;}
      else             {endexp = time2 + poststimdur;}
  }
  else                 {endexp = stimtime+stimdur+poststimdur;}
  // fprintf (stderr,"endexp %g\n", endexp);

  vclamp (ndn(gca,1,soma), vcl, simtime+predur, endexp);

  stim_backgr(minten);		// set background for run expt at simtime = -predur
  step (predur);
  step (endexp);
}

/*------------------------------------------------------*/
