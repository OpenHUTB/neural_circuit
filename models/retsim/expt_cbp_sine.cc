/* Experiment cbp_cclamp */
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
int no_recpnts;
int disp_c2;

int stimloc;
int axon_base;
int vcloc;
int itransducer;

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
int ampar;
int ampar2;
int hbp1_amh;
int hbp2_amh;

const char *cbptype;
const char *cbptype2;

double kdr_cond;
double hbp1_gtau;
double hbp2_gtau;
double amh_atau;
double amh2_atau;

double amh_sdur;
double amh2_sdur;
double amh_sfall;
double amh2_sfall;

double g_amh_hbp1;
double g_amh_hbp2;
double g_amh2_hbp1;
double g_amh2_hbp2;
double g_hbp1_amh;
double g_hbp1_amh2;
double g_hbp2_amh;
double g_hbp2_amh2;

double drmab;
double driab;
double dria;
double dvrevc;
double drmc;
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
double amna;
double amk;
double amrm;
double amvrev;
double hbp1_k1o;
double hbp2_k1o;
double cadist;
double cone_soma_z;
double hbp1_soma_z;
double hbp2_soma_z;
double spacing;
double c2rot;			// z rot for 2nd cell of stereo pair (sets disparity for side view)
double c2yrot;			// y rot for 2nd cell of stereo pair (sets disparity for bottom view)

double varicos_dia;

double predur;
double prestimdur;
double spotdur;
double stimdur;
double tailcurdur;
double poststimdur;
double stimtime;
double temp_freq;
double spotdia;
double spotloc;

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
double vcaxon;

double voltage;
double current;
double maxCurrent;
double pscal;
double elec_rs;
double elec_cap;




char savefile[30] = {0};

/*------------------------------------------------------*/

void defparams(void) { 

  setptr("cbptype",	&cbptype);
  setptr("cbptype2",	&cbptype2);
  setptr("ivplot",	&ivplot);
  setptr("cone_type",	&cone_type);
  setptr("bp_morph",	&bp_morph);
  setptr("amarr",	&amarr);
  setptr("hbp1_amh",	&hbp1_amh);
  setptr("hbp2_amh",	&hbp2_amh);
  setptr("no_tip_labels",&no_tip_labels);
  setptr("no_recpnts",  &no_recpnts);
  setptr("disp_c2",     &disp_c2);

  setptrn("tip1",	&tip1);		/* Use predefined values (801,802, etc) from above, */
  setptrn("tip2",	&tip2);		/*  but allow user to change on command line.       */
  setptrn("tip3",	&tip3);
  setptrn("tip4",	&tip4);
  setptrn("tip5",	&tip5);
  setptrn("tip6",	&tip6);
  setptrn("tip7",	&tip7);
  setptrn("tip8",	&tip8);
  setptrn("tip9",	&tip9);
  setptrn("tip10",	&tip10);
  setptrn("tip11",	&tip11);
  setptrn("tip12",	&tip12);
  setptrn("tip13",	&tip13);
  setptrn("tip14",	&tip14);
  setptrn("tip15",	&tip15);
  setptrn("tip16",	&tip16);
  setptrn("tip17",	&tip17);
  setptrn("tip18",	&tip18);
  setptrn("tip19",	&tip19);
  setptrn("tip20",	&tip20);
  setptrn("tip21",	&tip21);
  setptrn("tip22",	&tip22);
  setptrn("tip23",	&tip23);
  setptrn("tip24",	&tip24);
  setptrn("tip25",	&tip25);
  setptrn("tip26",	&tip26);
  setptrn("tip27",	&tip27);
  setptrn("tip28",	&tip28);
  setptrn("tip29",	&tip29);
  setptrn("tip30",	&tip30);

  setptr("stimloc",	&stimloc);
  setptr("axon_base",	&axon_base);
  setptr("vcloc",	&vcloc);
  setptr("recpnt1",	&recpnt1);
  setptr("recpnt2",	&recpnt2);
  setptr("recpnt3",	&recpnt3);
  setptr("recpnt4",	&recpnt4);
  setptr("recpnt5",	&recpnt5);
  setptr("recpnt6",	&recpnt6);
  setptr("recpnt7",	&recpnt7);
  setptr("recpnt8",	&recpnt8);
  setptr("recpnt9",	&recpnt9);
  setptr("presyn1",	&presyn1);
  setptr("presyn2",	&presyn2);
  setptr("itransducer",	&itransducer);

  setptr("drmab",       &drmab);
  setptr("driab",       &driab);
  setptr("dria",        &dria);
  setptr("dvrevc",      &dvrevc);
  setptr("drmc",        &drmc);
  setptr("dcmab",       &dcmab);
  setptr("axdia",       &axdia);
  setptr("axdiap",      &axdiap);
  setptr("ddia",        &ddia);
  setptr("naax",        &naax);
  setptr("naab",        &naab);
  setptr("nahd",        &nahd);
  setptr("ksoma",	&ksoma);
  setptr("kr6",		&kr6);
  setptr("kr7",		&kr7);
  setptr("amna",	&amna);
  setptr("amk",		&amk);
  setptr("amrm",	&amrm);
  setptr("amvrev",	&amvrev);
  setptr("hbp1_k1o",	&hbp1_k1o);
  setptr("hbp2_k1o",	&hbp2_k1o);
  setptr("cadist",	&cadist);
  setptr("kdr_cond",    &kdr_cond);
  setptr("g_amh_hbp1",   &g_amh_hbp1);
  setptr("g_amh_hbp2",   &g_amh_hbp2);
  setptr("g_amh2_hbp1",  &g_amh2_hbp1);
  setptr("g_amh2_hbp2",  &g_amh2_hbp2);
  setptr("g_hbp1_amh",   &g_hbp1_amh);
  setptr("g_hbp2_amh",   &g_hbp2_amh);
  setptr("g_hbp1_amh2",  &g_hbp1_amh2);
  setptr("g_hbp2_amh2",  &g_hbp2_amh2);
  setptr("elec_rs",	&elec_rs);
  setptr("elec_cap",	&elec_cap);

  setptr("varicos_dia", &varicos_dia);
  setptr("amh_sdur",     &amh_sdur);
  setptr("amh2_sdur",    &amh2_sdur);
  setptr("amh_sfall",    &amh_sfall);
  setptr("amh2_sfall",   &amh2_sfall);
  setptr("amh_atau",     &amh_atau);
  setptr("amh2_atau",    &amh2_atau);
  setptr("hbp1_gtau",   &hbp1_gtau);
  setptr("hbp2_gtau",   &hbp2_gtau);
  setptr("ampar",       &ampar);
  setptr("ampar2",      &ampar2);

  setptr("cone_soma_z",	&cone_soma_z);
  setptr("hbp1_soma_z",	&hbp1_soma_z);
  setptr("hbp2_soma_z",	&hbp2_soma_z);
  setptr("spacing",     &spacing);
  setptr("c2rot",       &c2rot);
  setptr("c2yrot",      &c2yrot);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("spotdur",	&spotdur);
  setptr("stimdur",	&stimdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststimdur",	&poststimdur);
  setptr("spotdia",	&spotdia);
  setptr("spotloc",	&spotloc);

  setptr("stimtime",	&stimtime);
  setptr("temp_freq",	&temp_freq);
  
  setptr("fstart",	&fstart);
  setptr("fincr",	&fincr);
  setptr("c1",		&c1);
  setptr("c2",		&c2);
  setptr("cmult",	&cmult);
  setptr("minten",	&minten);

  setptr("istart",	&istart);
  setptr("istop",	&istop);
  setptr("istep",	&istep);
  setptr("ipre",	&ipre);
  setptr("itail",	&itail);
  setptr("dci",		&dci);
  setptr("scontrast",	&scontrast);
  setptr("vcontrast",	&vcontrast);
  setptr("vcaxon",	&vcaxon);

  setptr("maxCurrent", &maxCurrent);
  setptr("pscal", &pscal);

  nvalfile = "nval_cbp_sine.n";
  chanparamsfile = "chanparams_cbp_sine";

  // hbp1_file = "morph_bp";			// default set in retsim.cc 
  // hbp1_file = "morph_DB_111005_12_db4";	// morphology file in retsim.cc 
  // hbp1_file = "cell504_t4";			// morphology file in retsim.cc 
  // hbp1_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 
  // hbp2_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 
  hbp1_file = "cbp_0572_t5t";			// morphology file in retsim.cc 
  hbp2_file = "cbp_0572_t5t";			// morphology file in retsim.cc 
  hbp1_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 

  hbp1_densfile  = "dens_cbp_sine.n";		// cbplamh, axdia, dria
  hbp1_densfile2 = "dens_cbp_sinec.n";		// cbplamh2, axdia, dria
  hbp2_densfile  = "dens_cbp_sinec.n";
  hbp2_densfile2 = "dens_cbp_sinec.n";
  hbp1_densfile  = "dens_cbp_sine.n";
  hbp1_densfile2 = "dens_cbp_sinec.n";
  amh_densfile   = "dens_amh_sine.n";
  amh2_densfile  = "dens_amh2_sine.n";

  if (notinit(g_amh_hbp1))  g_amh_hbp1  = 4e-10; 	/* am cond, feedback to hbp1, nval_cbp_sine.n */
  if (notinit(g_amh_hbp2))  g_amh_hbp2  = 4e-10; 	/* am cond, feedback to hbp2, nval_cbp_sine.n */
  if (notinit(g_amh2_hbp1)) g_amh2_hbp1 = 4e-10; 	/* amh2 cond, feedback to hbp1, nval_cbp_sine.n */
  if (notinit(g_amh2_hbp2)) g_amh2_hbp2 = 4e-10; 	/* amh2 cond, feedback to hbp2, nval_cbp_sine.n */
  if (notinit(g_hbp1_amh))  g_hbp1_amh  = 5e-10;	/* hbp1 cond, output to am, nval_cbp_sine.n */
  if (notinit(g_hbp1_amh2)) g_hbp1_amh2 = 5e-10; 	/* hbp1 cond, output to amh2, nval_cbp_sine.n */
  if (notinit(g_hbp2_amh))  g_hbp2_amh  = 5e-10; 	/* hbp2 cond, output to am, nval_cbp_sine.n */
  if (notinit(g_hbp2_amh2)) g_hbp2_amh2 = 5e-10; 	/* hbp2 cond, output to amh2, nval_cbp_sine.n */
  if (notinit(presyn1)) presyn1 = 0;  		/* presynaptic locus for synaptic output (label in nval.n) */
  if (notinit(presyn2)) presyn2 = 0;  		/* presynaptic locus for synaptic output (label in nval.n) */
  // if (notinit(recpnt1)) recpnt1 = 0; 	 	/* recording point (label) */
  // if (notinit(recpnt2)) recpnt2 = 0; 	 	/* recording point (label) */
  if (notinit(hbp1_gtau)) hbp1_gtau = 1; 	/* tau multiplier for GABA channel in hbp1 */
  if (notinit(hbp2_gtau)) hbp2_gtau = 1; 	/* tau multiplier for GABA channel in hbp2 */
  if (notinit(amh_atau))     amh_atau  = 1; 	/* tau multiplier for AMPA channel in am */
  if (notinit(amh2_atau))    amh2_atau = 1; 	/* tau multiplier for AMPA channel in amh2 */
  if (notinit(amh_sdur))     amh_sdur  = 2; 	/* sdur for am->hbp1 synapse */
  if (notinit(amh2_sdur))    amh2_sdur = 2; 	/* sdur for amh2->hbp1 synapse */
  if (notinit(amh_sfall))    amh_sfall  = 0; 	/* fall time const for am->hbp1 synapse */
  if (notinit(amh2_sfall))   amh2_sfall = 0; 	/* fall time const for amh2->hbp1 synapse */
  if (notinit(ampar))      ampar  = xampa5; 	/* AMPA response for hbp1->am synapse */
  if (notinit(ampar2))     ampar2 = xampa5; 	/* AMPA response for hbp2->am synapse */

  make_hbp1_amh = 1;
  make_hbp2_amh  = 1;
  make_hbp2_amh2  = 1;
  make_hbp2_gcaoff = 0;
  make_amh_hbp1 = 1;		/* lateral feedback from first amacrine type to the first bp type */
  make_amh2_hbp1 = 1;		/* lateral feedback from second amacrine type to the first bp type */
  make_amh2_hbp2 = 1;
  // make_amh2_amh = 1;
  make_amh_gcaoff  = 0;
}

/*------------------------------------------------------*/

int makrow (int ncells, double spacing, double rowtheta, double xloc, double yloc, double zrot, 
		double  *xarr, double *yarr, double *tharr, int center, int i, int s)

/* make one row of cells along an orientation */
/* i counts number of cells in xarr[], yarr[]. */

#define R 5
{
  int c,j,k,n;
  int r[20];
  // static int r[] = {R,R-1,R+1,R-2,R+2,R-3,R+3,R-4,R+4,R-5,R+5,R-6,R+6,R-7,R+7,R-8,R+8};

  c = center;
  r[0] = c;
  for (j=1,n=1;j<18; n++,j+=2) {
       r[j]   = c-n;				// alternating either side of center
       r[j+1] = c+n;
  }
  // for (j=0; j<ncells; j++) {
  //      fprintf (stderr,"j %d r %d\n",j,r[j]);
  // }
  for (j=0;j<(ncells-s); i++,j++) {
       k = r[j+s];  
       // fprintf (stderr,"c %d i %d j %d k %d\n",c,i,j,k);
       xarr[i]  = -cosdeg(rowtheta) * ((c-k)*spacing) + xloc;
       yarr[i]  =  sindeg(rowtheta) * ((c-k)*spacing) + yloc;
       tharr[i] = zrot-rowtheta;
  }
  return i;
}
#undef R

/* - - - - - - - - - - - - - - - - - - - */

int hbp1row (int ncells, double spacing, double rowtheta, double zrot, double xloc, double yloc, int center, int i, int s)
{
      double dist;		// radial distance from (xloc, yloc)

    return makrow (ncells, spacing, rowtheta, xloc, yloc, zrot, hbp1xarr, hbp1yarr, hbp1tharr, center, i, s);
}

/* - - - - - - - - - - - - - - - - - - - */

int hbp1row (int ncells, double spacing, double rowtheta, double zrot, int center, int i, int s)
{

    return makrow (ncells, spacing, rowtheta, 0, 0, zrot, hbp1xarr, hbp1yarr, hbp1tharr, center, i, s);
}

/* - - - - - - - - - - - - - - - - - - - */

int hbp2row (int ncells, double spacing, double rowtheta, double zrot, double xloc, double yloc, int center, int i, int s)
{
      double dist;		// radial distance from (xloc, yloc)

    return makrow (ncells, spacing, rowtheta, xloc, yloc, zrot, hbp2xarr, hbp2yarr, hbp2tharr, center, i,s);
}

/* - - - - - - - - - - - - - - - - - - - */

int hbp2row (int ncells, double spacing, double rowtheta, double zrot, int center, int i, int s)
{

    return makrow (ncells, spacing, rowtheta, 0, 0, zrot, hbp2xarr, hbp2yarr, hbp2tharr, center, i,s);
}

/* - - - - - - - - - - - - - - - - - - - */

int hbp1cent (int ncells, double theta, double zrot, int center, int i)
{
      int s;			// start
      double spacing;		// cell spacing along row (here, none)

    return makrow (ncells, spacing=0, theta, 0, 0, zrot, hbp1xarr, hbp1yarr, hbp1tharr, center, i, s=0);
}
/* - - - - - - - - - - - - - - - - - - - */

int hbp2cent (int ncells, double theta, double zrot, int center, int i)
{
      int s;			// start
      double spacing;		// cell spacing along row (here, none)

    return makrow (ncells, spacing=0, theta, 0, 0, zrot, hbp2xarr, hbp2yarr, hbp2tharr, center, i, s=0);
}
/* - - - - - - - - - - - - - - - - - - - */

int make_amh_cell (double theta, int i, double somadist)
{

   make_ct(amh);
   amhxarr[i]  = -cosdeg(theta) * somadist;
   amhyarr[i]  =  sindeg(theta) * somadist;
   amhtharr[i++] = theta;
   return (i);
}

/* - - - - - - - - - - - - - - - - - - - */

int make_amh2_cell (double theta, int i, double somadist)
{
   make_ct(amh2);
   amh2xarr[i]  = -cosdeg(theta) * somadist;
   amh2yarr[i]  =  sindeg(theta) * somadist;
   amh2tharr[i++] = theta;
   return (i);
}

/*------------------------------------------------------*/
   
void setparams(void)
{
#define CBPARR 100
#define AMARR  20

     int c,i,s, dist, cent;
     double theta;						/* orientation of amacrine cell to follow */
     double hbp1_zrot=0, hbp2_zrot=0; 				/* rotation of cell at amacrine cell */

  if (!notinit(cbptype)) {			/* primary cbp type */
       ct = find_ct(cbptype);
  } else {
       ct = hbp1;
  }

  if (!notinit(cbptype2)) {			/* secondary cbp type */
       ct2 = find_ct(cbptype2);
  } else {
       ct2 = hbp2;
  }

  if (!notinit(cone_soma_z)) setn(xcone,SOMAZ,cone_soma_z);
  if (!notinit(hbp1_soma_z)) setn(hbp1,SOMAZ,hbp1_soma_z);
  if (!notinit(hbp2_soma_z)) setn(hbp2,SOMAZ,hbp2_soma_z);
  if (notinit(spacing)) spacing = 40;

  make_ct(ct);
  if (ct2>0) make_ct(ct2);
  // make_ct(gcaoff);

  SOMA = R_3;			/* defined in retsim_var.cc, retsim.h */

  if (notinit(hbp1_amh)) hbp1_amh = 11;
  if (notinit(amarr)) amarr = 0;

  cent = int((hbp1_amh+1)/2);				    /* the central hbp */
  hbp1xarr   = (double *)emalloc(CBPARR*sizeof(double));   /* used by makrow() above */
  hbp1yarr   = (double *)emalloc(CBPARR*sizeof(double));
  hbp1tharr  = (double *)emalloc(CBPARR*sizeof(double));
  hbp1cent (1,      theta=  0,   hbp1_zrot,   cent, i=0);
  if (!hbp2xarr) { 
	hbp2xarr   = (double *)emalloc(CBPARR*sizeof(double));
	hbp2yarr   = (double *)emalloc(CBPARR*sizeof(double));
	hbp2tharr  = (double *)emalloc(CBPARR*sizeof(double));
  }

  switch (amarr) {						/* set up hbp1s */

     case 0: 
	if      (n_hbp1<=1) { } 
	else if (n_hbp1==2) {
	  hbp1xarr[0] = 0;
	  hbp1yarr[0] = 0;
	  hbp1tharr[0] = 0;
	  hbp1xarr[1] = spacing;
	  hbp1yarr[1] = 0;
	  if (notinit(c2rot)) c2rot = 5;
	  hbp1tharr[1] = c2rot;
	  if (c2yrot!=0) {		// if stereo pair seen from bottom
             hbp1ytharr  = (double *)emalloc(CBPARR*sizeof(double));
	     hbp1ytharr[0] = 0;
	     hbp1ytharr[1] = -c2yrot;
	  }
	}
	//remove_nconns = 0;
	break;

     case 1:	
	 // hbp1_amh = 6;
	 if (n_hbp1 < 0 || n_hbp1 > hbp1_amh) n_hbp1 = hbp1_amh;			/* limit to 6 hbp1s */
	 i = hbp1row (hbp1_amh, spacing, theta=  180, hbp1_zrot, cent,  i=1,s=0);
	break; 

     case 2:
	break;

     case 11:	
	  if (n_hbp1 < 0 || n_hbp1 > hbp1_amh) n_hbp1 = hbp1_amh;			/* limit to 11 hbp1s */
	  i = hbp1row (n_hbp1, spacing, theta=  0,   hbp1_zrot, cent,  i=1,s=1);
	  break;
 
     case 21:

	  if (n_hbp1 < 0 || n_hbp1 > 2*hbp1_amh) n_hbp1 = 2*hbp1_amh - 1; 	  /* limit to 21 hbp1s */
	  // i = hbp1cent (1,                  theta=  45, hbp1_zrot,  cent, i=0);	  /* set cell rot near both ams */
	  i = hbp1row (hbp1_amh, spacing, theta=  90, hbp1_zrot,  cent, i=1,s=1);
	  i = hbp1row (hbp1_amh, spacing, theta=   0, hbp1_zrot,  cent, i,  s=1);
	  setn(hbp1,MAXSDIST,15); 
	  break;
 
     case 22: 									/* use with hbp1, hbp2 */
	  if (n_hbp1 < 0 || n_hbp1 > hbp1_amh) n_hbp1 = hbp1_amh;			/* limit to 11 hbp1s */
	  i = hbp1cent (1,      theta=  0,   hbp1_zrot+315,         cent, i=0);
	  // i = hbp2row (hbp2_amh, spacing, theta= 0, hbp1_zrot+55, cent, i=1,s=1);
	  setn(hbp1,MAXSDIST,25); 
	  break;

     case 23: 									/* use with hbp1, hbp2 */
          if (n_hbp1 >= 0) {			/* random array */
             hbp1xarrsiz = hbp1_amh * spacing;
             hbp1yarrsiz = 30;
 	     hbp1_first_cent = 1;
          } 
	  break;

     case 3:
	  if (n_hbp1 < 0 || n_hbp1 > 3*hbp1_amh-2) n_hbp1 = 3*hbp1_amh-2;		/* limit to 31 hbp1s */
	  i = hbp1cent (1,                  theta=  45, hbp1_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = hbp1row (hbp1_amh, spacing, theta=  0,  hbp1_zrot, cent, i=1,s=1);
	  i = hbp1row (hbp1_amh, spacing, theta=  60, hbp1_zrot, cent, i,  s=1);
	  i = hbp1row (hbp1_amh, spacing, theta= -60, hbp1_zrot, cent, i,  s=1);
	  setn(hbp1,MAXSDIST,14);
	  break;

     case 4:
	  if (n_hbp1 < 0 || n_hbp1 > 4*hbp1_amh-3) n_hbp1 = 4*hbp1_amh-3;		/* limit to 41 hbp1s */
	  i = hbp1cent (1,                  theta=  45, hbp1_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = hbp1row (hbp1_amh, spacing, theta=   0, hbp1_zrot, cent, i=1,s=1);
	  i = hbp1row (hbp1_amh, spacing, theta=  45, hbp1_zrot, cent, i, s=1);
	  i = hbp1row (hbp1_amh, spacing, theta= -45, hbp1_zrot, cent, i, s=1);
	  i = hbp1row (hbp1_amh, spacing, theta=  90, hbp1_zrot, cent, i, s=1);
	  setn(hbp1,MAXSDIST,14.7);
	 break;

     default: break;
 
  }  /* switch (amarr) */

  setn(hbp1,DENS,pow(spacing*0.5,-2)*1e6);
 
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (notinit(hbp2_amh)) {
     if (!notinit(hbp1_amh)) hbp2_amh = hbp1_amh;
     else hbp2_amh = 11;
  }

  if (amarr > 0) {
     cent = int((hbp2_amh+1)/2);					/* skip the central hbp (when needed) */
     if (!hbp2xarr) { 
	hbp2xarr   = (double *)emalloc(CBPARR*sizeof(double));
	hbp2yarr   = (double *)emalloc(CBPARR*sizeof(double));
	hbp2tharr  = (double *)emalloc(CBPARR*sizeof(double));
     }
     // hbp2cent (1,      theta=  0,   hbp2_zrot,   cent, i=0);
  }

  if (n_hbp2 > 0) {
    switch (amarr) {						/* set up hbp2s */

     case 0:
	  n_hbp2 = 1;								/* just 1 hbp2 */
          hbp2cent (1,      theta=  0,   hbp2_zrot,   cent, i=0);
	  // remove_nconns = 0;
	  break;

     case 1:
	  hbp2_amh = 6;
	  if (n_hbp2 > hbp2_amh+1) n_hbp2 = hbp2_amh;				/* limit to 5 hbp2s */
	  i = hbp2row (hbp2_amh, spacing, theta=  180, hbp2_zrot, cent,  i=1,s=0);
	  break;

     case 11:
	  if (n_hbp2 > hbp2_amh) n_hbp2 = hbp2_amh;				/* limit to 11 hbp2s */
	  i = hbp2row (hbp2_amh, spacing, theta=  0,   hbp2_zrot, cent,  i=0,s=1);
	  break;

     case 21: 
	  if (n_hbp2 > 2*hbp2_amh) n_hbp2 = 2*hbp2_amh - 1; 			   /* limit to 21 hbp2s */
	  i = hbp2cent (1,                  theta=  45, hbp2_zrot,  cent,  i=0); /* set cell rot near both ams */
	  i = hbp2row (hbp2_amh, spacing, theta=  90, hbp2_zrot,  cent,  i=1,s=1);
	  i = hbp2row (hbp2_amh, spacing, theta=   0, hbp2_zrot,  cent,  i,  s=1);
	  setn(hbp2,MAXSDIST,15); 
	  break;

     case 22: 		 							/* use with hbp1, hbp2 */
	 if (n_hbp2 > hbp2_amh) n_hbp2 = hbp2_amh;				/* limit to 10 hbp2s */
	 i = hbp2row (hbp2_amh, spacing, theta= 0, hbp1_zrot+55, cent, i=0,s=1);
	 i = hbp2row (hbp2_amh, spacing, theta= 90,   hbp2_zrot, cent, i,  s=1);
	 break;

     case 23:	
	  if (n_hbp2 >= 0) {				/* random array */
             hbp2xarrsiz = 30;
             hbp2yarrsiz = hbp2_amh * spacing;
	     hbp2_first_cent = 0;
	     setn(hbp2,MAXSDIST,15); 
	  }
	 break;

     case 3:
	  if (n_hbp2 > 3*hbp2_amh+1) n_hbp2 = 3*hbp2_amh - 2;			   /* limit to 31 hbp2s */
	  i = hbp2cent (1,                  theta=  45, hbp2_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = hbp2row (hbp2_amh, spacing, theta=  0,  hbp2_zrot, cent, i=1,s=1);
	  i = hbp2row (hbp2_amh, spacing, theta=  60, hbp2_zrot, cent, i,  s=1);
	  i = hbp2row (hbp2_amh, spacing, theta= -60, hbp2_zrot, cent, i,  s=1);
	  setn(hbp2,MAXSDIST,14);
	  break;

     case 4:
	  if (n_hbp2 > 4*hbp2_amh+1) n_hbp2 = 4*hbp2_amh - 3;			   /* limit to 41 hbp2s */
	  i = hbp2cent (1,                  theta=  45, hbp2_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = hbp2row (hbp2_amh, spacing, theta=   0, hbp2_zrot, cent, i=1,s=1);
	  i = hbp2row (hbp2_amh, spacing, theta=  45, hbp2_zrot, cent, i,  s=1);
	  i = hbp2row (hbp2_amh, spacing, theta= -45, hbp2_zrot, cent, i,  s=1);
	  i = hbp2row (hbp2_amh, spacing, theta=  90, hbp2_zrot, cent, i,  s=1);
	  setn(hbp2,MAXSDIST,14.7);
          break;

    }       /* switch (amarr) */
  }       /* if (nhbp2 > 0) */

  setn(hbp2,DENS,pow(spacing*0.5,-2)*1e6);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (!notinit(amarr)) {
        int i,j;
        double theta;
	double dtreedia;
	double somadist;

   if (amarr > 0) make_ct(amh);

   setn(ct,AXARBT,BRANCHED);		/* connect to bipolar cell branched axonal arbor */
   if (ct2>0) setn(ct2,AXARBT,BRANCHED);	/* connect to bipolar cell branched axonal arbor */

   amhxarr  = (double *)emalloc(AMARR*sizeof(double));
   amhyarr  = (double *)emalloc(AMARR*sizeof(double));
   amhtharr = (double *)emalloc(AMARR*sizeof(double));

   amh2xarr  = (double *)emalloc(AMARR*sizeof(double));
   amh2yarr  = (double *)emalloc(AMARR*sizeof(double));
   amh2tharr = (double *)emalloc(AMARR*sizeof(double));

   // dtreedia = n_hbp1 * spacing;
   /*
    switch (amarr) {
	    default:
	    case 1:
	    case 11:
	    case 22:
	    case 21: spacing = 32; break;
	    case 3:  spacing = 35; break;
	    case 4:  spacing = 36; break;
   }
   */

   dtreedia = hbp1_amh * spacing * 2;
   if (dtreedia <= 0) dtreedia = 650;
   somadist = dtreedia * 0.3;

   if (ninfo > 2) ncfprintf (stderr,"# setparams(): amh dtreedia %g somadist %g\n",dtreedia,somadist);

   setn(amh, DTREEDIA,dtreedia);
   setn(amh2,DTREEDIA,dtreedia);

   switch (amarr) {

     case 1:
	  if (n_hbp1 <= 6 && n_hbp2 <=6) setn(amh,DTREEDIA,dtreedia);
          make_amh_cell (theta=180, i=0, somadist);
	  n_amh = 1;
          break;
   
    case 2: 
	  // if (n_hbp1 <= 6 && n_hbp2 <=6) setn(amh,DTREEDIA,dtreedia);
          // make_amh_cell (theta=180, i=0, somadist);
	  n_amh = 0;
          break;

    case 11:
	  if (n_hbp1 <= 11 && n_hbp2 <=11) setn(amh,DTREEDIA,dtreedia);
          make_amh_cell (theta=180, i=0, somadist);
	  n_amh = 1;
          break;

    case 21:
          i = make_amh_cell (theta=180, i=0, somadist);
          i = make_amh_cell (theta=90,  i,   somadist);
          n_amh = 2;
          break;

    case 22:
          make_amh_cell  (theta=180, i=0, somadist);
          make_amh2_cell (theta=90,  i=0, somadist);
	  n_amh  = 1;
	  n_amh2 = 1;
          break;

    case 23:
          i = make_amh_cell (theta=135, i=0, somadist);
          i = make_amh_cell (theta=45,  i,   somadist);
          n_amh = 2;
          break;

    case 3:
          i = make_amh_cell (theta=0,   i=0, somadist);
          i = make_amh_cell (theta=60,  i,   somadist);
          i = make_amh_cell (theta=-60, i,   somadist);
          n_amh = 3;
          break;
 
    case 4:
          i = make_amh_cell (theta=0,   i=0, somadist);
          i = make_amh_cell (theta=45,  i, somadist);
          i = make_amh_cell (theta=-45, i, somadist);
          i = make_amh_cell (theta=90,  i, somadist);
          n_amh = 4;
          break;
   }  /* switch (amarr) */

  }  /* if (!notinit(amarr)) */

/*
  if (n_amh > 0) {
	  make_ct(amh);
          setn(ct,AXARBT,BRANCHED);	// connect to bipolar cell branched axonal arbor 
	  amhxarr  = (double *)emalloc(AMARR*sizeof(double));
	  amhyarr  = (double *)emalloc(AMARR*sizeof(double));
	  amhtharr = (double *)emalloc(AMARR*sizeof(double));
	  
	  amhxarr[0]  = 150;
	  amhyarr[0]  = 0;
	  amhtharr[0] = 180;

	  amhxarr[1]  = 0;
	  amhyarr[1]  = 150;
	  amhtharr[1] = 90;

	  amhxarr[2]  = 100;
	  amhyarr[2]  = 100;
	  amhtharr[2] = 135;

	  amhxarr[3]  = -100;
	  amhyarr[3]  = 100;
	  amhtharr[3] = 45;
  } /* */

  if (make_cones) {
    if (notinit(cone_type)) setn(xcone,MORPH,cone_type=1); /* 1 => cone phototransduction */
    setn(xcone,MORPH,cone_type);
    if (cone_type == 2) itransducer = 0; 	/* cone -> vtransducer */
    else 
    if (cone_type  > 2) itransducer = 1;	/* cone -> itransducer */

  } 
  if (notinit(itransducer)) itransducer = 1;	/* use cclamp for transducer, not vclamp */

  if (notinit(bp_morph)) setn(ct,MORPH,0);		/* 0 => set cell morphology from file, default = "morph_bp" */
  else                   setn(ct,MORPH,bp_morph);
  setn(ct,BIOPHYS,1);		/* set cell biophys from file, default = "dens_default" */
  setn(ct,NCOLOR,RCOLOR);	/* set cell display color from region */

  // if (notinit(cone_timec)) cone_timec = 0.5;		/* time constant for cone transduction, ~100 ms */
  // if (notinit(cone_loopg)) cone_loopg = 0.8;		/* loop gain (Ca feedbac) for cone, set for no osc. */

  if (ct2>0) {
    if (notinit(bp_morph)) setn(ct2,MORPH,0);		/* set cell morphology from file, default = "morph_bp" */
    else                   setn(ct2,MORPH,bp_morph);
    setn(ct2,BIOPHYS,1);				/* set cell biophys from file, default = "dens_default" */
    // setn(ct2,NCOLOR,RCOLOR);				/* set cell display color from region */
  }
  if (notinit(no_tip_labels)) no_tip_labels = 0;	/* = 1 => display recording pts but not tip labels */
  if (notinit(no_recpnts))    no_recpnts = 0;	/* = 1 => display tip labels but not recording pts */
  if (notinit(disp_c2))       disp_c2 = 1;	/* = 1 => display plots for cell 2 */

  if (notinit(cbplam))   cbplam = 0.1; 		/* default complam for regions in density file */
  if (notinit(cbplam2)) cbplam2 = 0.1; 		/* default complam for regions in density file 2 (large comps) */
  // if (notinit(dispsize)) dispsize = 80; 	/* default display size */
  if (notinit(node_scale)) node_scale = -3.15;  /* 3: nodenum, 0.2 medium font, 0.05: small font */
  if (notinit(hbp1_nscale)) hbp1_nscale = -2.09;
  if (notinit(hbp2_nscale)) hbp2_nscale = -2.09;

  if (notinit(ddia))     ddia  = 1; 		/* multiplier for dendrite diameter */
  if (notinit(axdia))   axdia  = 1; 		/* multiplier for axon diameter */
  if (notinit(axdiap)) axdiap  = 1; 		/* multiplier for proximal axon diameter */
  if (notinit(dvrev))   dvrev  = -0.065;	/* Vrev for dens_cbp_sine.n */
  if (notinit(dvst))    dvst   = -0.07;		/* Vstart for dens_cbp_sine.n */
  if (notinit(dvrevc))  dvrevc = -0.065;	/* Vrev for hbp2, dens_cbp_sinec.n */
  if (notinit(cadist))  cadist = 1.0e-3;	/* Ca chan density in axonal tips, dens_cbp_sine.n */

  if (notinit(ivplot)) ivplot = 0; 		/* make I/V plot */
  if (notinit(pscal)) pscal   = 50e-12;		/* plot scale */

   if (notinit(drmc))   drmc = drm;              /* user set default Rm for hbp2 */
   if (notinit(drmab)) drmab = drm;              /* user set default Rm for axon branches */
   if (notinit(driab)) driab = dri;              /* user set default Ri for axon branches */
   if (notinit(dria))  dria  = dri;              /* user set default Ri for axon branches */
   if (notinit(dcmab)) dcmab = dcm;              /* user set default Cm for axon branches */
   if (notinit(naax))   naax = 0;                /* user set Na chan density in axon */
   if (notinit(naab))   naab = 0;                /* user set Na chan density in axon branches */
   if (notinit(nahd))   nahd = 0e-3;             /* user set Na chan high density region in axon */
   if (notinit(ksoma)) ksoma = 0e-3; 		 /* Kdr chan density in soma, dens_cbp_sine.n */
   if (notinit(kr6))     kr6 = 0e-3; 		 /* Kdr chan density in R6, dens_cbp_sine.n */
   if (notinit(kr7))     kr7 = 0e-3; 		 /* Kdr chan density in R7, dens_cbp_sine.n */
   if (notinit(amna))   amna = 0e-3; 		 /* Na chan density in amh, dens_amh_sine.n */
   if (notinit(amk))     amk = 0e-3; 		 /* Kdr chan density in amh, dens_amh_sine.n */
   if (notinit(amrm))   amrm = 20e3; 		 /* Rm for amh, dens_amh_sine.n */
   if (notinit(amvrev))  amvrev = -0.05;	 /* Vrev for amh in dens_amh_sine.n */
   if (notinit(hbp1_k1o)) hbp1_k1o = 0.010; 	 /* Kdr chan voltage offset, chanparams_cbp_chirp */
   if (notinit(hbp2_k1o)) hbp2_k1o = 0.010; 	 /* Kdr chan voltage offset, chanparams_cbp_chirp */

   if (notinit(axon_base)) axon_base = 701; 	 /* recording point at base of axon where it branches */
   if (make_cones==0) 
	if (notinit(stimloc)) stimloc = 0; 	 /* trandsucer stimulus location (default: soma) */

  // dicafrac = 1;				/* remove Ca flux from ICa */

  // Set channel types in density file.
  // To change, check manual, and uncomment and modify these:  
  //
  //  _NA  = _NA5;		// NaV1.1 channel from Clancy & Kass (2004)
  //  _KDR = _K1;
  //  _KA  = _K3;
  //  _KH  = _K4;
  //  _KIR = _K5;
  //  _CA_L = _CA0;

  if (!notinit(disp_zmax) || !notinit(disp_zmin)) {

  //  display_z(disp_zmax=-5, disp_zmin=-15);          /* display sublamina a */
  //  display_z(disp_zmax=-15, disp_zmin=-50);         /* display sublamina b */
      display_z(disp_zmax, disp_zmin);                 /* limit z display range */
  }
}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{
   int i, cn, maxnum;

    // if (!notinit(kdr_cond))    celdens[hbp1][0][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (hbp1,SCOND,1, 25e-10);
    //} 
    if (make_hbp1>0) {
       if (amarr==0 && n_hbp1==2) {
          ndens[hbp1][cn=1] = 0;              // set cn 1 to use hbp1_densfile
          ndens[hbp1][cn=2] = 0;              // set cn 2 to use hbp1_densfile
      } else {
        maxnum = getn(hbp1,MAXNUM);
        ndens[hbp1][cn=1] = 0;              // set cn 1 to use hbp1_densfile
        for (cn=2; cn<maxnum; cn++) {
            ndens[hbp1][cn] = 1;            // set other cns to use hbp1_densfile2
       }
      }
    } 
    else if (make_hbp2>0) {
        maxnum = getn(hbp2,MAXNUM);
        ndens[hbp2][cn=1] = 0;              // set cn 1 to use hbp2_densfile
        for (cn=2; cn<maxnum; cn++) {
            ndens[hbp2][cn] = 1;            // set other cns to use hbp2_densfile2
       }
    }
    else if (make_hbp1>0) {
        maxnum = getn(hbp1,MAXNUM);
        ndens[hbp1][cn=1] = 0;              // set cn 1 to use hbp1_densfile
        for (cn=2; cn<maxnum; cn++) {
            ndens[hbp1][cn] = 1;          // set other cns to use hbp1_densfile2
       }
    }
}

/*------------------------------------------------------*/

void addlabels(void)
{ 
     int i, c, ncells=0;

   if (ct==hbp1) {
     c = ct;
     for (i=1; i<=nhbp1; i++) {
       if (amarr > 0) { 
             label(nde(c,i,dendn_node(c,presyn1)),red,"S1");
	     continue;
       };
       if (no_tip_labels) continue;
       if (!notinit(tip1))    {label(nde(c,i,dendn_node(c,tip1)),red,"Tip1");}
       if (!notinit(tip2))    {label(nde(c,i,dendn_node(c,tip2)),red,"Tip2");}
       if (!notinit(tip3))    {label(nde(c,i,dendn_node(c,tip3)),red,"Tip3");}
       if (!notinit(tip4))    {label(nde(c,i,dendn_node(c,tip4)),red,"Tip4");}
       if (!notinit(tip5))    {label(nde(c,i,dendn_node(c,tip5)),red,"Tip5");}
       if (!notinit(tip6))    {label(nde(c,i,dendn_node(c,tip6)),red,"Tip6");}
       if (!notinit(tip7))    {label(nde(c,i,dendn_node(c,tip7)),red,"Tip7");}
       if (!notinit(tip8))    {label(nde(c,i,dendn_node(c,tip8)),red,"Tip8");}
       if (!notinit(tip9))    {label(nde(c,i,dendn_node(c,tip9)),red,"Tip9");}
       if (!notinit(tip10))   {label(nde(c,i,dendn_node(c,tip10)),red,"Tip10");}
       if (!notinit(tip11))   {label(nde(c,i,dendn_node(c,tip11)),red,"Tip11");}
       if (!notinit(tip12))   {label(nde(c,i,dendn_node(c,tip12)),red,"Tip12");}
       if (!notinit(tip13))   {label(nde(c,i,dendn_node(c,tip13)),red,"Tip13");}
       if (!notinit(tip14))   {label(nde(c,i,dendn_node(c,tip14)),red,"Tip14");}
       if (!notinit(tip15))   {label(nde(c,i,dendn_node(c,tip15)),red,"Tip15");}
       if (!notinit(tip16))   {label(nde(c,i,dendn_node(c,tip16)),red,"Tip16");}
       if (!notinit(tip17))   {label(nde(c,i,dendn_node(c,tip17)),red,"Tip17");}
       if (!notinit(tip18))   {label(nde(c,i,dendn_node(c,tip18)),red,"Tip18");}
       if (!notinit(tip19))   {label(nde(c,i,dendn_node(c,tip19)),red,"Tip19");}
       if (!notinit(tip20))   {label(nde(c,i,dendn_node(c,tip20)),red,"Tip20");}
       if (!notinit(tip21))   {label(nde(c,i,dendn_node(c,tip21)),red,"Tip21");}
       if (!notinit(tip22))   {label(nde(c,i,dendn_node(c,tip22)),red,"Tip22");}
       if (!notinit(tip23))   {label(nde(c,i,dendn_node(c,tip23)),red,"Tip23");}
       if (!notinit(tip24))   {label(nde(c,i,dendn_node(c,tip24)),red,"Tip24");}
       if (!notinit(tip25))   {label(nde(c,i,dendn_node(c,tip25)),red,"Tip25");}
       if (!notinit(tip26))   {label(nde(c,i,dendn_node(c,tip26)),red,"Tip26");}
       if (!notinit(tip27))   {label(nde(c,i,dendn_node(c,tip27)),red,"Tip27");}
       if (!notinit(tip28))   {label(nde(c,i,dendn_node(c,tip28)),red,"Tip28");}
       if (!notinit(tip29))   {label(nde(c,i,dendn_node(c,tip29)),red,"Tip29");}
       if (!notinit(tip30))   {label(nde(c,i,dendn_node(c,tip30)),red,"Tip30");}
     }
   }

   if (ct==hbp2) {
     c = ct;
     for (i=1; i<=nhbp2; i++) {
       if (amarr > 0) {
             label(nde(c,i,dendn_node(c,presyn1)),red,"S1");
             continue;
       };
       if (!notinit(tip1))    {label(nde(c,i,dendn_node(c,tip1)),red,"Tip1");}
       if (!notinit(tip2))    {label(nde(c,i,dendn_node(c,tip2)),red,"Tip2");}
       if (!notinit(tip3))    {label(nde(c,i,dendn_node(c,tip3)),red,"Tip3");}
       if (!notinit(tip4))    {label(nde(c,i,dendn_node(c,tip4)),red,"Tip4");}
       if (!notinit(tip5))    {label(nde(c,i,dendn_node(c,tip5)),red,"Tip5");}
       if (!notinit(tip6))    {label(nde(c,i,dendn_node(c,tip6)),red,"Tip6");}
       if (!notinit(tip7))    {label(nde(c,i,dendn_node(c,tip7)),red,"Tip7");}
       if (!notinit(tip8))    {label(nde(c,i,dendn_node(c,tip8)),red,"Tip8");}
       if (!notinit(tip9))    {label(nde(c,i,dendn_node(c,tip9)),red,"Tip9");}
       if (!notinit(tip10))   {label(nde(c,i,dendn_node(c,tip10)),red,"Tip10");}
       if (!notinit(tip11))   {label(nde(c,i,dendn_node(c,tip11)),red,"Tip11");}
     }
   }

   if (ct2==hbp2) {
     c = ct2;
     for (i=1; i<=nhbp2; i++) {
       if (amarr > 0) {
             label(nde(c,i,dendn_node(c,presyn2)),yellow,"S2");
             continue;
       };
       if (!notinit(tip1))    {label(nde(c,i,dendn_node(c,tip1)),blue,"Tip1");}
       if (!notinit(tip2))    {label(nde(c,i,dendn_node(c,tip2)),blue,"Tip2");}
       if (!notinit(tip3))    {label(nde(c,i,dendn_node(c,tip3)),blue,"Tip3");}
       if (!notinit(tip4))    {label(nde(c,i,dendn_node(c,tip4)),blue,"Tip4");}
       if (!notinit(tip5))    {label(nde(c,i,dendn_node(c,tip5)),blue,"Tip5");}
       if (!notinit(tip6))    {label(nde(c,i,dendn_node(c,tip6)),blue,"Tip6");}
       if (!notinit(tip7))    {label(nde(c,i,dendn_node(c,tip7)),blue,"Tip7");}
       if (!notinit(tip8))    {label(nde(c,i,dendn_node(c,tip8)),blue,"Tip8");}
       if (!notinit(tip9))    {label(nde(c,i,dendn_node(c,tip9)),blue,"Tip9");}
       if (!notinit(tip10))   {label(nde(c,i,dendn_node(c,tip10)),blue,"Tip10");}
       if (!notinit(tip11))   {label(nde(c,i,dendn_node(c,tip11)),blue,"Tip11");}
     }
   }

 if (ct>=0) {
   for (i=1; i<=nhbp1; i++) {
     if (no_recpnts) continue;
     if (!notinit(axon_base)) {label(nde(c,i,dendn_node(c,axon_base)),blue,"AB");}
     if (!notinit(stimloc)) {label(nde(ct,i,dendn_node(ct,stimloc)),yellow,"stim");}
     if (!notinit(vcloc))   {label(nde(ct,i,dendn_node(ct,vcloc)),blue,"VC");}
     if (!notinit(recpnt1)) {label(nde(ct,i,dendn_node(ct,recpnt1)),blue,"R1");}
     if (!notinit(recpnt2)) {label(nde(ct,i,dendn_node(ct,recpnt2)),blue,"R2");}
     if (!notinit(recpnt3)) {label(nde(ct,i,dendn_node(ct,recpnt3)),blue,"R3");}
     if (!notinit(recpnt4)) {label(nde(ct,i,dendn_node(ct,recpnt4)),blue,"R4");}
     if (!notinit(recpnt5)) {label(nde(ct,i,dendn_node(ct,recpnt5)),blue,"R5");}
     if (!notinit(recpnt6)) {label(nde(ct,i,dendn_node(ct,recpnt6)),blue,"R6");}
     if (!notinit(recpnt7)) {label(nde(ct,i,dendn_node(ct,recpnt7)),blue,"R7");}
     if (!notinit(recpnt8)) {label(nde(ct,i,dendn_node(ct,recpnt8)),blue,"R8");}
     if (!notinit(recpnt9)) {label(nde(ct,i,dendn_node(ct,recpnt9)),blue,"R9");}
     //label(ndn(ct,1,dendn_node(ct,0)),blue); 
     //label(ndn(ct,1,0),red);
   }
 }

  vna = 0.04;
}

/*------------------------------------------------------*/

void addcells() 
{
  if (make_gca) make_gc_comps(ct,1,R8,gca);			/* make postsynaptic compartments */
}

/*------------------------------------------------------*/


void runonexit (void)
{
  if (savefile[0]!=0)  unlink (savefile);
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

/*------------------------------------------------------*/

void runexpt(void) {

    int cn, i, n, plnum, electrode_node;
    int colr, midcone;
    int transdloc;
    double dst, t, fmax,fmin;
    double rmin, rmax, plsize;

    double xoffset = 0, yoffset = 0;
    double x,y;
    double Vmin, Vmax;
    double Imin, Imax;
    double cmin, cmax;
    double ipulse;
    double time2=0;
    photorec *p;
    node *npnt;
    char plabel[20] = {0};

  timinc = 1e-5;
  ploti = 1e-3;
  crit = 1e-10;

if (!make_cones) {
  for(npnt=nodepnt; npnt=foreach(npnt,ct,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {  	// add transducer to dendrite  
        if ((transdloc=dendn_node(ct,stimloc)) < 0) transdloc = soma;
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
  if (ct2>0)
   for(npnt=nodepnt; npnt=foreach(npnt,ct2,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {  	// add transducer to dendrite  
        if ((transdloc=dendn_node(ct2,stimloc)) < 0) transdloc = soma;
	if (itransducer) {
	   if ((p=(photorec*)make_itransducer(ndn(ct2,cn,transdloc))) != NULL) {
             p->xpos=npnt->xloc;                                                                     
             p->ypos=npnt->yloc;
	   }
	} else {
           if ((p=(photorec*)make_transducer(ndn(ct2,cn,transdloc))) != NULL) {
             p->xpos=npnt->xloc;                                                                     
             p->ypos=npnt->yloc;
	   }
	}
	if (ninfo >= 3) fprintf(stderr,"# ct2 %s cn %-2d transducer node %d\n",cname[ct2],cn, transdloc);
  }
}

  cn = 1;
  electrode_node = 5000;
  setonplot(onplot);

  if (notinit(temp_freq))     temp_freq = 10;
  // if (notinit(stimtime))       stimtime = (temp_freq>5 ? 0.05 : 0.1);
  if (notinit(stimtime))       stimtime = 0.1;
  if (notinit(fstart))           fstart = 1.0;
  if (notinit(fincr))             fincr = 0.005;
  if (notinit(c1))                   c1 = 0.01;
  if (notinit(c2))                   c2 = 1;
  if (notinit(vcontrast))     vcontrast = 0;
  if (notinit(spotdia))         spotdia = 300;
  if (notinit(spotloc))         spotloc = 0;
  if (notinit(vcaxon))           vcaxon = -0.06;

  if (notinit(prestimdur))   prestimdur = 0;
  if (notinit(spotdur))      spotdur    = 0;
  if (notinit(stimdur))      stimdur    = 5/temp_freq;
  if (notinit(poststimdur)) poststimdur = 0.05;

  if (notinit(itail)) itail  = 0e-12;
  if (notinit(dci))    dci   = 0e-12;
  if (notinit(ipre))   ipre  = dci;

  if (notinit(istart))      istart =   40e-12;
  if (notinit(istop))        istop =   120e-12;
  if (notinit(istep))        istep =   10e-12;

  if (itransducer==1) {
    if (notinit(minten))           minten = 10e-12;
    if (notinit(cmult))             cmult = 3e-12;
  } else {
    if (notinit(minten))           minten = -0.042;
    if (notinit(cmult))             cmult = 0.0125;
  }

  if (notinit(scontrast)) {
     if (cone_type > 1) scontrast = -1.0;
     else	        scontrast =  1.0;
  }
  if (notinit(elnode))     elnode  = soma;
  if (notinit(elec_rs))  elec_rs   = 400e6;
  if (notinit(elec_cap)) elec_cap  = 1e-14;

  if (elnode==electrode_node) {
      make_electrode  (nd(ct,cn,elnode), ndn(ct,cn,soma), elec_rs, elec_cap);
  }
	
  // midcbp  = findmid(ct,0,0);
   midcone  = findmid(xcone,0,0);

  if (ncones>0) {
      if (cone_type==1) {		// cone phototransduction
       	 plot_l_nod(xcone, midcone, soma, Vmin = 0, Vmax = 10000, colr=white, "", 22, 0.5);
	 plot_v_nod(xcone, midcone, soma, Vmin = -0.05, Vmax = -0.03, colr=blue, "", 21, 1.0);
      } else 
      if (cone_type==2) { 		// cone transducer
	 plot_v_nod(xcone, midcone, soma, Vmin = -0.05, Vmax = -0.03, colr=blue, "", 21, 0.5);
      } else
      if (cone_type==3) { 		// cone itransducer
	 plot_i_nod(xcone, midcone, soma, Vmin = -1e-10, Vmax = 1e-10, colr=blue, "", 21, 0.5);
      }
      plot_synrate(findsynloc(xcone,midcone,0,0),  fmin=0,fmax=400,  cyan,  20,"",0.5);
      // plot_syncond(findsynloc(xcone,midcone,0,0),  fmin=0,fmax=1e-10,  magenta,  20,"",0.5);
      if (vnoise > 0) 
	 plot_synves(findsynloc(xcone,midcone,0,0),   fmin=0,fmax=getn(xcone,STRCONC1), brown,  20,"",0.5);
      // plot_synrate_out(xcone,midcone,   rmin=0,rmax=400, fmax=1e-4, brown,"",0.5);
  }

  Vmin = -0.06; Vmax = -0.02;
  if (make_cones > 0) {
    plot_v_lnod(ct, cn, soma, Vmin,    Vmax,    colr=blue,   "", 3, 1.0);
    if (disp_c2) plot_v_lnod(ct2,findcell(ct2,0,0), soma,   Vmin,    Vmax,    colr=green,"", 3, 1.0);
  } else {
     if (itransducer) plot_i_nod(ct,cn,transdloc,Imin=minten-cmult, Imax=minten+cmult, colr=blue,  "Stim",  4,0.2);
     if (!notinit(stimloc)) plot_v_lnod(ct,cn,stimloc,Vmin, Vmax, colr=blue,  "Stim",  3,1.0);
     if (disp_c2) if (!notinit(stimloc)) plot_v_lnod(ct2,1,stimloc,Vmin, Vmax, colr=green,  "c2",  3,1.0);
  }
  if (!notinit(axon_base))plot_v_lnod(ct,cn,axon_base,Vmin, Vmax, colr=brown, "Axbase",3,1.0);
  if (!notinit(vcloc))   plot_v_lnod(ct, cn,vcloc,    Vmin, Vmax, colr=red,   "Cloc",  3,1.0);
  if (presyn1>0)         plot_v_lnod(ct, cn, presyn1, Vmin, Vmax, colr=red,   "psyn1", 3,1.0);
  if (!notinit(recpnt1)) plot_v_lnod(ct, cn, recpnt1, Vmin, Vmax, colr=magenta, "rp1", 3, 1.0);
  if (!notinit(recpnt2)) plot_v_lnod(ct, cn, recpnt2, Vmin, Vmax, colr=white, "rp2", 3, 1.0);
  if (!notinit(recpnt3)) plot_v_lnod(ct, cn, recpnt3, Vmin, Vmax, colr=yellow,"rp3", 3, 1.0);
  if (!notinit(recpnt4)) plot_v_lnod(ct, cn, recpnt4, Vmin, Vmax, colr=cyan,  "rp4", 3, 1.0);
  if (!notinit(recpnt5)) plot_v_lnod(ct, cn, recpnt5, Vmin, Vmax, colr=ltblue,"rp5", 3, 1.0);
  //if (!notinit(tip2))  plot_v_nod(ct, cn, dendn_node(ct,tip2),   Vmin, Vmax, colr=green,"", 3, 1.0);
  //if (!notinit(tip3))  plot_v_nod(ct, cn, dendn_node(ct,tip3),   Vmin, Vmax, colr=white,"", 3, 1.0);
  if (presyn2>0 && disp_c2) plot_v_lnod(ct2,cn, presyn2, Vmin, Vmax, colr=ltmag,"psyn2", 3, 1.0);

   // plot_chan(ct2, cn, 55, NA, 5, I, 5e-12, -10e-12);  // plot Na5 current in R6
   // plot_param("I chan",2,1.0);

  // plot_i_nod(ct, cn, soma, Imin = -3e-11, Imax = 3e-11, colr=magenta, "", 0, 0.5); 
  // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 

  //if (!notinit(tip1)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip1), 1, 10e-6, colr=red,   "Cabufb tip1", 6,0.5);
  //if (!notinit(tip2)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip2), 1, 10e-6, colr=green, "Cabufb tip2", 6,0.5);
  //if (!notinit(tip3)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip3), 1, 10e-6, colr=white, "Cabufb tip3", 6,0.5);
  if (presyn1>0) {
	plot_ca_nod(ct, cn, dendn_node(ct,presyn1),     1, 1.0e-6, cyan,       "Ca psyn1     ",    6,0.5); 
	// plot_ca_nod(ct, cn, dendn_node(ct,presyn1),     1, 1.0e-6, red,        "Ca psyn1filt ",    6,0.5); 
        // plotfilt(2,make_filt(0.01,0.2));
  	// plot_cabufb_nod(ct, cn, dendn_node(ct,presyn1), 1, 5.0e-6, colr=red,   "Cabufb psyn1", 6,0.5);

  }
  if (presyn2>0 && disp_c2) {
        if (n_hbp2==1) cn2 = 1; 
        else cn2 = 2;
        plot_ca_nod(hbp2, cn2, dendn_node(ct2,presyn2),  1, 1.0e-6, white,   "Ca psyn2     ",    6,0.5); 
        // plot_ca_nod(hbp2, cn2, dendn_node(ct2,presyn2),  1, 1.0e-6, brown,   "Ca psyn2filt ",    6,0.5); 
        //   plotfilt(2,make_filt(0.01,0.2));
  }

    plot_synrate(findsyn(hbp1,1,dendn_node(hbp1,presyn1),amh),   fmin=0,fmax=400,        red,  8,"psyn1",0.5);
    plot_synrate(findsyn(hbp2,1,dendn_node(hbp2,presyn2),amh),   fmin=0,fmax=400,      green,  8,"psyn2",0.5);
 if (vnoise > 0) 
    plot_synves(findsyn(hbp1,1,dendn_node(hbp1,presyn1),amh),    fmin=0,fmax=1e-4,     brown,  8,"psyn1",0.5);

    plot_synrate(findsyn(hbp1,1,dendn_node(hbp1,recpnt1),gca),   fmin=0,fmax=400,     magenta,8,"rp1",0.5);
    plot_synrate(findsyn(hbp1,1,dendn_node(hbp1,recpnt2),gca),   fmin=0,fmax=400,     white,8,"rp2",0.5);

    plot_syncond(findsyn(hbp1,1,dendn_node(hbp1,presyn1),amh),   cmin=0,cmax=100e-12,    red, 10,"psyn1",0.5);
    plot_syncond(findsyn(hbp2,1,dendn_node(hbp2,presyn2),amh),   cmin=0,cmax=100e-12,  green, 10,"psyn2",0.5);
    
    Vmin = -0.08; Vmax = 0.00; 
    plot_v_nod(amh,  1, findsynlocn(amh,1,0,0),    Vmin, Vmax, colr=red,     "", 12, 0.5);
    plot_v_nod(amh,  1, findnodloc (amh,1,100,0),  Vmin, Vmax, colr=blue,    "", 12, 0.5);
    plot_v_nod(amh,  1, findnodloc (amh,1,-100,0), Vmin, Vmax, colr=magenta, "", 12, 0.5);
    plot_v_nod(amh2, 1, findsynlocn(amh2,1,0,0),   Vmin, Vmax, colr=green,   "", 12, 0.5);
    // plot_v_nod(amh2, 1, 18,  Vmin, Vmax, colr=brown,   "", 12, 0.5);
    // plot_v_nod(amh2, 1, 20,  Vmin, Vmax, colr=yellow,   "", 12, 0.5);
    // plot_v_nod(amh, 1, findnodloc(amh,1,100,0), Vmin, Vmax, colr=white, "", 12, 0.5);
   
    plot_syncond(findsynloc(amh,1,0,0),    cmin=0,cmax=500e-12, red,     14,"",0.5);
    plot_syncond(findsynloc(amh2,1,0,0),   cmin=0,cmax=500e-12, green,   14,"",0.5);
    plot_syncond(findsyn(amh2,1,-1,amh),cmin=0,cmax=500e-12, blue,   14,"amh2_amh",0.5);

   if (disp) {                        // display the stimulus
        double t, dscale, starttime, disp_end;

      stim_backgr(minten);
      time2 = spot_sine(spotdia, x=0, y=0, temp_freq, cmult, scontrast, stimtime, stimdur);
      //time2 = spot_chirp(spotdia, 0, 0, fstart, fincr, cmult, scontrast, stimtime,stimdur);
      //time2 += spot_vcontrast(spotdia, 0, 0, temp_freq, 0, cmult, 0.005, 0.8,  stimtime,stimdur);

      display_size(500);
      disp_end = time2+0.05;
      for (starttime=simtime,t=stimtime; t<disp_end; starttime = t, t+= 0.001) {
           display_stim(starttime, t, dscale=4, -0.035, -0.055);
           //display_stim(t, dscale=4, -0.035, -0.045); 
           //display_stim(0+t, dscale=4); 
           simwait(0.10);
      }
      return;
   }

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 0.5;
  simtime = 0 - predur;
  if (vcontrast>0) {endexp = stimtime+2*stimdur+2*spotdur+2*prestimdur+poststimdur;}
  else             {endexp = stimtime+1*stimdur+2*spotdur+1*prestimdur+poststimdur;};
  stim_backgr(minten);

  //fprintf (stderr,"%d %g %g %g\n", vcloc, vcaxon, stimtime, endexp);

  if (make_cones) synaptau = 0.001;

  if (!notinit(vcloc)) vclamp (ndn(ct,cn,dendn_node(ct,vcloc)), vcaxon, stimtime, endexp);

  // cclamp(ndn(ct,cn,soma), ipre, simtime, predur+stimtime);
  step (predur);

  if (make_cones) synaptau = 1.0;

  if (dci>0) {
     cclamp (ndn(ct,cn,elnode),dci,stimtime,2*stimdur);
  }

  spot_sine(spotdia, x= spotloc, y=0, temp_freq, cmult, scontrast, stimtime, stimdur);

  /*
  stim_spot (spotdia, x=0, y=0,  cmult*scontrast, stimtime, spotdur);
  stim_spot (spotdia, x=0, y=0, -cmult*scontrast, stimtime+spotdur, spotdur);
  time2 = spot_chirp (spotdia, x=0, y=0, fstart, fincr, cmult, scontrast, stimtime+2*spotdur+prestimdur, stimdur);
  if (vcontrast>0) spot_vcontrast (spotdia, x=0, y=0, temp_freq, cmult, c1, c2, time2+prestimdur, stimdur);
  */

  step (endexp);
}


