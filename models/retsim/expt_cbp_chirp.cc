/* Experiment cbp_chirp */
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
int presyn3;
int ampar;
int ampar2;
int dbp1_am;
int dbp2_am;

const char *cbptype;
const char *cbptype2;

double kdr_cond;
double dbp1_gtau;
double dbp2_gtau;
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
double dbp1_mp;

double drmab;
double driab;
double dria;
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
double cadist;
double cadist2;
double catdist;
double catsoma;
double cone_soma_z;
double dbp1_soma_z;
double dbp2_soma_z;
double spacing;
double dtreedia;
double c2rot;			// z rot for 2nd cell of stereo pair (sets disparity for side view)
double c2yrot;			// y rot for 2nd cell of stereo pair (sets disparity for bottom view)

double varicos_dia;

double predur;
double prestimdur;
double spotdur;
double spotduri;
double stimdur;
double tailcurdur;
double poststimdur;
double stimtime;
double temp_freq;
double spotdia;

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
  setptr("dbp1_am",	&dbp1_am);
  setptr("dbp2_am",	&dbp2_am);
  setptr("no_tip_labels",&no_tip_labels);
  setptr("no_rpnt_labels",&no_rpnt_labels);
  setptr("disp_c2",     &disp_c2);
  setptr("nspots",      &nspots);
  setptr("stimtype",    &stimtype);

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
  setptr("presyn3",	&presyn3);
  setptr("itransducer",	&itransducer);

  setptr("drmab",       &drmab);
  setptr("driab",       &driab);
  setptr("dria",        &dria);
  setptr("dvreva",      &dvreva);
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
  setptr("kasoma",	&kasoma);
  setptr("kar6",	&kar6);
  setptr("kar7",	&kar7);
  setptr("dbp1_k1o",	&dbp1_k1o);
  setptr("khsoma",	&khsoma);
  setptr("cadist",	&cadist);
  setptr("cadist2",	&cadist2);
  setptr("catdist",	&catdist);
  setptr("catsoma",	&catsoma);
  setptr("kdr_cond",    &kdr_cond);
  setptr("g_am_dbp1",   &g_am_dbp1);
  setptr("g_am_dbp2",   &g_am_dbp2);
  setptr("g_am2_dbp1",  &g_am2_dbp1);
  setptr("g_am2_dbp2",  &g_am2_dbp2);
  setptr("g_dbp1_am",   &g_dbp1_am);
  setptr("g_dbp2_am",   &g_dbp2_am);
  setptr("g_dbp1_am2",  &g_dbp1_am2);
  setptr("g_dbp2_am2",  &g_dbp2_am2);
  setptr("elec_rs",	&elec_rs);
  setptr("elec_cap",	&elec_cap);
  setptr("dbp1_mp",	&dbp1_mp);

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
  setptr("ampar",       &ampar);
  setptr("ampar2",      &ampar2);

  setptr("cone_soma_z",	&cone_soma_z);
  setptr("dbp1_soma_z",	&dbp1_soma_z);
  setptr("dbp2_soma_z",	&dbp2_soma_z);
  setptr("spacing",     &spacing);
  setptr("dtreedia",    &dtreedia);
  setptr("c2rot",       &c2rot);
  setptr("c2yrot",      &c2yrot);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("spotdur",	&spotdur);
  setptr("spotduri",	&spotduri);
  setptr("stimdur",	&stimdur);
  setptr("tailcurdur",	&tailcurdur);
  setptr("poststimdur",	&poststimdur);
  setptr("spotdia",	&spotdia);

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

  nvalfile = "nval_cbp_chirp.n";
  chanparamsfile = "chanparams_cbp_chirp";

  // dbp1_file = "morph_bp";			// default set in retsim.cc 
  // dbp1_file = "morph_DB_111005_12_db4";	// morphology file in retsim.cc 
  // dbp1_file = "cell504_t4";			// morphology file in retsim.cc 
  // dbp1_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 
  // dbp2_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 
  dbp1_file = "cbp_0572_t5t";			// morphology file in retsim.cc 
  dbp2_file = "cbp_0572_t5t";			// morphology file in retsim.cc 
  hbp1_file = "morph_cbp_0572_t5r";		// morphology file in retsim.cc 

  dbp1_densfile  = "dens_cbp_chirp.n";		// cbplam, axdia, dria
  dbp1_densfile2 = "dens_cbp_chirpc.n";		// cbplam2, axdia, dria
  dbp2_densfile  = "dens_cbp2_chirp.n";
  dbp2_densfile2 = "dens_cbp2_chirp.n";
  hbp1_densfile  = "dens_cbp_chirp.n";
  hbp1_densfile2 = "dens_cbp_chirpc.n";

  g_am_dbp1  = 4e-10; 	/* am cond, feedback to dbp1, nval_cbp_chirp.n */
  g_am_dbp2  = 4e-10; 	/* am cond, feedback to dbp2, nval_cbp_chirp.n */
  g_am2_dbp1 = 4e-10; 	/* am2 cond, feedback to dbp1, nval_cbp_chirp.n */
  g_am2_dbp2 = 4e-10; 	/* am2 cond, feedback to dbp2, nval_cbp_chirp.n */
  g_dbp1_am  = 5e-10;	/* dbp1 cond, output to am, nval_cbp_chirp.n */
  g_dbp1_am2 = 5e-10; 	/* dbp1 cond, output to am2, nval_cbp_chirp.n */
  g_dbp2_am  = 5e-10; 	/* dbp2 cond, output to am, nval_cbp_chirp.n */
  g_dbp2_am2 = 5e-10; 	/* dbp2 cond, output to am2, nval_cbp_chirp.n */
  presyn1 = 0;  	/* presynaptic locus for synaptic output (label in nval.n) */
  presyn2 = 0;  	/* presynaptic locus for dbp2 synaptic output (label in nval.n) */
  presyn3 = 0;  	/* presynaptic locus for dbp1 2nd synaptic output (label in nval.n) */
  // recpnt1 = 0; 	 	/* recording point (label) */
  // recpnt2 = 0; 	 	/* recording point (label) */
  dbp1_gtau = 1; 	/* tau multiplier for GABA channel in dbp1 */
  dbp2_gtau = 1; 	/* tau multiplier for GABA channel in dbp2 */
  am_atau  = 1; 	/* tau multiplier for AMPA channel in am */
  am2_atau = 1; 	/* tau multiplier for AMPA channel in am2 */
  am_sdur  = 2; 	/* sdur for am->dbp1 synapse */
  am_sdur2 = 2; 	/* sdur for am->dbp2 synapse */
  am2_sdur = 2; 	/* sdur for am2->dbp1 synapse */
  am2_sdur2 = 2; 	/* sdur for am2->dbp2 synapse */
  am_sfall  = 0; 	/* fall time const for am->dbp1 synapse */
  am_sfall2 = 0; 	/* fall time const for am->dbp2 synapse */
  am2_sfall = 0; 	/* fall time const for am2->dbp1 synapse */
  am2_sfall2 = 0; 	/* fall time const for am2->dbp2 synapse */
  ampar  = xampa5; 	/* AMPA response for dbp1->am synapse */
  ampar2 = xampa5; 	/* AMPA response for dbp2->am synapse */
  dbp1_mp = 50;		/* max rrpool sustained rate */

  make_dbp1_am2 = 1;
  make_dbp2_am  = 1;
  make_dbp2_gca = 0;
  make_am2_dbp1 = 1;		/* lateral feedback from second amacrine type to the first bp type */
  make_am_gca  = 0;
  _CA_L = _CA1;			/* set type of L-type calcium channel for dens_cbp_chirp.n */
  _CA_T = _CA7;			/* set type of T-type calcium channel for dens_cbp_chirp.n */
  dscavg = 5e5;			/* Ca sensitivity for vesicle release */
}

/*------------------------------------------------------*/

int makrow (int ct, int ncells, double spacing, double rowtheta, double zrot, double xloc, double yloc, 
		int center, int i, int s)

/* make one row of cells along an orientation */
/* i counts number of cells in xarr[], yarr[]. */

#define CBPARR 100
#define R 5
{
  int c,j,k,n;
  int r[20];
  // static int r[] = {R,R-1,R+1,R-2,R+2,R-3,R+3,R-4,R+4,R-5,R+5,R-6,R+6,R-7,R+7,R-8,R+8};
  double *xarr, *yarr, *tharr;

  c = center;
  r[0] = c;
  for (j=1,n=1;j<18; n++,j+=2) {
       r[j]   = c-n;				// alternating either side of center
       r[j+1] = c+n;
  }
  // for (j=0; j<ncells; j++) {
  //      fprintf (stderr,"j %d r %d\n",j,r[j]);
  // }
  switch (ct) {
      case DBP1: xarr = dbp1xarr; yarr = dbp1yarr; tharr = dbp1tharr; break;
      case DBP2: xarr = dbp2xarr; yarr = dbp2yarr; tharr = dbp2tharr; break;
      case DBP3: xarr = dbp3xarr; yarr = dbp3yarr; tharr = dbp3tharr; break;
      case DBP4: xarr = dbp4xarr; yarr = dbp4yarr; tharr = dbp4tharr; break;
  }
  if (xarr==NULL) {
      xarr = (double *)emalloc(CBPARR*sizeof(double));
      yarr = (double *)emalloc(CBPARR*sizeof(double));
     tharr = (double *)emalloc(CBPARR*sizeof(double));
     switch (ct) {
       case DBP1: dbp1xarr = xarr; dbp1yarr = yarr; dbp1tharr = tharr; break;
       case DBP2: dbp2xarr = xarr; dbp2yarr = yarr; dbp2tharr = tharr; break;
       case DBP3: dbp3xarr = xarr; dbp3yarr = yarr; dbp3tharr = tharr; break;
       case DBP4: dbp4xarr = xarr; dbp4yarr = yarr; dbp4tharr = tharr; break;
     }
  }
  if (xarr!=NULL) {
     for (j=0;j<(ncells-s); i++,j++) {
         k = r[j+s];  
         // fprintf (stderr,"c %d i %d j %d k %d\n",c,i,j,k);
         xarr[i]  = -cosdeg(rowtheta) * ((c-k)*spacing) + xloc;
         yarr[i]  =  sindeg(rowtheta) * ((c-k)*spacing) + yloc;
         tharr[i] = zrot-rowtheta;
     }
  } else 
     fprintf (stderr,"expt_cbp_inh: makrow, can't allocate memory for %s\n",cname[ct]);

  return i;
}
#undef R

/* - - - - - - - - - - - - - - - - - - - */

int dbp_row (int ct, int ncells, double spacing, double rowtheta, double zrot, double xloc, double yloc, int center, int i, int s)
{
      double dist;		// radial distance from (xloc, yloc)

    return makrow (ct, ncells, spacing, rowtheta, zrot, xloc, yloc, center, i, s);
}

/* - - - - - - - - - - - - - - - - - - - */

int dbp_row (int ct, int ncells, double spacing, double rowtheta, double zrot, int center, int i, int s)
{

    return makrow (ct, ncells, spacing, rowtheta, zrot, 0, 0, center, i, s);
}

/* - - - - - - - - - - - - - - - - - - - */

int dbp_cent (int ct, int ncells, double theta, double zrot, int center, int i)
{
      int s;			// start
      double spacing;		// cell spacing along row (here, none)

    return makrow (ct, ncells, spacing=0, theta, zrot, 0, 0, center, i, s=0);
}
/* - - - - - - - - - - - - - - - - - - - */

int make_am_cell (int ct, double theta, double xloc, double yloc, int i, double somadist)
{
   int x;
   double *xarr,*yarr,*tharr;
#define AMARR  20

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
  if (xarr!=NULL) {
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

/*------------------------------------------------------*/
   
void setparams(void)
{
#define AMARR  20

     int c,i,s, dist, cent;
     double theta;						/* orientation of amacrine cell to follow */
     double dbp1_zrot=0, dbp2_zrot=0; 				/* rotation of cell at amacrine cell */

  if (!notinit(cbptype)) {			/* primary cbp type */
       ct = find_ct(cbptype);
  } else {
       ct = dbp1;
  }

  if (!notinit(cbptype2)) {			/* secondary cbp type */
       ct2 = find_ct(cbptype2);
  } else {
       ct2 = dbp2;
  }

  if (!notinit(cone_soma_z)) setn(xcone,SOMAZ,cone_soma_z);
  if (!notinit(dbp1_soma_z)) setn(dbp1,SOMAZ,dbp1_soma_z);
  if (!notinit(dbp2_soma_z)) setn(dbp2,SOMAZ,dbp2_soma_z);
  if (notinit(spacing)) spacing = 40;

  make_ct(ct);
  if (ct2>0) make_ct(ct2);
  // make_ct(gca);

  SOMA = R_3;			/* defined in retsim_var.cc, retsim.h */

  if (notinit(dbp1_am)) dbp1_am = 11;
  if (notinit(amarr)) amarr = 0;

  cent = int((dbp1_am+1)/2);				    /* the central dbp */
  dbp_cent (dbp1, 1,      theta=  0,   dbp1_zrot,   cent, i=0);

  switch (amarr) {						/* set up dbp1s */

     case 0: 
	if      (n_dbp1<=1) { } 
	else if (n_dbp1==2) {
	  dbp1xarr[0] = 0;
	  dbp1yarr[0] = 0;
	  dbp1tharr[0] = 0;
	  dbp1xarr[1] = spacing;
	  dbp1yarr[1] = 0;
	  if (notinit(c2rot)) c2rot = 5;
	  dbp1tharr[1] = c2rot;
	  if (c2yrot!=0) {		// if stereo pair seen from bottom
             dbp1ytharr  = (double *)emalloc(CBPARR*sizeof(double));
	     dbp1ytharr[0] = 0;
	     dbp1ytharr[1] = -c2yrot;
	  }
	}
	//remove_nconns = 0;
	break;

     case 1:	
	 // dbp1_am = 6;
	 if (n_dbp1 < 0 || n_dbp1 > dbp1_am) n_dbp1 = dbp1_am;			/* limit to 6 dbp1s */
	 i = dbp_row (dbp1, dbp1_am, spacing, theta=  180, dbp1_zrot, cent,  i=1,s=0);
	break; 

     case 2:
	break;

     case 11:	
	  if (n_dbp1 < 0 || n_dbp1 > dbp1_am) n_dbp1 = dbp1_am;			/* limit to 11 dbp1s */
	  i = dbp_row (dbp1, n_dbp1, spacing, theta=  0,   dbp1_zrot, cent,  i=1,s=1);
	  break;
 
     case 21:

	  if (n_dbp1 < 0 || n_dbp1 > 2*dbp1_am) n_dbp1 = 2*dbp1_am - 1; 	  /* limit to 21 dbp1s */
	  // i = dbp_cent (dbp1, 1,        theta=  45, dbp1_zrot,  cent, i=0);	  /* set cell rot near both ams */
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=  90, dbp1_zrot,  cent, i=1,s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=   0, dbp1_zrot,  cent, i,  s=1);
	  setn(dbp1,MAXSDIST,15); 
	  break;
 
     case 22: 									/* use with dbp1, dbp2 */
	  if (n_dbp1 < 0 || n_dbp1 > dbp1_am) n_dbp1 = dbp1_am;			/* limit to 11 dbp1s */
	  i = dbp_cent (dbp1, 1, theta=  0,   dbp1_zrot+315,         cent, i=0);
	  // i = dbp_row (dbp1, dbp1_am, spacing, theta= 0, dbp1_zrot+55, cent, i=1,s=1);
	  setn(dbp1,MAXSDIST,25); 
	  break;

     case 23: 									/* use with dbp1, dbp2 */
          if (n_dbp1 >= 0) {			/* random array */
             dbp1xarrsiz = dbp1_am * spacing;
             dbp1yarrsiz = 30;
 	     dbp1_first_cent = 1;
          } 
	  break;

     case 3:
	  if (n_dbp1 < 0 || n_dbp1 > 3*dbp1_am-2) n_dbp1 = 3*dbp1_am-2;		/* limit to 31 dbp1s */
	  i = dbp_cent (dbp1, 1,             theta=  45, dbp1_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=  0,  dbp1_zrot, cent, i=1,s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=  60, dbp1_zrot, cent, i,  s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta= -60, dbp1_zrot, cent, i,  s=1);
	  setn(dbp1,MAXSDIST,14);
	  break;

     case 4:
	  if (n_dbp1 < 0 || n_dbp1 > 4*dbp1_am-3) n_dbp1 = 4*dbp1_am-3;		/* limit to 41 dbp1s */
	  i = dbp_cent (dbp1, 1,        theta=  45, dbp1_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=   0, dbp1_zrot, cent, i=1,s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=  45, dbp1_zrot, cent, i, s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta= -45, dbp1_zrot, cent, i, s=1);
	  i = dbp_row (dbp1, dbp1_am, spacing, theta=  90, dbp1_zrot, cent, i, s=1);
	  setn(dbp1,MAXSDIST,14.7);
	 break;

     default: break;
 
  }  /* switch (amarr) */

  setn(dbp1,DENS,pow(spacing*0.5,-2)*1e6);	/* convert dbp1 density from spacing to mm2 */
 
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (notinit(dbp2_am)) {
     if (!notinit(dbp1_am)) dbp2_am = dbp1_am;
     else dbp2_am = 11;
  }

  if (amarr > 0) {
     cent = int((dbp2_am+1)/2);					/* skip the central dbp (when needed) */
     // dbp_cent (dbp2, 1,      theta=  0,   dbp2_zrot,   cent, i=0);
  }

  if (n_dbp2 > 0) {
    switch (amarr) {						/* set up dbp2s */

     case 0:
	  n_dbp2 = 1;								/* just 1 dbp2 */
          dbp_cent (dbp2, 1, theta=  0,   dbp2_zrot,   cent, i=0);
	  // remove_nconns = 0;
	  break;

     case 1:
	  dbp2_am = 6;
	  if (n_dbp2 > dbp2_am+1) n_dbp2 = dbp2_am;				/* limit to 5 dbp2s */
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  180, dbp2_zrot, cent,  i=1,s=0);
	  break;

     case 11:
	  if (n_dbp2 > dbp2_am) n_dbp2 = dbp2_am;				/* limit to 11 dbp2s */
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  0,   dbp2_zrot, cent,  i=0,s=1);
	  break;

     case 21: 
	  if (n_dbp2 > 2*dbp2_am) n_dbp2 = 2*dbp2_am - 1; 			   /* limit to 21 dbp2s */
	  i = dbp_cent (dbp2, 1,            theta=  45, dbp2_zrot,  cent,  i=0); /* set cell rot near both ams */
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  90, dbp2_zrot,  cent,  i=1,s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=   0, dbp2_zrot,  cent,  i,  s=1);
	  setn(dbp2,MAXSDIST,15); 
	  break;

     case 22: 		 							/* use with dbp1, dbp2 */
	 if (n_dbp2 >= dbp2_am) n_dbp2 = dbp2_am;				/* limit to 10 dbp2s */
	 i = dbp_row (dbp2, dbp2_am, spacing, theta=  0,  dbp2_zrot, cent, i=0, s=1);
	 i = dbp_row (dbp2, dbp2_am, spacing, theta=  90, dbp2_zrot, cent, i,   s=1);
	  setn(dbp2,MAXSDIST,15); 
	 break;

     case 23:	
	  if (n_dbp2 >= 0) {				/* random array */
             dbp2xarrsiz = 30;
             dbp2yarrsiz = dbp2_am * spacing;
	     dbp2_first_cent = 0;
	     setn(dbp2,MAXSDIST,15); 
	  }
	 break;

     case 3:
	  if (n_dbp2 > 3*dbp2_am+1) n_dbp2 = 3*dbp2_am - 2;			   /* limit to 31 dbp2s */
	  i = dbp_cent (dbp2, 1,             theta=  45, dbp2_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  0,  dbp2_zrot, cent, i=1,s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  60, dbp2_zrot, cent, i,  s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta= -60, dbp2_zrot, cent, i,  s=1);
	  setn(dbp2,MAXSDIST,14);
	  break;

     case 4:
	  if (n_dbp2 > 4*dbp2_am+1) n_dbp2 = 4*dbp2_am - 3;			   /* limit to 41 dbp2s */
	  i = dbp_cent (dbp2, 1,             theta=  45, dbp2_zrot, cent, i=0);  /* set cell rot near both ams */
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=   0, dbp2_zrot, cent, i=1,s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  45, dbp2_zrot, cent, i,  s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta= -45, dbp2_zrot, cent, i,  s=1);
	  i = dbp_row (dbp2, dbp2_am, spacing, theta=  90, dbp2_zrot, cent, i,  s=1);
	  setn(dbp2,MAXSDIST,14.7);
          break;

    }       /* switch (amarr) */
  }       /* if (ndbp2 > 0) */

  setn(dbp2,DENS,pow(spacing*0.5,-2)*1e6);	/* convert dbp2 density from spacing to mm2 */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (!notinit(amarr)) {
        int i,j;
        double theta;
	double somadist;

   if (amarr > 0) make_ct(am);

   setn(ct,AXARBT,BRANCHED);		/* connect to bipolar cell branched axonal arbor */
   if (ct2>0) setn(ct2,AXARBT,BRANCHED);	/* connect to bipolar cell branched axonal arbor */

   // dtreedia = n_dbp1 * spacing;
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

   if (notinit(dtreedia)) dtreedia = dbp1_am * spacing * 2;
   if (dtreedia <= 0) dtreedia = 650;
   somadist = dtreedia * 0.3;

   if (ninfo > 2) ncfprintf (stderr,"# setparams(): am dtreedia %g somadist %g\n",dtreedia,somadist);

   setn(am, DTREEDIA,dtreedia);
   setn(am2,DTREEDIA,dtreedia);

   switch (amarr) {

     case 1:
	  if (n_dbp1 <= 6 && n_dbp2 <=6) setn(am,DTREEDIA,dtreedia);
          make_am_cell (am, theta=180, i=0, somadist);
	  n_am = 1;
          break;
   
    case 2: 
	  // if (n_dbp1 <= 6 && n_dbp2 <=6) setn(am,DTREEDIA,dtreedia);
          // make_am_cell (am, theta=180, i=0, somadist);
	  n_am = 0;
          break;

    case 11:
	  if (n_dbp1 <= 11 && n_dbp2 <=11) setn(am,DTREEDIA,dtreedia);
          make_am_cell (am, theta=180, i=0, somadist);
	  n_am = 1;
          break;

    case 21:
          i = make_am_cell (am, theta=180, i=0, somadist);
          i = make_am_cell (am, theta=90,  i,   somadist);
          n_am = 2;
          break;

    case 22:
          make_am_cell (am,  theta = 180, i=0, somadist);
          make_am_cell (am2, theta = -90,  i=0, somadist);
	  n_am  = 1;
	  n_am2 = 1;
          break;

    case 23:
          i = make_am_cell (am, theta=135, i=0, somadist);
          i = make_am_cell (am, theta=45,  i,   somadist);
          n_am = 2;
          break;

    case 3:
          i = make_am_cell (am, theta=0,   i=0, somadist);
          i = make_am_cell (am, theta=60,  i,   somadist);
          i = make_am_cell (am, theta=-60, i,   somadist);
          n_am = 3;
          break;
 
    case 4:
          i = make_am_cell (am, theta=0,   i=0, somadist);
          i = make_am_cell (am, theta=45,  i, somadist);
          i = make_am_cell (am, theta=-45, i, somadist);
          i = make_am_cell (am, theta=90,  i, somadist);
          n_am = 4;
          break;
   }  /* switch (amarr) */

  }  /* if (!notinit(amarr)) */

  if (n_amh > 0) {
	  make_ct(amh);
          setn(ct,AXARBT,BRANCHED);	/* connect to bipolar cell branched axonal arbor */
	  
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
  }

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
  if (notinit(no_rpnt_labels)) no_rpnt_labels = 0;	/* = 1 => display tip labels but not recording pts */
  if (notinit(disp_c2))       disp_c2 = 1;	/* = 1 => display plots for cell 2 */

  if (notinit(cbplam))   cbplam = 0.02; 	/* default complam for regions in density file */
  if (notinit(cbplam2)) cbplam2 = 0.2; 		/* default complam for regions in density file 2 (large comps) */
  // if (notinit(dispsize)) dispsize = 80; 	/* default display size */
  if (notinit(node_scale)) node_scale = -3.15;  /* 3: nodenum, 0.2 medium font, 0.05: small font */
  if (notinit(dbp1_nscale)) dbp1_nscale = -2.09;
  if (notinit(dbp2_nscale)) dbp2_nscale = -2.09;

  if (notinit(ddia))     ddia  = 1; 		/* multiplier for dendrite diameter */
  if (notinit(axdia))   axdia  = 1; 		/* multiplier for axon diameter */
  if (notinit(axdiap)) axdiap  = 1; 		/* multiplier for proximal axon diameter */
  if (notinit(dvrev))   dvrev  = -0.065;	/* Vrev for dens_dbp1.n */
  if (notinit(dvst))    dvst   = -0.04;		/* Vstart for dens_dbp1.n */
  if (notinit(dvreva))  dvreva = -0.07;		/* Vrev for axon terminal in dens_dbp1.n */
  if (notinit(cadist))  cadist = 1.0e-3;	/* Ca chan density in axonal tips, dens_cbp_chirp.n */
  if (notinit(cadist2)) cadist2 = 1.0e-3;	/* Ca chan density in axonal tips, dens_cbp2_chirp.n */
  if (notinit(catdist)) catdist = 0e-3;		/* Ca T-type chan density in axonal tips, dens_cbp_chirp.n */
  if (notinit(catsoma)) catsoma = 0e-3;		/* Ca T-type chan density in soma, dens_cbp_chirp.n */

  if (notinit(ivplot)) ivplot = 0; 		/* make I/V plot */
  if (notinit(pscal)) pscal   = 50e-12;		/* plot scale */

   if (notinit(drmab)) drmab = drm;              /* user set default Rm for axon branches */
   if (notinit(driab)) driab = dri;              /* user set default Ri for axon branches */
   if (notinit(dria))  dria  = dri;              /* user set default Ri for axon branches */
   if (notinit(dcmab)) dcmab = dcm;              /* user set default Cm for axon branches */
   if (notinit(naax))   naax = 0;                /* user set Na density in axon */
   if (notinit(naab))   naab = 0;                /* user set Na density in axon branches */
   if (notinit(nahd))   nahd = 0e-3;             /* user set Na high density region in axon */
   if (notinit(ksoma)) ksoma = 0e-3; 		 /* Kdr chan density in soma, dens_cbp_chirp.n */
   if (notinit(kr6))     kr6 = 0e-3; 		 /* Kdr chan density in R6, dens_cbp_chirp.n */
   if (notinit(kr7))     kr7 = 0e-3; 		 /* Kdr chan density in R7, dens_cbp_chirp.n */
   if (notinit(dbp1_k1o)) dbp1_k1o = 0.010; 	 /* Kdr chan voltage offset, chanparams_cbp_chirp */
   if (notinit(kasoma)) kasoma = 0e-3; 		 /* KA chan density in soma, dens_cbp_chirp.n */
   if (notinit(kar6))     kar6 = 0e-3; 		 /* KA chan density in R6, dens_cbp_chirp.n */
   if (notinit(kar7))     kar7 = 0e-3; 		 /* KA chan density in R7, dens_cbp_chirp.n */
   if (notinit(khsoma)) khsoma = 0e-3; 		 /* Kh chan density in soma, dens_cbp_chirp.n */

   if (notinit(axon_base)) axon_base = 701; 	 /* recording point at base of axon where it branches */
   if (make_cones==0) 
	if (notinit(stimloc)) stimloc = 0; 	 /* trandsucer stimulus location (default: soma) */

  // dicafrac = 1;				/* remove Ca flux from ICa */

  // Set channel types in density file.
  // To change, check manual, and uncomment and modify these:  
  //
  //  _NA  = _NA5;		// NaV1.1 channel from Clancy & Kass (2004)
  //  _KA  = _K3;
  //  _KH  = _K4;
  //  _KIR = _K5;

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

    // if (!notinit(kdr_cond))    celdens[dbp1][0][_KDR][R_SOMA] = kdr_cond;  
    // setsv (ct,SCOND,1, 0);
    // if (arrsiz==100) {
    //  setsv (dbp1,SCOND,1, 25e-10);
    //} 
    if (make_dbp1>0) {
       if (amarr==0 && n_dbp1==2) {
          ndens[dbp1][cn=1] = 0;              // set cn 1 to use dbp1_densfile
          ndens[dbp1][cn=2] = 0;              // set cn 2 to use dbp1_densfile
      } else {
        maxnum = getn(dbp1,MAXNUM);
        ndens[dbp1][cn=1] = 0;              // set cn 1 to use dbp1_densfile
        for (cn=2; cn<maxnum; cn++) {
            ndens[dbp1][cn] = 1;            // set other cns to use dbp1_densfile2
       }
      }
    } 
    else if (make_dbp2>0) {
        maxnum = getn(dbp2,MAXNUM);
        ndens[dbp2][cn=1] = 0;              // set cn 1 to use dbp2_densfile
        for (cn=2; cn<maxnum; cn++) {
            ndens[dbp2][cn] = 1;            // set other cns to use dbp2_densfile2
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

   if (ct==dbp1) {
     c = ct;
     for (i=1; i<=ndbp1; i++) {
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

   if (ct==dbp2) {
     c = ct;
     for (i=1; i<=ndbp2; i++) {
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

   if (ct2==dbp2) {
     c = ct2;
     for (i=1; i<=ndbp2; i++) {
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
   for (i=1; i<=ndbp1; i++) {
     if (no_rpnt_labels) continue;
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
  // if (make_gca) make_gc_comps(ct,1,R8,gca);			/* make postsynaptic compartments */
  if (n_dbp1 != 2) make_gc_comps(ct,1,R8,gca);		/* make postsynaptic compartments */
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
    int colr, midcone, npixels, nframes;
    int transdloc;
    double dst, t, fmax,fmin;
    double rmin, rmax, plsize, msize;
    double mask,orient;
    double xoffset = 0, yoffset = 0;
    double x,y;
    double Vmin, Vmax;
    double Imin, Imax;
    double cmin, cmax;
    double ipulse;
    double time2=0;
    double *rndarr = NULL;
    char *checkerboard_file;
    FILE *fchk;
    photorec *p;
    node *npnt;
    char plabel[20] = {0};

  timinc = 2e-5;
  ploti = 1e-3;
  crit = 1e-9;

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

  if (notinit(prestimdur))   prestimdur = 0.5;
  if (notinit(spotdur))      spotdur    = 0.5;
  if (notinit(spotduri))     spotduri   = 0.2;
  if (notinit(stimdur))      stimdur    = 2.0;
  if (notinit(poststimdur)) poststimdur = 0.05;
  if (notinit(nspots))           nspots = 3;

  if (notinit(stimtime))       stimtime = 0.1;
  if (notinit(temp_freq))     temp_freq = 10;
  if (notinit(fstart))           fstart = 1.0;
  if (notinit(fincr))             fincr = 0.005;
  if (notinit(c1))                   c1 = 0.01;
  if (notinit(c2))                   c2 = 1;
  if (notinit(vcontrast))     vcontrast = 0;
  if (notinit(spotdia))         spotdia = 300;
  if (notinit(vcaxon))           vcaxon = -0.06;

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

  // if (elnode==electrode_node) {
  //     make_electrode  (nd(ct,cn,elnode), ndn(ct,cn,soma), elec_rs, elec_cap);
  // }
	
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
      plot_synrate(findsynloc(xcone,midcone,0,0),  fmin=0,fmax=400,  cyan,  20,"",0.3);
      if (vnoise > 0) 
	 plot_synves(findsynloc(xcone,midcone,0,0),   fmin=0,fmax=getn(xcone,STRCONC1), brown,  20,"",0.3);
      // plot_synrate_out(xcone,midcone,   rmin=0,rmax=400, fmax=1e-4, brown,"",0.5);
  }

  Vmin = -0.06; Vmax = -0.02;
  if (make_cones > 0) {
    plot_v_lnod(ct, cn, soma, Vmin,  Vmax,  colr=blue,   "", 3, 1.0);
    // plot_v_lnod(ct2,findcell(ct2,0,0), soma,  Vmin,    Vmax,    colr=green,"", 3, 1.0);
    plot_v_lnod(ct2, 1, soma, Vmin,  Vmax,  colr=white,"", 3, 1.0);
    plot_v_lnod(ct2, 2, soma, Vmin,  Vmax,  colr=gray,"", 3, 1.0);
  } else {
     if (itransducer) plot_i_nod(ct,cn,transdloc,Imin=minten-cmult, Imax=minten+cmult, colr=blue,  "Stim",  4,0.2);
     if (!notinit(stimloc)) plot_v_lnod(ct,cn,stimloc,Vmin, Vmax, colr=blue,  "Stim",  3,1.0);
     if (disp_c2) if (!notinit(stimloc)) plot_v_lnod(ct2,1,stimloc,Vmin, Vmax, colr=green,  "c2",  3,1.0);
  }
  if (!notinit(axon_base))plot_v_lnod(ct,cn, axon_base,Vmin, Vmax, colr=brown,"Axbase",3, 1.0);
  if (!notinit(vcloc))   plot_v_lnod(ct, cn, vcloc,    Vmin, Vmax, colr=red,  "Cloc", 3, 1.0);
  if (presyn1>0)         plot_v_lnod(ct, cn, presyn1, Vmin, Vmax, colr=red,   "psyn1", 3, 1.0);
  if (presyn2>0)         plot_v_lnod(ct, cn, presyn2, Vmin, Vmax, colr=green, "psyn2", 3, 1.0);
  if (presyn3>0 && disp_c2) {
	  plot_v_lnod(ct2,1, presyn3, Vmin, Vmax, colr=magenta,  "psyn3", 3, 1.0);
	  plot_v_lnod(ct2,2, presyn3, Vmin, Vmax, colr=ltgreen,"psyn4", 3, 1.0);
  }
  if (!notinit(recpnt1)) plot_v_lnod(ct, cn, recpnt1, Vmin, Vmax, colr=blue,   "rp1", 3, 1.0);
  if (!notinit(recpnt2)) plot_v_lnod(ct, cn, recpnt2, Vmin, Vmax, colr=ltgreen,"rp2", 3, 1.0);
  if (!notinit(recpnt3)) plot_v_lnod(ct, cn, recpnt3, Vmin, Vmax, colr=ltred,  "rp3", 3, 1.0);
  if (!notinit(recpnt4)) plot_v_lnod(ct, cn, recpnt4, Vmin, Vmax, colr=ltcyan, "rp4", 3, 1.0);
  // if (!notinit(recpnt5)) plot_v_lnod(ct, cn, recpnt5, Vmin, Vmax, colr=white,  "rp5", 3, 1.0);

  //if (!notinit(tip2))  plot_v_nod(ct, cn, dendn_node(ct,tip2),   Vmin, Vmax, colr=green,"", 3, 1.0);
  //if (!notinit(tip3))  plot_v_nod(ct, cn, dendn_node(ct,tip3),   Vmin, Vmax, colr=white,"", 3, 1.0);

  // plot_i_nod(ct, cn, soma, Imin = -3e-11, Imax = 3e-11, colr=magenta, "", 0, 0.5); 
  // plot_ca_nod(ct, cn, axtrm, 10e-6, cyan, "Ca axterm", 1, 0.5); 

  //if (!notinit(tip1)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip1), 1, 10e-6, colr=red,   "Cabufb tip1", 6,0.5);
  //if (!notinit(tip2)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip2), 1, 10e-6, colr=green, "Cabufb tip2", 6,0.5);
  //if (!notinit(tip3)) plot_cabufb_nod(ct, cn, dendn_node(ct,tip3), 1, 10e-6, colr=white, "Cabufb tip3", 6,0.5);
  if (presyn1>0) {
	plot_ca_nod(ct, cn, dendn_node(ct,presyn1),      1, 1.0e-6, red,       "Ca psyn1     ",    6,0.5); 
	// plot_ca_nod(ct, cn, dendn_node(ct,presyn1),     1, 1.0e-6, cyan,        "Ca psyn1filt ",    6,0.5); 
        // plotfilt(2,make_filt(0.02,0.1));
  	// plot_cabufb_nod(ct, cn, dendn_node(ct,presyn1), 1, 1.0e-6, colr=red,   "Cabufb psyn1", 6,0.5);

  } 
  if (presyn2>0) {
        plot_ca_nod(ct, cn, dendn_node(ct,presyn2),     1, 1.0e-6, green,     "Ca psyn2     ",    6,0.5); 
        // plot_ca_nod(dbp1, cn, dendn_node(ct,presyn2),  1, 1.0e-6, magenta,   "Ca psyn2filt ",    6,0.5); 
        //   plotfilt(2,make_filt(0.05,0.2));
  }
  // if (presyn3>0 && disp_c2) {
  //       //if (n_dbp2==1) cn2 = 1; 
  //       //else cn2 = 2;
  //       cn2 = 1;
  //       plot_ca_nod(dbp2, cn2, dendn_node(ct2,presyn3),  1, 1.0e-6, magenta, "Ca psyn3     ",    6,0.5); 
  //       plot_ca_nod(dbp2, cn2, dendn_node(ct2,presyn3),  1, 1.0e-6, brown,   "Ca psyn3filt ",    6,0.5); 
  //         plotfilt(2,make_filt(0.05,0.2));
  // }
  if (recpnt1>0) {
	plot_ca_nod(ct, cn, dendn_node(ct,recpnt1),     1, 1.0e-6, blue,       "Ca rp1     ",    6,0.5); 
	plot_ca_nod(ct, cn, dendn_node(ct,recpnt1),     1, 1.0e-6, ltblue,     "Ca rp1filt ",    6,0.5); 
        plotfilt(2,make_filt(0.02,0.1));
  }
	plot_ca_nod(ct, cn, dendn_node(ct,recpnt2),     1, 1.0e-6, ltgreen,      "Ca rp2     ",    6,0.5); 
	plot_ca_nod(ct, cn, dendn_node(ct,recpnt3),     1, 1.0e-6, ltred,       "Ca rp3     ",    6,0.5); 
	plot_ca_nod(ct, cn, dendn_node(ct,recpnt4),     1, 1.0e-6, ltcyan,     "Ca rp4     ",    6,0.5); 

  plot_synrate(findsyn(ct,1,dendn_node(ct,presyn1),gca),   fmin=0,fmax=400,        red,  8,"psyn1",0.8);
  plot_synrate(findsyn(ct,1,dendn_node(ct,presyn2),gca),  fmin=0,fmax=400,      green,  8,"psyn2",0.8);
  // plot_synrate(findsyn(ct,1,dendn_node(ct,presyn1),am),   fmin=0,fmax=400,        red,  8,"psyn1",0.8);
  // plot_synrate(findsyn(ct,1,dendn_node(ct,presyn2),am2),  fmin=0,fmax=400,      green,  8,"psyn2",0.8);
  // plot_synrate(findsyn(dbp2,1,dendn_node(dbp2,presyn3),am),   fmin=0,fmax=400,   magenta,  8,"psyn3",0.8);
  // plot_synrate(findsyn(dbp2,2,dendn_node(dbp2,presyn3),am2),  fmin=0,fmax=400,      blue,   8,"psyn4",0.8);

  if (vnoise > 0) 
    plot_synves(findsyn(ct,1,dendn_node(dbp1,presyn1),am),    fmin=0,fmax=1e-4,     brown,  8,"psyn1",0.5);

   plot_synrate(findsyn(ct,1,dendn_node(dbp1,recpnt1),gca),   fmin=0,fmax=400,     blue,8,"rp1",0.5);
   plot_synrate(findsyn(ct,1,dendn_node(dbp1,recpnt2),gca),   fmin=0,fmax=400,     ltgreen,8,"rp2",0.5);
   plot_synrate(findsyn(ct,1,dendn_node(dbp1,recpnt3),gca),   fmin=0,fmax=400,     ltred,8, "rp3",0.5);
   plot_synrate(findsyn(ct,1,dendn_node(dbp1,recpnt4),gca),   fmin=0,fmax=400,     ltcyan,8,"rp4",0.5);

   plot_syncond(findsyn(dbp1,1,dendn_node(dbp1,presyn1),am),   cmin=0,cmax=200e-12,  red,   10,"psyn1",0.8);
   plot_syncond(findsyn(dbp1,1,dendn_node(dbp1,presyn2),am2),  cmin=0,cmax=200e-12,  green, 10,"psyn2",0.8);
   plot_syncond(findsyn(dbp2,1,dendn_node(dbp2,presyn3),am),   cmin=0,cmax=200e-12,  magenta, 10,"psyn3",0.8);
   plot_syncond(findsyn(dbp2,2,dendn_node(dbp2,presyn3),am2),  cmin=0,cmax=200e-12,  ltgreen, 10,"psyn4",0.8);
    
   Vmin = -0.07; Vmax = -0.01; 
   plot_v_nod(am,  1, findsynlocn(am,1,0,0),   Vmin, Vmax, colr=red,     "", 12, 0.5);
   plot_v_nod(am,  1, findnodloc (am,1,100,0), Vmin, Vmax, colr=blue,    "", 12, 0.5);
   plot_v_nod(am2, 1, findsynlocn(am2,1,0,0),  Vmin, Vmax, colr=green,   "", 12, 0.5);
   // plot_v_nod(am2, 1, 18,  Vmin, Vmax, colr=brown,   "", 12, 0.5);
   // plot_v_nod(am2, 1, 20,  Vmin, Vmax, colr=yellow,   "", 12, 0.5);
   // plot_v_nod(am, 1, findnodloc(am,1,100,0), Vmin, Vmax, colr=white, "", 12, 0.5);
  
   plot_syncond(findsynloc(am,1,0,0),    cmin=0,cmax=500e-12, red,     14,"",0.5);
   plot_syncond(findsynloc(am2,1,0,0),   cmin=0,cmax=500e-12, green,   14,"",0.5);

   /* stimuli */

   if (notinit(stimtype))   stimtype = 1;      


   if (stimtype == 1) { 		// flashed spot and chirp
       s = 0;
       if (nspots==3) {
         stim_spot (spotdia, x=0, y=0, -cmult*scontrast, stimtime, spotduri);
       }
       stim_spot (spotdia, x=0, y=0,  cmult*scontrast, stimtime+spotduri+s*spotdur, spotdur); s++;
       stim_spot (spotdia, x=0, y=0, -cmult*scontrast, stimtime+spotduri+s*spotdur, spotdur); s++;
  
       time2 = spot_chirp (spotdia, x=0, y=0, fstart, fincr, cmult, scontrast, 
		       				stimtime+spotduri+s*spotdur+prestimdur, stimdur);
       if (vcontrast>0) time2 = spot_vcontrast (spotdia, x=0, y=0, temp_freq, 
        					cmult, c1, c2, time2+prestimdur, stimdur);
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
       stimdur = 0.1;
       npixels = 16;
       stim_checkerboard(200, 200, npixels, npixels, orient=0, xoffset=0, yoffset=0,
         	          temp_freq=100, 0, 0.03, stimtime, stimdur, &rndarr, &nframes, rseed);
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


   if (disp) {                        // display the stimulus
         double t, dscale, starttime, disp_end;

      stim_backgr(minten);		// set background for display at simtime=0
      display_size(500);
      disp_end = time2+0.05;
      for (starttime=simtime,t=stimtime; t<disp_end; starttime = t, t+= 0.001) {
           display_stim(starttime, t, dscale=4, -0.037, -0.043);
           //display_stim(t, dscale=4, -0.035, -0.045); 
           //display_stim(0+t, dscale=4); 
           simwait(0.10);
      }
      return;
   }

  if (notinit(setxmin)) setxmin = 0;			// set plot to start at 0
  if (notinit(predur)) predur = 0.5;
  simtime = 0 - predur;
  if (vcontrast>0) {endexp = stimtime+2*stimdur+nspots*spotdur+2*prestimdur+poststimdur;}
  else             {endexp = stimtime+1*stimdur+nspots*spotdur+1*prestimdur+poststimdur;};

  stim_backgr(minten);		// set background for run expt at simtime = -predur

  if (make_cones) synaptau = 0.001;		// set synapses to run with shorter time constant

  //fprintf (stderr,"%d %g %g %g\n", vcloc, vcaxon, stimtime, endexp);
 
  if (!notinit(vcloc)) { 
  	if (elec_rs > 5e6) {
	    make_electrode (nd(ct,cn,electrode_node), ndn(ct,cn,dendn_node(ct,vcloc)), elec_rs, elec_cap);
	    vclamp (ndn(ct,cn,electrode_node), vcaxon, simtime, predur+endexp);
	} else {
	    vclamp (ndn(ct,cn,dendn_node(ct,vcloc)), vcaxon, simtime, predur+endexp);
	}
  }

  // cclamp(ndn(ct,cn,soma), ipre, simtime, predur+stimtime);
  step (predur);

  if (make_cones) synaptau = 1.0;		// set synapses to run with normal time constant

  if (dci>0) { cclamp (ndn(ct,cn,elnode),dci,stimtime,2*stimdur); }

  step (endexp);
}


