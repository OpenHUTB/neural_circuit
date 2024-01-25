/* Experiment dsgc_cbp_bar for retsim */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "gprim.h"
#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int n_dsgc;


double theta;
double iroff;
int light_inhib;
int dsgc_prefdir;

int rec_ct;
int rec_cn;

int rec1;
int rec2;
int rec3;
int rec4;
int rec5;
int rec6;
int rec7;
int rec8;

int stim1;
int stim2;
int stim3;
int stim4;
int stim5;
int setscene;
int stim_loc;
int nodive;
int db2_morph;
int db2_biophys;
int set_cclamp;

double dendrm;
double dendcm;
double ddia;
double dvrev;
double dvst;
double dsgc_denddia;
double dcplam;
double g_dbp1_dsgc;

double barwidth;
double minten;
double sinten;
double spotdia;
double velocity;
double stimtime;
double disptime;
double stimdur;
double endwait;
double sblur;
double ioffset;
double istim;
double predur;
double sinc;
double sdur;
double dispsize;

const char *dend_text1 = NULL;
const char *dend_text2 = NULL;

int movein;

/*--------------------------------------------------------*/

void defparams(void)

{
  defparams_dsgc_movie();
  defparams_onplot_movie();

  setptr("theta", 	 &theta);
  setptr("iroff", 	 &iroff);
  setptr("light_inhib",  &light_inhib);
  setptr("rec_cn", 	 &rec_cn);
  setptr("rec_ct", 	 &rec_ct);
  setptr("dsgc_prefdir", &dsgc_prefdir);

  setptr("dendrm",   	&dendrm);
  setptr("dendcm",   	&dendcm);
  setptr("dvrev",   	&dvrev);
  setptr("dvst",   	&dvst);

  setptr("ddia",	&ddia);
  setptr("dsgc_denddia",&dsgc_denddia);
  setptr("stim_loc",	&stim_loc);
  setptr("nodive",	&nodive);
  setptr("db2_morph",	&db2_morph);
  setptr("db2_biophys",	&db2_biophys);
  setptr("dcplam",	&dcplam);
  setptr("set_cclamp",	&set_cclamp);
  setptr("spotdia",	&spotdia);

  setptr("barwidth",   &barwidth);
  setptr("minten",     &minten);
  setptr("sinten",     &sinten);
  setptr("velocity",   &velocity);
  setptr("stimtime",   &stimtime);
  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("stimdur",    &stimdur);
  setptr("endwait",    &endwait);
  setptr("sblur",      &sblur);
  setptr("movein",     &movein);
  setptr("predur",     &predur);
  nvalfile = "nval_dsgc_dend_dive.n";
  setptr("g_dbp1_dsgc",&g_dbp1_dsgc);

  setptr("rec1",     &rec1);
  setptr("rec2",     &rec2);
  setptr("rec3",     &rec3);
  setptr("rec4",     &rec4);
  setptr("rec5",     &rec5);
  setptr("rec6",     &rec6);

  setptr("stim1",    &stim1);
  setptr("stim2",    &stim2);
  setptr("stim3",    &stim3);
  setptr("stim4",    &stim4);
  setptr("stim5",    &stim5);
  setptr("setscene", &setscene);

  setptr("sinc",     &sinc);
  setptr("sdur",     &sdur);

  setvar();

  if (notinit(ddia)) ddia = 1.0;
  if (notinit(dsgc_denddia)) dsgc_denddia = 1.0;

}

/*--------------------------------------------------------*/

void setparams(void)

  /*  set up default configuration for sb expts */
  /* cones, cone bipolars, sb, dsgc */

{
   int i;
   double zmax, zmin;
   char tbuf1[200] = {0};
   char tbuf2[200] = {0};

  make_rods  = 0;
  make_cones = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_hbp1  = 0;
  make_hbp2  = 0;
  make_ams   = 1;
  make_sbac  = 0;
  make_dsgc  = 1;

  pickden[dsgc] = 0; //941;       	/* pick one dendrite to allow connections */
  if (notinit(ath_dia)) ath_dia = 0.4;	/* diameter of initial thin axon segment */

  if (n_dsgc<0) n_dsgc = 1;
  if(notinit(rec_ct)) rec_ct = dsgc;    /* type of cell to record from */
  if(notinit(rec_cn)) rec_cn=1;         /* cellnum to record from */

  if(notinit(dsgc_prefdir)) dsgc_prefdir=0;

  if (notinit(dendrm)) dendrm = drm;
  if (notinit(dendcm)) dendcm = dcm;
  if (notinit(dcplam)) dcplam = 0.01;
  if (notinit(dvrev))   dvrev = -0.07;
  if (notinit(dvst))     dvst = -0.07;
  if (notinit(db2_morph)) db2_morph = 2;     /* use simple artificial morphology */
  if (notinit(db2_biophys)) db2_biophys = 0; /* don't use density file */
  if (notinit(set_cclamp))   set_cclamp = 0; /* set cclamp instead of synaptic input */
  if (notinit(g_dbp1_dsgc)) g_dbp1_dsgc = 2e-9; /* set synpatic stimulus conductance */

  SOMA      = R_6;			/* set soma to be R6 */

  setn(dsgc,NCOLOR,RCOLOR);    		/* set cell display color from region */
  setn(dsgc,SOMAZ,-32);    		/* set cell display color from region */
  setn(dbp1,NCOLOR,NOCOLOR);   		/* set dbp1 to be transparent */

  // setn (dbp1,SDURH2,50);                /* set 50 ms transient excit input to dsgc */
  // setn (ams,SDURH1,50);                 /* set 50 ms transient inhib input to dsgc */
  setsv (dbp1,SCOND,2,g_dbp1_dsgc);        /* set dbp1 conductance */

  // disp_dsgc_zmax = 0;
  // disp_dsgc_zmin = -28;
  
 if (notinit(istim)) istim = 0; 

// if (!notinit(setscene)) {
//   switch (setscene) {
//	default:
// 	case 1:	_COL = C_COLOR;  dend_text = "dendrite region colors"; break;
// 	case 2:	_COL = C_COLOR2; dend_text = "all dendrites vcolor"; break;
// 	case 3:	_COL = C_COLOR3; dend_text = "Off layer"; break;
// 	case 4:	_COL = C_COLOR4; dend_text = "On layer"; break;
//   }
// }
// sprintf (tbuf1," denddia %g, %g pA: soma, S1 %d S2 %d S3 %d S4 %d",
//		   dsgc_denddia, istim*1e12,stim1,stim2,stim3,stim4);
//
// sprintf (tbuf1,"denddia %g, %g pA", dsgc_denddia, istim*1e12);
 plotlabel = "";

 if (notinit(nodive)) nodive = 0;
 if (!notinit(setscene)) {
   switch (setscene) {
	default:
 	case 1:	_COL = C_COLOR;  space_time = 1; dend_text1 = "dendrite region colors"; break;
 	case 2:	_COL = C_COLOR2; space_time = 1; dend_text1 = "all dendrites vcolor"; break;
 	case 3:	_COL = C_COLOR3; space_time = 1; dend_text1 = "Off layer"; break;
 	case 4:	_COL = C_COLOR4; space_time = 1; dend_text1 = "On layer"; break;
 	case 5: dend_layers = 1; space_time = 0; 
		denddens_color1 = C_COLOR3; 
		denddens_color1t = C_COLOR5; 
		if      (!notinit(stim1)) stim_loc = stim1;
		else if (!notinit(stim2)) stim_loc = stim2;
		if (!notinit(stim_loc)) {sprintf (tbuf2,"Off layer, stim_loc %d",stim_loc); dend_text1 = tbuf2;}
		else dend_text1 = "Off layer";
		denddens_color2 = C_COLOR4; 
		denddens_color2t = C_COLOR6; 
		dend_text2 = "On layer";
		if (set_cclamp==0) {
		   n_dbp1 = 1;
		   make_one_dbp1 = 1;
                   gcdistnod  = stim_loc;
		}
		break;
   }

   if (dend_layers) {
      movie_title1 = dend_text1;
      movie_title2 = dend_text2;
      if (set_cclamp > 0) 
        sprintf   (tbuf1,"denddia %g, %g pA", dsgc_denddia, istim*1e12);
      else
        sprintf   (tbuf1,"denddia %g, %g pS", dsgc_denddia, g_dbp1_dsgc*1e12);
      gframe ("/");
      gpen (white);
      gmove (0.7, 0.975);
      gtext (tbuf1);
   } else { 
      sprintf   (tbuf1,"%s, denddia %g, %g pA", dend_text1,dsgc_denddia, istim*1e12);
      movie_title1 = tbuf1;
   }
 }
 dend_vcolor1 = rcolor;
 dend_vcolor2 = rcolor;

 onplot_dsgc_movie_init();             /* initialize dsgc movie stuff */
 onplot_movie_init();                  /* initialize onplot_movie stuff */

}

/*--------------------------------------------------------*/

/* orig, for dens_dsgc_morph_ww.n */
/* 
int stimlist_green[] = {778,  2330,  2360, 144, 690,		// tips
			5007, 5727, 5743, 5788, 5715,
			5681, 4959, 5563, 3201, 3177, 
			1423, 673,  3306, 5269, 5941,
			5871, 5824, 5170, 5842, 3718,
			3760, 1959, 2417, 920,  2576,
			4271, 65,   
			
			750,  100,  5,    81,   122, 		// intermediate     
			14,   34,   4980, 5770, 4966, 
			5690, 4940, 4905, 3474, 3458, 
			3188, 3152, 1400, 1335, 1304, 
			1270, 3270, 3240, 3215, 650,  
			5960, 5240, 5226, 5915, 3724, 
			1916, 3740, 5155, 1940, 3705, 
			617,  800,  847,  2390, 4290, 
			2555, 2592,
			0};

int stimlist_red[] = {	1566, 1609, 6035, 5522, 1446, 		// diving dendrites, tips
			4475, 2756, 5395, 4588, 4206, 
			5087, 3692, 2086, 3781, 5369, 
			4096, 3800, 5552, 
		
			4500, 4520, 4540, 4565, 5384, 		// intermediate	 
			2670, 4428, 4444, 2734, 5420, 
			5460, 5490, 5510, 1502, 1520, 
			1540, 1555, 1585, 1596, 3655, 
			3675, 5792, 4157, 4177, 4193, 
			1980, 2025, 2065, 3772, 4080, 
			4068, 5354, 5535,
			0};

/* */

int stimlist_green[] = {3817, 6154, 6190, 1748, 3707,		// tips
			7199, 7891, 7922, 7973, 7881,
			7840, 7147, 7693, 5322, 5287, 
			3294, 654,  5499, 6411, 7573,
			7497, 7424, 6285, 7453, 3859,
			3923, 1797, 1337, 287,  1564,
			3663, 477,   
		
			3780, 1691, 1638, 1659, 1725,		// intermediate     
			421,  450,  7172, 7952, 7157, 
			7850, 7117, 7070, 5757, 5735, 
			5304, 5260, 3262, 3196, 3160, 
			3110, 5446, 5386, 5341, 624,  
			7590, 6372, 6350, 7543, 3868, 
			1758, 3892, 6263, 1778, 3839, 
			575,  113,  184,  1306, 3682, 
			1540, 1586,
			0};


int stimlist_red[] = {	2455, 2507, 8003, 7387, 2272, 		// diving dendrites, tips
			5989, 3455, 7236, 6126, 4566, 
			6647, 4864, 3954, 6535, 4419, 
			4009, 7677, 
		
			6019, 6036, 6057, 6091, 7217, 		// intermediate	 
			3344, 5930, 5954, 3421, 7264, 
			7300, 7332, 7368, 2346, 2378, 
			2410, 2434, 2480, 2493, 4810, 
			4844, 7651, 4486, 4521, 4545, 
			1830, 1885, 3940, 4394, 4377, 
			6520, 7661,
			0};

int stimlist_green_len = 0;
int stimlist_red_len = 0;


void addlabels(void)

/*  if (!notinit(stim1)) label(nde(c,1,dendn_node(c,stim1)),yellow," -S1"); */
/*   use dendn_node to refer to node on same line as label */

{
   int c=dsgc;
   int i;

  // if (!notinit(stim1)) label(nde(c,1,dendn_node(c,201)),yellow," -S1");  // label column 201
  
 // if (!notinit(stim_loc)) stim1 = stim_loc;
 
  // if (make_movie && setscene!=4) {
  if (setscene!=4) {
      // if (!notinit(stim1)) label(nde(c,1,stim1),red,    " -Stim_Off_arb");
      // if (!notinit(stim2)) label(nde(c,1,stim2),magenta," -Stim_Cross");
      // if (!notinit(stim1)) label(nde(c,1,stim1),red,    " -SO");
      // if (!notinit(stim2)) label(nde(c,1,stim2),magenta," -SC");
      // if (!notinit(stim3)) label(nde(c,1,stim3),green,  " -Off3");
      // if (!notinit(stim4)) label(nde(c,1,stim4),ltgreen," -Off4");
      if (!notinit(stim1)) label(nde(c,1,stim1),red,    " *");
      if (!notinit(stim2)) label(nde(c,1,stim2),magenta," *");
  }
  if (setscene!=3) {
      // if (!notinit(rec1))  label(nde(c,1,rec1),blue,  " -On1");
      // if (!notinit(rec2))  label(nde(c,1,rec2),ltblue," -On2");
      // if (!notinit(rec3))  label(nde(c,1,rec3),cyan,  " -On3");
      // if (!notinit(rec4))  label(nde(c,1,rec4),ltcyan," -On4");
   }
   if (notinit(setscene)) node_color = RCOLOR;

/* 
   if (notinit(stim_loc)) {
      if (!nodive) {
          for (i=0; i<100; i++) { if (stimlist_red[i]==0) break;}
          stimlist_red_len = i;
     }
     for (i=0; i<100; i++) { if (stimlist_green[i]==0) break;}
       stimlist_green_len = i;
     fprintf (stderr,"# stimlist red %d  green %d  total %d\n",stimlist_red_len, stimlist_green_len, 
			  stimlist_red_len + stimlist_green_len);
   }
/* */

}

/*--------------------------------------------------------*/

void runexpt(void)

{
    int c, ct, cn, cbp_cn, i, pl, t;
    int off_node1, off_node2, off_node3, off_node4;
    int on_node1, on_node2; 
    double start, dur, dscale;
    double ixoff, iyoff, disp_end, psize;
    double Vmax, Vmaxg, Vmin, plsiz;
    node *npnt;
    elem *e;
    photorec *p;
    chattrib *a;

  timinc = 1e-6;
  implicit = 1;

  cbp_cn    = 68;

  // e = at (ndn(dsgc,1,297), CHAN);
  // a = make_chan (e,NA,2);
  // chset(e);
  // xxx = e->elnum;
  // /* at [dsgc][1][297] chan Na type 2 chset ename xxx; */

  if (notinit(sblur)) sblur = 10;
  if (notinit(stimdur)) stimdur= 0.45;	/* used for non-moving stimuli */
  if (notinit(endwait)) endwait = 0.02;

  Vmax  = 0.00;
  Vmaxg = 0.00;
  Vmin = -0.08;

   /* add light transducer to each bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma)); 
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

     for(npnt=nodepnt; npnt=foreach(npnt,ams,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
       p = (photorec*)make_transducer(ndn(ams,cn,soma)); 
       p->xpos=npnt->xloc + ixoff; 
       p->ypos=npnt->yloc + iyoff;
     }
   }

  /*  - - - - - - - - - - - - - - - - - - - */

   if (notinit(barwidth))     barwidth = 100;
   if (notinit(spotdia))       spotdia = 500;
   if (notinit(minten))         minten = -0.05;
   if (notinit(sinten))         sinten =  0.02;
   if (notinit(velocity))     velocity =  2000; 
   if (notinit(stimtime))     stimtime =  0.01;
   if (notinit(disptime))     disptime =  0.15;

   /* stim_backgr(minten,start=0.02);				  	 /* background */

   if (notinit(ioffset)) ioffset = barwidth;
   if (notinit(movein))   movein = 1;

//   if (movein) {
//     stimdur = movebar (stimtime,0,0,300,-300,barwidth,theta,velocity,sinten);	 /* excitatory */
//     if (light_inhib)
//       stimdur = movebar (stimtime,ixoff,iyoff,300+ioffset,-300+ioffset,
//		     		barwidth,theta,velocity,sinten); /* inhib */ 
//   }
//   else {
//     stimdur = movebar (stimtime,0,0,-300,300,barwidth,theta,velocity,sinten);	 /* excitatory */
//     if (light_inhib)
//       stimdur = movebar (stimtime,ixoff,iyoff,-300+ioffset,300+ioffset,
//		     		barwidth,theta,velocity,sinten); /* inhib */ 
//   };


   if (disp & DSTIM) {
	double t;

      if (light_inhib) dispsize = 2500;
      else 	       dispsize = 300;
      set_disp_rot (mxrot,myrot,mzrot,arrcentx,arrcenty,arrcentz,0,0,0,dispsize);

      disp_end = 900/velocity;
      for (t=0; t<disp_end; t+= 0.005) {
           display_stim(stimtime, t, dscale=4, -0.02, -0.05);
	   simwait(0.25);
      }

      return;
    }

   // node in Off layer distal tip that dives to On layer: 2080
   // node in Off layer distal tip that dives to On layer: 4200
   // node in Off layer distal tip that dives to On layer: 4475
   // node in Off layer distal tip that dives to On layer: 6006

   // node in On layer ~10 um proximal from diving input: 2250
   // node in On layer distal tip from diving input: 4015

   if (notinit(sdur)) sdur = 0.01;
   if (notinit(sinc)) sinc = sdur + 0.02;
   t=0;
   // if (istim != 0) cclamp(ndn(dsgc,1,soma),  istim, start=stimtime+t++*sinc, sdur);
   // if (istim != 0) cclamp(ndn(dsgc,1,stim1), istim, start=stimtime+t++*sinc, sdur);
   // if (istim != 0) cclamp(ndn(dsgc,1,stim2), istim, start=stimtime+t++*sinc, sdur);
   // if (istim != 0) cclamp(ndn(dsgc,1,stim3), istim, start=stimtime+t++*sinc, sdur);
   // if (istim != 0) cclamp(ndn(dsgc,1,stim4), istim, start=stimtime+t++*sinc, sdur);

   /* stimulate the Off layer in multiple places */

   if (istim != 0) {
	  double totcur;

       t = 0;
       if (!notinit(stim_loc)) {
          if (set_cclamp>0) 
		 cclamp(ndn(dsgc,1,stim_loc), istim, start=stimtime+t*sinc, sdur);
	  else {
		 stim_backgr(minten);
		 stim_spot(spotdia, 0, 0, sinten, stimtime, sdur);
	  }

       } else {
         for (totcur=i=0; i<stimlist_red_len; i++) {		// the diving dendrites 
	  // if (ddia > 0.5) { 
          //   cclamp(ndn(dsgc,1,stimlist_red[i]), istim, start=stimtime+t*sinc, sdur);
	  //   totcur += istim;
	  // }
	  if (!nodive) {
               cclamp(ndn(dsgc,1,stimlist_red[i]), istim, start=stimtime+t*sinc, sdur);
	       totcur += istim;
	  }
         }
         for (i=0; i<stimlist_green_len; i++) {		// the regular Off denddrites
           cclamp(ndn(dsgc,1,stimlist_green[i]), istim, start=stimtime+t*sinc, sdur);
	   totcur += istim;
         }
         fprintf (stderr,"# total current injected %g\n",totcur);
      }
   }

   t++;
   stimdur = stimtime + t*(sdur+sinc);

   if (!noplots) {
     if      (setscene==5) plsiz = 0.25;
     else if (setscene<5)  plsiz = 0.3;
     else                  plsiz = 1;
     /*
     plot_v_nod(ct=dsgc,cn=1,stim1,Vmin,Vmax,c=blue, "S1",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim2,Vmin,Vmax,c=green,"S2",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim3,Vmin,Vmax,c=cyan, "S3",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim4,Vmin,Vmax,c=red,  "S4",pl=20,0.35*plsiz);
     */

     // plot_v_nod(ct=dbp1,cn=1,soma,Vmin,Vmax,c=blue,    "Vdbp1", pl=25, 1.0*plsiz);/* V at soma */

     plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmax,c=blue,    "Vsoma", pl=10, 1.0*plsiz);/* V at soma */
     if (!notinit(stim1)) plot_v_nod(ct=dsgc,cn=1,stim1,Vmin,Vmax,c=green,  "Stim_Off_arb", pl=10, 1.0*plsiz);
     if (!notinit(stim2)) plot_v_nod(ct=dsgc,cn=1,stim2,Vmin,Vmax,c=green,  "Stim_Cross",   pl=10, 1.0*plsiz);
     // plot_v_nod(ct=dsgc,cn=1,stim1,Vmin,Vmax,c=red,    "Off_dive1", pl=10, 1.0*plsiz);
     // plot_v_nod(ct=dsgc,cn=1,stim2,Vmin,Vmax,c=magenta,"Off_dive2", pl=10, 1.0*plsiz);
     // plot_v_nod(ct=dsgc,cn=1,stim3,Vmin,Vmax,c=green,  "Off3", pl=10, 1.0*plsiz);
     // plot_v_nod(ct=dsgc,cn=1,stim4,Vmin,Vmax,c=ltgreen,"Off4", pl=10, 1.0*plsiz);
     if (!notinit(rec1)) plot_v_nod(ct=dsgc,cn=1,rec1,Vmin,Vmax,c=red,    "On1", pl=10, 1.0*plsiz);
     if (!notinit(rec2)) plot_v_nod(ct=dsgc,cn=1,rec2,Vmin,Vmax,c=brown,  "On2", pl=10, 1.0*plsiz);
     if (!notinit(rec3)) plot_v_nod(ct=dsgc,cn=1,rec3,Vmin,Vmax,c=cyan,    "On3", pl=10, 1.0*plsiz);
     if (!notinit(rec4)) plot_v_nod(ct=dsgc,cn=1,rec4,Vmin,Vmax,c=ltcyan,  "On4", pl=10, 1.0*plsiz);

   //plot_ca_nod(dsgc,1,soma,cyan,1e-6,"Ca_soma", -1, -1);
   //plot_ca_nod(dsgc,1,gcdistnod,1e-6,yellow,"", -1, -1);
   //plot_v_nod(ct=cbp,cbp_cn,axtrm,   Vmin,Vmax,red,"", -1, 0.35);
   // plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmax,c=green,"Vtip1",pl=10,0.35);
   // plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmax,c=red,"Vtip2",pl=10,0.35);
   //plot_v_nod(ct=dsgc,cn=1,422,      Vmin,Vmax,red,"", 10, 0.35);
   //plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,c=red,"Vtip2", pl=10, .35);
   //plot_v_nod(ct=dsgc,cn=1,1336,      Vmin,Vmaxg,c=green,"Vtip1", pl=10, .35);
   //plot_v_nod(ct=dsgc,cn=1,1464,     Vmin,Vmaxg,c=blue,"", -1, -1);
   //plot_v_nod(ct=dsgc,cn=1,2328,     Vmin,Vmaxg,c=magenta,"",pl=10,1);
   //plot_synrate_out(cbp,cbp_cn,0,500,green);
   //plot_synrate_out(cbp,241,0,500,blue);
   //plot_currents(ct=dsgc,plgain=200e-12);

   // plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.3); 
   }

  if (notinit(predur)) predur=0.00;

  /* set movie plot routine */

    setonplot(onplot_movie);

  /* run experiment */

  setxmin=0;
  simtime=0-predur;
  endexp=stimtime+sdur+endwait;
  step(predur);
  step(stimtime+sdur+endwait);
}
