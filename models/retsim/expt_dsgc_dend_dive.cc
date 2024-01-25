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
int on_spot;
int off_spot;
int plot_cond;

int sreg1;
int sreg2;
int sreg3;
int sreg4;
int sreg5;

double dendrm;
double dendcm;
double ddia;
double dvrev;
double dvst;
double dsgc_denddia;
double dcplam;
double g_dbp1_dsgc;
double g_dbp1_nmda;
double g_hbp1_dsgc;
double g_hbp1_nmda;
double n_dbp1_dsgc;
double n_hbp1_dsgc;

double db1_ca6_offm;
double db2_ca6_offm;
double db1_ca7_offm;
double db2_ca7_offm;

double sb_ca6_offm;
double sb_ca6_tauc;
double sb_ca6_taud;

double rstim_theta;
double barwidth;
double minten;
double scontrast;
double spotdia;
double velocity;
double stimtime;
double disptime;
double endwait;
double ioffset;
double istim;
double predur;
double sinc;
double sdur;
double sdur2;
double dispsize;

double spot_offx;
double spot_offy;
double spot_onx;
double spot_ony;
double spot_delay;

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
  
  setptr("sreg1",   	&sreg1);
  setptr("sreg2",   	&sreg2);
  setptr("sreg3",   	&sreg3);
  setptr("sreg4",   	&sreg4);
  setptr("sreg5",   	&sreg5);

  setptr("ddia",	&ddia);
  setptr("dsgc_denddia",&dsgc_denddia);
  setptr("stim_loc",	&stim_loc);
  setptr("nodive",	&nodive);
  setptr("db2_morph",	&db2_morph);
  setptr("db2_biophys",	&db2_biophys);
  setptr("dcplam",	&dcplam);
  setptr("set_cclamp",	&set_cclamp);
  setptr("spotdia",	&spotdia);

  setptr("on_spot",   	&on_spot);
  setptr("off_spot",   	&off_spot);
  setptr("plot_cond", 	&plot_cond);

  setptr("barwidth",   &barwidth);
  setptr("minten",     &minten);
  setptr("scontrast",  &scontrast);
  setptr("velocity",   &velocity);
  setptr("stimtime",   &stimtime);
  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("endwait",    &endwait);
  setptr("movein",     &movein);
  setptr("predur",     &predur);
  nvalfile = "nval_dsgc_dend_dive.n";
  setptr("g_dbp1_dsgc",&g_dbp1_dsgc);
  setptr("g_dbp1_nmda",&g_dbp1_nmda);
  setptr("g_hbp1_dsgc",&g_hbp1_dsgc);
  setptr("g_hbp1_nmda",&g_hbp1_nmda);
  setptr("n_dbp1_dsgc",&n_dbp1_dsgc);
  setptr("n_hbp1_dsgc",&n_hbp1_dsgc);
  setptr("rstim_theta",&rstim_theta);

  setptr("db1_ca6_offm",&db1_ca6_offm);
  setptr("db2_ca6_offm",&db2_ca6_offm);
  setptr("db1_ca7_offm",&db1_ca7_offm);
  setptr("db2_ca7_offm",&db2_ca7_offm);
  setptr("sb_ca6_offm", &sb_ca6_offm);
  setptr("sb_ca6_tauc", &sb_ca6_tauc);
  setptr("sb_ca6_taud", &sb_ca6_taud);

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

  setptr("spot_offx", &spot_offx);
  setptr("spot_offy", &spot_offy);
  setptr("spot_onx",  &spot_onx);
  setptr("spot_ony",  &spot_ony);
  setptr("spot_delay",&spot_delay);

  setptr("sinc",     &sinc);
  setptr("sdur",     &sdur);
  setptr("sdur2",    &sdur2);

  setvar();

  if (notinit(ddia)) ddia = 1.0;
  if (notinit(dsgc_denddia)) dsgc_denddia = 1.0;

  if (notinit (db1_ca6_offm)) db1_ca6_offm = 0.028; // dbp1 Ca6 offsetm 
  if (notinit (db2_ca6_offm)) db2_ca6_offm = 0.028; // dbp1 Ca6 offsetm 
  if (notinit (db1_ca7_offm)) db1_ca7_offm = 0.027; // dbp1 Ca7 offsetm 
  if (notinit (db2_ca7_offm)) db2_ca7_offm = 0.027; // dbp1 Ca7 offsetm 

  if (notinit(sb_ca6_offm))  sb_ca6_offm = 0.026;      // offset for Ca6 chans in sbac dendrites
  if (notinit(sb_ca6_tauc))  sb_ca6_tauc = 8;          // tauc (inactivation) for Ca6 chans in sbac dendrites
  if (notinit(sb_ca6_taud))  sb_ca6_taud = 1;          // taud (reactivation) for Ca6 chans in sbac dendrites

  chanparamsfile = "chanparams_dsgc_dend_dive";
  
    sreg1 = 1001;         // _SREG1, Off-bp input to dsgc,  sets of regions in dens_dsgc_morph_ww.n
    sreg2 = 1002;         // _SREG2, On-bp input to dsgc,  sets of regions in dens_dsgc_morph_ww.n

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
  make_hbp1  = 1;
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
  if (notinit(dcplam)) dcplam = 0.05;
  if (notinit(dvrev))   dvrev = -0.07;
  if (notinit(dvst))     dvst = -0.07;
  if (notinit(db2_morph)) db2_morph = 2;     /* use simple artificial morphology */
  if (notinit(db2_biophys)) db2_biophys = 0; /* don't use density file */
  if (notinit(set_cclamp))   set_cclamp = 0; /* set cclamp instead of synaptic input */
  if (notinit(g_dbp1_dsgc)) g_dbp1_dsgc = 2e-9; /* set On synpatic input conductance */
  if (notinit(g_dbp1_nmda)) g_dbp1_nmda = 0;
  if (notinit(g_hbp1_dsgc)) g_hbp1_dsgc = 2e-9; /* set Off synpatic input conductance */
  if (notinit(g_hbp1_nmda)) g_hbp1_nmda = 0;
  if (notinit(n_dbp1_dsgc)) n_dbp1_dsgc = 0;
  if (notinit(n_hbp1_dsgc)) n_hbp1_dsgc = 0;

  if (notinit(plot_cond))  plot_cond = 0;

  SOMA      = R_6;			/* set soma to be R6 */

  setn(dsgc,NCOLOR,RCOLOR);    		/* set cell display color from region */
  setn(dsgc,SOMAZ,-32);    		/* set cell display color from region */
  // setn(dbp1,NCOLOR,NOCOLOR);   		/* set dbp1 to be transparent */
  setn(dbp1,NCOLOR,ltblue);   		/* set dbp1 to be transparent */
  setn(hbp1,NCOLOR,ltred);   		/* set dbp1 to be transparent */
  
  // setn (dbp1,SDURH2,50);                /* set 50 ms transient excit input to dsgc */
  // setn (ams,SDURH1,50);                 /* set 50 ms transient inhib input to dsgc */
  setsv (dbp1,SCOND,2,g_dbp1_dsgc);        /* set dbp1 conductance */
  setsv (dbp1,SCOND,8,g_dbp1_nmda);        /* set dbp1 nmda conductance */
  setsv (hbp1,SCOND,2,g_hbp1_dsgc);        /* set dbp1 conductance */
  setsv (hbp1,SCOND,8,g_hbp1_nmda);        /* set dbp1 nmda conductance */

  setsv (dbp1,SVNOISE,2,n_dbp1_dsgc);      /* set noise in dbp1 ampa conductance */
  setsv (hbp1,SVNOISE,2,n_hbp1_dsgc);      /* set noise in hbp1 ampa conductance */

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
//  plotlabel = "";

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
		//if (set_cclamp==0) {
		//   n_dbp1 = 1;
		//   make_one_dbp1 = 1;
                //   gcdistnod  = stim_loc;
		//}
		break;
   }

      if (dend_layers) {
         movie_title1 = dend_text1;
         movie_title2 = dend_text2;
        if (plot_cond) {
           if (set_cclamp > 0) 
            sprintf   (tbuf1,"denddia %g  %g pA", dsgc_denddia, istim*1e12);
          else {
          // sprintf   (tbuf1,"denddia %g  off  %g pS", dsgc_denddia, g_hbp1_dsgc*1e12);
          }
          gframe ("/");
          gpen (white);
          gmove (0.75, 0.93);
          sprintf   (tbuf1,"%d off  %g pS", off_spot, g_hbp1_dsgc*1e12);
          gtext (tbuf1);
          gmove (0.75, 0.90);
          sprintf   (tbuf1,"%d nmda %g pS", off_spot, g_hbp1_nmda*1e12);
          gtext (tbuf1);
          gmove (0.75, 0.87);
          sprintf   (tbuf1,"%d on   %g pS", on_spot, g_dbp1_dsgc*1e12);
          gtext (tbuf1);
          gmove (0.75, 0.84);
          sprintf   (tbuf1,"%d nmda %g pS", on_spot, g_dbp1_nmda*1e12);
          gtext (tbuf1);
	}
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

}

/*--------------------------------------------------------*/

void runexpt(void)

{
    int c, ct, cn, cbp_cn, i, pl, t;
    int off_node1, off_node2, off_node3, off_node4;
    int on_node1, on_node2; 
    int stimchan;
    double start, dur, wavel, mask;
    double dscale, cmin, cmax;
    double ixoff, iyoff, starttime, disp_end, psize;
    double Vmax, Vmaxg, Vmin, plsiz;
    node *npnt;
    elem *e;
    photorec *p;
    chattrib *a;
    synapse *s;
    char plbuf[50];
    double spots[4][2] = { 80,  20,
	                   20, -60, 
	                  -80,  20, 
	                  -60, -50 };

  timinc = 1e-6;
  implicit = 1;

  cbp_cn    = 68;

  // e = at (ndn(dsgc,1,297), CHAN);
  // a = make_chan (e,NA,2);
  // chset(e);
  // xxx = e->elnum;
  // /* at [dsgc][1][297] chan Na type 2 chset ename xxx; */

  if (notinit(endwait)) endwait = 0.005;

  Vmax  = 0.00;
  Vmaxg = 0.00;
  Vmin = -0.08;

   /* add light transducer to each bipolar cell */

   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma),stimchan=2); 
     // p = (photorec*)make_transducer(ndn(dbp1,cn,soma)); 
     p->xpos=npnt->xloc; 
     p->ypos=npnt->yloc;
   }
   for(npnt=nodepnt; npnt=foreach(npnt,hbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(hbp1,cn,soma),stimchan=1); 
     // p = (photorec*)make_transducer(ndn(hbp1,cn,soma)); 
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
   if (notinit(minten))         minten = -0.045;
   if (notinit(scontrast))   scontrast =  0.005;
   if (notinit(velocity))     velocity =  2000; 
   if (notinit(stimtime))     stimtime =  0.005;
   if (notinit(disptime))     disptime =  0.15;
   if (notinit(sdur))             sdur = 0.05;
   if (notinit(sdur2))           sdur2 = 0.01;

   if (notinit(spot_delay))   spot_delay =  0.005;

   if (notinit(rstim_theta)) rstim_theta =  0;


   if (notinit(ioffset))    ioffset = barwidth;
   if (notinit(movein))      movein = 0;
   if (notinit(predur))      predur = 0.01;

   simtime = -predur;
   setxmin = 0;

   stim_backgr(minten,stimchan=1,start=simtime);				  	 /* background */
   stim_backgr(minten,stimchan=2,start=simtime);				  	 /* background */
   // stim_spot(spotdia, -spotx, -spoty, scontrast, stimtime, sdur);
  
  /*  
   if (off_spot==1) {
       spot_offx = 30; spot_offy = -50;
   } else
   if (off_spot==2) {
       spot_offx = -50; spot_offy = 20;
   }

   if (on_spot==1) {
       spot_onx = 30; spot_ony = -50;
   } else 
   if (on_spot==2) {
       spot_onx = -50; spot_ony = 20;
   }
  /* */

   if (notinit(on_spot))    on_spot = 1;
   if (notinit(off_spot))  off_spot = 1;
   if (off_spot > 4) off_spot = 4;
   if (on_spot > 4)  on_spot = 4;

   if (notinit(spot_offx))  spot_offx = spots[off_spot-1][0];
   if (notinit(spot_offy))  spot_offy = spots[off_spot-1][1];
   if (notinit(spot_onx))   spot_onx  = spots[on_spot-1][0];
   if (notinit(spot_ony))   spot_ony  = spots[on_spot-1][1];

   stim_spot(spotdia, -spot_offx, -spot_offy, scontrast, stimtime, sdur, wavel=1, mask=0, stimchan=1);
   stim_spot(spotdia, -spot_onx,  -spot_ony,  scontrast, stimtime+spot_delay, sdur2, wavel=1, mask=0, stimchan=2);

   if (plot_cond) {
     plot_syncond(findsynloc(hbp1,-1,dsgc,1,-spot_offx,-spot_offy,0.0),    cmin=0,cmax=10e-10, red,     16,"",0.15);
     plot_syncond(findsynloc(hbp1,-1,dsgc,1,-spot_offx,-spot_offy,-0.002), cmin=0,cmax=10e-10, yellow, 16,"Ghbp1_nmda",0.15);

     plot_syncond(findsynloc(dbp1,-1,dsgc,1,-spot_onx, -spot_ony,0.0),   cmin=0,cmax=10e-10, green,   16,"",0.15);
     plot_syncond(findsynloc(dbp1,-1,dsgc,1,-spot_onx, -spot_ony,-0.002),cmin=0,cmax=10e-10, brown,   16,"Gdbp1_nmda",0.15);
   }

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

      // disp_end = 900/velocity;
      disp_end = 0.05;
      for (starttime=simtime,t=0; t<=disp_end; starttime=t, t+= 0.001) {
           display_stim(starttime, t, dscale=4, -0.04, -0.05);
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
       } 
   }

   if (!noplots) {
     if      (setscene==5) plsiz = 0.2;
     else if (setscene<5)  plsiz = 0.3;
     else                  plsiz = 1;
     if (plot_cond==0) plsiz = 0.25;

     /*
     plot_v_nod(ct=dsgc,cn=1,stim1,Vmin,Vmax,c=blue, "S1",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim2,Vmin,Vmax,c=green,"S2",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim3,Vmin,Vmax,c=cyan, "S3",pl=20,0.35*plsiz);
     plot_v_nod(ct=dsgc,cn=1,stim4,Vmin,Vmax,c=red,  "S4",pl=20,0.35*plsiz);
     */

     // plot_v_nod(ct=dbp1,cn=130,soma,-0.06,-0.04,c=blue,    "", pl=25, 1.0*plsiz);/* V at soma */
     // plot_v_nod(ct=dbp1,cn=130,1,   -0.06,-0.04,c=green,   "", pl=25, 1.0*plsiz);/* V at soma */

     cmin = 0; cmax = 2e-10;
     // plot_synrate(findsynloca(dbp1,dsgc,1,-100,rstim_theta),0,400, red,   17,"",0.5);
     // plot_syncond(findsynloca(dbp1,dsgc,1,-100,rstim_theta),cmin,cmax, red,   16,"",0.5);
     // s = findsynloca(dbp1,dsgc,1,-100,rstim_theta);

     off_node1 = findnodlocr(ct=dsgc,cn=1, -spot_offx, -spot_offy, sreg1);
     sprintf (plbuf,"Voff spot %d",off_spot);
     plot_v_nod(dsgc,1,off_node1, Vmin,Vmax,c=red, plbuf, pl=10,plsiz);

     on_node1 = findnodlocr(ct=dsgc,cn=1, -spot_onx, -spot_ony, sreg2);
     sprintf (plbuf,"Von  spot %d",on_spot);
     plot_v_nod(dsgc,1,on_node1, Vmin,Vmax,c=green, plbuf, pl=10,plsiz);

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
   //plot_synrate_out(cbp,241,0,500,blue); //plot_currents(ct=dsgc,plgain=200e-12);

   // plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.3); 
   }


  /* set movie plot routine */

    setonplot(onplot_movie);

  /* run experiment */

  endexp=stimtime+sdur+sdur2+endwait;
  step(predur);
  step(stimtime+sdur+sdur2+endwait);
}
