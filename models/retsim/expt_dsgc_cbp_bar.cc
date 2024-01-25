/* Experiment dsgc_cbp_bar for retsim */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int n_dsgc;


double theta;
double iroff;
int light_inhib;
int dsgc_prefdir;

int rec_ct;
int rec_cn;

double barwidth;
double minten;
double sinten;
double velocity;
double starttime;
double disptime;
double stimdur;
double endwait;
double sblur;
double ioffset;
double istim;
double predur;
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

  setptr("barwidth",   &barwidth);
  setptr("minten",     &minten);
  setptr("sinten",     &sinten);
  setptr("velocity",   &velocity);
  setptr("starttime",  &starttime);
  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("stimdur",    &stimdur);
  setptr("endwait",    &endwait);
  setptr("sblur",      &sblur);
  setptr("movein",     &movein);
  setptr("predur",     &predur);
  nvalfile = "nval_dsgc_cbp_bar.n";
  fprintf (stderr,"makestim %d\n",makestim);
}

/*--------------------------------------------------------*/

void setparams(void)

  /*  set up default configuration for sb expts */
  /* cones, cone bipolars, sb, dsgc */

{
   int i;
   double zmax, zmin;

  make_rods  = 0;
  make_cones = 0;
  make_dbp1  = 1;
  make_dbp2  = 0;
  make_hbp1  = 0;
  make_hbp2  = 0;
  make_ams   = 1;
  make_sbac  = 0;
  make_dsgc  = 1;

  onplot_dsgc_movie_init();		/* initialize dsgc movie stuff */
  onplot_movie_init();			/* initialize onplot_movie stuff */

  gcdistnod  = 582;
  pickden[dsgc] = 0; //941;       	/* pick one dendrite to allow connections */
  if (notinit(ath_dia)) ath_dia = 0.4;	/* diameter of initial thin axon segment */

  if (n_dsgc<0) n_dsgc = 1;
  if(notinit(rec_ct)) rec_ct = dsgc;    /* type of cell to record from */
  if(notinit(rec_cn)) rec_cn=1;         /* cellnum to record from */

  if(notinit(dsgc_prefdir)) dsgc_prefdir=0;

  setn (dbp1,SDURH2,50);                    /* set 50 ms transient excit input to dsgc */
  setn (ams,SDURH1,50);                     /* set 50 ms transient inhib input to dsgc */

  display_z(zmax=-23, zmin=-20);	    /* exclude dsgc Off-arborization layer */
}

/*--------------------------------------------------------*/

void runexpt(void)

{
    int c, ct, cn, cbp_cn, pl;
    double start, dur, dscale;
    double ixoff, iyoff, disp_end, psize;
    node *npnt;
    elem *e;
    photorec *p;
    chattrib *a;

  cbp_cn    = 68;
  gcdistnod = 582;

  // e = at (ndn(dsgc,1,297), CHAN);
  // a = make_chan (e,NA,2);
  // chset(e);
  // xxx = e->elnum;
  // /* at [dsgc][1][297] chan Na type 2 chset ename xxx; */

  if (notinit(sblur)) sblur = 10;
  if (notinit(stimdur)) stimdur=.45;	/* used for non-moving stimuli */
  if (notinit(endwait)) endwait = 0;

  Vmax  = 0.02;
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

   if (notinit(light_inhib)) light_inhib = 1; 
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
   if (notinit(minten))         minten = -0.045;
   if (notinit(sinten))         sinten =  0.006;
   if (notinit(velocity))     velocity =  2000; 
   if (notinit(starttime))   starttime =  0;
   if (notinit(disptime))     disptime =  0.15;

   stim_backgr(minten,start=0.02);				  	 /* background */

   if (notinit(ioffset)) ioffset = barwidth;
   if (notinit(movein))   movein = 1;
   if (movein) {
     stimdur = movebar (starttime,0,0,300,-300,barwidth,theta,velocity,sinten);	 /* excitatory */
     if (light_inhib)
       stimdur = movebar (starttime,ixoff,iyoff,300+ioffset,-300+ioffset,
		     		barwidth,theta,velocity,sinten); /* inhib */ 
   }
   else {
     stimdur = movebar (starttime,0,0,-300,300,barwidth,theta,velocity,sinten);	 /* excitatory */
     if (light_inhib)
       stimdur = movebar (starttime,ixoff,iyoff,-300+ioffset,300+ioffset,
		     		barwidth,theta,velocity,sinten); /* inhib */ 
   };

   if (disp) {if (light_inhib) display_size(2500); else display_size(500); 
	double t;

      disp_end = 900/velocity;
      for (t=0; t<disp_end; t+= 0.005) {
	   display_stim(starttime+disptime+t, dscale=4); 
	   simwait(0.25);
      }

      return;
    }

   if (notinit(istim)) istim = 0; 
   if (istim != 0) cclamp(ndn(dsgc,1,soma), istim, start=0.02, dur=1);

				/* for ds1d on-layer */
				/* bp  68, n 582   inside spot, at tip */
				/* bp 101, n 529   inside spot */
				/* bp 265, n 469   outside spot, near */
				/* bp 41,  n 422   on dend to soma */
				/* bp 40,  n 297   on dend to soma */
				/* bp 369, n 99    on dend to soma */
				/* bp 241, n 1464  middle top of cell */
				/* bp 488, n 2328  left middle top of cell */
				/* bp 269, n 1336  top of cell */
   if (make_movie) {
     if (space_time) {  /* movie */ 
      plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmaxg,c=1,"Vsoma",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmaxg,c=green,"Vtip1",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmaxg,c=red,"Vtip2",pl=10,0.35);
      //plot_na_inact(dsgc, 1, soma, red, "Na[soma]", 12, .35);
     };
   }
   else { 
    plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmax,c=1,"Vsoma",pl=10,.35);/* V at soma */
    //plot_ca_nod(dsgc,1,soma,cyan,1e-6,"Ca_soma", -1, -1);
    //plot_ca_nod(dsgc,1,gcdistnod,1e-6,yellow,"", -1, -1);
    //plot_v_nod(ct=cbp,cbp_cn,axtrm,   Vmin,Vmax,red,"", -1, 0.35);
    //plot_v_nod(ct=cbp,41,axtrm,      Vmin,Vmax,blue,"", -1, 0.35);
      plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmax,c=green,"Vtip1",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmax,c=red,"Vtip2",pl=10,0.35);
    //plot_v_nod(ct=dsgc,cn=1,422,      Vmin,Vmax,red,"", 10, 0.35);
    //plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,c=red,"Vtip2", pl=10, .35);
    //plot_v_nod(ct=dsgc,cn=1,1336,      Vmin,Vmaxg,c=green,"Vtip1", pl=10, .35);
    //plot_v_nod(ct=dsgc,cn=1,1464,     Vmin,Vmaxg,c=blue,"", -1, -1);
    //plot_v_nod(ct=dsgc,cn=1,2328,     Vmin,Vmaxg,c=magenta,"",pl=10,1);
    //plot_synrate_out(cbp,cbp_cn,0,500,green);
    //plot_synrate_out(cbp,241,0,500,blue);
    //plot_currents(ct=dsgc,plgain=200e-12);

      plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.3); 
    
   };

  if (notinit(predur)) predur=0.00;

  /* set movie plot routine */

  setonplot(onplot_movie);

  /* run experiment */

  setxmin=0;
  simtime=0-predur;
  endexp=starttime+stimdur+endwait;
  step(predur);
  step(starttime+stimdur+endwait);
}
