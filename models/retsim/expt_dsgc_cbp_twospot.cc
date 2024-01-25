/* Experiment dsgc_cbp_twospot for retsim */

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

double diameter;
double minten;
double sinten;
double starttime;
double disptime;
double stimdur;
double endwait;
double sblur;
double ioffset;
double istim;
double predur;
double dur;
double trialtime;
double r1,r2;
double maxtime;
double trialdur;
double setploti;
int movein;
int numtrials;
int spectimes;
int randtimes;
int trial;
int v_ind;
int curve_fit;
int doplots;
int nobiophys;
double *voltages;

#define MAXTRIALS 1000

double trialtimes[MAXTRIALS] = {0.05, 0.1};

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

  setptr("minten",     &minten);
  setptr("sinten",     &sinten);
  setptr("starttime",  &starttime);
  setptr("disptime",   &disptime);
  setptr("ioffset",    &ioffset);
  setptr("istim",      &istim);
  setptr("stimdur",    &stimdur);
  setptr("endwait",    &endwait);
  setptr("sblur",      &sblur);
  setptr("movein",     &movein);
  setptr("predur",     &predur);
  setptr("dur",        &dur);
  setptr("trialtime",   &trialtime);
  setptr("r1",         &r1);
  setptr("r2",         &r2);
  setptr("diameter",   &diameter);
  setptr("numtrials",  &numtrials);
  setptr("maxtime",    &maxtime);
  setptr("trialdur",   &trialdur);
  setptr("spectimes",  &spectimes);
  setptr("randtimes",  &randtimes);
  setptr("setploti",   &setploti);
  setptr("curve_fit",  &curve_fit);
  setptr("doplots",    &doplots);
  setptr("nobiophys",  &nobiophys);
  nvalfile = "nval_dsgc_cbp_bar.n";
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

//  display_z(zmax=-23, zmin=-20);	    /* exclude dsgc Off-arborization layer */
   if (notinit(nobiophys))	nobiophys = 0;

   if(nobiophys)
    setn(dsgc,BIOPHYS,0); //turn off biophysical properties

}

/*--------------------------------------------------------*/

double plottrialtime(double x, double y)

/* function to plot trialtime only at beginning of trial */

{
  static double retval;
  static int trialnum = 0;
  static int stim2 = 0; //1 if second spot hasn't come yet, 0 if it has
  if ((simtime+1e-10) >= trialnum*trialdur) {
     retval = trialtimes[trialnum++];
     stim2 = 1;
  }
  else if((simtime+1e-10) >= (trialnum - stim2)*trialdur + trialtimes[trialnum-stim2])
  {
    retval = trialtimes[trialnum-1];
    stim2 = 0;
  }
  else  
    retval = -1;
//retval = int((simtime+1e-10)/1e-6)*1e-6; //trialtime;  
  return retval;
}

/*--------------------------------------------------------*/
int maxvoltages; 
void storedata(void)
{
  int ct,cn;
  voltages[v_ind] = v(nde(ct=dsgc,cn=1,soma));
  if(v_ind < maxvoltages)
    v_ind++;
}

void initonplot(void)
{
  v_ind = 0;
  maxvoltages = (trialdur*numtrials)/ploti;
  voltages = (double *)emalloc(maxvoltages*sizeof(double));
  setonplot(storedata);
}

int findpeak(void)
{
  double largest = -1e6;
  int lgind = 0;
  for(int i=1; i < maxvoltages; i++)
  {
      if(voltages[i] > largest)
      {
	largest = voltages[i];
	lgind = i;
      }
  }

  return lgind;
}
/*--------------------------------------------------------*/

void runexpt(void)

{
    int c, ct, cn, cbp_cn, pl;
    double start, dscale;
    double ixoff, iyoff, disp_end, psize;
    double fmax, fmin;
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
  if (notinit(endwait)) endwait = 0.1;

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
   if(notinit(numtrials))	numtrials = 2;
   if (notinit(diameter))	diameter = 100;
   if (notinit(r1))		r1 = 0;
   if (notinit(r2))		r2 = 100;
   if (notinit(dur))		dur = 0.01;
   if (notinit(minten))		minten = -0.045;
   if (notinit(sinten))		sinten =  0.006;

   if (notinit(starttime))	starttime =  0;
   if (notinit(disptime))	disptime =  0;
   if (notinit(maxtime))	maxtime = 0.03;
   if (notinit(trialtime))	trialtime = 0.1;
   if (notinit(trialdur))	trialdur = maxtime+starttime+0.1;
   if (notinit(spectimes)) 	spectimes = 0;
   if (notinit(randtimes)) 	randtimes = 0;
   if (!notinit(setploti))	ploti = setploti;
   if (notinit(curve_fit))	curve_fit = 0;
   if (notinit(doplots))	doplots = 1;

   if(randtimes) {
      for(trial = 0; trial < numtrials; trial++) {
	  //make trialtime random in range [-maxtime, maxtime]
	  trialtime = rrange(-maxtime, maxtime);
	  if((trialtime <= dur && trialtime >=0) || (trialtime >= -dur && trialtime <=0)) 
	  {
	    trial--;
	    continue;
	  }
	  trialtimes[trial] = trialtime;
      }
   }

   stim_backgr(minten,start=0.0);				/* background */

   if (notinit(ioffset)) ioffset = diameter;

   for(trial = 0;  trial < numtrials; trial++)			/* set up stimuli */
   {
        double trialt;

      if (spectimes||randtimes) trialt = trialtimes[trial]; 
      else			trialt = trialtime; 

      stimdur = twospot (starttime+trial*trialdur,0,0,r1,r2,diameter,theta,sinten,dur,
		      			trialt);	 /* excitatory */
      if (light_inhib)
      stimdur = twospot (starttime+trial*trialdur,ixoff,iyoff,r1+ioffset,r2+ioffset,diameter,
				      theta,sinten,dur,trialt); /* inhib */
   }
 
   if (disp) {		/* display the stimulus (retsim ... -d 16) */
 	double t;

	if (light_inhib) display_size(2500); else display_size(500); 

        disp_end = trialdur * numtrials;
        for (t=0; t<disp_end; t+= 0.005) {
	   display_stim(t, dscale=4,fmax=-0.03,fmin=-0.05); 
	   if (ninfo>=2) fprintf (stderr,"stim time %-6.4g %6.4g\n",t,disp_end);
	   simwait(0.4);
        }
      return;
    }

   if (notinit(istim)) istim = 0; 
   if (istim != 0) cclamp(ndn(dsgc,1,soma), istim, start=0.2, dur=1);

   if (curve_fit) {
    initonplot();
   }
   if(doplots)
   {
    plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmax,c=1,"Vsoma",pl=10,.35);/* V at soma */
    //plot_ca_nod(dsgc,1,soma,cyan,1e-6,"Ca_soma", -1, -1);
    //plot_ca_nod(dsgc,1,gcdistnod,1e-6,yellow,"", -1, -1);
    //plot_v_nod(ct=cbp,cbp_cn,axtrm,   Vmin,Vmax,red,"", -1, 0.35);
    //plot_v_nod(ct=cbp,41,axtrm,      Vmin,Vmax,blue,"", -1, 0.35);
      plot_v_nod(ct=dsgc,cn=1,1336,Vmin,Vmax,c=green,"Vtip1",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn=1,582,Vmin,Vmax,c=red,"Vtip2",pl=10,0.35);
    //plot_synrate_out(cbp,cbp_cn,0,500,green);
    //plot_synrate_out(cbp,241,0,500,blue);
    //plot_currents(ct=dsgc,plgain=200e-12);

      plot_spike_rate(ct=dsgc, cn=1, soma, red, "spike rate", pl=6, psize=0.2); 

      plot_func(plottrialtime,0.2,0.5,-1); 
      plot_param("trialt", c=magenta, pl=1, 0.1);

   };

  if (notinit(predur)) predur=0.00;

  /* run experiment */

  setxmin=0;
  simtime=0-predur;
  endexp=numtrials*trialdur;
  //step(predur);
  step(endexp);
  
  if(curve_fit)
    fprintf(stdout,"%g",voltages[findpeak()]);

}


