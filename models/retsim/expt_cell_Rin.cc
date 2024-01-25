/* Experiment gc_Rin */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "stimfuncs.h"

int light_inhib;
int dsgc_prefdir;
int passive;

double velocity;
double starttime;
double disptime;
double stimdur;
double endwait;
double sblur;
double ioffset;
double istim;
double predur;
double dt;
double dv;
double dst;
double vstart;
double vstop;
double vstep;
double fracnodes;
double set_drm;

double axon_br_dia;
double varicos_dia;
const char *celltype;
int ct;
extern int cumcomp;

char savefile[30] = {0};

/*-------------------------------------------------*/

void defparams(void)

{
  setptr("dsgc_prefdir", &dsgc_prefdir);
  setptr("starttime",	&starttime);
  setptr("disptime",	&disptime);
  setptr("ioffset",	&ioffset);
  setptr("istim",	&istim);
  setptr("stimdur",	&stimdur);
  setptr("endwait",	&endwait);
  setptr("predur",	&predur);
  setptr("dt",		&dt);
  setptr("dv",		&dv);
  setptr("dst",		&dst);
  setptr("vstart",	&vstart);
  setptr("vstop",	&vstop);
  setptr("vstep",	&vstep);
  setptr("fracnodes",	&fracnodes);
  setptr("axon_br_dia",	&axon_br_dia);
  setptr("varicos_dia",	&varicos_dia);
  setptr("passive",	&passive);
  setptr("set_drm",	&set_drm);
  setptr("celltype",	&celltype);

  // nvalfile = "nval_dsgc_cbp_bar.n";
  nvalfile = "nval_cbp_flash.n";
  // passive = 1;

}

/*-------------------------------------------------*/


void setparams(void)

  /*  set up default configuration for sb expts */
  /* cones, cone bipolars, sb, dsgc */

{
   int i;
   double zmax, zmin;

  if (!notinit(celltype)) {
     ct = find_ct(celltype);
     make_ct(ct);                       /* make the cell type */
  }
  else {
     ct = gca;
     make_ct(ct);                       /* make ganglion cell only */
  }
  set_ncel(ct,1);			/* set one cell only */

   if (!notinit(set_drm)) {		/* user set default Rm */
      setn(ct,NRM,set_drm); /* set default Rm */
      drm = set_drm;
   }

   garrsiz = 2;
   if (!notinit(passive) && passive==1) {
       setn(ct,BIOPHYS,0);	/* turn off biophysical properties */
   } else {
      setn(ct,BIOPHYS,1);	/* turn on  biophysical properties */
   }
   setn(ct,MORPH,0);	/* ==0 -> set dbp1 morphology from file, def = "morph_bp" */
   setn(ct,NCOLOR,RCOLOR);	/* set cell display color from region */

   if (notinit(ath_dia)) ath_dia = 0.4;  /* diameter of initial thin axon segment */

   if (notinit(axon_br_dia)) axon_br_dia = 0.5;  /* default dia for axon branches */
   if (notinit(varicos_dia)) varicos_dia = 1.5;  /* default dia for varicosities */
   if (notinit(den_dia))     den_dia = 0.2;      /* default dia for dendrites */

   if (notinit(node_scale)) node_scale = -3.2;	 /* 3: nodenum, 0.04: small font */

   // Set default channel types in density file.
   // To change, check manual, and uncomment and modify these:
   //
   //   _NA  = _NA2;
   //   _KA  = _K3;
   //   _KH  = _K4;
   //   _KIR = _K5;


}

/*-------------------------------------------------*/

void setdens(void)
{
}

void runonexit (void)
{
  if (savefile[0]!=0)  unlink (savefile);
}

/*-------------------------------------------------*/

#define MAXNODES 10000

void runexpt(void)

{
     int n,cn,nruns,nod,nnodes,doflag,sign;
     double MHz,timeperrun,timepredur;
     double cpuspeed,eqtime;
     double i1, i2, v1, v2, Gin, Rin,rad;
     double x, y, z, xsoma, ysoma, zsoma;
     node *npnt;
     int nodenums[MAXNODES];

   timinc = 5e-6;
   ploti = 1e-4;
   crit = 1e-12;

   if (notinit(dst)) dst = 0.0001;
   if (notinit(dt))   dt = 0.5;
   if (notinit(dv))   dv = 0.001;

   if (notinit(vstart)) vstart = -0.07;
   if (notinit(vstop))   vstop = -0.07;
   if (notinit(vstep))   vstep = 0.005;

   cn = 1;		/* cell number */
   stimdur = dt;
   endexp = 3 * stimdur;

   for(nnodes=0,npnt=nodepnt; npnt=foreach(npnt,ct,cn,-1,NULL,NULL,&nod); npnt=npnt->next) {
        nodenums[nnodes++] = nod;
   }
   for(n=0; n<nnodes/2; n++) {			/* reverse order to make soma first */
	   int temp;
        temp = nodenums[n];
	nodenums[n] = nodenums[nnodes-n-1];
	nodenums[nnodes-n-1] = temp;
   }

   if (ninfo >= 3) {
         double vmin, vmax;
         double imin, imax, plsize;
         int c,plnum;
       vmax = -0.02;
       vmin = -0.09;
       //plot_v_nod(ct,cn,782,vmin,vmax,blue,"", -1, -1);
       //plot_v_nod(ct,cn,1514,vmin,vmax,green,"", -1, -1);
       plot_v_nod(ct,cn,soma,vmin,vmax,red,"", -1, -1);
       plot_i_nod(ct,cn,soma, imin=-10e-10,imax=1e-10, c=blue,"Isoma",plnum=1,plsize=1);
   };

   if (notinit(setxmin)) setxmin = 0;                    // set plot to start at 0
   if (notinit(predur)) predur = 0.5;
   simtime = 0 - predur;
   vclamp(nde(ct,cn,soma), vstart, simtime, predur);
   step (predur);
   //simtime = 0;

   printf ("#\n");
   printf ("# number of nodes = %d\n",nnodes);

   if (notinit(fracnodes))  fracnodes = (nnodes<100? 1.0 : 100.0/nnodes);
   printf ("# fracnodes       = %.3g\n",fracnodes);
  
   if (ninfo >= 1) {		// for large cells, calculate how long it will take
     nruns = abs(int((vstop-vstart)/vstep) + 1);
     timepredur = 0.25 * cumcomp * predur;  /* minutes timed by hand for dt=0.5 one run / cumcomp */
     timeperrun = 0.25 * cumcomp * dt;  /* minutes timed by hand for dt=0.5 one run / cumcomp */
     MHz = 3456;		     /*  for a cpu of this speed */
     cpuspeed = system_speed();
     if (cpuspeed>10e3) cpuspeed = 3200;
     printf ("# cpu speed = %.5g\n",cpuspeed);
     if (vstart==vstop) printf ("# Calculating Rin at Vm = %.3g\n",vstart);
     else               printf ("# Calculating Rin from Vm = %.3g to %.3g\n",vstart,vstop);
     printf ("# Est time per node per vstep: %8.4g minutes\n",
                                (timeperrun+timepredur)*MHz/cpuspeed);
     printf ("# Estimated total run time:   %8.4g minutes\n",
                                nnodes*fracnodes*(nruns*timeperrun+timepredur)*MHz/cpuspeed);
     printf ("#\n");
   };

	/* find Rin at nodes in the cell */

   if (ninfo >= 0) {
     printf ("#\n");
     printf ("# node V     Rin(ohms)  Gin(S)    dist(um)\n");
     printf ("#\n");
   }

   // save model voltages so they can be restored after each step
       
   set_run_on_exit(runonexit);                       // set to erase savefile on ^C
   sprintf (savefile,"cell_Rin%06d",getpid());       // add pid to file name
   savemodel (savefile);			     // save beginning state
   if (vstart <= vstop) sign = 1;
   else                 sign = -1;

   for (n= 0; n<nnodes; n++) {		/* for a fraction of nodes (including soma) */

     if (n<= 0) { nod = soma; doflag = 1; }
     else       { nod = nodenums[n]; doflag = (drand() < fracnodes); }

     if (doflag) {

       npnt = nde(ct,cn,nod);
       if (nod==soma) { xsoma = npnt->xloc; ysoma = npnt->yloc; zsoma = npnt->zloc; }

        // vclamp(npnt,vstart,simtime,stimdur*3);
        // step (stimdur*3);

       for (v1=vstart; (v1*sign)<=(vstop*sign+1e-6); v1 += vstep) {

         simtime = 0;
         v2 = v1 + dv;

         //nod = soma;
         vclamp(npnt,v1,simtime,stimdur);
         step (stimdur-dst);
         i1 = i(npnt);
         step (dst);
         vclamp(npnt,v2,simtime,stimdur);
         step (stimdur-dst);
         i2 = i(npnt);
         step (dst);
         Rin = (v2-v1) / (i2-i1);
         x = npnt->xloc - xsoma;
	 y = npnt->yloc - ysoma;
	 z = npnt->zloc - zsoma;
         rad = sqrt (x*x + y*y + z*z); 
         if (ninfo>=2) {
            printf ("v1 %-8.5g  v2 %-8.5g\n",v1,v2);
            printf ("i1 %-8.5g  i2 %-8.5g\n",i1,i2);
         }
	 if (Rin<=0) Gin = 1;
	 else        Gin = 1 / Rin;
         printf ("%4d %-6.3g  %-9.3g %-9.3g %-8.3g\n",nod,v1,Rin,Gin,rad);
      }
      restoremodel (savefile);
    }
  }
  unlink (savefile);
  savefile[0] = 0;
}
