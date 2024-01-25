/* Experiment gc_Rin */

#include <stdio.h>
#include <math.h>
#include "stimfuncs.h"
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

int light_inhib;
int dsgc_prefdir;

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
  nvalfile = "nval_dsgc_cbp_bar.n";
}

/*-------------------------------------------------*/


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
  make_ams   = 0;
  make_sbac  = 0;
  make_gca     = 0;
  make_gcb     = 0;
  make_gcaoff  = 0;
  make_gcboff  = 0;
  make_dsgc    = 0;

  setn(dsgc,BIOPHYS,0);			/* turn off biophysical properties */

  if (notinit(ath_dia)) ath_dia = 0.4;  /* diameter of initial thin axon segment */

  if (n_dsgc<0) n_dsgc = 1;
  
}

/*-------------------------------------------------*/

#define MAXNODES 10000

void runexpt(void)

{
     int n,cn,ct,nruns,nod,nnodes,doflag;
     double MHz,timeperrun;
     double cpuspeed,eqtime;
     double i1, i2, v1, v2, Rin,rad;
     double x, y, xsoma, ysoma;
     node *npnt;
     int nodenums[MAXNODES];

   if (notinit(dst)) dst = 0.0001;
   if (notinit(dt))   dt = 0.025;
   if (notinit(dv))   dv = 0.002;

   if (notinit(vstart)) vstart = -0.06;
   if (notinit(vstop))   vstop = -0.06;
   if (notinit(vstep))   vstep = 0.005;

   if (notinit(fracnodes))  fracnodes = 0.2;

   if       (n_dbp1 > 0)  ct = dbp1;
   else if  (n_dbp2 > 0)  ct = dbp2;
   else if  (n_hbp1 > 0)  ct = hbp1;
   else if  (n_hbp2 > 0)  ct = hbp2;
   else if  (n_gca > 0)   ct = gca;
   else if  (n_gcb > 0)   ct = gcb;
   else if  (n_gcaoff > 0) ct = gcaoff;
   else if  (n_gcboff > 0) ct = gcboff;
   else if  (n_dsgc > 0)   ct = dsgc;

   cn = 1;		/* cell number */
   stimdur = dt;
   endexp = 3 * stimdur;

   for(nnodes=0,npnt=nodepnt; npnt=foreach(npnt,ct,cn,-1,NULL,NULL,&nod); npnt=npnt->next) {
        if (nod != soma) nodenums[nnodes++] = nod;
   }

   fprintf (stderr,"# number of nodes = %d\n",nnodes);

   if (ninfo >= 0) {
     nruns = int((vstop-vstart)/vstep + 1);
     timeperrun = 0.05;			/* minutes timed by hand for 0.025 dt */
     MHz = 2600;			/*  for a cpu of this speed */
     cpuspeed = system_speed();
     fprintf (stderr,"# cpu speed = %8.5g\n",cpuspeed);
     fprintf (stderr,"# Calculating Rin from Vm = %.3g to %.3g\n",vstart,vstop);
     fprintf (stderr,"# Estimated run time: %8.5g minutes\n",
                                nnodes*nruns*timeperrun*MHz/cpuspeed);
     fprintf (stderr,"#\n");
   };

   if (ninfo >= 2) {
         double vmin, vmax;
         double imin, imax, plsize;
         int c,plnum;
       vmax = -0.02;
       vmin = -0.09;
       //plot_v_nod(ct,cn,782,vmin,vmax,blue,"", -1, -1);
       //plot_v_nod(ct,cn,1514,vmin,vmax,green,"", -1, -1);
       plot_v_nod(ct,cn,soma,vmin,vmax,red,"", -1, -1);
       plot_i_nod(ct,cn,soma, imin=-5e-9,imax=5e-9, c=blue,"Isoma",plnum=1,plsize=1);
   };

   eqtime = .005;
   vclamp(nde(ct,cn,soma), vstart, simtime, eqtime-dst);
   step (eqtime);
   //simtime = 0;

	/* find Rin at nodes in the cell */

   if (ninfo >= 0) {
     fprintf (stderr,"#\n");
     fprintf (stderr,"#  node   V   Rin(ohms) dist(um)\n");
     fprintf (stderr,"#\n");
   }

   for (n= -1; n<nnodes; n++) {		/* for a fraction of nodes (including soma) */

     if (n== -1)  { nod = soma; doflag = 1; }
     else         { nod = nodenums[n]; doflag = (drand() < fracnodes); }

     if (doflag) {

       npnt = nde(ct,cn,nod);
       if (nod==soma) { xsoma = npnt->xloc; ysoma = npnt->yloc; }

       for (v1=vstart; v1<=vstop+dv; v1+=vstep) {

         simtime = 0;
         v2 = v1 + dv;

         //nod = soma;
         vclamp(npnt,v1,simtime,stimdur);
         step (stimdur);
         i1 = i(npnt);
         step (dst);
         vclamp(npnt,v2,simtime,stimdur);
         step (stimdur);
         i2 = i(npnt);
         step (dst);
         Rin = (v1-v2) / (i1-i2);
         x = npnt->xloc - xsoma;
	 y = npnt->yloc - ysoma;
         rad = sqrt (x*x + y*y); 
         if (ninfo>=2) {
            fprintf (stderr,"v1 %-6.3g  v2 %-6.3g\n",v1,v2);
            fprintf (stderr,"i1 %-6.3g  i2 %-6.3g\n",i1,i2);
         }
         fprintf (stderr,"%4d %-6.3g  %-9.3g %-8.3g\n",nod,v1,Rin,rad);
      }
    }
  }
}
