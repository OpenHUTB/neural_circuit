/* Experiment dsgc_cbp_stim */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "stimfuncs.h"
#include "onplot_dsgc_movie.cc"

extern int make_one_dbp1;
extern int n_dsgc;

double light_inhib;
int dsgc_prefdir;
const char *stimtype;

int rec_ct;
int rec_cn;
int inhib_point;

double barwidth;
double minten;
double sinten;
double starttime;
double stimdur;
double endwait;
double sblur;
double istim;
double predur;
double currinj;
double sdia;
double sloc;
double iloc;
double iinten;
int ninhib;

/*--------------------------------------------------------*/

void defparams(void)
{
  defparams_dsgc_movie();
  defparams_onplot_movie();

  setptr("rec_cn",       &rec_cn);
  setptr("rec_ct",       &rec_ct);
  setptr("dsgc_prefdir", &dsgc_prefdir);
  setptr("light_inhib",  &light_inhib);
  setptr("inhib_point",  &inhib_point);

  setptr("minten",     &minten);
  setptr("sinten",     &sinten);
  setptr("starttime",  &starttime);
  setptr("istim",      &istim);
  setptr("stimdur",    &stimdur);
  setptr("endwait",    &endwait);
  setptr("sblur",      &sblur);
  setptr("predur",     &predur);
  setptr("stimtype",   &stimtype);
  setptr("currinj",    &currinj);
  setptr("sdia",       &sdia);
  setptr("sloc",       &sloc);
  setptr("iloc",       &iloc);
  setptr("iinten",     &iinten);
  setptr("ninhib",     &ninhib);
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
  gcdistnod  = 582;
  pickden[dsgc] = 0; //941;     /* pick one dendrite to allow connectio ns */

  if (notinit(make_one_dbp1)) make_one_dbp1 = 0;
  if (notinit(stimtype)) stimtype="cbp_spot";
  if (notinit(n_dsgc)) n_dsgc = 1;
  if(notinit(rec_ct)) rec_ct = dsgc;    /* type of cell to record from */
  if(notinit(rec_cn)) rec_cn=1;         /* cellnum to record from */

  if(notinit(dsgc_prefdir)) dsgc_prefdir=0;

  display_z(zmax=-23, zmin=-20);            /* exclude dsgc Off-arborization layer */

}

/*--------------------------------------------------------*/

#define CSIZ 80

void runexpt(void) 

{
    int c, ct, cn, cbp_cn, pl;
    int stimnod, inod;
    double start, dia, dur, dscale;
    double ixoff, iyoff;
    node *npnt;
    elem *e;
    synapse *s;
    photorec *p;
    chattrib *a;
    nattrib *napnt;
    char gbuf[CSIZ];

  if (!strcmp(dsgc_file,"ds1e")) {
    cbp_cn    = 68;
    gcdistnod = 582;
  }
  else if (dsgc_file=="beta8b") {
    gcdistnod = 775;
  };
  //fprintf(stderr,"stimtype = %s; cur = %g, stimnod=%d \n", stimtype, currinj, stimnod);

  if (notinit(currinj)) currinj=50e-12;
  if (notinit(sblur)) sblur = 10;
  if (notinit(starttime)) starttime = 0.02;
  if (notinit(stimdur)) stimdur=.25;	/* used for non-moving stimuli */
  if (notinit(endwait)) endwait = .05;

   Vmax = -0.00;
   Vmaxg = -0.00;
   Vmin = -0.08;
   stimnod = gcdistnod;

  if (!(strcmp(stimtype,"cclamp"))) {
    plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmaxg,blue,"", -1, -1);/* plot V at soma */
    plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,red,"", -1, -1);/* plot V at soma */
    plot_ca_nod(dsgc,1,stimnod,1e-6,cyan,"Ca_soma", -1, -1);
    if (notinit(predur)) predur=.00;
    cclamp(ndn(dsgc,1,stimnod), currinj, start=starttime, dur=stimdur);
    // vclamp(ndn(dsgc,1,stimnod), -0.06, start=starttime, dur=stimdur);
  }
  else if (!(strcmp(stimtype,"cbp_spot"))) {

   /* add light transducer to each bipolar cell */


   for(npnt=nodepnt; npnt=foreach(npnt,dbp1,-1,soma,NULL,&cn,NULL); npnt=npnt->next) {
     p = (photorec*)make_transducer(ndn(dbp1,cn,soma));
     p->xpos=npnt->xloc;
     p->ypos=npnt->yloc;
   }

   ixoff = 400;					/* x offset for inhibition */
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


   if (notinit(sdia))         sdia = 80;
   if (notinit(sloc))         sloc = 140;
   if (notinit(iloc))         iloc = 70;	/* location of inhibition */
   if (notinit(minten))     minten = -0.045;
   if (notinit(sinten))     sinten =  0.005;
   if (notinit(iinten))     iinten =  0.008;
 

   if (notinit(inhib_point)) inhib_point = 0; 
   inod = 1000;

   if (inhib_point) { 
     if (notinit(ninhib)) ninhib = 422;		/* on-path inhibition */
     //if (notinit(ninhib)) ninhib = 529;	/* on-path inhibition */
     //if (notinit(ninhib)) ninhib = 297;	/* on-path inhibition */
     //if (notinit(ninhib)) ninhib = 469;	/* outside spot inhibition */
     //if (notinit(ninhib)) ninhib = 582;	/* peripheral inhibition */
     make_sphere(loc(ndn(dbp1,inod), iloc,0,-25), dia=2);
     p = (photorec*)make_transducer(ndn(dbp1,inod));
     p->xpos = iloc;
     p->ypos = 0;
     s = make_synapse(ndn(dbp1,inod), ndn(dsgc,1,ninhib)); 
		s->maxcond=10e-10; s->ngain=3;
		s->nfilt2=1; s->timec2=makfiltarr(1,0,s->timec2,10);
		s->thresh=-0.045; s->vrev=-0.08;
		napnt=make_vesnoise(s);

     stim_spot (dia=10, iloc, 0, iinten,starttime,stimdur);
     gmove (.03, .96);
     gpen(blue);
     gcwidth(.02);
     sprintf (gbuf,"inhib: dsgc node %d\n",ninhib);
     gtext (gbuf);
   };

    stim_backgr(minten,start=0);		/* background */
    stim_spot(sdia, sloc,0,       sinten,starttime,stimdur);
    stim_spot(sdia, ixoff+iloc,0, sinten,starttime,stimdur);

    if (notinit(istim)) istim = 0; 
    if (istim != 0) cclamp (ndn(dsgc,1,soma), istim,start=0,dur=1);

				/* for ds1d on-layer */
				/* bp  68, n 582   inside spot, at tip */
				/* bp 101, n 529   inside spot */
				/* bp 265, n 469   outside spot, near */
				/* bp 41,  n 422   on dend to soma */
				/* bp 40,  n 297   on dend to soma */
				/* bp 369, n 99    on dend to soma */
				/* bp 241, n 1464  middle top of cell */
				/* bp 488, n 2328  left middle top of cell */
   if (disp & 32) {
     if (space_time) {  /* movie */ 
      plot_v_nod(ct=dsgc,cn=1,soma,Vmin,Vmaxg,c=1,"Vsoma",pl=10,0.35);
      plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,c=red,"Vtip",pl=10,0.35);
     };
   }
   else { 
    plot_v_nod(ct=dsgc,cn=1,soma,Vmin,0.03,c=1,"",pl=10,1);/* V at soma */
    plot_ca_nod(dsgc,1,soma,1e-6,cyan,"Ca_soma", -1, -1);
    plot_ca_nod(dsgc,1,gcdistnod,1e-6,yellow,"", -1, -1);
    plot_v_nod(ct=dbp1,cbp_cn,axtrm,   Vmin,Vmax,red,"", -1, -1);
    plot_v_nod(ct=dbp1,41,axtrm,       Vmin,Vmax,blue,"", -1, -1);
    plot_v_nod(ct=dsgc,cn=1,gcdistnod,Vmin,Vmaxg,c=red,"", -1, -1);
    plot_v_nod(ct=dsgc,cn=1,422,      Vmin,Vmaxg,c=green,"", -1, -1);
    plot_v_nod(ct=dsgc,cn=1,1464,     Vmin,Vmaxg,c=blue,"", -1, -1);
    plot_v_nod(ct=dsgc,cn=1,2328,     Vmin,Vmaxg,c=magenta,"",pl=10,1);
    plot_synrate_out(dbp1,cbp_cn,0,500,green);
    plot_synrate_out(dbp1,241,0,500,blue);
    plot_synrate_out(dbp1,inod,0,500,red);
   };
  }; /* cbp_spot */

  if (notinit(predur)) predur=0.00;


  /* run experiment */

  setxmin=0;
  simtime=0-predur;

  endexp=predur+starttime+stimdur+endwait;
  step (predur);
  step (starttime+stimdur+endwait);
}
