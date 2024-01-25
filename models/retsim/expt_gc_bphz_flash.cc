/* Experiment gc_bphz_flash */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"


double temp_freq;
double ntrials;
double dstim;
double sdia;
double stimtime;
double minten;
double scontrast;
double Rext; 
double gHemi;
double dvsha;
double dvsc;
double dcrm;
double coneha_cond;
double stimvolt;
double predur;

int vnoise;

int rec_ct;
int rec_cn;

/*------------------------------------------------------*/

#define EPHAP 1

double voffset_midcone;
double currentampa;
double currenthemi;
double currentca;
double currentclca;
double currentclcac;
double sumcurrent;

double voffset_surrcone;
double currentampas;
double currenthemis;
double currentcas;
double currentclcas;
double currentclcacs;
double sumcurrents;
int midcone, surrcone, htipcent, htipsurr;

// void ephapsynapnode(void)
// void readephapticfeedback(void)
// void clca_speed (int stype, double tau)


#ifdef EPHAP
  #include "conehz_ephap_inc.cc"
#endif

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("predur",    &predur);
  setptr("Rext",      &Rext);
  setptr("gHemi",     &gHemi);
  setptr("dvsha",     &dvsha);
  setptr("dvsc",      &dvsc);
  setptr("dcrm",      &dcrm);
  setptr("vnoise",    &vnoise);
  setptr("stimvolt", &stimvolt);
  setptr("coneha_cond", &coneha_cond);
  nvalfile = "nval_gc_bphz.n";
  #ifdef EPHAP
    ephap_init();
  #endif
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  make_cones= 1;        /* make cones, cbp, gc */
  make_ha   = 1;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 1;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 1;
  make_dsgc = 0;
  make_ha_dbp1 = 0;
  make_ha_hbp1 = 0;
  make_cone_hb = 0;

  make_ha_ha = 0;

  // set lamcrit smaller to keep comps presynaptic to cones separate 
  //  lamcrit = 0.01;
 
  // setn(ha,NCOLOR,RCOLOR);       /* set cell display color from region */
  setn(gca,NCOLOR,RCOLOR);       /* set cell display color from region */

  //Set Cone and Horizontal Cell Resting Potentials
  
  vk = -0.082;
  if (notinit(dvsc))  dvsc  = -0.035;
  if (notinit(dvsha)) dvsha = -0.04;
  if (notinit(dcrm))  dcrm  = 1.0e3;
  
  // Set Ephaptic Parameters
  //
  setsv (ha,SCOND,1, 0);        // turn off standard ha->cone synapse in nval file, conductance = 0

  // setsv (xcone, SCOND,5,getsv(xcone, SCOND, 5)); 
   if (!notinit(coneha_cond)) setsv (xcone, SCOND,5,coneha_cond); 

  if (!notinit(vnoise)) setsv (xcone, SVNOISE,5,vnoise); 
   
  if(notinit(rec_ct)) rec_ct = gca;
  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 50e4;      /* background light intensity */

  if(notinit(Rext))  Rext  = 75e6;
  if(notinit(gHemi)) gHemi = 5e-9;
}

/*------------------------------------------------------------------------------------*/


#ifdef EPHAP 
void addsyns(void)
{
   sethzconns();
}
#endif

/*------------------------------------------------------------------------------------*/

#ifdef EPHAP
void onplot(void) 
{
    onplot_ephap();
}
#endif

/*------------------------------------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, pct, pcn, plnum, cn1, pa, pb;
    int colr, pl;
    int midha, midcbp;
    double t, st, fmax,fmin, vol;
    double rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax;
    //double volinc = 0.01; //VCLAMP 1
    elem *epnt;

  stiminc = 2e-6;

  midcone = findmid(xcone,0,0);
  midha   = findmid(ha, 0,0);
  midcbp  = findmid(dbp1, 0,0);
  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d\n",  midcbp);
  if (ninfo >=1) fprintf (stderr,"# mid ha  # %d\n",  midha);  

  // if (ninfo >=1) fprintf (stderr,"# midcone %d htipcent %d\n", midcone, htipcent);  
  // if (ninfo >=1) fprintf (stderr,"# surrcone %d htipsurr %d\n", surrcone, htipsurr);  

  if (notinit(temp_freq)) temp_freq = 4;
  // fprintf(stderr, "# temp freq %d\n", temp_freq);
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  
  dtrial = 1 / temp_freq;
  exptdur = dtrial * ntrials;
  endexp  = exptdur;


  if (notinit(dstim))         dstim = .10;  /* stimulus duration */
  if (notinit(sdia))           sdia = 300;  /* spot diameter */
  if (notinit(stimtime))   stimtime = .10;  /* stimulus time */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0.5;   /* stimulus contrast */
  if (notinit(predur))       predur = 0.05;     /* equilib time */

  simtime = -predur;            /* start model before time 0 */
  setxmin = 0;                  /* start plotting at time 0 */

  #ifdef EPHAP 
    ephapsynapnode();
    setonplot(onplot);
  #endif
  
  plot_v_nod(ct=xcone,cn=midcone,n=axtrm,Vmin=-.055,Vmax =-.030,colr=cyan,"", -1, -1); 
  plot_v_nod(ct=xcone,cn=11,     n=axtrm,Vmin=-.055,Vmax =-.030,colr=green,"", -1, -1); 
  plot_synrate_out(ct=xcone,cn=midcone,ha,1,rmin=0,rmax=250,colr=magenta, 1); /* plot rate out */
  plot_synrate_out(ct=xcone,cn=11,     ha,1,rmin=0,rmax=250,colr=blue, 1); /* plot rate out */

  plot_var(&currentampa,1,1e-12,-10e-12);		        /* plot AMPA chan current */
  plot_param ("IAMPA", colr=7, pl=6);
   
  plot_var(&currenthemi,1,1e-12,-10e-12);		        /* plot hemichannel current */
  plot_param ("IHemi", colr=6, pl=6);
   
  plot_var(&currentclca,1,1e-12,-10e-12);		        /* plot ClCa current */
  plot_param ("IClCa", colr=5, pl=6);
   
  plot_var(&voffset_midcone,1,0, -0.02);		       	/* plot voffset_midcone */
  plot_param ("Vext", 4, 4);

  plot_v_nod(ct=ha,cn=midha,n=soma,Vmin=-.07,Vmax=0,colr=green,"", -1, 1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=htipcent,Vmin=-.07,Vmax=0,colr=ltgreen,"",-1,1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=10,Vmin=-.07,Vmax=0,colr=blue,"",    -1, 1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=11,Vmin=-.07,Vmax=0,colr=ltmag,"",   -1, 1.5); /* plot Vha*/
  // plot_i_nod(ct=ha, cn=1, n=soma, -2e-9, 1e-9, colr=blue, "",   -1,-1.5);
  // plot_synrate_out(ct=ha,cn=midha,pct=xcone,pcn=midcone,rmin=0,rmax=5000,colr=yellow,1); 
  // plot_v_nod(ct=dbp1,cn=midcbp,n=soma,Vmin=-.045,Vmax =-.040,colr=red,"", -1, -1);
  // plot_synrate_out(ct=dbp1,cn=midcbp,rmin=0,rmax=200,colr=magenta);
  // plot_v_nod(ct=gca,cn=1,n=soma,Vmin=-.075,Vmax =-.055,colr=blue,"", -1, -1);
  //if (getn(gca,BIOPHYS)) { plot(CA,1,ndn(gca,1,soma),fmax=0.5e-6,fmin=0); 
  //				plot_param("Cai", plnum=0,plsize=0.3);}
  if (getn(xcone,BIOPHYS)) plot(CA,1,ndn(xcone,midcone,1),fmax=2.0e-6,fmin=0); 
  				plot_param("Cai coneterm", plnum=0,plsize=0.5);

   if (notinit(stimvolt)) stimvolt = -0.045; 
   stim_backgr(minten);
    for(epnt=elempnt; epnt = foreach (epnt, CONE, xcone, -1, &pa, &pb); epnt=epnt->next){
      //fprintf (stderr,"cone # %d\n",pb);
      if (pb!=midcone) {
      //if (pb!=1000) {
        vclamp(ndn(xcone,pb,soma), -0.040, 0, stimtime );
        vclamp(ndn(xcone,pb,soma),  stimvolt, stimtime, dstim );
        vclamp(ndn(xcone,pb,soma), -0.040, stimtime+dstim,1 );
      }
      else vclamp(ndn(xcone,pb,soma), -0.040, 0, 2 );
    };
    //vclamp(ndn(ha,midha,soma), -40, 0.1, 0.1); //VCLAMP - 5 - HORZONTAL CELL

   for (t=0; t<exptdur; t+= dtrial){
        double start, dur;
     simtime = 0;
      stim_spot(sdia, 0, 0, minten*scontrast, start=t+stimtime,dur=dstim);
      //vclamp(ndn(ha,midha,soma), vol, 0.1, 0.1); //VCLAMP - 5 - HORZONTAL CELL
       step(dtrial);
   }
}


