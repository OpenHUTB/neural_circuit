/* Experiment cone_hz */
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
double stimdur;
double taildur;
double sdia;
double stimtime;
double minten;
double scontrast;
double Rext; 

double gHemi;
double dvsha;
double dvsc;
double coneha_cond;
double stimvolt;
double clca_cond;  
double clcac_cond;  
double coneca_cond;
double khz_cond;
double ca_pump;
double dcrm;
double predur;
double set_sarea;
double Imax;
double Imin;
double stimcur;
int cone_vc;
int cone_stim;
int light_stim;
int set_cclamp;
int vnoise;

int rec_ct;
int rec_cn;

void ephap_init(void);

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("stimdur",   &stimdur);
  setptr("taildur",   &taildur);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("Rext",      &Rext);
  setptr("dvsha",     &dvsha);
  setptr("dvsc",      &dvsc);
  setptr("vnoise",    &vnoise);
  setptr("stimvolt",  &stimvolt);
  setptr("coneha_cond", &coneha_cond);
  setptr("clca_cond",  &clca_cond);
  setptr("clcac_cond",  &clcac_cond);
  setptr("coneca_cond", &coneca_cond);
  setptr("khz_cond",    &khz_cond);
  setptr("ca_pump",     &ca_pump);
  setptr("dcrm",	&dcrm);
  setptr("predur",	&predur);
  setptr("set_sarea",   &set_sarea);
  setptr("Imax",        &Imax);
  setptr("Imin",        &Imin);
  setptr("cone_vc",     &cone_vc);
  setptr("cone_stim",   &cone_stim);
  setptr("light_stim",  &light_stim);
  setptr("set_cclamp",  &set_cclamp);
  setptr("stimcur",     &stimcur);
  nvalfile = "nval_bphz.n";
  ephap_init();
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  make_cones= 1;        /* make cones, cbp, gc */
  make_ha   = 1;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 0;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 0;
  make_dsgc = 0;
  make_ha_dbp1 = 0;
  make_ha_hbp1 = 0;
  make_cone_hb = 0;

  // set lamcrit smaller to keep comps presynaptic to cones separate 

  lamcrit = 0.01;
  setn(ha,NCOLOR,RCOLOR);       /* set cell display color from region */

  //Set Cone and Horizontal Cell Resting Potentials
  
  vk = -0.082;
  if (notinit(dvsc))  dvsc  = -0.04;		// for dens_cone.n
  if (notinit(dvsha)) dvsha = -0.04;		// for dens_ha.n
  if (notinit(dcrm))   dcrm = 0.5e3;  
  cone_maxcond = 880e-12;  

  if (!notinit(set_sarea)) dcasarea = set_sarea; // shell area, changes [Ca]i given ICa.
  if (notinit(cone_vc))       cone_vc = 1; 
  if (notinit(cone_stim))     cone_stim = 0; 
  if (notinit(light_stim)) light_stim = 0; 
  if (notinit(set_cclamp)) set_cclamp = 0; 
  if (notinit(stimcur))       stimcur = -10e-12; 
}

void setdens(void)
{
  // Set Ephaptic Parameters
  //
  setsv (ha,SCOND,1, 0);	//Set a value from the synaptic parameter table
  // setsv (xcone, SCOND,5,getsv(xcone, SCOND, 5)); 
   if (!notinit(coneha_cond)) setsv (xcone, SCOND,5,coneha_cond); 

  if (!notinit(vnoise)) setsv (xcone, SVNOISE,5,vnoise); 
 
 // set some channel conductances if wanted from the command line
 
  if (!notinit(clca_cond))   celdens[xcone][0][_CLCA][AXOND] = clca_cond;  
  if (!notinit(clcac_cond))  celdens[xcone][0][_CLCAC][AXOND]= clcac_cond;  
  if (!notinit(coneca_cond)) celdens[xcone][0][_CA]  [AXOND] = coneca_cond;  
  if (!notinit(ca_pump))     celdens[xcone][0][_CAP] [AXOND] = ca_pump;  
  //if (!notinit(khz_cond))    celdens[ha][0] [_KHZ] [SOMA] = khz_cond;  

  //setsv (hb,SCOND,1, 0);
  if(notinit(rec_ct)) rec_ct = gca;
  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 1e3;      /* background light intensity */
}

/*------------------------------------------------------*/
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

/*------------------------------------------------------*/

// double pnx_func(double atp) 
// void ephapsynapnode(void)
// void readephapticfeedback(void)
// void clca_speed (int stype, double tau)

#include "conehz_ephap_inc.cc"

/*------------------------------------------------------*/

void onplot(void)
{
    onplot_ephap();
}

/*------------------------------------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, pct, pcn, plnum, cn1, pa, pb;
    int colr, nc, pl, pves, prate;
    int midha, midcbp;
    double t, st, fmax,fmin, vol;
    double rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax, phmin, phmax;
    //double volinc = 0.01; //VCLAMP 1
    elem *epnt;

  midcone = findmid(xcone,0,0);
  midcone = 25;
  midha   = findmid(ha, 0,0);
  midcbp  = findmid(dbp1, 0,0);
  surrcone = 22;

  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d\n",  midcbp);
  if (ninfo >=1) fprintf (stderr,"# mid ha  # %d\n",  midha);  

  if(notinit(Rext))  Rext  = 75e6;
  if(notinit(gHemi)) gHemi = 5e-9;
  if(notinit(gjmod)) gjmod = 0;
 
  ephapsynapnode();
  setonplot(onplot);
    
  // if (ninfo >=1) fprintf (stderr,"# midcone %d htipcent %d\n",  midcone, htipcent);  
  // if (ninfo >=1) fprintf (stderr,"# surrcone %d htipsurr %d\n", surrcone, htipsurr);  

  if (notinit(temp_freq)) temp_freq = 4;
  // fprintf(stderr, "# temp freq %d\n", temp_freq);
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  

  if (notinit(stimtime))   stimtime = 0.20;  /* stimulus time */
  if (notinit(stimdur))     stimdur = 0.30;  /* stimulus duration */
  if (notinit(taildur))     taildur = 0.20;  /* tail current duration */
  if (notinit(sdia))           sdia = 300;  /* spot diameter */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0.9;   /* stimulus contrast */

  exptdur = stimtime + stimdur + taildur;
  endexp  = exptdur;


  plot_v_nod(ct=xcone,cn=midcone,n=axtrm,Vmin=-.055,Vmax =-.030,colr=red,"Vcone_cent", -1, -1); 
  plot_v_nod(ct=xcone,cn=surrcone,     n=axtrm,Vmin=-.055,Vmax =-.030,colr=green,"Vcone_surr", -1, -1); 
  plot_synrate_out(ct=xcone,cn=midcone,1,ha,1,rmin=0,rmax=500,fmax=1,colr=magenta,prate=1,pves=0,"cone_cent",1.0); 
  plot_synrate_out(ct=xcone,cn=surrcone,1,     ha,1,rmin=0,rmax=500,fmax=1,colr=blue,   prate=1,pves=0,"cone_surr",1.0); 

  Imax = 10e-12;
  Imin = -30e-12;

  plot_var(&currentampa,1,Imax,Imin);		        /* plot AMPA chan current */
  plot_param ("IAMPA_cent", colr=7, pl=8, plsize=1);
   
  plot_var(&currentampas,1,Imax,Imin);		        /* plot AMPA chan current */
  plot_param ("IAMPA_surr", colr=9, pl=8, plsize=1);
   
//  plot_var(&currentca,1,Imax,Imin);		        /* plot hemichannel current */
//  plot_param ("ICa_cent", colr=1, pl=8, plsize=1);
   
  plot_var(&currenthemi,1,Imax,Imin);		        /* plot hemichannel current */
  plot_param ("IHemi_cent", colr=6, pl=8, plsize=1);
   
//  plot_var(&currenthemis,1,Imax,Imin);		        /* plot hemichannel current */
//  plot_param ("IHemi_surr", colr=4, pl=8, plsize=1);
   
  plot_var(&currentclca,1,Imax/10,Imin/10);		        /* plot ClCa current */
  plot_param ("IClCa_cent", colr=5, pl=8, plsize=1);
   
  plot_var(&currentclcas,1,Imax/10,Imin/10);		        /* plot ClCa current */
  plot_param ("IClCa_surr", colr=2, pl=8, plsize=1);
   
  plot_var(&voffset_midcone,1,0, -0.02);	        	/* plot voffset */
  plot_param ("Vext_cent", colr=4, pl=6, plsize=1);

  plot_var(&voffset_surrcone,1,0, -0.02);		       	/* plot voffset */
  plot_param ("Vext_surr", colr=2, pl=6, plsize=1);

  plot_v_nod(ct=ha,cn=midha,n=soma,Vmin=-.08,Vmax=0,colr=green,"Vhz_soma",5, 1.5); /* plot Vha*/
  //plot_v_nod(ct=ha,cn=midha,n=htipcent,Vmin=-.07,Vmax=0,colr=ltgreen,"",5,1.5); /* plot Vha*/
  //plot_v_nod(ct=ha,cn=midha,n=10,Vmin=-.07,Vmax=0,colr=blue,"",         5, 1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=htipcent,Vmin=-.08,Vmax=0,colr=ltgreen,"Vhz_tip_cent", 5, 1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=htipsurr,Vmin=-.08,Vmax=0,colr=ltmag,"Vhz_tip_surr",   5, 1.5); /* plot Vha*/
  if (ph_gain != 0 || phh_gain != 0) {
   plot_v_nod(ct=ha,cn=midha,n=htipcent+nodeph,Vmin=-.01,Vmax=0.02,colr=ltgreen,"Vhz_tip_centph",   4, 1.5); /* plot Vha*/
   plot_v_nod(ct=ha,cn=midha,n=htipsurr+nodeph,Vmin=-.01,Vmax=0.02,colr=ltmag  ,"Vhz_tip_surrph",   4, 1.5); /* plot Vha*/

   plot_ph_nod(ct=ha,cn=midha,n=htipcent+node1,phmin=7,phmax=8,colr=ltgreen,"pH_tip_centph",   3, 1.5); /* plot Vha*/
   // plot_ph_nod(ct=ha,cn=midha,n=htipsurr+node1,phmin=7,phmax=8,colr=ltmag,  "pH_tip_surrph",   3, 1.5); /* plot Vha*/
  }
  // plot_v_nod(ct=ha,cn=midha,n=htipcent+nodeatp,Vmin=-.06,Vmax=0.02,colr=blue,"Vhz_tip_centatp",   2, 1.5); 
  plot (ATP,ndn(ha,midha,htipcent+nodeatp), 1e-6, 0); plot_param("ATP", red, 2, 1.);

  // plot_i_nod(ct=ha, cn=1, n=soma, -2e-9, 1e-9, colr=blue, "",   -1,-1.5);
  // plot_synrate_out(ct=ha,cn=midha,pct=xcone,pcn=midcone,rmin=0,rmax=5000,colr=yellow,1); 
  // plot_v_nod(ct=dbp1,cn=midcbp,n=soma,Vmin=-.045,Vmax =-.040,colr=red,"", -1, -1);
  // plot_synrate_out(ct=dbp1,cn=midcbp,rmin=0,rmax=200,colr=magenta);
  // plot_v_nod(ct=gca,cn=1,n=soma,Vmin=-.075,Vmax =-.055,colr=blue,"", -1, -1);
  //if (getn(gca,BIOPHYS)) { plot(CA,1,ndn(gca,1,soma),fmax=0.5e-6,fmin=0); 
  //				plot_param("Cai", plnum=0,plsize=0.3);}

 // if (getn(xcone,BIOPHYS)) plot(CA,1,ndn(xcone,midcone,1),fmax=3.0e-6,fmin=0); 
 //  				plot_param("Cai coneterm",magenta, plnum=0,plsize=1.0);
 // if (getn(xcone,BIOPHYS)) plot(CA,1,ndn(xcone,surrcone,1),fmax=3.0e-6,fmin=0); 
 //  				plot_param("Cai conesurr", blue,plnum=0,plsize=1.0);
 
 // if (getn(xcone,BIOPHYS)) plot(CA,100,ndn(xcone,midcone,1),fmax=3.0e-6,fmin=0); 
 // 				plot_param("Caicore", plnum=0,plsize=1.0);
  //if (getn(xcone,BIOPHYS)) plot(CA,-1,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
  // 				plot_param("Cao coneterm", plnum=0,plsize=1.0);
  //if (getn(xcone,BIOPHYS)) plot(CA,-100,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
  // 				plot_param("Caocore coneterm", plnum=0,plsize=1.0);

   stim_backgr(minten);
   if (notinit(stimvolt)) stimvolt = -0.045; 

//   if (notinit(setxmin)) setxmin = 0;                    // set plot to start at 0
   if (notinit(predur)) predur = 0.05;			// equilib time before simulation starts
//   simtime = 0 - predur;

  if (light_stim) {
    for (n=1; n<=ncones; n++) {
      if (n!=midcone) {
	  stim_cone (ndn(xcone,n,0), 1e5, stimtime, stimdur, 1);
      } 
    }
  } else {
    for(epnt=elempnt; epnt = foreach (epnt, CONE, xcone, -1, &pa, &pb); epnt=epnt->next){
      //fprintf (stderr,"cone # %d\n",pb);
      if (pb!=midcone) {
         if (set_cclamp) cclamp(ndn(xcone,pb,axtrm), stimcur, stimtime, stimdur );
         else            vclamp(ndn(xcone,pb,axtrm), stimvolt, stimtime, stimdur );
      }
      else {
	 if (cone_vc)   vclamp(ndn(xcone,pb,axtrm), -0.040, simtime, 2 );
	 else if (cone_stim) {
            if (set_cclamp) cclamp(ndn(xcone,pb,axtrm), stimcur, stimtime, stimdur );
            else            vclamp(ndn(xcone,pb,axtrm), stimvolt, stimtime, stimdur );
	 }
      }
    }
  }

   clca_speed (1,2);		// speed up channels in cleft
   clca_speed (2,3);		// set slow channels faster
   // for (st=0; st<predur; st+= stiminc){
   //     readephapticfeedback();
   //     step(stiminc);
   // } 
   step (predur);
   clca_speed (2,0.333);	// set slow channels slower again

   // for (st=0; st<exptdur; st+= stiminc){
   //     readephapticfeedback();
   //     step(stiminc);
   // } 
   step (exptdur);
}


