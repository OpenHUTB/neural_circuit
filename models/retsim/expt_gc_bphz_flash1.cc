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
double vtime; //vclamp
int vnoise;

int rec_ct;
int rec_cn;

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
  setptr("Rext",      &Rext);
  setptr("gHemi",     &gHemi);
  setptr("dvsha",     &dvsha);
  setptr("dvsc",      &dvsc);
  setptr("dcrm",      &dcrm);
  setptr("vnoise",    &vnoise);
  setptr("stimvolt", &stimvolt);
  setptr("coneha_cond", &coneha_cond);
  nvalfile = "nval_bphz.n";
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
 
  // set lamcrit smaller to keep comps presynaptic to cones separate 

  lamcrit = 0.01;
  setn(ha,NCOLOR,RCOLOR);       /* set cell display color from region */

  //Set Cone and Horizontal Cell Resting Potentials
  
  vk = -0.082;
  if (notinit(dvsc))  dvsc  = -0.035;
  if (notinit(dvsha)) dvsha = -0.04;
  if (notinit(dcrm))  dcrm  = 1.0e3;
  
  // Set Ephaptic Parameters
  //
  setsv (ha,SCOND,1, 0);	//Set a value from the synaptic parameter table
  // setsv (xcone, SCOND,5,getsv(xcone, SCOND, 5)); 
   if (!notinit(coneha_cond)) setsv (xcone, SCOND,5,coneha_cond); 

  if (!notinit(vnoise)) setsv (xcone, SVNOISE,5,vnoise); 
   
  //setsv (hb,SCOND,1, 0);
  if(notinit(rec_ct)) rec_ct = gca;
  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 50e4;      /* background light intensity */
}

/*------------------------------------------------------*/
double voffset;
double voffset_midcone;
double currentampa;
double currenthemi;
double currentclca;
int midcone, htip;

double plotvoffset(double val, double tim)
{
  return voffset_midcone;
}
/*------------------------------------------------------*/

double plotcurrentampa(double val, double tim)
{
  return currentampa;
}
/*------------------------------------------------------*/

double plotcurrenthemi(double val, double tim)
{
  return currenthemi;
}
/*------------------------------------------------------*/

double plotcurrentclca(double val, double tim)
{
  return currentclca;
}
/*------------------------------------------------------*/
int node1 = 5000;
int node2 = 6000;


void ephapsynapnode(void)
{
    int pa, pb, ct1, cn1, cd1, ct2, cn2, cd2;
    synapse *sepnt = NULL;
    elem *epnt = NULL;
    resistor *r1,*r2;
    loadelem *re = NULL;
    capac *ce = NULL;
    synap *spnt = NULL;
    conlst *c = NULL;
    chan *c1,*c2,*c3;
    comp *rext_pnt = NULL;
    node *nd = NULL;

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     if (sepnt->node2a==ha || sepnt->node2a==hb) {  /* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
        ct2 = sepnt->node2a;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	r1 = (resistor *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2, RESISTOR);   
	r1->r = 1/gHemi;
	r2 = (resistor *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2+node2, RESISTOR);   
	r2->r = Rext/2;
	re = (loadelem *)at(ct2, cn2, cd2+node2, LOAD);   
	re->r = Rext/2;
	re->vrev = re->vrest = 0;
	ce = (capac *)at(ct2, cn2, cd2+node1, GNDCAP);
	ce->c = 2e-12;
	ce = (capac *)at(ct2, cn2, cd2+node2, GNDCAP);
	ce->c = 2e-12;
	if (cn1==midcone) htip = cd2;
      }  /* if (sepnt->node2a==ha) */

   }  /* foreach (epnt, SYNAPSE) */

   initsim();		//  generate compartments

   /* add extracellular compartment from Rext to capture channel currents */

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     if (sepnt->node2a==ha || sepnt->node2a==hb) {	/* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
        ct2 = sepnt->node2a;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	if ((nd=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){		/* get conn list for presyn compartment */
	       if ((c1=findchan(c,CA,-1))!=NULL) {	/* find channel in cone terminal */
	         addchan_extern(rext_pnt, c1);		/* add external comp to channel */
	       }
	       if ((c2=findchan(c,ClCa,-1))!=NULL) {	/* find channel in cone terminal */
	         addchan_extern(rext_pnt, c2);		/* add external comp to channel */
	       }
	    }  /* if (comp1->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */
	if ((nd=ndn(ct2, cn2, cd2+node2))!=NULL) {	/* get Rext node for synapse */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp2->clst)!=NULL){		/* get conn list for postsyn compartment */
	       if ((c3=findchan(c,AMPA,-1))!=NULL) {	/* find channel in Hz dendr tip */
	         addchan_extern(rext_pnt, c3);		/* add external comp to channel */
	       }
	    }  /* if (comp2->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */
      }  /* if (sepnt->node2a==ha) */
   }  /* foreach (epnt, SYNAPSE) */
}

/*----------------------------------------------------------------------------------------*/

void readephapticfeedback(void)

/* connect cone (ct1) to horizontal cell (ct2) with ephaptic feedback */

{
  int pa, pb, pc, pd, ct1, cn1, cd1, ct2, cn2, cd2;
  synapse *sepnt = NULL;
  synap *spnt = NULL;
  elem *epnt = NULL;
  int pct, pcn;
  node *nd = NULL;
  capac *ce = NULL;
  conlst *c = NULL;
  chan *c1 = NULL;
  chan *c2 = NULL;
  cable *cb1 = NULL;
  attrib *apnt = NULL;
  double dtrial = 1;
  double t, st, fmax,fmin, vol;
  double gs,vs,rs, vcone, vhz;
  int midha, midcbp;
  double rmin, rmax, plsize;
  struct conn *r1pnt = NULL;
  struct conn *r2pnt = NULL;
  struct conn *r3pnt = NULL;
  double R1, R2, R;
  double gclca=0, gca=0;

       for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
	 sepnt = (synapse *)epnt;
	 if (sepnt->node2a==ha || sepnt->node2a==hb) {  /* if the postsynaptic cell is ha or hb */
            ct1 = sepnt->node1a;
	    cn1 = sepnt->node1b;
	    cd1 = sepnt->node1c;
            ct2 = sepnt->node2a;
	    cn2 = sepnt->node2b;
	    cd2 = sepnt->node2c;
	    if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* run-time synapse struct */
	      gs = record_chan(((elem*)sepnt), G, 0);
	      // fprintf(stderr, "ampacond= %g\n", gs);
	      vcone = spnt->comp1->v;
	      vhz   = spnt->comp2->v;
	      // fprintf(stderr, "vhz %g \n" , vhz);
	      if (spnt->spre) spnt = spnt->spre->sdyad;

	      if ((nd=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node */
	         voffset = nd->comptr->v; 
	           if ((c=nd->comptr->clst)!=NULL){
		      if ((c1=findchan(c,CA,-1))!=NULL) {
		        gca = c1->conduct;
		        // fprintf (stderr,"ctype1 %d %g\n",c1->ctype,c1->conduct);
		      }
		      if ((c2=findchan(c,ClCa,-1))!=NULL) {
		        gclca = c2->conduct;
		      }
	           } /* if (comp1->clst) */

		   if (cn1==midcone) {
		       // fprintf(stderr, "midcone %d  ha tip %d\n" , midcone, cd2);
		       // fprintf(stderr, "voltagedropext %g \n" , voffset);
		       currentampa = (vhz - voffset) * (gs);
		       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
		       currenthemi = (vhz - voffset) * gHemi;
		       // fprintf(stderr, "currenthemi %g \n" , currenthemi);
		       currentclca = (vcone - voffset - vcl) * (gclca);
		       // fprintf (stderr,"gclca %g\n",gclca);
		       // fprintf(stderr, "currentclca %g \n" , currentclca);
		       voffset_midcone = voffset;
		   } /* if (cn1==midcone) */

	        }  /* if (nd=ndn()) */
	    }  /* if ((spnt=sepnt->lptr) */
         }   /* if sepnt */
       }    /* for (epnt;;) */
}
       
/*------------------------------------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, pct, pcn, plnum, cn1, pa, pb;
    int colr, pl;
    int midha, midcbp;
    double t, st, fmax,fmin, vol;
    double rmin, rmax, plsize, voffset;
    double dtrial,exptdur;
    double Vmin, Vmax;
    //double volinc = 0.01; //VCLAMP 1
    elem *epnt;

  midcone = findmid(xcone,0,0);
  midha   = findmid(ha, 0,0);
  midcbp  = findmid(dbp1, 0,0);
  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d\n",  midcbp);
  if (ninfo >=1) fprintf (stderr,"# mid ha  # %d\n",  midha);  

  if(notinit(Rext))  Rext  = 75e6;
  if(notinit(gHemi)) gHemi = 5e-9;
 
  ephapsynapnode();
    
  // if (ninfo >=1) fprintf (stderr,"# midcone %d htip %d\n", midcone, htip);  

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
   vtime = 0.0; // VCLAMP 2



  plot_v_nod(ct=xcone,cn=midcone,n=axtrm,Vmin=-.055,Vmax =-.030,colr=cyan,"", -1, -1); 
  plot_v_nod(ct=xcone,cn=11,     n=axtrm,Vmin=-.055,Vmax =-.030,colr=green,"", -1, -1); 
  plot_synrate_out(ct=xcone,cn=midcone,ha,1,rmin=0,rmax=250,colr=magenta, 1); /* plot rate out */
  plot_synrate_out(ct=xcone,cn=11,     ha,1,rmin=0,rmax=250,colr=blue, 1); /* plot rate out */

  plot_func(plotcurrentampa,1,1e-12,-10e-12);		        /* plot AMPA chan current */
  plot_param ("IAMPA", colr=7, pl=6);
   
  plot_func(plotcurrenthemi,1,1e-12,-10e-12);		        /* plot hemichannel current */
  plot_param ("IHemi", colr=6, pl=6);
   
  plot_func(plotcurrentclca,1,1e-12,-10e-12);		        /* plot ClCa current */
  plot_param ("IClCa", colr=5, pl=6);
   
  plot_func(plotvoffset,1,0, -0.02);		        	/* plot voffset */
  plot_param ("Vext", 4, 4);

  plot_v_nod(ct=ha,cn=midha,n=soma,Vmin=-.07,Vmax=0,colr=green,"", -1, 1.5); /* plot Vha*/
  plot_v_nod(ct=ha,cn=midha,n=htip,Vmin=-.07,Vmax=0,colr=ltgreen,"",-1,1.5); /* plot Vha*/
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
      //vclamp(ndn(xcone,midcone,soma), vol, vtime, exptdur); // VCLAMP 4 - CONE
      //vclamp(ndn(ha,midha,soma), vol, 0.1, 0.1); //VCLAMP - 5 - HORZONTAL CELL
     for (st=0; st<dtrial; st+= stiminc){
      readephapticfeedback();
       step(stiminc);
     } 
   }
}


