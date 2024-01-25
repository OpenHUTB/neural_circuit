/* Experiment cone_hz */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "gprim.h"

double temp_freq;
double ntrials;
double stimdur;
double postdur;
double sdia;
double stimtime;
double fbtime;
double fbdur;
double minten;
double scontrast;
double Rext; 
double gHemi;
double dvsha;
double dvsc;
double coneha_cond;
double stimvolt;
double prevolt;
double startvolt;
double stopvolt;
double stepvolt;
double clca_cond;  
double clcac_cond;  
double coneca_cond;
double khz_cond;
double set_cone_v;
double set_ha_v;
double set_ha_hv;
double ca_pump;
double predur;
double set_sarea;
double frac_clca;
double scale_bf;
double Imax;
double Imin;
double cone_dark_v;

double dcrm;
int vnoise;
int reset_ca;
int ivplot;
int ivca;

int rec_ct;
int rec_cn;

/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("temp_freq", &temp_freq);
  setptr("ntrials",   &ntrials);
  setptr("stimdur",   &stimdur);
  setptr("postdur",   &postdur);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("fbtime",    &fbtime);
  setptr("fbdur",     &fbdur);
  setptr("minten",    &minten);
  setptr("scontrast", &scontrast);
  setptr("Rext",      &Rext);
  setptr("gHemi",     &gHemi);
  setptr("dvsha",     &dvsha);
  setptr("dvsc",      &dvsc);
  setptr("vnoise",    &vnoise);
  setptr("stimvolt",  &stimvolt);
  setptr("prevolt",   &prevolt);
  setptr("startvolt",   &startvolt);
  setptr("stopvolt",    &stopvolt);
  setptr("stepvolt",    &stepvolt);
  setptr("coneha_cond", &coneha_cond);
  setptr("clca_cond",   &clca_cond);
  setptr("clcac_cond",   &clcac_cond);
  setptr("coneca_cond", &coneca_cond);
  setptr("khz_cond",    &khz_cond);
  setptr("set_cone_v", 	&set_cone_v);
  setptr("set_ha_v", 	&set_ha_v);
  setptr("set_ha_hv", 	&set_ha_hv);
  setptr("ca_pump", 	&ca_pump);
  setptr("dcrm",	&dcrm);
  setptr("predur",	&predur);
  setptr("set_sarea",	&set_sarea);
  setptr("reset_ca",	&reset_ca);
  setptr("ivplot",	&ivplot);
  setptr("ivca",	&ivca);
  setptr("frac_clca",	&frac_clca);
  setptr("scale_bf",	&scale_bf);
  setptr("Imax",	&Imax);
  setptr("Imin",	&Imin);
  setptr("cone_dark_v",	&cone_dark_v);

  nvalfile = "nval_bphz.n";
  dcasarea = 1100;	// area (um2) for Ca comp shells 
  			//  relates Ca flux to [Ca]i
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

  //Set Cone and Horizontal Cell Resting Potentials
  
  vk = -0.082;
  if (notinit(dvsc))  dvsc  = -0.04;		// for dens_cone.n
  if (notinit(dvsha)) dvsha = -0.04;		// for dens_ha.n
  if (notinit(dcrm))   dcrm = 2e4;  
  if (!notinit(set_sarea)) dcasarea = set_sarea; // shell area, changes [Ca]i given ICa.
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
 
  if (notinit(cone_maxcond)) cone_maxcond  = 0;
  if (!notinit(clca_cond))   celdens[xcone][0][_CLCA] [AXOND] = clca_cond;  
  if (!notinit(clcac_cond))  celdens[xcone][0][_CLCAC][AXOND] = clcac_cond;  
  if (!notinit(coneca_cond)) celdens[xcone][0][_CA]   [AXOND] = coneca_cond;  
  if (!notinit(ca_pump))     celdens[xcone][0][_CAP]  [AXOND] = ca_pump;  
  if (!notinit(khz_cond))    celdens[ha]   [0][_KHZ]  [SOMA]  = khz_cond;  

  //setsv (hb,SCOND,1, 0);
  if(notinit(rec_ct)) rec_ct = gca;
  if (notinit(arrsiz)) arrsiz = 50;
  if (notinit(bg_inten)) bg_inten = 50e4;       /* background light intensity */
  if (notinit(ivplot)) ivplot = 0;      	/* make I/V plot */
  if (notinit(ivca)) ivca = 0;      		/* plot ICa in I/V plot */
  if (notinit(scale_bf)) scale_bf = 1;      	/* scale current before fb */
}

/*------------------------------------------------------*/
double voffset;
double voffset_midcone;
double voffset2_midcone;
double currentampa;
double currenthemi;
double currentca;
double currentclca;
double currentclcac;
double sumcurrent;
int midcone, htip;

/*------------------------------------------------------*/
int node1 = 5000;
int nres = 2;
int nodeclca  = 0;
#define NRESMUL 100

void ephapsynapnode(void)
{
    int i, pa, pb, ct1, cn1, cd1, ct2, cn2, cd2;
    synapse *sepnt = NULL;
    elem *epnt = NULL;
    resistor *r1;
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
	for (i=0; i<nres-1; i++) {
	  r1 = (resistor *)make_conn(ct2, cn2, cd2+node1+i*NRESMUL, ct2, cn2, cd2+node1+(i+1)*NRESMUL, RESISTOR);   
	  r1->r = Rext/nres;
	}
	re = (loadelem *)at(ct2, cn2, cd2+node1+(nres-1)*NRESMUL, LOAD);   
	re->r = Rext/nres;
	re->vrev = re->vrest = 0;
	for (i=0; i<nres; i++) {
	  ce = (capac *)at(ct2, cn2, cd2+node1+i*NRESMUL, GNDCAP);
	  ce->c = 0.5e-12;
	}
	if (cn1==midcone) htip = cd2;
      }  /* if (sepnt->node2a==ha) */

   }  /* foreach (epnt, SYNAPSE) */

   initsim();		//  generate compartments

   /* add extracellular compartment from Rext to capture channel currents */

   if (notinit(frac_clca)) frac_clca = 1;

   nodeclca = node1 + (nres - int(nres*frac_clca + 0.5)) * NRESMUL;
   nodeclca = max(nodeclca,node1);
   nodeclca = min(nodeclca,node1+(nres-1)*NRESMUL);

   if (ninfo >= 2) fprintf (stderr,"# nodeclca %d\n",nodeclca);
   if (ninfo >= 2) fprintf (stderr,"# nodeampa %d\n", node1+nres/2*NRESMUL);

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     if (sepnt->node2a==ha || sepnt->node2a==hb) {	/* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
        ct2 = sepnt->node2a;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	if ((nd=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node 1 */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){		/* get conn list for presyn compartment */
	       if ((c1=findchan(c,CA,-1))!=NULL) {	/* find channel in cone terminal */
	         addchan_extern_ca(rext_pnt, c1);	/* add external comp to channel */
	       }
	    }  /* if (comp1->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */
	if ((nd=ndn(ct2, cn2, cd2+nodeclca))!=NULL) {	/* get Rext node for ClCa */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){		/* get conn list for presyn compartment */
	       if ((c2=findchan(c,ClCa,1))!=NULL) {	/* find channel ClCa type 1 in cone terminal */
	         addchan_extern(rext_pnt, c2);		/* add external comp to channel */
	       }					/* Don't add external comp for ClCa type 2 chan */
	    }  /* if (comp1->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */
	if ((nd=ndn(ct2, cn2, cd2+node1+nres/2*NRESMUL))!=NULL) { /* get Rext node 2 for synapse */
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
  node *nd1 = NULL, *nd2 = NULL;
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
  double gclca=0, gclcac=0, gca=0;

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

	      if ((nd2=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node */
	         voffset = nd2->comptr->v; 
	         if ((nd1=ndn(ct1, cn1, cd1))!=NULL) {		/* get presyn node */
	           if ((c=nd1->comptr->clst)!=NULL){
		      if ((c1=findchan(c,CA,-1))!=NULL) {
		        gca = c1->conduct;
		        // fprintf (stderr,"ctype1 %d %g\n",c1->ctype,c1->conduct);
		      }
		      if ((c2=findchan(c,ClCa,1))!=NULL) {
		        gclca = c2->conduct;
		      }
		      if ((c2=findchan(c,ClCa,2))!=NULL) {
		        gclcac = c2->conduct;
		      }
		     }
	           } /* if (comp1->clst) */

		   if (cn1==midcone) {
		       // fprintf(stderr, "midcone %d  ha tip %d\n" , midcone, cd2);
		       // fprintf(stderr, "voltagedropext %g \n" , voffset);
		       currentampa = (vhz - voffset) * (gs);
		       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
		       currenthemi = (vhz - voffset) * gHemi;
		       // fprintf(stderr, "currenthemi %g \n" , currenthemi);
		       currentca   = (vcone - voffset - c1->gvrev) * (gca);
		       //fprintf(stderr, "currentca %g \n" , currentca);
		       // fprintf(stderr, "gca %g v %g voff %g gvrev %g\n" , gca, vcone, voffset, c1->gvrev);
	               nd2=ndn(ct2, cn2, cd2+nodeclca);	/* get Rext node for ClCa */
	               voffset2_midcone = nd2->comptr->v; 
		       currentclca = (vcone - voffset2_midcone - vcl) * (gclca);
		       currentclcac   = (vcone - vcl) * (gclcac);
		       sumcurrent = currentca + currentclca + currentclcac;
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

double Vhz=0;
double Ihz=0;

void onplot (void)

{
     Vhz = v(ndn(ha, 1, soma));
     Ihz = i(ndn(ha, 1, soma)) - Vhz * 4.95e-9 - 4.07e-10;
}

/*------------------------------------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, pct, pcn, plnum, cn1, pa, pb;
    int colr, pl;
    int midha, midcbp;
    double t, st, fmax,fmin, vol;
    double rmin, rmax, plsize, voffset;
    double dtrial,exptdur,postdur;
    double Vmin, Vmax;
    double Imin, Imax;
    double pscal;
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

  setonplot(onplot);
  ephapsynapnode();
    
  // if (ninfo >=1) fprintf (stderr,"# midcone %d htip %d\n", midcone, htip);  

  if (notinit(temp_freq)) temp_freq = 4;
  // fprintf(stderr, "# temp freq %d\n", temp_freq);
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  

  if (notinit(stimdur))     stimdur = 0.10;  /* stimulus duration */
  if (notinit(sdia))           sdia = 300;   /* spot diameter */
  if (notinit(stimtime))   stimtime = 0.10;  /* stimulus time */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0.5;   /* stimulus contrast */

  postdur = 0.1;
  dtrial = stimtime + stimdur + postdur;
  exptdur = 0;
  ploti = 1e-3;


  if (0) { 
    plot_v_nod(ct=xcone,cn=midcone,n=axtrm,Vmin=-.055,Vmax =-.030,colr=cyan,"Vcone_cent", -1, -1); 
    plot_v_nod(ct=xcone,cn=27,     n=axtrm,Vmin=-.055,Vmax =-.030,colr=green,"Vcone_surr", -1, -1); 
    plot_synrate_out(ct=xcone,cn=midcone,ha,1,rmin=0,rmax=500,fmax=1,colr=magenta,1,"cone_cent"); /* plot rate out */
    plot_synrate_out(ct=xcone,cn=27,     ha,1,rmin=0,rmax=500,fmax=1,colr=blue,1,"cone_surr"); /* plot rate out */

    plot_var(&currentampa,1,1e-12,-10e-12);		        /* plot AMPA chan current */
    plot_param ("IAMPA_cent", colr=7, pl=6);
   
    plot_var(&currenthemi,1,1e-12,-10e-12);		        /* plot hemichannel current */
    plot_param ("IHemi_cent", colr=6, pl=6);
   
    plot_var(&currentclca,1,1e-12,-10e-12);		        /* plot ClCa current */
    plot_param ("IClCa_cent", colr=5, pl=6);
   
    plot_var(&voffset,1,0, -0.02);		        	/* plot voffset */
    plot_param ("Vext_cent", 4, 4);

    plot_v_nod(ct=ha,cn=midha,n=soma,Vmin=-.07,Vmax=0,colr=green,"Vhz_soma", -1, 1.5); /* plot Vha*/
    //plot_v_nod(ct=ha,cn=midha,n=htip,Vmin=-.07,Vmax=0,colr=ltgreen,"",-1,1.5); /* plot Vha*/
    //plot_v_nod(ct=ha,cn=midha,n=10,Vmin=-.07,Vmax=0,colr=blue,"",    -1, 1.5); /* plot Vha*/
    plot_v_nod(ct=ha,cn=midha,n=44,Vmin=-.07,Vmax=0,colr=ltmag,"Vhz_tip",   -1, 1.5); /* plot Vha*/
    // plot_i_nod(ct=ha, cn=1, n=soma, -2e-9, 1e-9, colr=blue, "",   -1,-1.5);
    // plot_synrate_out(ct=ha,cn=midha,pct=xcone,pcn=midcone,rmin=0,rmax=5000,colr=yellow,1); 
    // plot_v_nod(ct=dbp1,cn=midcbp,n=soma,Vmin=-.045,Vmax =-.040,colr=red,"", -1, -1);
    // plot_synrate_out(ct=dbp1,cn=midcbp,rmin=0,rmax=200,colr=magenta);
    // plot_v_nod(ct=gca,cn=1,n=soma,Vmin=-.075,Vmax =-.055,colr=blue,"", -1, -1);
    //if (getn(gca,BIOPHYS)) { plot(CA,1,ndn(gca,1,soma),fmax=0.5e-6,fmin=0); 
    //				plot_param("Cai", plnum=0,plsize=0.3);}
   } 

   if (!ivplot) {
    plot_v_nod(ct=ha,cn=1,n=soma,Vmin=-.10,Vmax=-0.02,colr=green,"Vhz_soma", -1, .5); /* plot Vha*/
    plot_i_nod(ct=ha,cn=1,n=soma,Imin= -1e-9,Vmax = 1e-9,colr=yellow,"I_Hz", -1, 1); 
   }

    //plot_func(plotcurrentclca,1,10e-12,-1e-12);		        /* plot ClCa current */
    //plot_param ("IClCa_cent", colr=5, pl=6, plsize=0.5);

    //if (getn(xcone,BIOPHYS)) plot(CA,-100,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
    //				plot_param("Caocore", plnum=1,plsize=0.7);

   // if (getn(xcone,BIOPHYS)) plot(CA,-1,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
   //				plot_param("Cao coneterm", plnum=1,plsize=0.7);
   // if (getn(xcone,BIOPHYS)) plot(CA,1,ndn(xcone,midcone,1),fmax=2.0e-6,fmin=0); 
   //  				plot_param("Cai coneterm", plnum=0,plsize=0.7);

   stim_backgr(minten);


   if (notinit(prevolt))   prevolt   = -0.070; 
   if (notinit(startvolt)) startvolt = -0.12; 
   if (notinit(stopvolt))  stopvolt  = -0.02; 
   if (notinit(stepvolt))  stepvolt  =  0.005; 

   if (notinit(cone_dark_v))  cone_dark_v  =  -0.04; 

   endexp  = dtrial;
   pscal = 5e-10;

    if (ivplot) {
        graph_x(startvolt, stopvolt);
        graph_y(pscal, -pscal);
        graph_init();
     }
   				/* clamp cones at dark potential */

   for(epnt=elempnt; epnt = foreach (epnt, CONE, xcone, -1, &pa, &pb); epnt=epnt->next){
        vclamp(ndn(xcone,pb,axtrm), cone_dark_v, 0, 10);
   }

   for (stimvolt=startvolt; stimvolt<=stopvolt; stimvolt += stepvolt) {
      simtime = 0;
      vclamp(ndn(ha,1,soma), prevolt, 0, stimtime);
      step (stimtime);
      vclamp(ndn(ha,1,soma), stimvolt, stimtime, stimdur );
      step (0.005);
      if (ivplot) graph (Vhz, Ihz);
      step (stimdur-0.005);
      vclamp(ndn(ha,1,soma), prevolt, stimtime+stimdur,postdur );
      step (postdur);

     //for (st=0; st<dtrial; st+= stiminc){
     // readephapticfeedback();
     //  step(stiminc);
     //} 
   }
}


