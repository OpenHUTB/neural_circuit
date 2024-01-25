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
double ph_tau;
double ph_gain;
double ph_offset;
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
double set_cone_v_surr;
double set_ha_v;
double set_ha_hv;
double ca_pump;
double predur;
double set_sarea;
double frac_ampa;
double frac_clca;
double scale_bf;
double Imax;
double Imin;
int gjmod;

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
  setptr("ph_tau",    &ph_tau);
  setptr("ph_gain",   &ph_gain);
  setptr("ph_offset", &ph_offset);
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
  setptr("set_cone_v_surr",  &set_cone_v_surr);
  setptr("set_ha_v", 	&set_ha_v);
  setptr("set_ha_hv", 	&set_ha_hv);
  setptr("ca_pump", 	&ca_pump);
  setptr("dcrm",	&dcrm);
  setptr("predur",	&predur);
  setptr("set_sarea",	&set_sarea);
  setptr("reset_ca",	&reset_ca);
  setptr("ivplot",	&ivplot);
  setptr("ivca",	&ivca);
  setptr("frac_ampa",	&frac_ampa);
  setptr("frac_clca",	&frac_clca);
  setptr("scale_bf",	&scale_bf);
  setptr("Imax",	&Imax);
  setptr("Imin",	&Imin);
  setptr("gjmod",       &gjmod);

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
//  if (!notinit(khz_cond))    celdens[ha]   [0][_KHZ]  [SOMA]  = khz_cond;  

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
int nodeampa = 0;
int nodeph = 6000;
#define NRESMUL 100

void ephapsynapnode(void)
{
    int i, pa, pb, ct1, cn1, cd1, ct2, cn2, cd2;
    synapse *sepnt = NULL;
    elem *epnt = NULL;
    resistor *r1;
    gapjunc *gj1;
    loadelem *re = NULL;
    capac *ce = NULL;
    synap *spnt = NULL;
    conlst *c = NULL;
    chan *c1,*c2,*c3;
    comp *rext_pnt = NULL;
    node *nd = NULL;
    vbuf *vb = NULL;

   // add external resistors and capacitor
      
   if (notinit(ph_gain))   ph_gain   = 0;               // gain for pH effect, 0 => off
   if (notinit(ph_tau))    ph_tau    = 200e-3;          // pH low pass filter time constant
   if (notinit(ph_offset)) ph_offset = -0.035;          // offset for pH effect 
      
   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     if (sepnt->node2a==ha || sepnt->node2a==hb) {  /* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
        ct2 = sepnt->node2a;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	// r1 = (resistor *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2, RESISTOR);   
	// r1->r = 1/gHemi;
        gj1 = (gapjunc *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2, GJ);
        if (gjmod>0) {
           gj1->gmax = gHemi*1.5;
           gj1->gnv = 0.1;              // fraction non voltage-sensitive
           gj1->vgain = 0.9;            // mV for e-fold change in rate
           gj1->voff = 0.01;            // V (not mV) offset for e-fold change in rate
           gj1->taun = 100.;            // relative rate for opening and closing
        } else {
           gj1->gmax = gHemi;
           gj1->gnv = 1;                // fraction non voltage-sensitive
        }
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

        // add low-pass filter from horizontal cell to external compartment for pH feedback
		
	vb = (vbuf *)make_conn(ct2, cn2, cd2, ct2, cn2, cd2+nodeph, BUF);
	vb->offset = ph_offset;
	vb->gain = ph_gain;
	vb->tau = ph_tau;
	// vb->delay = 0.01;
	ce = (capac *)at(ct2, cn2, cd2+nodeph, GNDCAP);
	ce->c = 0.1e-12;

      }  /* if (sepnt->node2a==ha) */

   }  /* foreach (epnt, SYNAPSE) */

   initsim();		//  generate compartments

   /* add extracellular compartment from Rext to capture channel currents */

   if (notinit(frac_ampa)) frac_ampa = 0.5;
   if (notinit(frac_clca)) frac_clca = 1;

   nodeclca = node1 + (nres - int(nres*frac_clca + 0.5)) * NRESMUL;
   nodeclca = max(nodeclca,node1);
   nodeclca = min(nodeclca,node1+(nres-1)*NRESMUL);

   nodeampa = node1 + (nres - int(nres*frac_ampa + 0.5)) * NRESMUL;
   nodeampa = max(nodeampa,node1);
   nodeampa = min(nodeampa,node1+(nres-1)*NRESMUL);

   if (ninfo >= 2) fprintf (stderr,"# nodeclca %d\n",nodeclca);
   if (ninfo >= 2) fprintf (stderr,"# nodeampa %d\n",nodeampa);
   if (ph_gain != 0) if (ninfo >= 2) fprintf (stderr,"# nodeph   %d\n",nodeph);

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

        // connect ph node to calcium channel
	
        if ((nd=ndn(ct2, cn2, cd2+nodeph))!=NULL) {     /* get Rext node 1 */
          rext_pnt = nd->comptr;                        /* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {      /* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){           /* get conn list for presyn compartment */
               if ((c1=findchan(c,CA,-1))!=NULL) {      /* find channel in cone terminal */
                 addchan_extern_ph(rext_pnt, c1);       /* add external lowpass pH comp to Ca channel */
               }
             }  /* if (comp1->clst) */
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
  chan *c3 = NULL;
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
  double gclca=0, gclcac=0, gca=0, ggj=0;

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


                   if ((nd2=ndn(ct2, cn2, cd2))!=NULL) {        /* get gHemi node */
                       if ((c=nd2->comptr->clst)!=NULL){
                          if ((c3=findchan(c,GJ,-1))!=NULL) {
                             ggj = c3->conduct;                 /* hemichannel conductance */
                          }
                       }
                   }

		   if (cn1==midcone) {
		       // fprintf(stderr, "midcone %d  ha tip %d\n" , midcone, cd2);
		       // fprintf(stderr, "voltagedropext %g \n" , voffset);
		       // currenthemi = (vhz - voffset) * gHemi;
		       currenthemi = (vhz - voffset) * ggj;
		       // fprintf(stderr, "currenthemi %g \n" , currenthemi);
		       currentca   = (vcone - voffset - c1->gvrev) * (gca);
		       //fprintf(stderr, "currentca %g \n" , currentca);
		       // fprintf(stderr, "gca %g v %g voff %g gvrev %g\n" , gca, vcone, voffset, c1->gvrev);
	               nd2=ndn(ct2, cn2, cd2+nodeampa);	/* get Rext node for ClCa */
	               voffset2_midcone = nd2->comptr->v; 
		       currentampa = (vhz - voffset2_midcone) * (gs);
		       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
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

double ibefore;
double iduring;

double t_before;
double t_during;
int before = 0;
int during = 0;
int gflag  = 0;

void onplot(void)
{
  if (before && (simtime > t_before)) {
          if (ivca)   ibefore = currentca * scale_bf;
          else	      ibefore = i(xcone,midcone,axtrm) * scale_bf;
          before = 0;
  }
  if (during && (simtime > t_during)) {
          if (ivca) 	iduring = currentca;
          else	        iduring = i(xcone,midcone,axtrm);
          during = 0;
  }
  if (gflag && !before && !during) {
	  if (ivplot) graph (stimvolt,ibefore,iduring);
	  gflag = 0;
  }
}


/*------------------------------------------------------------------------------------*/


void clca_speed (int stype, double tau)
{
     int pa, pb;
     elem *epnt;
     chan *cp;

     for(epnt=elempnt; epnt = foreach (epnt, CHAN, xcone, -1, &pa, &pb); epnt=epnt->next) {
         //fprintf (stderr,"clca %d ctype %d\n",pb, epnt->ctype);
         if ((cp=(chan *)epnt->lptr)!= NULL) {
             if (cp->ctype==ClCa && cp->stype==stype) {
                //fprintf (stderr,"chan type %d arate %g\n",cp->ctype,cp->arate);
                cp->arate *= tau;
                cp->brate *= tau;
                // fprintf (stderr,"chan type %d arate %g brate %g\n",cp->ctype,cp->arate,cp->brate);
             }
         }
     }
}

/*------------------------------------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, i, n, pct, pcn, plnum, cn1, pa, pb;
    int colr, pl;
    int midha, midcbp;
    double t, st, fmax,fmin, vol;
    double rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax;
    //double volinc = 0.01; //VCLAMP 1
    elem *epnt;
    char corebuf[30];

  timinc = 2e-5;
  //timinc = 1e-6;
  ddcao = 2e-5;
  //euler = 1;
  setonplot(onplot);

  midcone = findmid(xcone,0,0);
  midha   = findmid(ha, 0,0);
  midcbp  = findmid(dbp1, 0,0);

  if (ninfo >=1) fprintf (stderr,"# mid cone # %d\n", midcone);
  if (ninfo >=1) fprintf (stderr,"# mid cbp  # %d\n",  midcbp);
  if (ninfo >=1) fprintf (stderr,"# mid ha  # %d\n",  midha);  

  if(notinit(Rext))  Rext  = 75e6;
  if(notinit(gHemi)) gHemi = 5e-9;
  if(notinit(gjmod)) gjmod = 0;

  ephapsynapnode(); 

  // if (ninfo >=1) fprintf (stderr,"# midcone %d htip %d\n", midcone, htip);  

  if (notinit(temp_freq)) temp_freq = 4;
  // fprintf(stderr, "# temp freq %d\n", temp_freq);
  if (notinit(ntrials)) ntrials = 1;
  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };
  

  if (notinit(stimdur))     stimdur = 1.5;  /* stimulus duration */
  if (notinit(postdur))     postdur = 1.0;  /* stimulus duration */
  if (notinit(sdia))           sdia = 300;  /* spot diameter */
  if (notinit(stimtime))   stimtime = 0.10;  /* stimulus time */
  if (notinit(fbtime))       fbtime = 1.0;   /* feedback pulse time */
  if (notinit(fbdur))         fbdur = 0.3;   /* feedback pulse time */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (notinit(scontrast)) scontrast = 0.5;   /* stimulus contrast */

  dtrial = stimtime + stimdur + postdur;
  exptdur = 0;
  ploti = 1e-3;
  t_before = stimtime+fbtime-0.01;
  t_during = stimtime+fbtime+0.01;


  if (!ivplot) {
    if (clca_cond==0) { Imax = 50e-12; Imin = -10e-12; }
    else              { Imax = 50e-12; Imin = -10e-12;}
    plot_i_nod(ct=xcone,cn=midcone,n=axtrm,Imin,Imax,colr=brown,"Icone_cent", 14, 1); 
    plot_v_nod(ct=ha,cn=midha,n=soma,Vmin=-.07,Vmax=0,colr=green,"Vhz_soma", 7, 0.5); /* plot Vha*/
    // if (ph_gain != 0)
    //  plot_v_nod(ct=ha,cn=midha,n=38+nodeph,Vmin=-.08,Vmax=0,colr=blue,"Vhz_tip_centph",   -1, 1.5); /* plot Vha*/

    //plot_var(&sumcurrent,1,50e-12,-10e-12);		        /* plot ClCa current */
    //plot_param ("Sum current", colr=yellow, pl=14, plsize=1);

    plot_var(&currentclca,1,50e-12,-10e-12);		        /* plot ClCa current */
    plot_param ("IClCa_cent", colr=5, pl=15, plsize=1);

    plot_var(&currentclcac,1,50e-12,-10e-12);		        /* plot ClCaC current */
    plot_param ("IClCaC_cent", colr=green, pl=15, plsize=1);

    // plot_func(plotcurrentca,1,1e-12,-10e-12);				/* plot Ca current */
    plot_var(&currentca,1,1e-12,-10e-12);				/* plot Ca current */
    plot_param ("ICa_cent", colr=1, pl=6,plsize=1);

    plot_var(&voffset_midcone,1, 0.01, -0.02);		        	/* plot voffset */
    plot_param ("Vext", red, 4, 0.5);

    //plot_var(&voffset2_midcone,1, 0.01, -0.02);		        	/* plot voffset */
    //plot_param ("Vext_clca", green, 4, 0.5);

    // if (getn(xcone,BIOPHYS)) plot(CA,-100,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
    // 				plot_param("Caocore", plnum=1,plsize=0.7);
    // if (getn(xcone,BIOPHYS)) plot(CA,-1,ndn(xcone,midcone,1),fmax=dcao,fmin=0); 
    // 				plot_param("Cao coneterm", plnum=1,plsize=0.7);
    if (getn(xcone,BIOPHYS)) plot(CA,1,ndn(xcone,midcone,1),fmax=3.0e-6,fmin=0); 
  				plot_param("Cai coneterm", magenta, plnum=1,plsize=0.3);
    //if (getn(xcone,BIOPHYS)) plot(CA,3,ndn(xcone,midcone,1),fmax=2.0e-6,fmin=0); 
    // 				plot_param("Cai coneterm", blue, plnum=1,plsize=0.3);
    sprintf(corebuf,"Cai core %g",dsclcac);
    if (getn(xcone,BIOPHYS)) plot(CA,dsclcac,ndn(xcone,midcone,1),fmax=3.0e-6,fmin=0); 
     				plot_param(corebuf, green,plnum=0,plsize=0.3);
  }

   stim_backgr(minten);
   // clca_speed (1,2);            // speed up channels in cleft

   if (notinit(set_cone_v))       set_cone_v  =  -0.04;
   if (notinit(set_cone_v_surr))  set_cone_v_surr  =  -0.055;

   if (notinit(prevolt))   prevolt   = -0.080; 
   if (notinit(startvolt)) startvolt = -0.050; 
   if (notinit(stopvolt))  stopvolt  = -0.000; 
   if (notinit(stepvolt))  stepvolt  =  0.005; 
   if (notinit(set_ha_hv)) set_ha_hv = -0.06; 

  if (ivplot) {
       if (notinit(Imax)) {
         if      (clca_cond <= 4e-3) Imax = 20e-12;
	 else if (clca_cond <= 8e-3) Imax = 20e-12;
         else                        Imax = 100e-12; 
       }
       graph_x(startvolt, stopvolt);
       graph_y(Imax, -10e-12);
       graph_set ("before fb", 1, 1.0);
       graph_y(Imax, -10e-12);
       graph_set ("during fb", 1, 1.0);
       graph_init();
       gpen (white);
       gmove (0.07, 0.58);
       gcwidth (0.026);
       if ((clca_cond==0 && clcac_cond==0) || (ivca >0)) {
	 gtext ("ICa");
       }
       else {
	 gtext ("Itot");
       }
       gmove (0.07, 0.58);
       gdraw (0.07, 0.58);
    }
    else {
      plot_v_nod(ct=xcone,cn=midcone,n=axtrm,Vmin=-.08,Vmax = stopvolt,colr=cyan,"Vcone_cent", 16, 0.3); 
    }

   endexp  = dtrial;


    if (notinit(setxmin)) setxmin = 0;                    // set plot to start at 0
    if (notinit(predur)) predur = 0.5;
    simtime = 0 - predur;
    vclamp(ndn(xcone,midcone,axtrm), prevolt, simtime, predur);
    step (predur);
    if (!notinit(reset_ca) && (reset_ca >0)) savemodel("conehz.m");

   for (i=0,stimvolt=startvolt; stimvolt<=stopvolt; i++, stimvolt += stepvolt) {
      simtime = 0;
      before = 1; during = 1; gflag = 1;
      // graph_pen(i+1,i+1,i+1,i+1,i+1,i+1);

     for(epnt=elempnt; epnt = foreach (epnt, CONE, xcone, -1, &pa, &pb); epnt=epnt->next){
      if (pb==midcone) {
        vclamp(ndn(xcone,pb,axtrm), prevolt, 0, stimtime );
        vclamp(ndn(xcone,pb,axtrm), stimvolt, stimtime, stimdur );
        vclamp(ndn(xcone,pb,axtrm), prevolt, stimtime+stimdur,postdur );
      }
      else if (!notinit(set_ha_v)) vclamp(ndn(xcone,pb,axtrm), set_cone_v, 0, dtrial );
      else {
	   vclamp(ndn(xcone,pb,soma), set_cone_v, 0, stimtime+stimdur-0.5 );
           vclamp(ndn(xcone,pb,soma), set_cone_v_surr, stimtime+fbtime, fbdur );
           vclamp(ndn(xcone,pb,soma), set_cone_v, stimtime+fbtime+fbdur, 0.5 );
      }
      if (!notinit(set_ha_v)) {
	   vclamp(ndn(ha,1,soma), set_ha_v, 0, stimtime+stimdur-0.5 );
           vclamp(ndn(ha,1,soma), set_ha_hv, stimtime+fbtime, fbdur );
           vclamp(ndn(ha,1,soma), set_ha_v, stimtime+fbtime+fbdur, 0.5 );
     };
    };
    for (st=0; st<dtrial; st+= stiminc){
      readephapticfeedback();
       step(stiminc);
    } 
    if (!notinit(reset_ca) && (reset_ca > 0)) restoremodel("conehz.m");
  }
}


