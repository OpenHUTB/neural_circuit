/* Ephaptic feedback from horizontal cells to cones */

/*------------------------------------------------------*/

// double pnx_func(double atp) 

#include "phbuf_inc.cc"

/*------------------------------------------------------*/

int node1 = 2000;		// allow for hz contacts with nres = 2
int nodeatp = 5000;
int nres = 2;
int nodeclca = 0;
int nodeampa = 0;
int nodeph = 4000;
int nodephl = 5000;
int nodephh = 6000;
#define NRESMUL 500		// allow for 500 hz contacts

double ph_tau;
double ph_gain;
double ph_offset;
double phh_tau;
double phh_gain;
double phh_offset;
double phn_gain;
double phn_offset;
double phn_ntoffset;

double frac_ampa;
double frac_clca;

int gjmod;
int make_pnx;

void setcelconn (int from_celltype, int from_cellnum, int to_celltype, int to_cellnum);
void set_ephaptic_pointers(void);

/* - - - - - - - - - - - - - - - - - - - */

void ephap_init()
{
    setptr("ph_tau",    &ph_tau);
    setptr("ph_gain",   &ph_gain);
    setptr("ph_offset", &ph_offset);
    setptr("phh_tau",   &phh_tau);
    setptr("phh_gain",  &phh_gain);
    setptr("phh_offset",&phh_offset);
    setptr("phn_gain",  &phn_gain);
    setptr("phn_offset",&phn_offset);
    setptr("phn_ntoffset",&phn_ntoffset);
    setptr("gHemi",     &gHemi);
    setptr("make_pnx",  &make_pnx);
    setptr("frac_ampa", &frac_ampa);
    setptr("frac_clca", &frac_clca);
}

/* - - - - - - - - - - - - - - - - - - - */

void sethzconns(void)
{
    int pa, pb, ct1, cn1, ct2, cn2;
    synapse *sepnt = NULL; 
    elem *epnt = NULL;

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
      sepnt = (synapse *)epnt;
      ct2 = sepnt->node2a;  
      if (ct2==ha || ct2==hb) {          /* if the postsynaptic cell is ha or hb */
             ct1 = sepnt->node1a;                                    
             cn1 = sepnt->node1b;
             cn2 = sepnt->node2b;
	     setcelconn(ct2,cn2,ct1,cn1);
      }
   }
}

/* - - - - - - - - - - - - - - - - - - - */

void ephapsynapnode(void)
{
    int i, pa, pb, ct1, cn1, cd1, ct2, cn2, cd2;
    synapse *sepnt = NULL;
    elem *epnt = NULL;
    resistor *r1;
    gapjunc *gj1;
    pannex *pnx1;
    loadelem *re = NULL;
    capac *ce = NULL;
    synap *spnt = NULL;
    conlst *c = NULL;
    chan *c1,*c2,*c3;
    comp *rext_pnt = NULL;
    node *nd = NULL;
    vbuf *vb = NULL;
    nbuf *nb = NULL;

  
   // add external resistors and capacitor
 
   if (notinit(ph_gain))   ph_gain   = 0;		// gain for pH effect, 0 => off
   if (notinit(ph_tau))    ph_tau    = 200e-3;		// pH low pass filter time constant
   if (notinit(ph_offset)) ph_offset = -0.035;		// offset for pH effect 
   if (notinit(phh_gain))  phh_gain  = 0;		// gain for pH effect, 0 => off
   if (notinit(phh_tau))   phh_tau   = 5e-3;		// pH high pass filter time constant
   if (notinit(phh_offset)) phh_offset = 0;		// offset for pH HP effect 
   if (notinit(phn_offset)) phn_offset = 0;		// volt offset for volt -> pH translation 
   if (notinit(phn_ntoffset)) phn_ntoffset = 7.4;	// pH offset for volt -> pH translation 
   if (notinit(phn_gain))    phn_gain  = 100;		// gain for volt -> pH translation 
   if (notinit(gjmod))         gjmod   = 100;
   if (notinit(make_pnx))   make_pnx   = 0;		// = 1 -> make pannexins, set pH feedback

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     ct2 = sepnt->node2a;
     if (ct2==ha || ct2==hb) {  	/* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	// r1 = (resistor *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2, RESISTOR);   
	// r1->r = 1/gHemi;
	gj1  = (gapjunc *)make_conn(ct2, cn2, cd2+node1, ct2, cn2, cd2, GJ);   
	if (gjmod>0) {
	   gj1->gmax = gHemi*1.5;
	   gj1->gnv = 0.1;		// fraction non voltage-sensitive
	   gj1->vgain = 0.9;		// mV for e-fold change in rate
	   gj1->voff = 0.01;		// V (not mV) offset for e-fold change in rate
	   gj1->taun = 100.;		// relative rate for opening and closing
	} else {
	   gj1->gmax = gHemi;
	   gj1->gnv = 1;		// fraction non voltage-sensitive
	}

	if (make_pnx) {
	    pnx1 = (pannex *)make_conn(ct2, cn2, cd2, ct2, cn2, cd2+nodeatp, PNX);   
	    pnx1->gmax = gHemi* 0.002;
	    pnx1->gnv = 0.1;		// fraction non voltage-sensitive
	    pnx1->vgain = 0.9;		// mV for e-fold change in rate
	    pnx1->voff = 0.01;		// V (not mV) offset for e-fold change in rate
	    pnx1->taun = 100.;		// relative rate for opening and closing
	    pnx1->rect = 1;			// rectifying pannexin
	    pnx1->atp_decr = 1 - 0.002;	// pannexin atp decrement
            pnx1->pnxfunc = (double(*)(pannex*,double))pnx_func;	// set pannexin function for atp, etc

	    // 	pnx1->atp  = atp;

	    ce = (capac *)at(ct2, cn2, cd2+nodeatp, GNDCAP);
	    ce->c = 0.5e-12;
	    re = (loadelem *)at(ct2, cn2, cd2+nodeatp, LOAD);   
	    re->r = 1e8;
	    re->vrev = re->vrest = -0.05;
	}  /* if (make-pnx) */

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
	if (cn1==midcone)  htipcent = cd2;
	if (cn1==surrcone) htipsurr = cd2;

	// add low-pass filter from horizontal cell to external compartment for lowp pH feedback

	vb = (vbuf *)make_conn(ct2, cn2, cd2, ct2, cn2, cd2+nodephl, BUF);
	vb->offset = ph_offset;
	vb->gain = ph_gain;
	vb->tau = ph_tau;
	// vb->delay = 0.01;
	ce = (capac *)at(ct2, cn2, cd2+nodephl, GNDCAP);
	ce->c = 1e-12;

	// add resistor to ph compartment
	
	r1 = (resistor *)make_conn(ct2, cn2, cd2+nodephl, ct2, cn2, cd2+nodeph, RESISTOR);   
	r1->r = 1e6;
	ce = (capac *)at(ct2, cn2, cd2+nodeph, GNDCAP);
	ce->c = 1e-12;

	// add high-pass filter from cone  to external compartment for highp pH feedback

	vb = (vbuf *)make_conn(ct1, cn1, cd1, ct2, cn2, cd2+nodephh, BUF);
	vb->offset = phh_offset;
	vb->gain = phh_gain;
	vb->tau = phh_tau;
	vb->lphp = HP;
	// vb->delay = 0.01;
	ce = (capac *)at(ct2, cn2, cd2+nodephh, GNDCAP);
	ce->c = 1e-12;

	// add resistor to ph compartment
	
	r1 = (resistor *)make_conn(ct2, cn2, cd2+nodephh, ct2, cn2, cd2+nodeph, RESISTOR);   
	r1->r = 1e6;

	// add ntrans buf to convert voltage to pH for ext comp 
	
	nb = (nbuf *)make_conn(ct2, cn2, cd2+nodeph, ct2, cn2, cd2+node1, NBUF);
	nb->offset = phn_offset;
	nb->ntoffset = phn_ntoffset;
	nb->gain = phn_gain;
	nb->ntrans = PH;

      }  /* if (sepnt->node2a==ha) */

   }  /* foreach (epnt, SYNAPSE) */

   initsim();		//  generate compartments

   /* add extracellular compartment from Rext to channels capture channel currents */

  if (notinit(frac_ampa))   frac_ampa = 0.5;
  if (notinit(frac_clca))   frac_clca = 1;


   nodeclca = node1 + (nres - int(nres*frac_clca + 0.5)) * NRESMUL;
   nodeclca = max(nodeclca,node1);
   nodeclca = min(nodeclca,node1+(nres-1)*NRESMUL);

   nodeampa = node1 + (nres - int(nres*frac_ampa + 0.5)) * NRESMUL;
   nodeampa = max(nodeampa,node1);
   nodeampa = min(nodeampa,node1+(nres-1)*NRESMUL);

   if (ninfo >= 2) fprintf (stderr,"# nodeampa %d\n",nodeampa);
   if (ninfo >= 2) fprintf (stderr,"# nodeclca %d\n",nodeclca);
   if (ph_gain  != 0) if (ninfo >= 2) fprintf (stderr,"# nodephl  %d\n",nodephl);
   if (phh_gain != 0) if (ninfo >= 2) fprintf (stderr,"# nodephh  %d\n",nodephh);

   for(epnt=elempnt; epnt = foreach (epnt, SYNAPSE, xcone, -1, &pa, &pb); epnt=epnt->next){
     sepnt = (synapse *)epnt;
     ct2 = sepnt->node2a;
     if (ct2==ha || ct2==hb) {  			/* if the postsynaptic cell is ha or hb */
        ct1 = sepnt->node1a;
        cn1 = sepnt->node1b;
	cd1 = sepnt->node1c;
	cn2 = sepnt->node2b;
	cd2 = sepnt->node2c;
	if ((nd=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node 1 */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){		/* get conn list for presyn compartment */
	       if ((c1=findchan(c,CA,-1))!=NULL) {	/* find channel in cone terminal */
	         addchan_extern_ca(rext_pnt, c1);	/* add external comp to channel */
		 // setnt(rext_pnt, PH, 7.4);		/* set external pH */
	       }
	    }  /* if (comp1->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */

	if ((nd=ndn(ct2, cn2, cd2+nodeclca))!=NULL) {	/* get Rext node for ClCa from frac_clca */
	  rext_pnt = nd->comptr;			/* get pointer to Rext comp */
          if ((spnt=(synap *)sepnt->lptr)!=NULL) {	/* get run-time synapse struct */
            if ((c=spnt->comp1->clst)!=NULL){		/* get conn list for presyn compartment */
	       if ((c2=findchan(c,ClCa,1))!=NULL) {	/* find channel in cone terminal */
	         addchan_extern(rext_pnt, c2);		/* add external comp to channel */
	       }
	    }  /* if (comp1->clst) */
	  }  /* if (spnt=sepnt->lptr) */
	}  /* if (nd=ndn()) */

	if ((nd=ndn(ct2, cn2, cd2+nodeampa))!=NULL) { /* get Rext node 2 for synapse */
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

   set_ephaptic_pointers();				/* save center and surround cone pointers */
}

/*----------------------------------------------------------------------------------------*/

synapse *conec_sepnt = NULL;
node *conec_rext = NULL;
node *conec_nodeampa = NULL;
node *conec_nodeclca = NULL;
chan *conec_ca_chan = NULL;
chan *conec_clca_chan = NULL;
chan *conec_clcac_chan = NULL;
chan *conec_hemi_chan = NULL;

synapse *cones_sepnt = NULL;
node *cones_rext = NULL;
node *cones_nodeampa = NULL;
node *cones_nodeclca = NULL;
chan *cones_ca_chan = NULL;
chan *cones_clca_chan = NULL;
chan *cones_clcac_chan = NULL;
chan *cones_hemi_chan = NULL;

void set_ephaptic_pointers(void)

{
     int pa, pb, pc, pd, ct1, cn1, cd1, ct2, cn2, cd2;
     synapse *sepnt = NULL;
     elem *epnt = NULL;
     node *nd1=NULL, *nd2=NULL, *nd3=NULL;
     conlst *c = NULL;
     chan *c1 = NULL;
     chan *c2 = NULL;
     chan *c3 = NULL;
     chan *c4 = NULL;
     double gs, vcone, vhz;
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
	  if (cn1!=midcone && cn1!=surrcone) continue;
	  if (sepnt->lptr!=NULL) {				/* low-lev pointer to synap */
	      if ((nd1=ndn(ct2, cn2, cd2+node1))!=NULL) {	/* get Rext node */
	           if ((nd2=ndn(ct1, cn1, cd1))!=NULL) {	/* get presyn node */
	             if ((c=nd2->comptr->clst)!=NULL){
		        c1=findchan(c,CA,-1);			/* save channel pointers */
		        c2=findchan(c,ClCa,1);
		        c3=findchan(c,ClCa,2);
		      }
	           } 
	           if ((nd3=ndn(ct2, cn2, cd2))!=NULL) {	/* find postsynaptic node */
	               if ((c=nd3->comptr->clst)!=NULL){
		          c4=findchan(c,GJ,-1); 		/* save hemichannel pointer */
	               }
	           }
		   if (cn1==midcone) {
			conec_sepnt = sepnt;
			conec_rext = nd1;
			conec_nodeampa  = ndn(ct2, cn2, cd2+nodeampa);	/* save Rext node for ClCa */
			conec_nodeclca  = ndn(ct2, cn2, cd2+nodeclca);	
			conec_ca_chan = c1;
			conec_clca_chan = c2;
			conec_clcac_chan = c3;
			conec_hemi_chan = c4;
		   }
		   else if (cn1==surrcone) {
			cones_sepnt = sepnt;
			cones_rext = nd1;
			cones_nodeampa  = ndn(ct2, cn2, cd2+nodeampa);	/* save Rext node for ClCa */
			cones_nodeclca  = ndn(ct2, cn2, cd2+nodeclca);	
			cones_ca_chan = c1;
			cones_clca_chan = c2;
			cones_clcac_chan = c3;
			cones_hemi_chan = c4;
		   }
	      }
            }     /* if ((sepnt= ...) */
       }     /* if (sepnt ...) */
   }      /* for(epnt=elempnt; ;) */
}

/* - - - - - - - - - - - - - - - - - - - - - */

void onplot_ephap(void)
{
	double gs, vcone, vhz, voffset;
	node *npnt;
	synap *spnt = NULL;
	synapse *sepnt = NULL;
        chan *c1, *c2, *c3;

    if ((npnt=conec_rext)!=NULL && (spnt=(synap*)conec_sepnt->lptr)!=NULL) {
         voffset = npnt->comptr->v; 
         vcone = spnt->comp1->v;
         vhz   = spnt->comp2->v;
	 gs = record_chan(spnt->resp1, G, 0);
         if (conec_hemi_chan!=NULL) {		/* hemichannel conductance */
            currenthemi = (vhz - voffset) * conec_hemi_chan->conduct;
	 }
         if ((c1=conec_ca_chan) != NULL) {
	    currentca = (vcone - voffset - c1->gvrev) * (c1->conduct);
         }
	 if ((npnt=conec_nodeampa) != NULL) {
            currentampa = (vhz - npnt->comptr->v) * gs;
	 }
	 if ((npnt=conec_nodeclca)!=NULL &&  (c2=conec_clca_chan) != NULL) {
            currentclca = (vcone - npnt->comptr->v - vcl) * (c2->conduct);
	    c3 = conec_clcac_chan; 
            currentclcac = (vcone - npnt->comptr->v - vcl) * (c3->conduct);
	 }
         sumcurrent = currentca + currentclca + currentclcac;
         voffset_midcone = voffset;
    }

    if ((npnt=cones_rext)!=NULL && (spnt=(synap*)cones_sepnt->lptr)!=NULL) {
         voffset = npnt->comptr->v; 
         vcone = spnt->comp1->v;
         vhz   = spnt->comp2->v;
	 gs = record_chan(spnt->resp1, G, 0);
         if (cones_hemi_chan!=NULL) {		/* hemichannel conductance */
            currenthemis = (vhz - voffset) * cones_hemi_chan->conduct;
	 }
         if ((c1=cones_ca_chan) != NULL) {
	    currentcas = (vcone - voffset - c1->gvrev) * (c1->conduct);
         }
	 if ((npnt=cones_nodeampa) != NULL) {
            currentampas = (vhz - npnt->comptr->v) * gs;
	 }
	 if ((npnt=cones_nodeclca)!=NULL &&  (c2=cones_clca_chan) != NULL) {
            currentclcas = (vcone - npnt->comptr->v - vcl) * (c2->conduct);
	    c3 = cones_clcac_chan; 
            currentclcacs = (vcone - npnt->comptr->v - vcl) * (c3->conduct);
	 }
         sumcurrents = currentcas + currentclcas + currentclcacs;
         voffset_surrcone = voffset;
    }
}

/* - - - - - - - - - - - - - - - - - - - - - */

void readephapticfeedback(void)

/* read out feedback signals in cone (ct1) to horizontal cell (ct2) with ephaptic feedback */

{
  int pa, pb, pc, pd, ct1, cn1, cd1, ct2, cn2, cd2;
  synapse *sepnt = NULL;
  synap *spnt = NULL;
  elem *epnt = NULL;
  int pct, pcn;
  node *nd1=NULL, *nd2=NULL, *nd3=NULL;
  conlst *c = NULL;
  chan *c1 = NULL;
  chan *c2 = NULL;
  chan *c3 = NULL;
  double t, st;
  double gs, vcone, vhz, voffset, voffset2;
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
		      if ((c3=findchan(c,ClCa,2))!=NULL) {
		        gclcac = c3->conduct;
		      }
		     }
	           } /* if (comp1->clst) */

	           if ((nd3=ndn(ct2, cn2, cd2))!=NULL) {	/* get gHemi node */
	               if ((c=nd3->comptr->clst)!=NULL){
		          if ((c3=findchan(c,GJ,-1))!=NULL) {
			     ggj = c3->conduct;			/* hemichannel conductance */
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
		       // fprintf(stderr, "currentca %g \n" , currentca);
		       // fprintf(stderr, "gca %g v %g voff %g gvrev %g\n" , gca, vcone, voffset, c1->gvrev);
		       
	               nd3=ndn(ct2, cn2, cd2+nodeampa);	/* get Rext node for ClCa */
	               voffset2 = nd3->comptr->v; 
		       currentampa = (vhz - voffset2) * (gs);
		       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
		       
	               nd3=ndn(ct2, cn2, cd2+nodeclca);	/* get Rext node for ClCa */
	               voffset2 = nd3->comptr->v; 
		       currentclca = (vcone - voffset2 - vcl) * (gclca);
		       currentclcac   = (vcone - vcl) * (gclcac);
		       sumcurrent = currentca + currentclca + currentclcac;
		       // fprintf (stderr,"gclca %g\n",gclca);
		       // fprintf(stderr, "currentclca %g \n" , currentclca);
		       voffset_midcone = voffset;
		   } /* if (cn1==midcone) */

                   if (cn1==surrcone) {
                       // fprintf(stderr, "midcone %d  ha tip %d\n" , midcone, cd2);
                       // fprintf(stderr, "voltagedropext %g \n" , voffset);
                       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
                       //currenthemis = (vhz - voffset) * gHemi;
                       currenthemis = (vhz - voffset) * ggj;
                       // fprintf(stderr, "currenthemi %g \n" , currenthemi);
		       
	               nd3=ndn(ct2, cn2, cd2+nodeampa);	/* get Rext node for ClCa */
	               voffset2 = nd3->comptr->v; 
                       currentampas = (vhz - voffset2) * (gs);
		       // fprintf(stderr, "currentampa %g %g\n" , currentampa,gs);
		       
	               nd3=ndn(ct2, cn2, cd2+nodeclca);	/* get Rext node for ClCa */
	               voffset2 = nd3->comptr->v; 
		       currentclcas = (vcone - voffset2 - vcl) * (gclca);
		       currentclcacs  = (vcone - vcl) * (gclcac);
                       // fprintf (stderr,"gclca %g\n",gclca);
                       // fprintf(stderr, "currentclca %g \n" , currentclca);
                       voffset_surrcone = voffset;
                   } /* if (cn1==surrcone) */

	        }  /* if (nd=ndn()) */
	    }  /* if ((spnt=sepnt->lptr) */
         }   /* if sepnt */
       }    /* for (epnt;;) */
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


