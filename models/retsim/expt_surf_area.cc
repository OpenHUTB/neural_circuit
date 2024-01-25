/* Experiment surf_area */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "colors.h"
#include "stimfuncs.h"

const char *celltype;
double vstep;	    /* example of how to set var */
double fracnodes;

double axon_br_dia;
double varicos_dia;
double set_drm;

int ct;

extern double chanunit[NCHANS];                    /* channel unitary conductances */
double qcond (double cond);			   /* convert from mS/cm2 to n/um2 */

/*----------------------------------------------*/

void defparams(void)

{
  setptr("celltype",    &celltype);
  setptr("vstep",       &vstep);     
  setptr("fracnodes",   &fracnodes);
  setptr("axon_br_dia", &axon_br_dia);
  setptr("varicos_dia", &varicos_dia);
  setptr("set_drm", 	&set_drm);

  nvalfile = "nval_cbp_flash.n";
}

/*----------------------------------------------*/

void setparams(void)

{
  if (!notinit(celltype)) {
     ct = find_ct(celltype);
     make_ct(ct);			/* make the cell type */
  }
  else {
     ct = gca;
     make_ct(ct);			/* make ganglion cell only */
  }
 
  set_ncel(ct,1);

  // dbp1_file = "morph_bp";
  //  if (notinit(arrsiz)) arrsiz = 50;
 

  setn(dbp1,MORPH,0);           /* set dbp1 morphology from file, def = "morph_bp" */

  if (notinit(axon_br_dia)) axon_br_dia = 0.5;  /* default dia for axon branches */
  if (notinit(varicos_dia)) varicos_dia = 1.5;  /* default dia for varicosities */
  if (notinit(den_dia)) den_dia = 0.2;          /* default dia for dendrites */
  if (!notinit(set_drm)) drm = set_drm;

  if (notinit(node_scale)) node_scale = -3.2;	/* 3: nodenum, 0.2: medium font */

}

/*----------------------------------------------*/

double region_sa[R_NREGIONS];
double region_cond[R_NREGIONS];			// 1/Rin
double region_chan[NCHANS][R_NREGIONS];
char rlabel[40];

void runexpt(void)

{
    double c,g,ri,rm,cpl,sa,dend_sa,soma_sa,vari_sa,axon_sa,totsa;
    double dend_rm,soma_rm,axon_rm,vari_rm,totrm;
    double dend_rin,soma_rin,axon_rin,vari_rin,totrin,totcap;
    double tot_sa, tot_g, tot_c;
    int n,ch,region,found;
    elem *epnt;


  step (0.0001);
  
  printf("#\n");
  printf("# Surface area and conductances of %s\n",cname[ct]);
  printf("#\n");
 
  for (n=0; n<NREGIONS; n++) {
    region_sa[n] = 0;
    region_cond[n] = 0;
    for (ch=0; ch<NCHANS; ch++) {
       region_chan[ch][n] = 0;
    }
  }

  for (epnt=elempnt; epnt=foreach(epnt,SPHERE,ct,1,-1); epnt=epnt->next) {
      double r;

    region = int(get_efield(epnt,REGION));
    r = ((sphere*)epnt)->dia/2;
    sa = 4*PI*r*r;
    if (region<R_NREGIONS) region_sa[region] += sa;

  }

  for (epnt=elempnt; epnt=foreach(epnt,CABLE,ct,1,-1); epnt=epnt->next) {
      double dia1,dia2,len, lenc;

    region = int(get_efield(epnt,REGION));
    dia1 = ((cable*)epnt)->dia;
    dia2 = ((cable*)epnt)->dia2;
    len = ((cable*)epnt)->length;
    sa = PI*(dia1+dia2)*0.5*len;
    //if (region==DENDD || region==DEND || region==DENDP) sa *= dend_dia_factor;
    //if (region==DENDD) sa *= dendd_dia_factor;
    //if (region==DENDP) sa *= dendp_dia_factor;
    if (region<R_NREGIONS) region_sa[region] += sa;
    // lenc = dist3d(epnt->nodp1,epnt->nodp2);   /* calc length */
    // fprintf(stderr,"node %d diam %6g; len %7g, calc length %7g \n", epnt->node1c,dia1, len, lenc);

  }
  
  if (!notinit(set_drm)) {			/* find Rm for regions */
     for (n=0; n<NREGIONS; n++) {
       region_chan[C_RM][n] = set_drm;
       region_chan[C_RI][n] = dri;
     }
     printf("# Rm set using set_drm\n");
  }
  else if (getn(ct,BIOPHYS)>=1) {		/* if biophysical properties set */
     for (n=0; n<NREGIONS; n++) {
       if ((rm=getdens(ct,0,C_RM,n))==0) rm = drm;  // get Rm from density file
       region_chan[_RM][n] = rm;
       if ((ri=getdens(ct,0,C_RI,n))==0) ri = dri;  // get Ri from density file
       region_chan[_RI][n] = ri;
       if ((c=getdens(ct,0,C_CM,n))==0) c = dcm;    // get Cm from density file
       region_chan[_CM][n] = c;
     }
     printf("# nval file          %s\n",nvalfile);
     printf("# density file       %s \n",densfil[ct][0]);
   } else {				/* otherwise get Rm from nval.n */
     printf("# Using Rms from nval file.\n");
     for (n=0; n<NREGIONS; n++) {
       region_chan[_RM][n] = getn(ct,NRM);
       region_chan[_RI][n] = dri;
       region_chan[_CM][n] = dcm;
     }
   }
  printf("#\n");

  for (n=0; n<NREGIONS; n++) {		// calculate region channel conductances, 
    if (region_sa[n] > 0) {
      region_cond[n] = region_sa[n] * 1e-8 / region_chan[_RM][n];
      for (ch=0; ch<C_CA7; ch++) {
        region_chan[ch][n] = region_sa[n] * 1e-8 * getdens(ct,0,ch,n) * qcond(chanunit[ch]);;
      }
    }
  }

  printf("# Region   ");
  for (n=0; n<NREGIONS; n++) 
     if (region_sa[n]>0) printf("R%-2d       ",n);
  printf("Tot \n");

  printf("# Label    ");
  for (n=0; n<NREGIONS; n++) 
     if (region_sa[n]>0) {
	  if (n==DENDD) sprintf(rlabel,"DendD"); else
	  if (n==DEND)  sprintf(rlabel,"Dend");  else
	  if (n==DENDP) sprintf(rlabel,"DendP"); else
	  if (n==SOMA)  sprintf(rlabel,"Soma");  else
	  if (n==HCK)   sprintf(rlabel,"HCK");   else
	  if (n==AXONT) sprintf(rlabel,"AxonT"); else
	  if (n==AXONP) sprintf(rlabel,"AxonP"); else
	  if (n==AXOND) sprintf(rlabel,"AxonD"); else
	  if (n==VARIC) sprintf(rlabel,"Varicosity");
	  printf("%-9s ",rlabel);
	}
  printf("Cell \n");

  printf("# %-5s    ", chname[_COL]);
  for (n=0; n<NREGIONS; n++) {
      if (region_sa[n]>0) printf("%-9s ",colornames[int(getdens(ct,0,_COL,n))]);
  }
  printf("\n");
  printf("#\n");

  tot_sa = tot_g = tot_c = 0;

  printf("# Area     ");
  for (n=0; n<NREGIONS; n++) {
      sa = region_sa[n];
      tot_sa += sa;
      if (sa>0) printf("%-9.4g ",sa);
  }
  printf("%-9.4g um2\n",tot_sa);

  printf("# CompLam  ");
  for (n=0; n<NREGIONS; n++) {
      if (region_sa[n]>0) {
        if ((cpl=getdens(ct,0,_CPL,n))<=0) cpl = getn(ct,COMPLAM);
        if (cpl<=0) cpl = complam;
        printf("%-9.4g ",cpl);
      }
  }
  printf("          Lambda/comp\n");
  printf("#\n");

  printf("# %-5s    ", chname[_RM]);
  for (n=0; n<NREGIONS; n++) {
      rm = region_chan[C_RM][n];
      g = region_cond[n];
      tot_g += g;
      if (region_sa[n]>0) printf("%-9.4g ",rm);
  }
  printf("%-9.4g Ohm-cm2\n",tot_sa/tot_g*1e-8);

  printf("# %-5s    ", chname[_RI]);
  for (n=0; n<NREGIONS; n++) {
      ri = region_chan[C_RI][n];
      if (region_sa[n]>0) printf("%-9.4g ",ri);
  }
  printf("          Ohm-cm\n");

  printf("# %-5s    ", chname[_CM]);
  for (n=0; n<NREGIONS; n++) {
      c = region_chan[C_CM][n];
      if (region_sa[n]>0) printf("%-9.4g ",c*1e6);
  }
  printf("          uF/cm2\n");
  printf("#\n");

  printf("# Rin      ");
  for (n=0; n<NREGIONS; n++) {
      g = region_cond[n];
      if (region_sa[n]>0) {
        if ( g>0) printf("%-9.3g ",1/g*1e-6);
        else      printf("inf        ");
      }
  }
  if (tot_g>0) printf("%-9.4g MOhm (from Rm only)\n",1/tot_g*1e-6);

  printf("# Cond     ");
  for (n=0; n<NREGIONS; n++) {
      g = region_cond[n];
      if (region_sa[n]>0) printf("%-9.4g ",g*1e12);
  }
  printf("%-9.4g pS   (from Rm only)\n",tot_g*1e12);

  printf("# %-5s    ", chname[_CM]);
  for (n=0; n<NREGIONS; n++) {
      c = region_sa[n]*getdens(ct,0,_CM,n)*1e-8;
      tot_c += c;
      if (region_sa[n]>0) printf("%-9.4g ",c*1e12);
  }
  printf("%-9.4g pF\n",tot_c*1e12);

  printf("# Cond/cap ");
  for (n=0; n<NREGIONS; n++) {
      g = region_cond[n];
      c = region_sa[n]*getdens(ct,0,_CM,n)*1e-8;
      if (region_sa[n]>0) printf("%-9.4g ",g/c);
  }
  printf("%-9.4g pS/pF\n",tot_g/tot_c);
  printf("#\n");

  for (ch=0; ch<C_CA7; ch++) {
    for (n=found=0; n<NREGIONS; n++) {  // check if channel nonzero
      if (region_chan[ch][n]>0) found=1;
    }
    if (found) {			// if channel nonzero, print region conductances 
      printf("# %-5s    ", chnamea[ch]);
      for (n=tot_g=0; n<NREGIONS; n++) {
        g = region_chan[ch][n];
        tot_g += g;
	if (region_sa[n]>0) printf("%-9.4g ", g/region_sa[n]*1e11);
      }
      printf("%-9.4g mS/cm2\n",tot_g/tot_sa*1e11);
    }
  }
  printf("#\n");

  for (ch=0; ch<C_CA7; ch++) {
    for (n=found=0; n<NREGIONS; n++) {  // check if channel nonzero
      if (region_chan[ch][n]>0) found=1;
    }
    if (found) {			// if channel nonzero, print region conductances 
      printf("# %-5s    ", chnamea[ch]);
      for (n=tot_g=0; n<NREGIONS; n++) {
        g = region_chan[ch][n];
        tot_g += g;
	if (region_sa[n]>0) printf("%-9.4g ",g*1e9);
      }
      printf("%-9.4g nS\n",tot_g*1e9);
    }
  }

}
