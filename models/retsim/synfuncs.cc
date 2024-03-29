/* Module synfuncs.cc */

/* Functions to connect cells, either by making synapses directly or 
   by growing dendrites and then making synapses.
   For use with retsim nc script.  
   2/16/04 J. Tukker
   8/02/05 R.G. Smith
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "ncio.h"
#include "retsim.h"
#include "retsim_var.h"

// #define DEBUG 1

struct syndistnod { 	// structure to hold node and its distance
	double dist;
	double prob;
	int node;
	elem *epnt;
    };

struct synpair { 	// structure to hold node pair and its distance
	double dist;
	double prob;
	int nodepre;
	int nodepost;
	elem *elempost;
    };

extern int *cellnums[NCELTYPES];      /* cell numbers indexed by [ct][cn] */

extern node **alt_hb_somas;

int trace_node (int ct, int cn, int sn, int tn);
int sbdsgc_conn(int prect, int precn, int postct, int postcn);    /* in sb_recfuncs.cc */
void conn_sbdsgc (int prect, int precn, int postct, int postcn);
void save_cbp_syns (int prect, int precn, int postct, int postcn, int postnode, int synout);

int noduniq(nodeint nodea, nodeint nodeb,nodeint nodec, nodeint noded, int seed);
char *makrand (int rndsiz, const char *mesg, int num);
void ncinitstate (unsigned long int seed, char *state, int n);
void drand_setstate (char *arg_state);
void restorstate(void);

/*------------------------------------------------------------*/

/* Arrays to store individual cell connections */
/* We use separate arrays to save space, because  */
/*  larger cells have more connections but are fewer. */

int *cell_conni[NCELTYPES];
int *cell_conno[NCELTYPES];
int ncells[NCELTYPES];
node **cell_somas[NCELTYPES];

void synfuncs_init(void)

{
   int ct, k, x, y;
   static int init=0;
   node *npnt;

  if (init) return;
  init = 1;

  for (ct=0; ct<nceltypes; ct++) {
    ncells[ct]	= int(getn(ct,NMADE));
    x = NCONNI+1;
    y = int(getn(ct,MAXSYNI)+1);
    cell_conni[ct] = (int *)emalloc((ncells[ct]+1)*x*y*sizeof(int));

    x = NCONNO+1;
    y = int(getn(ct,MAXSYNO)+1);
    cell_conno[ct] = (int *)emalloc((ncells[ct]+1)*x*y*sizeof(int));

    cell_somas[ct] = (node **)emalloc(ncells[ct]*sizeof(node*));
    if (ncells[ct] > 0) {
        for (k=0; k<ncells[ct]; k++) {		 	// remember somas
           npnt=findnode(ct,k+1,soma);
          *(cell_somas[ct] + k) = npnt; 
        }
    }
  } 
}

/*------------------------------------------------------------*/

void synfuncs_cleanup(void)

/* Remove cell connect arrays. Call after cell arrays have been pruned. */
{
   int ct;

  for (ct=0; ct<nceltypes; ct++) {
    efree (cell_conni[ct]);
    efree (cell_conno[ct]);
    efree (cell_somas[ct]);
  }
}

/*------------------------------------------------------------*/
    
#define cell_connin(ct,cellnum,n,cn)  (cell_conni[ct]+cellnum*x*y+n*y+cn)
#define cell_connout(ct,cellnum,n,cn) (cell_conno[ct]+cellnum*x*y+n*y+cn)

void set_cell_in(int ct, int cellnum, int connum, int cn, int val) 

{
   int x,y;

#ifdef DEBUG
  if (connum > NCONNI) 
    ncfprintf (stderr,"set_cell_int, connection number for %s too large %d\n",
						cname[ct],connum);
  if (cn > getn(ct,MAXSYNI)) 
      ncfprintf (stderr,"set_cell_in, too many connections %d\n",val);
#endif

  x = NCONNI+1;
  y = int(getn(ct,MAXSYNI)+1);
  if (ct<nceltypes && cellnum<=ncells[ct] && connum<x && cn < y) {
    *(cell_connin(ct,cellnum,connum,cn)) = val;
  }
}

/*------------------------------------------------------------*/

void set_cell_out(int ct, int cellnum, int connum, int cn, int val) 

{
   int x,y;

#ifdef DEBUG
  if (connum > NCONNO) 
    ncfprintf (stderr,"set_cell_out, connection number for %s too large %d\n",
						cname[ct],connum);
  if (cn > getn(ct,MAXSYNO)) 
    ncfprintf (stderr,"set_cell_out, too many connections %d\n",val);
#endif

  x = NCONNO+1;
  y = int(getn(ct,MAXSYNO)+1);
  if (ct<nceltypes && cellnum<=ncells[ct] && connum<x && cn < y) {
    *(cell_connout(ct,cellnum,connum,cn)) = val;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int get_cell_in(int ct, int cellnum, int connum, int cn)

{
   int x,y;
   int val;

#ifdef DEBUG
  if (connum > NCONNI) 
    ncfprintf (stderr,"get_cell_in, connection number for %s too large %d\n",
  						cname[ct],connum);
#endif
  x = NCONNI+1;
  val = 0;
  y = int(getn(ct,MAXSYNI)+1);
  if (ct<nceltypes && cellnum<=ncells[ct] && connum<x && cn < y) {
    val = *(cell_connin(ct,cellnum,connum,cn));
  }
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int get_cell_out(int ct, int cellnum, int connum, int cn) 

/* get number of output synapses from one cell type to another */

{
   int x,y;
   int val;

#ifdef DEBUG
  if (connum > NCONNO) 
    ncfprintf (stderr,"get_cell_out, connection number for %s too large %d\n",
  						cname[ct],connum);
#endif
  x = NCONNO+1;
  y = int(getn(ct,MAXSYNO)+1);
  val = 0;
  if (ct<nceltypes && cellnum<=ncells[ct] && connum<x && cn < y) {
    val = *(cell_connout(ct,cellnum,connum,cn));
  }
  return val;
}

/*------------------------------------------------------------*/

int ncel_in (int to_celltype, int to_cellnum, int from_celltype)

/* return number of presynaptic cells */

{
   int i,n,val;

  n = getconn(from_celltype,to_celltype);
  i = int(getsv(from_celltype,CONPOST,n));	/* presynaptic connection number */
  val = get_cell_in(to_celltype,to_cellnum,i,NCELLS);
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int tot_ncel_in (int to_celltype, int to_cellnum)

/* return number of presynaptic cells */

{
   int i,val,totval;

  for (totval=0,i=1; i<=NCONNI; i++) {
     totval += get_cell_in(to_celltype,to_cellnum,i,NCELLS);
  }
  return totval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int tot_ncel_ind (int to_celltype, int to_cellnum)

/* return number of presynaptic cells of different type */

{
   int i,val,totval;

  for (totval=0,i=1; i<=NCONNI; i++) {
     if (getcv(to_celltype,CELPRE,i)!=to_celltype) 
       totval += get_cell_in(to_celltype,to_cellnum,i,NCELLS);
  }
  return totval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int cel_in_num (int to_celltype, int to_cellnum, int from_celltype, int numcon)

/* return cell number of presynaptic cell */

{
   int connin,n,val;

  n = getconn(from_celltype,to_celltype);
  connin = int(getsv(from_celltype,CONPOST,n));	/* presynaptic connection number */
  val = get_cell_in(to_celltype,to_cellnum,connin,numcon);
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int ncel_out (int from_celltype, int from_cellnum, int to_celltype)

/* return number of postsynaptic cells of type "to_celltype" */

{
   int connum,val;

  connum = getconn(from_celltype,to_celltype);
  val = get_cell_out(from_celltype,from_cellnum,connum,NCELLS);
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int tot_ncel_out (int from_celltype, int from_cellnum)

/* return total number of postsynaptic cells of all types */

{
   int i,totval;

  for (totval=0,i=1; i<=NCONNO; i++) {
    totval += get_cell_out(from_celltype,from_cellnum,i,NCELLS);
  }
  //ncfprintf (stderr,"tot_ncel_out %d %d %d\n",from_celltype,from_cellnum,totval);
  return totval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int tot_ncel_outd (int from_celltype, int from_cellnum)

/* return total number of postsynaptic cells of different type */

{
   int i,totval;

  for (totval=0,i=1; i<=NCONNO; i++) {
    if (getsv(from_celltype,CELPOST,i)!=from_celltype)
       totval += get_cell_out(from_celltype,from_cellnum,i,NCELLS);
  }
  //ncfprintf (stderr,"tot_ncel_out %d %d %d\n",from_celltype,from_cellnum,totval);
  return totval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int cel_out_num (int from_celltype, int from_cellnum, int to_celltype, int numcon)

/* return cell number of postsynaptic cell */

{
   int n,val;

  n = getconn(from_celltype,to_celltype);
  val = get_cell_out(from_celltype,from_cellnum,n,numcon);
  return val;
}

/*------------------------------------------------------------*/

int connected (int from_celltype, int from_cellnum, int to_celltype, int to_cellnum)

/* determine whether one cell is connected to another */

{
    int n,i, ncells, found;

  found = 0;
  n = getconn(from_celltype,to_celltype);
  ncells = get_cell_out(from_celltype,from_cellnum,n,NCELLS);
  for (i=CELN; i<=ncells; i++) {
    if (to_cellnum==get_cell_out(from_celltype,from_cellnum,n,i)) {
      found = 1;
      break;
    }
  }
  return found;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int connected2 (int from_ct, int from_cn, int to_ct2, int to_ct3)

/* determine whether a cell converges from a second type to a third type */

{
    int n, i, cn, ncells, found;

  found = 0;
  n = getconn(from_ct,to_ct2);
  ncells = get_cell_out(from_ct,from_cn,n,NCELLS); 
  // ncfprintf(stderr,"connected2a: from %s %d to %s ncells %d\n", cname[from_ct],from_cn,cname[to_ct2],ncells);
  for (i=CELN; i<=ncells; i++) {
    cn = get_cell_out(from_ct,from_cn,n,i);   /* for each connected cell of type 2 */
    // ncfprintf(stderr,"connected2: %d %d %d\n", to_ct2,cn,ncel_out (to_ct2,cn,to_ct3));
    if (ncel_out(to_ct2,cn,to_ct3)) {       /* if it is connected to type 3 */
      found = 1;
      break;
    }
  }
  return found;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int connected3 (int from_ct, int from_cn, int to_ct2, int to_ct3, int to_ct4)

/* determine whether a cell converges through a second and third type to a fourth type */

{
    int n, i, cn, ncells, found;

  found = 0;
  n = getconn(from_ct,to_ct2);
  ncells = get_cell_out(from_ct,from_cn,n,NCELLS); 
  // ncfprintf(stderr,"connected3a: from %s %d to %s ncells %d\n", cname[from_ct],from_cn,cname[to_ct2],ncells);
  for (i=CELN; i<=ncells; i++) {
    cn = get_cell_out(from_ct,from_cn,n,i);   /* for each connected cell of type 2 */
    // ncfprintf(stderr,"connected3: %s %d %d\n", cname[to_ct2],cn,connected2(to_ct2,cn,to_ct3,to_ct4));
    if (connected2(to_ct2,cn,to_ct3,to_ct4)) {   /* if it is connected to type 4 */
      found = 1;
      break;
    }
  }
  return found;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void setcelconn (int from_celltype, int from_cellnum, int to_celltype, int to_cellnum)

/* set cell count and cell number of postsynaptic cell */

{
    int i, n, val;
						/* only register first syn */
  if (connected(from_celltype,from_cellnum,to_celltype,to_cellnum)) return;  

  n = getconn(from_celltype,to_celltype);	/* set output connection */

  val = get_cell_out(from_celltype,from_cellnum,n,NCELLS) + 1;
  if (val <= getn(from_celltype,MAXSYNO)) {
    set_cell_out(from_celltype,from_cellnum,n,NCELLS,val);
    set_cell_out(from_celltype,from_cellnum,n,val,to_cellnum);
  }
  else ncfprintf (stderr,"too many output conns %d from %s to %s\n",
			val, cname[from_celltype],cname[to_celltype]);

  i = int(getsv(from_celltype,CONPOST,n));	/* presynaptic connection number */
 					/* set input connection */

  val = get_cell_in(to_celltype,to_cellnum,i,NCELLS) + 1;
  if (val <= getn(to_celltype,MAXSYNI)) {
    set_cell_in(to_celltype,to_cellnum,i,NCELLS,val);
    set_cell_in(to_celltype,to_cellnum,i,val,from_cellnum);
  }
  else ncfprintf (stderr,"too many input conns %g from %s to %s\n",val,
			cname[from_celltype],cname[to_celltype],val);
  if (ninfo >= 4) 
    ncfprintf (stderr,"# setcelconn %s %d to %s %d\n",cname[from_celltype],from_cellnum,
						  cname[to_celltype],to_cellnum);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void rmconni (int to_celltype, int to_cellnum, int from_cellconn, int from_cellnum)

/* remove info about input connections to a cell */
{
    int i, k, n, nc, val, found;

  n = from_cellconn;
  nc = get_cell_in(to_celltype,to_cellnum,n,NCELLS);	/* get number of input cells */

  for (found=0,i=CELN; i<=nc; i++) {
    if (get_cell_in(to_celltype,to_cellnum,n,i)==from_cellnum) {
      found=1;
      set_cell_in(to_celltype,to_cellnum,n,NCELLS, nc-1);
      for (k=i; k<nc; k++) {
	 val = get_cell_in(to_celltype,to_cellnum,n,k+1); 
         set_cell_in      (to_celltype,to_cellnum,n,k,val);	  /* erase conn */ 
      }	   
      break;
    }
  }
  if (!found) ncfprintf 
	(stderr,"rmconni: connection not found to %s %d from conn %d %d\n",
			cname[to_celltype], to_cellnum, from_cellconn, from_cellnum);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void rmconno (int from_celltype, int from_cellnum, int to_cellconn, int to_cellnum)

/* remove info about output connection from a cell */
{
    int i, k, n, nc, val, found;

  n = to_cellconn;

  nc = get_cell_out(from_celltype,from_cellnum,n,NCELLS);  /* get number of output cells */
  for (found=0,i=CELN; i<=nc; i++) {
    if (get_cell_out(from_celltype,from_cellnum,n,i)==to_cellnum) {
      found=1;
      set_cell_out(from_celltype,from_cellnum,n,NCELLS,nc-1);
      for (k=i; k<nc; k++) {
         val = get_cell_out(from_celltype,from_cellnum,n,k+1);    /* erase conn */
	 set_cell_out      (from_celltype,from_cellnum,n,k,val); 
      }	   
      break;
    }
  }
  if (!found) ncfprintf 
     (stderr,"rmconno: connection not found from %s %d\n",cname[from_celltype],from_cellnum);
     (stderr,"rmconno: connection not found from %s %d %d %d\n",
			cname[from_celltype], from_cellnum, to_cellconn, to_cellnum);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void rmcelconn (int celltype,  int cellnum)

/* remove info about a cell's connections, and other cells' connections to it */

{
    int i, n, nc, cci, cti, cni, cco, cto, cno;

  for (n=1; n<=NCONNI; n++) {			 /* for each input connection type */
    nc = get_cell_in(celltype,cellnum,n,NCELLS); /* get number of presynaptic cells */
    cto = int(getcv(celltype,CELPRE,n));
    cco = int(getcv(celltype,CONPRE,n));
    for (i=CELN; i<=nc; i++) {
       cno = get_cell_in(celltype,cellnum,n,i);
       rmconno(cto,cno,cco,cellnum); 		/* remove other cell's output */
    }
    set_cell_in(celltype,cellnum,n,NCELLS,0);
  }

  for (n=1; n<=NCONNO; n++) {			  /* for each output connection type */
    nc = get_cell_out(celltype,cellnum,n,NCELLS); /* get number of postsynaptic cells */
    cti = int(getsv(celltype,CELPOST,n));	  /* get postsynaptic cell type */
    cci = int(getsv(celltype,CONPOST,n));	  /* get conn number for postsynaptic cell */
    for (i=CELN; i<=nc; i++) {			  /* for all postsynaptic cells */
       cni = get_cell_out(celltype,cellnum,n,i);  /* get the cell number */ 
       rmconni(cti,cni,cci,cellnum); 		  /* remove other cell's input */
    };
    set_cell_out(celltype,cellnum,n,NCELLS,0);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void print_connections (int ct)

/* Print connections for all cells of a given type. 
   Need to check all cells in connection arrays, even those that
   have been erased, because the cells' locations in array are not
   changed when others are erased.
*/
 
{
      int i, cn, l, n, nc, ncells, tci;
      int ctpre, conpre;
      char buf[50];
      const char *s,*sbuf;

  ncells = int(getn(ct,NMADE)) + ncell_erased[ct];
  if (ncells==1) s="";
  else           s="s";
  if (ncells>0) ncfprintf (stderr,"# cell type %s, %d cell%s:\n",cname[ct],ncells,s);
  for (cn=1; cn<=ncells; cn++) {
    if (tci=tot_ncel_in(ct,cn) > 0) {
      ncfprintf (stderr,"# cell %d: ",cn);

      for (n=1; n<=NCONNI; n++) {
        nc = get_cell_in(ct,cn,n,NCELLS);    /* get number of input cells */
        if (nc>0) {
          if (nc==1) s="";
          else       s="s";
          ctpre = int(getcv(ct,CELPRE,n));
          conpre = int(getcv(ct,CONPRE,n));
          ncfprintf(stderr,"%s from %-s, ",
		rname[int(getsv(ctpre,SRESP,conpre))],cname[int(getcv(ct,CELPRE,n))]);
          ncfprintf(stderr,"%d cell%s: ",nc,s);
          for (i=CELN; i<=nc; i++) {
            ncfprintf(stderr," %d",get_cell_in(ct,cn,n,i));
          };
          ncfprintf (stderr,"; ");
        };
      };
      ncfprintf (stderr,"\n");
    };

    if (tot_ncel_out(ct,cn) > 0) {
      sprintf (buf,"%d",cn);
      l = strlen(buf);
      if (l==1)      sbuf="#         ";
      else if (l==2) sbuf="#          ";
      else if (l==3) sbuf="#           ";
      else           sbuf="#            ";
      if (tci==0) ncfprintf (stderr,"# cell %d: ",cn);
      else        ncfprintf (stderr,sbuf);
      for (n=1; n<=NCONNO; n++) {
        nc = get_cell_out(ct,cn,n,NCELLS);    /* get number of output cells */
        if (nc>0) {
          if (nc==1) s="";
          else       s="s";		/* plural */
          ncfprintf(stderr,"%s to   %-s, ",
			rname[int(getsv(ct,SRESP,n))],cname[int(getsv(ct,CELPOST,n))]);
          ncfprintf(stderr,"%d cell%s: ", nc,s);
          for (i=CELN; i<=nc; i++) {
            ncfprintf(stderr," %d",get_cell_out(ct,cn,n,i));
          }
          ncfprintf (stderr,"; ");
        }
      }
      ncfprintf (stderr,"\n");
    }
  }  /* for (cn;;) */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void print_avg_connections (int ct)

/* for each cell type, print the average number of connections to other types */

{
      int ci, cn, n, nc, listlen, ncells, totcells, nci, nco;
      double tnci, tnco, tncisq, tncosq;
      double avgci, avgco, stdci, stdco;
      const char *s1, *s2;

  ncells = int(getn(ct,NMADE));
  totcells = ncells + ncell_erased[ct];
  if (ncells==1) { s1=""; s2=" ";}
  else           { s1="s"; s2="";}
  if (ncells>0) {
      ncfprintf (stderr,"# cell type %5s, %5d cell%s:%s ",cname[ct],ncells,s1,s2);

    stdci = stdco = 0;
    ci = 0;   				/* =1 -> has at least some inputs */
    for (listlen=0,n=1; n<=NCONNI; n++) {
      for (tnci=tncisq=0,cn=1; cn<=totcells; cn++) {
        nci = get_cell_in(ct,cn,n,NCELLS);    /* get number of input cells */
        tnci += nci;
	tncisq += nci * nci;
      }
      avgci = ((double)tnci)/ncells;
      if (ncells > 2) stdci = sqrt((tncisq - (tnci * avgci) ) / (ncells - 1));
      if (avgci>0) {
        ci = 1;
        if (n>1 && listlen>0) ncfprintf(stderr,", ");
        if (ninfo >= 3 && ncells > 2 && stdci > 0) 
		        ncfprintf (stderr,"conv from %-4s = %.3g, sd %.3g", cname[int(getcv(ct,CELPRE,n))],avgci,stdci);
	else            ncfprintf (stderr,"conv from %-4s = %.3g", cname[int(getcv(ct,CELPRE,n))],avgci);
	listlen++;
      }
    }

    for (listlen=0,n=1; n<=NCONNO; n++) {
      for (tnco=tncosq=0,cn=1; cn<=totcells; cn++) {
        nco = get_cell_out(ct,cn,n,NCELLS);    /* get number of output cells */
        tnco += nco;
	tncosq += nco * nco;
      }
      avgco = ((double)tnco)/ncells;
      if (ncells > 2) stdco = sqrt((tncosq - (tnco * avgco) ) / (ncells - 1));
      if (avgco>0) {
        if (ci>0 || listlen>0) ncfprintf (stderr,", ");
        if (ninfo>= 3 && ncells > 2 && stdco > 0) 
		ncfprintf(stderr,"div to %-4s = %.3g, sd %.3g", cname[int(getsv(ct,CELPOST,n))],avgco,stdco);
	else    ncfprintf(stderr,"div to %-4s = %.3g", cname[int(getsv(ct,CELPOST,n))],avgco);
	listlen++;
      }
    }
    ncfprintf(stderr,"\n");
  }  /* if (ncells>0) */
}

/*------------------------------------------------------------*/

#define HUGESYNDIST 1e10
#define MAXPOSTSYN 10

int sortdist(int len, syndistnod *list)
{
     int i, j, swapped;
     syndistnod tmp;
     
  for (swapped=1; swapped; ) {
    swapped = 0;
    for (i=0,j=1; j<len; i++,j++) {
       if (list[i].dist > list[j].dist) {
          tmp = list[j]; 
	  list[j] = list[i]; 
	  list[i] = tmp;
	  swapped = 1;
       }
    }
  }
  for (i=0; i<len; i++) {
     if (list[i].dist >= HUGESYNDIST) break;
  }
  return i; 		// return number of postsyn nodes within range
}

/*------------------------------------------------------------*/

void initpostsyn(syndistnod *list)
{
    int i;

 for (i=0; i<MAXPOSTSYN; i++) {
     list[i].node = -1;
     list[i].dist = HUGESYNDIST;
     list[i].epnt = NULL;
   }
}

/*------------------------------------------------------------*/

double setcondmul(double scmul, double scgrad, double segrad, double scond, int postct, int postcn, int postnode)

/* Compute synaptic conductance multiplier. 
   1. If scmul (from nval file) is set (non-zero), find _CMUL parameter for region in dens_xx.n file (default=0). 
   2. If scmul is not set, set cmul = 1.
   3. If cmul is zero (i.e. density file for postsyn cell has zero or CMUL row not included), set cmul = 1.
   4. Compute egrad = exp(somadist/segrad)
   5. If scgrad is set, compute cgrad parameter by multiplying scgrad and radial distance to soma. 
   6. Then multiply SCOND * cmul * egrad, or 1 if neither scmul and scgrad are set.
   7. Add cgrad to the product, return.
 */

{
     int region;
     double cmul, cgrad, egrad, condnew, somadist;

  if (scmul==0) cmul = 1;
  else { 
     region = nd(postct,postcn,postnode)->region;
     cmul = scmul * celdens[postct][ndens[postct][postcn]][C_CMUL][region];
     if (cmul==0) cmul = 1;
  }

  // compute exponential gradient
  
  if (segrad==0) egrad = 1;
  else {
	somadist = dist3d(ndn(postct,postcn,postnode),ndn(postct,postcn,soma));
 	egrad = exp(somadist/segrad); 
  }

  // compute linear gradient
  
  if (scgrad==0) cgrad = 0;
  else {
      cgrad = scgrad * somadist;
  }

  condnew = scond * cmul * egrad + cgrad;
  if (condnew < 0) condnew = 0;

  return condnew;
}

/*------------------------------------------------------------*/

FILE *save_syn_file = NULL;

void save_synapses ()

/* The idea here is to run the model with -d 1 and make all the synapses
   and save their info in the syn_savefile. Then for later model runs, just
   read syn_restorefile into an array, and when connect_cell() is called,
   call connect_synapse() using the saved synapse info as arguments. Since
   connect_synapse sometimes makes a new postsynaptic node, we save the 
   original node number instead of the new one. If both syn_savefile and 
   syn_restorefile are set, do not do a restore. 
*/

{
	const char *filnam;
	FILE *sfil=NULL;

  filnam = syn_savefile;			/* get save file from expt file or command line */
  if (filnam == NULL) return;			/* allow syn_restorefile to be set here */
 
  if      (strcmp(filnam,"stdout") ==0) sfil = stdout;
  else if (strcmp(filnam,"stderr") ==0) sfil = stderr;
  else if ((sfil=fopen(filnam,"w"))==NULL) {
     ncfprintf (stderr," %s: save_synapses: can't open file '%.50s'.\n",progname,filnam);
     save_syn_file = NULL;
     return;
  }
  save_syn_file = sfil;
  fprintf (sfil,"# prect precn prenode postct postcn postnode postden postepnt zdist mindist growpre growpost usedyad dyadtyp dyadc rcs rcr hcf\n");
}


/*------------------------------------------------------------*/

double *restore_syn_arr;
int syn_arr_rows;
int syn_arr_cols;
int syn_arr_indx;

void restore_synapses ()
{
     int rows, cols;
     const char *filnam;
     FILE *fil;

  filnam = syn_restorefile;			/* get restore file from expt file or command line */
  if (filnam == NULL || syn_savefile !=NULL) {  /* if syn_savefile set, don't run restore */
	  restore_syn_arr = NULL; 
	  return;
  }
  if ((fil=fopen(filnam,"r")) == NULL) {
	  ncfprintf (stderr,"# %s: restore_synapses: file %s not found\n",progname,filnam);
  } else fclose (fil);

  restore_syn_arr = fread (filnam, &syn_arr_rows, &syn_arr_cols);
  syn_arr_indx = 0;
}

/*------------------------------------------------------------*/

void save_syns (int prect, int precn, int prenode, int postct, int postcn, int postnode, int postden, 
		     elem *postepnt, double zdist, double mindist, 
		     int growpre, int growpost, int usedyad, int dyadtyp, int dyadc,
		     int rcs, int rcr, int hcf, int regp, int postsoma, FILE *filout)
{
   fprintf (filout," %-5d %-5d %-5d %-5d %-5d %-5d %-5d %8d %-8.4g %-8.4g %-3d %-3d %-3d %-3d %-3d %-3d %-3d %-3d %-3d %-5d\n",
		prect, precn, prenode, postct, postcn, postnode, postden, postepnt->elnum, 
		zdist, mindist, growpre, growpost, usedyad, dyadtyp, dyadc, rcs, rcr, hcf, regp, postsoma);
}

/*------------------------------------------------------------*/

int connect_synapse (int prect, int precn, int prenode, 
		     int postct, int postcn, int postnode, int postden, 
		     elem *postepnt,
		     double zdist, double mindist, 
		     int growpre, int growpost, int usedyad, int dyadtyp, int dyadc,
		     int rcs, int rcr, int hcf, int regp, int postsoma);
#define PRECT 0
#define PRECN 1
#define POSTCT 3
#define POSTCN 4
#define POSTELEM 7

int restore_syns (int prect, int precn, int postct, int postcn, double *restore_syn_arr)
{
   int i,j, connected;
   double *row;
   elem *postepnt;

   for (connected=0, i=syn_arr_indx; i<syn_arr_rows; i++) {
     row = restore_syn_arr + i*syn_arr_cols;
     if (row[PRECT] == prect &&
         row[PRECN] == precn &&
         row[POSTCT] == postct &&
         row[POSTCN] == postcn) {
	 if (!(postepnt=findelem(row[POSTELEM]))) {
         //      connected += connect_synapse (row[0], row[1], row[2], row[3], row[4], row[5], 
	 //	 row[6], row[7], postepnt, row[9], row[10], row[11], row[12], 
	 //	 row[13], row[14], row[15], row[16], row[17],row[18], row[19]);
            ncfprintf (stderr,"restore_syns: can't find elem %d\n",row[POSTELEM]);
	    return 0;
         }
         connected += connect_synapse (row[0], row[1], row[2], row[3], row[4], row[5], 
		 row[6], postepnt, row[8], row[9], row[10], row[11], 
		 row[12], row[13], row[14], row[15], row[16],row[17], row[18], row[19]);
     }
     else { syn_arr_indx = i; break; }
   }
   return connected;
}

#undef PRECT
#undef PRECN
#undef POSTCT
#undef POSTCN
#undef POSTELEM

/*------------------------------------------------------------*/

/* variables for "gs(param)" function, used for dyads */

/* synapse types */

// int rcs = 0;  // feedforward connection
// int rcr = 0;  // inhibitory feedback connection
// int hcf = 0;  // inhibitory feedforward connection

// #define gs(sparam) getsv(prect,sparam,rcs)	/* macro defined in retsim.h */

/*-------------------------------------------------------------*/

double get_prearbdia(int prect)
{
   if (prect >= a17 && prect <=amhs)
           return getn(prect,DTREEDIA);
   else    return getn(prect,AXARBDIA);
}

/*-------------------------------------------------------------*/

int get_arbpre(int prect)
{
   if (prect >= a17 && prect <=amhs) { 
      return  int(getn(prect,DENDARB));	/* arborization type: NBRANCHED, BRANCHED, etc. */
   } else {
      return  int(getn(prect,AXARBT));	/* arborization type: NBRANCHED, BRANCHED, etc. */
   }
}

/*-------------------------------------------------------------*/

double get_postarbdia(int prect, int postct)
{
   if ((prect >= a17  && prect <=amhs) && (postct >= rbp && postct <= hbp2)) 
	   return getn(postct,AXARBDIA);
   else    return getn(postct,DTREEDIA);
}

/*-------------------------------------------------------------*/

int get_arbpost(int prect, int postct)
{
   if ((prect >= a17 && prect <=amhs) && (postct >= rbp && postct <= hbp2)) { 
      return  int(getn(postct,AXARBT));	/* arborization type: NBRANCHED, BRANCHED, etc. */
   } else {
      return  int(getn(postct,DENDARB)); /* arborization type: NBRANCHED, BRANCHED, etc. */
   }
}


/*-------------------------------------------------------------*/

int connect_synapse (int prect, int precn, int prenode, 
		     int postct, int postcn, int postnode, int postden, 
		     elem *postepnt,
		     double zdist, double mindist, 
		     int growpre, int growpost, int usedyad, int dyadtyp, int dyadc,
		     int rcs, int rcr, int hcf, int regp, int postsoma)

{
    int k;
    int arbpre, arbpost;
    int do_connect,do_connect_rcr,do_connect_hcf;
    int midnod, nnod;
    int dfound, lsyn, nwodyad;
    int rand_cl, nfilt;
    int s_ct, s_cn;
    int dyad_syn;
    int make_new_syn;
    int rvsyn, synout;
    int isgj;
    int region, oldpostnode;
    double xn,yn,zn,FRACD;
    double xytol, ztol;
    double frac, cd1, cd2;
    double zmindist;
    double diabs,diaspc,dia_thin,cplam;
    double scond, caperm, vmax, pkm;
    node *npnt, *pre_npnt, *newpost_npnt, *newpost2_npnt;
    elem *epnt;
    synapse *spnt;
    nattrib *napnt;
    cattrib *capnt;
    chattrib *cpnt;
    gapjunc *gpnt;
    char* dstate;

#define MINZDIST 1.5

#define SYNSIZ 100
    static int dyadsyn[SYNSIZ][4];
    static int oksyn[SYNSIZ];

  if (save_syn_file != NULL && syn_restorefile == NULL) {
       save_syns (prect, precn, prenode, postct, postcn, postnode, postden, postepnt, 
            zdist, mindist, growpre, growpost, usedyad, dyadtyp, dyadc, rcs, rcr, hcf, regp, postsoma, save_syn_file);
  }

  // rcs  = getconn(prect,postct);		/* synapse type prect -> postct, defined in connect_cell */
  if (dyadc) hcf  = getconn(postct,dyadtyp);  /* synapse type postct -> dyad typ */
  else       hcf = 0;
  isgj = (getsv(prect,SRESP,rcs)==xgapj); /* = 1 -> conn is gap junction */

  arbpre  = get_arbpre(prect);		/* dend (bp) or axon (am) arborization type: NBR, BR, etc. */
  arbpost = get_arbpost(prect,postct);	/* arborization type: br, nbr, etc. */
  xytol   = getn(prect,MAXSDIST); 	/* arborization X,Y tolerance */
  ztol    = getn(postct,DENZDIST);	/* arborization Z tolerance */
  cplam = -1;
  oldpostnode = postnode;		/* save postnode in case it's too far away and a new one is made */
  if (regp==0) regp = DEND_MED;

  if (ninfo>=4) ncfprintf(stderr,"# connect_synapse arbpre %d arbpost %d\n",arbpre, arbpost);

/* Don't connect if outside of selected dendrite */

  if (postden>0)  if (!trace_node(postct,postcn,postnode,postden)) return 0;

     /* Don't connect if presynaptic cell is too far away in Z. */


  do_connect = do_connect_rcr = do_connect_hcf = 0;
  make_new_syn = 0;

/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

/* "celnode[][]" should be set to the pre-existing number of dendrites(nodes)*/

  if (growpost) {
      int newpost, newpost2;
      double tiplen;

    /* calculate position for new node */

    pre_npnt = ndn(prect,precn,prenode);
    xn = pre_npnt->xloc;
    yn = pre_npnt->yloc;
    zn = pre_npnt->zloc;

    /* Make new dendrite on postsynaptic cell */

    if (arbpost>=BRANCHED) { /* branched dendritic tree */

      newpost = ++celnode[postct][postcn];
      if ((tiplen=getn(postct,DTIPLEN)) > 0) {
	   newpost2 = newpost;
           loc(newpost2_npnt=  nd(postct,postcn,newpost2), xn,yn,zn-0.5); /* new tip end */
           newpost = ++celnode[postct][postcn];
           loc(newpost_npnt= nd(postct,postcn,newpost),xn,yn,zn-tiplen-0.5); /* near tip end */
      }
      else loc(newpost_npnt= nd(postct,postcn,newpost), xn,yn,zn-0.5); /* new tip end */

       if (postepnt->ctype != CABLE || postepnt->node2a!=postct || postepnt->node2b!=postcn) {
	  if (findelem(CABLE,postct,postcn,postnode)==NULL) {
	         ncfprintf (stderr,"connect_synapse: can't find cable connected to %d %d %d\n",
			      postct,postcn,postnode);
	      frac = 0;
	  }
          else frac=enfrac (postepnt,pre_npnt);   /* find closest place on cable */
      }
      else frac=enfrac (postepnt,pre_npnt);   /* find closest place on cable */
      midnod = 0;
      if (frac == 0.0)
        nnod = postepnt->node1c;
      else if (frac == 1.0)
        nnod = postepnt->node2c;
      else {                                /* make new node on cable */
        nnod = ++celnode[postct][postcn];
	at(postepnt, ndn(postct,postcn,postnode), frac*0.6, nd(postct,postcn,nnod));
	midnod = 1;
      }

  	/* ncfprintf (stderr,"N %d %d newpost %d nnod %d n %d frac %g %d %d\n", 
		postct, postcn, newpost, nnod, postepnt->elnum, frac,
		postepnt->node1c, postepnt->node2c); /* */


      /* if the new node is exactly at the old, then move it just a little */

      if (dist3d (newpost_npnt,ndn(postct,postcn,nnod)) < 0.1) {
        loc(newpost_npnt, newpost_npnt->xloc,newpost_npnt->yloc, newpost_npnt->zloc - 0.1);
      };


      do_connect = 1;
      diabs  = ((cable*)postepnt)->dia2;
      diaspc = getn(postct,TAPERSPC);
      // #define THIN_DEND 1.0 
      // if (diabs<THIN_DEND) {dia_thin = diabs; diaspc *= 5; } 
      // else		    dia_thin = THIN_DEND;
      #define THIN_DEND 0.5 
      if (diabs<THIN_DEND) {diaspc *= 10; } 
      dia_thin = diabs;

      // midnod=0;
      if (midnod) {   /* make tapering spine */
          // make_celseg(postct,postcn,nnod,newpost,diabs,dia_thin,regp,cplam);
          celnode[postct][postcn] = taperdistden(postct,postcn,nnod,newpost,dia_thin,
      					diaspc, 1.0, 0, regp, celnode[postct][postcn], postsoma);
      } else {
          celnode[postct][postcn] = taperdistden(postct,postcn,nnod,newpost,diabs,
      					diaspc, 1.0, 0, regp, celnode[postct][postcn], postsoma);
      }
      if (tiplen > 0) {
          make_celseg(postct,postcn,newpost2,newpost,cd1=getn(postct,DTIPDIA),cd2=getn(postct,DTIPDIA),regp,cplam);
          newpost = newpost2;
      }

      if (ninfo>=4) {
        ncfprintf(stderr,"#  conn br %s %d %d ",cname[postct],postcn,newpost);
        ncfprintf(stderr,"to %s %d 0 den %d\n", 
			cname[postct],postcn,celnode[postct][postcn]);
      }
    }	/* BRANCHED */

    else if (arbpost==NBRANCHED) {  /* unbranched dendritic tree*/
         int bnod;
         double dtol;

      if (isgj) {	

        FRACD = 0.45;  /* if gj, grow pre and postsyn branches half way */

        npnt = ndn(postct,postcn,soma);
        xn = xn*FRACD + npnt->xloc * (1-FRACD);
        yn = yn*FRACD + npnt->yloc * (1-FRACD);

        rcr = hcf = 0;			/* don't make reverse connection */
      }

	/* Check to see if an existing dendrite of postsynaptic cell is within a small dist
           of the presynaptic terminal. If so, then connect directly to it instead
           of growing a new dendrite */

      //  /* Commented out because mGluR6 receptors require separate postsynaptic compartment. */
      //  /* But OK to allow them to condense, if trconc is set to allow several to sum. */
      //
      //dtol = 1.5;			/* maximum distance to connect same terminal */
      //zmindist = 1e6;
      //bnod = pn = 0;
      //for (epnt=elempnt; epnt=foreach(epnt,CABLE,postct,postcn,-1,-1,NULL,NULL,&pn,NULL,
      //		dtol, dist2d, ndn(prect,precn,axtrm));
      //		epnt=epnt->next) {
      //
      //        zdist = endistzd(epnt,ndn(prect,precn,axtrm)); 
      //       if (zdist<zmindist) { 			/* closest? */
      //	    zmindist = zdist;
      //	    bnod = pn;				/* remenber this node */
      //        }
      //    }
	// fprintf (stderr,"cn %-3d nod %-3d dend zdist %g\n",postcn,pn,zmindist);
      zmindist = 100;		/* make new dendrite every time, required by mGluR6 ->cycG comp */
      if (zmindist < MINZDIST) {		/* if near enough, connect to the dendrite */
        newpost = bnod; 
      } 
      else {					/* otherwise, make new dendrite to soma */
        newpost = ++celnode[postct][postcn];
        loc(nd(postct,postcn,newpost), xn,yn,getn(prect,AXARBZ)-0.5);
        make_celseg (postct, postcn, newpost, soma,  cd1=getn(postct,DTIPDIA)*0.6,cd2=getn(postct,DTIPDIA),regp,cplam);
      }
      do_connect = 1;

      if (ninfo>=5) {
        printf("#  conn unbr %s %d %d\n",cname[postct],postcn,newpost);
      }

    }  /* if (arbpost==NBRANCHED) */

    //prenode  = axtrm;		/* set node for photorec to bip */
    postnode = newpost;		/* set node for synapse on bip */

  }   /* if (growpost) */

/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

  if (growpre) {
      int newpre;
      int growtype;
      double tipdia;

    newpre = ++celnode[prect][precn];

    if (prect>=a17 && prect <= amhs) {
	 growtype = DEND;
     	 tipdia = getn(prect,DTIPDIA);
    } else {
	 growtype = AXON;
	 tipdia = getn(prect,AXTIPDIA);
    }

    if (arbpre==BRANCHED) {	/* branched dendritic tree */
      if (mindist <= xytol) {

         if (postepnt->ctype != CABLE || postepnt->node2a!=postcn || postepnt->node2b!=postcn) {
	     if ((postepnt=findelem(CABLE,postct,postcn,postnode))==NULL) {
	         ncfprintf (stderr,"connect_synapse: can't find cable connected to %d %d %d\n",
			      postct,postcn,postnode);
	     frac = 0;
	     }
             else frac=enfrac (postepnt,ndn(prect,precn,axtrm));   /* find closest place on cable */
         }
         else frac=enfrac (postepnt,ndn(prect,precn,axtrm));   /* find closest place on cable */
         if (frac == 0.0)
           nnod = postepnt->node1c;
         else if (frac == 1.0)
           nnod = postepnt->node2c;
         else {                                /* make new node on cable */
           nnod = ++celnode[postct][postcn];
	   at(postepnt, ndn(postct,postcn,postnode), frac, nd(postct,postcn,nnod));
         }
	 npnt = ndn(postct,postcn,nnod);
	 loc(nd(prect,precn,newpre), npnt->xloc, npnt->yloc, npnt->zloc+0.2);
         make_celseg (prect, precn, newpre, axtrm, 
					cd1=tipdia,cd2=tipdia,growtype,cplam);

         prenode = newpre;	/* set node for synapse to gc */
         postnode = nnod;	/* set node for synapse from bip */
	 do_connect = 1;
      }
      else
        do_connect = 0;

    }	/* BRANCHED */

    else if (arbpre==NBRANCHED) { 	/* non-branched dendritic or axonal tree */
      if (isgj) {	

	npnt = ndn(prect,precn,soma);
        xn = npnt->xloc;
        yn = npnt->yloc;

	npnt = ndn(postct,postcn,soma);
        xn = xn*(1-FRACD) + npnt->xloc*FRACD;
        yn = yn*(1-FRACD) + npnt->yloc*FRACD;
      } else {
	xn = ndn(postct,postcn,postnode) -> xloc;
	yn = ndn(postct,postcn,postnode) -> yloc;
      }
      do_connect = 1;
      loc(nd(prect,precn,newpre), xn,yn,getn(prect,AXARBZ)-0.5);
      make_celseg (prect, precn, newpre, soma,  cd1=tipdia, cd2=tipdia, growtype, cplam);
      prenode = newpre;
    }
  }     /* if (growpre) */

/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

  if (!growpre && !growpost) {
int x;

    //if(prect==sbac && postct==dsgc) {
    //  if(sbdsgc_conn(prect,precn,postct,postcn))
    //    do_connect = 1;
    //}
    //else if (isgj) {	
    if (isgj) {	
        //prenode = soma;		/* set presynaptic node for gj */
        //postnode = soma;		/* set postsynaptic node for gj */
        do_connect = 1;
    }
    else { /* otherwise, make synapse if nearest elem is near enough */

      /* We know the elem is close enough, */
      /*   but if nearest node is too far away, then make new node on post elem */

// if (postct==6 && postcn==92 && precn==3179) {
// 	fprintf (stderr,"prect %d precn %d prenode %d postct %d postcn %d celnode %d\n",
// 			prect,precn,prenode,postct,postcn,celnode[postct][postcn]);
// }
      if (dist3d(ndn(prect,precn,prenode),ndn(postct,postcn,postnode)) > xytol) {

        if (postepnt->ctype != CABLE || postepnt->node2a!=postct || postepnt->node2b!=postcn) {
	     if ((postepnt=findelem(CABLE,postct,postcn,postnode))==NULL) {
	         ncfprintf (stderr,"connect_synapse: can't find cable connected to %d %d %d\n",
			      postct,postcn,postnode);
	     frac = 0;
	     }
             else frac=enfrac (postepnt,ndn(prect,precn,prenode));   /* find closest place on cable */
        }
        else frac=enfrac (postepnt,ndn(prect,precn,prenode));   /* find closest place on cable */
        midnod = 0;
        if (frac == 0.0)
          nnod = postepnt->node1c;
        else if (frac == 1.0)
          nnod = postepnt->node2c;
        else {                                /* make new node on cable */
          nnod = ++celnode[postct][postcn];
	  at(postepnt, ndn(postct,postcn,postnode), frac, nd(postct,postcn,nnod));
        }
        postnode = nnod;	/* set node for synapse from bip */
      }

      if (mindist<xytol)
        do_connect = 1;

      if (ninfo >= 4) ncfprintf (stderr,"# connect synapse, mindist %g xytol %g\n",mindist, xytol);

    }

   /*
    ncfprintf(stderr, "connect_synapse: prect axarbdia=%g mindist=%g\n", 
			getn(prect,AXARBDIA), mindist);
      ncfprintf(stderr, "postnode=%d postepnt=%d mindist=%g \n", 
			postnode, postepnt, mindist); /* */
  }

/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

  /* Now make synapse. */
  /* Find how many synapses from presyn -> postsyn dyadtype cell. */

  lsyn = 0;
  if (usedyad && dyadtyp>0) {
        int drseed;

    for (k=0; k<SYNSIZ; k++) {
      dyadsyn[k][0] = 0;	/* the synapse elem number */
      dyadsyn[k][1] = 0;	/* count of synapse's dyads */
      dyadsyn[k][2] = 0;	/* postsyn cell number */
      dyadsyn[k][3] = 0;	/* count of synapses going to same */
				/* postsynaptic cell */
    }

    /* for all synapses from presynaptic cell */
    lsyn = 0;
    for (epnt=elempnt; epnt=foreach(epnt,SYNAPSE,prect,precn,-1,-1); epnt=epnt->next) {
        int s;

      s = epnt->elnum; 
      s_ct = epnt->node2a;	/* postsyn cell type */
      s_cn = epnt->node2b;	/* postsyn cell number */

      if (epnt->node2a == dyadtyp) {

        /* If synapse connects to original (bipolar) dyad cell type, */
        /*  count number of synapses connected to each presyn cell */

	for (dfound=0,k=0; k<lsyn; k++) { /* look at all previous checked */
          if (dyadsyn[k][2] == s_cn) { /*  synapses, count postsyn */
	    dyadsyn[k][3]++;
	    dfound = 1;
	  }
        }
        if (!dfound) {	/* if cell num has not been saved yet */
          dyadsyn[lsyn][2] = s_cn; /* save postsyn cell num */
	  dyadsyn[lsyn][3]++;	 /* incr number of cells */
        }
	dyadsyn[lsyn++][0] = s;	 	 /* save synapse number */
      }
      else {	/* postsynaptic type != dyadtyp */
        if (s_ct == postct && get_efield(epnt,SDYAD)>0) {
	  /* If connects to this cell type && if it's a dyad */

	  for (dfound=0,k=0; k<lsyn; k++) {
	    if (dyadsyn[k][0]==get_efield(epnt,SDYAD)) {
	      dyadsyn[k][1]++;	/* incr # of dyads */
	      dfound = 1;
	    }
	  }
	  if (! dfound) {
	    printf ("connect_synapse: can't find synapse %d -> %d\n",
				s, (int)get_efield(epnt,SDYAD));
	  }
	}
      } 	/* type != dyadtyp */
    } 	/* foreach synapse ?s */

      drseed = noduniq(prect,precn,prenode,postcn, rseed);
      dstate = (char *)makrand (RNDSIZ,"connect_synapse",prect);
      if (dstate) ncinitstate(drseed,dstate,RNDSIZ);

  }       /* if (dyadtyp) */

  if (do_connect) {

    dyad_syn = 0;	/* element number for first element in dyad */

    if (rcs) {		/* standard feedforward synapse */

     scond = setcondmul(gs(SCMUL), gs(SCGRAD), gs(SEGRAD), gs(SCOND),  postct, postcn, postnode);
     if (scond==0) return (do_connect=0);

     // fprintf (stderr,"scond %g\n",scond);

      if (usedyad) {
        for (k=0; k<SYNSIZ; k++)
	  oksyn[k] = 0;
	for (nwodyad=0,k=0; k<lsyn; k++) {		/* find synapse dyads to Hz */
          if ( dyadsyn[k][1] < 1)			/* find synapses w/o dyads */
		oksyn[nwodyad++] = dyadsyn[k][0];	/* synapse elem # */
	} 

	//printf ("postct %d lsyn %d n %d\n",postct,lsyn,nwodyad);

	if (nwodyad>0) {	   /* if there are some synapses w/o dyads */
             if (dstate) drand_setstate(dstate);
	     rand_cl = int(drand()*nwodyad);   /* select one at random */
	     dyad_syn = oksyn[rand_cl];
	     if (dstate) restorstate();

	     /*
	     spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	     spnt->ntact  = OPEN;
	     spnt->casens = gs(SENSCA);
	     spnt->caegain = gs(CAEGAIN);
	     spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
	     spnt->dyadelem = dyad_syn;
	     spnt->ngain  = gs(SGAIN);
	     spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	     spnt->maxcond= scond;
	     spnt->thresh = gs(STHRESH);
	     spnt->vrev   = gs(SVREV);
	     nfilt        = (int)gs(SNFILT);
	     spnt->nfilt2 = (short int)nfilt;
	     spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	     spnt->tfall2 = gs(SFALL);
	     spnt->vsize  = gs(SVSIZ);
	     if (gs(SCNOISE)>0) {
	       napnt=make_chnoise(spnt);
	       napnt->unitary=22e-12;
	     }
	     if (gs(SVNOISE)>0) {
	       napnt=make_vesnoise(spnt);
	       napnt->vsize =gs(SVSIZ);
	       napnt->cov   =gs(SCOV);
	     }
	 /* */

            make_new_syn = 1;
	}			 /* Can bp have multiple dyads? */

	else 		         /* If no free dyads, */
          make_new_syn = 1;      /* make new synapse here. */

      }  /* if usedyad */

      else make_new_syn = 1;

      if (make_new_syn) {   /* Make new synapse from prect to postct cell */

	make_new_syn = 0;	/* reset */

	//if (prect==sbac && postct==dsgc) {		/* sb to dsgc synapse */
        //     conn_sbdsgc (prect,precn,postct,postcn);
	//}
        //else 

         scond = setcondmul(gs(SCMUL), gs(SCGRAD), gs(SEGRAD), gs(SCOND),  postct, postcn, postnode);
         if (scond==0) return (do_connect=0);

 	 if (gs(SRESP)==xmglur6) {	/* inverting cone-cbp synapse */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->maxcond= 0;
	  	spnt->ntact  = CLOSE;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		spnt->npow   = 1;
		spnt->vsize  = gs(SVSIZ);
	  	if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		nfilt        = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->nkd    = 1;
		spnt->cgain  = gs(SCGAIN);
		spnt->coff   = gs(SCOFF);
		nfilt 	     = (int)gs(SCNFILT);
		spnt->nfilt3 = (short int)nfilt;
		spnt->timec3 = makfiltarr(nfilt,0,spnt->timec3,gs(SCDUR));
		spnt->trconc = gs(STRCONC);		/* additional ves size factor */ 
		spnt->secmsg = 1;
		spnt->mesg2 = CGMP;			/* set second mesg type */
		ename (spnt,&synout);

	  /* make cGMP channel */

	  // epnt = at(nd(postct,postcn,postnode),CHAN);
	    cpnt = make_chan(spnt, CGMP, 1);
	  	cpnt->maxcond = scond;
	  	cpnt->vrev    = gs(SVREV);
	  	cpnt->taua    = 1; cpnt->taub    = 1;
	  	cpnt->tauc    = 1; cpnt->taud    = 1;
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  napnt->unitary=20e-12;
		}
		if ((caperm=gs(SPCA))!=0) {
		   capnt = (cattrib*)chanattr(spnt,CACOMP);
		   capnt->cashell = 3;
		   capnt->cai  = 10e-9;
		   capnt->caflg = 1;
		   capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcacgmp;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   else capnt->vmax = 2e-4;
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		   else capnt->pkm  = 10e-6;
		}
					/* set 2nd mesg conc factor from chan */ 
		spnt->mesgconc = get_chan_trconc(cpnt); 
	}
	else if (int(gs(SRESP))==xampa || int(gs(SRESP))==xampa1) {		/* deactivating ampa channel */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	  	spnt->maxcond= scond;
		spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	  	spnt->thresh = gs(STHRESH);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,AMPA);
		capnt->stype = 1;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in chanampa.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcaampa;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		}
		set_chan_offsets (capnt,postct,C_AMPA1);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		// capnt->tauc = 10; 		/* to prevent inactivation, or use xampa5 */
		ename (spnt,&synout);

	  }
	else if (gs(SRESP)==xampa2) {		/* deactivating ampa channel */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
		spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,AMPA);
		capnt->stype = 2;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in chanampa.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcaampa;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		}
		set_chan_offsets (capnt,postct,C_AMPA2);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);
	  }
	else if (gs(SRESP)==xampa5) {		/* non-deactivating ampa channel */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
		spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,AMPA);
		capnt->stype = 5;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in chanampa.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcaampa;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		}
		set_chan_offsets (capnt,postct,C_AMPA5);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  // napnt->unitary=22e-12;
		}
		ename (spnt,&synout);
	  }
	else if (gs(SRESP)==xach1) {		/* non-deactivating ampa channel */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
		spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,ACH);
		capnt->stype = 1;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in chanach.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcaampa;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		}
		set_chan_offsets (capnt,postct,C_ACH1);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  // napnt->unitary=80e-12;
		}
		ename (spnt,&synout);
	  }
	else if (gs(SRESP)==xnmda) {		/* simple nmda channel, Clements & Westbrook (1991) */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,NMDA);
		capnt->stype = 1;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in channmda.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcanmda;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump, default dcavmax */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;	/* sets Km for Ca pump, default dcapkm */
		   }
		}
		set_chan_offsets (capnt,postct,C_NMDA1);
		// fprintf (stderr,"nmda voffset %g\n",capnt->voffsm);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}

		ename (spnt,&synout);

	  }
	else if (gs(SRESP)==xnmda2) {		/* more complex nmda channel */

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,NMDA);
		capnt->stype = 2;
		if ((caperm=gs(SPCA))!=0) {
		  capnt->caflg = 1;		/* sets ca comp, but check pca in channmda.cc */
		  capnt->pump = 1;		/* sets Ca pump, can set vmax, km separately */
		   if (caperm < 0) caperm = dpcanmda;
		   capnt->caperm = caperm;	/* Ca perm for postsyn chan */
		   if ((vmax=gs(SCAVMAX))!=0) {
		      if (vmax > 0) capnt->vmax = vmax;	/* sets Vmax for Ca pump */
		   }
		   if ((pkm=gs(SCAKM))!=0) {
		      if (pkm > 0) capnt->pkm = pkm;		/* sets Km for Ca pump */
		   }
		}
		set_chan_offsets (capnt,postct,C_NMDA2);
		// fprintf (stderr,"nmda voffset %g\n",capnt->voffsm);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
	  else if (int(gs(SRESP))==xgaba || int(gs(SRESP)==xgaba1)) {

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,GABA);
		capnt->stype = 1;
		set_chan_offsets (capnt,postct,C_GABA1);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
	  else if (int(gs(SRESP))==xgaba2) {

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,GABA);
		capnt->stype = 2;
		set_chan_offsets (capnt,postct,C_GABA2);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
	  else if (int(gs(SRESP))==xgaba3) {

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,GABA);
		capnt->stype = 3;
		set_chan_offsets (capnt,postct,C_GABA3);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
	  else if (int(gs(SRESP))==xgaba4) {

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,GABA);
		capnt->stype = 4;
		set_chan_offsets (capnt,postct,C_GABA4);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
	  else if (gs(SRESP)==xgly) {

	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		spnt->trconc = gs(STRCONC);
		capnt = (cattrib*)chanattr(spnt,GLY);
		capnt->stype = 1;
		set_chan_offsets (capnt,postct,C_GLY1);
		if (gs(SCNOISE)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);

	  }
          else if (isgj) {
	    if (gs(SCOND) > 0)
	      //if (precn > postcn) {
	        gpnt = make_gj (nd(prect,precn,prenode), nd(postct,postcn,postnode), gs(SCOND));
	        // gpnt->rect = 1;			/* make the gj sharply rectifying */ 
	        // fprintf (stderr,"gj %d %d\n",precn,postcn);
	        // }
	  }
          else  { 	/* (resp == xglut) */


	  spnt = make_synapse(nd(prect,precn,prenode), nd(postct,postcn,postnode), rcs);
	  	spnt->ntact  = OPEN;
	  	spnt->maxcond= scond;
	  	spnt->casens = gs(SENSCA);
	        spnt->caegain = gs(CAEGAIN);
	        spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
		spnt->ngain  = gs(SGAIN);
	        spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = gs(STHRESH);
	  	spnt->vrev   = gs(SVREV);
		nfilt = (int)gs(SNFILTH);
		spnt->nfilt1h= (short int)nfilt;
		spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
		spnt->filt1hg= gs(SHGAIN);
		spnt->filt1ho= gs(SHOFFS);
		nfilt 	     = (int)gs(SNFILT);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
	        spnt->tfall2 = gs(SFALL);
		spnt->rrpool = gs(SRRPOOL);
		spnt->rrpoolg = gs(SRRPOOLG);
		spnt->mrrpool = gs(SMRRPOOL);
		spnt->maxsrate = gs(SMAXRATE);
		spnt->vsize  = gs(SVSIZ);
		if (gs(SVNOISE)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =gs(SVSIZ);
		  napnt->cov   =gs(SCOV);
		}
		if (gs(SCNOISE)>0) {
		  napnt=make_chsyn_noise(spnt);
		  //napnt->unitary=22e-12;
		}
		ename (spnt,&synout);
	  }

	  // if (prect==dbp1 || prect==dbp2 || prect==dbp3 || prect==dbp4) 
	  //	save_cbp_syns (prect,precn,postct,postcn,postnode,synout);

      }   /* if make_new_syn */

      if (dyad_syn) spnt->dyadelem = dyad_syn;

    } 	/* if (rcs) */


    if (rcr) {		/* make inhibitory feedback synapse at same nodes as above "rcs" */

         scond = setcondmul(getsv(postct,SCMUL,rcr), getsv(postct,SCGRAD,rcr), getsv(postct,SEGRAD,rcr),
			       getsv(postct,SCOND,rcr), postct, postcn, postnode);

         if (isgj) { }

	 else if (scond!=0) {
	   spnt = make_synapse(nd(postct,postcn,postnode), nd(prect,precn,prenode), rcr);
	  	spnt->ntact  = OPEN;
	  	spnt->casens = getsv(postct,SENSCA,rcr);
	  	spnt->caegain = getsv(postct,CAEGAIN,rcr);
	        spnt->curve  = (getsv(postct,SVGAIN,rcr)<=0 ? LINEAR: EXPON);
		spnt->ngain  = getsv(postct,SGAIN,rcr);
	        spnt->vgain  = (getsv(postct,SVGAIN,rcr)<=0 ? 1 : gs(SVGAIN));
	  	spnt->vrev   = getsv(postct,SVREV,rcr);
	  	spnt->maxcond= scond;
	  	spnt->thresh = getsv(postct,STHRESH,rcr);
		nfilt 	     = (int)getsv(postct,SNFILT,rcr);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,getsv(postct,SDUR,rcr));
	  	spnt->tfall2 = getsv(postct,SFALL,rcr);
	  	spnt->rrpool = getsv(postct,SRRPOOL,rcr);
	  	spnt->rrpoolg = getsv(postct,SRRPOOLG,rcr);
	  	spnt->mrrpool = getsv(postct,SMRRPOOL,rcr);
	  	spnt->maxsrate = getsv(postct,SMAXRATE,rcr);
		switch (int(getsv(postct,SRESP,rcr))) {
		   case XGABA:
		   case XGABA1:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 1;
			set_chan_offsets (capnt,prect,C_GABA1);
			break;
		   case XGABA2:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 2;
			set_chan_offsets (capnt,prect,C_GABA2);
			break;
		   case XGABA3:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 3;
			set_chan_offsets (capnt,prect,C_GABA3);
			break;
		   case XGABA4:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 4;
			set_chan_offsets (capnt,prect,C_GABA4);
			break;
		   case XGLY:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GLY);
			capnt->stype = 1;
			set_chan_offsets (capnt,prect,C_GLY1);
			break;
		   case XVFB:
			capnt = (cattrib*)chanattr(spnt,VFB);
			capnt->stype = 1;
	  		capnt->vfbg = getsv(postct,SVFBG,rcr); // voltage feedback gain 
			if (capnt->vfbg < 0) capnt->vfbg = dvfbg;
	  		capnt->vfbo = getsv(postct,SVFBO,rcr); // voltage feedback offset 
			break;
		}

		spnt->vsize = getsv(postct,SVSIZ,rcr);
		if (getsv(postct,SVNOISE,rcr)>0) {
		  napnt      = make_vesnoise(spnt);
		  napnt->vsize =getsv(postct,SVSIZ,rcr);
		  napnt->cov   =getsv(postct,SCOV,rcr);
		}
		if (getsv(postct,SCNOISE,rcr)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
		rvsyn = spnt->elnum;
	        do_connect_rcr = 1;  
	 }
	 else do_connect_rcr = 0;  
    } 

    if (hcf) {	/* inhibitory feedforward (from hz to bp, or amac to gc) */
        int d_n2,d_n3;	/* node #s for dyadtyp cel */

      //if (notinit(nwodyad)) printf ("pre %d post %d\n",prec,postct);

      scond = setcondmul(getsv(postct,SCMUL,hcf), getsv(postct,SCGRAD,hcf), getsv(postct,SEGRAD,hcf),
			       getsv(postct,SCOND,hcf), postct, postcn, postnode);

      if (scond!=0) {

        if (usedyad && (nwodyad>0)) { 		/* if a suitable synapse to bp exists */
	   epnt = get_elempnt(oksyn[rand_cl]);
	   d_n2 = epnt->node2b;
	   d_n3 = epnt->node2c;
        } else {
	  d_n2 = precn;				/* precn */
	  d_n3 = prenode;			/* prenode */
        }

	spnt = make_synapse(nd(postct,postcn,postnode), nd(prect,d_n2,d_n3), hcf);
	  	spnt->ntact  = OPEN;
	  	spnt->casens = getsv(postct,SENSCA,hcf);
	  	spnt->caegain = getsv(postct,CAEGAIN,hcf);
	        spnt->curve  = (getsv(postct,SVGAIN,hcf)<=0 ? LINEAR: EXPON);
		spnt->dyadelem  = rvsyn;
		spnt->ngain  = getsv(postct,SGAIN,hcf);
	        spnt->vgain  = (getsv(postct,SVGAIN,hcf)<=0 ? 1 : gs(SVGAIN));
	  	spnt->thresh = getsv(postct,STHRESH,hcf);
	  	spnt->vrev   = getsv(postct,SVREV,hcf);
	  	spnt->maxcond= scond;
		nfilt 	     = (int)getsv(postct,SNFILT,hcf);
		spnt->nfilt2 = (short int)nfilt;
		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,getsv(postct,SDUR,hcf));
	  	spnt->tfall2 = getsv(postct,SFALL,hcf);
	  	spnt->rrpool = getsv(postct,SRRPOOL,hcf);
	  	spnt->rrpoolg = getsv(postct,SRRPOOLG,hcf);
	  	spnt->mrrpool = getsv(postct,SMRRPOOL,hcf);
	  	spnt->maxsrate = getsv(postct,SMAXRATE,hcf);
		switch (int(getsv(postct,SRESP,rcr))) {
		   case XGABA:
		   case XGABA1:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 1;
			set_chan_offsets (capnt,prect,C_GABA1);
			break;
		   case XGABA2:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 2;
			set_chan_offsets (capnt,prect,C_GABA2);
			break;
		   case XGABA3:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 3;
			set_chan_offsets (capnt,prect,C_GABA3);
			break;
		   case XGABA4:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GABA);
			capnt->stype = 4;
			set_chan_offsets (capnt,prect,C_GABA4);
			break;
		   case XGLY:
			spnt->trconc = getsv(postct,STRCONC,rcr);
			capnt = (cattrib*)chanattr(spnt,GLY);
			capnt->stype = 1;
			set_chan_offsets (capnt,prect,C_GLY1);
			break;
		}
		spnt->vsize  = getsv(postct,SVSIZ,hcf);
		if (getsv(postct,SVNOISE,hcf)>0) {
		  napnt=make_vesnoise(spnt);
		  napnt->vsize =getsv(postct,SVSIZ,hcf);
		  napnt->cov   =getsv(postct,SCOV,hcf);
		}
		if (getsv(postct,SCNOISE,hcf)>0) {
		  napnt=make_chnoise(spnt);
		  //napnt->unitary=22e-12;
		}
                do_connect_hcf = 1;  
      }
//      else {		/* connect hz to bp without dyad */
//
//	for (k=0; k<SYNSIZ; k++)
//	  oksyn[k] = 0;
//	for (nwodyad=0,k=0; k<lsyn; k++) {	/* find synapse dyads to Hz */
//	   if ( dyadsyn[k][1] < 1)		/* find synapses w/o dyads */
//	     oksyn[nwodyad++] = dyadsyn[k][0];	/* synapse elem # */
//	};
//	if (nwodyad>0) {	   /* if there are appropriate bipolar cells  */
//          if (dstate) drand_setstate(dstate);
//          rand_cl = int(drand()*nwodyad);   /* select one at random */
//	   dyad_syn = oksyn[rand_cl];
//         if (dstate) restorstate();
//
//	  epnt = get_elempnt(dyad_syn);
//	  // d_n2 = epnt->node2b;		/* precn */
//	  // d_n3 = epnt->node2c;		/* prenode */
//	  d_n2 = precn;				/* precn */
//	  d_n3 = prenode;			/* prenode */
//
//       if (scond!=0) 
//	 spnt = make_synapse(nd(postct,postcn,postnode), nd(dyadtyp,d_n2,d_n3), hcf);
//	  	spnt->ntact  = OPEN;
//	  	spnt->casens = getsv(postct,SENSCA,hcf);
//	  	spnt->caegain = getsv(postct,CAEGAIN,hcf);
//	        spnt->curve  = (getsv(postct,SVGAIN,hcf)<=0 ? LINEAR: EXPON);
//		spnt->ngain  = getsv(postct,SGAIN,hcf);
//	        spnt->vgain  = (getsv(postct,SVGAIN,hcf)<=0 ? 1 : gs(SVGAIN));
//	  	spnt->vrev   = getsv(postct,SVREV,hcf);
//	  	spnt->maxcond= scond;
//	  	spnt->thresh = getsv(postct,STHRESH,hcf);
//		nfilt 	     = (int)getsv(postct,SNFILT,hcf);
//		spnt->nfilt2 = (short int)nfilt;
//		spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,getsv(postct,SDUR,hcf));
//	  	spnt->tfall2 = getsv(postct,SFALL,hcf);
//	  	spnt->rrpool = getsv(postct,SRRPOOL,hcf);
//	  	spnt->rrpoolg = getsv(postct,SRRPOOLG,hcf);
//	  	spnt->mrrpool = getsv(postct,SMRRPOOL,hcf);
//	  	spnt->maxsrate = getsv(postct,SMAXRATE,hcf);
//		spnt->vsize  = getsv(postct,SVSIZ,hcf);
//		if (getsv(postct,SVNOISE,hcf)>0) {
//		  napnt=make_vesnoise(spnt);
//		  napnt->vsize =getsv(postct,SVSIZ,hcf);
//		  napnt->cov   =getsv(postct,SCOV,hcf);
//		}
//		if (getsv(postct,SCNOISE,hcf)>0) {
//		  napnt=make_chnoise(spnt);
//		  //napnt->unitary=22e-12;
//		}
//                do_connect_hcf = 1;  
// 	  } 
//       }
      else do_connect_hcf = 0;  
    } 	/* hcf */

    /* store connections between individual cells */

    setcelconn (prect, precn, postct, postcn);
    if (rcr && do_connect_rcr) setcelconn (postct, postcn, prect, precn);
    if (isgj) setcelconn (postct, postcn, prect, precn);

  }   /* do_connect */

  if (ninfo>=4) {
    ncfprintf(stderr,"# connect_synapse end, do_connect %d\n", do_connect);
  } 

  return do_connect;

}   /* func connect_synapse() */

/*-------------------------------------------------------------*/

int connect_synapse (int prect, int precn, int prenode, 
		     int postct, int postcn, int postnode)

/* int connect_synapse (int prect, int precn, int prenode, 
		     int postct, int postcn, int postnode, int postden, elem *postepnt,
		     double zdist, double mindist, 
		     int growpre, int growpost, int usedyad, int dyadtyp, int dyadc,
		     int rcs, int rcr, int hcf, int regp)
*/

{
   int rcs;

    rcs  = getconn(prect,postct);		/* synapse type prect -> postct, defined in connect_cell */
    return connect_synapse (prect, precn, prenode, 
		     postct, postcn, postnode, 0, NULL, 
		     0, 0, 
		     0,0,0,0,0,
		     rcs, 0, 0, 0, soma);
}

/*-------------------------------------------------------------*/

int connect_synapse (node *npre, node *npost)
{
    return connect_synapse (npre->nodenm1, npre->nodenm2, npre->nodenm3,
		         npost->nodenm1, npost->nodenm2, npost->nodenm3);
}

/*-------------------------------------------------------------*/

static char *cstate = NULL;

double check_region(int ct, int cn, int region, int synreg, int pn, const char *prepost)

/* Check to see if region is included in synreg. */
/* If synreg is between 1 and 20, it is assumed to be a region */
/*  If synreg is between 21 and 1000, it is assumed to be a node label */
/*  If synreg is between 1001 and 1008 it is assumed to be a set of regions */
/*   defined in the dens_xxx.n file in one of the rows _SREG1 ... _SREG8 */
/*   with probability ranging from 0 -- 1 */

{
     double match_prob;
     int sreg;

  match_prob = 1;
  if (synreg>0 && synreg < R_NREGIONS) {
      if (synreg==region) {					/* by region */
          if (ninfo>=4) ncfprintf (stderr,"# on %s %scn %d %d region %d found match %s_reg %d\n",
					cname[ct], prepost, cn, pn, region, prepost, synreg);
      } else {
	  if (ninfo>=5) ncfprintf (stderr,"# %s %scn %d region %d differs %s_reg %d\n",
			 	        cname[ct], prepost, cn, region, prepost, synreg);
          match_prob = 0;
      }
  }

  else if (synreg>R_NREGIONS && synreg<R_CREGIONS) {
         if (pn==dendn_node(ct,synreg)) { 			/* by label in 8th column */

 	     if (ninfo>=4) ncfprintf (stderr,"# %s %scn %s %d found %snode %d label %d node %d\n",
					cname[ct], prepost, cn, pn, prepost, synreg, dendn_node(ct,synreg));
	 } else {
	     if (ninfo>=5) ncfprintf (stderr,"# %s %scn %d label %d differs %d node %d\n",
					cname[ct], prepost, cn, region, synreg, dendn_node(ct,synreg));
             match_prob = 0;
	 }
  }

  else if (synreg >= R_CREGIONS && synreg <= R_CREGE) {    /* check regions in collective region SREGx */
       sreg = synreg - R_CREGIONS + _SREG1;          	  //   find the SREG line in density file
       match_prob = celdens[ct][ndens[ct][cn]][sreg][region]; // 0 <= match_prob <= 1
       if (match_prob > 0) {
           if (ninfo>=4) ncfprintf (stderr,"# on %s %scn %d %d region %d found %s_reg %d\n",
		cname[ct], prepost, cn, pn, region, prepost, sreg);
       } else {
	   if (ninfo>=5) ncfprintf (stderr,"# %s %scn %d %d region %d not included in c%s_reg %d\n",
			    		cname[ct], prepost, cn, pn, region, prepost, synreg);
       }
    }
return (match_prob);
}

/*-------------------------------------------------------------*/

int connect_cell (int prect,int precn,int postct,int postcn, int postsoma, int usedyad,int dyadtyp,
		int dyadc, int makercr, int growpre,int growpost,int postden)

/* Connect pre and postsynaptic cell with synapse. */

/*  prect = presynaptic cell type    */
/*  precn = presynaptic cell number  */
/*  postct = postsynaptic cell type  */
/*  postcn = postsynaptic cell number */
/*  usedyad =1 -> use previously made synapse as dyad */
/*  dyadtyp = cell type to look for previous synapse for dyad */
/*  growpre  =1 -> grow presynaptic cell  */
/*  growpost =1 -> grow postsynaptic cell */
/*  dyadc = 1 -> connect postsynaptic cell to dyad cell */

{
    int c,i,j,k,n,np,st,npe;
    int arbpre,arbpost,match_region,region,region_pre,region_post,sreg;
    int pn, prend, prenode, postnode;
    int rcs, rcr, hcf, node1, node2, crseed;
    int sbnod,dsgcnod,conncount,ncables,postelemtyp;
    int connected, this_connect, no_connect, this_pair, npostsyn, postn;
    int index, npairs, nodepre, nodepost, prepos, synreg, synregp, totnp;
    double d1,d2,dist,xytol,zdist,ztol;
    double maxrc, isgj, match_prob, prob, prob_pre;
    double mindist,zmindist,syndist,somadist;
    double synanni,synanno,synanpi,synanpo,synang,synrng,synrange;
    double somangle1,somangle2,theta1,theta2;
    char *elemtype;
    elem *epnt, *postepnt;
    node *npnt, *pre_npnt, *pre_npnt_soma, *post_npnt_soma, *post_npnt;
    syndistnod *disttab=NULL;
    synpair *prepost=NULL;

#define MAXPOST   5000

  if (restore_syn_arr != NULL) {
      connected = restore_syns (prect, precn, postct, postcn, restore_syn_arr);
      return connected;
  }

syndistnod postsynnod[MAXPOSTSYN];	/* space for distances from pre to postsyn node */

//  rcs  = getconn(prect,postct);		/* synapse type prect -> postct */
  rcr  = getconn(postct,prect);		/* synapse type postct -> prect */
  if (prect==postct) rcr = 0;		/* don't make autosynapse */
  if (!makercr) rcr = 0;			/* don't make reciprocal synapse */

					/* check here for more than one synapse type */
					/* first synapse is st=0, next st=1, etc */
  for(connected=st=0; st<NSYNTYPES && (rcs=getconns(prect,postct,st)); st++) {	/* synapse type prect -> postct */
   this_connect = 0;

   if (dyadc) hcf  = getconn(postct,dyadtyp);  /* synapse type postct -> dyad typ */
   else       hcf = 0;
   isgj = (getsv(prect,SRESP,rcs)==xgapj); /* = 1 -> conn is gap junction */

   arbpre   = get_arbpre(prect);	/* dend (bp) or axon (am) arborization type: NBR, BR, etc. */
   arbpost  = get_arbpost(prect,postct); /* dend (bp) or axon (bp) from amx type: NBR, BR, etc. */
   xytol  = getn(prect,MAXSDIST); 	/* arborization X,Y tolerance */
   ztol    = getn(postct,DENZDIST);	/* arborization Z tolerance */
   if (ninfo>=4) ncfprintf(stderr,"# connect_cell arbpre %d arbpost %d\n",arbpre, arbpost);
   if (ninfo>=4) ncfprintf(stderr,"# connect_cell growpre %d growpost %d\n",growpre, growpost);
 
   if (ninfo>=4) ncfprintf(stderr,"# connect_cell %s %2d to %s %2d\n",
                           cname[prect],precn,cname[postct],postcn); /* */

   synreg = getsv(prect,SYNREG,rcs);	/* specify presyn region */
   synregp = getsv(prect,SYNREGP,rcs);  /* specify postsyn region */

//
//  First check to see if postct, postcn has any cables, otherwise look for spheres
//
    for (epnt=getepnt(postct,postcn),ncables=0; epnt=foreachn(epnt,CABLE,postct,postcn,&post_npnt);
		epnt = epnt->hnnext, i++) {
	      if (arbpost > 0) ncables++;
    }
    // ncfprintf(stderr,"connect_cell: prect=%d precn=%d postct=%d postcn=%d ncables %d\n",
    // 		      prect,precn,postct,postcn,ncables);
     
    if (ncables>0) postelemtyp = CABLE;
    else           postelemtyp = SPHERE;

    if (ninfo>=4) {
       if (postelemtyp==SPHERE) 
	       ncfprintf(stderr,"# connect_cell: no dendrites in postct=%s postcn=%d, switching to spheres\n",
		       cname[postct], postcn);
    }
/* - - - - - - - - - - - - - - - - - - - - - - - */

    /* set up random state for check_region() */

   crseed = noduniq(prect,precn,postct,postcn, rseed);
   cstate = (char *)makrand (RNDSIZ,"connect_cell",prect);
   if (cstate) ncinitstate(crseed,cstate,RNDSIZ);

/* - - - - - - - - - - - - - - - - - - - - - - - */

  /* first sort out the problem by eliminating branches that are too far away */

   if (arbpre==NBRANCHED && arbpost>=BRANCHED) {

    /* nonbranched presynaptic cells synapsing onto branched postsyn dendritic tree */

    mindist = zmindist = 1e10; /* now find closest dendrite on cell we're connecting to */
    c = -1;
    pn = prenode = postnode = 0;
    zdist = 0; postepnt = NULL;


    if ((pre_npnt=nde(prect,precn,axtrm))!=NULL) {   	/* if presynaptic cell has axon, use it */

	/* find closest axonal segment (by zdist) to postcn */

        for (prend=soma,npnt=nodepnt; npnt=foreach(npnt,prect,precn,-1,-1,NULL,NULL,&prenode,NULL);
		npnt = npnt->next) {
		zdist = distzd(npnt,nde(postct,postcn,soma)); 
		if (zdist < zmindist) {
			zmindist = zdist;
			prend = prenode;
		}
	        // fprintf (stderr,"zdist %g prenode %d\n",zdist,prenode);
        }
	pre_npnt = nde (prect,precn,prend);	
	// fprintf (stderr,"prend %d\n",prend);
    } else { 
	 prend = soma;
	 pre_npnt = nde (prect,precn,soma);	
    }

    zmindist = 1e10; /* now find closest dendrite on cell we're connecting to */

    /* ncfprintf(stderr,"foreach loop 1: prect=%d precn=%d prend %d\n",prect,precn,prend); /* */

//    for (epnt=elempnt; epnt=foreach(epnt,CABLE,postct,postcn,-1,-1,NULL,NULL,&pn,NULL);
//		epnt = epnt->next, i++) {

//    use foreach loop that uses hash table to access elems
//

      for (postn=0, epnt=setforeachn2(postelemtyp,postct,postcn,&post_npnt); 
		      epnt=foreachn2(epnt); epnt = epnt->hnnext, i++) {

        // if (ninfo>=4) ncfprintf(stderr,"# connect_cell pre %2d to post %2d\n",
        //                    epnt->nodp1->nodenm3,pre_npnt->nodenm3); /* */

	if (!growpost && endist2d(epnt,pre_npnt) > xytol+5) continue;

        pn = post_npnt->nodenm3;
	c = epnt->elnum;	
        // if (pn==soma && synregp >= 0) continue;
        elemtype = get_elabel(epnt,ELABL);
	region = epnt->region;				// postsynaptic region

	// if (growpost && region!=DEND_PROX && region!=DEND_MED) continue;
			//!(streq(elemtype,"dend_prox")||streq(elemtype,"dend_med"))) continue; 
	
	/* check postsyn region */
       if ((match_prob=check_region(postct, postcn, region, synregp, pn, "post"))==0) continue;

       // if (cstate) drand_setstate(cstate);
       // match_region = (drand() <= match_prob);
       // if (cstate) restorstate();
       // if (match_region==0) break;

       /*   if ((postct >= gca) && synregp >= 0 &&		// don't place synapses on axons of GCs
		((region==AXON) || 
		 (region==AXON_DIST) ||
                 (region==AXON_THIN) ||
		 (region==HILLOCK))) continue;
	*/

	if (growpost) {
	   dist  = endist3d(epnt,pre_npnt);		// find dist to nearest elem
           zdist = endistzd(epnt,pre_npnt); 
	} else {
           dist  = dist2d(post_npnt,pre_npnt);		// find dist to nearest node
           zdist = distzd(post_npnt,pre_npnt); 
	}

         if (ninfo>=4) ncfprintf(stderr,"# foreach loop 2: %s cable=%d, precn=%d, prend=%d dist=%g zdist %g\n", 
	 				elemtype,c,precn,prend,dist,zdist); /* */

        if (dist<=mindist && abs(zdist)<ztol) { /* closest? */

	/* ncfprintf (stderr,
		"foreach loop 3: distance from [%s][%d][%d] to region %d cable %d [%d]-[%d]= %g %g %g\n",
		cname[prect],precn,prend,region,c,epnt->node1c,epnt->node2c,dist,
		endist3d(epnt,pre_npnt),dist2d(post_npnt,pre_npnt)); /* */

	  mindist = dist;                   /* remember this dist */
	  zmindist = zdist;
	  postepnt = epnt;                       /* remember this cable */
	  postnode = pn;                         /* remenber this node */
	  prenode = prend;
	  prob = match_prob;
	  region_post = region;
          //ncfprintf(stderr, "pn=%d=postnode=%d; cable=%d\n",pn, postnode, postepnt->elnum);
        }
    }     /* foreach cable on postsyn cell */

    zdist = zmindist;
    if (ninfo>=4) ncfprintf(stderr,"# connect_cell 1 mindist %g zmindist %g\n",mindist, zmindist);
    if (postepnt==NULL) break;

    // If synanpi > 0, then eliminate this synapse if closer to soma in postsynaptic cell

    no_connect = 0;

    synanpi=getsv(prect,SYNANPI,rcs);
    post_npnt = ndn(postct,postcn,postnode); 		 // saved post node
    post_npnt_soma = ndn(postct,postcn,postsoma);

    if (synanpi>0) {					 // synanpi sets exclusion zone
        somadist = dist2d(post_npnt_soma,post_npnt);
        // fprintf (stderr,"precn %d  postcn %d synanpi %g somadist %g\n",
	//	   precn,postcn,synanpi,somadist); 
        if (somadist < synanpi) {			// If the node is too close,
	   no_connect = 1;
        } 
    }

    // If synanpo > 0, then eliminate all synapses farther from soma in postsynaptic cell
    
    synanpo=getsv(prect,SYNANPO,rcs);
    if (synanpo>0) { 					// synanpo sets exclusion zone
        somadist = dist2d(post_npnt_soma,post_npnt);
        // fprintf (stderr,"precn %d  postcn %d synanpo %g somadist %g\n",
	//	   precn,postcn,synanno,somadist); 
        if (somadist > synanpo) {			// If the node is too far away,
	   no_connect = 1;
        } 
    }

    if (cstate) drand_setstate(cstate);
    match_region = (drand() <= prob);
    if (cstate) restorstate();

    if (!no_connect && match_region)
    this_connect = connect_synapse (prect, precn, prenode, 
		     postct, postcn, postnode, postden, 
		     postepnt, zdist, mindist, 
		     growpre, growpost, usedyad, dyadtyp, dyadc, rcs, rcr, hcf, region_post, postsoma);
  
  }  /*  if (arbpre==NBRANCHED && arbpost>=BRANCHED) */

 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  else if (arbpre>=BRANCHED && arbpost>=BRANCHED) {

    // branched presynaptic cell synapsing onto branched postsyn dendritic tree 
    
				// allocate space for table of distances 
    disttab = (syndistnod*) emalloc (MAXPOST*MAXPOSTSYN*sizeof(syndistnod));  
    prepost = (synpair *)   emalloc(MAXPOST*sizeof(synpair));

    mindist = zmindist = HUGESYNDIST; /* now find closest dendrite on cell we're connecting to */
    c = -1;
    pn = postnode = 0;
    zdist = 0; 
    npostsyn = 0;
    for (pre_npnt=getnpnt(prect,precn); pre_npnt=foreachn(pre_npnt,prect,precn,&prenode); 
		    	   pre_npnt = pre_npnt->hcnext) {

       /* ncfprintf(stderr,"foreach loop 1: prect=%d precn=%d prenode=%d\n",prect,precn,prenode); /* */

      region_pre = pre_npnt->region;		// presynaptic region
      if (prenode==soma && synreg >= 0 && synreg != SOMA) continue;
      if ((match_prob=check_region(prect, precn, region_pre, synreg, prenode, "pre"))==0) continue;

      initpostsyn (postsynnod);			// initialize the list of postsyn nodes 
						
//      for (epnt=elempnt,postn=0; (epnt=foreach(epnt,CABLE,postct,postcn,-1,-1,NULL,NULL,&pn,NULL)) 
//		      		  && postn<MAXPOSTSYN; 
//				  epnt = epnt->next, i++) {

//    use foreach loop that uses hash table to access elems
//

//      for (epnt=getepnt(postct,postcn),postn=0; (epnt=foreachn(epnt,CABLE,postct,postcn,
//				xytol, dist2d,pre_npnt,&post_npnt)) && postn<MAXPOSTSYN; 
//				  epnt = epnt->hnnext) {

      prob_pre = match_prob;

      for (postn=0, epnt=setforeachn2(postelemtyp,postct,postcn, &post_npnt); 
		       (epnt=foreachn2(epnt))!=NULL && postn<MAXPOSTSYN; epnt = epnt->hnnext) {

        pn = post_npnt->nodenm3;
        dist = dist3d(post_npnt,pre_npnt);

	 if (ninfo >= 4) ncfprintf 
		(stderr,"prect %d precn %d prenode %-d postct %d postcn %d pn %-d, dist %-7.4g xytol %-7.3g\n",
	 					    prect, precn, prenode, postct, postcn, pn, dist, xytol);

        if (dist > xytol) continue;

	c = epnt->elnum;	

	match_region = 1;
        elemtype = get_elabel(epnt,ELABL);
	region = epnt->region;			// postsynaptic region
        if (pn==postsoma && synregp >= 0 && synregp != SOMA) continue;
        if ((match_prob=check_region(postct, postcn, region, synregp, pn, "post"))==0) continue;


	// if (ninfo>=4) fprintf (stderr,"connect cell postct bptype %d region %d %d\n",bptype(postct),region,AXON);

        // if ((synregp <= 0) || 
	//	gctype(postct) &&
	//	((region!=AXON) && 
	//	 (region!=AXON_DIST) &&
        //	 (region!=AXON_THIN) &&
	//	 (region!=HILLOCK))) 

	 {

        // dist  = endist2d(epnt,pre_npnt);
        // zdist = endistzd(epnt,pre_npnt); 
        // dist  = dist2d(post_npnt,pre_npnt);
        zdist = distzd(post_npnt,pre_npnt); 
	// fprintf (stderr," nod %d dist %g\n",pn,dist);

        /* ncfprintf(stderr,"foreach loop 2: %s cable=%d, precn=%d, dist=%g zdist %g\n", 
		elemtype,c,precn,dist,zdist); /* */


        if (abs(zdist)<ztol) { /* closest? */

	/*
	  if (epnt->node2c<0) { 
	     ncfprintf (stderr,
	    "foreach loop 3: distance from [%s][%d][%d] to [%s] sphere %d [%d]= %g\n",
	      cname[prect],precn,prenode,cname[postct],c,epnt->node1c,dist); 
	  } else {
	    ncfprintf (stderr,
	    "foreach loop 3: distance from [%s][%d][%d] to [%s] cable %d [%d]-[%d]= %g\n",
	      cname[prect],precn,prenode,cname[postct],c,epnt->node1c,epnt->node2c,dist); 
	  } /* */

          postsynnod[postn].node = pn;			// remember postsyn node and cable
          postsynnod[postn].dist = dist;
          postsynnod[postn].prob = match_prob * prob_pre;
	  postsynnod[postn].epnt = epnt;		
	  postnode = pn;                         /* remenber this node */
	  postn++;
	  if (postn > MAXPOSTSYN) { postn = MAXPOSTSYN - 1; break; }
          }
        }
      }     /* foreach postsyn cable */

      n=sortdist(MAXPOSTSYN,postsynnod);	// sort the postsyn nodes by distance
      if (n>0) {				// if some postsyn nodes are within dist,
	 disttab[npostsyn*MAXPOSTSYN].node = prenode;	// save presyn node
	 for (i=0; i<MAXPOSTSYN-1; i++) { 		// save distances to postsyn nodes.
	   disttab[npostsyn*MAXPOSTSYN+i+1] = postsynnod[i];
	 }
	 if ((++npostsyn) >= MAXPOST) npostsyn = MAXPOST - 1;
      }
    }     /* foreach presyn node */


      if (ninfo >= 4) {
	  if (npostsyn==0) 
     ncfprintf (stderr,"# connect_cell: no nodes on %s %d within %s %d synapse dist, check MAXSDIST and DENZDIST\n",
				cname[postct], postcn, cname[prect], precn);
      }
    // print out the synapse table 

    if (ninfo>=4) {
      for (j=0; j<MAXPOST; j++) {
        if (disttab[j*MAXPOSTSYN].node<=soma && synregp >=0 && synregp != SOMA) continue;
        if (disttab[j*MAXPOSTSYN+1].dist>=HUGESYNDIST) continue;
        ncfprintf (stderr,"# cell precn %d presyn %d distances to postcn node:\n",precn, disttab[j*MAXPOSTSYN].node);
        for (i=1; i<MAXPOSTSYN; i++) {
  	  if (disttab[j*MAXPOSTSYN+i].dist>=HUGESYNDIST) break;
	  ncfprintf (stderr,"#  %3d %.8g\n", disttab[j*MAXPOSTSYN+i].node, 
			               disttab[j*MAXPOSTSYN+i].dist);
        }
      }	
    } /* */

    //  At this point all the close (i.e. within allowable distance) pre and postsynaptic 
    //  nodes have been identified and for each prenode, the postsynaptic nodes are sorted 
    //  according to distance.
    //
    //  Now we do a merge sort over all the synapses: find the closest pair, save it for 
    //  a synapse, then eliminate that pair from the table. Repeat this until all the pairs 
    //  are saved.
    
    npairs = 0;
    for (mindist=0; mindist<HUGESYNDIST; ) {
      mindist = HUGESYNDIST;
      for (j=0; j<npostsyn; j++) {			// find the nearest pair overall
         index = j*MAXPOSTSYN;
         if (disttab[index].node<=soma && synregp >=0 && synregp != SOMA) continue; // ignore synapses to soma
         dist = disttab[index+1].dist;
         if (dist < mindist) {
            mindist = dist;
            nodepre  = disttab[index].node;
            nodepost = disttab[index+1].node;
	    postepnt = disttab[index+1].epnt;
	    prob = disttab[index+1].prob;
	    prepos = j;
         }
      }
      if (mindist >= HUGESYNDIST) continue;
      prepost[npairs].nodepre  = nodepre;		// save the pre - post pair
      prepost[npairs].nodepost = nodepost;
      prepost[npairs].elempost = postepnt;
      prepost[npairs].dist     = mindist;
      prepost[npairs].prob     = prob;
      npairs++;
      //  fprintf (stderr," nodpre %2d nodepost %2d dist %g mindist %g\n",
      //  			nodepre,nodepost,dist,mindist);

      index = prepos*MAXPOSTSYN+1;			// remove same presyn from table
      disttab[index].dist = HUGESYNDIST;

      for (j=0; j<npostsyn; j++) {			// scan to remove all same postsyn
         index = j*MAXPOSTSYN+1;
         if (disttab[index].node==nodepost){
             disttab[index].dist = HUGESYNDIST;
	 }
       } 
    }
    if (ninfo >= 4)
        ncfprintf (stderr,"# connect_cell precn %d postcn %d npairs %d \n",precn, postcn, npairs);

    // If synanni > 0, then eliminate all closer synapses to soma in presynaptic cell
    // If synanno > 0, then eliminate all further synapses from soma in presynaptic cell
    // If synanno > 0 and synanno < synanni, then eliminate all further synapses in exclusion annulus
  
    synanni=getsv(prect,SYNANNI,rcs);
    synanno=getsv(prect,SYNANNO,rcs);
    synanpi=getsv(prect,SYNANPI,rcs);
    synanpo=getsv(prect,SYNANPO,rcs);

    if (ninfo >= 4)
        ncfprintf (stderr,"# connect_cell synanni %g synanno %g synanpi %g synanpo %g\n",
			synanni,synanno,synanpi,synanpo);

    pre_npnt_soma = ndn(prect,precn,soma);

    if (synanni>0 && synanno>0 && synanno < synanni) { // synanno < synnani sets exclusion annulus
      for (i=0; i<npairs; i++) {			// 
        node1  = prepost[i].nodepre;
        somadist = dist2d(pre_npnt_soma,ndn(prect,precn,node1));
        // fprintf (stderr,"precn %d  postcn %d synanni %g somadist %g\n",
	//	   precn,postcn,synanni,somadist); 
        if (somadist < synanni && somadist > synanno) { // If the node is in exclusion annulus,
           for (k=i; k<(npairs-1); k++) {		//   erase it
	      prepost[k] = prepost[k+1];
           }
           i--; npairs--;
           // fprintf (stderr,"npairs2 %d\n",npairs);

        } 
      }
    }
    else {						// synnano >= synanni

      if (synanni>0) { 					// synanni sets inner exclusion zone
       for (i=0; i<npairs; i++) {			// 
         node1  = prepost[i].nodepre;
         somadist = dist2d(pre_npnt_soma,ndn(prect,precn,node1));
         // fprintf (stderr,"prect %d precn %d  postcn %d synanni %g somadist %g\n",
 	 //  	   prect,precn,postcn,synanni,somadist); 
          if (somadist < synanni) {			// If the node is too close,
            for (k=i; k<(npairs-1); k++) {		//   erase it
 	      prepost[k] = prepost[k+1];
            }
            i--; npairs--;
            // fprintf (stderr,"npairs2 %d\n",npairs);
 
         } 
       }
     }
    
    // If synanno > 0, then eliminate all more distant synapses to soma in presynaptic cell
  
     if (synanno>0) {					// synanno sets outer exclusion zone
       for (i=0; i<npairs; i++) {			// 
         node1  = prepost[i].nodepre;
         somadist = dist2d(pre_npnt_soma,ndn(prect,precn,node1));
         // fprintf (stderr,"precn %d  postcn %d synanno %g somadist %g\n",
 	//	   precn,postcn,synanno,somadist); 
         if (somadist > synanno) {			// If the node is too far away,
            for (k=i; k<(npairs-1); k++) {		//   erase it
 	      prepost[k] = prepost[k+1];
            }
            i--; npairs--;
            // fprintf (stderr,"npairs2 %d\n",npairs);
 
         } 
       }
      }
    } 

    post_npnt_soma = ndn(postct,postcn,soma);

    // If synanpi > 0, then eliminate all synapses closer to soma in postsynaptic cell

    synanpi=getsv(prect,SYNANPI,rcs);
    if (synanpi>0) { 					// synanpi sets exclusion zone
      for (i=0; i<npairs; i++) {			// 
        node1  = prepost[i].nodepost;
        somadist = dist2d(post_npnt_soma,ndn(postct,postcn,node1));
        // fprintf (stderr,"precn %d  postcn %d synanpi %g somadist %g\n",
	//	   precn,postcn,synanpi,somadist); 
        if (somadist < synanpi) {			// If the node is too close,
           for (k=i; k<(npairs-1); k++) {		//   erase it
	      prepost[k] = prepost[k+1];
           }
           i--; npairs--;
           // fprintf (stderr,"npairs2 %d\n",npairs);

        } 
      }
    }

    // If synanpo > 0, then eliminate all synapses outside farther from soma in postsynaptic cell
    

    synanpo=getsv(prect,SYNANPO,rcs);
    if (synanpo>0) {					// synanpo sets exclusion zone
      for (i=0; i<npairs; i++) {			// 
        node1  = prepost[i].nodepost;
        somadist = dist2d(post_npnt_soma,ndn(postct,postcn,node1));
         // fprintf (stderr,"precn %d  postcn %d synanpo %g somadist %g\n",
	 //	   precn,postcn,synanpo,somadist); 
        if (somadist > synanpo) {			// If the node is too far away,
           for (k=i; k<(npairs-1); k++) {		//   erase it
	      prepost[k] = prepost[k+1];
           }
           i--; npairs--;
           // fprintf (stderr,"npairs2 %d\n",npairs);

        } 
      }
    }

    // If SYNRNG > 0, then eliminate all synapses outside of an angular range from presyn soma,
    //   else if -1000 < SYNRNG < 0, then ignore SYNANG, only contact dends with same postsyn orient (within -synrng),
    //   else if SYNRNG <= -1000, contact dendrites with opposite orientation (within (-synrng) - 1000).
   
    synrange = 0; 
    synrng=getsv(prect,SYNRNG,rcs);
    if (synrng != 0) {					// if range is non-zero
      if (synrng>0) {
        synang=getsv(prect,SYNANG,rcs);			// predefined angle (deg)
	synrange = synrng;
      }
      for (i=0; i<npairs; i++) {			// 
        node1  = prepost[i].nodepre;
        somangle1 = node_angle(ndn(prect,precn,node1),pre_npnt_soma); // angle from presyn node to its soma
	if (synrng<0) {					// contact dends with same orient
           node2  = prepost[i].nodepost;
           synang = node_angle(ndn(postct,postcn,node2),post_npnt_soma);// angle from postsyn node to its soma
	   synrange = -synrng;
	   if (synrng <= -1000 && synrng > -2000) {	// contact dends with opposite orient
               synang = node_angle(post_npnt_soma, ndn(postct,postcn,node2));// angle from soma to node
	       synrange -= 1000;
	   }
	   else 		// if postsyn node has syn to dend of presyn ct of same orient as precn,node1
	     if (synrng <= -2000 && synrng > -3000) { 
		  int pcn, pnod;
		  double minang;
		  synapse *synp;
		  elem *lepnt;

	       synrange -= 2000;
	       for (minang=100,lepnt=elempnt;(synp=findsyn(postct, postcn, node2, prect, -1, &lepnt)) != NULL;) {
		    pcn = synp->node2b;			// check other prect,cn that is postsyn to presyn cell
		    pnod = synp->node2c;
		    if (pcn==precn) continue;		// only check soma angle for different precn
		    synang = node_angle(ndn(prect,pcn,pnod),ndn(prect,pcn,soma));// angle from soma to node
		    if (synang < minang) minang = synang;
		}
	        synang = minang;
		if (synang > 10) { 			// skip if no synapse postct,cn,node2 -> prect
		    synrange = 1;			// deg
		    synang = somangle1 + PI;		// make outside of range; rad
	        }
	   }
	   synang *= DEG;				// convert to degrees
	}
	theta1 = (synang - synrange*0.5) / DEG;		// convert back to radians
	theta2 = (synang + synrange*0.5) / DEG;
        // fprintf (stderr,"precn %d postcn %d node2 %d theta1 %g theta2 %g somangle1 %g\n",precn,postcn,node2,theta1,theta2,somangle1);

    //  if (inrange(theta1,theta2,somangle1)) 
    //      fprintf 
    // (stderr,"precn %d postcn %d postang %g synrng %g theta1 %g theta2 %g pre soma angle %g relangle %g in %d relangl %g\n",
    //  	   precn,postcn,synang,synrange,theta1*DEG,theta2*DEG,somangle1*DEG,(synang-somangle1*DEG),
    // 	       inrange(theta1,theta2,somangle1), 180 - abs(synang-somangle1*DEG)); 

	if (!inrange(theta1,theta2,somangle1)) {	// if node's angle is outside range (radians) 
           for (k=i; k<(npairs-1); k++) {		//   erase this synapse pair
	      prepost[k] = prepost[k+1];
           }
           i--; npairs--;
           // fprintf (stderr,"npairs2 %d\n",npairs);
        } 
	//else {
        // fprintf (stderr,"precn %d  postcn %d synang %g synrng %g soma angle %g OK\n",
	//	   precn,postcn,synang,synrng,somangle1); 
	//}
      }
    }


    // Now, we eliminate synapses closer to each other than the allowed presynaptic 
    // spacing in the presynaptic dendritic arbor, starting at the closest synapse 
    // (the first one in the prepost list).

    syndist = getsv(prect,SYNSPAC,rcs);		// get the synaptic distance

    totnp = npairs;
    for (npe=j=0; j<npairs; j++) {			// 
      for (i=j+1; i<npairs; i++) {
         node1  = prepost[j].nodepre;
         node2  = prepost[i].nodepre;
	 if (node1==node2) continue;
         dist = dist3d(ndn(prect,precn,node1),ndn(prect,precn,node2));
          // ncfprintf(stderr,"# precn %d j %d i %d node1 %d node2 %d dist %g syndist %g\n", 
	  // 		 precn, j, i, node1, node2, dist,syndist);
         if (dist < syndist) {			// If the i'th node is too close,
	   for (k=i; k<(npairs-1); k++) {	//   erase it
	      prepost[k] = prepost[k+1];
	   }
	   i--; npairs--;
	   npe++;
	 } 
      }
    } 

    if ((syndist==0) && (npairs > 1)) np = 1;	// SYNSPAC = 0 -> no multiple output synapses to postsyn cn
    else                              np = npairs;

     if (ninfo>=3) 
      ncfprintf(stderr,"# connect cell %s %d to cell %s %d, syndist %g, nsynapses to make %d totnp %d npe %d\n", 
			  			cname[prect],precn,cname[postct],postcn,syndist,np,totnp,npe);

    for (this_connect=j=0; j<np; j++) {			// connect pre and post nodes 

      if (cstate) drand_setstate(cstate);
      match_region = (drand() <= prepost[j].prob);
      if (cstate) restorstate();
      if (match_region==0) continue;

      prenode  = prepost[j].nodepre;
      postnode = prepost[j].nodepost;
      postepnt = prepost[j].elempost;
      mindist  = prepost[j].dist;

      if (postepnt==NULL) break;

      if (ninfo>=3) 
	ncfprintf(stderr,"# precn %d prenode %d postcn %d postnode %d mindist %g\n", 
		      	precn, prenode, postcn, postnode, mindist);

      this_pair = connect_synapse (prect, precn, prenode, 
		  postct, postcn, postnode, postden, 
		  postepnt, zdist, mindist, 
		  growpre, growpost, usedyad, dyadtyp, dyadc, rcs, rcr, hcf, postepnt->region, postsoma);
      this_connect += this_pair;
    }

    efree(disttab);				// erase the temporary tables
    efree(prepost);

    if (np<0) np = 0;

  } 	/* else if (arbpre>=BRANCHED && arbpost>=BRANCHED) */

 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  else {  /* arbpost == NBRANCHED */

    if (arbpre==NBRANCHED) {

      if (synreg==AXON && (pre_npnt=nde(prect,precn,axtrm))!=NULL) {   	/* if presynaptic cell has axon, use it */
        // if ((prenode=celnode[prect][precn]) >= 1) prenode = axtrm;
          prenode=celnode[prect][precn];
      }
      else prenode = soma;

      if (prect==postct && prenode>=axtrm) {		/* if same ct and it has axterm, use it in post cell */
	  postnode = prenode;
      }
      else if (synregp==AXON && (post_npnt=nde(postct,postcn,axtrm))!=NULL) {  /* if postsynaptic cell has axon, use it */
             //   if ((postnode=celnode[postct][postcn]) >= 1) postnode = axtrm;
                  postnode=celnode[postct][postcn];
      } 
      else postnode = soma;

      pre_npnt  = ndn(prect, precn, prenode);
      post_npnt = ndn(postct,postcn,postnode);

      mindist = dist3d(pre_npnt,post_npnt);
      zdist   = distzd(pre_npnt,post_npnt);
      prob = prob_pre = 1;
    } 
    else if (arbpre>=BRANCHED) {

    /* branched presynaptic cells synapsing onto unbranched cell */

      mindist = zmindist = 1e10; /* now find closest dendrite on cell we're connecting to */
      c = -1;
      pn = postnode = prenode = 0;
      zdist = 0; postepnt = NULL;

      for (pre_npnt=getnpnt(prect,precn); pre_npnt=foreachn(pre_npnt,prect,precn,&prend); pre_npnt = pre_npnt->hcnext) {

        // ncfprintf(stderr,"foreach loop 1: prect=%d precn=%d prend %d\n",prect,precn,prend);

        if (prend==soma && synreg >= 0 && synreg != SOMA) continue;

       region_pre = pre_npnt->region;            // presynaptic region
       if ((match_prob=check_region(prect, precn, region_pre, synreg, prend, "pre"))==0) continue;

       synanni=getsv(prect,SYNANNI,rcs);	  // check if prenode is within specified annulus
       synanno=getsv(prect,SYNANNO,rcs);
       if (synanni>0 || synanno>0) {
	      pre_npnt_soma = ndn(prect,precn,soma);
	      somadist = dist2d(pre_npnt_soma,pre_npnt);
	      // fprintf (stderr,"somadist %g synanni %g\n",somadist,synanni);
              if (synanni>0 && somadist < synanni) continue;         // If the node is too close, skip
              if (synanno>0 && somadist > synanno) continue;         // If the node is too far, skip
              if (synanni>0 && synanno>0 && synanno < synanni) 
		      if (somadist < synanni && somadist > synanno) continue; // syanno < synanni sets exclusion zone
       }

//        for (epnt=elempnt; epnt=foreach(epnt,ELEMENT,postct,postcn,-1,-1,NULL,NULL,&pn,NULL,
//		xytol+5, dist3d, pre_npnt);
//		epnt = epnt->next) {

       prob_pre = match_prob;

//    use foreach loop that uses hash table to access elems
//
      for (postn=0, epnt=setforeachn3(postelemtyp,postct,postcn, xytol+5, dist3d, pre_npnt, &post_npnt);
		     epnt = foreachn3(epnt); epnt = epnt->hnnext, i++) {

	  c = epnt->elnum;	
          pn = post_npnt->nodenm3;
          // if (pn==soma && synregp >= 0 && synregp != SOMA) continue;
          elemtype = get_elabel(epnt,ELABL);
	  region = epnt->region;
	  if (ninfo>=5) ncfprintf(stderr,"region %s\n",region);

          pn = post_npnt->nodenm3;
	  match_region = 1;
          elemtype = get_elabel(epnt,ELABL);
          // if (pn==soma && synregp >= 0 && synregp != SOMA) continue;
          if ((match_prob=check_region(postct, postcn, region, synregp, pn, "post"))==0) continue;

        // if (bptype (postct) || (synregp < 0) ||
	// 	((region!=AXON) && 
	// 	 (region!=AXON_DIST) &&
        //       (region!=AXON_THIN) &&
	// 	 (region!=HILLOCK))) 

	  {
          // dist  = endist2d(epnt,pre_npnt);
          // zdist = endistzd(epnt,pre_npnt); 

          dist  = dist2d(post_npnt,pre_npnt);
          zdist = distzd(post_npnt,pre_npnt); 

          /* ncfprintf(stderr,"foreach loop2: cable=%d, precn=%d, dist=%g zdist %g\n", 
		c,precn,dist,zdist); /* */
          /* ncfprintf(stderr,"foreach loop3: elemtype %s postct=%d postcn=%d postnd=%d\n",
	    	elemtype,postct,postcn,pn);  /* */

          if (dist<=mindist && abs(zdist)<ztol) { /* closest? */

	   /* ncfprintf (stderr,
	      "foreach loop: distance from [%s][%d][%d] to cable %d = %g\n",
	      cname[prect],precn,prend,c,dist); /* */

	    mindist = dist;                   /* remember this dist */
	    zmindist = zdist;
	    postepnt = epnt;                       /* remember this cable */
	    postnode = pn;                         /* remenber this node */
	    prenode = prend;
	    prob = match_prob;
            //ncfprintf(stderr, "pn=%d postnode=%d; cable=%d\n",pn, postnode, postepnt->elnum);
            }
          }
        }     /* foreach post cable */
      }    /* foreach pre */

      zdist = zmindist;
      if (ninfo>=4) ncfprintf(stderr,"# connect_cell 2 mindist %g zmindist %g\n",mindist, zmindist);
      if (postepnt==NULL) break;

    }
    if (mindist > xytol || abs(zdist) > ztol) break;
 
     if (cstate) drand_setstate(cstate);
     prob *= prob_pre;
     match_region = (drand() <= prob);
     if (cstate) restorstate();
     if (match_region==0) break;
      
    this_connect = connect_synapse (prect, precn, prenode, 
		     postct, postcn, postnode, postden, 
		     postepnt, zdist, mindist, 
		     growpre, growpost, usedyad, dyadtyp, dyadc, rcs, rcr, hcf, DEND, postsoma);
  
   }

   if (ninfo >= 4) ncfprintf (stderr,"# connect_cell zdist %g ztol %g connect %d \n",zdist,ztol,this_connect);

   if (ninfo>=3) {
	  ncfprintf(stderr,"# connect_cell prect %s precn %d prenode %d postcn %d postnode %d mindist %g zdist %g\n", 
			  cname[prect],precn,prenode,postcn,postnode,mindist,zdist);
   }
   if (this_connect>0) connected++;

  }  /* for (st;;) */

 return connected;

}	/* connect_cell */

/*-----------------------------------------*/

void findnearest (int prect, int postct, int nsomas, int *nearest, int nnearest)

/* Find closest postsynaptic cells based on the location of their somas. 
   This must work for both growing and non-growing cells, so can't use 
   existing dendritic trees. The nearest array is indexed by cell numbers
   and contains the cell number and the soma node number.
   When nsomas is non-zero, then use alt_hb_somas as the list of somas.
*/

{
    int i,j,k,p,n,found,ntrials;
    int prenum, postnum, postnn;
    int precn, postcn, prearb, postarb, *pnt;
    double nearrad, rad, dist;
    node *npre, *npost, *pre_npnt_soma;
    node **presoma=NULL, **postsoma=NULL;

  ntrials = 8;

  /* get somas of presynaptic, postsynaptic cells, set in synfuncs_init() */

  prenum  = ncells[prect];
  presoma  = cell_somas[prect];
  if (nsomas==0) {
      postsoma = cell_somas[postct];
      postnum = ncells[postct];
      // nearrad = (getn(prect,AXARBDIA) + getn(postct,DTREEDIA)) * 0.5;
      nearrad = (get_prearbdia(prect) + get_postarbdia(prect,postct)) * 0.5;
  }
  else {
      postsoma = alt_hb_somas;		// use list of alternate somas (hbat's)
      postnum = nsomas;
      nearrad = (getn(prect,AXARBDIA) + getn(hbat,DTREEDIA)) * 0.5;
  } 

  if (ninfo>=3) ncfprintf(stderr,"# findnearest: nearrad %g\n", nearrad);

//  for (npre=nodepnt; npre=foreach(npre,prect,-1,soma,-1,NULL,&precn,NULL,NULL);
//  		     npre=npre->next) {
  
  for (k=0; k<prenum; k++) {			// faster way to access all somas 
    if ((npre = *(presoma+k))==NULL) continue;
    precn = npre->nodenm2;			// get cell number 
    n = 1;
    pnt = nearest+precn*(nnearest+1)*2;
    pre_npnt_soma = npre;
    for (i=0; i<ntrials; i++) {
      rad = nearrad * ((double) i+1)/ntrials;
//      for(npost=nodepnt; npost=foreach(npost,postct,-1,soma,-1,NULL,&postcn,NULL,NULL,
//      		rad, dist2d, pre_npnt_soma);
//      		npost = npost->next) {
      for (p=0; p<postnum; p++) {
	 if ((npost = *(postsoma+p))==NULL) continue;
	 postcn = npost->nodenm2;
	 postnn = npost->nodenm3;		// save soma node number (could be alt)
         dist = dist2d(pre_npnt_soma, npost);
	 if (dist > rad) continue;

         for (found=0,j=1; j<n+1; j++) {
	    if (*(pnt+(j*2))==postcn) {
		found=1; break;
	    }
         }
	 if (found) continue;			/* don't add cells already counted */
         if (prect!=postct || precn!=postcn)
            if (n <= nnearest) { 		/* if this one is close, add it */
              *(pnt+n*2)   = postcn;	/* cell numbers */
              *(pnt+n*2+1) = postnn;	/* soma node num (could be alt) */
	      n++;
	    }
      } /* foreach (postcn) */

      *(pnt+NCELLS) = n-1;		/* number of cells */

      if (ninfo>=4) ncfprintf(stderr,"# findnearest: trial %d  rad %-6.4g ncells %d\n", i, rad, n-1);

     } /* for (i=0; i<ntrials; ) */
   } /* foreach (precn) */
}

/*-----------------------------------------*/

static char *sstate = NULL;

void shuffle(int *arr, int ncel)

/* make unique disordered index in array "arr", of size ncel+1 */

{
   int c,s,temp;

  for (c=0; c<ncel; c++)	/* initialize arr */
    arr[c] = c;	

  if (sstate) drand_setstate(sstate); 
  if(ncel>1){			/* if >1 cell, shuffle order of cells */
    for (c=0; c<ncel; c++) {
      for (s=c; c==s; s)  		/* make source diff than dest */
        s = int(drand() * ncel);
      temp = arr[c];			/* swap for shuffle */
      arr[c] = arr[s];
      arr[s] = temp;			/* */
      /* ncfprintf(stderr,"c=%d, arr[c]=%d\n", c, arr[c]);	/* */
    }
  }
  if (sstate) restorstate();
}

/*-----------------------------------------*/

double integ_circ (double x) 

/* return the fractional area of a circle outside of a limit = x */
/* x = limit up to radius of circle,  must range from 0 - 1 */
/* returns area normalized to area of whole circle */
/* ref:	https://math.stackexchange.com/questions/1469820/area-under-quarter-circle-by-integration */

{
    double z;
    z = 1.0 - x;
    return (0.5 - (z * sqrt(1-z*z) + asin(z)) / PI);
}

/*-----------------------------------------*/

double calc_conv_frac(node *np)

/* Scale the convergence from presynaptic cell according to
    how much of the dendritic tree is inside "arrsiz".  */

/* Uses single integral, need double integral, but approx correct */

{
     int ct2;
     double cellrad, retval,retval2;
     double xloc, yloc;
     double xlim, ylim;
     double xcent, ycent;
     double xinteg, yinteg;
     double xfrac_out, yfrac_out;

   ct2 = np->nodenm1;
   if ((ct2==ha) || (ct2==hb)) {

      if (!notinit(conearrsiz) || !notinit(conexarrsiz)) {
          xlim = ylim = conearrsiz / 2; 
	  if (notinit(conearrsiz)) {
              xlim = conexarrsiz / 2; 
              ylim = coneyarrsiz / 2; 
	  }
	  xcent = ycent = 0;
          if (!notinit(conearrcentx)) {		// offset for center of cone array
	     xcent = conearrcentx;
	     ycent = conearrcenty;
	  }
      } else {
          xlim = xarrsiz / 2;
          ylim = yarrsiz / 2;
	  xcent = ycent = 0;
      }
      xloc = np->xloc - xcent;
      yloc = np->yloc - ycent;
      if (xloc<0) xloc = -xloc;
      if (yloc<0) yloc = -yloc;
      cellrad = getn(ct2,DTREEDIA) / 2;
      if (cellrad == 0) cellrad = 1;
      xfrac_out = (xloc+cellrad-xlim) / (cellrad*2);	// fraction of dia outside limit
      yfrac_out = (yloc+cellrad-ylim) / (cellrad*2);
      if (xfrac_out > 0.5) xfrac_out = 0.5;
      if (yfrac_out > 0.5) yfrac_out = 0.5;
      // if (xfrac_out>0) retval =  1 - xfrac_out; else retval = 1.0;
      // if (yfrac_out>0) retval *= 1 - yfrac_out; 
      if (xfrac_out>0) xinteg = integ_circ(xfrac_out*2); else xinteg = 0;
      if (yfrac_out>0) yinteg = integ_circ(yfrac_out*2); else yinteg = 0;
      retval = (1-xinteg) * (1-yinteg);

       // fprintf (stderr,"x %g y %g cellrad %g ylim %g xfrac_out %g yfrac_out %g xinteg %g yinteg %g retval %g\n",
       // 	      xloc,yloc,cellrad,ylim,xfrac_out,yfrac_out,xinteg,yinteg,retval);
   } else retval = 1.0;
   return retval;
}

/*-----------------------------------------*/

void connect_types (int ctype1,int ctype2, int nsomas, int makercr, int dyadc) 

/* Connect one type of cell to another type, with
competition between cells of same type. Assume that
photoreceptors exist and skeleton bipolar cells already exist in
correct number and location.  The photoreceptors are checked in
random sequence.  For each photoreceptor, find all the bipolar
cells whose dendritic fields it lies in, then pick one of them at
random to connect to.  The probability of connecting to a cell is
proportional to its weighting function, so a stochastic gaussian
weighting function is generated to decide which cells connect.
If "growpost" is set, extend a dendrite from the postsynaptic cell
to the presynaptic cell.  

 If "usedyad" is set, don't make any new synapses, but find
  previously made synapses from the same type of photoreceptor
  to a cell of "dyadtyp" and make the new synaptic connections
  be dyads to the previously made synaptic connections.  Thus,
  photoreceptors -> bipolar cells are normal synapses, and
  photoreceptors -> HCs are dyads.

 If "makercr" is set, allow reciprocal feedback to be made at the 
 same node where a feedforward connection is made.

*/

{
    int s,i,ct1,ct2,cn,c1n,c2n,n,c1,c2,rcs,c2soma;
    int trials, maxtrials, maxcells, ncheck;
    int isgj,nconn,oldn,oldn2,ncel1,ncel2;
    int arbpre, arbpost, srseed;
    int growpre, growpost, usedyad, dyadtyp;
    int conv,div,postden,synnum;
    int t=1, conn_cnt1=0, conn_cnt2=0, conn_cnt3=0;
    double radmin, radweight;
    double rdist,ctrad;
    double dens,rad;
    double arb_scale1,arb_scale2;
    int *nearest;
    double *conv_frac=NULL;
    int *ckd1, *ckd2;
    node **postsoma;
    int postnum;

  ct1 = ctype1;			/* presynaptic cell type */
  ct2 = ctype2; 		/* postsynaptic cell type */
  				/* Don't make synaptic connections if the cells are missing. */
  if (!getn(ct1,MAKE) || !getn(ct2,MAKE)) return; 

  synfuncs_init();

  arbpre  = get_arbpre(ct1);		/* dend (bp) or axon (am) arborization type: NBR, BR, etc. */
  arbpost  = get_arbpost(ct1,ct2);	/* dend (bp) or axon (bp) from amx type: NBR, BR, etc. */
  // maxcells = (int)getn(ct1,MAXNUM);
  maxcells = ncells[ct1];

  /* Here we calculate how many postsynaptic cells to check */

  if (arbpre>NBRANCHED && arbpost==NBRANCHED) { /* area(pre) * density(post) */
     dens = getn(ct2,DENS) * 1e-6;
     rad = getn(ct1,DTREEDIA)*0.5;
     ncheck = int(dens * rad * rad * PI * 4) + 1;
     if (ninfo >=4) ncfprintf (stderr,"# 1 dens %g rad %g\n",dens,rad);
  }
  else if (arbpost > NBRANCHED) {		/* area(post) * density(post) */
     dens = getn(ct2,DENS) * 1e-6;
     // rad = getn(ct2,DTREEDIA)*0.5;
     rad = (get_prearbdia(ct1) + get_postarbdia(ct1,ct2))*0.5;
     ncheck = int(dens * rad * rad * PI * 2) + 1;
     if (ninfo >=4) ncfprintf (stderr,"# 2 dens %g rad %g\n",dens,rad);
  }
  else if (arbpre==NBRANCHED) {		/* area(post) * density(post) */
     dens = getn(ct2,DENS) * 1e-6;
     rad = getn(ct2,DTREEDIA)*0.5;
     ncheck = int (dens * rad * rad * PI * 2) + 1;
     if (ninfo >=4) ncfprintf (stderr,"# 3 dens %g rad %g\n",dens,rad);
  }
  else ncheck  = (int)(getn(ct1,MAXCOV)*getn(ct2,MAXCOV));
  if (ninfo >=4) ncfprintf (stderr,"# ct1 %s ct2 %s ncheck %d\n",cname[ct1],cname[ct2],ncheck);

  nearest = (int *)emalloc((maxcells+1)*(ncheck+1)*2*sizeof(int));
  ckd1    = (int *)emalloc((maxcells+1)*sizeof(int));
  ckd2    = (int *)emalloc((ncheck+1)*sizeof(int));
  for (i=1; i<=maxcells; i++) ckd1[i] = -1;
  for (i=1; i<=ncheck;   i++) ckd2[i] = -1;
  
  ncel1 = ncells[ct1];		/* number of cells of presynaptic type */
  ncel2 = ncells[ct2];		/* number of cells of presynaptic type */

  /* set up random state for shuffle() and prob_gaus */

  srseed = noduniq(ctype1,ncel1,ctype2,ncel2, rseed);
  sstate = (char *)makrand (RNDSIZ,"connect_types",ctype1);
  if (sstate) ncinitstate(srseed,sstate,RNDSIZ);

  findnearest(ct1,ct2,nsomas,nearest,ncheck);	/* find nearby postsyn cells */
  rcs = getconn(ct1,ct2);	/* synapse type prect -> postct */

  if (ninfo >= 2) {
    if (ncel1==1) {
        if (ncel2==1) ncfprintf(stderr,"# connecting %s to %s\n", cname[ct1], cname[ct2]);
        else          ncfprintf(stderr,"# connecting %s to %ss\n", cname[ct1], cname[ct2]);
    } else {
        if (ncel2==1) ncfprintf(stderr,"# connecting %ss to %s\n", cname[ct1], cname[ct2]);
        else          ncfprintf(stderr,"# connecting %ss to %ss\n", cname[ct1], cname[ct2]);
    }
  } else 

  if (ninfo >= 3) {
    if (ncel1==1) {   
        if (ncel2==1) ncfprintf(stderr,"# connecting %d %s to %d %s\n", ncel1, cname[ct1], ncel2, cname[ct2]);
        else          ncfprintf(stderr,"# connecting %d %s to %d %s's\n", ncel1, cname[ct1], ncel2, cname[ct2]);
    } else {
        if (ncel2==1) ncfprintf(stderr,"# connecting %d %ss to %d %s\n", ncel1, cname[ct1], ncel2, cname[ct2]);
        else          ncfprintf(stderr,"# connecting %d %ss to %d %ss\n", ncel1, cname[ct1], ncel2, cname[ct2]);
    }
  }

  growpre  = (int)getsv(ct1,GROWPRE,rcs); 	/* grow presynaptic cell ? */
  usedyad  = (int)getsv(ct1,USEDYAD,rcs); 	/* use dyad when presynaptic cell ? */
  dyadtyp  = (int)getsv(ct1,DYADTYP,rcs); 	/* dyad type to attach */
  growpost = (int)getcv(ct2,GROWPOST,(int)getsv(ct1,CONPOST,rcs));/* grow postsyn cell?*/

  conv   = (int)getcv(ct2,CELCONV,(int)getsv(ct1,CONPOST,rcs));
  div    = (int)getsv(ct1,CELDIV,rcs);
  isgj   = ((int)getsv(ct1,SRESP,rcs)==xgapj);
  synnum = (int)getsv(ct1,SYNNUM,rcs);

  if (isgj) {	/* for gj connections, start at 0.5 radius */
    radmin = 0.5;
  }
  else {	/* for all other interactions, start a little closer */
    radmin = 0.2;
  };

   if (growpost) {	// calculate area of dendritic trees inside array
      if (nsomas==0) {
          postsoma = cell_somas[ct2];
          postnum = ncells[ct2];
      }
      else if (ct2==hb) {
          postsoma = alt_hb_somas;		// use list of alternate somas (hbat's)
          postnum = nsomas;
      } else fprintf (stderr,"alt soma for non-hb cell\n");
      // fprintf (stderr,"ct2 %d n %d\n",ct2,postnum);
      conv_frac = (double *)emalloc((postnum+1)*sizeof(double));
      for (i=0; i<=postnum; i++) { conv_frac[i] = 1.0; }
      if (ct2==ha || ct2==hb) {	// calculate area of dendritic trees inside array
         for (i=0; i<postnum; i++) {
	    conv_frac[i+1] = calc_conv_frac(postsoma[i]);	// when ha, hb outside arrsiz, scale conv
	    // fprintf(stderr,"ct2 %d cn2 %d x %g y %g conv_frac %g\n",
	    // 		ct2,i,postsoma[i]->xloc,postsoma[i]->yloc,conv_frac[i+1]);
         }
      }
   }

  	/* Eliminate cells that are too far away, before we check */
	/* carefully. Note that if DTREEDIA does not control size */
  	/* of dendritic tree (i.e. with morpology file), its */
	/* only use is to set this criterion for selection.*/

   arb_scale1 = arb_scale2 = 1.0;
   //arb_scale1 = getn(ct1,ARBSCALE);
   //if (arb_scale1==0) arb_scale1 = 1.0;
   //arb_scale2 = getn(ct2,ARBSCALE);
   //if (arb_scale2==0) arb_scale2 = 1.0;
   ctrad = (getn(ct1,DTREEDIA)*arb_scale1 + getn(ct2,DTREEDIA)*arb_scale2 + 
	    getn (ct1,MAXSDIST) + getn(ct2,MAXSDIST)) * 0.5; /* radius to check */
  if (ninfo==2) ncfprintf (stderr,"# ");
  oldn = oldn2 = 0;
  if (!growpost) maxtrials = 1;		/* no random synapse selection, only 1 trial */
  else		 maxtrials = 15;
  for (trials=0; trials<maxtrials; trials++) { /* trials to connect */

    /* The idea here is to find the nearest connections first,
       which will tend to reliably connect the cells near the center 
       of the dendritic field, then stepwise increase the radius of 
       the Gaussian weighting function to allow the remainder of the
       connections to be made more randomly.
    */

   radweight = ((double)(trials+1))/maxtrials * (1-radmin) + radmin;

   shuffle (ckd1,ncel1);

   nconn = 0;
   for (c1=0; c1<ncel1; c1++) {		/* For each presyn cell */

     c1n = cellnums[ct1][ckd1[c1]+1];		/* shuffle number for presyn cell*/

     ncel2 = *(nearest+(c1n*(ncheck+1)+NCELLS)*2); // indexed by cell numbers 
     shuffle (ckd2,ncel2);

      /* Find all postsyn cells that overlap with presynaptic cell */
      /*  and therefore could potentially connect. */

     for (c2=0; c2<ncel2; c2++) {		// index into number of cells, not cell numbers
	  int ok,nsyn;

       c2n = *(nearest+(c1n*(ncheck+1)+ckd2[c2]+CELN)*2);     /* index into nearby cell2s */
       if (ct1==ct2 && c1n==c2n) continue;   /* don't connect cell to itself */

       c2soma = *(nearest+((c1n*(ncheck+1)+ckd2[c2]+CELN)*2)+1);   /* get soma node num (could be alt) */

       rdist = dist2d(ndn(ct1,c1n,soma),ndn(ct2,c2n,c2soma)); /* find radial dist*/

       if (ninfo>=4) {
     ncfprintf (stderr,
	"# connect_types %s %d to %s %d rdist %4.3g, ctrad %g radweight %g radcrit %g growpre %d growpost %d\n",
          cname[ct1],c1n,cname[ct2],c2n,rdist,ctrad,radweight,ctrad*radweight,growpre,growpost);
	 if (connected(ct1,c1n,ct2,c2n))
           ncfprintf (stderr,"# %s %-3d and %s %-3d are already connected\n", 
			cname[ct1],c1n,cname[ct2],c2n);
       }

       if (growpost) {
	     int prob_gaus;

         if (rdist <= ctrad && !connected(ct1,c1n,ct2,c2n) && 
           (ncel_in(ct2,c2n,ct1) < conv * conv_frac[c2n]) && (ncel_out(ct1,c1n,ct2) < div)) {

	    if (sstate) drand_setstate(sstate);
            prob_gaus = (drand() < gauss(rdist,ctrad*radweight));
	    if (sstate) restorstate();

            if (prob_gaus) {
            for (nsyn=0; nsyn<synnum; nsyn++) {
                 ok = connect_cell (ct1,c1n,ct2,c2n,c2soma,usedyad,dyadtyp,dyadc,makercr,
					growpre,growpost,postden=pickden[ctype2]);
            }
            if (ok) { nconn++; 
		    if (ninfo>=2) { 
		        if (ninfo==2) {
			   if (++conn_cnt1<10) ncfprintf (stderr,".");
			   conn_cnt2 += t;
		    	   if (conn_cnt2>=10)  { conn_cnt2=0; if (conn_cnt3 < 99) ncfprintf (stderr,"+");}
			   if (++conn_cnt3>=100) { t=0; conn_cnt3=0; ncfprintf (stderr,"*");}
			} else ncfprintf (stderr,".");
		    } 
	    }

           } /* if (drand()<gauss) */
         } /* if (rdist<=ctrad) */
       }  /* if (growpost) */

       else { /* if (!growpost) */
         if (rdist <= ctrad*radweight && !connected(ct1,c1n,ct2,c2n) && 
           (ncel_in(ct2,c2n,ct1) < conv) && (ncel_out(ct1,c1n,ct2) < div)) {

            for (nsyn=0; nsyn<synnum; nsyn++) {
              ok = connect_cell (ct1,c1n,ct2,c2n,soma,usedyad,dyadtyp,dyadc,makercr,
					growpre,growpost,postden=pickden[ctype2]);
            }
            if (ok) { nconn++; 
		    if (ninfo>=2) { 
			if (ninfo==2) {
		           if (++conn_cnt1<10) ncfprintf (stderr,".");
			   conn_cnt2 += t;
		    	   if (conn_cnt2>=10)  { conn_cnt2=0; if (conn_cnt3 < 99) ncfprintf (stderr,"+");}
			   if (++conn_cnt3>=100) {t=0; conn_cnt3=0; ncfprintf (stderr,"*");}
			} else ncfprintf (stderr,".");
		    } 
	    }
         }
       }
       if (ncel_out(ct1,c1n,ct2) >= div) break;

      }; /* for (c2;;)  For all postsyn cells */
    };  /* for (c1;;) */
   if (ninfo >= 3) ncfprintf (stderr,"# trial %d nconn %d\n",trials,nconn);
   if (nconn==0 && oldn==0 && oldn2==0 && trials>9) break;
   oldn2 = oldn;
   oldn = nconn;
  }  /* for (trials;;) */

  efree(nearest);		/* free the local arrays */
  efree (ckd1);
  efree (ckd2);
  if (growpost) efree (conv_frac);
  if (ninfo>=2) ncfprintf (stderr,"\n");
  
  if (save_syn_file!=NULL) fflush (save_syn_file);

} /* proc connect_types() */

/*-----------------------------------------*/
void connect_types (int ctype1, int ctype2, int makercr=0, int dyadc=0)

{
      int nsomas;

    connect_types (ctype1, ctype2, nsomas=0, makercr=0, dyadc=0);

}
/*-----------------------------------------*/

void connect_types (int ctype1,int ctype2) 

{
   int dyadc, makercr;

  connect_types (ctype1, ctype2, makercr=0, dyadc=0);
}

/*-----------------------------------------*/

void make_synapse_type (node *prenod, node *postnod, int syntype) {

	int rcs,stype,nfilt,ntrans;
        int prect, postct, postcn, postnode;
        cattrib *capnt;
        nattrib *napnt;
	double scond;
        synapse *spnt;

   if (prenod==NULL || postnod==NULL) return;

   prect  = prenod->nodenm1;
   postct = postnod->nodenm1;
   postcn = postnod->nodenm2;
   postnode = postnod->nodenm3;
   rcs=getconns(prect,postct,syntype);	// get connection number */
   scond = setcondmul(gs(SCMUL), gs(SCGRAD), gs(SEGRAD), gs(SCOND),  postct, postcn, postnode);
    
   spnt = make_synapse(prenod, postnod, rcs);
   spnt->ntact  = OPEN;
   spnt->maxcond= scond;
   spnt->casens = gs(SENSCA);
   spnt->caegain = gs(CAEGAIN);
   spnt->curve  = (gs(SVGAIN)<=0 ? LINEAR: EXPON);
   spnt->ngain  = gs(SGAIN);
   spnt->vgain  = (gs(SVGAIN)<=0 ? 1 : gs(SVGAIN));
   spnt->thresh = gs(STHRESH);
   spnt->vrev   = gs(SVREV);
   nfilt = (int)gs(SNFILTH);
   spnt->nfilt1h= (short int)nfilt;
   spnt->timec1h= makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
   spnt->filt1hg= gs(SHGAIN);
   spnt->filt1ho= gs(SHOFFS);
   nfilt        = (int)gs(SNFILT);
   spnt->nfilt2 = (short int)nfilt;
   spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
   spnt->tfall2 = gs(SFALL);
   spnt->rrpool = gs(SRRPOOL);
   spnt->rrpoolg = gs(SRRPOOLG);
   spnt->mrrpool = gs(SMRRPOOL);
   spnt->maxsrate = gs(SMAXRATE);
   spnt->vsize  = gs(SVSIZ);
   if (gs(SVNOISE)>0) {
     napnt=make_vesnoise(spnt);
     napnt->vsize =gs(SVSIZ);
     napnt->cov   =gs(SCOV);
   }
   spnt->trconc = gs(STRCONC);
   switch (int(gs(SRESP))) {
       case XGABA:  ntrans = GABA; stype = 1; break;
       case XGABA2: ntrans = GABA; stype = 2; break;
       case XGABA3: ntrans = GABA; stype = 3; break;
       case XGABA4: ntrans = GABA; stype = 4; break;
       case XGLY:   ntrans = GLY;  stype = 1; break;
       case XAMPA:  ntrans = AMPA; stype = 1; break;
       case XAMPA1: ntrans = AMPA; stype = 1; break;
       case XAMPA2: ntrans = AMPA; stype = 2; break;
       case XAMPA5: ntrans = AMPA; stype = 5; break;
       case XKAINATE: ntrans = KAINATE; stype = 1; break;
       case XMGLUR6: ntrans = GLU; stype = 1; break;
   }
   capnt = (cattrib*)chanattr(spnt,ntrans);
   capnt->stype = stype;
   if (gs(SCNOISE)>0) {
     napnt=make_chnoise(spnt);
     //napnt->unitary=22e-12;
   }
}

