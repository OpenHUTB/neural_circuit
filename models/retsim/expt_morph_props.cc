/* Experiment morph_props for program retsim */

#include <stdlib.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "morphfuncs.h"

const char *expt_desc;			//description of experiment, used for file naming
const char *celltype; 			//name of cell type
int cellnum; 				//cell number
const char *morphfile;			//path to morphology file
const char *outname;			//description for output file
const char *delim;			//delimiter to separate fields


void defparams()
{
	expt = "morph_props";

	setptr("expt_desc", &expt_desc);
	setptr("celltype", &celltype);
	setptr("cellnum", &cellnum);
	setptr("morphfile", &morphfile);
	setptr("outname", &outname);
	setptr("delim", &delim);
}


void setparams() 
{

	make_rods  = 0;
	make_cones = 0;

	make_ha    = 0;
	make_hb    = 0;
	make_hbat  = 0;
	make_dbp1  = 0;
	make_dbp2  = 0;
	make_rbp   = 0;
	make_gca   = 0;
	make_gcb   = 0;
	make_dsgc  = 1;

	timinc = 5e-6;
}


void runexpt()
{
	if (notinit(expt_desc)) expt_desc = "_default";
	if (notinit(celltype))  celltype  = "dsgc";
	if (notinit(cellnum))   cellnum   = 1;
	if (notinit(morphfile)) morphfile = "runconf/morph_ds1e";
	if (notinit(outname))   outname   = "ds1e";
	if (notinit(delim))     delim     = " ";

	int ctnum = find_ct(celltype);
	if (ctnum == -1) {
		fprintf(stderr, "Cell type '%s' not found, exiting...\n", celltype);
		exit(1);
	}

	//obtain the morphological properties of each node for the
	//cell specified and write to a file
	vector<morph_sample *> *samps = get_morph_props(ctnum, cellnum);
	fprintf(stderr, "# Sampled from %d nodes\n", (int)samps->size());

	char *fname = new char[1024];
	sprintf(fname, "morph_props_%s.data", outname);

	ofstream outfile(fname);  

	node *n;

	for (int k = 0; k < samps->size(); k++) {

		//fprintf(stderr, "k=%d: system=%d, order=%d, dist=%g, edist=%g, cdm=%g, dia=%g, sa=%g, lambda=%g, dist_sa=%g\n",
		//    k, (*samps)[k]->system, (*samps)[k]->order, (*samps)[k]->node2soma_dist, (*samps)[k]->node2soma_edist,
		//    (*samps)[k]->cdm,(*samps)[k]->dia, (*samps)[k]->surf_area, (*samps)[k]->space_const,
		//    (*samps)[k]->dist_surf_area);

		n = nd(ctnum, cellnum, k);


		outfile << k << delim << n->xloc << delim << n->yloc << delim << n->zloc << delim
		<< (*samps)[k]->system << delim << (*samps)[k]->order << delim
		<< (*samps)[k]->node2soma_dist << delim << (*samps)[k]->node2soma_edist << delim
		<< (*samps)[k]->cdm << delim << (*samps)[k]->dia << delim << (*samps)[k]->surf_area
		<< delim << (*samps)[k]->dist_surf_area << delim << (*samps)[k]->space_const << endl;
	}

	fprintf(stderr, "# Wrote morphology properties to %s\n", fname);
	outfile.close();
	delete fname;  

	/*       
  // Compute node-to-node distances for a given morphology,
  //  a fairly time-consuming process for large morphologies.
 
  vector <vector<double>*> *n2n_dists = get_node2node_dists(ctnum, cellnum);
  char *fname = new char[300];
  sprintf(fname, "node2node_%s.data", morphfile);

  ofstream outfile(fname);  

  int numnodes = n2n_dists->size();

  for (int k = 0; k < numnodes; k++) {
    for (int j = 0; j < numnodes; j++) {
      outfile << (*(*n2n_dists)[k])[j];
      if (j != (numnodes-1)) outfile << " ";
    }
    outfile << endl;
  }

  outfile.close();
  delete fname;
	 */
}

