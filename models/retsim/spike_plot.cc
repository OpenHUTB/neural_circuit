/* module spike_plot.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

/* Sets up a spike frequency plot that can be accessed by calling
   the "spike_plot()" procedure:

   plot_freq = 1;
   if (plot_freq) plot(spike_plot, max=500, min=0, vpen=freq_color);
*/

/* Process intracellular voltage to find spike times, */
/* inter-spike intervals, and instantaneous spike */
/* frequencies. */

#define MAX_SPIKE_PLOTS 5

int spikyet[MAX_SPIKE_PLOTS]	    = {0, 0, 0, 0, 0};
double oldvh[MAX_SPIKE_PLOTS]	    = {0, 0, 0, 0, 0};
double oldspiktim[MAX_SPIKE_PLOTS]	    = {0, 0, 0, 0, 0};
double max_spike_rate[MAX_SPIKE_PLOTS] = {0, 0, 0, 0, 0};
int spike_node_map[MAX_SPIKE_PLOTS][3] = {-1};
int num_spike_nodes = 0;

#define SPIKTHRESH (-0.04)
#define _CT 0 
#define _CN 1 
#define _N  2 


double freq_color (int nplot, double xval, double yval) 

{
   double retval;

  if (yval > 5) retval = 15;
  else if (yval > 20)  retval = 14;
  else                 retval = 12;
  if (yval < 1) retval = -1;
  return (retval);
}

int add_spike_node(int ct, int cn, int n)
{
  if (num_spike_nodes < MAX_SPIKE_PLOTS) {
    spike_node_map[num_spike_nodes][_CT] = ct;
    spike_node_map[num_spike_nodes][_CN] = cn;
    spike_node_map[num_spike_nodes][_N]   = n;
    //fprintf(stderr, "# Added node %d to map with index %d\n", n, num_spike_nodes);
    num_spike_nodes++;
    return (num_spike_nodes - 1);
  } else {
    fprintf(stderr, "# Could not add spike node %d\n", n);
  }
  return -1;
}

int get_spike_node_index(int ct, int cn, int n)
{ 
   int i;

  for (i = 0; i < num_spike_nodes; i++) {
    if ((spike_node_map[i][_CT] == ct) &&
        (spike_node_map[i][_CN] == cn) &&
        (spike_node_map[i][_N] == n)) {
      return i;
    }
  }
  return -1;
}

double get_max_spike_rate(int ct, int cn, int n)
{
  int i;
  i = get_spike_node_index(ct,cn,n);
  if (i != -1) {
    return max_spike_rate[i];
  }
  return 0;
}

double spike_plot(double cct, double ccn, double nodn, double thresh, double time)

/* calculate instantaneous frequency from inter-spike interval */

{
  double spikerate, spikint, spiktim;
  double vh;
  int index, ct, cn, nodenum;

  ct = int(cct);
  cn = int(ccn);
  nodenum = int(nodn);
  //fprintf(stderr, "#Plotting spikerate of node %g\n", nodenum);
  spikerate = 0;

  if (nodenum != -1) {
    spikerate = 0;
    index = get_spike_node_index(ct,cn,nodenum);
    if (index == -1) {
      index = add_spike_node(ct,cn,nodenum);
    };
    //fprintf(stderr, "#Got index for node %g as %g\n", nodenum, index);
    if (index > -1) {
      vh = (v(ndn(ct,cn,nodenum)) > SPIKTHRESH);
      if (vh && !oldvh[index]) {           /* spike here */
        spiktim = time;

        if (spikyet[index]) {
        
          spikint = spiktim - oldspiktim[index] + 1e-12;
          spikerate = 1 / spikint;
        } else {
          spikyet[index] = 1;
          spikerate = 0;
        };

        oldspiktim[index] = spiktim;
      }
      else { spikerate = 0; }
      oldvh[index] = vh;

      if (spikerate > max_spike_rate[index]) {
        max_spike_rate[index] = spikerate;
      }
    }    
  }
  return spikerate;
}

#undef _CT 
#undef _CN 
#undef _N 

