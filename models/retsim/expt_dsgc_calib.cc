/* Experiment dsgc_calib for program retsim */
 
#include <cstring>
#include <string>
using namespace std;

#include "rectask.h"
//#include "densfuncs.h"

#include "onplot_dsgc_movie.cc"


const char* expt_type;			//type of experiment
const char* stimtype;			//type of stimulus ("vclamp" | "cclamp")
int    samprate;			//sample rate (in Hz)
double stimamp;				//amplitude of stimulus, in Volts or Amps
double stimdur;				//duration of stimulus
double starttime;			//start time of stimulus
double endwait;				//time to wait after stimulus before exit
const char* expt_desc;			//description of experiment, used for file naming

int rec_ct;
int rec_cn;

void defparams()
{
  expt = "dsgc_node";

  setptr("expt_type", &expt_type);
  setptr("stimtype", &stimtype);
  setptr("stimamp", &stimamp);
  setptr("stimdur", &stimdur);
  setptr("starttime", &starttime);
  setptr("endwait", &endwait);
  setptr("expt_desc", &expt_desc);
  setptr("ath_dia", &ath_dia);
  setptr("samprate", &samprate);

  defparams_dsgc_movie();
  defparams_onplot_movie();
}


void setparams() 
{
  make_rods = 0;
  make_cones= 0;
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 0;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 0;
  make_gcb  = 0;
  make_dsgc = 1;

  rec_ct = dsgc;  /* type of cell to record from */
  rec_cn = 1;     /* cellnum to record from */

  if (notinit(ath_dia)) ath_dia = .5;
  
  /*
  setn(dsgc, MAKE_DEND, 0);
  setn(dsgc, MAKE_AXON, 0);
  setn(dsgc, MAKE_LONG, 0);
	*/

  onplot_dsgc_movie_init();
  onplot_movie_init(); 

  timinc = 5e-6;
  
  
}


void runexpt()
{
  if (notinit(stimtype))  stimtype   = "cclamp";
  if (notinit(stimamp))   stimamp    = .3e-9;
  if (notinit(stimdur))   stimdur    = .05;
  if (notinit(starttime)) starttime  = .005;
  if (notinit(endwait))   endwait    = .01;
  if (notinit(space_time)) space_time  = 1;
  if (notinit(make_movie)) make_movie = 1;
  if (notinit(expt_desc))  expt_desc = "_default";

  //apply_regparams(dsgc, 1, "runconf/regparams_dsgc_DS060825");

  if (space_time) {    
    if (strcmp(stimtype, "cclamp") == 0) {      
      plot_v_nod(dsgc, 1, soma, -.11, .02, blue, "V[soma]", 1, 0.3); 
    } else if (strcmp(stimtype, "vclamp") == 0) {
      plot_i_nod(dsgc, 1, soma, -1e-9, 1e-9, blue, "I[soma]", 1, 0.5); 
    }
  }

  node *somaNode = nd(dsgc, 1, soma);
  
  //schedule the movie display to run every 25us
  sched.addTask(25e-6, onplot_movie);

  //create and schedule a task to record voltage every 50us
  char fname[1024];
  sprintf(fname, "vtrace_%s.data", expt_desc);
  fprintf(stderr, "# Writing voltage trace to '%s'\n", fname);

  rectask rt(fname, ' ');
  rt.add(somaNode, "v");
  double sampint = (double) 1.0 / samprate;
  fprintf(stderr, "# Sample Rate: %gkHz, Sample Interval: %gus\n", samprate*1e-3, sampint*1e6);
  sched.addTask(sampint, &rt);  
  
  setxmin = 0;

  //run experiment
  endexp= starttime + stimdur + endwait;

  step(starttime);

  fprintf(stderr, "############\n");
  fprintf(stderr, "# Model Info\n");
  fprintf(stderr, "############\n");
  fprintf(stderr, "# Stepsize: %gus, expt_desc=\"%s\"\n", timinc*1e6, expt_desc);
  fprintf(stderr, "# Temperature=%g, dqm=%g, dqh=%g, vna=%g, vk=%g\n", tempcel, dqm, dqh, vna, vk);  
  fprintf(stderr, "# Experiment duration: %gms\n", endexp*1000);
  fprintf(stderr, "############\n");

  if (strcmp(stimtype, "cclamp") == 0) {
    fprintf(stderr, "# Current injection of %gpA for %gms\n", stimamp*1e12, stimdur*1e3);
    cclamp(somaNode, stimamp, starttime, stimdur);
  } else if (strcmp(stimtype, "vclamp") == 0) {
    fprintf(stderr, "# Voltage clamp at %gmV for %gms\n", stimamp*1e3, stimdur*1e3);
    vclamp(somaNode, stimamp, starttime, stimdur);
  }

  step(stimdur + endwait);
};
