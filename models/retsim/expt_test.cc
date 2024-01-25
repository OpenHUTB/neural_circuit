#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ncfuncs.h>
#include "retsim.h"
#include "retsim_var.h"

double testvar1;
double testvar2;
const char *testvar3;

void defparams(void) 
{
  setptr("testvar1", &testvar1);
  setptr("testvar2", &testvar2);
  setptr("testvar3", &testvar3); 
}


void setparams(void)
{
  make_rods = 0;
  make_cones= 0;        /* make cones, cbp, gc */
  make_ha   = 0;
  make_hb   = 0;
  make_hbat = 0;
  make_dbp1 = 0;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_gca  = 0;
  make_gcb  = 0;
  make_dsgc = 0;
}



void runexpt(void)
{
  if (notinit(testvar1)) testvar1 = -1;
  if (notinit(testvar2)) testvar2 = -1;
  if (notinit(testvar3)) testvar3 = "uninit";

  fprintf(stderr, "testvar1=%g, testvar2=%g, testvar3=%s\n", testvar1, testvar2, testvar3);
}

