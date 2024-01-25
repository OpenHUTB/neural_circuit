/* 	Program ttest for t-test of 2 lists of numbers */

/*   	Latest mod	23-Feb-18	R.G.Smith

        Accepts 2 columns of numbers for pair-wise t-test
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *pictin;

int indep_flg=0;
int paired_flg=1;
static char *filnam = 0;

double sqrt(double arg);
double atof();
void run();

/****************************************/

int main(int argc, char **argv)

{
   int i;
   FILE *temp,*freopen();
   char *cptr=0;
	 
 pictin = stdin;
 if (argc==1)			/* if user needs help */
   run();
 else
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {

          case 'i':
		indep_flg = 1; paired_flg = 0;
		break;
     
          case 'p':
		paired_flg = 1;
		break;
     
	  default:
		fprintf (stderr,"var: unknown switch '%s'\n",*argv);
		exit(2);

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"var: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         }
       else {
         filnam = cptr;
         run();
       }
       if (argc <= 1) break;
     }
    }
    else run();
  }
 while (argc > 0);
}

/****************************************/

void run()

/* skip blank lines and read list of
   x,y values until either a blank line or the
   end of text */

#define STRSIZ 80
 
{
    static int  i;
    static int n,nfields;
    static char str[STRSIZ]={0};
    static int found;
    static double x,y;

    double k = 1.0;
    double data1 = 0.0, data2 = 0.0;	/* data entry */
    double diff=0.0, diffsq=0.0;
    double sum1=0.0, sum1sq=0.0;	/* sum of 1st data */
    double sum2=0.0, sum2sq=0.0;	/* sum of 2nd data */
    double tscore = 0.0;
    double var=0.0, stdev=0.0,mean1=0.0, mean2=0.0;
    double var2=0.0, stdev2=0.0;
    double stdevx=0.0, meanstd=0.0;
    double se;  /* standard error of the mean */

 n = 0;
 x = y = 0;
 for(found=0; !found; )			/* ignore leading blank lines */
  {					/* ignore lines starting with '#' */
   if (! fgets (str,STRSIZ,pictin)) return;
   if (*str == '#') continue;
   if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {found=1; break;}
  }

if (paired_flg) {
   if (found) {
      switch (nfields) {
        case 1: data1 = x; break;
        case 2: data1 = x; data2 = y; break;
      }
      diff = data2 - data1;
      diffsq = diff * diff;
      sum1 += diff;
      sum1sq += diffsq;
      n++;
  }
  for(found=0; !found; ) {
    if (! fgets (str,STRSIZ,pictin)) break;
    if (*str == '#') continue;
    if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {
       switch (nfields) {
         case 1: data1 = x; break;
         case 2: data1 = x; data2 = y; break;
       }
       diff = data2 - data1;
       diffsq = diff * diff;
       sum1 += diff;
       sum1sq += diffsq;
       if (++n >= 1000000) break;
    }
    else found = 1;
   }

   if (n == 0) n = 1;

   tscore = (sum1 / n) /
	 sqrt ( (sum1sq - (sum1*sum1 / n)) / ((n-1) * n));

   printf ("paired tscore %g n %d\n", tscore, n);

//  mean1 = sum1 / n;
//  var = (sum1sq - (sum1 * mean1)) / n;
//  if (n==1) n = 2;
//  var2 = (sum1sq - (sum1 * mean)) / (n-1);
//  stdev = sqrt (var);
//  stdev2 = sqrt (var2);
//  stdevx = stdev;
//  if (stdevx == 0) stdevx = 1;
//  meanstd = mean1 / stdevx;
//  se = stdev2 / sqrt(n);
//
//  printf ("%9.3g %9.3g %9.3g %9.3g\n",mean1,stdev,meanstd,se);

 }   /* if paired_flg */

 else if (indep_flg) {
   if (found) {
      switch (nfields) {
        case 1: data1 = x; break;
        case 2: data1 = x; data2 = y; break;
      }
      sum1 += data1;
      sum2 += data2;
      sum1sq += data1*data1;
      sum2sq += data2*data2;
      n++;
  }
  for(found=0; !found; ) {
    if (! fgets (str,STRSIZ,pictin)) break;
    if (*str == '#') continue;
    if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {
       switch (nfields) {
         case 1: data1 = x; break;
         case 2: data1 = x; data2 = y; break;
       }
      sum1 += data1;
      sum2 += data2;
      sum1sq += data1*data1;
      sum2sq += data2*data2;
      if (++n >= 1000000) break;
    }
    else found = 1;
   }

   if (n == 0) n = 1;
   mean1 = sum1 / n;
   mean2 = sum2 / n;

   tscore = (mean2 - mean1) /
	 sqrt ( ((sum2sq - (sum2*sum2 / n)) + (sum1sq - (sum1*sum1)/n)) / (n+n-2) * (1.0/n + 1.0/n) );


printf ("indep  tscore %g n %d\n", tscore, n);

}


}

