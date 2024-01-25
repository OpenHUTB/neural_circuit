 
/*    stimaddc: add up stimuli at each time to shorten stimfile.  */

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define INBUFSIZ 100

main(int argc, char **argv) 

{
	int i, row;
	double t, n1, n2, n3, n4, val, w, m, st, se;
	double ot, on1, on2, on3, on4, oval, ow, om, ost, ose;
	char a, oa;
	FILE *fd;
	static char inbuf[INBUFSIZ];

ot = on1 = on2 = on3 = on4 = oval = ow = om = ost = ose = 0; 
oa = 0;

// fprintf (stderr,"%s",argv[1]);

// fd = fopen (argv[1],"r");

printf ("# stimfile sorted, added.\n");

for (row=0; fgets(inbuf,INBUFSIZ,stdin); ) {

    switch (inbuf[0]) {

	case '#': printf ("%s",inbuf); 
		  break;

	default:
		sscanf(inbuf,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %c %lg",&t,&n1,&n2,&n3,&n4,&val,&w,&m,&st,&a,&se);

		// printf ("%-.4g %-g %-g %-g %-g %-12.6g %-g %-g %-g %c %g\n", 
		// 		  t,n1,n2,n3,n4,val,w,m,st,a,se);

	  if (				// if same time, node, action, etc
	       t != ot  || 
              n1 != on1 ||
              n2 != on2 ||
              n3 != on3 ||
              n4 != on4 ||
              w  != ow  ||
              m  != om  ||
              st != ost ||
              a  != oa  ||
              se != ose )  
	   {   				// print row, load next
		if (row++ > 0) {
		  printf ("%-.4g %-g %-g %-g %-g %-12.6g %-g %-g %-g %c %g\n", 
				  ot,on1,on2,on3,on4,oval,ow,om,ost,oa,ose);
		};
		ot  = t; 
		on1 = n1; 
		on2 = n2;
		on3 = n3;
		on4 = n4;
		oval = val;
		ow  = w;
		om  = m;
		ost = st;
		oa  = a;
		ose = se;
	    }
	    else oval += val; 
	    break;

    }   /* */
 }

    printf ("%-.4g %-g %-g %-g %-g %-12.6g %-g %-g %-g %c %g\n",        // print the last line
			  ot,on1,on2,on3,on4,oval,ow,om,ost,oa,ose);
}
