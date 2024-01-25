#!/bin/sh
#
#  script to fix SVSIZ lines in nval.n files
#
#

awk ' BEGIN {first = 0}

$1=="xcone" { if (!first)  {ncols = NF; first = 1;} }
   // { if (index($0,"SVSIZ")) { 
	   printf (" ");
	   for (i=1; i<ncols; i++) {
	      printf ("     1.0");
           }
	   printf ("   %s    # synaptic vesicle size\n",$i);
         } 
         else  print $0; 
      } 
   ' $1
