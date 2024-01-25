#!/bin/sh
#
#  script to fix converted Neuromantic files for standard anatomy format for Neuron-C 
#
# format of neuron-C morphology file:
#
#   node  parent    dia     xbio     ybio    zbio    region   dendr;
#
#

awk ' BEGIN {	diascale = 0.3;
		xscale = 0.2;
		yscale = 0.5;
		zscale = -0.2;
	    }
      /^#/ { print } 
     !/^#/ { printf ("  %4d     %4d     %-7s   %-8.4g   %-8.4g  %-5.3g     %s      %d\n",
		$1, $2, $3*diascale, $4*xscale, $6*yscale, $5*zscale, $7, $8); 
	} 
   ' $1
