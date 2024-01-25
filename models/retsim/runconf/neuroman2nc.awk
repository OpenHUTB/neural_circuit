#!/bin/sh
#
#  script to convert Neuromantic files to standard anatomy format for Neuron-C 
#
# x,y scale is 1024 for 260 um, or about 0.25;
# z scale is 1 um
#
# Neuromantic format: n T x y z R P
#
# n is an integer label that identifies the current point and 
#  increments by one from one line to the next.
#
# T is an integer representing the type of neuronal segment, 
#    such as soma, axon, apical dendrite, etc. 
#
#  0 = undefined
#  1 = soma
#  2 = axon
#  3 = dendrite
#  4 = apical dendrite
#  5 = fork point
#  6 = end point
#  7 = custom
#
#   x, y, z gives the cartesian coordinates of each node.
#
#   R is the radius at that node.
#
#   P indicates the parent (the integer label) of the current point or -1 
#      to indicate an origin (soma). 
#

awk 'BEGIN {
     xyscale = 1.0; 
     zscale = 1.0;
     par   = 0;
     dia   = 0;
     reg   = "DEND";
     dend = 0;
     first = 1;
     printf ("#   node  parent    dia       xbio       ybio      zbio      region   dendr\n#\n");
     }

     !/^#/ { 
	    n = $1-1; par = $7-1; dia = $6; x = $3/xyscale; y = $4/xyscale; z = $5/zscale; 
	    if (par<0) { par = n; reg = "SOMA";} 
	    else if ($2==1) {reg = "SOMA"} 
	    else if ($2==2) {reg = "AXON";} 
	    else {reg = "DEND"};
	    if (first==1) { first=0; xs=x; ys=y; zs=z; }
	    printf ("  %4s     %4s     %-7.4g   %-8.4g   %-8.4g  %-5.3g     %s      %d\n",
		n, par, dia, x-xs, y-ys, z-zs, reg, dend); 
	} 
   ' $1
