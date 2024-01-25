#!/bin/sh
#
#  Script to convert regions in retsim files.
#   Also re-zeroes the soma location.
#

awk 'BEGIN {
     xyscale  = 1.0
     zscale   = 1.0;
     diascale = 1.0;
     par   = 0;
     dia   = 0;
     reg   = "DEND";
     dend = 0;
     first = 1;
    
     radincr = 10; 
     reg1 = radincr;
     reg2 = reg1+radincr;
     reg3 = reg2+radincr;
     reg4 = reg3+radincr;
     reg5 = reg4+radincr;
     reg6 = reg5+radincr;
     reg7 = reg6+radincr;
     reg8 = reg7+radincr;
     reg9 = 1000;
     }

     /^#/ && $2!="node" { print }

     !/^#/ { 
	    n = $1; par = $2; dia = $3; x = $4*xyscale; y = $5*xyscale; z = $6*zscale; dend = $7
	    if (first==1) {
                 printf ("#   node  parent    dia       xbio       ybio      zbio      region   dendr\n#\n");
		 first=0; xs=x; ys=y; zs=z; 
	    }
	    if (diascale != 1.0) dia *= diascale;
	    dx = x - xs;
	    dy = y - ys;
	    dz = z - zs;
	    dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if (dist > 0) {
	    if (dist < reg1) 	    reg = "R1";
	      else if (dist < reg2) reg = "R2";
	      else if (dist < reg3) reg = "R3";
	      else if (dist < reg4) reg = "R5";
	      else if (dist < reg5) reg = "R6";
	      else if (dist < reg6) reg = "R7";
	      else if (dist < reg7) reg = "R8";
	      else if (dist < reg8) reg = "R9";
	      else		    reg = "R9";
	     }
	     else		    reg = "SOMA";
	    printf ("  %4d     %4d     %-8.4g   %-8.4g   %-8.4g  %-5.3g    %-4s      %d\n",
		n, par, dia, x-xs, y-ys, z-zs, reg, dend); 
	} 
   ' $1
