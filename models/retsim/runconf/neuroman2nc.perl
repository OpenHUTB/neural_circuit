#!/usr/bin/perl
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

$par = 0;
$dia = 0;
$reg = 'DEND';
$dend = 0;
$first = 1;

$diascale = 1.0;
$xscale   = 1.0;
$yscale   = 1.0;
$zscale   = 1.0;
$xyscale  = 0;
$xyzscale = 0;


use Getopt::Long;

&GetOptions ("diascale=f"   => \$diascale,
             "xscale=f"     => \$xscale,
             "yscale=f"     => \$yscale,
             "zscale=f"     => \$zscale,
             "xyzscale=f"   => \$xyzscale,
             "xyscale=f"    => \$xyscale
             );

if ($xyzscale > 0) { $xscale *= $xyzscale; $yscale *= $xyzscale; $zscale *= $xyzscale; }
if ($xyscale > 0)  { $xscale  *= $xyscale;  $yscale *= $xyscale; }

printf ("#   node  parent    dia       xbio       ybio      zbio      region   dendr\n#\n");

while (<>) {
    @Fld = split(' ', $_, -1);
    print $_ if $awk;if (!/^#/) {
	$n = $Fld[0] - 1;
	$par = $Fld[6] - 1;
	$dia = $Fld[5] * $diascale;
	$X = $Fld[2] * $xscale;
	$Y = $Fld[3] * $yscale;
	$z = $Fld[4] * $zscale;

	if ($par < 0) { $par = $n; $reg = 'SOMA'; }
	elsif ($Fld[1] == 1) { $reg = 'SOMA'; }
	elsif ($Fld[1] == 2) { $reg = 'AXON'; }
	else { $reg = 'DEND'; }

	if ($first == 1) {
	    $first = 0;
	    $xs = $X;
	    $ys = $Y;
	    $zs = $z;
	}
	printf "  %4d     %4d     %-7.4g   %-8.4g   %-8.4g  %-5.3g     %s      %d\n",
	$n, $par, $dia, $X - $xs, $Y - $ys, $z - $zs, $reg, $dend;
    }
}
