#!/usr/bin/perl
#
#  script to fix converted Neuromantic files for standard anatomy format for Neuron-C 
#
# format of neuron-C morphology file:
#
#   node  parent    dia     xbio     ybio    zbio    region   dendr;
#
#

$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

$dia      = 0;
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

if ($xyzscale > 0) { $xyscale = $xyzscale; $zscale = $xyzscale; }
if ($xyscale > 0)  { $xscale  = $xyscale;  $yscale = $xyscale; }

while (<>) {
    chomp;	# strip record separator
    @Fld = split(' ', $_, -1);
    if (/^#/) { print $_; }
    if (!/^#/) {

       $dia = $Fld[2];
       if ($diascale != 1.0) {
           $dia *= $diascale;		# this changes dia into a number if it was a char string
       }

        printf "  %4d     %4d     %-7s    %-8.4g   %-8.4g  %-8.4g     %s      %d\n",
          $Fld[0], $Fld[1], $dia, $Fld[3] * $xscale, $Fld[4] * $yscale, $Fld[5] * $zscale, $Fld[6], $Fld[7];

    }
}

