#!/usr/bin/perl
#
#  Script to print labels from column 8 in retsim files.
#

$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

use Getopt::Long;

&GetOptions ("diascale=f"   => \$diascale,
             "xscale=f"     => \$xscale,
             "yscale=f"     => \$yscale,
             "zscale=f"     => \$zscale,
             "xyzscale=f"   => \$xyzscale,
             "xyscale=f"    => \$xyscale,
             "radincr=f"    => \$radincr
             );

while (<>) {
    chomp;	# strip record separator
    @Fld = split(' ', $_, -1);
    if (!/^#/ && $Fld[7] ne '0') { print $_; }

}
