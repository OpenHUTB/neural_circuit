#!/usr/bin/perl
#
#  Script to convert regions in retsim files.
#   Also re-zeroes the soma location.
#

$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

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
$axincr   = 10;
$dendincr = 10;
$soma_thresh = 0;

use Getopt::Long;

&GetOptions ("diascale=f"   => \$diascale,
             "xscale=f"     => \$xscale,
             "yscale=f"     => \$yscale,
             "zscale=f"     => \$zscale,
             "xyzscale=f"   => \$xyzscale,
             "xyscale=f"    => \$xyscale,
             "dendincr=f"   => \$dendincr,
             "axincr=f"     => \$axincr,
             "soma_thresh=f" => \$soma_thresh
             );

if ($xyzscale > 0) { $xyscale = $xyzscale; $zscale = $xyzscale; }
if ($xyscale > 0)  { $xscale  = $xyscale;  $yscale = $xyscale; }

$reg2 = $dendincr + $soma_thresh;
$reg1 = $reg2 + $dendincr;
# $reg3 = $reg2 + $denddincr;
$reg4 = $axincr + $soma_thresh;
$reg5 = $reg4 + $axincr;
$reg6 = $reg5 + $axincr;
$reg7 = $reg6 + $axincr;
$reg8 = $reg7 + $axincr;
$reg9 = 1000;

while (<>) {
    chomp;	# strip record separator
    @Fld = split(' ', $_, -1);
    if (/^#/ && $Fld[1] ne 'node') { print $_; }

    if (!/^#/) {
	$n = $Fld[0];
	$par = $Fld[1];
	$dia = $Fld[2];
	$X = $Fld[3] * $xscale;
	$Y = $Fld[4] * $yscale;
	$Z = $Fld[5] * $zscale;
	$oreg = $Fld[6];
	$dend = $Fld[7];
	if ($first == 1) {		# soma is first node
	    printf

	      (("#   node  parent    dia       xbio       ybio      zbio      region   dendr\n#\n"));
	    $first = 0;
	    $xs = $X;
	    $ys = $Y;
	    $zs = $Z;
	}
	if ($diascale != 1.0) {
	    $dia *= $diascale;		# this changes dia into a number if it was a char string
	}
	$dx = $X - $xs;
	$dy = $Y - $ys;
	$dz = $Z - $zs;
#	$dist = sqrt($dx * $dx + $dy * $dy + $dz * $dz);
	$dist = $dz;
#
#	Look for regions in axon (negative z dist)
#
	if ($oreg eq 'R8') { $reg = 'R8'; }
	else {
	  if ($dist > 0) {
	      if    ($dist < $reg2)  { $reg = 'R2'; }
	      elsif ($dist < $reg1)  { $reg = 'R1'; }
	      else                   { $reg = 'R1'; }
	  }
	  elsif ($dist < 0) {
	      if    (-$dist < $reg4) { $reg = 'R4'; }
	      elsif (-$dist < $reg5) { $reg = 'R5'; }
	      elsif (-$dist < $reg6) { $reg = 'R6'; }
	      elsif (-$dist < $reg7) { $reg = 'R7'; }
	      # elsif (-$dist < $reg8) { $reg = 'R8'; }
	      # elsif (-$dist < $reg9) { $reg = 'R9'; }
	      else                   { $reg = 'R7'; }
	  }
	  elsif ($oreg eq 'PSOMA') { $reg = 'PSOMA'; }
	  else                     { $reg = 'SOMA'; }
	}
        
	printf "  %4d     %4d     %-8s   %-8.4g   %-8.4g  %-5.3g    %-4s      %d\n",
	        $n, $par, $dia, $X - $xs, $Y - $ys, $Z - $zs, $reg, $dend;
    }
}

