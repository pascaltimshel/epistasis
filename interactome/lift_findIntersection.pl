#!/usr/bin/env perl

# USAGE
# perl /Users/pascaltimshel/git/epistasis/lift_findIntersection.pl supFile3_K562_interactingLoci_clusters.loci1.hg19.bed supFile3_K562_interactingLoci_clusters.loci2.hg19.bed
# paste lift_findItersection.intersection.file1 lift_findItersection.intersection.file2 > lift_findItersection.intersection.paste

### EXAMPLE INPUT file
# chr12   67280788        67280789        interaction1
# chr3    150566213       150566214       interaction2
# chr10   28973034        28973035        interaction3
# chr10   28978176        28978177        interaction4

use strict;

die "ERROR: wrong number of arguments passed to program" unless @ARGV == 2;
open(fhIN1, "<", @ARGV[0]) or die "$!";
open(fhIN2, "<", @ARGV[1]) or die "$!";

open(fhOUT1, ">", "lift_findItersection.intersection.file1") or die "$!";
open(fhOUT2, ">", "lift_findItersection.intersection.file2") or die "$!";

my %set1;
my %set2;

## Loci1
while (<fhIN1>) {
	my @fields = split();
	$set1{$fields[-1]} = 1;
}
## Loci2
while (<fhIN2>) {
	my @fields = split();
	$set2{$fields[-1]} = 1;
	if ( exists($set1{$fields[-1]}) ) {
		print fhOUT2;
	}
}
## Loci1 again
seek(fhIN1, 0, 0); #  rewind to the beginning of the file
while (<fhIN1>) {
	my @fields = split();
	if ( exists($set2{$fields[-1]}) ) {
		print fhOUT1;
	}
}

close(fhIN1);
close(fhIN2);



