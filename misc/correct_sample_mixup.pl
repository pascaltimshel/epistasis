#!/usr/bin/env perl

use strict;

my %mixup_pcodes_hash = (
				"27998" => "76535",
				"76535" => "27998",
				"1083" => "24644",
				"24644" => "1083",
				"31900" => "41164",
				"41164" => "31900",
				"85682" => "2964",
				"2964" => "85682",
				"55067" => "28579",
				"28579" => "55067"
				);

my $outfile = @ARGV > 0 ? $ARGV[0] . "_no_mixup" : die "ERROR: no arguments passed to program";
#my $outfile = $ARGV[0] . "_no_mixup";

open(OUT, '>', $outfile) or die "$!";
while (<>) {
	my @fields = split();
	if ( exists($mixup_pcodes_hash{$fields[0]}) ) {
		my $n_subs = s/$fields[0]/$mixup_pcodes_hash{$fields[0]}/g;
		print "Line $. substitution (n = $n_subs): $fields[0] ==> $mixup_pcodes_hash{$fields[0]}\n";
	}
	print OUT;
}
close(OUT);

