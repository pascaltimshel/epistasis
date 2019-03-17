#!/usr/bin/env perl

### USAGE
#perl /cvar/jhlab/timshel/egcut/GTypes_hapmap2_expr/infer_chromosome_for_gen_file.pl <inputfile> <outputfile>

### Remember: this script does not accept gzipped files

use strict;

my $sttime = time;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

die "ERROR: wrong number of arguments passed to program" unless @ARGV == 2;

my $chromosome = 1;
my $position_old = 0;

### EXAMPLE .gen/inpute input file
# --- rs2905035 765522 A G 0.0200 0.2500 0.7300 0 0 1 0.0400 0.3600 0.5900 0 0.1800 0.8200 0 0.2600 0.7400 0.0800 0.4400 0.4700 0 0.1300 0.8700 0 0.1800 0.8100 0 0.1800 0
# 1 rs12124819 766409 A G 1 0 0 0 0 1 1 0 0 0 1 0 0 1 0 1 0 0 0 1 0 1 0 0 0.3300 0.3300 0.3300 0 1 0
# --- rs2980319 766985 A T 0.0200 0.2500 0.7300 0 0 1 0.0400 0.3600 0.5900 0 0.1800 0.8200 0 0.2600 0.7400 0.0800 0.4400 0.4700 0 0.1300 0.8700 0 0.1800 0.8100 0 0.1800 0

my $n_with_chr;
my $n_without_chr;
my $n_total;

open(IN, '<', $infile) or die "$!";
open(OUT, '>', $outfile) or die "$!";
while (defined(my $line = <IN>)) {
	my ($chromosome_new, $snpID, $position_new) = split(' ', $line);
	if ( $position_new < $position_old ) { # the position is now smaller than what we have observed before. That means we are now at a new chromosome.
		# note that this also works in the initial case because position_old is 0
		$chromosome++;
		print "[#Line $.; Pos_new=$position_new; Pos_old=$position_old]: New chromosome detected: $chromosome\n"
	}
	if ( $chromosome_new == '---' ) {
		$line =~ s/---/$chromosome/;
		$n_without_chr++;
	} else { # the chromosome already have a chromosome annotation. Let's make sure that this chromosome is the same as we think it is!
		$n_with_chr++;
		unless ( $chromosome_new == $chromosome ) {
			print "chromosome_new ($chromosome_new) is not equal to chromosome ($chromosome):\t $line";
		}
	}
	$n_total++;

	print OUT $line;

	$position_old = $position_new; # updateing position
}
close(IN);
close(OUT);

print "n_with_chr: $n_with_chr (" . ($n_with_chr/$n_total*100) . "%)\n";
print "n_without_chr: $n_without_chr (" . ($n_without_chr/$n_total*100) . "%)\n";
print "n_total: $n_total\n";



my $entime = time;
my $elapse = $entime - $sttime;
print "Elapsed time: ", $elapse/60, " min\n";


# while (<IN>) {
# 	chomp;
# 	#DO SOMETHING with $_
# }
### EQUIVALENT OF 
# while (defined(my $line=<IN>) ) {
# 	chomp $line;
# 	#DO SOMETHING!
# }

## Basic 
#http://docstore.mik.ua/orelly/perl3/lperl/ch06_01.htm#FOOTNOTE-142
## Advances
#http://docstore.mik.ua/orelly/perl3/prog/ch02_11.htm
