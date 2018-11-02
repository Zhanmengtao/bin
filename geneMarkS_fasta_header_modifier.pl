#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl $0 genome.fasta.faa Locus_tag > tmp; mv tmp genome.fasta.faa
	perl $0 genome.fasta.fnn Locus_tag > tmp; mv tmp genome.fasta.fnn

	if no Locus_tag given, then Locust_tag = gm

USAGE
if (@ARGV==0) {die $usage}

my $locus_tag;
if ($ARGV[1]) { $locus_tag = $ARGV[1] }
else { $locus_tag = "gm" }

open IN, $ARGV[0] or die $!;
while (<IN>) {
	if (/^>gene_(\d+)/) {
		my $code = $1;
		if (length $code == 1) { $code = $locus_tag . "_00000" . $code; }
		if (length $code == 2) { $code = $locus_tag . "_0000" . $code; }
		if (length $code == 3) { $code = $locus_tag . "_000" . $code; }
		if (length $code == 4) { $code = $locus_tag . "_00" . $code; }
		if (length $code == 5) { $code = $locus_tag . "_0" . $code; }

		print ">$code\n";
	}
	else {
		print;
	}
}
