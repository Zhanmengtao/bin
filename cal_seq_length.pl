#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl cal_seq_length.pl fastaFile

USAGE
if (@ARGV==0) {die $usage}

my ($geneID, @geneID, %length);
while (<>) {
	chomp;
	if (/^>(.*)/) { $geneID = $1; push @geneID, $geneID }
	else          { $length{$geneID} += length }
}

foreach (@geneID) {
	print "$_\t$length{$_}\n";
}
