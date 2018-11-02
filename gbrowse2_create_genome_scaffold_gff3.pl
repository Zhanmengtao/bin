#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta > genome.gff

USAGE
if (@ARGV == 0) { die $usage }

my $seqID;
my @seqID;
my %seq;
while (<>) {
	chomp;
	if (/^>(.*)/) { $seqID = $1; push @seqID, $seqID }
	else          { $seq{$seqID} .= $_ }
}

foreach (@seqID) {
	my $seq = $seq{$_};
	my $length = length $seq;
	print "$_\tGENOME\tscaffold\t1\t$length\t\.\t\.\t\.\tID=$_;Name=$_\n";
}
