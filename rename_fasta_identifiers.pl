#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 FastaFile Seq_identity[default:scaffold]

USAGE
if (@ARGV==0) {die $usage}

open FASTA, '<', $ARGV[0] or die "Cannot open the fasta file $ARGV[0]! $!";
my $fasta_identity;
if ($ARGV[1]) { $fasta_identity = $ARGV[1] }
else          { $fasta_identity = "scaffold" }

my $num = 0;
while (<FASTA>) {
	if (/^>/) { $num ++; print ">$fasta_identity"."_$num\n" }
	else      { print }
}
print STDERR "$num Sequences in total\n";
