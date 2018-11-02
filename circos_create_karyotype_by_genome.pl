#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta > karyotype.txt

USAGE
if (@ARGV==0){die $usage }

open IN, $ARGV[0] or die $!;
my ($seq_id, @seq_id, %length);
while (<IN>) {
	chomp;
	if (m/^>(\S+)/) { $seq_id = $1; push @seq_id, $1; }
	else            { $length{$seq_id} += length; }
}

my $num = 1;
foreach (@seq_id) {
	my $length = $length{$_} - 1;
	$num = $num % 22;
	print "chr\t-\t$_\t$num\t0\t$length\tchr$num\n";
	$num ++;
}
