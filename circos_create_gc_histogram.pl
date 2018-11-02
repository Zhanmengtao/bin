#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta 100000 > gc.histogram.txt

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my ($seq_id, @seq_id, %seq);
while (<IN>) {
	chomp;
	if (m/^>(\S+)/) { $seq_id = $1; push @seq_id, $1; }
	else            { $seq{$seq_id} .= $_; }
}
close IN;

foreach (@seq_id) {
	my $seq = $seq{$_};
	my $number = 0;
	my $length = length $seq;
	while ($length > $ARGV[1]) {
		substr($seq, 0, $ARGV[1]) =~ s/(.*)//;
		my $subseq = $1;
		my @subseq = split //, $subseq;
		my ($gc_number, $atgc_number);
		foreach (@subseq) {
			if (/[GCgc]/) {
				$gc_number ++; $atgc_number ++;
			}
			elsif (/[ATCGatcg]/) {
				$atgc_number ++;
			}
		}
		my $gc_rate = 0.5;
		if ($atgc_number != 0) {
			$gc_rate = int(($gc_number /$atgc_number) * 10000) / 10000;
		}
		my $start = $ARGV[1] * $number;
		my $end = $ARGV[1] * $number + $ARGV[1] - 1;
		print "$_\t$start\t$end\t$gc_rate\n";
		$number ++;
		$length = length $seq;
	}
	if ($seq) {
		my $subseq = $seq;
		my @subseq = split //, $subseq;
		my ($gc_number, $atgc_number);
		foreach (@subseq) {
			if (/[GCgc]/) {
				$gc_number ++; $atgc_number ++;
			}
			elsif (/[ATCGatcg]/) {
				$atgc_number ++;
			}
		}
		my $gc_rate = 0.5;
		if ($atgc_number != 0) {
			$gc_rate = int(($gc_number /$atgc_number) * 10000) / 10000;
		}
		my $start = $ARGV[1] * $number;
		my $end = $ARGV[1] * $number + $length - 1;
		print "$_\t$start\t$end\t$gc_rate\n";
	}
}
