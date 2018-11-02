#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 genome.fasta > contig.fasta 2> genome.agp

USAGE
if (@ARGV==0){die $usage}

print STDERR "##agp-version\t2.0\n";
open IN, $ARGV[0] or die $!;
my $contig_code = 1;
my ($seqID, $seq);
while (<IN>) {
	chomp;
	if (/>(\S+)/) {
		my $start = 1;
		my $part_num = 1;
		my $end;
		while ($seq =~ s/([^Nn]+)//) {
			my $length = length $1;
			$end = $length + $start - 1;
			print STDERR "$seqID\t$start\t$end\t$part_num\tW\tcontig_$contig_code\t1\t$length\t+\n";
			my $contig_seq = $1;
			$contig_seq =~ s/(\w{60})/$1\n/g;
			$contig_seq =~ s/\n$//;
			print ">contig_$contig_code\n$contig_seq\n";
			$part_num ++;
			$contig_code ++;

			if ($seq =~ s/([Nn]+)//) {
				$length = length $1;
				$start = $end + 1;
				$end = $length + $start - 1;
				print STDERR "$seqID\t$start\t$end\t$part_num\tN\t$length\tscaffold\tyes\tpaired-ends\n";
				$start = $end + 1;
				$part_num ++;
			}
		}
		
		$seqID = $1;
		$seq = "";
	}
	else {
		$seq .= $_;
	}
}

my $start = 1;
my $part_num = 1;
my $end;
while ($seq =~ s/([^Nn]+)//) {
	my $length = length $1;
	$end = $length + $start - 1;
	print STDERR "$seqID\t$start\t$end\t$part_num\tW\tcontig_$contig_code\t1\t$length\t+\n";
	my $contig_seq = $1;
	$contig_seq =~ s/(\w{60})/$1\n/g;
	$contig_seq =~ s/\n$//;
	print ">contig_$contig_code\n$contig_seq\n";
	$part_num ++;
	$contig_code ++;

	if ($seq =~ s/([Nn]+)//) {
		$length = length $1;
		$start = $end + 1;
		$end = $length + $start - 1;
		print STDERR "$seqID\t$start\t$end\t$part_num\tN\t$length\tscaffold\tyes\tpaired-ends\n";
		$start = $end + 1;
		$part_num ++;
	}
}
