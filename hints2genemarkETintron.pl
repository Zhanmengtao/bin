#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta bam2hints.gff > genemartET.intron.gff

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%seq, $seq_id);
while (<IN>) {
	chomp;
	if (m/^>(\S+)/) { $seq_id = $1; }
	else { $seq{$seq_id} .= $_; }
}
close IN;

my %split_site = (
	"GT\tAG" => "GT_AG\t+",
	"GC\tAG" => "GC_AG\t+",
	"AT\tAC" => "AT_AC\t+",
	"CT\tAC" => "GT_AG\t-",
	"CT\tGC" => "GC_AG\t-",
	"GT\tAT" => "AT_AC\t-",
);
open IN, $ARGV[1] or die $!;
while (<IN>) {
	if (m/mult=(\d+)/) {
		next if $1 < 4;
		@_ = split /\t/;
		my $start = substr($seq{$_[0]}, $_[3] - 1, 2);
		my $end = substr($seq{$_[0]}, $_[4] - 2, 2);
		my $split_site = "$start\t$end";
		if (exists $split_site{$split_site}) {
			my @ss = split /\t/, $split_site{$split_site};
			print "$_[0]\t$_[1]\tintron\t$_[3]\t$_[4]\t$1\t$ss[1]\t\.\t$ss[0]\n";
		}
		else {
			#print STDERR;
		}
	}
}
close IN;
