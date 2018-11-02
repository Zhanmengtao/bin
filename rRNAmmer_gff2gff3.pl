#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 rRNA.gff2 > rRNA.gff3

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my %rRNA_annot;
while (<IN>) {
	next if /^#/;
	next if /^\s*$/;
	my @fileds = split /\t/, $_;
	$rRNA_annot{$_}{"pos"} = $fileds[3];
	$rRNA_annot{$_}{"seq"} = $fileds[0];
}
close IN;

my @rRNA_annot = keys %rRNA_annot;
@rRNA_annot = sort { $rRNA_annot{$a}{"seq"} cmp $rRNA_annot{$b}{"seq"} or $rRNA_annot{$a}{"pos"} <=> $rRNA_annot{$b}{"pos"} } @rRNA_annot;

my $number = 0;
foreach (@rRNA_annot) {
	$number ++;
	if (length $number == 1) { $number = "00$number"; }
	if (length $number == 2) { $number = "0$number"; }
	if (length $number == 3) { $number = "$number"; }
	s/\t(\d+s_rRNA)/\tID=rRNA_${number}_$1;Name=$1;/;
	print;
}

