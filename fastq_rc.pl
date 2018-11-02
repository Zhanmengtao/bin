#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 reads.fastq > out.fastq

USAGE
if (@ARGV==0) {die $usage}

while (<>) {
	chomp;
	print "$_\n";
	$_ = <>;
	chomp;
	$_ = reverse $_;
	$_ =~ tr/ATCG/TAGC/;
	print "$_\n";
	$_ = <>;
	chomp;
	print "$_\n";
	$_ = <>;
	chomp;
	$_ = reverse $_;
	print "$_\n";
}
