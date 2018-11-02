#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 Dir | sh

USAGE
if (@ARGV==0) { die $usage }

my $dir = shift @ARGV;
my $xmls = `ls $dir`;
my @xmls = split /\n/, $xmls;

foreach ( @xmls ) {
	chomp;
	s/\|/\\|/g;
	my $first = $_;
	s/\.xml$//;
	my $second = $_;
	print "mv $dir/$first $dir/$second\n";
#	system `mv $dir/$first $dir/$second`;
}
