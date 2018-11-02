#!/usr/bin/perl
use strict;
use Statistics::Descriptive;

my $usage = <<USAGE;
Usage:
    perl $0 file.sam min_insertSize max_insertSize

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (@insertSize, $read_name);
while (<IN>){
    @_ = split /\t/, $_;
    $_ = abs($_[8]);
    push @insertSize, $_ if ($_ >= $ARGV[1] && $_ <= $ARGV[2] && $read_name ne $_[0]);
    $read_name = $_[0];
}
close IN;

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@insertSize);
my $mean = $stat->mean();
my $SD = $stat->standard_deviation();
print "$mean\t$SD\n";
