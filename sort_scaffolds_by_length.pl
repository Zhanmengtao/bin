#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 file.fasta > file.sort.fasta

USAGE
if (@ARGV==0) { die $usage }

open IN, $ARGV[0] or die $!;
my $sequences = join "", <IN>;
$sequences =~ s/^>//;
my @sequences = split />/, $sequences;
my %sequences;

foreach ( @sequences ) {
	my $sequence = $_;
	s/.*\n//;
	s/\n//g;
	my $length = length $_;
	$sequences{$sequence} = $length;
}

my @sorted_sequences = sort { $sequences{$b} <=> $sequences{$a} } keys %sequences;

foreach ( @sorted_sequences ) {
	print ">$_";
}
