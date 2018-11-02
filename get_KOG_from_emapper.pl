#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 emapper.annotations > emapper.KOG

USAGE
    if (@ARGV==0) { die $usage }

    open IN, '<', $ARGV[0] or die $!;
    while (<IN>) {
        my $gene_id = $1 if m/^(\S+)/;
        my $descrip = $1 if m/\@euNOG\s+(.+\s(\S+))/;
        if (my $kog = m/(KOG\d+)/) {
            print "$gene_id\t$1\t$descrip\n";
        }

}
close IN;
