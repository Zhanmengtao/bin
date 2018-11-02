#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 Trinity.fasta > unigene.fasta

USAGE
if (@ARGV==0) {die $usage}
my ($fastaHead, @comp, %seq, %comp, %comp_length, @bestcomp);
while (<>) {
    chomp;
    if (/^>/) {
        /^>((TRINITY_DN\d+_c\d+_g\d+)_i\d+) len=(\d+)/;
        $fastaHead = $1;
        $comp_length{$2}{$1} = $3;
        $comp{$2} = 1;
    }
    else { $seq{$fastaHead} .= $_ }
}

@comp = keys %comp;
my (%sort1, %sort2, %sort3);
foreach (@comp) {
    if (/TR(\d+)\|_c(\d+)_g(\d+)/) {
        $sort1{$_} = $1;
        $sort2{$_} = $2;
        $sort3{$_} = $3;
    }
}
@comp = sort { $sort1{$a} <=> $sort1{$b} or $sort2{$a} <=> $sort2{$b} or $sort3{$a} <=> $sort3{$b} } @comp;

foreach (@comp) {
    my %hash = %{$comp_length{$_}};
    my @hash = keys %hash;
    @hash = sort { $hash{$b} <=> $hash{$a} } @hash;
    my $bestComp = $hash[0];
    push @bestcomp, $bestComp;
}

foreach (@bestcomp) {
    print ">$_\n$seq{$_}\n"
}
