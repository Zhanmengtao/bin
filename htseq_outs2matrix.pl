#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 L1,L2,...,LN htseq_out1.txt htseq_out2.txt ... htseq_outn.txt

USAGE
if (@ARGV==0){die $usage}

my $labels = shift @ARGV;
my @labels = split /,/, $labels;
my (%count, @label);
foreach (@ARGV) {
    my $label = shift @labels;
    push @label, $label;
    warn "Warning: lacking the lable of $_\n" unless $label;
    open IN, $_ or die $!;
    while (<IN>) {
        last if m/^__/;
        chomp;
        @_ = split /\t/;
        $count{$_[0]}{$label} = $_[1];
    }
    close IN;
}

my $out = "\t" . (join "\t", @label);
print "$out\n";
foreach my $gene_id (sort keys %count) {
    $out = "$gene_id\t";
    foreach (@label) {
        $out .= "$count{$gene_id}{$_}\t";
    }
    $out =~ s/\t$/\n/;
    print $out;
}
