#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 blast2go.annot interpro.tsv > go.annot

USAGE
if (@ARGV==0) { die $usage }

open IN, '<', $ARGV[0] or die $!;
my (%go, %go_blast2go, %go_interpro);
while (<IN>) {
    if (m/(\S+)\s+(GO:\d+)/) {
        $go{$1}{$2} = 1;
        $go_blast2go{$1}{$2} = 1;
    }
}
close IN;

open IN, '<', $ARGV[1] or die $!;
while (<IN>) {
    my $gene_id = $1 if m/^(\S+)/;
    if (my @go = m/(GO:\d+)/g) {
        foreach (@go) {
            $go{$gene_id}{$_} = 1;
            $go_interpro{$gene_id}{$_} = 1;
        }
    }
}
close IN;

my ($num_blast2go, $num_interpro, $num_common, $num_total);
my ($number_blast2go, $number_interpro, $number_common, $number_total);
foreach my $gene_id (sort keys %go) {
    $number_total ++;
    $number_blast2go ++ if exists $go_blast2go{$gene_id};
    $number_interpro ++ if exists $go_interpro{$gene_id};
    $number_common ++ if (exists $go_blast2go{$gene_id} && exists $go_interpro{$gene_id});
    foreach (sort keys %{$go{$gene_id}}) {
        $num_total ++;
        $num_blast2go ++ if $go_blast2go{$gene_id}{$_};
        $num_interpro ++ if $go_interpro{$gene_id}{$_};
        $num_common ++ if ($go_blast2go{$gene_id}{$_} && $go_interpro{$gene_id}{$_});
        print "$gene_id\t$_\n";
    }
}
print STDERR "\tTotal\tBlast2GO\tInterpro\tCommon\n";
print STDERR "GO:\t$num_total\t$num_blast2go\t$num_interpro\t$num_common\n";
print STDERR "Gene:\t$number_total\t$number_blast2go\t$number_interpro\t$number_common\n";
