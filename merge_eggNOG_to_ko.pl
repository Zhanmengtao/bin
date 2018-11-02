#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 KAAS.ko eggNOG.ko > eggNOG_KAAS_merge.ko

USAGE
if (@ARGV==0) { die $usage }

open IN, '<', $ARGV[0] or die $!;
my (%ko, %ko_KAAS, %ko_eggNOG);
while (<IN>) {
    if (m/(\S+)\s+(K\d+)/) {
        $ko{$1}{$2} = 1;
        $ko_KAAS{$1}{$2} = 1;
    }
}
close IN;

open IN, '<', $ARGV[1] or die $!;
while (<IN>) {
    my $gene_id = $1 if m/^(\S+)/;
    if (my @ko = m/(K\d+)/g) {
        foreach (@ko) {
            $ko{$gene_id}{$_} = 1;
            $ko_eggNOG{$gene_id}{$_} = 1;
        }
    }
}
close IN;

my ($num_KAAS, $num_eggNOG, $num_common, $num_total);
my ($number_KAAS, $number_eggNOG, $number_common, $number_total);
foreach my $gene_id (sort keys %ko) {
    $number_total ++;
    $number_KAAS ++ if exists $ko_KAAS{$gene_id};
    $number_eggNOG ++ if exists $ko_eggNOG{$gene_id};
    $number_common ++ if (exists $ko_KAAS{$gene_id} && exists $ko_eggNOG{$gene_id});
    foreach (sort keys %{$ko{$gene_id}}) {
        $num_total ++;
        $num_KAAS ++ if $ko_KAAS{$gene_id}{$_};
        $num_eggNOG ++ if $ko_eggNOG{$gene_id}{$_};
        $num_common ++ if ($ko_KAAS{$gene_id}{$_} && $ko_eggNOG{$gene_id}{$_});
        print "$gene_id\t$_\n";
    }
}
print STDERR "\tTotal\tKAAS\teggNOG\tCommon\n";
print STDERR "KO:\t$num_total\t$num_KAAS\t$num_eggNOG\t$num_common\n";
print STDERR "Gene:\t$number_total\t$number_KAAS\t$number_eggNOG\t$number_common\n";
