#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 Trinity.fasta isoforms.counts.matrix min_total_counts_number > unigene.fasta

For Example:
    perl $0 Trinity.fasta isoforms.counts.matrix 50 > unigene.most_abundant.fasta
    perl $0 Trinity.fasta isoforms.TMM.fpkm.matrix 10 > unigene.most_abundant.fasta

USAGE
if (@ARGV==0) {die $usage}

open IN, $ARGV[0] or die $!;
my ($seq_id, %seq);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

open IN, $ARGV[1] or die $!;
my (%abundance, %gene_count);
while (<IN>) {
    chomp;
    if (s/^((TRINITY_DN\d+_c\d+_g\d+)_i\d+)\s+//) {
        my ($isoform_id, $gene_id) = ($1, $2);
        
        @_ = split /\s+/;
        my $count_number;
        foreach (@_) { $count_number += $_; }

        $abundance{$gene_id}{$isoform_id} = $count_number;
        $gene_count{$gene_id} += $count_number;
    }
}
close IN;

my @gene_id = sort { $gene_count{$b} <=> $gene_count{$a} } keys %gene_count;
foreach my $gene_id (@gene_id) {
    if ($gene_count{$gene_id} >= $ARGV[2]) {
        my @isoform_id = sort { $abundance{$gene_id}{$b} <=> $abundance{$gene_id}{$a} } keys %{$abundance{$gene_id}};
        print ">$isoform_id[0]    total_count_or_FPKM_of_gene:$gene_count{$gene_id}\n";

        my $seq = $seq{$isoform_id[0]};
        $seq =~ s/(\w{60})/$1\n/g;
        $seq =~ s/\n*$/\n/;
        print $seq;
    }
}
