#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 alternative_splicing_statisctis.txt > alternative_splicing_statisctis.stats

USAGE

open IN, @ARGV[0] or die $!;
my (%stats, %gene, %as_event);
my @as_event = ("ES", "MEE", "IR", "A5SS", "A3SS", "AFE", "ALE", "AFE1", "ALE1");
my %as_desc = ("ES", "Exon Skipping", "MEE", "Mutually Exclusive Exons", "IR", "Intron Retention", "A5SS", "Alternative 5' Splicing Sites", "A3SS", "Alternative 3' Splicing Sites", "AFE", "Alternative First Exon of both isoform", "ALE", "Alternative Last Exon of both isoform", "AFE1", "Alternative First Exon of one isoform", "ALE1", "Alternative Last Exon of one  isoform");
while (<IN>) {
	chomp;
	split;
	my $gene_id = shift @_;
	my $chr_id = shift @_;
	my $strand = shift @_;
	my $as_type = shift @_;
	my $data = join "\t", @_;
	$gene{$gene_id} = 1;

	$stats{$as_type}{"gene"}{$gene_id} = 1;
	#print "$_\n" if exists $stats{$as_type}{"position"}{"$chr_id\t$strand\t$data"};
	$stats{$as_type}{"position"}{"$chr_id\t$strand\t$data"} = 1;
	$as_event{"$chr_id\t$strand\t$data"} = 1;
}
close IN;

my $gene = keys %gene;
my $as_event = keys %as_event;
print "\n$gene genes undergo as event\n";
print "$as_event as event in total\n\n";

print "AS_event\tAS_description\tGene_number\tAS_event_number\n";
foreach (@as_event) {
	my $gene_number = keys %{$stats{$_}{"gene"}};
	my $number = keys %{$stats{$_}{"position"}};
	print "$_\t$as_desc{$_}\t$gene_number\t$number\n";
}

