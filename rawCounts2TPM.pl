#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 raw_counts.matrix genome.gtf > TPM.matrix

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
$_ = <IN>;
print;
s/^\s+//;
chomp;
my @sample = split /\t/, $_;

my (%raw_count, @gene_id);
while (<IN>) {
	s/([^\t]+)\t//;
	my $gene_id = $1;
	push @gene_id, $gene_id;
	chomp;
	@_ = split /\t/, $_;
	foreach (@sample) {
		my $count = shift @_;
		$raw_count{$_}{$gene_id} = $count;
	}
}
close IN;

open IN, $ARGV[1] or die $!;
my %length;
while (<IN>) {
    if (m/\texon\t/) {
        @_ = split /\t/;
        my ($gene_id, $transcript_id) = ($1, $2) if $_[8] =~ m/gene_id\s+\"(.*?)\".*transcript_id\s+\"(.*?)\"/;
        $length{$gene_id}{$transcript_id} += abs($_[4] - $_[3]);
    }
}
close IN;

my %cDNA_length;
foreach my $gene_id (sort keys %length) {
    my @transcript = sort {$length{$gene_id}{$b} <=> $length{$gene_id}{$a}} keys %{$length{$gene_id}};
	$cDNA_length{$gene_id} = $length{$gene_id}{$transcript[0]};
	#print "$gene_id\t$length{$gene_id}{$transcript[0]}\n";
}

my (%norm_by_cDNA_length, %sum);
foreach my $gene_id (@gene_id) {
	foreach (@sample) {
		my $value = $raw_count{$_}{$gene_id} * 1000 / $cDNA_length{$gene_id};
		$norm_by_cDNA_length{$_}{$gene_id} = $value;
		$sum{$_} += $value;
		#print "$gene_id\t$_\t$raw_count{$_}{$gene_id}\t$value\n";
	}
}

my %coefficient;
foreach (@sample) {
	my $coefficient = 1000000 / $sum{$_};
	#print "$_\t$sum{$_}\t$coefficient\n";
	$coefficient{$_} = 1000000 / $sum{$_};
}

foreach my $gene_id (@gene_id) {
	my $out .= "$gene_id\t";
	my @out;
	foreach (@sample) {
		my $value = $norm_by_cDNA_length{$_}{$gene_id} * $coefficient{$_};
		$value =~ s/(\.\d\d\d\d).*/$1/;
		push @out, $value;
	}
	$out .= join "\t", @out;
	print "$out\n";
}
