#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 prefix.isoforms.results > prefix.genes.results

USAGE
if (@ARGV==0){die $usage}

my %hash;
open IN, $ARGV[0] or die $!;
print "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n";
<IN>;
while (<IN>) {
	@_ = split /\t/;
	my $gene_id = $_[1];
	$gene_id =~ s/_i\d+//;
	$hash{$gene_id}{$_[0]}{"length"} = $_[2];
	$hash{$gene_id}{$_[0]}{"effective_length"} = $_[3];
	$hash{$gene_id}{$_[0]}{"expected_count"} = $_[4];
	$hash{$gene_id}{$_[0]}{"TPM"} = $_[5];
	$hash{$gene_id}{$_[0]}{"FPKM"} = $_[6];
}
close IN;

foreach my $gene (sort keys %hash) {
	my ($len, $el, $ec, $tpm, $fpkm, $count, @isoforms);
	foreach (sort keys %{$hash{$gene}}) {
		$len += $hash{$gene}{$_}{"length"};
		$el += $hash{$gene}{$_}{"effective_length"};
		$ec += $hash{$gene}{$_}{"expected_count"};
		$tpm += $hash{$gene}{$_}{"TPM"};
		$fpkm += $hash{$gene}{$_}{"FPKM"};
		$count ++;
		push @isoforms, $_;
	}
	$len = int($len / $count);
	$el = (int(($el / $count) * 100 + 0.5)) / 100;
	my $isoforms = join ',', @isoforms;
	print "$gene\t$isoforms\t$len\t$el\t$ec\t$tpm\t$fpkm\n";
}
