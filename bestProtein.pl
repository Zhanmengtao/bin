#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 protein.fasta 1> bestProtein.fasta 2> protein_statistic.txt

USAGE
if (@ARGV==0){die $usage}

my ($proteinID, %protein, %length);
while (<>) {
	chomp;
	if (/^>(.*)/) { $proteinID = $1 }
	else          { $protein{$proteinID} .= $_; $length{$proteinID} += length }
}

my @protein = keys %protein;

my %gene;
foreach (@protein) {
	if (/(.+)\.t(\d+)/) { $gene{$1} .= "$2\t" }
}

my @gene = keys %gene;
my %geneSort;
foreach (@gene) {
	/(\d{5})$/;
	$geneSort{$_} = "$1$2";
}
@gene = sort { $geneSort{$a} <=> $geneSort{$b} } @gene;

my ($geneNum, $transcriptNum, $completeGeneNum, $noStopGenes, $noStartM, $noStopGenesNum, $noStartM_Num);
foreach (@gene) {
	$geneNum ++;
	my $geneID = $_;
	my @proteins = split /\t/, $gene{$_};
	my ($length,$gene);
	foreach (@proteins) {
		$transcriptNum ++;
		my $protein = "$geneID.t$_";
		my $length1 = $length{$protein};
#		print "$geneID.t$_\t$length1\n";
		if ($length1 > $length) { $gene = $protein; $length = $length1 }
	}
	my $sequence = $protein{$gene};
	if    ($sequence =~ s/^(M.*)\*$/$1/) { $completeGeneNum ++ }
	elsif ($sequence =~ s/^(M.*)/$1/)    { $noStopGenes .= "$_\t"; $noStopGenesNum ++ }
	elsif ($sequence =~ s/(.*)\*$/$1/)   { $noStartM .= "$_\t"; $noStartM_Num ++ }
	print ">$_\n$sequence\n";
}
print STDERR "Gene Num\t$geneNum\nComplete Gene Num\t$completeGeneNum\nTranscript Num\t$transcriptNum\n";
print STDERR ">Uncomplete Genes without stop codon\t$noStopGenesNum\n$noStopGenes\n>Uncomplete Genes without M in start condon\t$noStartM_Num\n$noStartM\n";
