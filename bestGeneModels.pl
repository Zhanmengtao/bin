#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 gff3File 1> bestGeneModels.gff3 2> geneModelsStatistic

USAGE
if (@ARGV==0) {die $usage}

my ($geneID, @geneID, %annot);
while (<>) {
	next if /^\s*$/;
	next if /^#/;
	if (/\tgene\t.*ID=([^;]+);/) { $geneID=$1; push @geneID, $geneID; $annot{$geneID} .= $_ }
	else                         { $annot{$geneID} .= $_ }
}

my ($total_geneNum, $noAS_geneNum, $AS_geneNum, %AS_geneNum, %AS_genes);
foreach (@geneID) {
	$total_geneNum ++;
	my $annot = $annot{$_};
	unless ($annot =~ m/\tmRNA\t.*\tmRNA\t/s) { 
		$noAS_geneNum ++;
		print "$annot\n"
	}else {
		$AS_geneNum ++;
		my $AS_num;
		my @annot = split /\n/, $annot;
		my $row_gene = shift @annot;
		my ($mRNAID, @mRNAID, %mRNAAnnot);
		foreach (@annot) { 
			if (/\tmRNA\t.*ID=([^;]+);/) { $mRNAID=$1; push @mRNAID,$mRNAID; $mRNAAnnot{$mRNAID} .= "$_\n" }
			else                         { $mRNAAnnot{$mRNAID} .= "$_\n" }
		}

		my ($cdsLength,$choose);
		foreach (@mRNAID) {
			$AS_num ++;
			my $mRNAAnnot = $mRNAAnnot{$_};
			my @mRNAAnnot = split /\n/, $mRNAAnnot;
			my $mRNAcdsLength;
			foreach (@mRNAAnnot) {
				if (/\tCDS\t(\d+)\t(\d+)\t/) {
					my $acdsLength = $2 - $1;
					if ($acdsLength < 0) { $acdsLength = 0 -  $acdsLength }
					$mRNAcdsLength += $acdsLength;
				}
			}
			if ($mRNAcdsLength > $cdsLength) { $choose = $_; $cdsLength = $mRNAcdsLength }
		}
		$AS_geneNum{$AS_num} ++;
		$AS_genes{$AS_num} .= "$_\t";
		print "$row_gene\n$mRNAAnnot{$choose}\n";
	}
}

print STDERR "Total Genes Num is                  : $total_geneNum\n";
print STDERR "No Alternative Spliced Genes Num is : $noAS_geneNum\n";
print STDERR "Alternative Spliced Genes Num is    : $AS_geneNum\n";

my @AS_geneNum = keys %AS_geneNum;
@AS_geneNum = sort { $a <=> $b } @AS_geneNum;

foreach (@AS_geneNum) {
	print STDERR "$_\t$AS_geneNum{$_}\t$AS_genes{$_}\n";
}
