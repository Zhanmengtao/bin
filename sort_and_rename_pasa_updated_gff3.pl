#!/usr/bin/perl
use strict;
# 2013-11-02	chenllianfu@foxmail.com
# 2014-06-22	chenllianfu@foxmail.com
# 2014-09-11	chenllianfu@foxmail.com
my $usage = <<USAGE;
Usage:
	perl $0 pasa_updated.gff3 genePrefix > genome.gff3
USAGE
if (@ARGV == 0) { die $usage }
my $genePrefix;
if ($ARGV[1]) { $genePrefix = $ARGV[1] }
else          { $genePrefix = "g" }

open GFF, '<', $ARGV[0] or die $!;
my ($geneSym, %geneAnnot, %lines);
while ( <GFF> ) {
	next if /^#/;
	next if /^\s/;
	next if exists $lines{$_};
	$lines{$_} = 1;
	if (/([^\t]+).*?\tgene\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t/) { $geneSym = "$1\t$2\t$3\t$5\t$4"; $geneAnnot{$geneSym} .= $_;  }
	else                         { $geneAnnot{$geneSym} .= $_ }
}

my @geneSym = keys %geneAnnot;
my (%geneSymSort1, %geneSymSort2, %geneSymSort3, %geneSymSort4, %geneSymSort5);
foreach ( @geneSym ) { $geneSymSort1{$_} = $1 if /(\d+)\t/ }
foreach ( @geneSym ) { $geneSymSort2{$_} = $1 if /\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort3{$_} = $1 if /[^\t]+\t[^\t]+\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort4{$_} = $1 if /(\+|\-)/ }
foreach ( @geneSym ) { $geneSymSort5{$_} = $1 if /.*\t(\d+)/ }
@geneSym = sort { $geneSymSort1{$a} <=> $geneSymSort1{$b} or $geneSymSort2{$a} <=> $geneSymSort2{$b} or $geneSymSort3{$a} <=> $geneSymSort3{$b} or $geneSymSort4{$a} cmp $geneSymSort4{$b} or $geneSymSort5{$a} <=> $geneSymSort5{$b} } @geneSym;

my ($geneNum,$scaffoldName,$scaffoldGeneNum);
foreach (@geneSym) {
	my $geneAnnot = $geneAnnot{$_};
	$geneNum ++;
	$geneNum = '0000'."$geneNum" if length $geneNum == 1;
	$geneNum = '000'."$geneNum" if length $geneNum == 2;
	$geneNum = '00'."$geneNum" if length $geneNum == 3;
	$geneNum = '0'."$geneNum" if length $geneNum == 4;

	my @geneAnnot = split /\n/, $geneAnnot;
	my ($transcriptIdNum, $fiveUTR_Num, $threeUTR_Num, $exonNum, $cdsNum);
	foreach (@geneAnnot) {
		if (/\tgene\t/i) { s/(.*\t).*/$1ID=$genePrefix$geneNum;Name=$genePrefix$geneNum;/; print "$_\n" }
		if (/\tmRNA\t/i) { my $aed; $aed="AED=$1;" if /_AED=(\S+?);/; $transcriptIdNum ++; $fiveUTR_Num=0; $threeUTR_Num=0; $exonNum=0; $cdsNum=0; s/(.*\t).*/$1ID=$genePrefix$geneNum.t$transcriptIdNum;Parent=$genePrefix$geneNum;Name=$genePrefix$geneNum.t$transcriptIdNum;$aed/; print "$_\n" }
		if (/\tfive_prime_UTR\t/i)  { $fiveUTR_Num ++; s/(.*\t).*/$1ID=$genePrefix$geneNum.t$transcriptIdNum.utr5p$fiveUTR_Num;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
		if (/\tthree_prime_UTR\t/i) { $threeUTR_Num ++; s/(.*\t).*/$1ID=$genePrefix$geneNum\.t$transcriptIdNum.utr3p$threeUTR_Num;Parent=$genePrefix$geneNum\.t$transcriptIdNum;/; print "$_\n" }
		if (/\texon\t/i) { $exonNum ++; s/(.*\t).*/$1ID=$genePrefix$geneNum.t$transcriptIdNum.exon$exonNum;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
		if (/\tCDS\t/i)  { $cdsNum ++; s/(.*\t).*/$1ID=$genePrefix$geneNum.t$transcriptIdNum.cds$cdsNum;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
	}
		print "\n";
}
