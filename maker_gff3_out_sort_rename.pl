#!/usr/bin/perl
use strict;
use Getopt::Long;
my $usage = <<USAGE;
Usage:
	perl $0 [options] genome.maker.gff3

	--genePrefix    set gene prefix, default: g
	--maxAED        max AED value allowed, default: 1

USAGE
if (@ARGV==0){die $usage}

my ($genePrefix, $maxAED);
GetOptions(
	"genePrefix:s"=>\$genePrefix,
	"maxAED:s"=>\$maxAED
);
$genePrefix ||= "g";
$maxAED ||= 1;

print STDERR "GenePrefix was setted to $genePrefix, MaxAED value was setted to $maxAED\n";

open IN, $ARGV[0] or die $!;
my ($geneSym, %geneAnnot, %lines, %aed);
while (<IN>) {
	next if /^#/;
	next if /^\s/;
	next if exists $lines{$_};
	$lines{$_} = 1;
	if (/([^\t]+).*?\tgene\t(\d+)\t(\d+)\t[^\t]+\t([^\t]+)\t/) {
		$geneSym = "$1\t$2\t$3\t$4"; $geneAnnot{$geneSym} .= $_;
	}
	else {
		$aed{$geneSym} = $1 if /AED=(\S+?);/;
		$geneAnnot{$geneSym} .= $_
	}
}

my @geneSym = keys %geneAnnot;
my (%geneSymSort1, %geneSymSort2, %geneSymSort3, %geneSymSort4);
foreach ( @geneSym ) { $geneSymSort1{$_} = $1 if /(\d+)\t/ }
foreach ( @geneSym ) { $geneSymSort2{$_} = $1 if /\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort3{$_} = $1 if /.*\t(\d+)/ }
foreach ( @geneSym ) { $geneSymSort4{$_} = $1 if /(\+|\-)/ }
@geneSym = sort { $geneSymSort1{$a} <=> $geneSymSort1{$b} or $geneSymSort2{$a} <=> $geneSymSort2{$b} or $geneSymSort3{$a} <=> $geneSymSort3{$b} or $geneSymSort4{$a} cmp $geneSymSort4{$b} } @geneSym;

my ($geneNum,$scaffoldName,$scaffoldGeneNum);
foreach (@geneSym) {
	my $geneAnnot = $geneAnnot{$_};
	my $aed = $aed{$_};

	if ($aed <= $maxAED) {
		$geneNum ++;
		$geneNum = '0000'."$geneNum" if length $geneNum == 1;
		$geneNum = '000'."$geneNum" if length $geneNum == 2;
		$geneNum = '00'."$geneNum" if length $geneNum == 3;
		$geneNum = '0'."$geneNum" if length $geneNum == 4;

		my @geneAnnot = split /\n/, $geneAnnot;
		my ($transcriptIdNum, $fiveUTR_Num, $threeUTR_Num, $exonNum, $cdsNum);
		foreach (@geneAnnot) {
			if (/\tgene\t/) { s/ID=.*/ID=$genePrefix$geneNum;Name=$genePrefix$geneNum;/; print "$_\n" }
			if (/\tmRNA\t/) { $transcriptIdNum ++; $fiveUTR_Num=0; $threeUTR_Num=0; $exonNum=0; $cdsNum=0; s/ID=.*/ID=$genePrefix$geneNum.t$transcriptIdNum;Parent=$genePrefix$geneNum;Name=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
			if (/\tfive_prime_UTR\t/)  { $fiveUTR_Num ++; s/ID=.*/ID=$genePrefix$geneNum.t$transcriptIdNum.utr5p$fiveUTR_Num;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
			if (/\tthree_prime_UTR\t/) { $threeUTR_Num ++; s/ID=.*/ID=$genePrefix$geneNum\.t$transcriptIdNum.utr3p$threeUTR_Num;Parent=$genePrefix$geneNum\.t$transcriptIdNum;/; print "$_\n" }
			if (/\texon\t/) { $exonNum ++; s/ID=.*/ID=$genePrefix$geneNum.t$transcriptIdNum.exon$exonNum;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
			if (/\tCDS\t/)  { $cdsNum ++; s/ID=.*/ID=$genePrefix$geneNum.t$transcriptIdNum.cds$cdsNum;Parent=$genePrefix$geneNum.t$transcriptIdNum;/; print "$_\n" }
		}
		print "\n";
	}
}
$geneNum =~ s/^0*//;
print STDERR "Total gene number: $geneNum\n";
