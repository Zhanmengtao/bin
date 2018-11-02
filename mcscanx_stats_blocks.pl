#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 collinearityFile

USAGE
if (@ARGV==0){die $usage}

my (%stats, %speices1, %speices2);
open IN, $ARGV[0] or die $!;
while (<IN>) {
	next if m/^#/;
	if (/(\d+)-\s*(\d+):\s*(\S+)\s+(\S+)/) {
		$stats{$1}{$2} = 1;
		$speices1{$1}{$3} = 1;
		$speices2{$1}{$4} = 1;
	}
}

my (%blockCode2blockGeneNum, $total_num_blocks);
foreach (keys %stats) {
	my @num = keys %{$stats{$_}};
	my $num = @num;
	$blockCode2blockGeneNum{$_} = $num;
}

my (%blockGeneNum2blockCodes, %species1_5, %species1_10, %species2_5, %species2_10);
my $num_of_blocks_with_no_less_than_10_genes = 0;
foreach (keys %blockCode2blockGeneNum) {
	$total_num_blocks ++;
	my $blockGeneNum = $blockCode2blockGeneNum{$_};
	$num_of_blocks_with_no_less_than_10_genes ++ if $blockGeneNum >= 10;
	
	if ($blockGeneNum >= 10) {
		foreach (keys %{$speices1{$_}}) { $species1_10{$_} = 1; }
		foreach (keys %{$speices2{$_}}) { $species2_10{$_} = 1; }
	}
	foreach (keys %{$speices1{$_}}) { $species1_5{$_} = 1; }
	foreach (keys %{$speices2{$_}}) { $species2_5{$_} = 1; }

	$blockGeneNum2blockCodes{$blockGeneNum}{"blocksNum"} ++;
	$blockGeneNum2blockCodes{$blockGeneNum}{"blocks"} .= "$_,";
}

print "genes per block\t# of blocks\tblocks names\n";
foreach (sort {$b <=> $a} keys %blockGeneNum2blockCodes) {
	my $out = "$_\t$blockGeneNum2blockCodes{$_}{\"blocksNum\"}\t$blockGeneNum2blockCodes{$_}{\"blocks\"}";
	$out =~ s/,$/\n/;
	print $out;
}
print "\n\n$total_num_blocks blocks in all\n$num_of_blocks_with_no_less_than_10_genes blocks with 10 or more genes\n";

print "Blocks with >= 10 genes:\n";
my $num1 = keys %species1_10;
my $num2 = keys %species2_10;
print "Genes number of Speices1: $num1\n";
print "Genes number of Speices2: $num2\n";
print "Blocks with >= 5 genes:\n";
my $num1 = keys %species1_5;
my $num2 = keys %species2_5;
print "Genes number of Speices1: $num1\n";
print "Genes number of Speices2: $num2\n";
