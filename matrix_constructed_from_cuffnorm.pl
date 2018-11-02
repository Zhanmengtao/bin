#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 cuffnorm1/genes.count_table cuffnorm2/genes.count_table ... > genes.count_table.matrix

	当需要使用cufflinks软件对不同的样品进行表达量计算，而同时使用其它软件进行差异表达分析的时候。可以按如下步骤进行分析：
	1. 使用cuffquant和cuffnorm分别对每个样品进行表达量分析。由于cuffnorm输入的数据至少有2个，因此输入2个相同的cuffquant结果即可得到该样品的raw count。
	2. 使用本程序根据多个样品的cuffnorm结果，得到count或fpkm的矩阵文件。
	3. 使用count或fpkm的矩阵文件，利用其它差异表达分需软件进行差异表达基因分析。

USAGE
if (@ARGV==0){die $usage}

my (%matrix, @sample);
foreach (@ARGV) {
	open IN, $_ or die $!;
	$_ = <IN>;
	my $sample = $1 if m/tracking_id\t(.*?)_0\t/;
	push @sample, $sample;

	while (<IN>) {
		split;
		$matrix{$_[0]}{$sample} = $_[1];
	}
	close IN;
}

my $sample = join "\t", @sample;
print "tracking_id\t$sample\n";

foreach my $gene_id (sort keys %matrix) {
	my @out;
	foreach (@sample) {
		push @out, $matrix{$gene_id}{$_},
	}
	my $out = join "\t", @out;
	print "$gene_id\t$out\n";
}
