#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 gene_exp.diff out_dir

USAGE
if (@ARGV==0) {die $usage}

my (%stats, %sort_by_q_value, %fpkm, %samples, %deg);
open IN, $ARGV[0] or die $!;
<IN>;
while (<IN>) {
	chomp;
	@_ = split /\t/;

	$fpkm{$_[0]}{$_[4]} = $_[7];
	$fpkm{$_[0]}{$_[5]} = $_[8];
	$samples{$_[4]} = 1;
	$samples{$_[5]} = 1;

	if ($_[13] eq "yes") {
		$stats{"$_[4]_vs_$_[5]"}{"$_[4]-up"}{$_[0]} = "$_[0]\t$_[9]\t$_[11]\t$_[12]" if $_[9] < 0;
		$stats{"$_[4]_vs_$_[5]"}{"$_[5]-up"}{$_[0]} = "$_[0]\t$_[9]\t$_[11]\t$_[12]" if $_[9] > 0;
		$sort_by_q_value{"$_[4]_vs_$_[5]"}{$_[0]} = $_[12];

		$deg{$_[0]} = 1;
	}
}
close IN;

mkdir $ARGV[1] unless -e $ARGV[1];
foreach my $vs (sort keys %stats) {
	foreach my $up (sort keys %{$stats{$vs}}) {
		open OUT, '>', "$ARGV[1]/$vs.cuffdiff_DE_results.$up.subset" or die $!;

		print OUT "ID\tlog2FC\tp_value\tq_value";
		foreach (sort keys %samples) { print OUT "\t$_";} 
		print OUT "\n";

		my $diff_num = 0;
		foreach my $id (sort {$sort_by_q_value{$vs}{$a} <=> $sort_by_q_value{$vs}{$b}} keys %{$stats{$vs}{$up}}) {
			$diff_num ++;
			print OUT "$stats{$vs}{$up}{$id}";

			foreach (sort keys %samples) { print OUT "\t$fpkm{$id}{$_}"; }
			print OUT "\n";
		}
		print "$vs.cuffdiff_DE_results.$up.subset\t$diff_num\n";
		close OUT;
	}
}

open OUT, '>', "$ARGV[1]/DEG.matrix" or die $!;
print OUT "ID";
foreach (sort keys %samples) { print OUT "\t$_";}
print OUT "\n";

my $deg_num;
foreach my $gene (sort keys %deg) {
	$deg_num ++;
	print OUT "$gene";
	foreach (sort keys %samples) { print OUT "\t$fpkm{$gene}{$_}"; }
	print OUT "\n";
}
close OUT;
print "Total DEG num: $deg_num\n";
