#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 Min_count_per_1000bp Min_samples_number Trinity.fasta out_prefix RSEM_out/A.isoforms.results RSEM_out/B.isoforms.results RSEM_out/C.isoforms.results ...

For Example:
    perl $0 10 2 Trinity.fasta out RSEM_out/A.isoforms.results RSEM_out/B.isoforms.results RSEM_out/C.isoforms.results RSEM_out/D.isoforms.results

    程序输入Trinity denovo组装的Trinity.fasta文件，和Trinity附带脚本align_and_estimate_abundance.pl使用RSEM计算得到的多个isoforms.results文件，结果得到:
	out.gene_rawCounts.matrix          输出gene水平的raw count矩阵结果
	out.gene_TPM.not_cross.matrix      输出gene水平的TPM矩阵结果
	out.Trinity.fasta                  输出gene水平上表达量达到指定阈值的Trinity组装序列
	out.Unigene.fasta                  输出gene水平上表达量达到指定阈值的Trinity gene中表达量最大的isoform序列作为该Trinity gene的代表性序列

    示例中表示，输入了4个样品的表达量结果；程序计算有效长度平均值，将expected_count转换成每1000bp的count，得到基因水平的expected_count和count_per_1000bp，认定count_per_1000bp>=10代表基因有表达，确保基因在4个样品中有>=2个样品有表达或者4个样品总的count_per_1000bp之和>=10*4；程序最终给出这些基因的raw count矩阵文件结果、这些基因所有isoforms的Trinity组装结果、这些基因表达量最高的isoform的Trinity组装结果。

USAGE
if (@ARGV==0){die $usage}

my $min_count = shift @ARGV;
my $min_sample = shift @ARGV;
my $trinity_out = shift @ARGV;
my $prefix = shift @ARGV;

my (%seq, $seq_id);
open IN, $trinity_out or die $!;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) { $seq_id = $1; }
	else { $seq{$seq_id} .= $_; }
}
close IN;

my (%effective_length, %expected_count, %tpm, %isoform, @sample);
foreach (@ARGV) {
	open IN, $_ or die $!;
	s/.*\///;
	s/\.isoforms\.results//;
	my $sample_name = $_;
	push @sample, $sample_name;

	<IN>;
	while (<IN>) {
		@_ = split /\t/;
		$expected_count{$_[0]}{$sample_name} = $_[4];
		$effective_length{$_[0]}{$sample_name} = $_[3];
		$tpm{$_[0]}{$sample_name} = $_[5];
		$isoform{$1}{$_[0]} = 1 if $_[0] =~ m/(.*)_i\d+/;;
	}
	close IN;
}

my (%keep, %ex_count, %gene_tpm);
foreach my $gene_id (sort keys %isoform) {
	my (%count_1000bp, %tpm_gene);
	foreach my $isoform_id (keys %{$isoform{$gene_id}}) {
		my $effective_length_total;
		foreach my $sample (@sample) {
			$effective_length_total += $effective_length{$isoform_id}{$sample};
		}
		my $effective_length = $effective_length_total / @sample;

		foreach my $sample (@sample) {
			$ex_count{$gene_id}{$sample} += $expected_count{$isoform_id}{$sample};
			$gene_tpm{$gene_id}{$sample} += $tpm{$isoform_id}{$sample};
			$count_1000bp{$gene_id}{$sample} += $expected_count{$isoform_id}{$sample} * 1000 / $effective_length;
			$tpm_gene{$isoform_id} += $tpm{$isoform_id}{$sample};
		}
	}

	my @tpm_gene = sort {$tpm_gene{$b} <=> $tpm_gene{$a}} keys %tpm_gene;

	my ($number, $total);
	foreach my $isoform_id (keys %count_1000bp) {
		foreach my $sample (@sample) {
			$number ++ if $count_1000bp{$isoform_id}{$sample} >= $min_count;
			$total += $count_1000bp{$isoform_id}{$sample};
		}
	}
	if (($number >= $min_sample) or ($total >= $min_count * @sample)) {
		$keep{$gene_id} = $tpm_gene[0];
		#print "$gene_id\t$tpm_gene[0]\t$number\t$total\n";
	}
}

open OUT1, ">", "$prefix.Trinity.fasta" or die $!;
open OUT2, ">", "$prefix.Unigene.fasta" or die $!;
open OUT3, ">", "$prefix.gene_rawCounts.matrix" or die $!;
open OUT4, ">", "$prefix.gene_TPM.not_cross.matrix" or die $!;
my $sample = join "\t", @sample;
print OUT3 "\t$sample\n";
print OUT4 "\t$sample\n";
foreach my $gene_id (sort keys %keep) {
	my $isoform_id = $keep{$gene_id};
	foreach (keys %{$isoform{$gene_id}}) {
		print OUT1 ">$_\n$seq{$_}\n";
	}
	print OUT2 ">$gene_id\n$seq{$isoform_id}\n";
	my (@out1, @out2);
	foreach (@sample) {
		push @out1, $ex_count{$gene_id}{$_};
		push @out2, $gene_tpm{$gene_id}{$_};
	}
	my $out1 = join "\t", @out1;
	my $out2 = join "\t", @out2;
	print OUT3 "$gene_id\t$out1\n";
	print OUT4 "$gene_id\t$out2\n";
}
