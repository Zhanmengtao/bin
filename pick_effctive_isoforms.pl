#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 predicted_spliceGraphs.gtf cufflinks_isoforms.count_table > reads_supporting.gtf

    程序用于根据isoforms的counts数(本程序中的counts数全都表示为: cDNA序列上每1000bp的counts数)来挑选有效的isoforms信息，以用于下游的可变剪接统计。
    --min_gene_count
        单个基因的counts数需要>=该参素设定的值。若基因水平上的counts数低于此值，则表示该基因无表达，从而不需要对该基因进行可变剪接分析。若不设置该参数，则统计所有基因的counts数，将>=10的counts值按从大到小排列，取 (50%分位数/10) 作为该参数的值。
    --min_rate
        该参数的默认值为 0.1 。默认设置下，若一个基因有 n 个 counts >= 3 的 isoforms，则有效的isoform的counts数 >= gene的counts数 * 0.1 / n OR 有效的isoform的counts数 >= 10 。

USAGE
if (@ARGV==0){die $usage}

my (%opts, $min_gene_count, $min_rate);
GetOptions(
    \%opts,
    "min_gene_count:s"=>\$min_gene_count,
    "min_rate:f"=>\$min_rate
);

open IN, $ARGV[0] or die $!;
my (%genes, %gtf, %isoform_length);
while (<IN>) {
    $genes{$1}{$2} = 1 if m/gene_id\s+?\"(\S+?)\".*transcript_id\s+\"(\S+?)\"/;
    $gtf{$2} .= $_;

    @_ = split /\t/;
    $isoform_length{$2} += abs($_[4] - $_[3]) + 1;
}
close IN;

open IN, $ARGV[1] or die $!;
my %counts_isoform;
<IN>;
while (<IN>){
    @_ = split /\t/;
    $counts_isoform{$_[0]} = $_[1] * 1000 / $isoform_length{$_[0]};
}
close IN;

my %effective_isoform;
my %counts_gene;
foreach my $gene_id (keys %genes) {
    my @isoform_id = keys %{$genes{$gene_id}};
    my $counts_gene;
    foreach (@isoform_id) {
        $counts_gene += $counts_isoform{$_};
    }
    $counts_gene{$gene_id} = $counts_gene;
}

my @counts_gene = values %counts_gene;
my @counts_gene_gt10;
foreach (@counts_gene) { push @counts_gene_gt10, $_ if $_ > 10; }
@counts_gene_gt10 = sort {$b <=> $a} @counts_gene_gt10;
my $quantile = $counts_gene_gt10[@counts_gene_gt10*0.5]; 
print STDERR "Gene per 1000bp counts values' 50% quantile is: $quantile\n";
$min_gene_count ||= int($quantile / 10);
print STDERR "Min gene count was set to $min_gene_count\n";
$min_rate ||= 0.1;

my (%effective_isoform, $effective_genes_num);
foreach my $gene_id (keys %genes) {
    if ($counts_gene{$gene_id} >= $min_gene_count) {
        $effective_genes_num ++;

        my @isoform_id = keys %{$genes{$gene_id}};
        my $isoform_number;
        foreach (@isoform_id) {
            $isoform_number ++ if $counts_isoform{$_} >= 3;
        }
        foreach (@isoform_id) {
            $effective_isoform{$_} = 1 if $isoform_number > 0 && ($counts_isoform{$_} >= $counts_gene{$gene_id} * $min_rate / $isoform_number or $counts_isoform{$_} >= 10);
        }
    }
}

my @effective_isoform_num = sort keys %effective_isoform;
my $effective_isoform_num = @effective_isoform_num;
print STDERR "Effective genes number:   $effective_genes_num\n";
print STDERR "Effective isoform number: $effective_isoform_num\n";

foreach (sort keys %effective_isoform) {
    print $gtf{$_};
}
