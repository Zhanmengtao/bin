#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 annot_kaas.txt genelist.txt

	--fdr    set the false discovery rate, default: 0.05
    --pdn    When a kegg term was perpared for Fisher's exact test, the i(the num of this kegg terms have so many differential expression genes) should >= this value, default: 3

    annot_kaas.txt: the first column should be the gene ID; Kd+ can be one or more in one row; the same geneID can be in one or more row.

    This program is designed to perform pathway enrichment analysis using the kaas result. These files will be produced: enrichment.keggOrtholog, enrichment.keggPathway, KEGGpathway.annot, K_codes_not_in_pathway_map.list, 4value_ko.tab, 4value_K.tab and ko00001.keg.

USAGE
if (@ARGV==0){die $usage}

my ($fdr_threshold,$pathway2degs_num_threshold);
GetOptions(
        "fdr:s"=>\$fdr_threshold,
        "pdn:s"=>\$pathway2degs_num_threshold,
);
$fdr_threshold ||= 0.05;
$pathway2degs_num_threshold ||= 3;

# 下载并读取 ko00001.keg 信息
system `wget "http://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=" -O ko00001.keg` unless -e "ko00001.keg";
open KO,'<',"ko00001.keg" or die $!;
my ($A,$B,$C,$D,$path_DESC,$KO_DESC);
my %K2ko; my %K2DESC; my %ko2desc;
while (<KO>) {
	if (/^A<b>(.+)<\/b>/) {$A=$1;}
	elsif (/^B\s+<b>(.+)<\/b>/) {$B=$1;}
	elsif (/^C\s+\d+\s+(.+)\s+\[PATH:(ko\d+)\]/) {
		$path_DESC=$1;
		$C=$2;
		$ko2desc{$C} = "$A\t$B\t$path_DESC";
	}
	elsif (/^D\s+(K\d+)\s+(.*)/) {
		$D=$1;
		$KO_DESC=$2;
		$K2ko{$D}{$C} = 1;
		$K2DESC{$D} = $KO_DESC;
	}
}
close KO;

# 读取 KAAS 注释结果，得到 gene => K 的 Hash。
open K_ANNOT, '<', $ARGV[0] or die $!;
my (%K_annot,  @geneID, %geneID);
while (<K_ANNOT>) {
	if (/^(\S+)/) {
		my $geneID = $1;
		push @geneID, $geneID if ( !exists $geneID{$geneID} );
		$geneID{$geneID} = 1;
		my @K_annot = $_ =~ m/(K\d+)/g;
		foreach (@K_annot) {
			$K_annot{$geneID}{$_} = 1;
		}
	}
}
close K_ANNOT;

# 读取 DEG 信息。
open DEGLIST, '<', $ARGV[1] or die $!;
my %deg;
while (<DEGLIST>) {
	if (/^(\S+)/) { $deg{$1} = 1 }
}
close DEGLIST;
my $total_deg_num = keys %deg;

# 对 KAAS 注释结果进行整理，得到 ko => gene 和 K => gene 的 Hash。 并输出 pathway 注释结果，和不能 map 到 pahtway 的 K 编号，分别输出到文件 'KEGGpathway.annot' 和 'K_codes_not_in_pathway_map.list'。
open PATHWAY, '>', 'KEGGpathway.annot' or die $!;
open K_NOT_IN_MAP, '>', 'K_codes_not_in_pathway_map.list' or die $!;
my ($total_gene_num, $total_gene_own_K_num, $total_gene_own_ko_num, $deg_own_K_num, $deg_own_ko_num);
my (%ko2gene, %K2gene);
foreach (@geneID) {
	my $geneID = $_;
	$total_gene_num ++;
	if (exists $K_annot{$_}) {
		$total_gene_own_K_num ++; 
		if (exists $deg{$geneID}) { $deg_own_K_num ++ }
	}
	my $gene_own_ko = 0;
	my @K_annot = sort keys %{$K_annot{$_}};
	foreach (@K_annot) {
		my $K_annot = $_;
		$K2gene{$_}{$geneID} = 1;
		if (exists $K2ko{$_}) {
			$gene_own_ko = 1;
			foreach (keys %{$K2ko{$_}}) {
				$ko2gene{$_}{$geneID} = 1;
				my $ko_desc = $ko2desc{$_};
				my $k_desc = $K2DESC{$K_annot};
				print PATHWAY "$geneID\t$K_annot\t$k_desc\t$_\t$ko_desc\n";
			}
		}else {
			print K_NOT_IN_MAP "$_\n";
			print PATHWAY "$geneID\t$K_annot\n";
		}
	}
	if ($gene_own_ko == 1) {
		$total_gene_own_ko_num ++;
		if (exists $deg{$geneID}) { $deg_own_ko_num ++ }
	}
}
print STDERR "$total_gene_num genes in total, $total_deg_num differential expression genes in total.\n$total_gene_own_K_num genes own KEGG Ortholog annotation, $total_gene_own_ko_num genes own pathway annotation.\n$deg_own_K_num DEGs own KEGG Ortholog annotation, $deg_own_ko_num DEGs own pathway annotation.\n";
close PATHWAY;
close K_NOT_IN_MAP;

# 对 K 编号 进行富集分析
# 计算富集分析的输入文件，4 个 值
open K_R, '>', '4value_K.tab' or die $!;
my (%K2deg, %K2nondeg, $K_fisher_exact_test_num);
foreach my $K (sort keys %K2gene) {
	my ($x, $m, $n, $k);
	# 使用 R 的 dhyper 函数进行超几何分布计算： p = dhyper(x, m, n, k) . 其中， m 是桶里面白球的个数， n 是黑球的个数， k 是从桶中随机取出的球数， x 是取出球中白球的个数。 该 p 值表示，取出的 k 个球中白球有 x 个的概率。
	# 若 p(X>=x) 的值低于 0.05, 则表示，取出的球 >= x 的事件是小概率事件。如果实际中取出的球为 x ， 则表示发生了富集现象。
	# 因此，进行富集分析，则要计算 1 - p(X<=(x-1)) 的概率。该值低于 0.05 ， 则表示有富集。
	my $K_R_out = 0;
	foreach (keys %{$K2gene{$K}}) {
		$m ++;
		if (exists $deg{$_}) { $x ++; $K_R_out = 1; $K2deg{$K}{$_} = 1; }
		else { $K2nondeg{$K}{$_} = 1 }
	}
	$n = $total_gene_own_K_num - $m;
	$k = $deg_own_K_num;
	$x = $x - 1;         # 这一步很重要
	if ($K_R_out == 1 && $x >= ($pathway2degs_num_threshold - 1)) {
		print K_R "$K\t$x\t$m\t$n\t$k\n";
		$K_fisher_exact_test_num ++;
	}
}
print STDERR "$K_fisher_exact_test_num KEGG Ortholog terms is perpared for Fisher's exact test.\n";
close K_R;
# 进行超几何分布计算
open R_SCRIPT, '>', "phyper_adjust.R" or die $!;
print R_SCRIPT 
'phy <- read.table(file="4value_K.tab")
pvalue <- phyper(phy$V2,phy$V3,phy$V4,phy$V5,lower.tail=FALSE)
qvalue <- p.adjust(pvalue,method=\'fdr\')
value <- cbind ( pvalue, qvalue )
rownames(value)=phy$V1
value
';
close R_SCRIPT;
my @fdr_K = `cat phyper_adjust.R | R --vanilla --slave`;
`rm phyper_adjust.R`;
# 输出富集分析结果
open OUT, '>', 'enrichment.keggOrtholog' or die $!;
my $enrichment_K;
foreach (@fdr_K) {
	next unless /(\S+)\s+(\S+)\s+(\S+)/;
	my $K = $1; my $pvalue = $2; my $fdr = $3;
	if ($fdr <= $fdr_threshold) {
		my $deg_num = keys %{$K2deg{$K}};
		my $background_num = keys %{$K2gene{$K}};
		my $K_desc;
		if ($K2DESC{$K}) { $K_desc = $K2DESC{$K} } else { $K_desc = "---NA---" }
		$enrichment_K .= "$K\t$K_desc\t$deg_num\t$deg_own_K_num\t$background_num\t$total_gene_own_K_num\t$pvalue\t$fdr\t";
		foreach (keys %{$K2deg{$K}}) { $enrichment_K .= "$_\," }
		$enrichment_K =~ s/,$/\t/;
		foreach (keys %{$K2nondeg{$K}}) { $enrichment_K.= "$_\," }
		$enrichment_K =~ s/,$/\n/;
	}
}
# 对富集结果按 FDR 值进行排序
my (%for_sort, @enrichment_K);
@enrichment_K = split /\n/, $enrichment_K;
foreach (@enrichment_K) { if (/^[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t([^\t]*?)\t/) { $for_sort{$_} = $1 } }
@enrichment_K = sort { $for_sort{$a} <=> $for_sort{$b} } @enrichment_K;
foreach (@enrichment_K) { print OUT "$_\n" }
close OUT;

# 对 ko 编号 进行富集分析
open KO_R, '>', '4value_ko.tab' or die $!;
my (%ko2deg,%ko2nondeg,$ko_fisher_exact_test_num);
foreach my $ko (sort keys %ko2gene) {
	my ($x, $m, $n, $k);
	# 使用 R 的 dhyper 函数进行超几何分布计算： p = dhyper(x, m, n, k) . 其中， m 是桶里面白球的个数， n 是黑球的个数， k 是从桶中随机取出的球数， x 是取出球中白球的个数。 该 p 值表示，取出的 k 个球中白球有 x 个的概率。
	# 若 p(X>=x) 的值低于 0.05, 则表示，取出的球 >= x 的事件是小概率事件。如果实际中取出的球为 x ， 则表示发生了富集现象。
	# 因此，进行富集分析，则要计算 1 - p(X<=(x-1)) 的概率。该值低于 0.05 ， 则表示有富集。
	my $ko_R_out = 0;
	foreach (keys %{$ko2gene{$ko}}) {
		$m ++;
		if (exists $deg{$_}) { $x ++; $ko_R_out = 1; $ko2deg{$ko}{$_} = 1; }
		else { $ko2nondeg{$ko}{$_} = 1 }
	}
	$n = $total_gene_own_ko_num - $m;
	$k = $deg_own_ko_num;
	$x = $x - 1;         # 这一步很重要
	if ($ko_R_out == 1 && $x >= ($pathway2degs_num_threshold - 1)) {
		print KO_R  "$ko\t$x\t$m\t$n\t$k\n";
		$ko_fisher_exact_test_num ++;
	}
}
print STDERR "$ko_fisher_exact_test_num KEGG Pathway terms is perpared for Fisher's exact test.\n";
close KO_R;
# 进行超几何分布计算
open R_SCRIPT, '>', "phyper_adjust.R" or die $!;
print R_SCRIPT 
'phy <- read.table(file="4value_ko.tab")
pvalue <- phyper(phy$V2,phy$V3,phy$V4,phy$V5,lower.tail=FALSE)
qvalue <- p.adjust(pvalue,method=\'fdr\')
value <- cbind ( pvalue, qvalue )
rownames(value)=phy$V1
value
';
close R_SCRIPT;
my @fdr = `cat phyper_adjust.R | R --vanilla --slave`;
`rm phyper_adjust.R`;
# 输出富集分析结果
open OUT, '>', 'enrichment.keggPathway' or die $!;
my $enrichment_result;
foreach (@fdr) {
	next unless /(\S+)\s+(\S+)\s+(\S+)/;
	my $ko = $1; my $pvalue = $2; my $fdr = $3;
	if ($fdr <= $fdr_threshold) {
		my $deg_num = keys %{$ko2deg{$ko}};
		my $background_num = keys %{$ko2gene{$ko}};
		$enrichment_result .= "$ko\t$ko2desc{$ko}\t$deg_num\t$deg_own_ko_num\t$background_num\t$total_gene_own_ko_num\t$pvalue\t$fdr\t";
		foreach (keys %{$ko2deg{$ko}}) { $enrichment_result .= "$_\," }
		$enrichment_result =~ s/,$/\t/;
		foreach (keys %{$ko2nondeg{$ko}}) { $enrichment_result .= "$_\," }
		$enrichment_result =~ s/,$/\n/;
	}
}
# 对富集结果按 FDR 值进行排序
my (%for_sort, @enrichment_result);
@enrichment_result = split /\n/, $enrichment_result;
foreach (@enrichment_result) { if (/^[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t[^\t]*?\t([^\t]*?)\t/) { $for_sort{$_} = $1 } }
@enrichment_result = sort { $for_sort{$a} <=> $for_sort{$b} } @enrichment_result;
foreach (@enrichment_result) { print OUT "$_\n" }
close OUT;

=head1 Name

  kaas2pathwayEnrichmentAnalysis.pl

=head1 Description

	This program is designed to perform pathway enrichment analysis using the kaas result. These files will be produced: enrichment.keggOrtholog, enrichment.keggPathway, KEGGpathway.annot, K_codes_not_in_pathway_map.list, 4value_ko.tab, 4value_K.tab and ko00001.keg .

=head1 Version

  Author: Chen Lianfu, chenllianfu@foxmail.com
  Version: 1.1,  Date: Thu May  8 09:25:17 2014
  Version: 1.0,  Date: 2013年 07月 04日 星期四 14:39:29 CST

=head1 Usage

	--fdr         set the false discovery rate, default: 0.05
	--pdn         When a kegg term was perpared for Fisher's exact test, the i(the num of this kegg terms have so many differential expression genes) should >= this value, default: 3
	--help        output help information to screen

=head1 Example

	perl kaas2pathwayEnrichmentAnalysis.pl annot_kaas.txt genelist.txt
	annot_kaas.txt: the first column should be the gene ID; K\d+ can be one or more in one row; the same geneID can be in one or more row.

=cut
