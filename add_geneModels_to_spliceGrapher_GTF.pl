#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 generate_putative_sequences.out.gtf GeneModels.gtf > sample.gtf

    使用SpliceGrapher的命令generate_putative_sequences.py得到了GTF格式的文件。由于有些基因没有表达，或者其可变剪接序列过多（默认设置下超过1000，基本都是unresolved-node transcripts），则GTF文件中不包含这些基因。为了能对这些基因进行表达量和作图分析，需要使用本程序补齐GTF中的基因信息。

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my %predicted_graph;
while (<IN>) {
    if (m/gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/) {
        $predicted_graph{$1}{$2} .= $_;
    }
}
close IN;

my %detect = %predicted_graph;
my %add;
open IN, $ARGV[1] or die $!;
while (<IN>) {
    if (m/gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/ && !exists $detect{$1}) {
        $add{$1} = 1;
        $predicted_graph{$1}{$2} .= $_;
    }
}
close IN;
my $add_number = keys %add;
print STDERR "$add_number gene models were added\n";

foreach my $gene_id (sort keys %predicted_graph) {
    foreach (sort keys %{$predicted_graph{$gene_id}}) {
        print $predicted_graph{$gene_id}{$_};
    }
}
