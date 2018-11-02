#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta augustus_species CPUs > augustus.gff3

    程序将文件genome.fasta中的序列分割成单条后。对每条序列写出一个augustus命令，再进行并行化计算。最终合并结果。
    由于结果文件是合并了各条序列的比对结果，基因会有重名的问题。因此，需要执行如下命令来解决该问题：
    join_aug_pred.pl <augustus.out > augustus.gff3

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%seq, $id);
while (<IN>) {
    chomp;
    if (/>(\S+)/) { $id = $1; }
    else { $seq{$id} .= $_; }
}
close IN;

mkdir "aug_para.tmp" unless -e "aug_para.tmp";
open COM, ">", "command.augustus.list" or die $!;
foreach (sort keys %seq) {
    open OUT, ">", "aug_para.tmp/$_.fasta" or die $!;
    print OUT ">$_\n$seq{$_}\n";
    print COM "augustus --gff3=on --species=$ARGV[1] aug_para.tmp/$_.fasta > aug_para.tmp/$_.out\n";
    close OUT;
}
close COM;

my $cmdString = "ParaFly -c command.augustus.list -CPU $ARGV[2] &> command.augustus.log";
system ($cmdString) == 0 or die "Failed to execute: $cmdString\n$!\n";

foreach (sort keys %seq) {
    open IN, "aug_para.tmp/$_.out" or die $!;
    print <IN>;
    close IN;
}
