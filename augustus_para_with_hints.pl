#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 genome.fasta hints.gff augustus_species CPUs > augustus.out

    程序将文件genome.fasta和hints.gff中的序列和信息分割成单条后。对每条序列写出一个augustus命令，再进行并行化计算。最终合并结果。
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

open IN, $ARGV[1] or die $!;
my %hints;
while (<IN>) {
    m/(^\S+)/;
    $hints{$1} .= $_;
}
close IN;

mkdir "aug_para_with_hints.tmp" unless -e "aug_para_with_hints.tmp";
open COM, ">", "command.augustus.list" or die $!;
my $AUGUSTUS_CONFIG_PATH = `echo \$AUGUSTUS_CONFIG_PATH`;
chomp($AUGUSTUS_CONFIG_PATH);
foreach (sort keys %seq) {
    open OUT1, ">", "aug_para_with_hints.tmp/$_.fasta" or die $!;
    open OUT2, ">", "aug_para_with_hints.tmp/$_.hint.gff" or die $!;
    print OUT1 ">$_\n$seq{$_}\n";
    print OUT2 "$hints{$_}";
    close OUT1;
    close OUT2;

    print COM "augustus --gff3=on --species=$ARGV[2] --hintsfile=$ARGV[1] --extrinsicCfgFile=$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg --alternatives-from-evidence=true --allow_hinted_splicesites=atac aug_para_with_hints.tmp/$_.fasta > aug_para_with_hints.tmp/$_.out\n";
}
close COM;

my $cmdString = "ParaFly -c command.augustus.list -CPU $ARGV[3] &> command.augustus.log";
system ($cmdString) == 0 or die "Failed to execute: $cmdString\n$!\n";

foreach (sort keys %seq) {
    open IN, "aug_para_with_hints.tmp/$_.out" or die $!;
    print <IN>;
    close IN;
}
