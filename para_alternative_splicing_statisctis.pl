#!/usr/bin/perl
use strict;
use File::Basename;

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf CPUs > alternative_splicing_statisctis.txt

    本程序运行需要该程序同目录下有alternative_splicing_statisctis.pl程序支持，同时需要Paraly命令能直接使用。
    本程序将gtf文件分割成一个个基因，对每个基因都使用alternative_splicing_statisctis.pl命令进行AS分析。同时调用Paraly进行并行化计算。

USAGE
if (@ARGV==0) {die $usage}

my $dirname = dirname $0;

open IN, $ARGV[0] or die $!;
mkdir "para_alternative_splicing_statisctis" unless -e "para_alternative_splicing_statisctis";
my %gtf;
while (<IN>) {
    $gtf{$1} .= $_ if m/gene_id\s+?\"(\S+?)\"/;
}
close IN;

open COM, ">", "command.alternative_splicing_statisctis.list" or die $!;
foreach (keys %gtf) {
    open OUT, ">", "para_alternative_splicing_statisctis/$_.gtf" or die $!;
    print OUT $gtf{$_};
    close OUT;
    print COM "$dirname/alternative_splicing_statisctis.pl para_alternative_splicing_statisctis/$_.gtf > para_alternative_splicing_statisctis/$_.out\n";
}
close COM;

my $cmdString = "ParaFly -c command.alternative_splicing_statisctis.list -CPU $ARGV[1] > /dev/null";
(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

foreach (sort keys %gtf) {
    open IN, "para_alternative_splicing_statisctis/$_.out" or die $!;
    print <IN>;
    close IN;
}
