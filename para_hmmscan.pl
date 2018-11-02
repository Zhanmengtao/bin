#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] protein.fasta

    --out <string>    default: out
    设置输出文件前缀。程序会输出3个文件: out.txt, out.tbl, out.domtbl。

    --chunk <in>    default: 100
    设置每个数据块的protein序列条数。程序将protein.fasta序列从头到尾的分割成多份，每100条相邻的序列分配到一个fasta文件中；每个fasta文件写出一条hmmscan命令。

    --cpu <int>    default: 4
    程序调用ParaFly对hmmscan命令进行并行化计算，此参数传递给ParaFly，表示并行运行数目。

    --hmmscan <string>    default: " --cpu 1 -E 1e-5 --domE 1e-5"
    传递参数给hmmscan，注意该字符串最前面要有一个空格。

    --hmm_db <string>    default: "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm"
    设置hmmscan比对的数据库路径。

    --tmp_prefix <string>    default: "hmmscan"
    设置临时文件或文件夹前缀。默认设置下，程序生成hmmscan.command.list, hmmscan.tmp/等临时文件或目录。

USAGE
if (@ARGV==0){die $usage}

my ($out, $chunk, $cpu, $hmmscan, $hmm_db, $tmp_prefix);
GetOptions(
    "out:s" => \$out,
    "chunk:i" => \$chunk,
    "cpu:i" => \$cpu,
    "hmmscan:s" => \$hmmscan,
    "hmm_db:s" => \$hmm_db,
    "tmp_prefix:s" => \$tmp_prefix,
);
$out ||= "out";
$chunk ||= 100;
$cpu ||= 4;
$hmmscan ||= " --cpu 1 -E 1e-5 --domE 1e-5";
$hmm_db ||= "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm";
$tmp_prefix ||= "hmmscan";

mkdir "$tmp_prefix.tmp" unless -e "$tmp_prefix.tmp";
open IN, $ARGV[0] or die "Can not open the input file: $ARGV[0]\n$!\n";
open CMD, ">", "$tmp_prefix.command.list"  or die "Cannot create file $tmp_prefix.command.list\n$!\n";
my ($fasta, $number, $chunk_number, @chunk);
while (<IN>) {
    if (m/^>/) {
        $number ++;
        if ($number > $chunk) {
            $chunk_number ++;
            push @chunk, "$tmp_prefix.tmp/chunk.$chunk_number";
            open OUT, ">", "$tmp_prefix.tmp/chunk.$chunk_number.fasta" or die "Can not create file $tmp_prefix.tmp/chunk.$chunk_number.fasta\n$!\n";
            print OUT $fasta;
            print CMD "hmmscan $hmmscan -o $tmp_prefix.tmp/chunk.$chunk_number.txt --tblout $tmp_prefix.tmp/chunk.$chunk_number.tbl --domtblout $tmp_prefix.tmp/chunk.$chunk_number.domtbl --pfamtblout $tmp_prefix.tmp/chunk.$chunk_number.pfamtbl $hmm_db $tmp_prefix.tmp/chunk.$chunk_number.fasta\n";
            close OUT;
            $number  = 1;
            $fasta = "";
        }
    }
    $fasta .= $_;
}
close IN;
if ($fasta) {
    $chunk_number ++;
    push @chunk, "$tmp_prefix.tmp/chunk.$chunk_number";
    open OUT, ">", "$tmp_prefix.tmp/chunk.$chunk_number.fasta" or die "Can not create file $tmp_prefix.tmp/chunk.$chunk_number.fasta\n$!\n";
    print CMD "hmmscan $hmmscan -o $tmp_prefix.tmp/chunk.$chunk_number.txt --tblout $tmp_prefix.tmp/chunk.$chunk_number.tbl --domtblout $tmp_prefix.tmp/chunk.$chunk_number.domtbl --pfamtblout $tmp_prefix.tmp/chunk.$chunk_number.pfamtbl $hmm_db $tmp_prefix.tmp/chunk.$chunk_number.fasta\n";
    print OUT $fasta;
    close OUT;
}
close CMD;

my $cmdString = "ParaFly -c $tmp_prefix.command.list -CPU $cpu &> /dev/null";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

open OUT1, ">", "$out.txt" or die "Can not create file $out.txt\n$!\n";
open OUT2, ">", "$out.tbl" or die "Can not create file $out.tbl\n$!\n";
open OUT3, ">", "$out.domtbl" or die "Can not create file $out.domtbl\n$!\n";
open OUT4, ">", "$out.pfamtbl" or die "Can not create file $out.pfamtbl\n$!\n";
foreach (@chunk) {
	open IN, "$_.txt" or die "Can not open file $_.txt\n$!\n";
	print OUT1 <IN>;
	close IN;

	open IN, "$_.tbl" or die "Can not open file $_.tbl\n$!\n";
	print OUT2 <IN>;
	close IN;

	open IN, "$_.domtbl" or die "Can not open file $_.domtbl\n$!\n";
	print OUT3 <IN>;
	close IN;

	open IN, "$_.pfamtbl" or die "Can not open file $_.pfamtbl\n$!\n";
	print OUT4 <IN>;
	close IN;
}
