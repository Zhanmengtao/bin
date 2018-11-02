#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage
    perl $0 input.sam output_prefix.tab > unpaired_match.readsID.list

    本程序用于从sam文件提取目标比对结果，转换成SSPACE的输入tab格式文件。意义：当双端测序，比如BAC end测序数据使用bowtie2比对到基因组序列上时，需要根据比对结果得到scaffold连接信息。需要利用双端reads都unique比对到基因组上的比对结果。同时将不能匹配或单端匹配到基因组上的reads ID输出到标准错误输出。
    本程序输入文件是 bowtie2 的比对结果文件（必须是paired数据比对结果）。（1）首先，若reads对有任何一个不是unique比对，则剔除该数据，再进行后续分析。（2）若reads对都能匹配，则根据匹配结果得到SSPACE的输入TAB格式文件（当前目录生成文件output_prefix.tab）；（3）若reads对单端匹配或不能匹配，则输出reads ID到标准输出；（4）标准错误输出是统计信息。

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
open OUT, ">", "$ARGV[1]" or die $!;
my ($multi_align_num, $paired_align_num, $non_paired_align_num);
while (<IN>) {
	next if m/^\@/;
	
	my $first = $_;
	my $second = <IN>;
	if ($first =~ m/XS:/ or $second =~ m/XS:/) { $multi_align_num ++ ; next; }

	my @first = split /\t/, $first;
	my @second = split /\t/, $second;
	my $flag1 = sprintf("%b",$first[1])+0;
	my $flag2 = sprintf("%b",$second[1])+0;

	#print  "$first[0]\t$flag1\n";
	if ($flag1 =~ m/(\d\d)\d\d$/ && $1 eq "00") {
		$paired_align_num ++;
		my $length1 = length $first[9];
		my $length2 = length $second[9];

		my ($chr1,$start1,$end1,$chr2,$start2,$end2);

		$flag1 =~ m/(\d)\d\d\d\d$/;
		my $strand1 = $1;
		if ($strand1 == 0) {
			$chr1 = $first[2];
			$start1 = $first[3];
			$end1 = $start1 + $length1 - 1;
		}
		elsif ($strand1 == 1) {
			$chr1 = $first[2];
			$end1 = $first[3];
			$start1 = $end1 + $length1 - 1;
		}

		$flag2 =~ m/(\d)\d\d\d\d$/;
		my $strand2 = $1;
		if ($strand2 == 0) {
			$chr2 = $second[2];
			$start2 = $second[3];
			$end2 = $start2 + $length2 - 1;
		}
		elsif ($strand2 == 1) {
			$chr2 = $second[2];
			$end2 = $second[3];
			$start2 = $end2 + $length2 - 1;
		}

		print OUT "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\n";
		#print OUT "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\n" if $chr1 ne $chr2;
	}
	else {
		$non_paired_align_num ++;
		print "$first[0]\n";
	}
}

my $total = $multi_align_num + $paired_align_num + $non_paired_align_num;
print STDERR "multi_align_num:\t$multi_align_num\npaired_align_num:\t$paired_align_num\nnon_paired_align_num:\t$non_paired_align_num\ntotal read pairs num:\t$total\n";
