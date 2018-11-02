#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 gatk.vcf [filter_ratio, default: 0.1]

	程序对gatk.vcf中Variants数据的各项得分(QUAL, QD, MQ, FS, MQRankSum, ReadPosRankSum)进行统计，找到相应的分位数对应的得分。以利于后续对Variants进行过滤。
	程序将6项得分分别按从小到大排序（除了FS得分是从大到小排序），并将指定的分位数值打印到屏幕上。

USAGE
if (@ARGV==0){die $usage}

my $filter_ratio = 0.1;
$filter_ratio = $ARGV[1] if $ARGV[1];
$filter_ratio = $filter_ratio * 100;

open IN, $ARGV[0] or die $!;
my (@qual, @qd, @mq, @fs, @mqranksum, @readposranksum);
while (<IN>) {
	next if m/^#/;
	@_ = split /\t/;
	push @qual, $_[5];

	my %hash = $_[7] =~ m/(\w+?)=([^;]+)/g;
	push @qd, $hash{"QD"} if exists $hash{"QD"};
	push @mq, $hash{"MQ"} if exists $hash{"MQ"};
	push @fs, $hash{"FS"} if exists $hash{"FS"};
	push @mqranksum, $hash{"MQRankSum"} if exists $hash{"MQRankSum"};
	push @readposranksum, $hash{"ReadPosRankSum"} if exists $hash{"ReadPosRankSum"};
}
close IN;

@qual = sort {$a <=> $b} @qual;
@qd = sort {$a <=> $b} @qd;
@mq = sort {$a <=> $b} @mq;
@fs = sort {$b <=> $a} @fs;
@mqranksum = sort {$a <=> $b} @mqranksum;
@readposranksum = sort {$a <=> $b} @readposranksum;

print "将各项得分按从小达到排序(除了FS得分是从大到小排序)后，其 $filter_ratio \% 处的得分分别是：\n";
print "QUAL          :\t$qual[@qual * $filter_ratio / 100]\n";
print "QD            :\t$qd[@qd * $filter_ratio / 100]\n";
print "MQ            :\t$mq[@mq * $filter_ratio / 100]\n";
print "FS            :\t$fs[@fs * $filter_ratio / 100]\n";
print "MQRankSum     :\t$mqranksum[@mqranksum * $filter_ratio / 100]\n";
print "ReadPosRankSum:\t$readposranksum[@readposranksum * $filter_ratio / 100]\n";
