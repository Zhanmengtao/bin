#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 num1 num2 num3 num4

    num1 : 小范围 的 目标个数
    num2 : 小范围 的 总体个数
    num3 : 大背景 的 目标个数
    num4 : 大背景 的 总体个数

结果给出富集分析的 P 值。
num1 值推荐 >= 3 才进行富集分析。

USAGE
if (@ARGV==0){die $usage}

my $m = $ARGV[2];
my $n = $ARGV[3] - $ARGV[2];
my $k = $ARGV[1];
my $x = $ARGV[0];
# 使用 R 的 dhyper 函数进行超几何分布计算： p = dhyper(x, m, n, k) . 其中， m 是桶里面白球的个数， n 是黑球的个数， k 是从桶中随机取出的球数， x 是取出球中白球的个数。 该 p 值表示，取出的 k 个球中白球有 x 个的概率。
# 因此，进行富集分析，则要计算 1 - p(X<=(x-1)) 的概率。该值低于 0.05 ， 则表示有富集。
$x --;
my $CmdString = "echo \"phyper($x, $m, $n, $k, lower.tail=FALSE)\" | Rscript -";
my $p_value = `$CmdString`;
print STDERR "CMD: $CmdString\n";
$p_value =~ s/\[1\] //;
print "Enrichment P value is : $p_value";
