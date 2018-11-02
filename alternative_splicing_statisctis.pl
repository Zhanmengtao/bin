#!/usr/bin/perl
use strict;
use Algorithm::Combinatorics qw(combinations permutations);

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf > alternative_splicing_statisctis.txt

USAGE
if (@ARGV==0){die $usage}

# 对每个基因进行AS分析，统计的AS类型有：
# 可变剪接类型                        缩写      对应的子程序
# 1. exon skipping                    ES        get_es
# 2. mutually exclusice exons         MEE       get_es
# 3. intron retention                 IR        get_ir
# 4. alternative 5' splicing sites    A5SS      get_a5ss,get_a3ss
# 5. alternative 3' splicing sites    A3SS      get_a5ss,get_a3ss
# 6. alternative first exon           AFE       get_afe,get_ale
# 7. alternatie last exon             ALE       get_afe,get_ale
#
# 本程序中AS分析的原理：
# 首先统计各个基因的可变剪接的exon信息，然后对可变剪接之间进行两两比较，从而得到可变剪接类型。
# 调用子程序进行两两比较。不同的可变剪接类型使用不同的子程序。以下是对这些子程序的介绍。
# 1. get_es
#    该子程序可以鉴定ES和MEE两种类型的可变剪接。
#    ES(exon skipping)                   MEE(mutually exclusice exons)
#    |||||-------------|||||             |||||----|||||-------------|||||
#    |||||----|||||----|||||             |||||-------------|||||----|||||
#
#    分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序，找2个可变剪接exon结束位置（方向为从小到大的方向，不是exon的正负义链方向）一致的点，和距离该点最近的exon起始位位置一致的点，作为目标区域。（2）然后，对目标区域进行exon数目分析。若2个可变剪接的目标区域中都有>=1个exons，且这些exons的边界位置都不一致，则为MEE；若仅有1个可变剪接的目标区域中有>=1个exons，则为ES。
#    每个ES信息由3个数据组成，这3个数据用制表符分割：共有的exon边界位点、Skipped exons首尾位置信息（多个exons信息使用分号分割）、共有的exon边界位点。这3个数据按数值从小到大排序。
#    每个MEE信息由4个数据组成，这4个数据用制表符分割：共有的exon边界位点、Skipped exons首尾位置信息（多个exons信息使用分号分割）、另一个可变剪接的Skipped exons首尾位置信息（多个exons信息使用分号分割））、共有的exon边界位点。这4个数据按数值从小到大排序。
#
# 2. get_ir
#    该子程序可以鉴定IR类型的可变剪接。
#    IR(intron retention)
#    |||||----|||||
#    ||||||||||||||
#
#    分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序。（2）寻找第一个可变剪接的intron信息，在第二个可变剪接的exon信息中寻找是否存在能覆盖该intron的exon，从而鉴定在第一个可变剪接上的IR。（3）同理，寻找第二个可变剪接的intron信息，在第一个可变剪接的exon信息中寻找是否存在能覆盖该intron的exon，从而鉴定在第二个可变剪接上的IR。
#    每个IR信息由1个数据组成，即IR的intron位置信息。
#
# 3. get_a5ss和get_a3ss
#    这两个子程序共同用于鉴定A5SS和A3SS类型的可变剪接。
#    A5SS(alternative 5' splicing sites)    A3SS(alternative 3' splicing sites)
#    |||||--------|||||                     |||||--------|||||
#    ||||||||-----|||||                     |||||-----||||||||
#    get_a5ss的分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序。（2）得到第一个可变剪接的相邻的两个exon信息，对其使用第二个可变剪接的相邻exon信息进行分析，若它们的第二个exon具有相同的起始位点（从小到大的位置排序方向，而不是正负链方向），且第二个可变剪接的第一个exon的结束位点 > 第一个可变剪接的第一个exon的结束位点，且第二个可变剪接的第一个exon的起始位点 < 第一个可变剪接的第一个exon的结束位点，则认为该位点属于A5SS和A3SS类型的可变剪接。（3）同理，得到第二个可变剪接的相邻exon信息，对其使用第一个可变剪接的相邻exon信息进行分析... （4）get_a5ss的鉴定结果若是在正义链上的基因，则是A5SS类型的可变剪接，若是在负义链上基因，则为A3SS的可变剪接。
#    get_a3ss的分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序。（2）得到第一个可变剪接的相邻的两个exon信息，对其使用第二个可变剪接的相邻exon信息进行分析，若它们的第一个exon具有相同的结束位点（从小到大的位置排序方向，而不是正负链方向），且第二个可变剪接的第二个exon的起始位点 < 第一个可变剪接的第二个exon的起始位点，且第二个可变剪接的第二个exon的结束位点 > 第一个可变剪接的第二个exon的起始位点，则认为该位点属于A5SS和A3SS类型的可变剪接。（3）同理，得到第二个可变剪接的相邻exon信息，对其使用第一个可变剪接的相邻exon信息进行分析... （4）get_a3ss的鉴定结果若是在正义链上的基因，则是A3SS类型的可变剪接，若是在负义链上基因，则为A5SS的可变剪接。
#    每个A5SS信息由3个数据组成，这3个数据用制表符分割：共同的第二个exon起始位点、一个可变剪接的第一个exon的结束位点、另一个可变剪接的第一个exon的结束位点。后2个数据按数值从小到大排序。
#    每个A3SS信息由3个数据组成，这3个数据用制表符分割：共同的第一个exon结束位点、一个可变剪接的第二个exon的起始位点、另一个可变剪接的第二个exon的起始位点。后2个数据按数值从小到大排序。
#
# 4. get_afe和get_ale
#    这两个子程序共同用于鉴定AFE和ALE类型的可变剪接。
#    AFE(alternative first exon)        ALE(alternative last exon)
#    |||||-------------|||||            |||||-------------|||||
#             |||||----|||||            |||||----|||||
#
#    get_afe的分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序。（2）分别得到这两个可变剪接的第一个和第二个exon信息，进行比较从而得到是否存在AFE或ALE。若具有相同的第二个exon起始位点，且第一个exon的起始位点与结束位点都不一致，则认为存在AFE或ALE类型的可变剪接。（3）若基因位于正义链上，则得到AFE类型的可变剪接；若基因位于负义链，则得到ALE类型的可变剪接。
#    get_ale的分析原理：（1）首先，将两个可变剪接的exon信息都按数值从小到大排序。（2）分别得到这两个可变剪接的倒数第一个和倒数第二个exon信息，进行比较从而得到是否存在AFE或ALE。若具有相同的倒数第二个exon起始位点，且倒数第一个exon的起始位点与结束位点都不一致，则认为存在AFE或ALE类型的可变剪接。（3）若基因位于正义链上，则得到ALE类型的可变剪接；若基因位于负义链，则得到AFE类型的可变剪接。
#    每个AFE信息由3个数据组成，这3个数据用制表符分割：共同的第二个exon起始位点、一个可变剪接的第一个exon的位置信息、另一个可变剪接的第一个exon的位置信息。后2个数据按数值从小到大排序。
#    每个ALE信息由3个数据组成，这3个数据用制表符分割：共同的倒数第二个exon的结束位点、一个可变剪接的倒数第一个exon的位置信息、另一个可变剪接的倒数第一个exon的位置信息。后2个数据按数值从小到大排序。

# 5. get_afe1和get_ale1
#    这两个子程序共同用于鉴定AFE1和ALE1类型的可变剪接。
#    AFE1                              ALE1
#    |||||----|||||----|||||           |||||----|||||----|||||
#           |||||||----|||||           |||||----|||||||
#    AFE1和AFE类型的可变剪接类似。不一样在于：AFE1中只有一个可变剪接对另外一个可变剪接是Alternative First exon的，而AFE则是两个可变剪接是相互都是Alternative First exon的。
#    ALE1和ALE类型的可变剪接类似。不一样在于：ALE1中只有一个可变剪接对另外一个可变剪接是Alternative Last exon的，而ALE则是两个可变剪接是相互都是Alternative Last exon的。
#    get_afe1的分析原理：（1）首先，将两个可变剪接（分别用A和B表示）的exon信息都按数值从小到大排序（不设定基因的方向，起始位点永远是小于结束位点）。（2）得到可变剪接A的前两个exon信息，在另外一个可变剪接B的exons信息中寻找：先找和A的第二个exon起始位点相同的exon、B的上一个exon的结束位点和A的第一个exon结束位点一致、B的上一个exon的起始位点和A的第一个exon起始位点不一致，则认为存在AFE1或ALE1类型的可变剪接。（3）若基因位于正义链上，则得到AFE1类型的可变剪接；若基因位于负义链，则得到ALE1类型的可变剪接。
#    get_ale1的分析原理：（1）首先，将两个可变剪接（分别用A和B表示）的exon信息都按数值从大到小排序（不设定基因的方向，起始位点永远是小于结束位点）。（2）得到可变剪接A的前两个exon信息，在另外一个可变剪接B的exons信息中寻找：先找和A的第二个exon结束位点相同的exon、B的上一个exon的起始位点和A的第一个exon起始位点一致、B的上一个exon的结束位点和A的第一个exon结束位点不一致，则认为存在AFE1或ALE1类型的可变剪接。（3）若基因位于正义链上，则得到ALE1类型的可变剪接；若基因位于负义链，则得到AFE1类型的可变剪接。
#    每个AFE1信息由3个数据组成，这3个数据用制表符分割：AFE1可变剪接的第二个exon的起始位点、AFE1可变剪接的第一个exon的位置信息、另一个可变剪接中和AFE1可变剪接的第一个exon具有相同结束位点的exon信息。
#    每个ALE1信息由3个数据组成，这3个数据用制表符分割：ALE1可变剪接的倒数第二个exon的结束位点、ALE1可变剪接的最后一个exon的位置信息、另一个可变剪接中和ALE1可变剪接的最后一个exon具有相同起始位点的exon信息。

open IN, $ARGV[0] or die $!;
my (%gene, %gene_loci);
while (<IN>) {
    if (m/\texon\t.*gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/) {
        @_ = split /\t/;
        $gene{$1}{$2} .= "$_[3]\t$_[4]\t$_[6]\n";
        $gene_loci{$1}{"scaffold"} = $_[0];
        $gene_loci{$1}{"strand"} = $_[6];
    }
}
close IN;

my %statistics;
foreach my $gene_id (sort keys %gene) {
    #print "$gene_id\n";
    my @isoform = keys %{$gene{$gene_id}};

    if (@isoform > 100) {
        my ($num, $last_as_num);
        while ($num < 5000) {
            my @two_exon_info;
            push @two_exon_info, $gene{$gene_id}{$isoform[rand(@isoform)]};
            push @two_exon_info, $gene{$gene_id}{$isoform[rand(@isoform)]};

            my $as_num;
            foreach (keys %{$statistics{$gene_id}}) {
                $as_num += keys %{$statistics{$gene_id}{$_}};
            }
            if ($as_num > $last_as_num) {
                $num = 0;
            }
            else {
                $num ++;
            }
            $last_as_num = $as_num;

            my %es = &get_es(@two_exon_info);
            foreach (keys %es) {
                $statistics{$gene_id}{"es"}{$_} = 1;
            }

            my %ir = &get_ir(@two_exon_info);
            foreach (keys %ir) {
                $statistics{$gene_id}{"ir"}{$_} = 1;
            }

            my %a5ss = &get_a5ss(@two_exon_info);
            foreach (keys %a5ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $a53ss_info = "A5SS\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $a53ss_info = "A3SS\t$_"; }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }

            my %a3ss = &get_a3ss(@two_exon_info);
            foreach (keys %a3ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $a53ss_info = "A3SS\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $a53ss_info = "A5SS\t$_"; }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }

            my %afe = &get_afe(@two_exon_info);
            foreach (keys %afe) {
                my $afle_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle_info = "AFE\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle_info = "ALE\t$_"; }
                $statistics{$gene_id}{"afle"}{$afle_info} = 1;
            }

            my %ale = &get_ale(@two_exon_info);
            foreach (keys %ale) {
                my $afle_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle_info = "ALE\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle_info = "AFE\t$_"; }
                $statistics{$gene_id}{"afle"}{$afle_info} = 1;
            }

            my %afe1 = &get_afe1(@two_exon_info);
            foreach (keys %afe1) {
                my $afle1_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle1_info = "AFE1\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle1_info = "ALE1\t$_"; }
                $statistics{$gene_id}{"afle1"}{$afle1_info} = 1;
            }

            my %ale1 = &get_ale1(@two_exon_info);
            foreach (keys %ale1) {
                my $afle1_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle1_info = "ALE1\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle1_info = "AFE1\t$_"; }
                $statistics{$gene_id}{"afle1"}{$afle1_info} = 1;
            }
        }
    }
    elsif (@isoform >= 2) {
        my $iter = combinations(\@isoform, 2);
        while (my $c = $iter->next) {
            my @two_exon_info;
            foreach (@$c) {
                push @two_exon_info, $gene{$gene_id}{$_};
            }

            my %es = &get_es(@two_exon_info);
            foreach (keys %es) {
                $statistics{$gene_id}{"es"}{$_} = 1;
            }

            my %ir = &get_ir(@two_exon_info);
            foreach (keys %ir) {
                $statistics{$gene_id}{"ir"}{$_} = 1;
            }

            my %a5ss = &get_a5ss(@two_exon_info);
            foreach (keys %a5ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $a53ss_info = "A5SS\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $a53ss_info = "A3SS\t$_"; }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }

            my %a3ss = &get_a3ss(@two_exon_info);
            foreach (keys %a3ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $a53ss_info = "A3SS\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $a53ss_info = "A5SS\t$_"; }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }

            my %afe = &get_afe(@two_exon_info);
            foreach (keys %afe) {
                my $afle_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle_info = "AFE\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle_info = "ALE\t$_"; }
                $statistics{$gene_id}{"afle"}{$afle_info} = 1;
            }

            my %ale = &get_ale(@two_exon_info);
            foreach (keys %ale) {
                my $afle_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle_info = "ALE\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle_info = "AFE\t$_"; }
                $statistics{$gene_id}{"afle"}{$afle_info} = 1;
            }

            my %afe1 = &get_afe1(@two_exon_info);
            foreach (keys %afe1) {
                my $afle1_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle1_info = "AFE1\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle1_info = "ALE1\t$_"; }
                $statistics{$gene_id}{"afle1"}{$afle1_info} = 1;
            }

            my %ale1 = &get_ale1(@two_exon_info);
            foreach (keys %ale1) {
                my $afle1_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") { $afle1_info = "ALE1\t$_"; }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") { $afle1_info = "AFE1\t$_"; }
                $statistics{$gene_id}{"afle1"}{$afle1_info} = 1;
            }
        }
    }
}

foreach  my $gene_id (sort keys %statistics) {
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"es"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"ir"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"a53ss"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"afle"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"afle1"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
}

sub get_es {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my @exons_1 = @first;
    my @exons_2 = @second;
    my %target_regions;
    while (my $exon1_1 = shift @exons_1) {
        my ($target_region_start, $target_region_end);
        my $exon1_1_end = $1 if $exon1_1 =~ m/\t(\d+)/;
        my @exons_2_tmp = @exons_2;
        while (my $exon2_1 = shift @exons_2_tmp) {
            my ($start, $end) = ($1, $2) if $exon2_1 =~ m/(\d+)\t(\d+)/;
            if ($end == $exon1_1_end) {
                $target_region_start = $exon1_1_end;
                @exons_2 = @exons_2_tmp;
                last;
            }
        }

T1:     foreach my $exon1_2 (@exons_1) {
            my $exon1_2_start = $1 if $exon1_2 =~ m/(\d+)/;
            foreach my $exon2_2 (@exons_2) {
                my $exon2_2_start = $1 if $exon2_2 =~ m/(\d+)/;
                if ($exon1_2_start == $exon2_2_start) {
                    $target_region_end = $exon1_2_start;
                    last T1;
                }
            }
        }

        $target_regions{"$target_region_start\t$target_region_end"} = 1 if ($target_region_start && $target_region_end);
    }

    foreach (sort {$a <=> $b} keys %target_regions) {
        my ($target_region_start, $target_region_end) = ($1, $2) if m/(\d+)\t(\d+)/;

        my @exons_1 = @first;
        my @exons_2 = @second;
        my ($exon_skipping_number1, $exon_skipping_number2, $exon_skipping_info1, $exon_skipping_info2, %positions1, %positions2);
T2:     while (my $exon1_1 = shift @exons_1) {
            my $exon1_1_end = $1 if $exon1_1 =~ m/\t(\d+)/;
            if ($exon1_1_end == $target_region_start) {
                while (my $exon1_2 = shift @exons_1) {
                    my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
                    if ($exon1_2_start != $target_region_end) {
                        $positions1{"$exon1_2_start\t$exon1_2_end"} = 1;
                        $exon_skipping_number1 ++;
                        $exon_skipping_info1 .= "$exon1_2_start-$exon1_2_end;";
                    }
                    else {
                        last T2;
                    }
                }
            }
        }

T3:     while (my $exon2_1 = shift @exons_2) {
            my $exon2_1_end = $1 if $exon2_1 =~ m/\t(\d+)/;
            if ($exon2_1_end == $target_region_start) {
                while (my $exon2_2 = shift @exons_2) {
                    my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;
                    if ($exon2_2_start != $target_region_end) {
                        $positions2{"$exon2_2_start\t$exon2_2_end"} = 1;
                        $exon_skipping_number2 ++;
                        $exon_skipping_info2 .= "$exon2_2_start-$exon2_2_end;";
                    }
                    else {
                        last T3;
                    }
                }
            }
        }

        if ($exon_skipping_number1 >= 1 && $exon_skipping_number2 >= 1) {
            my $keep = 1;
            foreach (keys %positions1) {
                my ($exon1_start, $exon1_end) = split /\t/;
                foreach (keys %positions2) {
                    my ($exon2_start, $exon2_end) = split /\t/;
                    $keep = 0 if ($exon1_start <= $exon2_end && $exon2_start <= $exon1_end);
                }
            }
            my $info;
            if ($exon_skipping_info2 > $exon_skipping_info1) { $info = "$exon_skipping_info1\t$exon_skipping_info2"; }
            else { $info = "$exon_skipping_info2\t$exon_skipping_info1"; }
            $return{"MEE\t$target_region_start\t$info\t$target_region_end"} = 1 if $keep == 1;
        }
        elsif ($exon_skipping_number1 >= 1) {
            $return{"ES\t$target_region_start\t$exon_skipping_info1\t$target_region_end"} = 1;
        }
        elsif ($exon_skipping_number2 >= 1) {
            $return{"ES\t$target_region_start\t$exon_skipping_info2\t$target_region_end"} = 1;
        }
    }
    return %return;
}

sub get_ir {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my @exons_1 = @first;
    my @exons_2 = @second;

    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my $region_start = $1 + 1 if $exon =~ m/\t(\d+)/;
            my $region_end = $1 - 1 if $next_exon =~ m/(\d+)/;
            #print "$region_start\t$region_end\n";

            foreach (@exons_2) {
                #print "$record\t";
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                next if $end < $region_start;
                last if $start > $region_end;
                $return{"IR\t$region_start-$region_end"} = 1 if ($start < $region_start && $end > $region_end);
            }
        }
    }

    my @exons_1 = @second;
    my @exons_2 = @first;

    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my $region_start = $1 + 1 if $exon =~ m/\t(\d+)/;
            my $region_end = $1 - 1 if $next_exon =~ m/(\d+)/;
            #print "$region_start\t$region_end\n";

            foreach (@exons_2) {
                #print "$record\t";
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                next if $end < $region_start;
                last if $start > $region_end;
                $return{"IR\t$region_start-$region_end"} = 1 if ($start < $region_start && $end > $region_end);
            }
        }
    }

    return %return;
}

sub get_a5ss {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my @exons_1 = @first;
    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my ($region_start_before, $region_start) = ($1, $2) if $exon =~ m/(\d+)\t(\d+)/;
            my $region_end = $1 if $next_exon =~ m/(\d+)/;

            my @exons_2 = @second;
            my $exon5p = shift @exons_2;
            my ($exon5p_start, $exon5p_end) = ($1, $2) if $exon5p =~ m/(\d+)\t(\d+)/;
            foreach (@exons_2) {
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                if ($start == $region_end) {
                    $return{"$start\t$region_start\t$exon5p_end"} = 1 if ($exon5p_end > $region_start && $exon5p_start < $region_start);
                    last;
                }
                $exon5p_start = $start;
                $exon5p_end = $end;
            }
        }
    }

    my @exons_1 = @second;
    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my ($region_start_before, $region_start) = ($1, $2) if $exon =~ m/(\d+)\t(\d+)/;
            my $region_end = $1 if $next_exon =~ m/(\d+)/;

            my @exons_2 = @first;
            my $exon5p = shift @exons_2;
            my ($exon5p_start, $exon5p_end) = ($1, $2) if $exon5p =~ m/(\d+)\t(\d+)/;
            foreach (@exons_2) {
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                if ($start == $region_end) {
                    $return{"$start\t$region_start\t$exon5p_end"} = 1 if ($exon5p_end > $region_start && $exon5p_start < $region_start);
                    last;
                }
                $exon5p_start = $start;
                $exon5p_end = $end;
            }
        }
    }

    return %return;
}

sub get_a3ss {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my @exons_1 = @first;
    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my $region_start = $1 if $exon =~ m/\t(\d+)/;
            my ($region_end, $region_end_after) = ($1, $2) if $next_exon =~ m/(\d+)\t(\d+)/;

            my @exons_2 = @second;
            my $exon5p = shift @exons_2;
            my ($exon5p_start, $exon5p_end) = ($1, $2) if $exon5p =~ m/(\d+)\t(\d+)/;
            foreach (@exons_2) {
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                if ($exon5p_end == $region_start) {
                    $return{"$exon5p_end\t$start\t$region_end"} = 1 if ($start < $region_end && $end > $region_end);
                    last;
                }
                $exon5p_start = $start;
                $exon5p_end = $end;
            }
        }
    }

    my @exons_1 = @second;
    while (my $exon = shift @exons_1) {
        if (@exons_1) {
            my $next_exon = $exons_1[0];
            my $region_start = $1 if $exon =~ m/\t(\d+)/;
            my ($region_end, $region_end_after) = ($1, $2) if $next_exon =~ m/(\d+)\t(\d+)/;

            my @exons_2 = @first;
            my $exon5p = shift @exons_2;
            my ($exon5p_start, $exon5p_end) = ($1, $2) if $exon5p =~ m/(\d+)\t(\d+)/;
            foreach (@exons_2) {
                my ($start, $end) = ($1, $2) if m/(\d+)\t(\d+)/;
                if ($exon5p_end == $region_start) {
                    $return{"$exon5p_end\t$start\t$region_end"} = 1 if ($start < $region_end && $end > $region_end);
                    last;
                }
                $exon5p_start = $start;
                $exon5p_end = $end;
            }
        }
    }

    return %return;
}
sub get_afe {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my $exon1_1 = shift @first;
    my $exon1_2 = shift @first;
    my $exon2_1 = shift @second;
    my $exon2_2 = shift @second;

    my ($exon1_1_start, $exon1_1_end) = ($1, $2) if $exon1_1 =~ m/(\d+)\t(\d+)/;
    my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
    my ($exon2_1_start, $exon2_1_end) = ($1, $2) if $exon2_1 =~ m/(\d+)\t(\d+)/;
    my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;

    if ($exon1_2_start && $exon2_2_start && $exon1_2_start == $exon2_2_start) {
        if ($exon1_1_start != $exon2_1_start && $exon1_1_end != $exon2_1_end) {
            if ($exon1_1_start < $exon2_1_start) {
                $return{"$exon1_2_start\t$exon1_1_start-$exon1_1_end\t$exon2_1_start-$exon2_1_end"} = 1;
            }
            else {
                $return{"$exon1_2_start\t$exon2_1_start-$exon2_1_end\t$exon1_1_start-$exon1_1_end"} = 1;
            }
        }
    }
    
    return %return;
}

sub get_ale {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;
    
    my $exon1_1 = pop @first;
    my $exon1_2 = pop @first;
    my $exon2_1 = pop @second;
    my $exon2_2 = pop @second;

    my ($exon1_1_start, $exon1_1_end) = ($1, $2) if $exon1_1 =~ m/(\d+)\t(\d+)/;
    my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
    my ($exon2_1_start, $exon2_1_end) = ($1, $2) if $exon2_1 =~ m/(\d+)\t(\d+)/;
    my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;

    if ($exon1_2_end && $exon2_2_end && $exon1_2_end == $exon2_2_end) {
        if ($exon1_1_start != $exon2_1_start && $exon1_1_end != $exon2_1_end) {
            if ($exon1_1_start < $exon2_1_start) {
                $return{"$exon1_2_end\t$exon1_1_start-$exon1_1_end\t$exon2_1_start-$exon2_1_end"} = 1;
            }
            else {
                $return{"$exon1_2_end\t$exon2_1_start-$exon2_1_end\t$exon1_1_start-$exon1_1_end"} = 1;
            }
        }
    }

    return %return;
}

sub get_afe1 {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$a <=> $b} split /\n/, $first;
    my @second = sort {$a <=> $b} split /\n/, $second;

    my $exon1_1 = $first[0];
    my $exon1_2 = $first[1];
    my $exon2_1 = $second[0];
    my $exon2_2 = $second[1];

    if ($exon1_2 && $exon2_2) {
        my ($exon1_1_start, $exon1_1_end) = ($1, $2) if $exon1_1 =~ m/(\d+)\t(\d+)/;
        my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
        my ($exon2_1_start, $exon2_1_end) = ($1, $2) if $exon2_1 =~ m/(\d+)\t(\d+)/;
        my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;

        my ($last_exon_start, $last_exon_end) = ($exon2_1_start, $exon2_1_end);
        shift @second;
        foreach (@second) {
            @_ = split /\t/;
            if ($_[0] == $exon1_2_start && $last_exon_end == $exon1_1_end && $last_exon_start != $exon1_1_start) {
                $return{"$exon1_2_start\t$exon1_1_start-$exon1_1_end\t$last_exon_start-$last_exon_end"} = 1;
                last;
            }
            ($last_exon_start, $last_exon_end) = @_;
        }
        my ($last_exon_start, $last_exon_end) = ($exon1_1_start, $exon1_1_end);
        shift @first;
        foreach (@first) {
            @_ = split /\t/;
            if ($_[0] == $exon2_2_start && $last_exon_end == $exon2_1_end && $last_exon_start != $exon2_1_start) {
                $return{"$exon2_2_start\t$exon2_1_start-$exon2_1_end\t$last_exon_start-$last_exon_end"} = 1;
                last;
            }
            ($last_exon_start, $last_exon_end) = @_;
        }
    }

    return %return;
}

sub get_ale1 {
    my %return;
    my ($first, $second) = @_;
    $first =~ s/\n$//; $second =~ s/\n$//;
    my @first = sort {$b <=> $a} split /\n/, $first;
    my @second = sort {$b <=> $a} split /\n/, $second;

    my $exon1_1 = $first[0];
    my $exon1_2 = $first[1];
    my $exon2_1 = $second[0];
    my $exon2_2 = $second[1];

    if ($exon1_2 && $exon2_2) {
        my ($exon1_1_start, $exon1_1_end) = ($1, $2) if $exon1_1 =~ m/(\d+)\t(\d+)/;
        my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
        my ($exon2_1_start, $exon2_1_end) = ($1, $2) if $exon2_1 =~ m/(\d+)\t(\d+)/;
        my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;

        my ($last_exon_start, $last_exon_end) = ($exon2_1_start, $exon2_1_end);
        shift @second;
        foreach (@second) {
            @_ = split /\t/;
            if ($_[1] == $exon1_2_end && $last_exon_start == $exon1_1_start && $last_exon_end != $exon1_1_end) {
                $return{"$exon1_2_end\t$exon1_1_start-$exon1_1_end\t$last_exon_start-$last_exon_end"} = 1;
                last;
            }
            ($last_exon_start, $last_exon_end) = @_;
        }
        my ($last_exon_start, $last_exon_end) = ($exon1_1_start, $exon1_1_end);
        shift @first;
        foreach (@first) {
            @_ = split /\t/;
            if ($_[1] == $exon2_2_end && $last_exon_start == $exon2_1_start && $last_exon_end != $exon2_1_end) {
                $return{"$exon2_2_end\t$exon2_1_start-$exon2_1_end\t$last_exon_start-$last_exon_end"} = 1;
                last;
            }
            ($last_exon_start, $last_exon_end) = @_;
        }
    }

    return %return;
}
