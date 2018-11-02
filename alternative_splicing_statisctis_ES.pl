#!/usr/bin/perl
use strict;
use Algorithm::Combinatorics qw(combinations permutations);

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf > alternative_splicing_statisctis_MEE.txt

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%gene, %gene_loci);
while (<IN>) {
    if (m/gene_id\s+?\"(\S+?)\".*transcript_id\s+?\"(\S+?)\"/) {
        @_ = split /\t/;
        $gene{$1}{$2} .= "$_[3]\t$_[4]\t$_[6]\n";
        $gene_loci{$1}{"scaffold"} = $_[0];
        $gene_loci{$1}{"strand"} = $_[6];
    }
}
close IN;

# 对每个基因进行AS分析，统计的AS类型有：
# 1. exon skipping                es
# 2. intron retention            ir
# 3. alternative 5' splicing    a5s
# 4. alternative 3' splicing    a3s
# 5. alternative first exon        afe
# 6. alternatie last exon        ale
# 7. mutually exclusice exons    mee

my %statistics;

foreach my $gene_id (sort keys %gene) {
    #print "$gene_id\n";
    my @isoform = keys %{$gene{$gene_id}};
    if (@isoform >= 2) {
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
        }
    }
}

foreach  my $gene_id (sort keys %statistics) {
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"es"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
}

# 子程序用于比较2个可变剪接的exon信息，得到ES和MEE信息。
# 每个ES信息由3个数据组成，这3个数据用制表符分割：共有的exon边界位点、Skipped exons首尾位置信息（多个exons信息使用分号分割）、共有的exon边界位点。这3个数据按数值从小到大排序。
# 每个MEE信息由4个数据组成，这4个数据用制表符分割：共有的exon边界位点、Skipped exons首尾位置信息（多个exons信息使用分号分割）、另一个可变剪接的Skipped exons首尾位置信息（多个exons信息使用分号分割））、共有的exon边界位点。这4个数据按数值从小到大排序。
# 找ES和MEE原理：
# （1）首先，将两个可变剪接的exon信息都按数值从小到大排序，找2个可变剪接exon结束位置（方向为从小到大的方向，不是exon的正负义链方向）一致的点，和距离该点最近的exon起始位位置一致的点，作为目标区域。
# （2）对目标区域进行exon数目分析。若2个可变剪接的目标区域中都有>=1个exons，且这些exons的边界位置都不一致，则为MEE；若仅有1个可变剪接的目标区域中有>=1个exons，则为ES。
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

T1:        foreach my $exon1_2 (@exons_1) {
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
        my ($exon_skipping_number1, $exon_skipping_number2, $exon_skipping_info1, $exon_skipping_info2, %positions);
T2:        while (my $exon1_1 = shift @exons_1) {
            my $exon1_1_end = $1 if $exon1_1 =~ m/\t(\d+)/;
            if ($exon1_1_end == $target_region_start) {
                while (my $exon1_2 = shift @exons_1) {
                    my ($exon1_2_start, $exon1_2_end) = ($1, $2) if $exon1_2 =~ m/(\d+)\t(\d+)/;
                    if ($exon1_2_start != $target_region_end) {
                        $positions{$exon1_2_start} = 1;
                        $positions{$exon1_2_end} = 1;
                        $exon_skipping_number1 ++;
                        $exon_skipping_info1 .= "$exon1_2_start-$exon1_2_end;";
                    }
                    else {
                        last T2;
                    }
                }
            }
        }

T3:        while (my $exon2_1 = shift @exons_2) {
            my $exon2_1_end = $1 if $exon2_1 =~ m/\t(\d+)/;
            if ($exon2_1_end == $target_region_start) {
                while (my $exon2_2 = shift @exons_2) {
                    my ($exon2_2_start, $exon2_2_end) = ($1, $2) if $exon2_2 =~ m/(\d+)\t(\d+)/;
                    if ($exon2_2_start != $target_region_end) {
                        $positions{$exon2_2_start} = 1;
                        $positions{$exon2_2_end} = 1;
                        $exon_skipping_number2 ++;
                        $exon_skipping_info2 .= "$exon2_2_start-$exon2_2_end;";
                    }
                    else {
                        last T3;
                    }
                }
            }
        }

        my @positions = keys %positions;
        my $positions_number = @positions;
        if ($positions_number > 0 && ($positions_number == ($exon_skipping_number1 + $exon_skipping_number2) * 2)) {
            $exon_skipping_info1 =~ s/;$//;
            $exon_skipping_info2 =~ s/;$//;
            if ($exon_skipping_number1 >= 1 && $exon_skipping_number2 >= 1) {
                my $info;
                if ($exon_skipping_info2 > $exon_skipping_info1) { $info = "$exon_skipping_info1\t$exon_skipping_info2"; }
                else { $info = "$exon_skipping_info2\t$exon_skipping_info1"; }
                $return{"MEE\t$target_region_start\t$info\t$target_region_end"} = 1;
            }
            elsif ($exon_skipping_number1 >= 1) {
                $return{"ES\t$target_region_start\t$exon_skipping_info1\t$target_region_end"} = 1;
            }
            elsif ($exon_skipping_number2 >= 1) {
                $return{"ES\t$target_region_start\t$exon_skipping_info2\t$target_region_end"} = 1;
            }
        }
    }
    return %return;
}
