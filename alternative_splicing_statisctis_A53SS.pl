#!/usr/bin/perl
use strict;
use Algorithm::Combinatorics qw(combinations permutations);

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf > alternative_splicing_statisctis_A53SS.txt

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
                #print "$gene{$gene_id}{$_}\n";
            }

            my %a5ss = &get_a5ss(@two_exon_info);
            foreach (keys %a5ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") {
                    $a53ss_info = "A5SS\t$_";
                }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") {
                    $a53ss_info = "A3SS\t$_";
                }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }

            my %a3ss = &get_a3ss(@two_exon_info);
            foreach (keys %a3ss) {
                my $a53ss_info;
                if ($gene_loci{$gene_id}{"strand"} eq "+") {
                    $a53ss_info = "A3SS\t$_";
                }
                elsif ($gene_loci{$gene_id}{"strand"} eq "-") {
                    $a53ss_info = "A5SS\t$_";
                }
                $statistics{$gene_id}{"a53ss"}{$a53ss_info} = 1;
            }
        }
    }
}

foreach  my $gene_id (sort keys %statistics) {
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"a53ss"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
}

# 子程序用于比较2个可变剪接的exon信息，得到A5SS和A3SS信息。
# 每个A5SS信息由3个数据组成，这3个数据用制表符分割：intron的3'端位点、intron的5'端位点、intron的5'端可变剪接位点。后2个数据按数值从小到大排序。
# 每个A3SS信息由3个数据组成，这3个数据用制表符分割：intron的5'端位点、intron的3'端位点、intron的3'端可变剪接位点。后2个数据按数值从小到大排序。
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
