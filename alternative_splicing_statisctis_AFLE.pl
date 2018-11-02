#!/usr/bin/perl
use strict;
use Algorithm::Combinatorics qw(combinations permutations);

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf > alternative_splicing_statisctis_AFLE.txt

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
        }
    }
}

foreach  my $gene_id (sort keys %statistics) {
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"afle"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
}

# 子程序用于比较2个可变剪接的exon信息，得到AFE和ALE信息。
# AFE信息由3个数据组成，这3个数据用制表符分割：共同的intron的3'端位点、第一个exon的首尾位点、可变剪接的第一个exon的首尾位点。后2个数据按数值从小到大排序。
# ALE信息由3个数据组成，这3个数据用制表符分割：共同的intron的5'端位点、最后一个exon的首尾位点、可变剪接的最后一个exon的首尾位点。后2个数据按数值从小到大排序。
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

    if ($exon1_2_start == $exon2_2_start) {
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

    if ($exon1_2_end == $exon2_2_end) {
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
