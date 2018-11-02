#!/usr/bin/perl
use strict;
use Algorithm::Combinatorics qw(combinations permutations);

my $usage = <<USAGE;
Usage:
    perl $0 spliceGrapher.gtf > alternative_splicing_statisctis_IR.txt

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

            my %ir = &get_ir(@two_exon_info);
            foreach (keys %ir) {
                $statistics{$gene_id}{"ir"}{$_} = 1;
            }
        }
    }
}

foreach  my $gene_id (sort keys %statistics) {
    foreach (sort {$a <=> $b} keys %{$statistics{$gene_id}{"ir"}}) {
        print "$gene_id\t$gene_loci{$gene_id}{\"scaffold\"}\t$gene_loci{$gene_id}{\"strand\"}\t$_\n";
    }
}

# 子程序用于比较2个可变剪接的exon信息，得到IR信息。
# 每个IR信息由1个数据组成，即IR的intron位置信息。
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
