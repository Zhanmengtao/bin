#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 assemblies.fasta.transdecoder.genome.gff3 pasa_bamboo.assemblies.fasta.transdecoder.cds Min_exon_number_of_CDS Min_cds_length Min_cds_ratio_vs_cDNA Max_intron_length_ratio Max_cds_length_ratio

For Example:
    perl $0 assemblies.fasta.transdecoder.genome.gff3 pasa_bamboo.assemblies.fasta.transdecoder.cds 3 900 0.60 0.95 0.95 > best_candidates.gff3

    1. 按照先后顺序对assemblies.fasta.transdecoder.genome.gff3中的基因模型进行过滤。
    2. 首先，根据pasa_bamboo.assemblies.fasta.transdecoder.cds文件中的complete关键词过滤不完整基因模型。
    3. 再过滤cds数目<3的基因模型。
    4. 再过滤cds长度<900的基因模型。
    5. 再过滤cds region占exon region比例<0.6的基因模型。
    6. 再过滤含有过长intron的基因模型。根据全部数据的cds(而不是exon)信息得到的intron长度信息；再将intron长度按从小到大进行排序，选取95%分位数为阈值。
    7. 最后过滤cds长度过长的基因模型。根据全部数据得到的cds长度信息；将cds长度长度从小到大排序，选取95%分为数为阈值。

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[1] or die $!;
my %complete;
while (<IN>) {
    $complete{$1} = 1 if m/>(\S+).*type:complete/;
}
close IN;

open IN, $ARGV[0] or die $!;
$/ = "\n\n";
my (@intron_length, @cds_length, %validate);
my ($gene_number_total, $gene_nuber_validate, $gene_number_uncomplete, $gene_number_exonLack, $gene_number_short, $gene_number_cdsRatio, $gene_number_intron, $gene_number_long);
while (<IN>) {
    $gene_number_total ++;

    if (m/mRNA\t.*?ID=(\S+?);/ && exists $complete{$1}) {
        my @cds_line = m/(CDS\t\d+?\t\d+?)\t/g;
        my @exon_line = m/(exon\t\d+?\t\d+?)\t/g;
        push @intron_length, &get_intron_length(@cds_line);

        my $cds_line = @cds_line;
        if ($cds_line < $ARGV[2]) {
            $gene_number_exonLack ++;
            next;
        }

        my ($cds_length, $exon_length);
        foreach (@cds_line) {
            m/CDS\t(\d+)\t(\d+)/;
            $cds_length += ($2 - $1 + 1);
        }
        push @cds_length, $cds_length;
        foreach (@exon_line) {
            m/exon\t(\d+)\t(\d+)/;
            $exon_length += ($2 - $1 + 1);
        }

        if ($cds_length >= $ARGV[3]) {
            my $cdsRation = $cds_length / $exon_length;
            if ($cdsRation >= $ARGV[4]) {
                $validate{$_} = 1;
            }
            else {
                $gene_number_cdsRatio ++;
            }
        }
        else {
            $gene_number_short ++;
        }
    }
    else {
        $gene_number_uncomplete ++;
    }
}

my @intron_length = sort {$a <=> $b} @intron_length;
my $intron_length_threshold = $intron_length[@intron_length*$ARGV[5]];
print STDERR "Max inron length threshold was set to: $intron_length_threshold\n";
my @cds_length = sort {$a <=> $b} @cds_length;
my $cds_length_threshold = $cds_length[@cds_length*$ARGV[6]];
print STDERR "Max CDS length threshold was set to:   $cds_length_threshold\n";

foreach (keys %validate) {
    my @cds_region = m/\tCDS\t(\d+?\t\d+?)\t/g;
    my $cds_length;
    foreach (@cds_region) {
        @_ = split /\t/;
        $cds_length += ($_[1] - $_[0] + 1);
    }

    my @intron_length = &get_intron_length(@cds_region);
    @intron_length = sort {$b <=> $a} @intron_length;

    if ($intron_length[0] <= $intron_length_threshold) {
        if ($cds_length <= $cds_length_threshold) {
            $gene_nuber_validate ++;
            print;
        }
        else {
            $gene_number_long ++;
        }
    }
    else {
        $gene_number_intron ++;
    }
}
print STDERR "
Total gene number:                  $gene_number_total
Validate gene number:               $gene_nuber_validate
Uncomplete gene number:             $gene_number_uncomplete
Exon < $ARGV[2] gene number:               $gene_number_exonLack
CDS length < $ARGV[3] gene number:       $gene_number_short
CDS region ratio < $ARGV[4] gene number: $gene_number_cdsRatio
Intron length > $intron_length_threshold gene number:   $gene_number_intron
CDS length > $cds_length_threshold gene number:      $gene_number_long\n\n";

sub get_intron_length {
    my @exon;
    foreach (@_) {
        m/(\d+)\t(\d+)$/;
        push @exon, "$1\t$2";
    }
    @exon = sort {$a <=> $b} @exon;
    my $first = shift @exon;
    my $end = $1 if $first =~ m/(\d+)$/;
    my @intron_length;
    foreach (@exon) {
        m/^(\d+)\t(\d+)$/;
        my $intron_length = $1 - $end;
        $end = $2;
        push @intron_length, $intron_length;
    }
    return @intron_length;
}
