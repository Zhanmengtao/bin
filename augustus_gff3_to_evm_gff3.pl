#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 augustus.gff3 > augustus_for_evm.gff3

    当Augustus的结果有可变剪接信息的时候，使用EVM软件自带的程序augustus_GFF3_to_EVM_GFF3.pl的结果有问题。
    若有可变剪接，使用本程序则只给出得分最高的可变剪接。
    此外，程序根据CDS信息转换出exon信息，并去掉intron、start_codon和stop_codon信息。
    不适合有UTR预测的结果文件。

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my (%transcript, %gene, @gene);
while (<IN>) {
    next if m/^#/;
    next if m/^\s*$/;

    if (m/\ttranscript\t\d+\t\d+\t(\S+)\t.*ID=([^;]+);Parent=([^;\s]+)/) {
        $transcript{$2}{"info"} .= $_;
        $transcript{$2}{"score"} = $1;

        push @gene, $3 unless exists $gene{$3};
        $gene{$3}{$2} = 1;
    }
    elsif (m/Parent=([^;\s]+)/) {
        $transcript{$1}{"info"} .= $_;
    }
}
close IN;

foreach my $gene_id (@gene) {
    my @transcript_id = sort {$transcript{$b}{"score"} <=> $transcript{$a}{"score"}} keys %{$gene{$gene_id}};
    my $out = $transcript{$transcript_id[0]}{"info"};

    $out =~ s/(.*)\ttranscript\t(\d+?)\t(\d+?)\t(.*)\tID/$1\tgene\t$2\t$3\t$4\tID=$gene_id\n$1\tmRNA\t$2\t$3\t$4\tID/;
    my @out = split /\n/, $out;
    print "$out[0]\n$out[1]\n";

    my @cds = grep /\tCDS\t/, @out;
    my $strand = $1 if $cds[0] =~ m/\t([+-])\t/;
    my %sort;
    foreach (@cds) { @_ = split /\t/; $sort{$_} = $_[3]; }

    if ($strand eq "+") {
        @cds = sort {$sort{$a} <=> $sort{$b}} @cds;
    }
    elsif ($strand eq "-") {
        @cds = sort {$sort{$b} <=> $sort{$a}} @cds;
    }

    my $exon_number = 0;
    foreach (@cds) {
        $exon_number ++;
        @_ = split /\t/;
        print "$_[0]\t$_[1]\texon\t$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]\tID=$transcript_id[0].exon$exon_number;Parent=$transcript_id[0]\n";
        print "$_\n";
    }
    print "\n";
}
