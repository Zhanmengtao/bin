#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] pacbio_subreads0.fasta pacbio_subreads1.fasta ...

For example:
    perl $0 --genome_size 4600000 --coverage 25 --out out.fasta pacbio_subreads0.fasta pacbio_subreads1.fasta

    --out <string>    default:out.fasta
    设置输出文件。

    --genome_size <int>
    设置基因组大小。

    --coverage <float>
    设置提取相应覆盖深度的数据量。

    --size <int>
    设置提取相应大小的数据量。若设置了该参数，则以上两个参数失效。
USAGE
if (@ARGV==0){die $usage}

my ($out, $genome_size, $coverage, $size);
GetOptions(
    "out:s" => \$out,
    "genome_size:i" => \$genome_size,
    "coverage:f" => \$coverage,
    "size:i" => \$size,
);
$out ||= "out.fasta";
die "Wrong parameters\n" unless (($genome_size && $coverage) || $size);

my (%seq, $seq_id, %len, $total);
foreach (<>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else {
        $seq{$seq_id} .= $_;
        $len{$seq_id} += length($_);
        $total += length($_)
    }
}
print STDERR "Input $total base in total\n";

$size = $genome_size * $coverage unless $size;
open OUT, ">", $out or die $!;
my ($read_number, $bp);
foreach (sort {$len{$b} <=> $len{$a}} keys %seq) {
    $read_number ++;
    $bp += $len{$_};
    print OUT ">$_\n$seq{$_}\n";
    last if $bp >= $size;
}
close OUT;
print STDERR "$read_number reads were picked out and $bp base in total\n";
