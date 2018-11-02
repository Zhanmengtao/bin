#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 genome.vcf genome.fasta > VCF_InDel_primer3.out

    --InDel_distance <INT>  Default: 5
        设置InDel相差碱基数大于此值的时候，才进行引物设计。

    --flanking_length <INT>  Default: 300
        设计引物时提取InDel两侧翼该长度的序列作为Primer3输入的模板序列。

    --min_product_length <INT>  Default: 100
        引物得到的最小产物长度。

    --max_product_length <INT>  Default: 250
        引物得到的最大产物长度。

    --CPU <INT>  Default:1
        程序调用ParaFly命令进行并行化计算，该参数设置并行数。

    --p3_setting_file <STRING>
        程序使用Primer3进行引物批量设计，该参数设置所使用的Primer3配置文件。若不输入该参数，Primer3默认得到的结果会很差。
        
    --gff3_out <STRING>
        可以选择输出GFF3格式的结果。

USAGE
if (@ARGV==0){die $usage}

my ($InDel_distance, $flanking_length, $min_product_length, $max_product_length, $CPU, $p3_setting_file, $gff3_out);
GetOptions(
    "InDel_distance:i" => \$InDel_distance,
    "flanking_length:i" => \$flanking_length,
    "min_product_length:i" => \$min_product_length,
    "max_product_length:i" => \$max_product_length,
    "CPU:i" => \$CPU,
    "p3_setting_file:s" => \$p3_setting_file,
    "gff3_out:s" => \$gff3_out,
);
$InDel_distance ||= 5;
$flanking_length ||= 300;
$min_product_length ||= 100;
$max_product_length ||= 250;
$CPU ||= 1;
warn "Warning: No Primer3 config file\n" unless $p3_setting_file;

open IN, $ARGV[0] or die $!;
my %indel;
while (<IN>) {
    last if m/^#CHROM/;
}
while (<IN>) {
    @_ = split /\t/;
    my $distance = abs(length($_[3]) - length($_[4]));
    #print "$_[0]\t$_[1]\t$_[3]\t$_[4]\t$distance\n" if $distance > 2;
    
    if ($distance >= $InDel_distance) {
        $indel{$_[0]}{$_[1]} = "$_[0]\t$_[1]\t$_[2]\t$_[3]\t$_[4]";
    }
}
close IN;

open IN, $ARGV[1] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (m/>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

mkdir "VCF_InDel_primer3.tmp" unless -e "VCF_InDel_primer3.tmp";
my @ID;
open COM, ">", "command.primer3.list" or die $!;
foreach my $chr (sort {$a cmp $b} keys %indel) {
    my $seq = $seq{$chr};
    foreach my $pos (sort {$a <=> $b} keys %{$indel{$chr}}) {
        my $indel_info = $indel{$chr}{$pos};
        my @indel_info = split /\t/, $indel_info;
        my $length = abs(length($indel_info[-1]) - length($indel_info[-2])) + 1;
        push @ID, "$chr\t$pos\t$length";
        my $start = $pos - $flanking_length - 1;
        my $target_pos;
        if ($start < 0) {
            $start = 0;
            $target_pos = $pos;
        }
        else {
            $target_pos = $flanking_length + 1;
        }
        my $length1 = $flanking_length * 2 + $length;
        my $subSeq = substr($seq, $start, $length1);
        open OUT, ">", "VCF_InDel_primer3.tmp/$chr\_$pos" or die $!;
        print OUT "SEQUENCE_ID=$chr\_$pos\nSEQUENCE_TEMPLATE=$subSeq\nSEQUENCE_TARGET=$target_pos,$length\nPRIMER_PRODUCT_SIZE_RANGE=$min_product_length-$max_product_length\nSEQUENCE_INTERNAL_EXCLUDED_REGION=$target_pos,$length\n=\n";
        close OUT;
        if ($p3_setting_file) {
            print COM "primer3_core -p3_settings_file $p3_setting_file -strict_tags VCF_InDel_primer3.tmp/$chr\_$pos > VCF_InDel_primer3.tmp/$chr\_$pos.out\n";
        }
        else {
            print COM "primer3_core -strict_tags VCF_InDel_primer3.tmp/$chr\_$pos > VCF_InDel_primer3.tmp/$chr\_$pos.out\n";
        }
    }
}
close COM;

my $cmdString = "ParaFly -c command.primer3.list -CPU $CPU &> /dev/null";
(system $cmdString) == 0 or die "Excute Failed: $cmdString\n";

print "CHROM\tPOS\tID\tREF\tALT\tLeftPrimer1\tTm\tSize\tRightPrimer1\tTm\tSize\tProductSize\tLeftPrimer2\tTm\tSize\tRightPrimer2\tTm\tSize\tProductSize\tLeftPrimer3\tTm\tSize\tRightPrimer3\tTm\tSize\tProductSize\tLeftPrimer4\tTm\tSize\tRightPrimer4\tTm\tSize\tProductSize\tLeftPrimer5\tTm\tSize\tRightPrimer5\tTm\tSize\tProductSize\n";
if ($gff3_out) {
    open GFF3, ">", $gff3_out or die $!; 
}
foreach (@ID) {
    @_ = split /\t/;
    my $out = "$indel{$_[0]}{$_[1]}\t";
    my $end = $_[1] + $_[2] - 1;
    my $gff3_out = "$_[0]\tLianfu_Chen\tInDel\t$_[1]\t$end\t\.\t\.\t\.\tID=$_[0]_$_[1]_$_[2];";
    open IN, "VCF_InDel_primer3.tmp/$_[0]_$_[1].out" or die $!;
    my $primer_annotation = join "", <IN>;
    my $num = 0;
    foreach my $code (0 .. 4) {
        $num ++;
        if ($primer_annotation =~ /PRIMER_LEFT_${code}_SEQUENCE=(\w+)/) {
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_left_seq=$1;";

            $primer_annotation =~ m/PRIMER_LEFT_${code}_TM=(.*)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_left_tm=$1;";

            $primer_annotation =~ m/PRIMER_LEFT_${code}=\d+,(\d+)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_left_size=$1;";

            $primer_annotation =~ m/PRIMER_RIGHT_${code}_SEQUENCE=(\w+)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_right_seq=$1;";

            $primer_annotation =~ m/PRIMER_RIGHT_${code}_TM=(.*)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_right_tm=$1;";

            $primer_annotation =~ m/PRIMER_RIGHT_${code}=\d+,(\d+)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_right_size=$1;";

            $primer_annotation =~ m/PRIMER_PAIR_${code}_PRODUCT_SIZE=(\d+)/;
            $out .= "$1\t";
            $gff3_out .= "Primer_${num}_product_size=$1;";
        }
    }
    close IN;
    print GFF3 "$gff3_out\n" if $gff3_out;
    print "$out\n";
}
