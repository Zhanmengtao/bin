#!/usr/bin/perl
use strict;
use FileHandle;

my $usage = <<USAGE;
Usage:
    perl $0 base_num file1.fastq file2.fastq ... > out.fastq

    依次从输入的多个fastq文件中读取一条read，存放到out.fastq文件，再继续依次从多个fastq文件中读取下一条read，存放到out.fastq文件中... 
    若out.fastq文件中的碱基数量超过设定的base_num后，则停止读取数据。

USAGE
if (@ARGV==0){die $usage}

my $target_num = shift @ARGV;
my %fh;
foreach (@ARGV) {
    $fh{$_} = FileHandle->new($_, "r");
}

my $num = 0;
my $reads_num = 0;
my $last = 1;
while ($last) {
    my $over = 1;
    foreach my $file (@ARGV) {
        my $fh = $fh{$file};

        $_ = <$fh>; print;
        if ($_) {
            $reads_num ++;
            $over = 0;
        }

        $_ = <$fh>; print;
        $num += length($_) - 1;
        if ($num >= $target_num) {
            $last = 0;
        }

        $_ = <$fh>; print;
        $_ = <$fh>; print;

        last if $last == 0;
    }
    $last = 0 if $over == 1;
}
print STDERR "$reads_num Reads were selected with $num nuclotide bases\n";
