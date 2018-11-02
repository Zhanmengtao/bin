#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    Perl $0 proteins.fasta > new_proteins.fasta

    --no_sort
        设置不对序列按长度进行排序。默认设置下，程序会对输入的序列按从长到短进行排序。添加该参数则不进行排序。
    --no_rename
        设置不对序列进行重命名。若添加该参数，表示不会对序列名称进行重命名，参数--seq_prefix不会生效。
    --seq_prefix <String>  default: "proteins"
        设置--no_rename参数后，该参数失效。重命名后的序列名称以该指定的参数为前缀，后接逐一递增的数字编号，编号前用数字0补齐以使所有序列数字编号的字符数一致。
    --min_length <Int>  default: 50
        设置最短的序列长度。丢弃长度低于此阈值的序列。
    --line_length <Int>  default: 60
        设置输出的fasta文件中，序列在每行的最大字符长度。若该值 < 0，则表示不对序列进行换行处理。
    
USAGE
if (@ARGV==0){die $usage}

my ($no_sort, $no_rename, $seq_prefix, $min_length, $line_length);
GetOptions(
    "no_sort" => \$no_sort,
    "no_rename" => \$no_rename,
    "seq_prefix:s" => \$seq_prefix,
    "min_length:i" => \$min_length,
    "line_length:i" => \$line_length,
);
$seq_prefix ||= "proteins";
$min_length ||= 50;
$line_length ||= 60;

open IN, @ARGV[0] or die $!;
my (%seq, $seq_id, @seq_id, %seq_id, %seq_length);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) {
        $seq_id = $1;
        if (exists $seq_id{$seq_id}) {
            my $num = $seq_id{$seq_id};
            $seq_id{$seq_id} ++;
            print STDERR "Warning: $seq_id appears $seq_id{$seq_id} times! forcibly rename this sequence id to $seq_id\_$num\n";
            $seq_id = "$seq_id\_$num";
        }
        $seq_id{$seq_id} ++;
        push @seq_id, $seq_id;
    }
    else {
        $seq{$seq_id} .= $_;
        $seq_length{$seq_id} += length;
    }
}
close IN;

unless ($no_sort) {
    @seq_id = sort {$seq_length{$b} <=> $seq_length{$a}} @seq_id;
}

my @old_seq_id = @seq_id;
@seq_id = ();
foreach (@old_seq_id) {
	push @seq_id, $_ if $seq_length{$_} >= $min_length;
}

my $number = 0;
foreach my $id (@seq_id) {
    if ($seq_length{$id} >= $min_length) {
        $number ++;
        my $seq_name = $id;
        unless ($no_rename) {
            $seq_name = $seq_prefix . "0" x ((length @seq_id) - (length $number)) . $number;
        }
        my $seq = $seq{$id};
        if ($line_length > 0) {
            $seq =~ s/(\w{$line_length})/$1\n/g;
            $seq =~ s/\n$//;
        }
        print ">$seq_name\n$seq\n";
    }
    else {
        next;
    }
}
