#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    diff2file.pl file1.txt file2.txt

	--out     output 3 files: file_common, file1_only, file2_only. default: false
	--help        output help information to screen

USAGE
if (@ARGV==0){die $usage}

my ($out,$Help);

GetOptions(
	"out"=>\$out,
	"help"=>\$Help
);

die `pod2text $0` if (@ARGV == 0 || $Help);

open INPUT1, $ARGV[0] or die $!;
open INPUT2, $ARGV[1] or die $!;

my @file1 = <INPUT1>;
my @file2 = <INPUT2>;

my (%file1,%file1_only,%file2,%file2_only,%file_common);

foreach (@file1) {
	$file1{$_} = 1;
}
foreach (@file2) {
	$file2{$_} = 1;
}

foreach (@file1) {
	if (exists $file2{$_}) {
		$file_common{$_} = 1;
	}else {
		$file1_only{$_} = 1;
	}
}

foreach (@file2) {
	unless (exists $file_common{$_}) {
		$file2_only{$_} = 1;
	}
}

my @file1_only_num = keys %file1_only;
my @file2_only_num = keys %file2_only;
my @file_common = keys %file_common;

my $file1_only_num = @file1_only_num;
my $file2_only_num = @file2_only_num;
my $file_common = @file_common;

print "file_common\tfile1_only_num\tfile2_only_num\n$file_common\t$file1_only_num\t$file2_only_num\n";

if ($out) {
	open OUTPUT1, '>', "file1_only" or die $!;
	open OUTPUT2, '>', "file2_only" or die $!;
	open OUTPUT3, '>', "file_common" or die $!;
	print OUTPUT1 @file1_only_num;
	print OUTPUT2 @file2_only_num;
	print OUTPUT3 @file_common;
}
