#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 reads.fastq > reads_in_phred33.fastq
	script will detect the format of fastq file, and convert format to phred33 whether the origion format is phred33 or not.

USAGE
if (@ARGV==0) {die $usage}

my $convert = 1;
open IN, $ARGV[0] or die $!;
my $line_num = 1;
while ($line_num < 400000) {
	$_ = <IN>;$_ = <IN>;$_ = <IN>;$_ = <IN>;
	if (/=/) {
		$convert = 0;
		my $line_num = $line_num * 4;
		print STDERR "character '=' was found, so the origional format is Phred33, and not needed to transform\nyou can check line $line_num of $ARGV[0] :\n$_";
		last;
	}
	$line_num ++;
}
close IN;

if ($convert eq 1) {
	print STDERR "character '=' was not found in the first 100000 reads of $ARGV[0], so the origional format was evaluated to be Phred64, and now transforming begin ...\n";
	open IN, $ARGV[0] or die $!;
	while (<IN>) {
		s/ .*//;
		s/#.*//;
		print;
		$_ = <IN>;
		print;
		$_ = <IN>;
		print "+\n";
		$_ = <IN>;
		chomp;
		my $out = "";
		my @chrs = split //;
		foreach (@chrs) {
			$out .= chr((ord($_) - 64) + 33);
		}
		print $out . "\n";
	}
	close IN;
}
else {
	open IN, $ARGV[0] or die $!;
	while (<IN>) {
		s/ .*//;
		s/#.*//;
		print;
		$_ = <IN>;
		print;
		$_ = <IN>;
		print "+\n";
		$_ = <IN>;
		print;
	}
	close IN;
}
