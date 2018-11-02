#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
	perl $0 [options] genome.fasta hints.gff

	--minsize   genome scaffolds are extracted to a split fasta, until the total sequence length
	            of this fasta > this threashold.
	--output    the output file and directory

USAGE
if (@ARGV==0){die $usage}

my ($minsize, $output);
GetOptions(
	"minsize:s"=>\$minsize,
	"output:s"=>\$output
);
$minsize ||= 1000000;
$output ||= 'split';

open IN, $ARGV[0] or die "Cannot open the gneome fasta file!\n$!";
my ($id, @id, %seq, %length);
while (<IN>) {
	chomp;
	if (/>(\S+)/) { $id = $1; push @id, $id; }
	else          { $seq{$id} .= $_; $length{$id} += length; }
}
close IN;
@id = sort { $length{$b} <=> $length{$a} } @id;

open IN, $ARGV[1] or die "Cannot open the hint file~\n$!";
my %hints;
while (<IN>) {
	if (/^(\S+)/) {
		$hints{$1} .= $_;
	}
}
close IN;

mkdir $output unless -e $output;
my $split_number = 1;
my $split_length = 0;
foreach (@id) {
	open OUT, ">>", "$output/genome.split.$split_number.fa" or die $!;
	print OUT ">$_\n$seq{$_}\n";
	close OUT;
	open OUT, ">>", "$output/genome.split.$split_number.hints" or die $!;
	print OUT $hints{$_};
	close OUT;

	$split_length += $length{$_};
	if ($split_length >= $minsize) {
		$split_length = 0;
		$split_number ++;
	}
}
