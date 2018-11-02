#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 isolasso.pred.gtf > custom_alignments.gff3

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my %align;
while (<IN>) {
	if (/FPKM\s+?\"(.*?)\"/ && $1 < 1) {
	   next;
	}	   
	if (/\texon\t.*transcript_id \"(.*?)\"/) {
		$align{$1} .= $_;
	}
}
close IN;

my $align_num = 0;
foreach my $mRNA_id (sort keys %align) {
	$align_num ++;
	my @content = split /\n/, $align{$mRNA_id};

	my @field = split /\t/, $content[0];
	my $strand = $field[6];

	my %sort_content;
	foreach (@content) {
		@_ = split /\t/;
		$sort_content{$_} = $_[3];
	}
	
	my $out;
	if ($strand eq "+") {
		@content = sort { $sort_content{$a} <=> $sort_content{$b} } @content;
		my $start = 1;
		foreach (@content) {
			@_ = split /\t/;
			my $end = $start + ($_[4] - $_[3]);
			$out .= "$_[0]\tSGgtf2alignedGff3\tcDNA_match\t$_[3]\t$_[4]\t100\t$strand\t\.\tID=align_$align_num;Target=$mRNA_id $start $end \+\n";
			$start = $end + 1;
		}
	}
	if ($strand eq "-") {
		@content = sort { $sort_content{$b} <=> $sort_content{$a} } @content;
		my $start = 1;
		foreach (@content) {
			@_ = split /\t/;
			my $end = $start + ($_[4] - $_[3]);
			$out = "$_[0]\tSGgtf2alignedGff3\tcDNA_match\t$_[3]\t$_[4]\t100\t$strand\t\.\tID=align_$align_num;Target=$mRNA_id $start $end \+\n" . $out;
			$start = $end + 1;
		}
	}
	print $out;
}
