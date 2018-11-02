#!/usr/bin/perl
use strict;
# 2014.06.30	Lianfu Chen		chenllianfu@foxmail.com

my $usage = <<USAGE;
Usage:
	perl $0 tRNA.out tRNA.ss > tRNA.gff3

USAGE
if (@ARGV==0){die $usage}

my (%second_structure, $ss_seqID, $ss_start, $ss_seq);
open IN, $ARGV[1] or die $!;
while (<IN>) {
	if (/^(\S+)\.trna\d+ \((\d+)-(\d+)\)/) {
		$ss_seqID = $1;
		if ($3 > $2) { $ss_start = $2 } else { $ss_start = $3 }
	}
	elsif (/^Seq: (\w+)/) { $second_structure{$ss_seqID}{$ss_start}{"seq"} = $1 }
	elsif (/^Str: (\S+)/) { $second_structure{$ss_seqID}{$ss_start}{"str"} = $1 }
}

my %gff3;
open IN, $ARGV[0] or die $!;
while (<IN>) {
	my @annot = split /\s+/, $_;
	if ($annot[2] =~ m/\d+/) {
		my ($seqID,$start,$end,$intron_start,$intron_end,$type,$score,$anti_codon,$strand);
		if ($annot[3] > $annot[2]) { $start = $annot[2]; $end = $annot[3]; $strand = "+"; $intron_start = $annot[6] - 1; $intron_end = $annot[7] + 1; }
		else { $start = $annot[3]; $end = $annot[2]; $strand = "-"; $intron_start = $annot[7] - 1; $intron_end = $annot[6] + 1; }
		$seqID = $annot[0]; $type = $annot[4]; $score = $annot[8]; $anti_codon = $annot[5];
		my $attribute = "Name=tRNA-$type;Anti-codon=$anti_codon;Sequence=$second_structure{$seqID}{$start}{\"seq\"};Structure=$second_structure{$seqID}{$start}{\"str\"};";

		if ($type eq "Pseudo") {
			if ($intron_start == -1) {
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$end\t$score\t$strand\t.\tID=pseudo_number_tRNA-$type;$attribute\n";
			}
			else {
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$intron_start\t$score\t$strand\t.\tID=pseudo_number_tRNA-$type;$attribute\n";
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$intron_end\t$end\t$score\t$strand\t.\tID=pseudo_number_tRNA-$type;$attribute\n";
			}
		}
		else {
			if ($intron_start == -1) {
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$end\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
			}
			else {
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$intron_start\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
				$gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$intron_end\t$end\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
			}
		}
	}
}
close IN;

my @seqID = keys %gff3;
my %sort_seqID;
foreach (@seqID) { $sort_seqID{$_} = $1 if /(\d+)\D*$/ }
@seqID = sort {$sort_seqID{$a} <=> $sort_seqID{$b}} @seqID;

my ($pseudo_num, $tRNA_num) = (0, 0);
foreach my $seqID (@seqID) {
	foreach (sort {$a <=> $b} keys %{$gff3{$seqID}}) {
		my $out = $gff3{$seqID}{$_};
		if ($out =~ m/pseudo/) {
			$pseudo_num ++;
			$out =~ s/number/$pseudo_num/g;
			print $out;
		}
		else {
			$tRNA_num ++;
			$out =~ s/number/$tRNA_num/g;
			print $out;
		}
	}
}

print STDERR "$pseudo_num Pseudo tRNAs\n$tRNA_num tRNAs\n";
