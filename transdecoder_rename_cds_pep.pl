#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 trinity.fasta.transdecoder.cds trinity.fasta.transdecoder.pep

USAGE
if (@ARGV==0){die $usage}

open CDS, '<', $ARGV[0] or die $!;
open PEP, '<', $ARGV[1]  or die $!;

my ($cdsID, $pepID, %cds, %pep, %orf);

while (<CDS>) {
	chomp;
	if (/^>/) {
		if (/^>((\S+)\.p\d+)\s.*type:(\S+?)\slen:(\d+?)\s/) {
			$cdsID = $1;
			$orf{$2}{$1}{"type"} = $3;
			$orf{$2}{$1}{"length"} = $4;
		}
	}else{
		$cds{$cdsID} .= $_
	}
}

print STDERR "Read CDS file OVER\n";

while (<PEP>) {
	chomp;
	if (/^>((\S+)\.p\d+)\s.*type:(\S+?)\slen:(\d+?)\s/) { $pepID = $1 }
	else           { $pep{$pepID} .= $_ }
}
print STDERR "Read PEP file OVER\n";

close CDS;
close PEP;

open CDS, '>', 'cds.fasta' or die $!;
open PEP, '>', 'protein.fasta' or die $!;

my @titles = keys %orf;
my (%sort1, %sort2, %sort3, %sort4);
foreach (@titles) {
	if (/(\d+?)_c(\d+?)_g(\d+?)_i(\d+)/) {
		$sort1{$_} = $1;
		$sort2{$_} = $2;
		$sort3{$_} = $3;
		$sort4{$_} = $4;
	}
}
@titles = sort { $sort1{$a} <=> $sort1{$b} or $sort2{$a} <=> $sort2{$b} or $sort3{$a} <=> $sort3{$b} or $sort4{$a} <=> $sort4{$b} } @titles;

foreach my $title (@titles) {
	my $title1 = $title;
	$title1 =~ s/\|/_/;
	my %idfeature = %{$orf{$title}};
	my @ids = keys %idfeature;
	@ids = sort { $idfeature{$b}{"length"} <=> $idfeature{$a}{"length"} } @ids;

	my $orf_order = 1;
	foreach (@ids) {
		print CDS ">$title1\_orf$orf_order type:$idfeature{$_}{\"type\"} len:$idfeature{$_}{\"length\"}\n$cds{$_}\n";
		my $pepSeq = $pep{$_};
		$pepSeq =~ s/\*$//;
		print PEP ">$title1\_orf$orf_order type:$idfeature{$_}{\"type\"} len:$idfeature{$_}{\"length\"}\n$pepSeq\n";
		$orf_order ++;
	}
}

print STDERR "Results has been written to cds.fasta and protein.fasta\n";
