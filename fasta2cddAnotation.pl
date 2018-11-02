#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 proteins.fasta CPUs

USAGE
if (@ARGV==0) {die $usage}

my ($fasta, $cpus) = @ARGV;
my ($seqID, @seqID, %seq);
open IN, '<', $fasta or die $!;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) { $seqID = $1; $seqID =~ s/\|/_/; push @seqID, $seqID }
	else           { $seq{$seqID} .= $_ }
}
close IN;

mkdir "CDD.tmp" unless -e "CDD.tmp";
open COM, '>', "CDD.command" or die $!;
foreach my $query (@seqID) {
	open FASTA, '>', "CDD.tmp/$query.fa" or die $!;
	print FASTA ">$query\n$seq{$query}\n";
	print COM "bwrpsb.pl < CDD.tmp/$query.fa > CDD.tmp/$query.fa.out\n";
}
close COM;

system `ParaFly -c CDD.command -CPU $cpus`;

unless ( -e "FailedCommands" ) {
	open OUT, '>', "CDD.out" or die $!;
	print OUT "Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tIncomplete\tSuperfamily\n";

	foreach (@seqID) {
		open IN, "CDD.tmp/$_.fa.out" or die $!;
		while (<IN>) {
			if (s/^Q#\d+ - >//) { print OUT }
		}
	}

	system `rm -rf CDD.tmp`;
	unlink "CDD.command";
	unlink "CDD.command.completed";
}
