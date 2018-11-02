#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 fastaFile piecesNum[10]

USAGE
if (@ARGV==0) { die $usage }

my $pieces;
if ($ARGV[1]) { $pieces = $ARGV[1] - 1 }
else          { $pieces = 10 }

open FASTA, '<', $ARGV[0] or die $!;
my ($seqID, @seqID, %seq);
while (<FASTA>) {
	chomp;
	if (/^>(.*)/) { $seqID = $1 ; push @seqID, $seqID }
	else          { $seq{$seqID} .= $_ }
}
close FASTA;

while (@seqID) {
	foreach (0..$pieces) {
		open OUTPUT, '>>', "$_.fasta" or die $!;
		$seqID = shift @seqID;
		my $seq = $seq{$seqID};
		my $handle="FASTA$_";
		print OUTPUT ">$seqID\n$seq\n";
		last unless @seqID;
		close OUTPUT;
	 }
}
