#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl locus_of_NNN_in_genome.pl genome.fasta > nnn.gff3

USAGE
if (@ARGV==0) { die $usage }

open IN, '<', $ARGV[0] or die $!;
my ($seqID, @seqID, %seq);
while (<IN>) {
	chomp;
	if (/^>(\S+)/) { $seqID = $1; push @seqID, $1 }
	else           { $seq{$seqID} .= $_ }
}
close IN;
print STDERR "Fasta File Read Over!\n";

my $num = 0;
foreach my $seqId(@seqID) {
	my $seq = $seq{$seqId};
	my $position = 1;
	while ($seq) {
		if ( $seq =~ s/^([^Nn]+)// ) {
			$position += length $1;
		}
		elsif ( $seq =~ s/^([Nn]+)// ) {
			my $next_position = $position + ( length $1 ) - 1 ;
			$num ++;
			my $length = length $1;
			print "$seqId\tGAP\tgap\t$position\t$next_position\t.\t.\t.\tID=$num;Length=$length\n";
			$position = $next_position + 1 ;
		}
	}
}

